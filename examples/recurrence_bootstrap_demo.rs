//! Compare a bounded recurrence with and without explicit bootstrap.
//!
//! Bootstrap uses `BootstrapConfig::reference(N)`, and the recurrence runs
//! at its message scale. Three runs share the inputs: a short chain (the
//! bootstrap depth plus `RECURRENCE_LEVELS`) with refresh, the same chain
//! without refresh, and a chain long enough for all rounds without refresh.
//! Expected round and refresh counts are derived from the plans' levels.
//!
//! ```text
//! cargo run --release --example recurrence_bootstrap_demo -- [--degree 16|32|64|128]
//! ```
use std::sync::Arc;
use std::time::Instant;

use num_bigint::BigUint;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
use toy_heaan_ckks::ckks::algorithms::recurrence::RecurrencePlan;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

type Error = Box<dyn std::error::Error>;

const ROUNDS: usize = 32;
/// Chain levels above the bootstrap depth: room for the recurrence's
/// rounds between refreshes (18 primes at N=16/64, as before).
const RECURRENCE_LEVELS: usize = 7;

fn chain<const N: usize>(primes_count: usize) -> Result<Arc<RnsBasis<N>>, Error> {
    Ok(Arc::new(RnsBasis::<N>::new(generate_primes(
        61,
        primes_count,
        N as u64,
    ))?))
}

/// Run the recurrence and return the plan's refresh level.
fn run<const N: usize>(
    config: &BootstrapConfig,
    scale: Scale,
    primes_count: usize,
    refresh: bool,
) -> Result<usize, Error> {
    let started = Instant::now();
    let context = chain::<N>(primes_count)?;
    let modulus_bits = context
        .moduli()
        .iter()
        .fold(BigUint::from(1u8), |product, prime| product * prime)
        .bits();
    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        context.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: config.max_secret_weight,
            scale,
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(2026);
    let secret = engine.generate_secret_key(&mut rng)?;
    let public = engine.generate_public_key(&secret, &mut rng)?;
    let recurrence = RecurrencePlan::new(context.clone(), scale)?;
    let keys = recurrence.generate_keys(&secret, 3.2f64.sqrt(), &mut rng)?;
    let bootstrap = if refresh {
        Some(BootstrapPlan::new(context.clone(), config.clone())?)
    } else {
        None
    };
    let bootstrap_keys = bootstrap
        .as_ref()
        .map(|p| p.generate_keys(&secret, 3.2f64.sqrt(), &mut rng))
        .transpose()?;
    let encoder = CkksEncoder::new(scale);
    let values: Vec<_> = (0..N / 2)
        .map(|j| -1.0 + 2.0 * j as f64 / (N / 2 - 1) as f64)
        .collect();
    let driver = engine.encrypt(
        &encoder.encode(&values, context.clone()),
        &public,
        &mut rng,
    );
    let mut state = driver.clone();
    let mut expected = values.clone();
    println!(
        "\nN={N} primes={primes_count} modulus_bits={} scale=2^{:.0} \
         refresh_level={} bootstrap={refresh} setup={:?}",
        modulus_bits,
        scale.log2(),
        recurrence.refresh_level(),
        started.elapsed()
    );
    let mut compute = std::time::Duration::ZERO;
    let mut diagnostics = std::time::Duration::ZERO;
    let mut refreshes = 0;
    let mut max_error = 0f64;
    let mut completed = 0;
    for round in 0..ROUNDS {
        if recurrence.needs_refresh(&state) {
            let Some(plan) = &bootstrap else {
                println!(
                    "stop before round {}: insufficient levels ({})",
                    round + 1,
                    state.level()
                );
                break;
            };
            let started = Instant::now();
            state = plan.bootstrap(
                &state,
                bootstrap_keys.as_ref().ok_or("missing bootstrap keys")?,
            )?;
            compute += started.elapsed();
            refreshes += 1;
            println!("  refresh {refreshes}: level={}", state.level());
            let started = Instant::now();
            let actual = encoder.decode(&CkksEngine::decrypt(
                &state,
                &secret.mod_switch_to(state.level())?,
            ));
            for (a, b) in actual.iter().zip(&expected) {
                let error = (a - b).abs();
                if !error.is_finite() || error > 1e-3 {
                    return Err("refresh accuracy failed".into());
                }
                max_error = max_error.max(error);
            }
            diagnostics += started.elapsed();
        }
        let started = Instant::now();
        state = recurrence.step(&state, &driver, round, &keys)?;
        compute += started.elapsed();
        for (x, u) in expected.iter_mut().zip(&values) {
            *x = 0.5 * *x
                + 0.125 * *x * *x
                + if round % 2 == 0 { 0.25 * u } else { -0.25 * u };
        }
        let started = Instant::now();
        let actual = encoder.decode(&CkksEngine::decrypt(
            &state,
            &secret.mod_switch_to(state.level())?,
        ));
        let error = actual
            .iter()
            .zip(&expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0f64, f64::max);
        if actual.iter().any(|x| !x.is_finite()) || error > 1e-3 {
            return Err("recurrence accuracy failed".into());
        }
        max_error = max_error.max(error);
        diagnostics += started.elapsed();
        completed += 1;
        println!(
            "round={:2} level={:2} max_error={error:.3e}",
            round + 1,
            state.level()
        );
    }
    // A step needs `refresh_level` and drops one level; a refresh returns
    // to the bootstrap output level.
    let lowest = recurrence.refresh_level();
    let first = primes_count + 1 - lowest;
    let (want_rounds, want_refreshes) = match &bootstrap {
        _ if first >= ROUNDS => (ROUNDS, 0),
        None => (first, 0),
        Some(plan) => {
            let per_cycle = plan.schedule().output_level() + 1 - lowest;
            (ROUNDS, (ROUNDS - first).div_ceil(per_cycle))
        }
    };
    if completed != want_rounds || refreshes != want_refreshes {
        return Err(format!(
            "unexpected schedule: {completed} rounds, {refreshes} refreshes; \
             expected {want_rounds}, {want_refreshes}"
        )
        .into());
    }
    println!(
        "completed={completed}/{ROUNDS} refreshes={refreshes} max_error={max_error:.3e} evaluation={compute:?} diagnostics={diagnostics:?}"
    );
    Ok(lowest)
}

fn compare<const N: usize>() -> Result<(), Error> {
    let config = BootstrapConfig::reference(N)?;
    // Bootstrap depth from the config: CtS + EvalMod polynomial + double
    // angles + StC (the schedule derives the same count).
    let depth = config.coeffs_to_slots_levels
        + config.polynomial_degree.ilog2() as usize
        + 1
        + config.double_angles as usize
        + config.slots_to_coeffs_levels;
    let short = depth + RECURRENCE_LEVELS;
    let scale = config.message_scale(chain::<N>(short)?.moduli()[0])?;
    run::<N>(&config, scale, short, true)?;
    let lowest = run::<N>(&config, scale, short, false)?;
    run::<N>(&config, scale, ROUNDS + lowest, false)?;
    Ok(())
}

fn main() -> Result<(), Error> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    match args
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .as_slice()
    {
        [] | ["--degree", "16"] => compare::<16>(),
        ["--degree", "32"] => compare::<32>(),
        ["--degree", "64"] => compare::<64>(),
        ["--degree", "128"] => compare::<128>(),
        _ => Err("usage: recurrence_bootstrap_demo [--degree 16|32|64|128]".into()),
    }
}
