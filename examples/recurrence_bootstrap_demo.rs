//! A 32-round recurrence, with and without bootstrapping.
//!
//! ```bash
//! cargo run --release --example recurrence_bootstrap_demo -- [--degree 16|32|64|128]
//! ```
//!
//! A narrated comparison: a recurrence that spends one level per round
//! runs on a short chain with refreshes, on the same chain without them
//! (it stops), and on a chain long enough for all rounds. Expected round and
//! refresh counts are derived from the plans' levels and checked.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;
use std::time::{Duration, Instant};

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

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{
    bad, check, confirm, dim, elapsed, formula, heading, ok, param, say, step,
    subheading,
};

type Error = Box<dyn std::error::Error>;

const ROUNDS: usize = 32;
/// Chain levels above the bootstrap depth: room for the recurrence's
/// rounds between refreshes (18 primes at N=16/64).
const RECURRENCE_LEVELS: usize = 7;
const BOUND: f64 = 1e-3;
const SEED: u64 = 2026;

fn chain<const N: usize>(primes_count: usize) -> Result<Arc<RnsBasis<N>>, Error> {
    Ok(Arc::new(RnsBasis::<N>::new(generate_primes(
        61,
        primes_count,
        N as u64,
    ))?))
}

/// What one run did.
struct Outcome {
    lowest: usize,
    passed: bool,
}

/// Runs the recurrence on `primes_count` primes, with or without refresh.
fn run<const N: usize>(
    title: &str,
    config: &BootstrapConfig,
    scale: Scale,
    primes_count: usize,
    refresh: bool,
) -> Result<Outcome, Error> {
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
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
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

    println!("\n{}", subheading(title));
    step(
        "chain",
        format!(
            "{primes_count} primes, {modulus_bits} bits, refresh below level {}",
            recurrence.refresh_level()
        ),
    );
    println!(
        "  {}",
        dim(&format!("setup {}", elapsed(started.elapsed())))
    );

    let decode_error = |state: &_, expected: &[f64]| -> Result<f64, Error> {
        let actual = encoder.decode(&CkksEngine::decrypt(
            state,
            &secret.mod_switch_to(state.level())?,
        ));
        if actual.iter().any(|x| !x.is_finite()) {
            return Ok(f64::INFINITY);
        }
        Ok(actual
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max))
    };
    let mut compute = Duration::ZERO;
    let mut refreshes = 0;
    let mut max_error = 0f64;
    let mut completed = 0;
    // One entry per round: the level after it, "R" before a refresh.
    let mut trace = Vec::new();
    let mut stopped = None;
    for round in 0..ROUNDS {
        if recurrence.needs_refresh(&state) {
            let Some(plan) = &bootstrap else {
                stopped = Some((round + 1, state.level()));
                break;
            };
            let started = Instant::now();
            state = plan.bootstrap(
                &state,
                bootstrap_keys.as_ref().ok_or("missing bootstrap keys")?,
            )?;
            compute += started.elapsed();
            refreshes += 1;
            trace.push(format!("R->{}", state.level()));
            max_error = max_error.max(decode_error(&state, &expected)?);
        }
        let started = Instant::now();
        state = recurrence.step(&state, &driver, round, &keys)?;
        compute += started.elapsed();
        for (x, u) in expected.iter_mut().zip(&values) {
            *x = 0.5 * *x
                + 0.125 * *x * *x
                + if round % 2 == 0 { 0.25 * u } else { -0.25 * u };
        }
        max_error = max_error.max(decode_error(&state, &expected)?);
        completed += 1;
        trace.push(state.level().to_string());
    }
    for (index, part) in trace.chunks(16).enumerate() {
        let label = if index == 0 { "levels:" } else { "       " };
        println!("  {} {}", formula(label), part.join(" "));
    }
    if let Some((round, level)) = stopped {
        println!(
            "  {}",
            bad(&format!(
                "stopped before round {round}: level {level} cannot hold another step"
            ))
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
    step(
        "rounds, refreshes",
        format!(
            "{completed}/{ROUNDS}, {refreshes}   evaluation {}",
            elapsed(compute)
        ),
    );
    let mut passed = confirm(
        &format!(
            "matches the derived {want_rounds} rounds, {want_refreshes} refreshes"
        ),
        completed == want_rounds && refreshes == want_refreshes,
    );
    passed &= check("every completed round", max_error, BOUND);
    Ok(Outcome { lowest, passed })
}

fn compare<const N: usize>() -> Result<bool, Error> {
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

    println!("{}", heading("Parameters"));
    param("N              ", N, &format!("{} real slots", N / 2));
    param(
        "message scale  ",
        format!("2^{:.0}", scale.log2()),
        "the bootstrap's; the recurrence runs at it",
    );
    param("bootstrap depth", depth, "levels one refresh spends");
    param(
        "short chain    ",
        short,
        &format!("depth + {RECURRENCE_LEVELS} levels for rounds"),
    );
    param(
        "seed           ",
        SEED,
        "fixed, so every run prints the same numbers",
    );

    let mut passed =
        run::<N>("1. short chain, with refresh", &config, scale, short, true)?
            .passed;
    let without =
        run::<N>("2. short chain, no refresh", &config, scale, short, false)?;
    passed &= without.passed;
    let long = ROUNDS + without.lowest;
    passed &=
        run::<N>("3. long chain, no refresh", &config, scale, long, false)?.passed;
    Ok(passed)
}

fn main() -> Result<(), Error> {
    println!(
        "{}\n",
        heading("A long computation: refresh or a long chain")
    );
    say(&format!(
        "Each round computes\n  {}\n\
         (u is a fixed encrypted input), one multiplication, so one level\n\
         per round. 32 rounds need 32 levels. Two ways to get them:",
        formula("x -> x/2 + x^2/8 +- u/4")
    ));
    say(
        "  - a chain with 32+ primes: simple, but every prime makes every\n\
         \x20   operation slower and needs a larger ring for security;\n  \
         - a short chain plus bootstrapping: refresh whenever the levels\n\
         \x20   run low. The levels trace below shows the sawtooth (R = refresh).\n",
    );
    let args: Vec<_> = std::env::args().skip(1).collect();
    let passed = match args
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
    }?;

    println!("\n{}", heading("Takeaways"));
    say(
        "Without refresh the short chain stops early; with it the same chain\n\
         finishes all 32 rounds, trading bootstrap time for a small modulus.\n\
         The long chain finishes too, but every operation on it is slower.",
    );
    if !passed {
        return Err("a run missed its schedule or accuracy target".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
