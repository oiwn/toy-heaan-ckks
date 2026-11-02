//! Compare a bounded recurrence with and without explicit bootstrap.
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

fn run<const N: usize>(
    primes_count: usize,
    refresh: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let started = Instant::now();
    let primes = generate_primes(61, primes_count, N as u64);
    let context = Arc::new(RnsBasis::<N>::new(primes)?);
    let modulus_bits = context
        .moduli()
        .iter()
        .fold(BigUint::from(1u8), |product, prime| product * prime)
        .bits();
    let scale = Scale::power_of_two(44)?;
    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        context.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale,
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(2026);
    let secret = engine.generate_secret_key(&mut rng)?;
    let public = engine.generate_public_key(&secret, &mut rng)?;
    let recurrence = RecurrencePlan::new(context.clone())?;
    let keys = recurrence.generate_keys(&secret, 3.2f64.sqrt(), &mut rng)?;
    let bootstrap = if refresh {
        Some(BootstrapPlan::new(
            context.clone(),
            BootstrapConfig::reference(N)?,
        )?)
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
        "\nN={N} primes={primes_count} modulus_bits={} bootstrap={refresh} setup={:?}",
        modulus_bits,
        started.elapsed()
    );
    let mut compute = std::time::Duration::ZERO;
    let mut diagnostics = std::time::Duration::ZERO;
    let mut refreshes = 0;
    let mut max_error = 0f64;
    let mut completed = 0;
    for round in 0..32 {
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
    let target = if primes_count == 18 && !refresh {
        17
    } else {
        32
    };
    if completed != target || refreshes != if refresh { 3 } else { 0 } {
        return Err("unexpected depth or refresh schedule".into());
    }
    println!(
        "completed={completed}/32 refreshes={refreshes} max_error={max_error:.3e} evaluation={compute:?} diagnostics={diagnostics:?}"
    );
    Ok(())
}

fn compare<const N: usize>() -> Result<(), Box<dyn std::error::Error>> {
    run::<N>(18, true)?;
    run::<N>(18, false)?;
    run::<N>(34, false)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    match args
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .as_slice()
    {
        [] | ["--degree", "16"] => compare::<16>(),
        ["--degree", "64"] => compare::<64>(),
        _ => Err("usage: recurrence_bootstrap_demo [--degree 16|64]".into()),
    }
}
