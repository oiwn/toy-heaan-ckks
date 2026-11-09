//! Bootstrap, multiply, bootstrap again: a ciphertext that never runs out.
//!
//! ```bash
//! cargo run --release --example bootstrap_demo
//! ```
//!
//! A narrated walk-through of bootstrapping as a black box at N=16 and
//! N=64: encrypt complex slots, square them (spending a level), refresh,
//! and repeat. Each refresh must keep the slots within 1e-3.
//! `bootstrap_stages` opens the box stage by stage.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use num_complex::Complex64;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::RelinearizationKey;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, dim, formula, heading, ok, param, say, subheading};

const PRIMES: usize = 18;
const CYCLES: usize = 2;
const BOUND: f64 = 1e-3;
const SEED: u64 = 0;

/// Runs the refresh cycles at degree `N`; returns whether all passed.
fn run<const N: usize>() -> Result<bool, Box<dyn std::error::Error>> {
    let basis =
        Arc::new(RnsBasis::<N>::new(generate_primes(61, PRIMES, N as u64))?);
    let config = BootstrapConfig::reference(N)?;
    let plan = BootstrapPlan::new(basis.clone(), config.clone())?;
    let scale = plan.schedule().message_scale();
    println!(
        "\n{}",
        subheading(&format!("N = {N}, {} complex slots", N / 2))
    );
    param(
        "message scale",
        format!("2^{:.0}", scale.log2()),
        "q0 / 2^17: keeps m / q0 small",
    );
    param(
        "secret weight",
        config.max_secret_weight,
        "nonzeros in s; bounds the ModRaise overflow I",
    );
    param(
        "depth        ",
        plan.schedule().depth(),
        "levels one refresh spends",
    );
    param(
        "output level ",
        plan.schedule().output_level(),
        "where a refreshed ciphertext lands",
    );
    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: config.max_secret_weight,
            scale,
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let secret = engine.generate_secret_key(&mut rng)?;
    let public = engine.generate_public_key(&secret, &mut rng)?;
    let keys = plan.generate_keys(&secret, 3.2f64.sqrt(), &mut rng)?;
    let encoder = CkksEncoder::new(scale);
    let mut expected: Vec<_> = (0..N / 2)
        .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
        .collect();
    let initial = engine.encrypt(
        &encoder.encode_complex(&expected, basis),
        &public,
        &mut rng,
    );
    let relin = RelinearizationKey::generate(&secret, 3.2f64.sqrt(), &mut rng);
    let mut squared = initial.mul(&initial, &relin).rescale()?;
    for x in &mut expected {
        *x *= *x;
    }
    let mut passed = true;
    for cycle in 0..CYCLES {
        let before = squared.level();
        let refreshed = plan.bootstrap(&squared, &keys)?;
        let sk = secret.mod_switch_to(refreshed.level())?;
        let actual = encoder.decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
        let error = actual
            .iter()
            .zip(&expected)
            .map(|(a, b)| (*a - *b).norm())
            .fold(0.0, f64::max);
        println!(
            "  {}",
            formula(&format!(
                "cycle {cycle}: z^{} at level {before} -> bootstrap -> level {}, scale 2^{:.1}",
                1 << (cycle + 1),
                refreshed.level(),
                refreshed.scale.log2()
            ))
        );
        passed &= check(&format!("  cycle {cycle} slots"), error, BOUND);
        let relin = RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
        squared = refreshed.mul(&refreshed, &relin).rescale()?;
        for x in &mut expected {
            *x *= *x;
        }
    }
    Ok(passed)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Bootstrapping: refresh a spent ciphertext"));
    say(
        "Every multiplication spends a prime of the modulus chain. Without a\n\
         way to get levels back, the chain length caps the computation.\n\
         Bootstrapping evaluates the decryption circuit homomorphically and\n\
         returns the same slots at a high level, using only public keys.",
    );
    say(&format!(
        "\nThis demo squares the slots and refreshes, twice:\n  {}\n\
         |z| = 0.8, so z^2 and z^4 stay below 1, as bootstrapping requires.\n",
        formula("z -> z^2 -> bootstrap -> z^4 -> bootstrap")
    ));
    println!("{}", heading("Parameters"));
    param(
        "primes       ",
        "18 x 61 bits",
        "11 levels for the refresh, the rest to compute",
    );
    param(
        "config       ",
        "reference(N)",
        "polynomial 31, double angles, 1-level transforms",
    );
    param(
        "seed         ",
        SEED,
        "fixed, so every run prints the same numbers",
    );

    let mut passed = run::<16>()?;
    passed &= run::<64>()?;

    println!("\n{}", heading("Takeaways"));
    say(
        "A refresh costs 11 levels at these sizes and leaves 6 for new work,\n\
         so the computation can go on indefinitely. The refreshed error\n\
         (~1e-5 to 1e-4) comes from the sine approximation and key-switch\n\
         noise. Next: bootstrap_stages for the inside,\n\
         recurrence_bootstrap_demo for a long computation.",
    );
    println!(
        "{}",
        dim("References: CKKS bootstrapping (eprint 2018/153); OpenFHE, Lattigo.")
    );
    if !passed {
        return Err("bootstrap accuracy target failed".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
