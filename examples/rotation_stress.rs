//! How rotation noise grows over hundreds of rotations.
//!
//! ```bash
//! cargo run --example rotation_stress
//! ```
//!
//! A narrated measurement: rotate one ciphertext by +1 again and again
//! (800 times) and watch the error at checkpoints. Each rotation is one key
//! switch, which adds a small, fresh noise; the noise adds up roughly
//! linearly, while the level and scale never change.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, formula, heading, ok, param, say, table};

const N: usize = 32;
const SLOTS: usize = N / 2;
const SCALE_BITS: i32 = 58;
const SEED: u64 = 42;
const BOUND: f64 = 1e-3;
const CHECKPOINTS: [usize; 8] = [1, 5, 10, 50, 100, 200, 400, 800];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Rotation noise over many rotations"));
    say(&format!(
        "A rotation keeps level and scale, so nothing stops you from rotating\n\
         forever, except noise. Every key switch adds an independent error of\n\
         about digit_base * N * std_dev in the coefficients, so after k\n\
         rotations:\n  {}\n\
         The error grows slowly, roughly linearly in k. Multiplications,\n\
         not rotations, are what exhaust a ciphertext.\n",
        formula("error ~ k * (per-rotation noise) / scale")
    ));

    println!("{}", heading("Parameters"));
    param("N     ", N, "ring degree; 16 slots");
    param("primes", "3 x 30 bits", "rotation spends no level");
    param(
        "scale ",
        format!("2^{SCALE_BITS}"),
        "keeps 800 rotations of noise small",
    );
    param(
        "seed  ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();
    let basis = Arc::new(RnsBasis::<N>::new(generate_primes(30, 3, N as u64))?);
    let scale = Scale::power_of_two(SCALE_BITS)?;
    let engine = CkksEngine::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale,
        },
    );
    let encoder = CkksEncoder::<N>::new(scale);
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let sk = engine.generate_secret_key(&mut rng)?;
    let pk = engine.generate_public_key(&sk, &mut rng)?;
    // One key for offset +1, reused for every rotation.
    let key = engine.generate_rotation_key(&sk, 1, &mut rng);

    let input: Vec<f64> = (1..=SLOTS).map(|i| i as f64).collect();
    let mut ct = engine.encrypt(&encoder.encode(&input, basis), &pk, &mut rng);
    let mut expected = input;

    println!("{}", heading("Rotate by +1, 800 times"));
    let mut rows = Vec::new();
    let mut done = 0;
    let mut worst: f64 = 0.0;
    for checkpoint in CHECKPOINTS {
        for _ in done..checkpoint {
            ct = ct.rotate(&key);
            expected.rotate_left(1);
        }
        done = checkpoint;
        let decoded =
            encoder.decode(&CkksEngine::<RnsPoly<N>, N>::decrypt(&ct, &sk));
        let error = expected
            .iter()
            .zip(&decoded)
            .map(|(e, g)| (e - g).abs())
            .fold(0.0, f64::max);
        worst = worst.max(error);
        rows.push(vec![
            checkpoint.to_string(),
            format!("{error:.3e}"),
            format!("{:.3e}", error / checkpoint as f64),
            format!("{:.6}", decoded[0]),
            format!("{:.0}", expected[0]),
            format!("{}, 2^{:.0}", ct.level(), ct.scale.log2()),
        ]);
    }
    table(
        &[
            "rotations",
            "max error",
            "error / k",
            "slot 0",
            "expected",
            "level, scale",
        ],
        &rows,
    );
    say(
        "\nAt first the encryption noise (~1e-9) dominates, so error / k\n\
         falls. From about 100 rotations it settles near 4e-11: each\n\
         rotation adds about the same noise, and the total grows linearly.\n\
         Level and scale never move.\n",
    );
    let passed = check("after 800 rotations", worst, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(
        "Rotations are cheap in depth but not free in noise: a large scale\n\
         buys room for many of them. Next: horner_chain, where\n\
         multiplications spend levels.",
    );
    if !passed {
        return Err("rotation noise above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
