//! Rotating the slots of a ciphertext.
//!
//! ```bash
//! cargo run --example rotation_demo
//! ```
//!
//! A narrated walk-through: why slots rotate under the ring map
//! `X -> X^(5^k)`, why that needs a rotation key, and a small pipeline
//! (rotate, add, rotate) checked slot by slot. Rotations sum slots, move
//! data between them, and drive the bootstrap's linear transforms.
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
use narrate::{check, dim, formula, heading, ok, param, say, step, table};

const N: usize = 32;
const SLOTS: usize = N / 2;
const SCALE_BITS: i32 = 58;
const PRIMES: usize = 3;
const SEED: u64 = 42;
const BOUND: f64 = 1e-4;
const FIRST: i32 = 1;
const SECOND: i32 = 2;

fn rotate_vec(v: &[f64], k: i32) -> Vec<f64> {
    let k = k.rem_euclid(v.len() as i32) as usize;
    v[k..].iter().chain(&v[..k]).copied().collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Rotation: shifting all slots at once"));
    say(&format!(
        "The slots of a plaintext are values of m(X) at the roots\n\
         zeta^(5^j). Substituting {} moves the value at root j\n\
         to root j - k: every slot shifts left by k, cyclically. The same\n\
         map applied to (c0, c1) rotates an encrypted vector.",
        formula("X -> X^(5^k)")
    ));
    say(
        "\nThe catch: the result decrypts under s(X^(5^k)), not s. A rotation\n\
         key (one per offset k) key-switches it back to s. Rotation consumes\n\
         no level and keeps the scale; it only adds key-switch noise.\n",
    );

    println!("{}", heading("Parameters"));
    let primes = generate_primes(30, PRIMES, N as u64);
    param("N     ", N, "ring degree; 16 slots");
    param(
        "primes",
        "3 x 30 bits",
        "no multiplication, so no level is spent",
    );
    param(
        "scale ",
        format!("2^{SCALE_BITS}"),
        "far above the key-switch noise",
    );
    param(
        "seed  ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();
    let basis = Arc::new(RnsBasis::<N>::new(primes)?);
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
    let key_first = engine.generate_rotation_key(&sk, FIRST, &mut rng);
    let key_second = engine.generate_rotation_key(&sk, SECOND, &mut rng);

    let input: Vec<f64> = (1..=SLOTS).map(|i| i as f64).collect();
    let ct = engine.encrypt(&encoder.encode(&input, basis), &pk, &mut rng);
    println!("{}", heading("Pipeline"));
    say(&format!(
        "  {}\n  {}\n  {}",
        formula("A: rotate(x, +1)"),
        formula("B: x + A"),
        formula("C: rotate(B, +2)")
    ));
    step("x", format!("{input:?}"));

    let rotated = ct.rotate(&key_first);
    let added = ct.add(&rotated)?;
    let result = added.rotate(&key_second);
    step(
        "level, scale after C",
        format!(
            "{}, 2^{:.0}  (same as x: rotation is level-free)",
            result.level(),
            result.scale.log2()
        ),
    );

    let decoded =
        encoder.decode(&CkksEngine::<RnsPoly<N>, N>::decrypt(&result, &sk));
    let ref_rotated = rotate_vec(&input, FIRST);
    let ref_added: Vec<f64> =
        input.iter().zip(&ref_rotated).map(|(a, b)| a + b).collect();
    let expected = rotate_vec(&ref_added, SECOND);

    println!("\n{}", heading("Slot by slot"));
    let rows: Vec<Vec<String>> = (0..SLOTS)
        .map(|j| {
            vec![
                j.to_string(),
                format!("{:.0}", input[j]),
                format!("{:.0}", ref_rotated[j]),
                format!("{:.0}", ref_added[j]),
                format!("{:.0}", expected[j]),
                format!("{:.6}", decoded[j]),
            ]
        })
        .collect();
    table(&["slot", "x", "A", "B", "C expected", "C decrypted"], &rows);
    let error = expected
        .iter()
        .zip(&decoded)
        .map(|(e, g)| (e - g).abs())
        .fold(0.0, f64::max);
    println!();
    let passed = check("rotate, add, rotate", error, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(
        "Each offset needs its own key; real systems keep keys for powers\n\
         of two and compose them. Rotate + add repeatedly sums all slots in\n\
         log2(slots) steps. Next: rotation_stress (noise over many\n\
         rotations), conjugation (the other automorphism).",
    );
    println!(
        "{}",
        dim("Reference: fhetextbook.github.io, homomorphic rotation.")
    );
    if !passed {
        return Err("rotation error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
