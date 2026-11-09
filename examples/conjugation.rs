//! Complex conjugation of every slot, and the real/imaginary split.
//!
//! ```bash
//! cargo run --example conjugation
//! ```
//!
//! A narrated walk-through: why conjugation needs a key switch, what it does
//! to each slot, and how `(z + conj(z)) / 2` and `(z - conj(z)) / 2` split a
//! complex slot into its real and imaginary parts (bootstrap uses this).
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use num_complex::Complex64;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::ConjugationKey;
use toy_heaan_ckks::ckks::{Ciphertext, Scale};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{
    check, dim, formula, heading, ok, param, say, step, subheading, value,
};

const N: usize = 16;
const SCALE_BITS: i32 = 50;
const SEED: u64 = 1;
const BOUND: f64 = 1e-5;

type Ct = Ciphertext<RnsPoly<N>, N>;

fn max_error(actual: &[Complex64], expected: &[Complex64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Conjugation: conj(z) in every slot"));
    say(
        "CKKS packs N/2 complex numbers into the slots of one ciphertext. Some\n\
         algorithms need the conjugate conj(a + bi) = a - bi of every slot, or\n\
         the real and imaginary parts separately: bootstrapping, for example,\n\
         reduces the real and imaginary parts of its slots one at a time.",
    );
    say(&format!(
        "\nThe ring map {} conjugates every slot at once.\n\
         Applied to a ciphertext (c0, c1) it gives a valid encryption of\n\
         conj(z), but under the key s(X^(2N-1)) instead of s. A conjugation\n\
         key switches it back to s, exactly like a rotation key does for\n\
         rotations. No prime is used: the level and scale stay.\n",
        formula("X -> X^(2N-1)")
    ));

    println!("{}", heading("Parameters"));
    let basis = Arc::new(RnsBasis::<N>::new(generate_primes(61, 2, N as u64))?);
    let scale = Scale::power_of_two(SCALE_BITS)?;
    param("N     ", N, "ring degree; N/2 = 8 complex slots");
    param(
        "primes",
        "2 x 61 bits",
        "no multiplication, so no level is used",
    );
    param(
        "scale ",
        format!("2^{SCALE_BITS}"),
        "large: key switch noise (~2^24) vanishes after dividing",
    );
    param(
        "seed  ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();

    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale,
        },
    );
    let encoder = CkksEncoder::<N>::new(scale);
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let secret = engine.generate_secret_key(&mut rng)?;
    let public = engine.generate_public_key(&secret, &mut rng)?;
    let key = ConjugationKey::generate(&secret, 3.2f64.sqrt(), &mut rng);
    let decrypt = |c: &Ct| encoder.decode_complex(&CkksEngine::decrypt(c, &secret));

    let z: Vec<_> = (0..N / 2)
        .map(|j| Complex64::from_polar(0.9, 0.8 * j as f64))
        .collect();
    let ct = engine.encrypt(&encoder.encode_complex(&z, basis), &public, &mut rng);

    println!("{}", heading("Step 1: conjugate"));
    say("Encrypt z, apply the map, switch the key back to s, then decrypt.");
    let conj = ct.conjugate(&key)?;
    let decoded = decrypt(&conj);
    println!(
        "\n  {:>4} | {:>20} | {:>20}",
        "slot", "z", "decrypted conj(z)"
    );
    println!("  {}", "-".repeat(50));
    for (j, (v, d)) in z.iter().zip(&decoded).enumerate() {
        let left = format!("{v:>20.6}");
        let right = format!("{d:>20.6}");
        println!("  {j:>4} | {} | {}", left, value(right));
    }
    println!(
        "  {}",
        dim(&format!(
            "level {} -> {}, scale 2^{:.0} -> 2^{:.0}: unchanged",
            ct.level(),
            conj.level(),
            ct.scale.log2(),
            conj.scale.log2()
        ))
    );
    let conj_expected: Vec<_> = z.iter().map(|v| v.conj()).collect();
    let mut passed = check("conj(z)", max_error(&decoded, &conj_expected), BOUND);

    println!(
        "\n{}",
        heading("Step 2: split into real and imaginary parts")
    );
    say(&format!(
        "Adding and subtracting the conjugate cancels one part:\n  {}\n  {}",
        formula("(z + conj(z)) / 2 = Re(z)"),
        formula("(z - conj(z)) / 2 = i * Im(z)")
    ));
    say(&format!(
        "\nThe division by 2 costs no level: {} encodes\n\
         the integer 1 at scale 2, so the slots halve while the scale doubles.",
        formula("mul_const(0.5, 2^1)")
    ));
    let half = Scale::power_of_two(1)?;
    let real = ct.add(&conj)?.mul_const(0.5, half);
    let imag = ct.sub(&conj)?.mul_const(0.5, half);
    for (name, c, expected) in [
        (
            "Re(z)    ",
            &real,
            z.iter()
                .map(|v| Complex64::new(v.re, 0.0))
                .collect::<Vec<_>>(),
        ),
        (
            "i*Im(z)  ",
            &imag,
            z.iter().map(|v| Complex64::new(0.0, v.im)).collect(),
        ),
    ] {
        println!(
            "\n{}",
            subheading(&format!(
                "{} at level {}, scale 2^{:.0}",
                name.trim(),
                c.level(),
                c.scale.log2()
            ))
        );
        let decoded = decrypt(c);
        step(&format!("{name}slot 1"), format!("{:.6}", decoded[1]));
        passed &= check(name.trim(), max_error(&decoded, &expected), BOUND);
    }

    println!("\n{}", heading("Takeaways"));
    say(
        "Conjugation is one more automorphism plus key switch, the same\n\
         machinery as rotation (see key_switch, rotation_demo). Next:\n\
         bootstrap_stages uses this split before EvalMod.",
    );
    if !passed {
        return Err("conjugation error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
