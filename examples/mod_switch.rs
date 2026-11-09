//! Modulus switching versus rescaling, and aligning levels before adding.
//!
//! ```bash
//! cargo run --example mod_switch
//! ```
//!
//! A narrated walk-through: what a level is, how `mod_switch_to` drops
//! primes without touching the message, how `rescale` also divides the
//! scale after a multiplication, what goes wrong when a product is only
//! mod-switched, and why ciphertexts must share a level before `add`.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::SecretKey;
use toy_heaan_ckks::ckks::{Ciphertext, Scale};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{bad, check, dim, formula, heading, ok, param, say, value};

const N: usize = 16;
const SCALE_BITS: i32 = 40;
const SEED: u64 = 2;
const BOUND: f64 = 1e-6;

type Ct = Ciphertext<RnsPoly<N>, N>;
type Sk = SecretKey<RnsPoly<N>, N>;
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn decode(encoder: &CkksEncoder<N>, ct: &Ct, sk: &Sk) -> Result<Vec<f64>> {
    let sk = sk.mod_switch_to(ct.level())?;
    Ok(encoder.decode(&CkksEngine::decrypt(ct, &sk)))
}

fn max_error(actual: &[f64], expected: &[f64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

/// One row: name, level, logq, scale, error (green within the bound,
/// red above it).
fn row(name: &str, ct: &Ct, error: f64) {
    let error_text = format!("{error:>9.2e}");
    let error_text = if error <= BOUND {
        ok(&error_text)
    } else {
        bad(&error_text)
    };
    println!(
        "  {} | {:>5} | {:>4} | {:>6} | {}",
        formula(&format!("{name:<26}")),
        ct.level(),
        ct.logq(),
        value(format!("2^{:<4.1}", ct.scale.log2())),
        error_text
    );
}

fn header() {
    println!(
        "  {:<26} | {:>5} | {:>4} | {:>6} | {:>9}",
        "ciphertext", "level", "logq", "scale", "max error"
    );
    println!("  {}", "-".repeat(64));
}

fn main() -> Result<()> {
    println!("{}\n", heading("Levels: mod switch, rescale, match_level"));
    say(
        "A ciphertext lives modulo Q = q0 * q1 * ... * qL, a product of\n\
         primes. Its level is how many primes it still carries. Two\n\
         operations drop a prime, and they are easy to confuse:",
    );
    say(&format!(
        "  {}  forget the last prime; message and scale unchanged\n  \
         {}     divide the message by the dropped prime:\n              \
         {}",
        formula("mod switch"),
        formula("rescale"),
        formula("scale -> scale / q_last")
    ));
    say(
        "\nMultiplication squares the scale; rescale brings it back. Mod switch\n\
         is for aligning: two ciphertexts must sit at the same level (and\n\
         scale) before they can be added.\n",
    );

    println!("{}", heading("Parameters"));
    let mut primes = generate_primes(60, 1, N as u64);
    primes.extend(generate_primes(40, 4, N as u64));
    let basis = Arc::new(RnsBasis::<N>::new(primes)?);
    param("N     ", N, "ring degree; 8 real slots");
    param(
        "q0    ",
        "60 bits",
        "the last level must hold a 2^40 message",
    );
    param(
        "q1..q4",
        "40 bits",
        "close to the scale, so rescale returns to ~2^40",
    );
    param(
        "scale ",
        format!("2^{SCALE_BITS}"),
        "about 12 decimal digits",
    );
    param(
        "seed  ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();

    let scale = Scale::power_of_two(SCALE_BITS)?;
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
    let relin = engine.generate_relin_key(&secret, &mut rng);

    let x: Vec<f64> = (0..N / 2).map(|j| 0.9 - 0.25 * j as f64).collect();
    let square: Vec<f64> = x.iter().map(|v| v * v).collect();
    let ct = engine.encrypt(&encoder.encode(&x, basis.clone()), &public, &mut rng);
    // Largest error over the rows that are meant to be correct.
    let mut worst: f64 = 0.0;

    println!(
        "{}",
        heading("Step 1: mod switch keeps the message and scale")
    );
    say(
        "Drop primes one by one. The same secret key (cut to the level)\n\
         still decrypts, and the scale never changes.\n",
    );
    header();
    for level in (1..=ct.level()).rev() {
        let switched = ct.mod_switch_to(level)?;
        let error = max_error(&decode(&encoder, &switched, &secret)?, &x);
        worst = worst.max(error);
        row(&format!("x.mod_switch_to({level})"), &switched, error);
    }

    println!("\n{}", heading("Step 2: after a multiplication"));
    say(&format!(
        "x * x has scale 2^80. {} divides by the 40-bit q4 and\n\
         returns to about 2^40; {} keeps 2^80.\n",
        formula("rescale()"),
        formula("mod_switch_to(4)")
    ));
    let product = ct.mul(&ct, &relin);
    let rescaled = product.rescale()?;
    let switched = product.mod_switch_to(rescaled.level())?;
    let low = rescaled.mod_switch_to(1)?;
    let overflow = product.mod_switch_to(1)?;
    header();
    for (name, c) in [
        ("x*x", &product),
        ("(x*x).rescale()", &rescaled),
        ("(x*x).mod_switch_to(4)", &switched),
        ("rescaled, mod_switch_to(1)", &low),
    ] {
        let error = max_error(&decode(&encoder, c, &secret)?, &square);
        worst = worst.max(error);
        row(name, c, error);
    }
    row(
        "x*x, mod_switch_to(1)",
        &overflow,
        max_error(&decode(&encoder, &overflow, &secret)?, &square),
    );
    say(&format!(
        "\n{} the unrescaled product needs 80 bits plus a\n\
         sign, but level 1 has a 60-bit modulus: the message wraps around\n\
         and decodes to garbage. Rescale first, then mod switch.",
        bad("Shown on purpose:")
    ));

    println!("\n{}", heading("Step 3: align levels before adding"));
    say(
        "y is a fresh input at the top level, encoded at the scale x*x\n\
         reached. Adding it to the rescaled product fails: they sit at\n\
         different levels, and levels are never aligned silently.\n",
    );
    let y_encoder = CkksEncoder::<N>::new(rescaled.scale);
    let y = engine.encrypt(&y_encoder.encode(&x, basis), &public, &mut rng);
    match y.add(&rescaled) {
        Err(e) => println!("  {} {}", formula("y.add(x*x):"), bad(&e.to_string())),
        Ok(_) => return Err("add across levels should fail".into()),
    }
    say(&format!(
        "\n{} mod-switches y down to the level of x*x:\n",
        formula("y.match_level(&x*x)")
    ));
    let sum = y.match_level(&rescaled)?.add(&rescaled)?;
    let expected: Vec<f64> = x.iter().zip(&square).map(|(a, b)| a + b).collect();
    header();
    let error = max_error(&decode(&encoder, &sum, &secret)?, &expected);
    row("y.match_level(x*x) + x*x", &sum, error);
    worst = worst.max(error);
    println!();
    let passed = check("all rows except the one shown on purpose", worst, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(
        "Rescale after every multiplication; mod switch only to line levels\n\
         up. Each prime is spent once, so the chain length is the depth\n\
         budget. Next: encrypt_mul (relinearize + rescale), key_switch.",
    );
    println!(
        "{}",
        dim("Reference: fhetextbook.github.io, GLWE modulus switch chapter.")
    );
    if !passed {
        return Err("mod switch error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
