//! Encrypt two vectors, add them without decrypting, decrypt the sum.
//!
//! ```bash
//! cargo run --example encrypt_add
//! ```
//!
//! A narrated walk-through of the full CKKS pipeline for the simplest
//! operation: encode -> encrypt -> add -> decrypt -> decode. Shows what
//! encryption does to a plaintext, why adding ciphertexts adds the slots,
//! and where the small error comes from.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::PolyRing;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, dim, formula, heading, ok, param, say, step, value};

const N: usize = 16;
const SCALE_BITS: i32 = 30;
const SEED: u64 = 42;
const ERROR_VARIANCE: f64 = 3.2;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Homomorphic addition"));
    say(&format!(
        "Encryption hides an encoded message m behind a random mask that\n\
         only the secret key s removes:\n  {}\n  {}\n\
         c1 is random-looking; c0 hides m. Adding two ciphertexts adds\n\
         both parts, so the sum decrypts to m_a + m_b plus both errors.\n\
         Addition needs no key and consumes no level.\n",
        formula("encrypt: (c0, c1) with c0 + c1 * s = m + e"),
        formula("decrypt: m + e = c0 + c1 * s")
    ));

    println!("{}", heading("Parameters"));
    let scale = Scale::power_of_two(SCALE_BITS)?;
    param("N     ", N, "ring degree; 8 slots, 4 used");
    param(
        "scale ",
        format!("2^{SCALE_BITS}"),
        "error about 1e-7 after decoding",
    );
    param(
        "primes",
        "3 x 31 bits",
        "Q ~ 2^93; addition needs no spare level",
    );
    param(
        "seed  ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();
    let basis = Arc::new(RnsBasis::<N>::new(generate_primes(31, 3, N as u64))?);
    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        basis.clone(),
        CkksParams {
            error_variance: ERROR_VARIANCE,
            hamming_weight: N / 2,
            scale,
        },
    );
    let encoder = CkksEncoder::<N>::new(scale);
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let sk = engine.generate_secret_key(&mut rng)?;
    let pk = engine.generate_public_key(&sk, &mut rng)?;

    println!("{}", heading("Step 1: encode and encrypt"));
    let a = [1.0, 2.0, 3.0, 4.0];
    let b = [0.5, -1.5, 0.25, -0.75];
    step("a", format!("{a:?}"));
    step("b", format!("{b:?}"));
    let pt_a = encoder.encode(&a, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);
    let ct_b = engine.encrypt(&encoder.encode(&b, basis), &pk, &mut rng);
    step("plaintext  m_a[0]", pt_a.poly.to_coeffs()[0]);
    step(
        "ciphertext c0[0] ",
        format!("~2^{:.0}", ct_a.c0.to_centered_f64()[0].abs().log2()),
    );
    println!(
        "  {}",
        dim("m is about 2^30; c0 looks like a random number mod Q ~ 2^93")
    );
    step(
        "level, scale",
        format!("{}, 2^{:.0}", ct_a.level(), ct_a.scale.log2()),
    );

    println!("\n{}", heading("Step 2: add the ciphertexts"));
    say(&format!(
        "{}  (no key used)",
        formula("(c0_a + c0_b, c1_a + c1_b)")
    ));
    let ct_sum = ct_a.add(&ct_b)?;
    step(
        "level, scale",
        format!(
            "{}, 2^{:.0}  unchanged",
            ct_sum.level(),
            ct_sum.scale.log2()
        ),
    );

    println!("\n{}", heading("Step 3: decrypt and decode"));
    let decoded = encoder.decode(&CkksEngine::decrypt(&ct_sum, &sk));
    println!(
        "  {:>4} | {:>8} | {:>8} | {:>9} | {:>12} | {:>9}",
        "slot", "a", "b", "a + b", "decrypted", "error"
    );
    println!("  {}", "-".repeat(64));
    let mut worst: f64 = 0.0;
    for (j, ((x, y), got)) in a.iter().zip(&b).zip(&decoded).enumerate() {
        let error = (x + y - got).abs();
        worst = worst.max(error);
        println!(
            "  {j:>4} | {x:>8.4} | {y:>8.4} | {:>9.4} | {} | {error:>9.2e}",
            x + y,
            value(format!("{got:>12.8}"))
        );
    }
    // The e1 * s term dominates: about std_dev * sqrt(weight * N) per
    // coefficient; 10x margin for 4 slots, plus encoding rounding.
    // `encrypt` samples with 3.2 as the std_dev (a known quirk, see
    // specs/ctx.md), so the bound uses 3.2 directly.
    let std_dev = ERROR_VARIANCE;
    let bound = 10.0 * std_dev * ((N / 2 * N) as f64).sqrt() / scale.value()
        + 4.0 / scale.value();
    say(&format!(
        "\nThe error is the encryption noise divided by the scale:\n  {}",
        formula("~ 10 * std_dev * sqrt(weight * N) / scale")
    ));
    let passed = check("a + b", worst, bound);

    println!("\n{}", heading("Takeaways"));
    say(
        "Addition is free: no key, no level, and the noise only adds up.\n\
         Every result carries a small error, so CKKS is approximate\n\
         arithmetic. Next: encrypt_mul, where noise and scale grow.",
    );
    if !passed {
        return Err("addition error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
