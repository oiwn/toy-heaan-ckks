//! Multiply two ciphertexts: tensor, relinearize, rescale.
//!
//! ```bash
//! cargo run --example encrypt_mul
//! ```
//!
//! A narrated walk-through of ciphertext multiplication: why the product
//! decrypts under s^2, how the relinearization key brings it back to s,
//! and why the scale doubles and a rescale divides it back down.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::RelinearizationKey;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, formula, heading, ok, param, say, step, value};

const N: usize = 16;
const SCALE_BITS: i32 = 30;
const SEED: u64 = 42;
const BOUND: f64 = 1e-4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Homomorphic multiplication"));
    say(&format!(
        "Multiplying (c0 + c1*s)(d0 + d1*s) gives three terms:\n  {}\n\
         The result decrypts with 1, s and s^2: a ciphertext of degree 2.\n\
         Left alone, every product would add a degree.",
        formula("c0*d0 + (c0*d1 + c1*d0)*s + c1*d1*s^2")
    ));
    say(&format!(
        "\nTwo fixes keep ciphertexts small:\n  {}  key-switch the s^2 term back to s\n  \
         {}      the scale became scale^2; divide by the last\n               \
         prime (about scale) and drop it: one level\n",
        formula("relinearize"),
        formula("rescale")
    ));

    println!("{}", heading("Parameters"));
    let scale = Scale::power_of_two(SCALE_BITS)?;
    param("N     ", N, "ring degree; 8 slots, 4 used");
    param(
        "scale ",
        format!("2^{SCALE_BITS}"),
        "the product sits at 2^60 before rescale",
    );
    param("primes", "4 x 31 bits", "one is dropped by the rescale");
    param(
        "seed  ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();
    let basis = Arc::new(RnsBasis::<N>::new(generate_primes(31, 4, N as u64))?);
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
    let sk = engine.generate_secret_key(&mut rng)?;
    let pk = engine.generate_public_key(&sk, &mut rng)?;
    let rlk: RelinearizationKey<RnsPoly<N>, N> =
        engine.generate_relin_key(&sk, &mut rng);

    let a = [1.0, 2.0, 3.0, 4.0];
    let b = [0.5, 1.0, 1.5, 2.0];
    println!("{}", heading("Inputs"));
    step("a", format!("{a:?}"));
    step("b", format!("{b:?}"));
    let ct_a = engine.encrypt(&encoder.encode(&a, basis.clone()), &pk, &mut rng);
    let ct_b = engine.encrypt(&encoder.encode(&b, basis), &pk, &mut rng);
    let report = |label: &str, level: usize, logq: u32, scale_bits: f64| {
        step(
            label,
            format!("level {level}, logq {logq}, scale 2^{scale_bits:.1}"),
        );
    };
    report(
        "encrypted          ",
        ct_a.level(),
        ct_a.logq(),
        ct_a.scale.log2(),
    );

    println!("\n{}", heading("Step 1: multiply and relinearize"));
    say(
        "`mul` forms the three terms and key-switches c1*d1 with the\n\
         relinearization key, so the result is again a pair (c0, c1).",
    );
    let product = ct_a.mul(&ct_b, &rlk);
    report(
        "after mul          ",
        product.level(),
        product.logq(),
        product.scale.log2(),
    );

    println!("\n{}", heading("Step 2: rescale"));
    say(
        "Divide by the last prime and drop it. The 31-bit primes are a bit\n\
         above 2^30, so the scale lands just below 2^30. CKKS tracks the\n\
         exact value, so decoding stays correct.",
    );
    let rescaled = product.rescale()?;
    report(
        "after rescale      ",
        rescaled.level(),
        rescaled.logq(),
        rescaled.scale.log2(),
    );

    println!("\n{}", heading("Step 3: decrypt and decode"));
    say("The key is cut to the ciphertext's level before decrypting.");
    let sk_low = sk.mod_switch_to(rescaled.level())?;
    let decoded = encoder.decode(&CkksEngine::decrypt(&rescaled, &sk_low));
    println!(
        "\n  {:>4} | {:>6} | {:>6} | {:>7} | {:>12} | {:>9}",
        "slot", "a", "b", "a * b", "decrypted", "error"
    );
    println!("  {}", "-".repeat(58));
    let mut worst: f64 = 0.0;
    for (j, ((x, y), got)) in a.iter().zip(&b).zip(&decoded).enumerate() {
        let error = (x * y - got).abs();
        worst = worst.max(error);
        println!(
            "  {j:>4} | {x:>6.2} | {y:>6.2} | {:>7.3} | {} | {error:>9.2e}",
            x * y,
            value(format!("{got:>12.8}"))
        );
    }
    let passed = check("a * b", worst, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(
        "Each multiplication costs one prime (one level); the chain length\n\
         is the multiplicative depth. Noise grows with each product, and\n\
         relinearization adds key-switch noise. Next: mod_switch for levels,\n\
         key_switch for what relinearization does inside.",
    );
    if !passed {
        return Err("multiplication error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
