//! Key switching: move a ciphertext from one secret key to another.
//!
//! ```bash
//! cargo run --example key_switch
//! ```
//!
//! A narrated walk-through: why key switching exists, how a gadget
//! decomposition keeps its noise small (with a decimal toy), one switch at
//! N=16 with the noise it adds, and what decrypting with the wrong key
//! looks like. Relinearization, rotation and conjugation are this one
//! operation with different old keys.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::PolyRing;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::{DIGIT_BITS, KeySwitchingKey};
use toy_heaan_ckks::ckks::{Ciphertext, Plaintext, Scale};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{
    bad, check, dim, formula, heading, ok, param, say, step, subheading,
};

const N: usize = 16;
const SCALE_BITS: i32 = 50;
const SEED: u64 = 3;
const BOUND: f64 = 1e-6;

fn max_error(actual: &[f64], expected: &[f64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Key switching: from s_old to s_new"));
    say(&format!(
        "A ciphertext (c0, c1) decrypts as {}. Several\n\
         operations leave it under a different key: multiplication produces\n\
         a term in s^2, rotation one in s(X^(5^k)), conjugation one in\n\
         s(X^(2N-1)). Key switching turns c1 * s_old into a pair that\n\
         decrypts under s_new, without ever decrypting.",
        formula("m ~ c0 + c1 * s_old")
    ));
    say(&format!(
        "\nThe switching key encrypts s_old under s_new. Multiplying c1 by it\n\
         directly would multiply its noise by c1, a number as large as the\n\
         modulus. So c1 is first split into small digits, and the key holds\n\
         one encryption of {} per digit position.\n",
        formula("s_old * base^j")
    ));

    println!(
        "{}",
        subheading("Toy in decimal: base 100, key noise e = 3")
    );
    let d: u64 = 734_512;
    let digits = [d % 100, d / 100 % 100, d / 10_000];
    step("d                              ", d);
    step(
        "digits d_j (base 100)          ",
        format!("{} + {}*100 + {}*100^2", digits[0], digits[1], digits[2]),
    );
    step("noise without digits:  d * e   ", d * 3);
    step(
        "noise with digits:  sum(d_j) * e",
        digits.iter().sum::<u64>() * 3,
    );
    println!(
        "  {}\n",
        dim("each digit is below the base, so the noise stays near base * e")
    );

    println!("{}", heading("Parameters"));
    let basis = Arc::new(RnsBasis::<N>::new(generate_primes(61, 3, N as u64))?);
    let scale = Scale::power_of_two(SCALE_BITS)?;
    param("N         ", N, "ring degree; 8 real slots");
    param("primes    ", "3 x 61 bits", "logq 180; no level is used");
    param(
        "DIGIT_BITS",
        DIGIT_BITS,
        "base 2^20: each 61-bit residue gives 4 digits",
    );
    param(
        "scale     ",
        format!("2^{SCALE_BITS}"),
        "well above the switch noise",
    );
    param(
        "seed      ",
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
    let s_old = engine.generate_secret_key(&mut rng)?;
    let s_new = engine.generate_secret_key(&mut rng)?;
    let public_old = engine.generate_public_key(&s_old, &mut rng)?;
    let ksk =
        KeySwitchingKey::generate(&s_new, &s_old.poly, 3.2f64.sqrt(), &mut rng);

    let x: Vec<f64> = (0..N / 2).map(|j| (j as f64 * 0.9).sin()).collect();
    let ct = engine.encrypt(&encoder.encode(&x, basis), &public_old, &mut rng);

    println!("{}", heading("Step 1: switch"));
    say(&format!(
        "Split c1 into digits and multiply each by its key entry (a_k, b_k):\n  \
         {}\n  {}\n\
         Then c0' + c1' * s_new ~ c0 + c1 * s_old: the same message.",
        formula("c0' = c0 + sum_k digit_k * b_k"),
        formula("c1' =      sum_k digit_k * a_k")
    ));
    step("gadget terms (3 primes x 4 digits)", ksk.a.len());
    let (extra, c1) = ksk.switch(&ct.c1);
    let mut c0 = ct.c0.clone();
    c0 += &extra;
    let switched = Ciphertext {
        c0,
        c1,
        scale: ct.scale,
    };

    println!("\n{}", heading("Step 2: decrypt with each key"));
    let before = CkksEngine::decrypt(&ct, &s_old);
    let after = CkksEngine::decrypt(&switched, &s_new);
    let wrong_new = CkksEngine::decrypt(&ct, &s_new);
    let wrong_old = CkksEngine::decrypt(&switched, &s_old);
    let rows: [(&str, &Plaintext<RnsPoly<N>, N>, bool); 4] = [
        ("original  with s_old", &before, true),
        ("original  with s_new", &wrong_new, false),
        ("switched  with s_new", &after, true),
        ("switched  with s_old", &wrong_old, false),
    ];
    for (name, pt, right_key) in rows {
        let error = format!("{:>10.3e}", max_error(&encoder.decode(pt), &x));
        let error = if right_key { ok(&error) } else { bad(&error) };
        println!("  {} max error {error}", formula(name));
    }
    say(&format!(
        "\n{} a wrong key leaves c1 * (s_old - s_new), a random\n\
         polynomial as large as the modulus, so the decoded values are\n\
         about modulus / scale = 2^130.",
        bad("Shown on purpose:")
    ));

    println!("\n{}", heading("Step 3: the noise the switch adds"));
    let mut noise = after.poly.clone();
    noise -= &before.poly;
    let added = noise
        .to_centered_f64()
        .iter()
        .fold(0.0f64, |m, c| m.max(c.abs()));
    step(
        "max |added noise| per coefficient",
        format!("2^{:.1}", added.log2()),
    );
    step(
        "digit base                       ",
        format!("2^{DIGIT_BITS}"),
    );
    step(
        "after dividing by the scale      ",
        format!("2^{:.1}", added.log2() - scale.log2()),
    );
    say(
        "\nThe noise follows the digit base (2^20, times the number of digits,\n\
         N and the key noise), not the 180-bit modulus. Smaller digits mean\n\
         less noise but more key entries and more multiplications.\n",
    );
    let passed = check(
        "switched under s_new",
        max_error(&encoder.decode(&after), &x),
        BOUND,
    );

    println!("\n{}", heading("Takeaways"));
    say(
        "One operation, many uses: relinearization switches from s^2,\n\
         rotation from s(X^(5^k)), conjugation from s(X^(2N-1)). Next:\n\
         rotation_demo and conjugation use it to move slots.",
    );
    println!(
        "{}",
        dim("Reference: fhetextbook.github.io, key switching (D-3.8).")
    );
    if !passed {
        return Err("key switch error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
