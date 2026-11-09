//! Key generation: secret, public and relinearization keys.
//!
//! ```bash
//! cargo run --example keys
//! ```
//!
//! A narrated walk-through of the three CKKS key types and the relation
//! each one satisfies, checked on real keys: the secret is a small ternary
//! polynomial, the public key is an RLWE sample `b + a*s = e`, and the
//! relinearization key is a gadget encryption of `s^2`.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::keys::DIGIT_BITS;
use toy_heaan_ckks::ckks::keys::public_key::{PublicKey, PublicKeyParams};
use toy_heaan_ckks::ckks::keys::relin_key::RelinearizationKey;
use toy_heaan_ckks::ckks::keys::secret_key::{SecretKey, SecretKeyParams};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
use toy_heaan_ckks::rings::traits::PolyRing;

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{confirm, dim, formula, heading, ok, param, say, step};

const N: usize = 16;
const HAMMING_WEIGHT: usize = N / 2;
const ERROR_STD: f64 = 3.2;
const SEED: u64 = 42;
/// A Gaussian with std_dev 3.2 stays far below this over 16 coefficients.
const SMALL: u64 = 50;

fn max_abs(coefficients: &[i64]) -> u64 {
    coefficients
        .iter()
        .map(|c| c.unsigned_abs())
        .max()
        .unwrap_or(0)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Keys: secret, public, relinearization"));
    say(
        "CKKS security rests on Ring-LWE: given a random polynomial a and\n\
         b = -a*s + e with a small error e, finding s is hard. Every key is\n\
         built from that one idea. The checks below use the secret key to\n\
         confirm each relation; a real evaluator never sees it.\n",
    );

    println!("{}", heading("Parameters"));
    let primes = generate_primes(20, 2, N as u64);
    param("N             ", N, "ring degree; toy size, no security");
    param(
        "secret weight ",
        HAMMING_WEIGHT,
        "nonzero coefficients of s",
    );
    param("error std_dev ", ERROR_STD, "width of the Gaussian error e");
    param(
        "primes        ",
        "2 x 20 bits",
        "small, so the residues print",
    );
    param(
        "seed          ",
        SEED,
        "fixed, so every run prints the same keys",
    );
    println!();
    let basis: Arc<RnsBasis<N>> = Arc::new(RnsBasis::new(primes)?);
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let mut passed = true;

    println!("{}", heading("Step 1: secret key s"));
    say(
        "A ternary polynomial: each coefficient is -1, 0 or 1, with exactly\n\
         `weight` nonzeros. Small keys keep products like c1 * s small.",
    );
    let sk_params = SecretKeyParams::<N>::new(HAMMING_WEIGHT)?;
    let sk = SecretKey::<RnsPoly<N>, N>::generate(&sk_params, &basis, &mut rng)?;
    let coefficients = sk.poly.to_coeffs();
    step("s", format!("{coefficients:?}"));
    let weight = coefficients.iter().filter(|&&c| c != 0).count();
    passed &= confirm(
        "every coefficient in {-1, 0, 1}",
        coefficients.iter().all(|c| (-1..=1).contains(c)),
    );
    passed &= confirm(
        &format!("weight {weight} == {HAMMING_WEIGHT}"),
        weight == HAMMING_WEIGHT,
    );

    println!("\n{}", heading("Step 2: public key (a, b)"));
    say(&format!(
        "Sample a uniformly and a small error e, then set\n  {}\n\
         so that {} is small. Encryption multiplies (a, b) by\n\
         a fresh random u; only s can strip the a*s part off again.",
        formula("b = -(a * s) + e"),
        formula("b + a * s = e")
    ));
    let pk_params = PublicKeyParams::<N>::new(ERROR_STD)?;
    let pk =
        PublicKey::<RnsPoly<N>, N>::generate(&sk, &pk_params, &basis, &mut rng)?;
    let mut residual = pk.a.clone();
    residual *= &sk.poly;
    residual += &pk.b;
    let residual = residual.to_coeffs();
    step("a[0] (uniform, ~2^40)", pk.a.to_coeffs()[0]);
    step(
        "b + a*s, first 8 coefficients",
        format!("{:?}", &residual[..8]),
    );
    passed &= confirm(
        &format!(
            "max |b + a*s| = {} is small (< {SMALL})",
            max_abs(&residual)
        ),
        max_abs(&residual) < SMALL,
    );

    println!("\n{}", heading("Step 3: relinearization key"));
    say(&format!(
        "Multiplying two ciphertexts leaves a term d2 * s^2. The\n\
         relinearization key encrypts s^2 under s, one entry per gadget\n\
         digit (base 2^{DIGIT_BITS}), so key switching can turn d2 * s^2 into\n\
         a pair under s (see key_switch):\n  {}",
        formula("b[k] + a[k] * s = s^2 * g_k + e_k")
    ));
    let rk =
        RelinearizationKey::<RnsPoly<N>, N>::generate(&sk, ERROR_STD, &mut rng);
    step("gadget terms (2 primes x 1 digit)", rk.ksk.a.len());
    let mut s_squared = sk.poly.clone();
    s_squared *= &sk.poly;
    let mut worst = 0;
    for (k, (a, b)) in rk.ksk.a.iter().zip(&rk.ksk.b).enumerate() {
        let mut residual = a.clone();
        residual *= &sk.poly;
        residual += b;
        residual -= &s_squared.gadget_term(k, DIGIT_BITS);
        let residual = residual.to_coeffs();
        if k == 0 {
            step(
                "b[0] + a[0]*s - s^2*g_0, first 8",
                format!("{:?}", &residual[..8]),
            );
        }
        worst = worst.max(max_abs(&residual));
    }
    passed &= confirm(
        &format!("max error over all terms = {worst} is small (< {SMALL})"),
        worst < SMALL,
    );

    println!("\n{}", heading("Takeaways"));
    say("s stays with the owner; (a, b) lets anyone encrypt; the\n\
         relinearization key (and rotation/conjugation keys, built the same\n\
         way) lets a server compute without decrypting. Next: encrypt_add.");
    println!(
        "{}",
        dim("Reference: fhetextbook.github.io/KeyGeneration.html")
    );
    if !passed {
        return Err("a key relation failed".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
