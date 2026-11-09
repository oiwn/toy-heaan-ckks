//! A chain of five multiply-add steps over 4096 slots.
//!
//! ```bash
//! cargo run --release --example horner_chain   # ~30 s in debug
//! ```
//!
//! A narrated walk-through of multiplicative depth: apply `x -> alpha*x +
//! beta` five times at N=8192, one prime spent per multiplication, and watch
//! level, scale and modulus shrink. Also shows why the last level needs
//! headroom above the scale.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::{PublicKey, RelinearizationKey};
use toy_heaan_ckks::ckks::types::Ciphertext;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, dim, formula, heading, ok, param, say, step, table};

const N: usize = 8192;
const SLOTS: usize = N / 2;
const SCALE_BITS: i32 = 61;
const ITERS: usize = 5;
/// One prime per multiplication, plus two so the last level holds the
/// message with room to spare.
const PRIMES: usize = ITERS + 2;
const ALPHA: f64 = 0.8;
const BETA: f64 = 0.1;
const SEED: u64 = 42;
const BOUND: f64 = 1e-5;

type Ct = Ciphertext<RnsPoly<N>, N>;
type Basis = Arc<RnsBasis<N>>;
type Engine = CkksEngine<RnsPoly<N>, N>;

fn make_engine(basis: Basis, scale: Scale) -> Engine {
    CkksEngine::new(
        basis,
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale,
        },
    )
}

/// Encrypt `value` in every slot, at `encoder`'s scale on `basis`.
fn encrypt_constant(
    value: f64,
    basis: Basis,
    encoder: &CkksEncoder<N>,
    engine: &Engine,
    pk: &PublicKey<RnsPoly<N>, N>,
    rng: &mut ChaCha20Rng,
) -> Ct {
    engine.encrypt(&encoder.encode(&vec![value; SLOTS], basis), pk, rng)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Depth: five multiply-add steps"));
    say(&format!(
        "Each step computes {}. Unrolled, five steps give\n  {}\n\
         Every multiplication squares the scale and needs a rescale, which\n\
         spends one prime. The number of primes is the depth budget.",
        formula("x -> alpha * x + beta"),
        formula("x5 = alpha^5 * x0 + beta * (alpha^4 + ... + alpha + 1)")
    ));
    say(
        "\nalpha and beta are encrypted here (as a client would send them),\n\
         each freshly at the accumulator's current level. Keys for a lower\n\
         level come from the secret key cut down with mod_switch_to.\n",
    );

    println!("{}", heading("Parameters"));
    let scale = Scale::power_of_two(SCALE_BITS)?;
    param(
        "N          ",
        N,
        "4096 slots: one operation acts on all of them",
    );
    param(
        "primes     ",
        "7 x 61 bits",
        "5 multiplications + 2 left at the end",
    );
    param(
        "scale      ",
        format!("2^{SCALE_BITS}"),
        "equals the prime size, so rescale returns ~2^61",
    );
    param(
        "alpha, beta",
        format!("{ALPHA}, {BETA}"),
        "contracts to beta/(1-alpha) = 0.5",
    );
    param(
        "seed       ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!(
        "  {}\n",
        dim(
            "Why 2 primes at the end: the mean slot (~0.5) puts the constant\n  \
             coefficient near scale/2 = 2^60, above q0/2 with one 61-bit prime."
        )
    );

    let primes = generate_primes(SCALE_BITS as usize, PRIMES, N as u64);
    let basis_top = Arc::new(RnsBasis::<N>::new(primes)?);
    let engine_top = make_engine(basis_top.clone(), scale);
    let encoder = CkksEncoder::<N>::new(scale);
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let sk_top = engine_top.generate_secret_key(&mut rng)?;
    let pk_top = engine_top.generate_public_key(&sk_top, &mut rng)?;
    let mut rlk: RelinearizationKey<RnsPoly<N>, N> =
        engine_top.generate_relin_key(&sk_top, &mut rng);

    let x0: Vec<f64> = (0..SLOTS).map(|i| (i + 1) as f64 / SLOTS as f64).collect();
    let mut ct =
        engine_top.encrypt(&encoder.encode(&x0, basis_top), &pk_top, &mut rng);
    let mut expected = x0.clone();
    let mut sk = sk_top.clone();
    let mut pk = pk_top;
    let mut engine = engine_top;

    println!("{}", heading("Steps"));
    let mut rows = vec![vec![
        "start".to_string(),
        ct.level().to_string(),
        ct.logq().to_string(),
        format!("2^{:.1}", ct.scale.log2()),
    ]];
    for iter in 1..=ITERS {
        let alpha = encrypt_constant(
            ALPHA,
            ct.c0.basis().clone(),
            &encoder,
            &engine,
            &pk,
            &mut rng,
        );
        let product = ct.mul(&alpha, &rlk);
        rows.push(vec![
            format!("{iter}: mul alpha"),
            product.level().to_string(),
            product.logq().to_string(),
            format!("2^{:.1}", product.scale.log2()),
        ]);
        ct = product.rescale()?;
        rows.push(vec![
            format!("{iter}: rescale"),
            ct.level().to_string(),
            ct.logq().to_string(),
            format!("2^{:.1}", ct.scale.log2()),
        ]);

        // Keys for the new, lower level.
        let basis = ct.c0.basis().clone();
        sk = sk_top.mod_switch_to(ct.level())?;
        engine = make_engine(basis.clone(), scale);
        pk = engine.generate_public_key(&sk, &mut rng)?;
        // beta at ct's exact scale (rescale divided by the actual prime).
        let beta = encrypt_constant(
            BETA,
            basis,
            &CkksEncoder::<N>::new(ct.scale),
            &engine,
            &pk,
            &mut rng,
        );
        ct = ct.add(&beta)?;
        for v in &mut expected {
            *v = *v * ALPHA + BETA;
        }
        if iter < ITERS {
            rlk = engine.generate_relin_key(&sk, &mut rng);
        }
    }
    table(&["step", "level", "logq", "scale"], &rows);
    say(
        "\nmul doubles the scale bits, rescale halves them and removes ~61\n\
         bits of modulus; add changes nothing.\n",
    );

    println!("{}", heading("Result"));
    let decoded = encoder.decode(&CkksEngine::<RnsPoly<N>, N>::decrypt(&ct, &sk));
    let mut rows = Vec::new();
    for j in [0, 1, 2, SLOTS / 2, SLOTS - 1] {
        rows.push(vec![
            j.to_string(),
            format!("{:.6}", x0[j]),
            format!("{:.6}", expected[j]),
            format!("{:.6}", decoded[j]),
        ]);
    }
    table(&["slot", "x0", "expected x5", "decrypted x5"], &rows);
    println!(
        "  {}",
        dim(&format!("(5 of {SLOTS} slots shown; all are checked)"))
    );
    let error = expected
        .iter()
        .zip(&decoded)
        .map(|(e, g)| (e - g).abs())
        .fold(0.0, f64::max);
    step(
        "fixed point beta / (1 - alpha)",
        format!("{:.6}", BETA / (1.0 - ALPHA)),
    );
    println!();
    let passed = check(&format!("all {SLOTS} slots"), error, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(
        "Five multiplications cost five primes; the circuit's depth, not\n\
         its size, decides the chain length. When the primes run out,\n\
         bootstrapping refreshes the ciphertext (bootstrap_demo).",
    );
    if !passed {
        return Err("Horner chain error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
