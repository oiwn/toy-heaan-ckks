//! Bootstrap one stage at a time, checking each against its plain model.
//!
//! ```bash
//! cargo run --release --example bootstrap_stages
//! ```
//!
//! A narrated walk-through of the five public stages of
//! `BootstrapPlan::bootstrap`: prepare, ModRaise, CoeffsToSlots, EvalMod
//! and SlotsToCoeffs. After each stage the example decrypts (a diagnostic
//! only; bootstrapping itself never decrypts) and compares with the same
//! stage computed on plain numbers.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::f64::consts::TAU;
use std::sync::Arc;

use num_complex::Complex64;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::PolyRing;
use toy_heaan_ckks::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::SecretKey;
use toy_heaan_ckks::ckks::{Ciphertext, CkksResult};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, dim, formula, heading, ok, param, say, step};

const N: usize = 16;
const PRIMES: usize = 18;
const SEED: u64 = 4;
/// Every stage must match its plain model this closely, relative to the
/// largest value of that stage.
const STAGE_BOUND: f64 = 1e-3;
/// The refreshed slots must match the input this closely.
const BOUND: f64 = 1e-3;

type Ct = Ciphertext<RnsPoly<N>, N>;
type Sk = SecretKey<RnsPoly<N>, N>;

fn decrypt(ct: &Ct, secret: &Sk) -> CkksResult<Vec<Complex64>> {
    let pt = CkksEngine::decrypt(ct, &secret.mod_switch_to(ct.level())?);
    Ok(CkksEncoder::<N>::new(ct.scale).decode_complex(&pt))
}

fn coefficients(ct: &Ct, secret: &Sk) -> CkksResult<[f64; N]> {
    let pt = CkksEngine::decrypt(ct, &secret.mod_switch_to(ct.level())?);
    Ok(pt.poly.to_centered_f64())
}

fn max_error(actual: &[Complex64], expected: &[Complex64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max)
}

fn max_norm(values: &[Complex64]) -> f64 {
    values.iter().map(|v| v.norm()).fold(0.0, f64::max)
}

/// Level and scale of a stage output.
fn where_now(ct: &Ct) {
    step(
        "level, scale",
        format!("{}, 2^{:.1}", ct.level(), ct.scale.log2()),
    );
}

/// Compares a stage output with its plain model, relative to the model's
/// largest value.
fn compare(
    name: &str,
    ct: &Ct,
    expected: &[Complex64],
    secret: &Sk,
) -> CkksResult<bool> {
    let size = max_norm(expected);
    let error = max_error(&decrypt(ct, secret)?, expected);
    step(&format!("{name} max |value|"), format!("{size:.3e}"));
    Ok(check(
        &format!("{name} relative"),
        error / size,
        STAGE_BOUND,
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Bootstrapping, one stage at a time"));
    say(
        "Every rescale spends a prime. At level 1 none is left, and the\n\
         ciphertext cannot be multiplied again. Bootstrapping refreshes it:\n\
         it returns the same slots at a high level, homomorphically, using\n\
         only public evaluation keys.",
    );
    say(&format!(
        "\nThe idea: read the level-1 ciphertext modulo the full chain\n\
         (ModRaise). That is free, but every coefficient gains an unknown\n\
         multiple of the bottom prime: {}. The other stages\n\
         remove q0 * I by evaluating a function with period q0 on the\n\
         coefficients, which first have to be moved into slots.\n",
        formula("m -> m + q0 * I")
    ));

    println!("{}", heading("Parameters"));
    let basis =
        Arc::new(RnsBasis::<N>::new(generate_primes(61, PRIMES, N as u64))?);
    let config = BootstrapConfig::reference(N)?;
    let plan = BootstrapPlan::new(basis.clone(), config.clone())?;
    let scale = plan.schedule().message_scale();
    let q0 = plan.bottom_modulus() as f64;
    let radius = config.radius();
    param(
        "N            ",
        N,
        "small, so every slot fits on screen; 8 complex slots",
    );
    param(
        "primes       ",
        "18 x 61 bits",
        "bootstrap uses 11 levels, 7 remain",
    );
    param(
        "message scale",
        format!("2^{:.0}", scale.log2()),
        "q0 / message = 2^17 keeps m / q0 small",
    );
    param(
        "secret weight",
        config.max_secret_weight,
        "bounds I: |I| <= (weight + 2) / 2",
    );
    param(
        "radius       ",
        radius,
        "EvalMod input interval [-radius, radius]",
    );
    param(
        "seed         ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();

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
    let z: Vec<_> = (0..N / 2)
        .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
        .collect();
    let input = engine.encrypt(
        &CkksEncoder::<N>::new(scale).encode_complex(&z, basis),
        &public,
        &mut rng,
    );
    let mut passed = true;
    println!(
        "{}\n",
        dim("Decrypting after every stage is a diagnostic only.")
    );

    println!("{}", heading("Stage 1: prepare"));
    say(
        "Lift the scale to the message scale by an exact power of two and\n\
         drop to level 1. The slots do not change.",
    );
    let prepared = plan.prepare(&input)?;
    where_now(&prepared);
    passed &= compare("slots", &prepared, &z, &secret)?;

    println!("\n{}", heading("Stage 2: ModRaise"));
    say(&format!(
        "Read (c0, c1) modulo the full chain. Decryption now gives\n\
         {}, with I an integer polynomial.",
        formula("m + q0 * I")
    ));
    let raised = plan.mod_raise(&prepared)?;
    where_now(&raised);
    let m = coefficients(&prepared, &secret)?;
    let c = coefficients(&raised, &secret)?;
    let overflow: Vec<f64> = c.iter().zip(&m).map(|(c, m)| (c - m) / q0).collect();
    let largest = overflow.iter().fold(0.0f64, |a, i| a.max(i.abs()));
    let off_integer = overflow
        .iter()
        .fold(0.0f64, |a, i| a.max((i - i.round()).abs()));
    step(
        "max |I|",
        format!("{largest:.0} (bound {})", config.interval_bound()),
    );
    passed &= check("(c - m) / q0 is an integer", off_integer, 1e-9);

    println!("\n{}", heading("Stage 3: CoeffsToSlots"));
    say(&format!(
        "A linear transform (rotations and plaintext diagonals) moves\n\
         coefficient j into the real part and j + N/2 into the imaginary\n\
         part of a slot, in bit-reversed order. A conjugation splits them\n\
         into two ciphertexts with real slots:\n  {}",
        formula("slot j = c[rev(j)] / (radius * q0) = (I + m/q0) / radius")
    ));
    let (real, imaginary) = plan.coeffs_to_slots(&raised, &keys)?;
    where_now(&real);
    let index = |j: usize| j.reverse_bits() >> (usize::BITS - (N / 2).ilog2());
    let model = |values: &[f64; N], offset: usize, f: &dyn Fn(f64) -> f64| {
        (0..N / 2)
            .map(|j| Complex64::new(f(values[index(j) + offset]), 0.0))
            .collect::<Vec<_>>()
    };
    let to_slot = |coefficient: f64| coefficient / (radius * q0);
    passed &= compare("real half", &real, &model(&c, 0, &to_slot), &secret)?;
    passed &= compare(
        "imaginary half",
        &imaginary,
        &model(&c, N / 2, &to_slot),
        &secret,
    )?;

    println!("\n{}", heading("Stage 4: EvalMod"));
    say(&format!(
        "A Chebyshev polynomial of a scaled cosine, then double angles,\n\
         evaluates {} on each half. The sine has\n\
         period 1 in I + m/q0, so the integer I drops out, and since\n\
         m << q0 the result is about 2 * pi * m / q0. This stage uses most\n\
         of the levels.",
        formula("sin(2 * pi * (I + m/q0))")
    ));
    let real = plan.eval_mod(&real, &keys)?;
    let imaginary = plan.eval_mod(&imaginary, &keys)?;
    where_now(&real);
    let reduce = |coefficient: f64| (TAU * coefficient / q0).sin();
    passed &= compare("real half", &real, &model(&m, 0, &reduce), &secret)?;
    passed &= compare(
        "imaginary half",
        &imaginary,
        &model(&m, N / 2, &reduce),
        &secret,
    )?;

    println!("\n{}", heading("Stage 5: SlotsToCoeffs"));
    say(
        "Recombine real + i * imaginary and apply the inverse transform,\n\
         with the factor q0 / (2 * pi * message scale) folded into its\n\
         diagonals. The coefficients are m again: the original slots, now\n\
         at a high level.",
    );
    let output = plan.slots_to_coeffs(&real, &imaginary, prepared.scale, &keys)?;
    where_now(&output);
    let error = max_error(&decrypt(&output, &secret)?, &z);
    passed &= check("refreshed slots against the input", error, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(&format!(
        "Level 1 -> {}: {} levels spent on the refresh, {} free for\n\
         new multiplications. The approximation (sin x ~ x, the Chebyshev\n\
         fit) bounds the accuracy; bootstrap_demo runs whole refresh\n\
         cycles, recurrence_bootstrap_demo runs 32 rounds with refreshes.",
        output.level(),
        plan.schedule().depth(),
        output.level() - 1
    ));
    println!(
        "{}",
        dim("References: CKKS bootstrapping (eprint 2018/153);\n\
             OpenFHE ckksrns-fhe; Lattigo ckks/bootstrapping.")
    );
    if !passed {
        return Err("a bootstrap stage missed its bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
