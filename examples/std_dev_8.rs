//! Encrypted variance of 8 numbers, one operation at a time.
//!
//! ```bash
//! cargo run --example std_dev_8
//! ```
//!
//! A narrated walk-through of the statistics pipeline at toy size (N=16,
//! 8 slots): sum the slots with a rotation tree, divide by 8, center,
//! square, sum again and divide. Each step is decrypted (a diagnostic only)
//! and compared with plain f64. The square root comes in
//! `sqrt_chebyshev_plain` and, encrypted at full size, `sensor_std_dev`.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::algorithms::sum_slots::sum_slots;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::evaluator::generate_rotation_key_set;
use toy_heaan_ckks::ckks::keys::SecretKey;
use toy_heaan_ckks::ckks::types::Ciphertext;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, formula, heading, ok, param, say, step, table};

const N: usize = 16;
const SLOTS: usize = 8;
const SCALE_BITS: i32 = 61;
const PRIMES: usize = 5;
const SEED: u64 = 42;
const BOUND: f64 = 1e-6;
const INPUT: [f64; SLOTS] = [0.12, 0.18, 0.21, 0.27, 0.31, 0.38, 0.42, 0.49];

type Ct = Ciphertext<RnsPoly<N>, N>;
type Sk = SecretKey<RnsPoly<N>, N>;
type Engine = CkksEngine<RnsPoly<N>, N>;

fn make_engine(basis: Arc<RnsBasis<N>>, scale: Scale) -> Engine {
    CkksEngine::new(
        basis,
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale,
        },
    )
}

/// Decrypts every slot (diagnostic only).
fn peek(
    ct: &Ct,
    sk: &Sk,
    encoder: &CkksEncoder<N>,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let sk = sk.mod_switch_to(ct.level())?;
    Ok(encoder.decode(&CkksEngine::decrypt(ct, &sk))[..SLOTS].to_vec())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Encrypted variance of 8 slots"));
    say(&format!(
        "Goal: from Enc([x0, ..., x7]) compute Enc([v, ..., v]) with\n  {}\n  {}\n\
         using only CKKS operations. Two tricks do most of the work:",
        formula("mean = (x0 + ... + x7) / 8"),
        formula("v    = ((x0 - mean)^2 + ... + (x7 - mean)^2) / 8")
    ));
    say(&format!(
        "  {}  x + rot(x,1), then + rot(.,2), then + rot(.,4):\n\
         \x20                after log2(8) = 3 steps every slot holds the sum\n  \
         {}  multiply by the constant 1/8 and rescale (one level)\n",
        formula("rotation tree"),
        formula("division     ")
    ));

    println!("{}", heading("Parameters"));
    let scale = Scale::power_of_two(SCALE_BITS)?;
    param("N     ", N, "8 slots, all used");
    param(
        "primes",
        "5 x 61 bits",
        "mean, square, variance: 3 levels, 2 to spare",
    );
    param("scale ", format!("2^{SCALE_BITS}"), "matches the primes");
    param(
        "seed  ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!();
    let basis = Arc::new(RnsBasis::<N>::new(generate_primes(
        SCALE_BITS as usize,
        PRIMES,
        N as u64,
    ))?);
    let engine = make_engine(basis.clone(), scale);
    let encoder = CkksEncoder::<N>::new(scale);
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let sk = engine.generate_secret_key(&mut rng)?;
    let pk = engine.generate_public_key(&sk, &mut rng)?;
    let rotation_keys = generate_rotation_key_set(&engine, &sk, SLOTS, &mut rng);

    let sum: f64 = INPUT.iter().sum();
    let mean = sum / SLOTS as f64;
    let centered: Vec<f64> = INPUT.iter().map(|x| x - mean).collect();
    let squares: Vec<f64> = centered.iter().map(|c| c * c).collect();
    let variance = squares.iter().sum::<f64>() / SLOTS as f64;
    step("x", format!("{INPUT:?}"));

    println!("\n{}", heading("Steps"));
    let input = engine.encrypt(&encoder.encode(&INPUT, basis), &pk, &mut rng);
    let mut rows = Vec::new();
    let mut worst: f64 = 0.0;
    let mut record = |name: &str,
                      ct: &Ct,
                      expected: &[f64]|
     -> Result<(), Box<dyn std::error::Error>> {
        let got = peek(ct, &sk, &encoder)?;
        let error = got
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        worst = worst.max(error);
        rows.push(vec![
            name.to_string(),
            format!("{}, 2^{:.1}", ct.level(), ct.scale.log2()),
            format!("{:.8}", got[0]),
            format!("{:.8}", expected[0]),
            format!("{error:.2e}"),
        ]);
        Ok(())
    };
    record("input x", &input, &INPUT)?;

    let total = sum_slots(&input, SLOTS, &rotation_keys);
    record("sum (rotation tree)", &total, &[sum; SLOTS])?;

    // Land on the input's exact scale so the subtraction is compatible.
    let mean_ct = total.mul_const_to_scale(1.0 / SLOTS as f64, input.scale)?;
    record("mean = sum * 1/8", &mean_ct, &[mean; SLOTS])?;

    let centered_ct = input.match_level(&mean_ct)?.sub(&mean_ct)?;
    record("x - mean", &centered_ct, &centered)?;

    // Keys for the lower level the centered values live at.
    let sk_low = sk.mod_switch_to(centered_ct.level())?;
    let engine_low = make_engine(centered_ct.c0.basis().clone(), scale);
    let relin = engine_low.generate_relin_key(&sk_low, &mut rng);
    let squared = centered_ct.mul(&centered_ct, &relin).rescale()?;
    record("(x - mean)^2", &squared, &squares)?;

    let sk_sq = sk.mod_switch_to(squared.level())?;
    let engine_sq = make_engine(squared.c0.basis().clone(), scale);
    let rotation_keys_sq =
        generate_rotation_key_set(&engine_sq, &sk_sq, SLOTS, &mut rng);
    let sum_squares = sum_slots(&squared, SLOTS, &rotation_keys_sq);
    let variance_ct =
        sum_squares.mul_const_to_scale(1.0 / SLOTS as f64, sum_squares.scale)?;
    record("variance", &variance_ct, &[variance; SLOTS])?;

    table(
        &["step", "level, scale", "slot 0", "expected", "max error"],
        &rows,
    );
    say(
        "\nRotations and additions keep the level; each division and the\n\
         square spend one. Keys are made for each level the data reaches.\n",
    );
    step(
        "std_dev = sqrt(variance)",
        format!("{:.8} (plain f64)", variance.sqrt()),
    );
    println!();
    let passed = check("every step, every slot", worst, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(
        "Statistics over packed slots cost log2(slots) rotations per sum\n\
         and one level per multiplication. sqrt is not a CKKS operation:\n\
         sqrt_chebyshev_plain fits it as a polynomial, sensor_std_dev runs\n\
         the full pipeline encrypted.",
    );
    if !passed {
        return Err("variance error above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
