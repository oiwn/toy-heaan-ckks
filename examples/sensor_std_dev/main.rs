//! Encrypted standard deviation of a day of temperature readings.
//!
//! ```bash
//! cargo run --release --example sensor_std_dev
//! ```
//!
//! A narrated application: a client normalizes and encrypts 24 h of sensor
//! readings, a server computes their standard deviation on the ciphertext
//! (slot sums, a square, and a Chebyshev square root), and the client
//! decrypts and compares with plain f64 math. N = 8192 (32 of 4096 slots
//! used), 13 x 61-bit primes; educational, reduced-security parameters
//! (see `specs/demo.md`).
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.

mod contract;
#[path = "../support/narrate.rs"]
mod narrate;
mod std_dev;

use std::sync::Arc;
use std::time::Instant;

use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::CkksParams;
use toy_heaan_ckks::ckks::evaluator::LeveledEvaluator;
use toy_heaan_ckks::ckks::{CkksEngine, Scale};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

use contract::SensorContract;
use narrate::{
    check, confirm, dim, elapsed, formula, heading, ok, param, say, step, table,
};
use std_dev::encrypted_std_dev;

const N: usize = 8192;
const READINGS: usize = 32;
const SCALE_BITS: i32 = 61;
const CHEBYSHEV_DEGREE: usize = 7;
const DATA_SEED: u64 = 2026;
const CRYPTO_SEED: u64 = 4321;

/// Std error budget in degC (degree-7 fit on the [1, 8] degC envelope).
const STD_BUDGET: f64 = 0.1;
/// Intermediate variance error budget in degC^2: CKKS key-switch noise at
/// N = 8192 (measured ~5e-4 degC^2), far below the sqrt approximation error.
const VARIANCE_BUDGET: f64 = 1e-2;

/// One day of readings: a diurnal curve (coldest at dawn, warmest mid
/// afternoon) of 15 +- 4 degC plus up to 0.5 degC of seeded jitter.
fn diurnal_readings(seed: u64) -> [f64; READINGS] {
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    std::array::from_fn(|i| {
        let hour = i as f64 * 24.0 / READINGS as f64;
        let phase = (hour - 9.0) / 24.0 * std::f64::consts::TAU;
        15.0 + 4.0 * phase.sin() + rng.random_range(-0.5..0.5)
    })
}

/// Population variance (divide by n), matching the encrypted pipeline.
fn population_variance(values: &[f64]) -> f64 {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / values.len() as f64
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "{}\n",
        heading("Encrypted standard deviation of sensor data")
    );
    say(&format!(
        "A sensor owner wants the day's temperature spread computed by a\n\
         server that must not see the readings. The server evaluates\n  {}\n  {}\n\
         on encrypted slots. Every step is a CKKS operation: rotations sum\n\
         the slots, one multiplication squares, and the square root, which\n\
         CKKS cannot compute exactly, is a degree-{CHEBYSHEV_DEGREE} Chebyshev polynomial.",
        formula("variance = mean((x - mean(x))^2)"),
        formula("std_dev  = sqrt(variance)")
    ));
    say(
        "\nA client-side contract fixes the reading range and the expected\n\
         std_dev range. Readings are normalized to [0, 1], and the\n\
         polynomial is fitted only on that std_dev range, which keeps its\n\
         degree low.\n",
    );

    let contract = SensorContract::demo();
    let readings = diurnal_readings(DATA_SEED);
    let variance = population_variance(&readings);
    let std_dev = variance.sqrt();

    println!("{}", heading("Parameters"));
    param("readings        ", READINGS, "one per 45 min over 24 h");
    param(
        "reading range   ",
        format!("[{}, {}] degC", contract.lo, contract.hi),
        "the contract clamps to it",
    );
    param(
        "std_dev range   ",
        format!("[{}, {}] degC", contract.std_min, contract.std_max),
        "the sqrt is fitted on it",
    );
    param(
        "N               ",
        N,
        "4096 slots; the 32 readings are tiled",
    );
    param(
        "primes          ",
        format!("{} x 61 bits", CHEBYSHEV_DEGREE + 6),
        "6 levels for statistics + sqrt depth",
    );
    param(
        "scale           ",
        format!("2^{SCALE_BITS}"),
        "matches the primes",
    );
    param(
        "Chebyshev degree",
        CHEBYSHEV_DEGREE,
        "error budget 0.1 degC",
    );
    println!();

    println!("{}", heading("Client: check, normalize, encrypt"));
    step(
        "readings",
        format!(
            "{:.1} .. {:.1} degC",
            readings.iter().copied().fold(f64::INFINITY, f64::min),
            readings.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        ),
    );
    let mut passed = confirm(
        &format!("std_dev {std_dev:.3} degC lies in the contract range"),
        contract.contains_std(std_dev),
    );
    // Tile the readings over all N/2 slots: rotations wrap over every slot,
    // so each 32-slot window must hold one full copy of the day.
    let normalized: Vec<f64> = (0..N / 2)
        .map(|slot| contract.normalize(readings[slot % READINGS]))
        .collect();
    let primes =
        generate_primes(SCALE_BITS as usize, CHEBYSHEV_DEGREE + 6, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes)?);
    let scale = Scale::power_of_two(SCALE_BITS)?;
    let params = CkksParams {
        error_variance: 3.2,
        hamming_weight: N / 2,
        scale,
    };
    let engine = CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
    let mut rng = ChaCha20Rng::seed_from_u64(CRYPTO_SEED);
    let started = Instant::now();
    let sk = engine.generate_secret_key(&mut rng)?;
    let pk = engine.generate_public_key(&sk, &mut rng)?;
    let encoder = CkksEncoder::<N>::new(scale);
    let ct =
        engine.encrypt(&encoder.encode(&normalized, basis.clone()), &pk, &mut rng);
    println!(
        "  {}",
        dim(&format!(
            "keys + encode + encrypt {}",
            elapsed(started.elapsed())
        ))
    );

    println!("\n{}", heading("Server: compute on ciphertexts"));
    say("Each stage spends levels (channels) and keeps the scale near 2^61.\n");
    let mut evaluator = LeveledEvaluator::new(params, sk);
    let started = Instant::now();
    let run = encrypted_std_dev(
        &ct,
        READINGS,
        &contract,
        CHEBYSHEV_DEGREE,
        &mut evaluator,
        &mut rng,
    );
    let compute_time = started.elapsed();
    let mut rows = vec![vec![
        "fresh ciphertext".to_string(),
        basis.channel_count().to_string(),
        format!("2^{:.1}", ct.scale.log2()),
        "-".to_string(),
    ]];
    for stage in &run.stages {
        rows.push(vec![
            stage.name.to_string(),
            stage.channels.to_string(),
            format!("2^{:.1}", stage.log2_scale),
            elapsed(stage.elapsed),
        ]);
    }
    table(&["stage", "channels", "scale", "time"], &rows);
    println!(
        "  {}",
        dim(&format!(
            "total, including per-level keys: {}",
            elapsed(compute_time)
        ))
    );

    println!("\n{}", heading("Client: decrypt and compare"));
    let variance_got = encoder.decode(&evaluator.decrypt(&run.variance))[0]
        * contract.width().powi(2);
    let std_got = contract
        .restore_std(encoder.decode(&evaluator.decrypt(&run.std_normalized))[0]);
    table(
        &["", "plain f64", "encrypted", "error"],
        &[
            vec![
                "variance (degC^2)".to_string(),
                format!("{variance:.9}"),
                format!("{variance_got:.9}"),
                format!("{:.2e}", (variance_got - variance).abs()),
            ],
            vec![
                "std_dev (degC)".to_string(),
                format!("{std_dev:.9}"),
                format!("{std_got:.9}"),
                format!("{:.2e}", (std_got - std_dev).abs()),
            ],
        ],
    );
    say(
        "\nThe variance error is CKKS noise; the std_dev error is mostly the\n\
         degree-7 square-root approximation.\n",
    );
    passed &= check(
        "variance (degC^2)",
        (variance_got - variance).abs(),
        VARIANCE_BUDGET,
    );
    passed &= check("std_dev (degC)", (std_got - std_dev).abs(), STD_BUDGET);

    println!("\n{}", heading("Takeaways"));
    say(
        "Statistics need only additions, rotations and one multiplication;\n\
         nonlinear steps like sqrt become polynomials, accurate only on a\n\
         range agreed in advance. sensor_std_plain shows the same pipeline\n\
         in plain f64, step by step.",
    );
    if !passed {
        return Err("sensor std_dev missed its budget".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
