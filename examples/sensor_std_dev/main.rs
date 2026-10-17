//! Sensor std-dev demo: encrypted standard deviation of temperature readings.
//!
//! A sensor collected 24 h of readings (one per 45 min). The client clamps
//! and normalizes them with a [`SensorContract`] and encrypts them; the
//! "server" computes the std on the ciphertext with library ops only; the
//! client decrypts and restores degC, then compares with plain f64 math.
//!
//! Demo parameters: N = 8192 (32 of 4096 slots used), 13 x 61-bit primes,
//! degree-7 Chebyshev sqrt. These are educational, reduced-security
//! parameters (see `specs/demo.md`).
//!
//! Run with: `cargo run --release --example sensor_std_dev`

mod contract;
mod std_dev;

use std::sync::Arc;
use std::time::Instant;

use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::CkksEngine;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::CkksParams;
use toy_heaan_ckks::ckks::evaluator::LeveledEvaluator;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

use contract::SensorContract;
use std_dev::encrypted_std_dev;

const N: usize = 8192;
const READINGS: usize = 32;
const SCALE: u32 = 61;
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

fn main() {
    let contract = SensorContract::demo();
    let readings = diurnal_readings(DATA_SEED);

    // ── Client: plaintext reference + envelope check ─────────────────────────
    let variance = population_variance(&readings);
    let std_dev = variance.sqrt();
    assert!(
        contract.contains_std(std_dev),
        "std {std_dev} is outside the contract envelope"
    );
    // Tile the readings over all N/2 slots: rotations wrap over every slot,
    // so each 32-slot window must hold one full copy of the day.
    let normalized: Vec<f64> = (0..N / 2)
        .map(|slot| contract.normalize(readings[slot % READINGS]))
        .collect();

    println!("== Sensor std-dev demo ==");
    println!(
        "Readings:   {READINGS} x 45 min, {:.1}..{:.1} degC",
        readings.iter().copied().fold(f64::INFINITY, f64::min),
        readings.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    );
    println!(
        "Contract:   range [{}, {}] degC, std envelope [{}, {}] degC",
        contract.lo, contract.hi, contract.std_min, contract.std_max
    );

    // ── Parameters and keys ──────────────────────────────────────────────────
    let primes = generate_primes(SCALE as usize, CHEBYSHEV_DEGREE + 6, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();
    println!(
        "Params:     N = {N}, {} x {SCALE}-bit primes (logq = {logq}), \
         scale 2^{SCALE}, Chebyshev degree {CHEBYSHEV_DEGREE}",
        basis.channel_count()
    );

    let params = CkksParams {
        error_variance: 3.2,
        hamming_weight: N / 2,
        scale_bits: SCALE,
    };
    let engine = CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
    let mut rng = ChaCha20Rng::seed_from_u64(CRYPTO_SEED);

    let started = Instant::now();
    let sk = engine.generate_secret_key(&mut rng).expect("secret key");
    let pk = engine
        .generate_public_key(&sk, &mut rng)
        .expect("public key");
    println!("\nKeygen (sk, pk):        {:>9.2?}", started.elapsed());

    // ── Client: encode + encrypt ─────────────────────────────────────────────
    let started = Instant::now();
    let encoder = CkksEncoder::<N>::new(SCALE);
    let pt = encoder.encode(&normalized, basis.clone());
    let ct = engine.encrypt(&pt, &pk, &mut rng);
    println!("Encode + encrypt:       {:>9.2?}", started.elapsed());

    // ── Server: encrypted std-dev (gadget keys generated per level) ──────────
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

    println!("\nStage               channels  logp       time");
    println!(
        "fresh ciphertext    {:>8}  {:>4}          -",
        basis.channel_count(),
        ct.logp
    );
    for stage in &run.stages {
        println!(
            "{:<18}  {:>8}  {:>4}  {:>9.2?}",
            stage.name, stage.channels, stage.logp, stage.elapsed
        );
    }
    println!("total (incl. keys)                {compute_time:>9.2?}");

    // ── Client: decrypt + restore ────────────────────────────────────────────
    let started = Instant::now();
    let variance_got = encoder.decode(&evaluator.decrypt(&run.variance))[0]
        * contract.width().powi(2);
    let std_got = contract
        .restore_std(encoder.decode(&evaluator.decrypt(&run.std_normalized))[0]);
    println!("\nDecrypt + decode (x2):  {:>9.2?}", started.elapsed());

    let variance_error = (variance_got - variance).abs();
    let std_error = (std_got - std_dev).abs();

    println!("\n                     plaintext f64       encrypted      abs error");
    println!(
        "variance (degC^2) {variance:>16.9} {variance_got:>15.9} {variance_error:>14.2e}"
    );
    println!(
        "std-dev (degC)    {std_dev:>16.9} {std_got:>15.9} {std_error:>14.2e}"
    );
    println!(
        "\nVariance error is CKKS noise; std error is the degree-{CHEBYSHEV_DEGREE} \
         sqrt approximation."
    );

    assert!(
        variance_error < VARIANCE_BUDGET,
        "variance error {variance_error:e} exceeds {VARIANCE_BUDGET:e} degC^2"
    );
    assert!(
        std_error < STD_BUDGET,
        "std error {std_error} exceeds the {STD_BUDGET} degC budget"
    );
}
