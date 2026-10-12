//! Sensor std-dev demo: encrypted standard deviation of temperature readings.
//!
//! The client clamps and normalizes readings with a [`SensorContract`],
//! encrypts them, and the "server" computes the std on the ciphertext with
//! library ops only. The client decrypts and restores degC.
//!
//! Toy parameters (N = 16, 8 readings, degree 7, 13 x 61-bit primes) for now;
//! the N = 8192 demo parameters come next.
//!
//! Run with: `cargo run --release --example sensor_std_dev`

mod contract;
mod std_dev;

use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::CkksEngine;
use toy_heaan_ckks::crypto::engine::CkksParams;
use toy_heaan_ckks::crypto::pipeline::{LeveledEvaluator, decode_wide};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

use contract::SensorContract;
use std_dev::encrypted_std_dev;

const N: usize = 16;
const SLOTS: usize = 8;
const SCALE: u32 = 61;
const DEGREE: usize = 7;

fn main() {
    let contract = SensorContract::demo();
    let readings: [f64; SLOTS] = [12.0, 14.0, 17.0, 20.0, 21.0, 19.0, 16.0, 13.0];

    // ── Client: plaintext reference + envelope check ─────────────────────────
    let mean = readings.iter().sum::<f64>() / SLOTS as f64;
    let std_dev = (readings.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
        / SLOTS as f64)
        .sqrt();
    assert!(
        contract.contains_std(std_dev),
        "std {std_dev} is outside the contract envelope"
    );
    let normalized: Vec<f64> =
        readings.iter().map(|&x| contract.normalize(x)).collect();

    println!("Readings (degC):  {readings:?}");
    println!(
        "Contract:         range [{}, {}] degC, std envelope [{}, {}] degC",
        contract.lo, contract.hi, contract.std_min, contract.std_max
    );

    // ── Keys and encryption ──────────────────────────────────────────────────
    let primes = generate_primes(SCALE as usize, DEGREE + 6, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();
    let params = CkksParams {
        error_variance: 3.2,
        hamming_weight: N / 2,
        scale_bits: SCALE,
    };
    let engine = CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
    let mut rng = ChaCha20Rng::seed_from_u64(4321);
    let sk = engine.generate_secret_key(&mut rng).expect("secret key");
    let pk = engine
        .generate_public_key(&sk, &mut rng)
        .expect("public key");

    let pt = CkksEncoder::<N>::new(SCALE).encode(&normalized, basis.clone());
    let ct = engine.encrypt(&pt, &pk, logq, &mut rng);
    println!(
        "Encrypted:        {} RNS channels, logq = {logq}",
        basis.channel_count()
    );

    // ── Server: encrypted std-dev ────────────────────────────────────────────
    let mut evaluator = LeveledEvaluator::new(params, sk);
    let result =
        encrypted_std_dev(&ct, SLOTS, &contract, DEGREE, &mut evaluator, &mut rng);
    println!(
        "Computed:         {} RNS channels left",
        result.c0.basis().channel_count()
    );

    // ── Client: decrypt and restore ──────────────────────────────────────────
    let decoded = decode_wide(&evaluator.decrypt(&result));
    let encrypted_std = contract.restore_std(decoded[0]);
    let error = (encrypted_std - std_dev).abs();

    println!("Plaintext std:    {std_dev:.6} degC");
    println!("Encrypted std:    {encrypted_std:.6} degC");
    println!("Absolute error:   {error:.2e} degC");
    assert!(error < 0.1, "std error {error} exceeds the 0.1 degC budget");
}
