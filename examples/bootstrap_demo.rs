//! Toy complex bootstrap followed by multiplication and a second refresh.
use std::sync::Arc;

use num_complex::Complex64;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::RelinearizationKey;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

fn run<const N: usize>() -> Result<(), Box<dyn std::error::Error>> {
    let primes = generate_primes(61, 18, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes)?);
    let config = BootstrapConfig::reference(N)?;
    let plan = BootstrapPlan::new(basis.clone(), config.clone())?;
    // Encrypt at the bootstrap's message scale (2^44 at N=16/64).
    let scale = plan.schedule().message_scale();
    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: config.max_secret_weight,
            scale,
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let secret = engine.generate_secret_key(&mut rng)?;
    let public = engine.generate_public_key(&secret, &mut rng)?;
    println!(
        "N={N} bootstrap depth={} output_level={} message_scale=2^{:.0}",
        plan.schedule().depth(),
        plan.schedule().output_level(),
        scale.log2()
    );
    let keys = plan.generate_keys(&secret, 3.2f64.sqrt(), &mut rng)?;
    let encoder = CkksEncoder::new(scale);
    let mut expected: Vec<_> = (0..N / 2)
        .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
        .collect();
    let initial = engine.encrypt(
        &encoder.encode_complex(&expected, basis),
        &public,
        &mut rng,
    );
    let relin = RelinearizationKey::generate(&secret, 3.2f64.sqrt(), &mut rng);
    let mut squared = initial.mul(&initial, &relin).rescale()?;
    for x in &mut expected {
        *x *= *x;
    }
    for cycle in 0..2 {
        let refreshed = plan.bootstrap(&squared, &keys)?;
        let sk = secret.mod_switch_to(refreshed.level())?;
        let actual = encoder.decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
        let error = actual
            .iter()
            .zip(&expected)
            .map(|(a, b)| (*a - *b).norm())
            .fold(0.0, f64::max);
        println!(
            "N={N} cycle={cycle} level={} log2(scale)={:.1} max_error={error:e}",
            refreshed.level(),
            refreshed.scale.log2()
        );
        if error > 1e-3 || !error.is_finite() {
            return Err("bootstrap accuracy target failed".into());
        }
        let relin = RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
        squared = refreshed.mul(&refreshed, &relin).rescale()?;
        for x in &mut expected {
            *x *= *x;
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    run::<16>()?;
    run::<64>()
}
