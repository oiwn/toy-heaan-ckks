use std::sync::Arc;

use num_complex::Complex64;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::PolyRing;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::{
    ConjugationKey, RelinearizationKey, RotationKey, SecretKey,
};
use toy_heaan_ckks::ckks::types::Ciphertext;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

fn pow2(exponent: i32) -> Scale {
    Scale::power_of_two(exponent).unwrap()
}

fn assert_same<const N: usize>(
    left: &Ciphertext<RnsPoly<N>, N>,
    right: &Ciphertext<RnsPoly<N>, N>,
) {
    assert_eq!((left.level(), left.scale), (right.level(), right.scale));
    for (a, b) in [(&left.c0, &right.c0), (&left.c1, &right.c1)] {
        let mut difference = a.clone();
        difference -= b;
        assert_eq!(difference.to_coeffs(), [0; N]);
    }
}

fn suite<const N: usize>() {
    let primes = generate_primes(61, 18, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes.clone()).unwrap());
    let plan = BootstrapPlan::<RnsPoly<N>, N>::new(
        basis.clone(),
        BootstrapConfig::reference(N).unwrap(),
    )
    .unwrap();
    let encoder = CkksEncoder::new(pow2(44));
    let engine = CkksEngine::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale: pow2(44),
        },
    );
    let other = Arc::new(RnsBasis::<N>::new(primes.clone()).unwrap());
    for seed in 0..8 {
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let secret = engine.generate_secret_key(&mut rng).unwrap();
        let public = engine.generate_public_key(&secret, &mut rng).unwrap();
        let keys = plan
            .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
            .unwrap();

        for scenario in [
            "zero",
            "constant",
            "asymmetric",
            "unit_circle",
            "random_disk",
            "alternating",
        ]
        .iter()
        {
            let values: Vec<_> = (0..N / 2)
                .map(|j| match *scenario {
                    "zero" => Complex64::default(),
                    "constant" => Complex64::new(0.5, -0.25),
                    "asymmetric" => Complex64::from_polar(0.8, j as f64 * 0.7),
                    "unit_circle" => Complex64::from_polar(1.0, j as f64 * 0.9),
                    "random_disk" => Complex64::from_polar(
                        rng.random::<f64>().sqrt(),
                        rng.random::<f64>() * std::f64::consts::TAU,
                    ),
                    _ => Complex64::new(if j % 2 == 0 { 1.0 } else { -1.0 }, 0.0),
                })
                .collect();
            let initial = engine.encrypt(
                &encoder.encode_complex(&values, basis.clone()),
                &public,
                &mut rng,
            );
            if *scenario == "asymmetric" {
                check_automorphisms(&initial, &values, &secret, &mut rng);
            }
            let relin =
                RelinearizationKey::generate(&secret, 3.2f64.sqrt(), &mut rng);
            let mut squared = initial.mul(&initial, &relin).rescale().unwrap();
            let mut expected: Vec<_> = values.iter().map(|x| x * x).collect();
            let not_power_of_two = Scale::new(1.5 * 2f64.powi(44)).unwrap();
            for scale in [pow2(0), pow2(45), not_power_of_two] {
                let mut wrong_scale = squared.clone();
                wrong_scale.scale = scale;
                assert!(plan.bootstrap(&wrong_scale, &keys).is_err());
            }
            if seed == 0 && (*scenario == "zero" || *scenario == "asymmetric") {
                let snapshot = initial.clone();
                let fresh = plan.bootstrap(&initial, &keys).unwrap();
                let decoded = encoder.decode_complex(&CkksEngine::decrypt(
                    &fresh,
                    &secret.mod_switch_to(fresh.level()).unwrap(),
                ));
                for (a, b) in decoded.iter().zip(&values) {
                    assert!((*a - *b).norm() < 1e-3);
                }
                assert_same(&snapshot, &initial);

                let mut wrong_context = squared.clone();
                wrong_context.c0 =
                    RnsPoly::zero(&other.at_level(squared.level()).unwrap());
                wrong_context.c1 = wrong_context.c0.clone();
                assert!(plan.bootstrap(&wrong_context, &keys).is_err());
                wrong_context.c1 = squared.c1.clone();
                assert!(plan.bootstrap(&wrong_context, &keys).is_err());
                let other_plan = BootstrapPlan::<RnsPoly<N>, N>::new(
                    other.clone(),
                    BootstrapConfig::reference(N).unwrap(),
                )
                .unwrap();
                assert!(other_plan.bootstrap(&squared, &keys).is_err());
            }

            for cycle in 0..2 {
                let snapshot = squared.clone();
                let explicit =
                    squared.mul_const(1.0, pow2(17)).mod_switch_to(1).unwrap();

                let refreshed = plan.bootstrap(&squared, &keys).unwrap();
                assert_same(&squared, &snapshot);
                if seed == 0 && *scenario == "zero" && cycle == 0 {
                    assert_same(
                        &refreshed,
                        &plan.bootstrap(&explicit, &keys).unwrap(),
                    );
                }
                // Output properties come from the schedule; the scale is
                // tracked, not snapped, and the input is lifted by 2^17.
                let schedule = plan.schedule();
                assert_eq!(refreshed.level(), schedule.output_level());
                assert!(
                    refreshed.scale.is_compatible(
                        schedule.output_scale(explicit.scale).unwrap()
                    )
                );
                assert!((refreshed.scale.log2() - 44.0).abs() < 1e-9);
                let sk = secret.mod_switch_to(refreshed.level()).unwrap();
                let actual =
                    encoder.decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
                let before = encoder.decode_complex(&CkksEngine::decrypt(
                    &squared,
                    &secret.mod_switch_to(squared.level()).unwrap(),
                ));
                for ((a, want), prior) in actual.iter().zip(&expected).zip(before) {
                    assert!(
                        (*a - *want).norm() < 1e-3,
                        "N={N} seed={seed} case={scenario} cycle={cycle} got={a} expected={want}"
                    );
                    assert!((*a - prior).norm() < 1e-3);
                }
                let key =
                    RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
                squared = refreshed.mul(&refreshed, &key).rescale().unwrap();
                for x in &mut expected {
                    *x *= *x;
                }
            }
        }
    }
}

#[test]
fn complex_bootstrap_16() {
    suite::<16>();
}
#[test]
fn complex_bootstrap_64() {
    suite::<64>();
}

fn check_automorphisms<const N: usize>(
    initial: &Ciphertext<RnsPoly<N>, N>,
    values: &[Complex64],
    secret: &SecretKey<RnsPoly<N>, N>,
    rng: &mut ChaCha20Rng,
) {
    let encoder = CkksEncoder::new(pow2(44));

    for offset in [-9, -1, 0, 1, 9] {
        let key = RotationKey::generate(secret, offset, 3.2f64.sqrt(), rng);
        let actual = encoder
            .decode_complex(&CkksEngine::decrypt(&initial.rotate(&key), secret));
        for (j, a) in actual.iter().enumerate() {
            assert!(
                (*a - values
                    [(j as i32 + offset).rem_euclid((N / 2) as i32) as usize])
                    .norm()
                    < 1e-3,
                "N={N} offset={offset} slot={j} actual={a}"
            );
        }
    }
    let key = ConjugationKey::generate(secret, 3.2f64.sqrt(), rng);
    let actual = encoder.decode_complex(&CkksEngine::decrypt(
        &initial.conjugate(&key).unwrap(),
        secret,
    ));
    for (a, v) in actual.iter().zip(values) {
        assert!((*a - v.conj()).norm() < 1e-3);
    }
}
