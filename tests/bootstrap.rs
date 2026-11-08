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
    let config = BootstrapConfig::reference(N).unwrap();
    let plan =
        BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), config.clone()).unwrap();
    let message = plan.schedule().message_scale();
    let encoder = CkksEncoder::new(message);
    let engine = CkksEngine::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: config.max_secret_weight,
            scale: message,
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
            let not_power_of_two = Scale::new(1.5 * message.value()).unwrap();
            let above = message.checked_mul(pow2(1)).unwrap();
            for scale in [pow2(0), above, not_power_of_two] {
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
                let lift = (message.log2() - squared.scale.log2()).round() as i32;
                let explicit =
                    squared.mul_const(1.0, pow2(lift)).mod_switch_to(1).unwrap();

                let refreshed = plan.bootstrap(&squared, &keys).unwrap();
                assert_same(&squared, &snapshot);
                if seed == 0 && *scenario == "zero" && cycle == 0 {
                    assert_same(
                        &refreshed,
                        &plan.bootstrap(&explicit, &keys).unwrap(),
                    );
                }
                // Output properties come from the schedule; the scale is
                // tracked, not snapped, and the input is lifted by 2^lift.
                let schedule = plan.schedule();
                assert_eq!(refreshed.level(), schedule.output_level());
                assert!(
                    refreshed.scale.is_compatible(
                        schedule.output_scale(explicit.scale).unwrap()
                    )
                );
                assert!((refreshed.scale.log2() - message.log2()).abs() < 1e-9);
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

/// Multi-group transforms: several (CtS, StC) budgets, one seed, two
/// scenarios, two cycles each, against the same 1e-3 gate.
fn multi_group<const N: usize>() {
    let primes = generate_primes(61, 18, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).unwrap());
    let reference = BootstrapConfig::reference(N).unwrap();
    let reference_plan =
        BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), reference.clone())
            .unwrap();
    let message = reference_plan.schedule().message_scale();
    let encoder = CkksEncoder::new(message);
    let engine = CkksEngine::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: reference.max_secret_weight,
            scale: message,
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let secret = engine.generate_secret_key(&mut rng).unwrap();
    let public = engine.generate_public_key(&secret, &mut rng).unwrap();
    let reference_keys = reference_plan
        .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
        .unwrap();

    for (cts, stc) in [(2, 2), (3, 3), (1, 3), (3, 1)] {
        let config = BootstrapConfig {
            coeffs_to_slots_levels: cts,
            slots_to_coeffs_levels: stc,
            ..reference.clone()
        };
        let plan =
            BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), config).unwrap();
        let schedule = plan.schedule();
        let polynomial_depth = reference.polynomial_degree.ilog2() as usize + 1;
        assert_eq!(
            schedule.depth(),
            cts + polynomial_depth + reference.double_angles as usize + stc
        );
        let keys = plan
            .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
            .unwrap();

        for scenario in ["asymmetric", "random_disk"] {
            let values: Vec<_> = (0..N / 2)
                .map(|j| match scenario {
                    "asymmetric" => Complex64::from_polar(0.8, j as f64 * 0.7),
                    _ => Complex64::from_polar(
                        rng.random::<f64>().sqrt(),
                        rng.random::<f64>() * std::f64::consts::TAU,
                    ),
                })
                .collect();
            let initial = engine.encrypt(
                &encoder.encode_complex(&values, basis.clone()),
                &public,
                &mut rng,
            );
            // Keys for another schedule miss this plan's levels (deeper CtS
            // shifts every later stage); the check runs before any
            // evaluation and names the kind and level of the first gap.
            let error = plan
                .bootstrap(&initial, &reference_keys)
                .expect_err("reference keys lack this schedule's keys")
                .to_string();
            assert!(
                error.contains("missing") && error.contains("at level"),
                "{error}"
            );

            let mut current = initial;
            let mut expected = values;
            for cycle in 0..2 {
                let refreshed = plan.bootstrap(&current, &keys).unwrap();
                assert_eq!(refreshed.level(), schedule.output_level());
                // Bootstrap lifts the input by 2^k to within half a bit of
                // the message scale.
                let lift = (message.log2() - current.scale.log2()).round() as i32;
                let prepared = current.scale.checked_mul(pow2(lift)).unwrap();
                assert!(
                    refreshed
                        .scale
                        .is_compatible(schedule.output_scale(prepared).unwrap())
                );
                let sk = secret.mod_switch_to(refreshed.level()).unwrap();
                let actual =
                    encoder.decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
                for (a, want) in actual.iter().zip(&expected) {
                    assert!(
                        (*a - *want).norm() < 1e-3,
                        "N={N} budgets=({cts},{stc}) case={scenario} cycle={cycle} \
                         got={a} expected={want}"
                    );
                }
                let relin =
                    RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
                current = refreshed.mul(&refreshed, &relin).rescale().unwrap();
                for x in &mut expected {
                    *x *= *x;
                }
            }
        }
    }
}

#[test]
fn multi_group_bootstrap_16() {
    multi_group::<16>();
}
#[test]
fn multi_group_bootstrap_64() {
    multi_group::<64>();
}

/// `count` primes of `bits` bits.
fn uniform<const N: usize>(bits: usize, count: usize) -> Vec<u64> {
    generate_primes(bits, count, N as u64)
}

/// Settings outside the reference: one seed, the asymmetric case and two
/// cycles against the same 1e-3 gate.
fn smoke<const N: usize>(
    primes: Vec<u64>,
    config: BootstrapConfig,
    message_bits: i32,
) {
    let bits: Vec<_> = primes.iter().map(|p| 64 - p.leading_zeros()).collect();
    let basis = Arc::new(RnsBasis::<N>::new(primes).unwrap());
    let plan = BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), config.clone())
        .unwrap_or_else(|e| panic!("N={N} bits={bits:?} {config:?}: {e}"));
    let schedule = plan.schedule();
    let message_scale = pow2(message_bits);
    assert_eq!(schedule.message_scale(), message_scale);
    let encoder = CkksEncoder::new(message_scale);
    let engine = CkksEngine::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: config.max_secret_weight,
            scale: message_scale,
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let secret = engine.generate_secret_key(&mut rng).unwrap();
    let public = engine.generate_public_key(&secret, &mut rng).unwrap();
    let keys = plan
        .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
        .unwrap();
    let mut expected: Vec<_> = (0..N / 2)
        .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
        .collect();
    let mut current = engine.encrypt(
        &encoder.encode_complex(&expected, basis.clone()),
        &public,
        &mut rng,
    );
    for cycle in 0..2 {
        let refreshed = plan.bootstrap(&current, &keys).unwrap();
        assert_eq!(refreshed.level(), schedule.output_level());
        let lift = (f64::from(message_bits) - current.scale.log2()).round() as i32;
        let prepared = current.scale.checked_mul(pow2(lift)).unwrap();
        assert!(
            refreshed
                .scale
                .is_compatible(schedule.output_scale(prepared).unwrap())
        );
        let sk = secret.mod_switch_to(refreshed.level()).unwrap();
        let actual = CkksEncoder::new(refreshed.scale)
            .decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
        let error = actual
            .iter()
            .zip(&expected)
            .map(|(a, b)| (*a - *b).norm())
            .fold(0f64, f64::max);
        // Shown with `--nocapture`; the validation matrix records these.
        println!(
            "smoke N={N} primes={} weight={} poly={} angles={} cycle={cycle} \
             max_error={error:.1e}",
            bits.len(),
            config.max_secret_weight,
            config.polynomial_degree,
            config.double_angles
        );
        assert!(
            error < 1e-3,
            "N={N} bits={bits:?} {config:?} cycle={cycle} max_error={error}"
        );
        let relin = RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
        current = refreshed.mul(&refreshed, &relin).rescale().unwrap();
        for x in &mut expected {
            *x *= *x;
        }
    }
}

#[test]
fn general_configs() {
    let reference = BootstrapConfig::reference(16).unwrap();
    smoke::<32>(
        uniform::<32>(61, 18),
        BootstrapConfig::reference(32).unwrap(),
        44,
    );
    smoke::<16>(
        uniform::<16>(60, 18),
        BootstrapConfig {
            working_scale: pow2(60),
            ..reference.clone()
        },
        43,
    );
    smoke::<16>(
        uniform::<16>(61, 18),
        BootstrapConfig {
            polynomial_degree: 63,
            double_angles: 3,
            ..reference.clone()
        },
        44,
    );
    // Secret weights below and at the ring degree.
    for weight in [4, 16] {
        smoke::<16>(
            uniform::<16>(61, 18),
            BootstrapConfig {
                max_secret_weight: weight,
                ..reference.clone()
            },
            44,
        );
    }
    // A longer chain: the input enters four levels above the minimum.
    smoke::<16>(uniform::<16>(61, 22), reference.clone(), 44);
    // Mixed widths: a 61-bit bottom prime under 60-bit primes, working
    // scale 2^60; the message scale follows the bottom prime.
    let mut mixed = uniform::<16>(61, 1);
    mixed.extend(uniform::<16>(60, 17));
    smoke::<16>(
        mixed,
        BootstrapConfig {
            working_scale: pow2(60),
            ..reference
        },
        44,
    );
}

/// The degree-dependent reference above N=64: ratio 10, message 2^51.
#[test]
fn reference_128() {
    smoke::<128>(
        uniform::<128>(61, 19),
        BootstrapConfig::reference(128).unwrap(),
        51,
    );
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
