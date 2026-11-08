//! Scheme operations and bootstrap per backend and degree.
//!
//! Every group checks its results once, outside timing, so a fast but wrong
//! path cannot produce numbers.

#[macro_use]
mod common;

use std::hint::black_box;
use std::time::Duration;

use common::Backend;
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use num_complex::Complex64;
use toy_heaan_ckks::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::{ConjugationKey, RelinearizationKey, RotationKey};
use toy_heaan_ckks::ckks::{Plaintext, Scale};
use toy_heaan_ckks::rings::backends::rns_ntt::RnsPoly;

const ERROR_VARIANCE: f64 = 3.2;
const SCHEME_PRIMES: usize = 8;
/// Sanity bound for key-switched results: a wrong path is off by ~1, while
/// BV key-switch noise grows with N (1.8e-4 at N=2048, scale 2^50; see
/// `specs/research/large-degree-bootstrap.md`). Not a precision gate.
const KEY_SWITCH_SANITY: f64 = 1e-2;

fn values(slots: usize) -> Vec<Complex64> {
    (0..slots)
        .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
        .collect()
}

fn max_error(actual: &[Complex64], expected: &[Complex64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max)
}

/// Encode/encrypt, arithmetic, key switching and key generation on an
/// 8-prime 61-bit chain at scale 2^50.
fn scheme_ops<P: Backend<N>, const N: usize>(c: &mut Criterion) {
    let top = P::chain(&common::primes(61, SCHEME_PRIMES, N));
    let scale = Scale::power_of_two(50).expect("scale");
    let engine = CkksEngine::<P, N>::new(
        top.clone(),
        CkksParams {
            error_variance: ERROR_VARIANCE,
            hamming_weight: N / 2,
            scale,
        },
    );
    let std_dev = ERROR_VARIANCE.sqrt();
    let mut rng = common::rng(N as u64);
    let secret = engine.generate_secret_key(&mut rng).expect("secret");
    let public = engine
        .generate_public_key(&secret, &mut rng)
        .expect("public");
    let relin = engine.generate_relin_key(&secret, &mut rng);
    let rotation = engine.generate_rotation_key(&secret, 1, &mut rng);
    let conjugation = ConjugationKey::generate(&secret, std_dev, &mut rng);

    let encoder = CkksEncoder::new(scale);
    let input = values(N / 2);
    let plain = encoder.encode_complex(&input, top.clone());
    let ct = engine.encrypt(&plain, &public, &mut rng);
    let other = engine.encrypt(&plain, &public, &mut rng);
    let product = ct.mul(&other, &relin);

    // Correctness once, outside timing.
    let decrypt = |ct: &toy_heaan_ckks::ckks::Ciphertext<P, N>| {
        let sk = secret.mod_switch_to(ct.level()).expect("level");
        CkksEncoder::new(ct.scale).decode_complex(&CkksEngine::decrypt(ct, &sk))
    };
    assert!(max_error(&decrypt(&ct), &input) < 1e-6, "encrypt/decrypt");
    let squares: Vec<_> = input.iter().map(|x| x * x).collect();
    let rescaled = product.rescale().expect("rescale");
    let mul_error = max_error(&decrypt(&rescaled), &squares);
    assert!(
        mul_error < KEY_SWITCH_SANITY,
        "mul N={N} error {mul_error:e}"
    );
    let rotated_plain = Plaintext {
        poly: plain.poly.rotate_slots(1),
        scale: plain.scale,
        slots: plain.slots,
    };
    let rotated = encoder.decode_complex(&rotated_plain);
    let rotate_error = max_error(&decrypt(&ct.rotate(&rotation)), &rotated);
    assert!(
        rotate_error < KEY_SWITCH_SANITY,
        "rotate N={N} error {rotate_error:e}"
    );
    let conjugates: Vec<_> = input.iter().map(Complex64::conj).collect();
    let conjugated = ct.conjugate(&conjugation).expect("conjugate");
    let conjugate_error = max_error(&decrypt(&conjugated), &conjugates);
    assert!(
        conjugate_error < KEY_SWITCH_SANITY,
        "conjugate N={N} error {conjugate_error:e}"
    );

    let mut group = c.benchmark_group("scheme");
    let id = |op: &str| common::id::<P, N>(op, 61, SCHEME_PRIMES);
    group.bench_function(id("encode"), |b| {
        b.iter(|| black_box(encoder.encode_complex::<P>(&input, top.clone())));
    });
    group.bench_function(id("decode"), |b| {
        b.iter(|| black_box(encoder.decode_complex(&plain)));
    });
    group.bench_function(id("encrypt"), |b| {
        b.iter(|| black_box(engine.encrypt(&plain, &public, &mut rng)));
    });
    group.bench_function(id("decrypt"), |b| {
        b.iter(|| black_box(CkksEngine::decrypt(&ct, &secret)));
    });
    group.bench_function(id("add"), |b| {
        b.iter(|| black_box(ct.add(&other).expect("add")));
    });
    group.bench_function(id("mul_relin"), |b| {
        b.iter(|| black_box(ct.mul(&other, &relin)));
    });
    group.bench_function(id("rescale"), |b| {
        b.iter(|| black_box(product.rescale().expect("rescale")));
    });
    group.bench_function(id("rotate"), |b| {
        b.iter(|| black_box(ct.rotate(&rotation)));
    });
    group.bench_function(id("conjugate"), |b| {
        b.iter(|| black_box(ct.conjugate(&conjugation).expect("conjugate")));
    });
    group.bench_function(id("keygen_relin"), |b| {
        b.iter(|| {
            black_box(RelinearizationKey::generate(&secret, std_dev, &mut rng))
        });
    });
    group.bench_function(id("keygen_rotation"), |b| {
        b.iter(|| black_box(RotationKey::generate(&secret, 1, std_dev, &mut rng)));
    });
    group.finish();
}

/// Bootstrap config for degree `N`, matching `examples/bootstrap_large.rs`:
/// the reference (weight, ratio and double angles chosen from N) with
/// grouped transforms.
fn bootstrap_config<const N: usize>(
    cts: Option<usize>,
    stc: Option<usize>,
) -> BootstrapConfig {
    let stages = (N / 2).ilog2() as usize;
    BootstrapConfig {
        coeffs_to_slots_levels: cts.unwrap_or(stages.div_ceil(2)),
        slots_to_coeffs_levels: stc.unwrap_or(stages.div_ceil(3)),
        ..BootstrapConfig::reference(N).expect("reference")
    }
}

/// Primes for `config`: its depth, the bottom prime and one spare level.
fn bootstrap_primes(config: &BootstrapConfig) -> usize {
    let polynomial_depth = config.polynomial_degree.ilog2() as usize + 1;
    config.coeffs_to_slots_levels
        + polynomial_depth
        + config.double_angles as usize
        + config.slots_to_coeffs_levels
        + 2
}

/// Plan construction, key generation and one bootstrap.
fn bootstrap<P: Backend<N>, const N: usize>(
    c: &mut Criterion,
    group_name: &str,
    cts: Option<usize>,
    stc: Option<usize>,
    full_pipeline: bool,
) {
    let config = bootstrap_config::<N>(cts, stc);
    let top = P::chain(&common::primes(61, bootstrap_primes(&config), N));
    let plan =
        BootstrapPlan::<P, N>::new(top.clone(), config.clone()).expect("plan");
    let message_scale = plan.schedule().message_scale();
    let engine = CkksEngine::<P, N>::new(
        top.clone(),
        CkksParams {
            error_variance: ERROR_VARIANCE,
            hamming_weight: config.max_secret_weight,
            scale: message_scale,
        },
    );
    let std_dev = ERROR_VARIANCE.sqrt();
    let mut rng = common::rng(N as u64);
    let secret = engine.generate_secret_key(&mut rng).expect("secret");
    let public = engine
        .generate_public_key(&secret, &mut rng)
        .expect("public");
    let keys = plan
        .generate_keys(&secret, std_dev, &mut rng)
        .expect("keys");
    let input = values(N / 2);
    let plain = CkksEncoder::new(message_scale).encode_complex(&input, top.clone());
    let ct = engine.encrypt(&plain, &public, &mut rng);

    // Correctness once, outside timing: within the 1e-3 bootstrap gate.
    let refreshed = plan.bootstrap(&ct, &keys).expect("bootstrap");
    let sk = secret.mod_switch_to(refreshed.level()).expect("level");
    let decoded = CkksEncoder::new(refreshed.scale)
        .decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
    let error = max_error(&decoded, &input);
    assert!(error < 1e-3, "bootstrap N={N} error {error:e}");

    let label = format!(
        "{}/N{N}/cts{}/stc{}",
        P::NAME,
        config.coeffs_to_slots_levels,
        config.slots_to_coeffs_levels
    );
    let mut group = c.benchmark_group(group_name);
    group.sample_size(10);
    // Ten samples of the slowest step (N=1024 bootstrap ~3 s).
    let seconds = match (N >= 1024, common::full()) {
        (true, _) => 60,
        (false, true) => 15,
        (false, false) => 6,
    };
    group.measurement_time(Duration::from_secs(seconds));
    if full_pipeline {
        group.bench_function(BenchmarkId::new("plan", &label), |b| {
            b.iter(|| {
                black_box(BootstrapPlan::<P, N>::new(top.clone(), config.clone()))
            });
        });
        group.bench_function(BenchmarkId::new("keygen", &label), |b| {
            b.iter(|| black_box(plan.generate_keys(&secret, std_dev, &mut rng)));
        });
    }
    group.bench_function(BenchmarkId::new("bootstrap", &label), |b| {
        b.iter(|| black_box(plan.bootstrap(&ct, &keys).expect("bootstrap")));
    });
    group.finish();
}

fn bootstrap_default<P: Backend<N>, const N: usize>(c: &mut Criterion) {
    bootstrap::<P, N>(c, "bootstrap", None, None, true);
}

/// Level budgets at N=256 (7 butterfly stages): one, two and seven levels.
fn bootstrap_budgets<P: Backend<N>, const N: usize>(c: &mut Criterion) {
    let stages = (N / 2).ilog2() as usize;
    for budget in [1, 2, stages] {
        bootstrap::<P, N>(
            c,
            "bootstrap_budgets",
            Some(budget),
            Some(budget),
            false,
        );
    }
}

fn benches(c: &mut Criterion) {
    for_degrees!(c, scheme_ops, RnsPoly, [1024, 2048]);
    for_degrees!(c, bootstrap_default, RnsPoly, [64, 256]);
    if common::full() {
        for_degrees!(c, scheme_ops, RnsPoly, [4096, 8192]);
        for_degrees!(c, bootstrap_default, RnsPoly, [1024]);
        for_degrees!(c, bootstrap_budgets, RnsPoly, [256]);
    }
}

criterion_group! {
    name = scheme;
    config = common::config();
    targets = benches
}
criterion_main!(scheme);
