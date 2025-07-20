use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use crypto_bigint::{NonZero, U256};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::hint::black_box;
use std::sync::Arc;
use toy_heaan_ckks::rings::backends::rns::{RnsBasisBuilder, RnsPolyRing};
use toy_heaan_ckks::{
    CkksEngine, NaivePolyRing, Plaintext, PolyRing, PolyRingU256,
};

// Benchmark parameters
// const DEGREES: &[usize] = &[256, 512, 1024, 2048];
const SCALE_BITS: u32 = 40;

// Test moduli
const NAIVE_MODULUS: u64 = 741507920154517877; // ~60 bits
const U256_MODULUS_WORDS: [u64; 4] = [0x1, 0x0, 0x0, 0x8000000000000000]; // Large prime for U256

fn get_u256_modulus() -> NonZero<U256> {
    NonZero::new(U256::from_words(U256_MODULUS_WORDS)).unwrap()
}

// Benchmark encrypt -> add -> decrypt cycle for Naive backend
fn bench_naive_cycle<const DEGREE: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("naive_cycle_degree_{}", DEGREE));

    let engine = CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::builder()
        .build_naive(NAIVE_MODULUS, SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);

    // Generate keys once
    let sk = engine.generate_secret_key(&mut rng).unwrap();
    let pk = engine.generate_public_key(&sk, &mut rng).unwrap();

    // Create test plaintexts with proper coefficient arrays
    let mut coeffs1 = vec![0i64; DEGREE];
    let mut coeffs2 = vec![0i64; DEGREE];
    coeffs1[0..4].copy_from_slice(&[1, 2, 3, 4]);
    coeffs2[0..4].copy_from_slice(&[5, 6, 7, 8]);

    let plaintext1 = Plaintext {
        poly: NaivePolyRing::from_coeffs(&coeffs1, engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };
    let plaintext2 = Plaintext {
        poly: NaivePolyRing::from_coeffs(&coeffs2, engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    group.bench_function("encrypt_add_decrypt", |b| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);

            // Encrypt
            let ct1 = engine.encrypt(black_box(&plaintext1), &pk, &mut rng);
            let ct2 = engine.encrypt(black_box(&plaintext2), &pk, &mut rng);

            // Homomorphic addition
            let ct_sum =
                CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::add_ciphertexts(
                    &ct1, &ct2,
                );

            // Decrypt
            let result =
                CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::decrypt(&ct_sum, &sk);

            black_box(result);
        });
    });

    group.finish();
}

// Benchmark encrypt -> add -> decrypt cycle for U256 backend
fn bench_u256_cycle<const DEGREE: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("u256_cycle_degree_{}", DEGREE));

    let engine = CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::builder()
        .build_bigint_u256(get_u256_modulus(), SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);

    // Generate keys once
    let sk = engine.generate_secret_key(&mut rng).unwrap();
    let pk = engine.generate_public_key(&sk, &mut rng).unwrap();

    // Create test plaintexts with proper coefficient arrays
    let mut coeffs1 = vec![0i64; DEGREE];
    let mut coeffs2 = vec![0i64; DEGREE];
    coeffs1[0..4].copy_from_slice(&[1, 2, 3, 4]);
    coeffs2[0..4].copy_from_slice(&[5, 6, 7, 8]);

    let plaintext1 = Plaintext {
        poly: PolyRingU256::from_coeffs(&coeffs1, engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };
    let plaintext2 = Plaintext {
        poly: PolyRingU256::from_coeffs(&coeffs2, engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    group.bench_function("encrypt_add_decrypt", |b| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);

            // Encrypt
            let ct1 = engine.encrypt(black_box(&plaintext1), &pk, &mut rng);
            let ct2 = engine.encrypt(black_box(&plaintext2), &pk, &mut rng);

            // Homomorphic addition
            let ct_sum =
                CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::add_ciphertexts(
                    &ct1, &ct2,
                );

            // Decrypt
            let result =
                CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::decrypt(&ct_sum, &sk);

            black_box(result);
        });
    });

    group.finish();
}

// Benchmark encrypt -> add -> decrypt cycle for RNS backend
fn bench_rns_cycle<const DEGREE: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("rns_cycle_degree_{}", DEGREE));

    // Create RNS basis with multiple primes
    let rns_basis = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_prime_bits(vec![17, 19, 23]) // Three primes of different sizes
            .build()
            .unwrap(),
    );

    let engine = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_rns(rns_basis.clone(), SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);

    let sk = engine.generate_secret_key(&mut rng).unwrap();
    let pk = engine.generate_public_key(&sk, &mut rng).unwrap();

    let mut coeffs = vec![0i64; DEGREE];
    coeffs[0..4].copy_from_slice(&[1, 2, 3, 4]);
    let plaintext = Plaintext {
        poly: RnsPolyRing::from_i64_slice(&coeffs, rns_basis.clone()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    group.bench_with_input(BenchmarkId::new("full_cycle", "rns"), &(), |b, _| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);
            let ct1 = engine.encrypt(black_box(&plaintext), &pk, &mut rng);
            let ct2 = engine.encrypt(black_box(&plaintext), &pk, &mut rng);
            let ct_sum = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::add_ciphertexts(
                &ct1, &ct2,
            );
            let result =
                CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::decrypt(&ct_sum, &sk);
            black_box(result);
        });
    });

    group.finish();
}

// Individual operation benchmarks for Naive backend
fn bench_naive_operations<const DEGREE: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("naive_ops_degree_{}", DEGREE));

    let engine = CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::builder()
        .build_naive(NAIVE_MODULUS, SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);

    let sk = engine.generate_secret_key(&mut rng).unwrap();
    let pk = engine.generate_public_key(&sk, &mut rng).unwrap();

    let mut coeffs = vec![0i64; DEGREE];
    coeffs[0..4].copy_from_slice(&[1, 2, 3, 4]);

    let plaintext = Plaintext {
        poly: NaivePolyRing::from_coeffs(&coeffs, engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    // Pre-encrypt ciphertext for decrypt benchmark
    let ciphertext = engine.encrypt(&plaintext, &pk, &mut rng);

    group.bench_function("encrypt_only", |b| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);
            let ct = engine.encrypt(black_box(&plaintext), &pk, &mut rng);
            black_box(ct);
        });
    });

    group.bench_function("decrypt_only", |b| {
        b.iter(|| {
            let result = CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::decrypt(
                black_box(&ciphertext),
                &sk,
            );
            black_box(result);
        });
    });

    group.bench_function("homomorphic_add", |b| {
        b.iter(|| {
            let ct_sum =
                CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::add_ciphertexts(
                    black_box(&ciphertext),
                    black_box(&ciphertext),
                );
            black_box(ct_sum);
        });
    });

    group.finish();
}

// Individual operation benchmarks for U256 backend
fn bench_u256_operations<const DEGREE: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("u256_ops_degree_{}", DEGREE));

    let engine = CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::builder()
        .build_bigint_u256(get_u256_modulus(), SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);

    let sk = engine.generate_secret_key(&mut rng).unwrap();
    let pk = engine.generate_public_key(&sk, &mut rng).unwrap();

    let mut coeffs = vec![0i64; DEGREE];
    coeffs[0..4].copy_from_slice(&[1, 2, 3, 4]);

    let plaintext = Plaintext {
        poly: PolyRingU256::from_coeffs(&coeffs, engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    // Pre-encrypt ciphertext for decrypt benchmark
    let ciphertext = engine.encrypt(&plaintext, &pk, &mut rng);

    group.bench_function("encrypt_only", |b| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);
            let ct = engine.encrypt(black_box(&plaintext), &pk, &mut rng);
            black_box(ct);
        });
    });

    group.bench_function("decrypt_only", |b| {
        b.iter(|| {
            let result = CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::decrypt(
                black_box(&ciphertext),
                &sk,
            );
            black_box(result);
        });
    });

    group.bench_function("homomorphic_add", |b| {
        b.iter(|| {
            let ct_sum =
                CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::add_ciphertexts(
                    black_box(&ciphertext),
                    black_box(&ciphertext),
                );
            black_box(ct_sum);
        });
    });

    group.finish();
}

// Individual operation benchmarks for RNS backend
fn bench_rns_operations<const DEGREE: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("rns_ops_degree_{}", DEGREE));

    // Create RNS basis with multiple primes
    let rns_basis = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_prime_bits(vec![17, 19, 23])
            .build()
            .unwrap(),
    );

    let engine = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_rns(rns_basis.clone(), SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);

    let sk = engine.generate_secret_key(&mut rng).unwrap();
    let pk = engine.generate_public_key(&sk, &mut rng).unwrap();

    let mut coeffs = vec![0i64; DEGREE];
    coeffs[0..4].copy_from_slice(&[1, 2, 3, 4]);

    let plaintext = Plaintext {
        poly: RnsPolyRing::from_i64_slice(&coeffs, rns_basis.clone()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    // Pre-encrypt ciphertext for decrypt benchmark
    let ciphertext = engine.encrypt(&plaintext, &pk, &mut rng);

    group.bench_function("encrypt_only", |b| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);
            let ct = engine.encrypt(black_box(&plaintext), &pk, &mut rng);
            black_box(ct);
        });
    });

    group.bench_function("decrypt_only", |b| {
        b.iter(|| {
            let result = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::decrypt(
                black_box(&ciphertext),
                &sk,
            );
            black_box(result);
        });
    });

    group.bench_function("homomorphic_add", |b| {
        b.iter(|| {
            let ct_sum = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::add_ciphertexts(
                black_box(&ciphertext),
                black_box(&ciphertext),
            );
            black_box(ct_sum);
        });
    });

    group.finish();
}

// Comparison benchmark - same test for both backends
fn bench_backend_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("backend_comparison_degree_1024");

    const DEGREE: usize = 1024;

    // Setup Naive
    let naive_engine = CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::builder()
        .build_naive(NAIVE_MODULUS, SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);
    let naive_sk = naive_engine.generate_secret_key(&mut rng).unwrap();
    let naive_pk = naive_engine
        .generate_public_key(&naive_sk, &mut rng)
        .unwrap();
    let mut naive_coeffs = vec![0i64; DEGREE];
    naive_coeffs[0..4].copy_from_slice(&[1, 2, 3, 4]);
    let naive_plaintext = Plaintext {
        poly: NaivePolyRing::from_coeffs(&naive_coeffs, naive_engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    // Setup U256
    let u256_engine = CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::builder()
        .build_bigint_u256(get_u256_modulus(), SCALE_BITS)
        .unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(42);
    let u256_sk = u256_engine.generate_secret_key(&mut rng).unwrap();
    let u256_pk = u256_engine.generate_public_key(&u256_sk, &mut rng).unwrap();
    let mut u256_coeffs = vec![0i64; DEGREE];
    u256_coeffs[0..4].copy_from_slice(&[1, 2, 3, 4]);
    let u256_plaintext = Plaintext {
        poly: PolyRingU256::from_coeffs(&u256_coeffs, u256_engine.context()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    group.bench_with_input(BenchmarkId::new("full_cycle", "naive"), &(), |b, _| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);
            let ct1 = naive_engine.encrypt(
                black_box(&naive_plaintext),
                &naive_pk,
                &mut rng,
            );
            let ct2 = naive_engine.encrypt(
                black_box(&naive_plaintext),
                &naive_pk,
                &mut rng,
            );
            let ct_sum =
                CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::add_ciphertexts(
                    &ct1, &ct2,
                );
            let result = CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::decrypt(
                &ct_sum, &naive_sk,
            );
            black_box(result);
        });
    });

    group.bench_with_input(BenchmarkId::new("full_cycle", "u256"), &(), |b, _| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);
            let ct1 =
                u256_engine.encrypt(black_box(&u256_plaintext), &u256_pk, &mut rng);
            let ct2 =
                u256_engine.encrypt(black_box(&u256_plaintext), &u256_pk, &mut rng);
            let ct_sum =
                CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::add_ciphertexts(
                    &ct1, &ct2,
                );
            let result = CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::decrypt(
                &ct_sum, &u256_sk,
            );
            black_box(result);
        });
    });

    let rns_basis = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_prime_bits(vec![17, 19, 23])
            .build()
            .unwrap(),
    );

    let rns_engine = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_rns(rns_basis.clone(), SCALE_BITS)
        .unwrap();

    let rns_sk = rns_engine.generate_secret_key(&mut rng).unwrap();
    let rns_pk = rns_engine.generate_public_key(&rns_sk, &mut rng).unwrap();

    let mut rns_coeffs = vec![0i64; DEGREE];
    rns_coeffs[0..4].copy_from_slice(&[1, 2, 3, 4]);
    let rns_plaintext = Plaintext {
        poly: RnsPolyRing::from_i64_slice(&rns_coeffs, rns_basis.clone()),
        scale: 2.0_f64.powi(SCALE_BITS as i32),
    };

    group.bench_with_input(BenchmarkId::new("full_cycle", "rns"), &(), |b, _| {
        b.iter(|| {
            let mut rng = ChaCha8Rng::seed_from_u64(42);
            let ct1 =
                rns_engine.encrypt(black_box(&rns_plaintext), &rns_pk, &mut rng);
            let ct2 =
                rns_engine.encrypt(black_box(&rns_plaintext), &rns_pk, &mut rng);
            let ct_sum = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::add_ciphertexts(
                &ct1, &ct2,
            );
            let result = CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::decrypt(
                &ct_sum, &rns_sk,
            );
            black_box(result);
        });
    });

    group.finish();
}

// Macro to generate benchmark functions for all degrees
macro_rules! bench_all_degrees {
    ($bench_fn:ident, $($degree:expr),+) => {
        $(
            paste::paste! {
                #[allow(dead_code)]
                fn [<$bench_fn _ $degree>](c: &mut Criterion) {
                    $bench_fn::<$degree>(c);
                }
            }
        )+
    };
}

// Generate benchmark functions for all degrees
bench_all_degrees!(bench_naive_cycle, 256, 512, 1024, 2048);
bench_all_degrees!(bench_u256_cycle, 256, 512, 1024, 2048);
bench_all_degrees!(bench_naive_operations, 256, 512, 1024, 2048);
bench_all_degrees!(bench_u256_operations, 256, 512, 1024, 2048);
bench_all_degrees!(bench_rns_cycle, 256, 512, 1024, 2048);
bench_all_degrees!(bench_rns_operations, 256, 512, 1024, 2048);

criterion_group!(
    benches,
    // Full cycle benchmarks for Naive
    bench_naive_cycle_256,
    bench_naive_cycle_512,
    bench_naive_cycle_1024,
    bench_naive_cycle_2048,
    // Full cycle benchmarks for U256
    bench_u256_cycle_256,
    bench_u256_cycle_512,
    // bench_u256_cycle_1024,
    // bench_u256_cycle_2048,
    // Full cycle benchmarks for RNS
    bench_rns_cycle_256,
    bench_rns_cycle_512,
    bench_rns_cycle_1024,
    bench_rns_cycle_2048,
    // Individual operation benchmarks for Naive
    bench_naive_operations_256,
    bench_naive_operations_512,
    bench_naive_operations_1024,
    bench_naive_operations_2048,
    // Individual operation benchmarks for U256
    bench_u256_operations_256,
    bench_u256_operations_512,
    // bench_u256_operations_1024,
    // bench_u256_operations_2048,
    // Individual operation benchmarks for RNS
    bench_rns_operations_256,
    bench_rns_operations_512,
    bench_rns_operations_1024,
    bench_rns_operations_2048,
    // Comparison benchmark
    bench_backend_comparison
);

criterion_main!(benches);
