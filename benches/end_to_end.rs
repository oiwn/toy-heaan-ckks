use criterion::{Criterion, black_box, criterion_group, criterion_main};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    PolyRing, PublicKey, PublicKeyParams, SecretKey, SecretKeyParams, decode,
    decrypt, encode, encoding, encrypt,
};

fn bench_ckks_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("ckks_pipeline");

    // Parameters
    let ring_dim = 8192; // Ensure this is large enough for our data
    let scale_bits = 40;
    let modulus = (1u64 << 60) - 1;

    // Benchmark with different data sizes
    for &size in &[1000, 5000, 10000] {
        group.bench_function(format!("full_pipeline_{}_elements", size), |b| {
            // Setup
            let mut rng = ChaCha20Rng::seed_from_u64(123);

            // Create keys
            let sk_params = SecretKeyParams {
                ring_dim,
                modulus,
                hamming_weight: 100,
            };
            let secret_key = SecretKey::generate(&sk_params, &mut rng);

            let pk_params = PublicKeyParams {
                ring_dim,
                modulus,
                error_variance: 3.0,
            };
            let public_key =
                PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng);

            // Generate test data
            let values: Vec<f64> = (0..size).map(|i| (i as f64) * 0.01).collect();

            let encoding_params =
                encoding::EncodingParams::new(ring_dim, scale_bits)
                    .expect("Failed to create encoding parameters");

            // Measure the complete pipeline
            b.iter(|| {
                // 1. Encode
                let coeffs = encode(black_box(&values), &encoding_params)
                    .expect("Encoding failed");

                // 2. Create polynomial and encrypt
                let poly = PolyRing::from_coeffs(&coeffs, modulus, ring_dim);
                let scale = (1u64 << scale_bits) as f64;
                let ciphertext = encrypt(&poly, &public_key, scale, &mut rng);

                // 3. Decrypt
                let decrypted_poly = decrypt(&ciphertext, &secret_key);

                // 4. Convert to coefficients and decode
                let decrypted_coeffs = poly_to_coeffs(&decrypted_poly);
                let result = decode(&decrypted_coeffs, &encoding_params)
                    .expect("Decoding failed");

                // Return the result to prevent the compiler from optimizing it away
                black_box(result)
            });
        });
    }

    group.finish();
}

// Helper function to convert polynomial to signed coefficients
fn poly_to_coeffs(poly: &PolyRing) -> Vec<i64> {
    let modulus = poly.modulus();
    let half_modulus = modulus / 2;

    poly.into_iter()
        .map(|&c| {
            if c > half_modulus {
                -((modulus - c) as i64)
            } else {
                c as i64
            }
        })
        .collect()
}

criterion_group!(benches, bench_ckks_pipeline);
criterion_main!(benches);
