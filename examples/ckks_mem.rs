#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    PolyRing, PublicKey, PublicKeyParams, SecretKey, SecretKeyParams, decrypt,
    encoding, encrypt,
};

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _dhat = dhat::Profiler::new_heap();

    println!("Beginning memory profiling for large plaintext operations");

    // Parameters setup
    let ring_degree = 8192 * 2;
    let scale_bits = 40;
    let modulus = (1u64 << 60) - 1;
    let mut rng = ChaCha20Rng::seed_from_u64(123);

    // Create keys
    let sk_params = SecretKeyParams {
        ring_degree,
        modulus,
        hamming_weight: 100,
    };
    println!("Generating secret key...");
    let secret_key = SecretKey::generate(&sk_params, &mut rng);

    let pk_params = PublicKeyParams {
        poly_len: ring_degree,
        modulus,
        error_variance: 3.0,
    };
    println!("Generating public key...");
    let public_key = PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng);

    // Generate large test vector (10K elements)
    println!("Generating 10K test values...");
    let large_values: Vec<f64> = (0..8000).map(|i| (i as f64) * 0.01).collect();

    println!("Secret key ring_dim: {}", secret_key.s.ring_dim());
    println!("Public key a ring_dim: {}", public_key.a.ring_dim());
    println!("Public key b ring_dim: {}", public_key.b.ring_dim());

    // Encoding parameters
    let encoding_params = encoding::EncodingParams::new(ring_degree, scale_bits)
        .expect("Failed to create encoding parameters");

    // Profile encoding
    println!("Profiling encoding...");
    let coeffs =
        encoding::encode(&large_values, &encoding_params).expect("Encoding failed");

    // Profile polynomial creation
    println!("Creating polynomial...");
    let poly = PolyRing::from_coeffs(&coeffs, modulus, ring_degree);
    let scale = (1u64 << scale_bits) as f64;

    // Profile encryption
    println!("Profiling encryption...");
    let ciphertext = encrypt(&poly, &public_key, scale, &mut rng);

    // Profile decryption
    println!("Profiling decryption...");
    let decrypted_poly = decrypt(&ciphertext, &secret_key);

    // Transform back to coefficients
    println!("Converting to coefficients...");
    let decrypted_coeffs = poly_to_coeffs(&decrypted_poly);

    // Profile decoding
    println!("Profiling decoding...");
    let _result = encoding::decode(&decrypted_coeffs, &encoding_params)
        .expect("Decoding failed");

    println!("Memory profiling complete!");
}

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
