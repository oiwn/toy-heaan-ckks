/* use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    PolyRing, PublicKey, PublicKeyParams, RelinearizationKey,
    RelinearizationKeyParams, SecretKey, SecretKeyParams, decrypt, encoding,
    encrypt,
};

fn main() -> Result<(), String> {
    // Parameters setup
    let ring_dim = 8; // Small for example
    let scale_bits = 40;
    let modulus = (1u64 << 60) - 1; // Large prime-like modulus
    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    // Create Secret-Key
    let sk_params = SecretKeyParams {
        ring_dim,
        modulus,
        hamming_weight: 2,
    };
    let secret_key = SecretKey::generate(&sk_params, &mut rng);

    let pk_params = PublicKeyParams {
        ring_dim,
        modulus,
        error_variance: 3.0,
    };
    let public_key = PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng);

    let relin_params: RelinearizationKeyParams = (&pk_params).into();
    let relin_key =
        RelinearizationKey::from_secret_key(&secret_key, &relin_params, &mut rng);

    // Original values
    let values1 = vec![1.5, 2.5, 3.5, 4.5];
    let values2 = vec![0.5, 1.0, 1.5, 2.0];
    let values3 = vec![1.0, 2.0, 3.0, 4.0];

    // Parameters for encoding
    let encoding_params = encoding::EncodingParams::new(ring_dim, scale_bits)?;

    // Encode
    let coeffs1 = encoding::encode(&values1, &encoding_params)?;
    let coeffs2 = encoding::encode(&values2, &encoding_params)?;
    let coeffs3 = encoding::encode(&values3, &encoding_params)?;

    // Convert to polynomial (you might need to add a helper function)
    let poly1 = PolyRing::from_coeffs(&coeffs1, modulus, ring_dim);
    let poly2 = PolyRing::from_coeffs(&coeffs2, modulus, ring_dim);
    let poly3 = PolyRing::from_coeffs(&coeffs3, modulus, ring_dim);

    // Encrypt
    let scale = (1u64 << scale_bits) as f64;
    let ct1 = encrypt(&poly1, &public_key, scale, &mut rng);
    let ct2 = encrypt(&poly2, &public_key, scale, &mut rng);
    let ct3 = encrypt(&poly3, &public_key, scale, &mut rng);

    // Homomorphic Addition and multiplication
    let ct_sum = ct1.add(&ct2);
    let ct_mul = ct_sum.mul(&ct3, &relin_key);

    // Decrypt
    let decrypted_poly_sum = decrypt(&ct_sum, &secret_key);
    let decrypted_poly_mul = decrypt(&ct_mul, &secret_key);

    // Convert back to coefficients
    let decrypted_coeffs = poly_to_coeffs(&decrypted_poly_mul);

    // Decode
    let result = encoding::decode(&decrypted_coeffs, &encoding_params)?;

    // Print results
    println!("values_1: {:?} poly: {}", values1, poly1);
    println!("values_2: {:?} poly: {}", values2, poly2);
    println!("values_3: {:?} poly: {}", values3, poly3);
    println!("Result of (values1 + values2) * values3:");
    println!("values1 + values2 = {}", decrypted_poly_sum);
    println!("(values1 + values2) * values3 = {}", decrypted_poly_mul);
    println!("Decrypted result: {:?}", result);
    println!("Reference i64 max: {}", i64::MAX);

    Ok(())
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
} */

fn main() {
    todo!();
}
