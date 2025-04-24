use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    PolyRing, PublicKey, PublicKeyParams, SecretKey, SecretKeyParams, decrypt,
    encoding, encrypt,
};

fn main() -> Result<(), String> {
    // Parameters setup
    let ring_degree = 8;
    let scale_bits = 30;
    let modulus = (1u64 << 60) - 1; // Large prime-like modulus
    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    // Create Secret-Key
    let sk_params = SecretKeyParams {
        ring_degree,
        modulus,
        hamming_weight: 4,
    };
    let secret_key = SecretKey::generate(&sk_params, &mut rng);

    // Create Public-Key
    let pk_params = PublicKeyParams {
        poly_len: ring_degree,
        modulus,
        error_variance: 3.0,
    };
    let public_key = PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng);

    // Original values
    let values1 = vec![1.5, 2.5, 3.5, 4.5];
    let values2 = vec![0.5, 1.0, 1.5, 2.0];

    // Parameters for encoding
    let encoding_params = encoding::EncodingParams::new(ring_degree, scale_bits)?;

    // Encode to polynomial
    let coeffs1 = encoding::encode(&values1, &encoding_params)?;
    let coeffs2 = encoding::encode(&values2, &encoding_params)?;

    // Convert to polynomial (you might need to add a helper function)
    let poly1 = PolyRing::from_coeffs(&coeffs1, modulus, ring_degree);
    let poly2 = PolyRing::from_coeffs(&coeffs2, modulus, ring_degree);

    // Encrypt
    let scale = (1u64 << scale_bits) as f64;
    let ct1 = encrypt(&poly1, &public_key, scale, &mut rng);
    let ct2 = encrypt(&poly2, &public_key, scale, &mut rng);

    // Homomorphic Addition
    let ct_sum = ct1.add(&ct2);

    // Decrypt
    let decrypted_poly = decrypt(&ct_sum, &secret_key);

    // Convert back to coefficients
    let decrypted_coeffs = poly_to_coeffs(&decrypted_poly);

    // Decode
    let result = encoding::decode(&decrypted_coeffs, &encoding_params)?;

    // Print results
    println!("values_1: {:?}", values1);
    println!("values_2: {:?}", values2);
    println!(
        "Expected sum: {:?}",
        values1
            .iter()
            .zip(values2.iter())
            .map(|(a, b)| a + b)
            .collect::<Vec<_>>()
    );
    println!("Decrypted result: {:?}", result);

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
}
