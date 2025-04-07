use toy_heaan_ckks::{
    PolyRing, PublicKey, PublicKeyParams, SecretKey, SecretKeyParams, decrypt,
    encoding, encrypt,
};

fn main() -> Result<(), String> {
    // Parameters setup
    let ring_degree = 4; // Small for example
    let scale_bits = 30;
    let modulus = (1u64 << 60) - 1; // Large prime-like modulus

    // Create Secret-Key
    let sk_params = SecretKeyParams {
        ring_degree,
        modulus,
        sparsity: 0.5,
    };
    let secret_key = SecretKey::generate(&sk_params);

    let pk_params = PublicKeyParams {
        n: ring_degree,
        modulus,
        error_variance: 3.0,
        relin_base: 3, // Base for relinearization key decomposition
        relin_components: 1, // Number of components for our toy example
    };
    let public_key = PublicKey::from_secret_key(&secret_key, &pk_params);

    // Original values
    let values1 = vec![3.5, 3.0, 4.0, 3.5];
    let values2 = vec![1.5, 3.0, 3.0, 4.5];
    let values3 = vec![1.0, 2.0, 3.0, 4.0];

    // Parameters for encoding
    let encoding_params =
        encoding::EncodingParams::new(ring_degree * 2, scale_bits)?;

    // Encode
    let coeffs1 = encoding::encode(&values1, &encoding_params)?;
    let coeffs2 = encoding::encode(&values2, &encoding_params)?;
    let coeffs3 = encoding::encode(&values3, &encoding_params)?;

    // Convert to polynomial (you might need to add a helper function)
    let poly1 = PolyRing::from_signed_coeffs(&coeffs1, modulus, ring_degree as u64);
    let poly2 = PolyRing::from_signed_coeffs(&coeffs2, modulus, ring_degree as u64);
    let poly3 = PolyRing::from_signed_coeffs(&coeffs3, modulus, ring_degree as u64);

    // Encrypt
    let scale = (1u64 << scale_bits) as f64;
    let ct1 = encrypt(&poly1, &public_key, scale);
    let ct2 = encrypt(&poly2, &public_key, scale);
    let ct3 = encrypt(&poly3, &public_key, scale);

    // Homomorphic Addition and multiplication
    let ct_sum = ct1.add(&ct2).mul(&ct3, &public_key.relin_key);

    // Decrypt
    let decrypted_poly = decrypt(&ct_sum, &secret_key);

    // Convert back to coefficients
    let decrypted_coeffs = poly_to_coeffs(&decrypted_poly);

    // Decode
    let result = encoding::decode(&decrypted_coeffs, &encoding_params)?;

    // Print results
    println!("values_1: {:?}", values1);
    println!("values_2: {:?}", values2);
    println!("values_3: {:?}", values3);
    println!("Expected result of (values1 + values2) * values3: [56, 36, 2, 60]");
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
