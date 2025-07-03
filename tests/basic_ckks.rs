/* #[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use toy_heaan_ckks::{
        PolyRing, PublicKey, PublicKeyParams, RelinearizationKey,
        RelinearizationKeyParams, SecretKey, SecretKeyParams, decrypt, encoding,
        encrypt,
    };

    #[test]
    fn test_basic_ckks_operations() {
        // Parameters setup - small for testing
        let ring_dim = 8;
        let scale_bits = 20;
        let modulus = (1u64 << 61) - 1; // Large prime-like modulus
        let mut rng = ChaCha20Rng::seed_from_u64(42); // Fixed seed for reproducibility

        // Create Secret Key
        let sk_params = SecretKeyParams {
            ring_dim,
            modulus,
            hamming_weight: 3,
        };
        let secret_key = SecretKey::generate(&sk_params, &mut rng);

        // Create Public Key
        let pk_params = PublicKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.0,
        };
        let public_key =
            PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng);

        // Create relinearization key
        let relin_params: RelinearizationKeyParams = (&pk_params).into();
        let relin_key = RelinearizationKey::from_secret_key(
            &secret_key,
            &relin_params,
            &mut rng,
        );

        // Test vectors
        let values1 = vec![1.0, 2.0, 3.0, 4.0];
        let values2 = vec![5.0, 6.0, 7.0, 8.0];

        // Expected results
        let values_add: Vec<f64> =
            values1.iter().zip(&values2).map(|(&a, &b)| a + b).collect();

        let values_mult: Vec<f64> =
            values1.iter().zip(&values2).map(|(&a, &b)| a * b).collect();

        // Parameters for encoding
        let encoding_params =
            encoding::EncodingParams::new(ring_dim, scale_bits).unwrap();

        // Encode to polynomials
        let coeffs1 = encoding::encode(&values1, &encoding_params).unwrap();
        let coeffs2 = encoding::encode(&values2, &encoding_params).unwrap();

        // Convert to polynomial
        let poly1 = PolyRing::from_coeffs(&coeffs1, modulus, ring_dim);
        let poly2 = PolyRing::from_coeffs(&coeffs2, modulus, ring_dim);

        println!("Original values1: {:?}", values1);
        println!("Original values2: {:?}", values2);

        // Encrypt
        let scale = (1u64 << scale_bits) as f64;
        let ct1 = encrypt(&poly1, &public_key, scale, &mut rng);
        let ct2 = encrypt(&poly2, &public_key, scale, &mut rng);

        // Homomorphic addition
        println!("--- Performing homomorphic addition ---");
        let ct_add = ct1.add(&ct2);

        // Homomorphic multiplication and rescaling
        println!("--- Performing homomorphic multiplication ---");
        let ct_mult = ct1.mul(&ct2, &relin_key);

        // Decrypt results
        let poly_add = decrypt(&ct_add, &secret_key);
        let poly_mult = decrypt(&ct_mult, &secret_key);

        // Convert to coefficients
        let decrypted_coeffs_add = poly_to_coeffs(&poly_add);
        let decrypted_coeffs_mult = poly_to_coeffs(&poly_mult);

        // Create params with the actual scales
        let add_params =
            encoding::EncodingParams::new(ring_dim, scale_bits).unwrap();
        let mult_params =
            encoding::EncodingParams::new(ring_dim, scale_bits).unwrap();

        // Decode
        let result_add =
            encoding::decode(&decrypted_coeffs_add, &add_params).unwrap();
        let result_mult =
            encoding::decode(&decrypted_coeffs_mult, &mult_params).unwrap();

        // Print and compare results
        println!("Expected addition: {:?}", values_add);
        println!("Decrypted addition: {:?}", result_add);

        println!("Expected multiplication: {:?}", values_mult);
        println!("Decrypted multiplication: {:?}", result_mult);

        // Verify results are close to expected (within some error margin)
        for (expected, actual) in values_add.iter().zip(&result_add) {
            assert!(
                (expected - actual).abs() < 0.1,
                "Addition error too large: expected {}, got {}",
                expected,
                actual
            );
        }

        for (expected, actual) in values_mult.iter().zip(&result_mult) {
            assert!(
                (expected - actual).abs() < 0.5,
                "Multiplication error too large: expected {}, got {}",
                expected,
                actual
            );
        }

        println!("!!! END TEST BASIC !!!");
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

    #[test]
    fn test_formula_sanity_check() {
        let modulus = 65537u64;
        let ring_dim = 8;

        // Simple test values
        let c0 = PolyRing::from_coeffs(&[1, 0, 0, 0], modulus, ring_dim);
        let c1 = PolyRing::from_coeffs(&[2, 0, 0, 0], modulus, ring_dim);
        let d0 = PolyRing::from_coeffs(&[3, 0, 0, 0], modulus, ring_dim);
        let d1 = PolyRing::from_coeffs(&[4, 0, 0, 0], modulus, ring_dim);

        // Expected result:
        // c0*d0 + c0*d1 + c1*d0 + c1*d1 =
        // 1*3 + 1*4 + 2*3 + 2*4 = 3 + 4 + 6 + 8 = 21

        // Using formula
        let c0_plus_c1 = c0.clone() + c1.clone();
        let d0_plus_d1 = d0.clone() + d1.clone();
        let sum_prod = c0_plus_c1.clone() * d0_plus_d1.clone();

        let c0d0 = c0.clone() * d0.clone();
        let c1d1 = c1.clone() * d1.clone();

        println!("c0+c1: {}", c0_plus_c1); // Should be 3
        println!("d0+d1: {}", d0_plus_d1); // Should be 7
        println!("sum_prod: {}", sum_prod); // Should be 21
        println!("c0d0: {}", c0d0); // Should be 3
        println!("c1d1: {}", c1d1); // Should be 8

        // Verify basic arithmetic
        assert_eq!(sum_prod.into_iter().next().unwrap(), &21u64);
    }
} */

fn main() {
    todo!();
}
