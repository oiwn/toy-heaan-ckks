use crate::{PolyRing, RelinearizationKey};

/// Represents a ciphertext in the CKKS scheme
#[derive(Clone, Debug)]
pub struct Ciphertext {
    /// c0 component of ciphertext
    pub c0: PolyRing,
    /// c1 component of ciphertext
    pub c1: PolyRing,
    /// c2 component of ciphertext (coefficient of s^2),
    /// present only after multiplication before relinearization
    pub c2: Option<PolyRing>,
    /// Tracks the scaling factor used
    pub scale: f64,
}

impl Ciphertext {
    /// Homomorphic addition of ciphertexts
    pub fn add(&self, other: &Self) -> Self {
        // Check scaling factors and moduli
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale"
        );

        // Add corresponding polynomials component-wise
        let c0 = self.c0.clone() + other.c0.clone();
        let c1 = self.c1.clone() + other.c1.clone();

        Self {
            c0,
            c1,
            c2: None,
            scale: self.scale,
        }
    }

    pub fn mul(&self, other: &Self, relin_key: &RelinearizationKey) -> Self {
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale for multiplication"
        );

        // Same steps for creating d0, d1, d2...
        let d0 = self.c0.clone() * other.c0.clone();
        let d1 = (self.c0.clone() * other.c1.clone())
            + (self.c1.clone() * other.c0.clone());
        let d2 = self.c1.clone() * other.c1.clone();

        let intermediate = Self {
            c0: d0,
            c1: d1,
            c2: Some(d2),
            scale: self.scale * other.scale,
        };

        // Apply relinearization to eliminate the s^2 term
        let relinearized = intermediate.relinearize(relin_key);

        // Rescale back down to the original scale
        relinearized.rescale(self.scale)
    }

    pub fn relinearize(&self, relin_key: &RelinearizationKey) -> Self {
        if self.c2.is_none() {
            // No s^2 term, nothing to do
            return self.clone();
        }

        let c2 = self.c2.as_ref().unwrap();

        // Use the b and a components directly from the relin_key
        let b0 = &relin_key.b;
        let a0 = &relin_key.a;

        // Compute c2 * (b0 + a0*s)
        let c0_new = self.c0.clone() + (c2.clone() * b0.clone());
        let c1_new = self.c1.clone() + (c2.clone() * a0.clone());

        Self {
            c0: c0_new,
            c1: c1_new,
            c2: None,
            scale: self.scale,
        }
    }

    pub fn rescale(&self, target_scale: f64) -> Self {
        // Calculate how much we need to scale down by
        let scale_factor = self.scale / target_scale;
        println!(
            "Rescaling from {} to {} (factor: {})",
            self.scale, target_scale, scale_factor
        );

        // Rescale the polynomials by dividing coefficients
        let c0_scaled = rescale_poly(&self.c0, scale_factor);
        let c1_scaled = rescale_poly(&self.c1, scale_factor);

        Self {
            c0: c0_scaled,
            c1: c1_scaled,
            c2: None,
            scale: target_scale,
        }
    }
}

pub fn rescale_poly(poly: &PolyRing, scale_factor: f64) -> PolyRing {
    println!("rescale_poly: scale_factor = {}", scale_factor);
    let shift_bits = scale_factor.log2().round() as u32;
    println!("rescale_poly: shift_bits = {}", shift_bits);

    let modulus = poly.modulus();
    let half_modulus = modulus / 2;
    let ring_dim = poly.ring_dim();

    // Convert scale_factor to power-of-2 shift amount
    let shift_bits = scale_factor.log2().round() as u32;
    println!("Bit shifts: {}", shift_bits);

    let mut new_coeffs = Vec::with_capacity(poly.len());

    for &coeff in poly {
        // Convert to centered representation (-q/2, q/2)
        let centered_coeff = if coeff > half_modulus {
            (coeff as i128) - (modulus as i128)
        } else {
            coeff as i128
        };

        // Apply bit shifting with proper rounding
        let round_bit = 1i128 << (shift_bits - 1);
        let scaled_coeff = if centered_coeff >= 0 {
            (centered_coeff + round_bit) >> shift_bits
        } else {
            (centered_coeff - round_bit) >> shift_bits
        };

        // Convert back to modular representation
        let scaled = if scaled_coeff < 0 {
            (modulus as i128 + scaled_coeff) as u64
        } else {
            scaled_coeff as u64
        } % modulus;

        new_coeffs.push(scaled);
    }

    // Print first few coefficients before/after
    println!("Before rescaling: {:?}", &poly.coefficients[..4]);
    println!("After rescaling: {:?}", &new_coeffs[..4]);

    PolyRing::from_coeffs(&new_coeffs, modulus, ring_dim)
}

#[cfg(test)]
mod tests {
    use crate::{
        Ciphertext, PolyRing, PublicKey, PublicKeyParams, RelinearizationKey,
        RelinearizationKeyParams, SecretKey, SecretKeyParams, decrypt, encoding,
        encrypt, rescale_poly,
    };
    use rand::rng;

    #[test]
    fn test_rescale_function_simple() {
        let modulus = (1u64 << 61) - 1;
        let ring_degree = 8;

        // Create a simple polynomial with known coefficients
        let coeffs = vec![
            1000000000u64,
            2000000000u64,
            3000000000u64,
            4000000000u64,
            0,
            0,
            0,
            0,
        ];
        let poly = PolyRing::from_coeffs(&coeffs, modulus, ring_degree);

        println!("Original poly coeffs: {:?}", &coeffs[..4]);

        let scale_factor = 1000000.0; // Simple scale factor

        println!("About to call rescale_poly with factor: {}", scale_factor);

        // This is where it's likely hanging
        let rescaled = rescale_poly(&poly, scale_factor);

        println!("Rescaling completed!");
        println!("Rescaled coeffs: {:?}", &rescaled.coefficients[..4]);
    }

    #[test]
    fn test_multiplication_pipeline_debug() {
        let mut rng = rng();
        let ring_dim = 8;
        let scale_bits = 20;
        let modulus = (1u64 << 61) - 1;

        // Setup keys
        let secret_key_params = SecretKeyParams {
            ring_dim,
            modulus,
            hamming_weight: 4,
        };
        let secret_key = SecretKey::generate(&secret_key_params, &mut rng);

        let public_key_params = PublicKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.2,
        };
        let public_key =
            PublicKey::from_secret_key(&secret_key, &public_key_params, &mut rng);

        let relin_params: RelinearizationKeyParams = (&public_key_params).into();
        let relin_key = RelinearizationKey::from_secret_key(
            &secret_key,
            &relin_params,
            &mut rng,
        );

        // Create test vectors
        let values1 = vec![1.0, 2.0, 3.0, 4.0];
        let values2 = vec![5.0, 6.0, 7.0, 8.0];
        let expected: Vec<f64> =
            values1.iter().zip(&values2).map(|(&a, &b)| a * b).collect();

        println!("=== MULTIPLICATION PIPELINE DEBUG ===");
        println!("Values1: {:?}", values1);
        println!("Values2: {:?}", values2);
        println!("Expected product: {:?}", expected);

        // Encode and encrypt
        let enc_params =
            encoding::EncodingParams::new(ring_dim, scale_bits).unwrap();
        let coeffs1 = encoding::encode(&values1, &enc_params).unwrap();
        let coeffs2 = encoding::encode(&values2, &enc_params).unwrap();
        let poly1 = PolyRing::from_coeffs(&coeffs1, modulus, ring_dim);
        let poly2 = PolyRing::from_coeffs(&coeffs2, modulus, ring_dim);

        let base_scale = (1u64 << scale_bits) as f64;
        let ct1 = encrypt(&poly1, &public_key, base_scale, &mut rng);
        let ct2 = encrypt(&poly2, &public_key, base_scale, &mut rng);

        println!("\n--- Step 1: Initial ciphertexts ---");
        println!("ct1.scale: {}", ct1.scale);
        println!("ct2.scale: {}", ct2.scale);

        // Manual multiplication step by step
        println!("\n--- Step 2: Tensor multiplication ---");
        let d0 = ct1.c0.clone() * ct2.c0.clone();
        let d1 =
            (ct1.c0.clone() * ct2.c1.clone()) + (ct1.c1.clone() * ct2.c0.clone());
        let d2 = ct1.c1.clone() * ct2.c1.clone();

        let intermediate = Ciphertext {
            c0: d0,
            c1: d1,
            c2: Some(d2),
            scale: ct1.scale * ct2.scale,
        };
        println!(
            "Intermediate scale (should be 2^40): {}",
            intermediate.scale
        );

        // Test decryption at this point
        println!("\n--- Step 3: Decrypt BEFORE relinearization ---");
        // Create a temporary ciphertext without c2 for decryption test
        let temp_ct = Ciphertext {
            c0: intermediate.c0.clone(),
            c1: intermediate.c1.clone(),
            c2: None,
            scale: intermediate.scale,
        };

        let decrypted_before_relin = decrypt(&temp_ct, &secret_key);
        let coeffs_before_relin = poly_to_coeffs(&decrypted_before_relin);

        // We need to create encoding params for the high scale
        let high_scale_bits = 40; // since scale is 2^40
        let high_enc_params =
            encoding::EncodingParams::new(ring_dim, high_scale_bits).unwrap();
        let decoded_before_relin =
            encoding::decode(&coeffs_before_relin, &high_enc_params).unwrap();
        println!("Decoded BEFORE relinearization: {:?}", decoded_before_relin);

        // Relinearization
        println!("\n--- Step 4: Relinearization ---");
        let relinearized = intermediate.relinearize(&relin_key);
        println!("After relinearization scale: {}", relinearized.scale);

        // Test decryption after relinearization
        let decrypted_after_relin = decrypt(&relinearized, &secret_key);
        let coeffs_after_relin = poly_to_coeffs(&decrypted_after_relin);
        let decoded_after_relin =
            encoding::decode(&coeffs_after_relin, &high_enc_params).unwrap();
        println!("Decoded AFTER relinearization: {:?}", decoded_after_relin);

        // Rescaling
        println!("\n--- Step 5: Rescaling ---");
        let final_result = relinearized.rescale(base_scale);
        println!("Final scale: {}", final_result.scale);

        // Final decryption
        let final_decrypted = decrypt(&final_result, &secret_key);
        let final_coeffs = poly_to_coeffs(&final_decrypted);
        let final_decoded = encoding::decode(&final_coeffs, &enc_params).unwrap();
        println!("FINAL decoded: {:?}", final_decoded);
        println!("Expected: {:?}", expected);
    }

    #[test]
    fn test_rescale_poly_debug() {
        let ring_degree = 8;
        let scale_bits = 20; // Base scale 2^20
        let modulus = (1u64 << 61) - 1;
        let base_scale = (1u64 << scale_bits) as f64; // 2^20 = 1048576
        let high_scale = base_scale * base_scale; // 2^40 = 1099511627776

        println!("=== RESCALE DEBUG TEST ===");
        println!("Base scale (2^20): {}", base_scale);
        println!("High scale (2^40): {}", high_scale);
        println!("Scale factor (high/base): {}", high_scale / base_scale);

        // Step 1: Encode [5, 12, 21, 32] with base scale
        let values = vec![5.0, 12.0, 21.0, 32.0];
        let enc_params =
            encoding::EncodingParams::new(ring_degree, scale_bits).unwrap();
        let coeffs_base = encoding::encode(&values, &enc_params).unwrap();
        let _poly_base = PolyRing::from_coeffs(&coeffs_base, modulus, ring_degree);

        println!("\n--- Step 1: Base encoding ---");
        println!("Original values: {:?}", values);
        println!("Encoded coeffs (first 4): {:?}", &coeffs_base[..4]);

        // Step 2: Manually simulate what happens after multiplication
        // When we multiply two ciphertexts with scale 2^20, we get scale 2^40
        // The polynomial coefficients should be multiplied by the extra scale factor
        let scale_multiplier = (high_scale / base_scale) as u64; // Should be 2^20
        let mut coeffs_high = coeffs_base.clone();
        for coeff in &mut coeffs_high {
            *coeff = (*coeff * (scale_multiplier as i64)) % modulus as i64;
        }
        let poly_high = PolyRing::from_coeffs(&coeffs_high, modulus, ring_degree);

        println!("\n--- Step 2: High scale simulation ---");
        println!("Scale multiplier: {}", scale_multiplier);
        println!("High scale coeffs (first 4): {:?}", &coeffs_high[..4]);

        // Step 3: Use rescale_poly to scale back down
        println!("\n--- Step 3: Rescaling back down ---");
        println!(
            "About to call rescale_poly with factor: {}",
            scale_multiplier as f64
        );
        let poly_rescaled = rescale_poly(&poly_high, scale_multiplier as f64);

        // Step 4: Check if we got back the original coefficients
        println!("\n--- Step 4: Comparison ---");
        println!("Original coeffs (first 4): {:?}", &coeffs_base[..4]);
        println!(
            "Rescaled coeffs (first 4): {:?}",
            &poly_rescaled.coefficients[..4]
        );

        // Step 5: Decode the rescaled result
        let rescaled_coeffs_signed = poly_to_coeffs(&poly_rescaled);
        let decoded_values =
            encoding::decode(&rescaled_coeffs_signed, &enc_params).unwrap();

        println!("\n--- Step 5: Final decode ---");
        println!("Expected values: {:?}", values);
        println!("Decoded values: {:?}", decoded_values);

        // Verify the results
        for (expected, actual) in values.iter().zip(&decoded_values) {
            let error = (expected - actual).abs();
            println!("Expected: {}, Got: {}, Error: {}", expected, actual, error);
            assert!(
                error < 0.1,
                "Rescaling failed: expected {}, got {}",
                expected,
                actual
            );
        }
    }

    #[test]
    fn test_key_generation_isolated() {
        let mut rng = rng();
        let ring_dim = 8;
        let modulus = (1u64 << 61) - 1;

        println!("1a. Creating secret key params...");
        let secret_key_params = SecretKeyParams {
            ring_dim,
            modulus,
            hamming_weight: 4,
        };
        println!("Secret key params created!");

        println!("1b. Generating secret key...");
        let secret_key = SecretKey::generate(&secret_key_params, &mut rng);
        println!("Secret key generated!");

        println!("1c. Creating public key params...");
        let public_key_params = PublicKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.2,
        };
        println!("Public key params created!");

        println!("1d. Generating public key...");
        let _public_key =
            PublicKey::from_secret_key(&secret_key, &public_key_params, &mut rng);
        println!("Public key generated!");
    }

    #[test]
    fn test_rescaling_isolated() {
        let mut rng = rng();
        let ring_dim = 8;
        let scale_bits = 20;
        let modulus = (1u64 << 61) - 1;

        println!("1. Setting up keys...");
        let secret_key_params = SecretKeyParams {
            ring_dim,
            modulus,
            hamming_weight: 4,
        };
        let secret_key = SecretKey::generate(&secret_key_params, &mut rng);
        let public_key_params = PublicKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.2,
        };
        let public_key =
            PublicKey::from_secret_key(&secret_key, &public_key_params, &mut rng);
        println!("Keys generated!");

        let values = vec![5.0, 12.0, 21.0, 32.0];

        println!("2. Encoding...");
        let enc_params =
            encoding::EncodingParams::new(ring_dim, scale_bits).unwrap();
        let coeffs = encoding::encode(&values, &enc_params).unwrap();
        let poly = PolyRing::from_coeffs(&coeffs, modulus, ring_dim);
        println!("Encoding complete!");

        println!("3. Encrypting...");
        let base_scale = (1u64 << scale_bits) as f64;
        let ct = encrypt(&poly, &public_key, base_scale, &mut rng);
        println!("Encryption complete! Scale: {}", ct.scale);

        println!("4. Creating high scale ciphertext...");
        let high_scale_ct = Ciphertext {
            c0: ct.c0.clone(),
            c1: ct.c1.clone(),
            c2: None,
            scale: ct.scale * ct.scale,
        };
        println!("High scale: {}", high_scale_ct.scale);

        println!("5. Rescaling...");
        let rescaled_ct = high_scale_ct.rescale(base_scale);
        println!("Rescaling complete! New scale: {}", rescaled_ct.scale);

        println!("6. Decrypting...");
        let decrypted_poly = decrypt(&rescaled_ct, &secret_key);
        println!("Decryption complete!");

        println!("7. Converting to coeffs...");
        let decrypted_coeffs = poly_to_coeffs(&decrypted_poly);
        println!("Coeffs: {:?}", &decrypted_coeffs[..4]);

        println!("8. Decoding...");
        let decrypted_values =
            encoding::decode(&decrypted_coeffs, &enc_params).unwrap();
        println!("Decoding complete!");

        println!("Expected: {:?}", values);
        println!("Got: {:?}", decrypted_values);
    }

    /* #[test]
    fn test_rescaling_isolated() {
        let mut rng = rng();
        let ring_degree = 8;
        let scale_bits = 20; // 2^20 base scale
        let modulus = (1u64 << 61) - 1;

        // Create keys
        let secret_key_params = SecretKeyParams {
            ring_degree,
            modulus,
            hamming_weight: 100,
        };
        let secret_key = SecretKey::generate(&secret_key_params, &mut rng);
        let public_key_params = PublicKeyParams {
            poly_len: ring_degree,
            modulus,
            error_variance: 3.6,
        };
        let public_key =
            PublicKey::from_secret_key(&secret_key, &public_key_params, &mut rng);

        // Create test values
        let values = vec![5.0, 12.0, 21.0, 32.0]; // Expected multiplication result

        // Encode with base scale
        let enc_params =
            encoding::EncodingParams::new(ring_degree, scale_bits).unwrap();
        let coeffs = encoding::encode(&values, &enc_params).unwrap();
        let poly = PolyRing::from_coeffs(&coeffs, modulus, ring_degree);

        // Encrypt with base scale
        let base_scale = (1u64 << scale_bits) as f64;
        let ct = encrypt(&poly, &public_key, base_scale, &mut rng);

        println!("Original ciphertext scale: {}", ct.scale);

        // Simulate what happens after multiplication - create a ciphertext with doubled scale
        let high_scale_ct = Ciphertext {
            c0: ct.c0.clone(),
            c1: ct.c1.clone(),
            c2: None,
            scale: ct.scale * ct.scale, // This simulates post-multiplication scale (2^40)
        };

        println!("High scale (post-multiplication): {}", high_scale_ct.scale);

        // Now test rescaling back down
        let rescaled_ct = high_scale_ct.rescale(base_scale);

        println!("After rescaling: {}", rescaled_ct.scale);

        // Decrypt and check if we get back the original values
        let decrypted_poly = decrypt(&rescaled_ct, &secret_key);
        let decrypted_coeffs = poly_to_coeffs(&decrypted_poly);
        let decrypted_values =
            encoding::decode(&decrypted_coeffs, &enc_params).unwrap();

        println!("Expected values: {:?}", values);
        println!("Decrypted after rescale: {:?}", decrypted_values);

        // Check if rescaling preserved the values
        for (expected, actual) in values.iter().zip(&decrypted_values) {
            assert!(
                (expected - actual).abs() < 0.1,
                "Rescaling error: expected {}, got {}",
                expected,
                actual
            );
        }
    } */

    // Helper function you'll need if not already present
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
}
