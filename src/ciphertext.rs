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
    let modulus = poly.modulus();
    let half_modulus = modulus / 2;
    let ring_dim = poly.ring_dim();

    // Convert scale_factor to power-of-2 shift amount
    let shift_bits = scale_factor.log2().round() as u32;

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

    PolyRing::from_coeffs(&new_coeffs, modulus, ring_dim)
}

#[cfg(test)]
mod tests {
    use crate::{
        Ciphertext, PolyRing, PublicKey, PublicKeyParams, SecretKey,
        SecretKeyParams, decrypt, encoding, encrypt, rescale_poly,
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
    fn test_key_generation_isolated() {
        let mut rng = rng();
        let ring_degree = 8;
        let modulus = (1u64 << 61) - 1;

        println!("1a. Creating secret key params...");
        let secret_key_params = SecretKeyParams {
            ring_degree,
            modulus,
            hamming_weight: 4,
        };
        println!("Secret key params created!");

        println!("1b. Generating secret key...");
        let secret_key = SecretKey::generate(&secret_key_params, &mut rng);
        println!("Secret key generated!");

        println!("1c. Creating public key params...");
        let public_key_params = PublicKeyParams {
            poly_len: ring_degree,
            modulus,
            error_variance: 3.2,
        };
        println!("Public key params created!");

        println!("1d. Generating public key...");
        let public_key =
            PublicKey::from_secret_key(&secret_key, &public_key_params, &mut rng);
        println!("Public key generated!");
    }

    #[test]
    fn test_rescaling_isolated() {
        let mut rng = rng();
        let ring_degree = 8;
        let scale_bits = 20;
        let modulus = (1u64 << 61) - 1;

        println!("1. Setting up keys...");
        let secret_key_params = SecretKeyParams {
            ring_degree,
            modulus,
            hamming_weight: 4,
        };
        let secret_key = SecretKey::generate(&secret_key_params, &mut rng);
        let public_key_params = PublicKeyParams {
            poly_len: ring_degree,
            modulus,
            error_variance: 3.2,
        };
        let public_key =
            PublicKey::from_secret_key(&secret_key, &public_key_params, &mut rng);
        println!("Keys generated!");

        let values = vec![5.0, 12.0, 21.0, 32.0];

        println!("2. Encoding...");
        let enc_params =
            encoding::EncodingParams::new(ring_degree, scale_bits).unwrap();
        let coeffs = encoding::encode(&values, &enc_params).unwrap();
        let poly = PolyRing::from_coeffs(&coeffs, modulus, ring_degree);
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
