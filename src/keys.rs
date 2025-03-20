use crate::PolyRing;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

/// A secret key in the CKKS scheme
pub struct SecretKey {
    /// The secret polynomial s
    pub s: PolyRing,
}

/// A public key in the CKKS scheme
pub struct PublicKey {
    /// The first component of the public key (b = -a·s + e)
    pub b: PolyRing,
    /// The second component of the public key (random polynomial a)
    pub a: PolyRing,
}

/// Parameters for key generation
pub struct KeyGenParams {
    /// Polynomial degree (must be power of 2)
    pub n: usize,
    /// Modulus for the polynomial coefficients
    pub modulus: u64,
    /// Variance for error distribution
    pub error_variance: f64,
}

impl SecretKey {
    /// Generate a new secret key with binary coefficients
    pub fn generate(params: &KeyGenParams) -> Self {
        // Create RNG with proper seeding
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]); // In production need to change it to entropy gen

        // Generate binary polynomial for secret key
        let mut coeffs = Vec::with_capacity(params.n);
        for _ in 0..params.n {
            // Generate 0 or 1 coefficients - use random() instead of gen_bool
            let coeff: u64 = if rng.random::<bool>() { 1 } else { 0 };
            coeffs.push(coeff);
        }

        // Create the polynomial
        let s = PolyRing::from_coeffs(&coeffs, params.modulus);

        Self { s }
    }

    /// Generate the corresponding public key
    pub fn public_key(&self, params: &KeyGenParams) -> PublicKey {
        // Generate random polynomial a
        let mut rng = ChaCha20Rng::from_seed([1u8; 32]); // Different seed from secret key
        let a = generate_random_poly(params.n, params.modulus, &mut rng);

        // Generate small error polynomial e
        let e = generate_error_poly(
            params.n,
            params.modulus,
            params.error_variance,
            &mut rng,
        );

        // Compute b = -(a·s) + e
        let modulus = params.modulus;
        let a_times_s = a.clone() * self.s.clone();
        // Create negative of (a·s)
        let neg_a_s = negate_poly(&a_times_s, modulus);

        // b = -(a·s) + e
        let b = neg_a_s + e;

        PublicKey { a, b }
    }
}

/// Generate a polynomial with random coefficients in [0, modulus)
fn generate_random_poly<T: Rng>(n: usize, modulus: u64, rng: &mut T) -> PolyRing {
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Generate a random value modulo q - updated to use random_range
        let random_u64 = rng.random::<u64>();
        let random_val = random_u64 % modulus;
        coeffs.push(random_val);
    }

    PolyRing::from_coeffs(&coeffs, modulus)
}

/// Generate a polynomial with small random coefficients (error polynomial)
fn generate_error_poly<T: Rng>(
    n: usize,
    modulus: u64,
    variance: f64,
    rng: &mut T,
) -> PolyRing {
    // Simplified discrete Gaussian sampler
    let std_dev = variance.sqrt();
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Generate a simple approximation of Gaussian noise - using random() instead of gen
        let noise =
            (rng.random::<f64>() + rng.random::<f64>() + rng.random::<f64>() - 1.5)
                * std_dev;
        let noise_int = noise.round() as i64;

        // Convert to unsigned and reduce modulo q
        let coeff = if noise_int < 0 {
            let abs_noise = noise_int.unsigned_abs();
            modulus - (abs_noise % modulus)
        } else {
            noise_int as u64 % modulus
        };

        coeffs.push(coeff);
    }

    PolyRing::from_coeffs(&coeffs, modulus)
}

/// Compute the negative of a polynomial in the ring
fn negate_poly(poly: &PolyRing, modulus: u64) -> PolyRing {
    // Create a polynomial with coefficients (q - a_i) for each coefficient a_i
    let mut neg_coeffs = Vec::with_capacity(poly.degree() + 1);

    for coeff in poly {
        if *coeff == 0 {
            neg_coeffs.push(0);
        } else {
            neg_coeffs.push(modulus - *coeff);
        }
    }

    PolyRing::from_coeffs(&neg_coeffs, modulus)
}

/// A ciphertext in the CKKS scheme
pub struct Ciphertext {
    /// First component of the ciphertext
    pub c0: PolyRing,
    /// Second component of the ciphertext
    pub c1: PolyRing,
    /// Scaling factor to track precision
    pub scale: f64,
}

impl PublicKey {
    /// Encrypt a plaintext polynomial
    pub fn encrypt(
        &self,
        plaintext: &PolyRing,
        params: &KeyGenParams,
        scale: f64,
    ) -> Ciphertext {
        let mut rng = ChaCha20Rng::from_seed([2u8; 32]); // Different seed

        // Generate small random polynomial u
        let u = generate_random_poly(params.n, params.modulus, &mut rng);

        // Generate small error polynomials
        let e1 = generate_error_poly(
            params.n,
            params.modulus,
            params.error_variance,
            &mut rng,
        );
        let e2 = generate_error_poly(
            params.n,
            params.modulus,
            params.error_variance,
            &mut rng,
        );

        // Compute ciphertext components
        // c0 = b·u + e1 + m
        let b_times_u = self.b.clone() * u.clone();
        let c0 = b_times_u + e1 + plaintext.clone();

        // c1 = a·u + e2
        let a_times_u = self.a.clone() * u;
        let c1 = a_times_u + e2;

        Ciphertext { c0, c1, scale }
    }
}

impl SecretKey {
    /// Decrypt a ciphertext
    pub fn decrypt(&self, ciphertext: &Ciphertext) -> PolyRing {
        // m' = c0 + c1·s
        let c1_times_s = ciphertext.c1.clone() * self.s.clone();
        ciphertext.c0.clone() + c1_times_s
    }
}

#[cfg(test)]
mod encryption_tests {
    use super::*;
    use crate::encoding::{EncodingParams, decode, encode};

    #[test]
    fn test_key_generation() {
        let modulus = 1152921504606748673u64; // A large prime near 2^60
        let params = KeyGenParams {
            n: 8, // Small value for testing
            modulus,
            error_variance: 3.2,
        };

        // Generate secret key
        let sk = SecretKey::generate(&params);
        assert_eq!(sk.s.degree() + 1, params.n);

        // Generate public key
        let pk = sk.public_key(&params);
        assert_eq!(pk.a.degree() + 1, params.n);
        assert_eq!(pk.b.degree() + 1, params.n);
    }

    #[test]
    fn test_encrypt_decrypt() {
        // Setup parameters
        let modulus = 1152921504606748673u64;
        let params = KeyGenParams {
            n: 8,
            modulus,
            error_variance: 3.2,
        };

        // Generate keys
        let sk = SecretKey::generate(&params);
        let pk = sk.public_key(&params);

        // Create a test polynomial
        let plaintext_coeffs =
            (0..params.n).map(|i| (i as u64) % 10).collect::<Vec<_>>();
        let plaintext = PolyRing::from_coeffs(&plaintext_coeffs, modulus);

        // Encrypt
        let ciphertext = pk.encrypt(&plaintext, &params, 1.0);

        // Decrypt
        let decrypted = sk.decrypt(&ciphertext);

        // Compare
        for (original, recovered) in
            plaintext.into_iter().zip(decrypted.into_iter())
        {
            assert_eq!(original, recovered);
        }
    }

    #[test]
    fn test_end_to_end_with_encoding() {
        // Setup parameters
        let modulus = 1152921504606748673u64;
        let key_params = KeyGenParams {
            n: 8,
            modulus,
            error_variance: 3.2,
        };
        let encoding_params = EncodingParams::new(8, 30).unwrap();

        // Generate keys
        let sk = SecretKey::generate(&key_params);
        let pk = sk.public_key(&key_params);

        // Original data
        let original = vec![1.23, 4.56, -7.89, 0.12];

        // Encode
        let encoded_coeffs = encode(&original, &encoding_params).unwrap();

        // Convert to polynomial
        let encoded_poly_coeffs = encoded_coeffs
            .into_iter()
            .map(|x| {
                if x < 0 {
                    modulus - ((-x) as u64)
                } else {
                    x as u64
                }
            })
            .collect::<Vec<_>>();
        let plaintext = PolyRing::from_coeffs(&encoded_poly_coeffs, modulus);

        // Encrypt
        let scale = 2.0f64.powi(30); // Must match encoding_params.scale_bits
        let ciphertext = pk.encrypt(&plaintext, &key_params, scale);

        // Decrypt
        let decrypted_poly = sk.decrypt(&ciphertext);

        let decrypted_coeffs = decrypted_poly
            .into_iter()
            .map(|x| {
                let half_modulus = modulus / 2;
                if *x > half_modulus {
                    -((modulus - *x) as i64)
                } else {
                    *x as i64
                }
            })
            .collect::<Vec<_>>();

        // Decode
        let recovered = decode(&decrypted_coeffs, &encoding_params).unwrap();

        // Compare (with some tolerance due to encryption noise)
        for (orig, rec) in original.iter().zip(recovered.iter()) {
            assert!(
                (orig - rec).abs() < 0.1,
                "Values differ too much: {} vs {}",
                orig,
                rec
            );
        }
    }
}
