//! Secret Key Generation for CKKS Scheme
//!
//! This module implements secret key generation for the CKKS homomorphic encryption scheme.
//! Secret keys are "small" polynomials where coefficients are sampled from the ternary
//! distribution {-1, 0, 1} with controlled sparsity (hamming weight).
//!
//! # Security Properties
//!
//! The secret key polynomial `s(X) ∈ R` satisfies:
//! - Coefficients are from `{-1, 0, 1}` (ternary distribution)
//! - Exactly `hamming_weight` coefficients are non-zero
//! - Non-zero positions are uniformly random
//! - +1/-1 values are uniformly distributed among non-zero coefficients
//!
//! # Example
//!
//! ```rust
//! use toy_heaan_ckks::{SecretKey, SecretKeyParams, NaivePolyRing};
//! use rand::rng;
//!
//! const DEGREE: usize = 1024;
//!
//! // Create parameters with half the coefficients non-zero
//! let params = SecretKeyParams::<DEGREE>::new(DEGREE).unwrap();
//! let modulus = 1125899906842679u64; // Example prime
//!
//! // Generate secret key
//! let mut rng = rng();
//! let secret_key: SecretKey<NaivePolyRing<DEGREE>, DEGREE> =
//!     SecretKey::generate(&params, &modulus, &mut rng).unwrap();
//! ```
use crate::{PolyRing, PolySampler};
use rand::Rng;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SecretKeyError {
    #[error("Hamming weight {0} exceeds ring dimension {1}")]
    InvalidHammingWeight(usize, usize),
    #[error("Invalid polynomial context: {0}")]
    InvalidContext(String),
}

/// Parameters for secret key generation
#[derive(Debug)]
pub struct SecretKeyParams<const DEGREE: usize> {
    /// Number of non-zero coefficients in the secret polynomial
    ///
    /// Must satisfy: 0 ≤ hamming_weight ≤ DEGREE
    ///
    /// # Security vs Performance Trade-off
    /// - Lower values: Potentially less secure but faster operations
    /// - Higher values: More secure but slower homomorphic operations
    /// - Sane choice: DEGREE/2 for balanced security/performance
    pub hamming_weight: usize,
}

impl<const DEGREE: usize> SecretKeyParams<DEGREE> {
    // Create new secret key parameters with validation
    pub fn new(hamming_weight: usize) -> Result<Self, SecretKeyError> {
        if hamming_weight > DEGREE {
            return Err(SecretKeyError::InvalidHammingWeight(
                hamming_weight,
                DEGREE,
            ));
        }
        Ok(Self { hamming_weight })
    }
}

/// Generic secret key for CKKS scheme
#[derive(Debug, Clone)]
pub struct SecretKey<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub poly: P,
}

impl<P, const DEGREE: usize> SecretKey<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    /// Generate a new secret key using ternary sampling
    ///
    /// Creates a polynomial with exactly `params.hamming_weight` non-zero coefficients,
    /// where each non-zero coefficient is ±1 with equal probability.
    ///
    /// # Arguments
    /// * `params` - Secret key parameters (validated automatically)
    /// * `context` - Polynomial ring context (modulus for naive, basis for RNS, etc.)
    /// * `rng` - Cryptographically secure random number generator
    ///
    /// # Returns
    /// * `Ok(SecretKey)` containing the generated polynomial
    /// * `Err(SecretKeyError)` if parameters are invalid
    ///
    /// # Security Requirements
    /// * Use a cryptographically secure RNG (e.g., `ChaCha20Rng`)
    /// * Keep the secret key confidential and secure
    /// * Use appropriate hamming weight for security level
    ///
    /// # Example
    /// ```rust
    /// # use toy_heaan_ckks::{SecretKey, SecretKeyParams, NaivePolyRing, PolyRing};
    /// # use rand::SeedableRng;
    /// # use rand_chacha::ChaCha20Rng;
    /// # const DEGREE: usize = 16;
    /// let params = SecretKeyParams::<DEGREE>::new(8)?;
    /// let modulus = 97u64;
    /// let mut rng = ChaCha20Rng::seed_from_u64(42);
    ///
    /// let secret_key: SecretKey<NaivePolyRing<DEGREE>, DEGREE> =
    ///     SecretKey::generate(&params, &modulus, &mut rng)?;
    ///
    /// // Verify hamming weight
    /// let coeffs = secret_key.poly.to_coeffs();
    /// let nonzero_count = coeffs.iter().filter(|&&c| c != 0).count();
    /// assert_eq!(nonzero_count, 8);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn generate<R: Rng>(
        params: &SecretKeyParams<DEGREE>,
        context: &P::Context,
        rng: &mut R,
    ) -> Result<Self, SecretKeyError> {
        // Use the polynomial backend's ternary sampling
        let poly = P::sample_tribits(params.hamming_weight, context, rng);

        Ok(SecretKey { poly })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::backends::naive::NaivePolyRing;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    const TEST_DEGREE: usize = 16;
    const TEST_MODULUS: u64 = 97; // Small prime for testing

    type TestSecretKey = SecretKey<NaivePolyRing<TEST_DEGREE>, TEST_DEGREE>;

    #[test]
    fn test_secret_key_params_creation() {
        // Valid hamming weight
        let params = SecretKeyParams::<TEST_DEGREE>::new(8).unwrap();
        assert_eq!(params.hamming_weight, 8);

        // Hamming weight equal to degree should be valid
        let params = SecretKeyParams::<TEST_DEGREE>::new(TEST_DEGREE).unwrap();
        assert_eq!(params.hamming_weight, TEST_DEGREE);
    }

    #[test]
    fn test_secret_key_params_invalid_hamming_weight() {
        // Hamming weight exceeding degree should fail
        let result = SecretKeyParams::<TEST_DEGREE>::new(TEST_DEGREE + 1);
        assert!(result.is_err());

        match result.unwrap_err() {
            SecretKeyError::InvalidHammingWeight(hw, deg) => {
                assert_eq!(hw, TEST_DEGREE + 1);
                assert_eq!(deg, TEST_DEGREE);
            }
            _ => panic!("Expected InvalidHammingWeight error"),
        }
    }

    #[test]
    fn test_secret_key_params_validation() {
        let _ = SecretKeyParams::<TEST_DEGREE>::new(8).is_ok();

        // Test with manually created invalid params (bypassing constructor)
        let invalid_params: Result<SecretKeyParams<TEST_DEGREE>, SecretKeyError> =
            SecretKeyParams::<TEST_DEGREE>::new(TEST_DEGREE + 5);
        assert!(invalid_params.is_err());
    }

    #[test]
    fn test_secret_key_generation() {
        let params = SecretKeyParams::<TEST_DEGREE>::new(8).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let secret_key = TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng);
        assert!(secret_key.is_ok());
    }

    #[test]
    fn test_secret_key_hamming_weight() {
        let hamming_weight = 6;
        let params = SecretKeyParams::<TEST_DEGREE>::new(hamming_weight).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(123);

        let secret_key =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng).unwrap();

        // Check that the polynomial has exactly the requested hamming weight
        let coeffs = secret_key.poly.to_coeffs();
        let actual_hamming_weight = coeffs.iter().filter(|&&c| c != 0).count();

        assert_eq!(
            actual_hamming_weight, hamming_weight,
            "Expected hamming weight {}, got {}",
            hamming_weight, actual_hamming_weight
        );
    }

    #[test]
    fn test_secret_key_ternary_coefficients() {
        let params = SecretKeyParams::<TEST_DEGREE>::new(8).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(789);

        let secret_key =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng).unwrap();
        let coeffs = secret_key.poly.to_coeffs();

        // All coefficients should be in {-1, 0, 1}
        for &coeff in &coeffs {
            assert!(
                coeff == -1 || coeff == 0 || coeff == 1,
                "Coefficient {} is not ternary",
                coeff
            );
        }
    }

    #[test]
    fn test_secret_key_zero_hamming_weight() {
        let params = SecretKeyParams::<TEST_DEGREE>::new(0).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(999);

        let secret_key =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng).unwrap();
        let coeffs = secret_key.poly.to_coeffs();

        // All coefficients should be zero
        for &coeff in &coeffs {
            assert_eq!(coeff, 0, "Expected zero coefficient, got {}", coeff);
        }
    }

    #[test]
    fn test_secret_key_full_hamming_weight() {
        let params = SecretKeyParams::<TEST_DEGREE>::new(TEST_DEGREE).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(555);

        let secret_key =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng).unwrap();
        let coeffs = secret_key.poly.to_coeffs();

        // All coefficients should be non-zero (±1)
        for &coeff in &coeffs {
            assert!(
                coeff == -1 || coeff == 1,
                "Expected ±1 coefficient, got {}",
                coeff
            );
        }
    }

    #[test]
    fn test_secret_key_deterministic_generation() {
        let params = SecretKeyParams::<TEST_DEGREE>::new(6).unwrap();

        // Same seed should produce same secret key
        let mut rng1 = ChaCha20Rng::seed_from_u64(12345);
        let mut rng2 = ChaCha20Rng::seed_from_u64(12345);

        let sk1 =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng1).unwrap();
        let sk2 =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng2).unwrap();

        let coeffs1 = sk1.poly.to_coeffs();
        let coeffs2 = sk2.poly.to_coeffs();

        assert_eq!(
            coeffs1, coeffs2,
            "Same seed should produce identical secret keys"
        );
    }

    #[test]
    fn test_secret_key_different_seeds() {
        let params = SecretKeyParams::<TEST_DEGREE>::new(8).unwrap();

        // Different seeds should produce different secret keys (with high probability)
        let mut rng1 = ChaCha20Rng::seed_from_u64(11111);
        let mut rng2 = ChaCha20Rng::seed_from_u64(22222);

        let sk1 =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng1).unwrap();
        let sk2 =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng2).unwrap();

        let coeffs1 = sk1.poly.to_coeffs();
        let coeffs2 = sk2.poly.to_coeffs();

        assert_ne!(
            coeffs1, coeffs2,
            "Different seeds should produce different secret keys"
        );
    }

    #[test]
    fn test_secret_key_coefficient_distribution() {
        let hamming_weight = 10;
        let params = SecretKeyParams::<TEST_DEGREE>::new(hamming_weight).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(777);

        // Generate multiple secret keys to test distribution
        let num_samples = 100;
        let mut plus_one_count = 0;
        let mut minus_one_count = 0;
        let mut zero_count = 0;

        for _ in 0..num_samples {
            let secret_key =
                TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng).unwrap();
            let coeffs = secret_key.poly.to_coeffs();

            for &coeff in &coeffs {
                match coeff {
                    1 => plus_one_count += 1,
                    -1 => minus_one_count += 1,
                    0 => zero_count += 1,
                    _ => panic!("Invalid coefficient: {}", coeff),
                }
            }
        }

        let total_coeffs = num_samples * TEST_DEGREE;
        let expected_nonzero = num_samples * hamming_weight;
        let expected_zero = total_coeffs - expected_nonzero;

        // Check that we have the expected number of zeros
        assert_eq!(zero_count, expected_zero);

        // Check that ±1 are roughly balanced (allow some variance)
        let nonzero_total = plus_one_count + minus_one_count;
        assert_eq!(nonzero_total, expected_nonzero);

        // Both should be > 0 and within reasonable bounds
        let min_expected = expected_nonzero / 3; // Allow significant variance for small sample
        assert!(plus_one_count > min_expected, "Too few +1 coefficients");
        assert!(minus_one_count > min_expected, "Too few -1 coefficients");
    }

    #[test]
    fn test_error_display() {
        let error = SecretKeyError::InvalidHammingWeight(100, 64);
        let display = format!("{}", error);
        assert!(display.contains("100"));
        assert!(display.contains("64"));
        assert!(display.contains("exceeds"));

        let error = SecretKeyError::InvalidContext("test error".to_string());
        let display = format!("{}", error);
        assert!(display.contains("test error"));
    }

    #[test]
    fn test_secret_key_clone() {
        let params = SecretKeyParams::<TEST_DEGREE>::new(4).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(333);

        let original =
            TestSecretKey::generate(&params, &TEST_MODULUS, &mut rng).unwrap();
        let cloned = original.clone();

        let orig_coeffs = original.poly.to_coeffs();
        let clone_coeffs = cloned.poly.to_coeffs();

        assert_eq!(orig_coeffs, clone_coeffs);
    }

    // Integration test with different ring degrees
    #[test]
    fn test_different_ring_degrees() {
        // Test with power-of-2 degrees
        test_degree_helper::<8>();
        test_degree_helper::<32>();
        test_degree_helper::<128>();
    }

    fn test_degree_helper<const DEGREE: usize>() {
        let hamming_weight = DEGREE / 4;
        let params = SecretKeyParams::<DEGREE>::new(hamming_weight).unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let secret_key: SecretKey<NaivePolyRing<DEGREE>, DEGREE> =
            SecretKey::generate(&params, &TEST_MODULUS, &mut rng).unwrap();

        let coeffs = secret_key.poly.to_coeffs();
        let actual_hw = coeffs.iter().filter(|&&c| c != 0).count();

        assert_eq!(actual_hw, hamming_weight);

        for &coeff in &coeffs {
            assert!(coeff >= -1 && coeff <= 1);
        }
    }
}
