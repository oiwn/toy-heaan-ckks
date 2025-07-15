//! Secret Key (sk): Sample a "small" polynomial s(X) from R.
//! "Small" means its coefficients are small (e.g., chosen from {-1, 0, 1})
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
pub struct SecretKeyParams<const DEGREE: usize> {
    /// Number of non-zero coefficients in the secret polynomial
    pub hamming_weight: usize,
}

impl<const DEGREE: usize> SecretKeyParams<DEGREE> {
    /// Create new secret key parameters
    pub fn new(hamming_weight: usize) -> Result<Self, SecretKeyError> {
        if hamming_weight > DEGREE {
            return Err(SecretKeyError::InvalidHammingWeight(
                hamming_weight,
                DEGREE,
            ));
        }
        Ok(Self { hamming_weight })
    }

    /// Validate parameters
    pub fn validate(&self) -> Result<(), SecretKeyError> {
        if self.hamming_weight > DEGREE {
            Err(SecretKeyError::InvalidHammingWeight(
                self.hamming_weight,
                DEGREE,
            ))
        } else {
            Ok(())
        }
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
    pub fn generate<R: Rng>(
        params: &SecretKeyParams<DEGREE>,
        context: &P::Context,
        rng: &mut R,
    ) -> Result<Self, SecretKeyError> {
        // Validate parameters first
        params.validate()?;

        // Use the polynomial backend's ternary sampling
        let poly = P::sample_tribits(rng, params.hamming_weight, context);

        Ok(SecretKey { poly })
    }
}

/* #[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::RnsBasisBuilder;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use std::sync::Arc;

    /// Helper to create a test RNS basis with small primes
    fn create_test_basis<const DEGREE: usize>() -> Arc<RnsBasis> {
        Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19, 23]) // Small test primes
                .build()
                .unwrap(),
        )
    }

    #[test]
    fn test_secret_key_hamming_weight() {
        const DEGREE: usize = 128;
        let basis = create_test_basis::<DEGREE>();
        let params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 40,
        };

        let mut rng = StdRng::seed_from_u64(0);
        let sk = SecretKey::generate(&params, &mut rng).unwrap();

        // Count non-zero coefficients by reconstructing to integers
        let mut non_zero_count = 0;
        for i in 0..DEGREE {
            let coeff = sk.s.coefficient_to_u64(i);
            if coeff != 0 {
                non_zero_count += 1;
            }
        }

        assert_eq!(non_zero_count, params.hamming_weight);
    }

    #[test]
    fn test_coefficient_values() {
        const DEGREE: usize = 64;
        let basis = create_test_basis::<DEGREE>();
        let params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 20,
        };

        let mut rng = StdRng::seed_from_u64(0);
        let sk = SecretKey::generate(&params, &mut rng).unwrap();

        // Calculate modulus product for the test basis
        let modulus_product: u64 = basis.primes().iter().product();

        for i in 0..DEGREE {
            let coeff = sk.s.coefficient_to_u64(i);
            // In RNS, coefficients should be 0, 1, or modulus_product-1 (representing -1)
            assert!(
                coeff == 0 || coeff == 1 || coeff == modulus_product - 1,
                "Coefficient {} at position {} should be 0, 1, or {} (representing -1), got {}",
                i,
                i,
                modulus_product - 1,
                coeff
            );
        }
    }

    #[test]
    fn test_value_distribution() {
        const DEGREE: usize = 1024;
        let basis = create_test_basis::<DEGREE>();
        let params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 500,
        };

        let modulus_product: u64 = basis.primes().iter().product();
        let mut ones_count = 0;
        let mut neg_ones_count = 0;

        for seed in 0..10 {
            // Generate 10 keys
            let mut rng = StdRng::seed_from_u64(seed);
            let sk = SecretKey::generate(&params, &mut rng).unwrap();

            for i in 0..DEGREE {
                let coeff = sk.s.coefficient_to_u64(i);
                if coeff == 1 {
                    ones_count += 1;
                } else if coeff == modulus_product - 1 {
                    neg_ones_count += 1;
                }
            }
        }

        // We expect roughly equal numbers of 1 and -1
        let ratio = ones_count as f64 / (ones_count + neg_ones_count) as f64;
        assert!(
            (ratio - 0.5).abs() < 0.1,
            "1 and -1 should be equally distributed, got ratio: {}, ones: {}, neg_ones: {}",
            ratio,
            ones_count,
            neg_ones_count
        );
    }

    #[test]
    fn test_secret_key_reproducible() {
        const DEGREE: usize = 32;
        let basis = create_test_basis::<DEGREE>();
        let params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 16,
        };

        // Generate two keys with same seed
        let mut rng1 = StdRng::seed_from_u64(42);
        let mut rng2 = StdRng::seed_from_u64(42);

        let sk1 = SecretKey::generate(&params, &mut rng1).unwrap();
        let sk2 = SecretKey::generate(&params, &mut rng2).unwrap();

        // Should be identical
        for i in 0..DEGREE {
            assert_eq!(
                sk1.s.coefficient_to_u64(i),
                sk2.s.coefficient_to_u64(i),
                "Coefficients should be identical with same seed"
            );
        }
    }

    #[test]
    fn test_secret_key_params_validation() {
        const DEGREE: usize = 8;
        let basis = create_test_basis::<DEGREE>();

        // Test valid params
        let valid_params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 4,
        };
        assert!(valid_params.validate().is_ok());

        // Test invalid params (hamming weight > degree)
        let invalid_params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 10, // > DEGREE
        };

        let result = invalid_params.validate();
        assert!(result.is_err());

        if let Err(SecretKeyError::InvalidHammingWeight(weight, degree)) = result {
            assert_eq!(weight, 10);
            assert_eq!(degree, DEGREE);
        } else {
            panic!("Expected InvalidHammingWeight error");
        }
    }

    #[test]
    fn test_secret_key_generation_fails_invalid_params() {
        const DEGREE: usize = 8;
        let basis = create_test_basis::<DEGREE>();
        let params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 10, // > DEGREE
        };

        let mut rng = StdRng::seed_from_u64(0);
        let result = SecretKey::generate(&params, &mut rng);

        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SecretKeyError::InvalidHammingWeight(10, 8)
        ));
    }

    #[test]
    fn test_secret_key_zero_hamming_weight() {
        const DEGREE: usize = 16;
        let basis = create_test_basis::<DEGREE>();
        let params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 0,
        };

        let mut rng = StdRng::seed_from_u64(0);
        let sk = SecretKey::generate(&params, &mut rng).unwrap();

        // All coefficients should be zero
        for i in 0..DEGREE {
            assert_eq!(
                sk.s.coefficient_to_u64(i),
                0,
                "All coefficients should be zero"
            );
        }
    }

    #[test]
    fn test_secret_key_full_hamming_weight() {
        const DEGREE: usize = 8;
        let basis = create_test_basis::<DEGREE>();
        let params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: DEGREE, // All coefficients non-zero
        };

        let mut rng = StdRng::seed_from_u64(0);
        let sk = SecretKey::generate(&params, &mut rng).unwrap();

        let modulus_product: u64 = basis.primes().iter().product();
        let mut zero_count = 0;

        for i in 0..DEGREE {
            let coeff = sk.s.coefficient_to_u64(i);
            if coeff == 0 {
                zero_count += 1;
            } else {
                // Should be ±1
                assert!(
                    coeff == 1 || coeff == modulus_product - 1,
                    "Non-zero coefficient should be ±1, got {}",
                    coeff
                );
            }
        }

        assert_eq!(
            zero_count, 0,
            "No coefficients should be zero with full hamming weight"
        );
    }

    #[test]
    fn test_secret_key_with_different_bases() {
        const DEGREE: usize = 16;

        // Create two different bases
        let basis1 = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19])
                .build()
                .unwrap(),
        );

        let basis2 = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![23, 29, 31])
                .build()
                .unwrap(),
        );

        let params1: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis1.clone(),
            hamming_weight: 8,
        };

        let params2: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis2.clone(),
            hamming_weight: 8,
        };

        let mut rng = StdRng::seed_from_u64(123);

        let sk1 = SecretKey::generate(&params1, &mut rng).unwrap();
        let sk2 = SecretKey::generate(&params2, &mut rng).unwrap();

        // Verify they have different bases
        assert_ne!(sk1.s.basis.primes(), sk2.s.basis.primes());

        // But same hamming weight
        let count1 = (0..DEGREE)
            .filter(|&i| sk1.s.coefficient_to_u64(i) != 0)
            .count();
        let count2 = (0..DEGREE)
            .filter(|&i| sk2.s.coefficient_to_u64(i) != 0)
            .count();

        assert_eq!(count1, 8);
        assert_eq!(count2, 8);
    }
} */
