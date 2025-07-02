//! Secret Key (sk): Sample a "small" polynomial s(X) from R.
//! "Small" means its coefficients are small (e.g., chosen from {-1, 0, 1})
use crate::{RnsBasis, RnsPolyRing, sample_ternary};
use rand::Rng;
use std::sync::Arc;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SecretKeyError {
    #[error("Hamming weight {0} exceeds ring dimension {1}")]
    InvalidHammingWeight(usize, usize),
}

/// Parameters for generating an RNS-encoded secret key.
pub struct SecretKeyParams<const DEGREE: usize> {
    pub basis: Arc<RnsBasis>,
    pub hamming_weight: usize,
}

impl<const DEGREE: usize> SecretKeyParams<DEGREE> {
    fn validate(&self) -> Result<(), SecretKeyError> {
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

/// RNS-encoded secret key wrapper.
pub struct SecretKey<const DEGREE: usize> {
    pub s: RnsPolyRing<DEGREE>,
}

impl<const DEGREE: usize> SecretKey<DEGREE> {
    /// Generate a new sparse ternary secret key in the RNS domain.
    pub fn generate<R: Rng + ?Sized>(
        params: &SecretKeyParams<DEGREE>,
        rng: &mut R,
    ) -> Result<Self, SecretKeyError> {
        params.validate()?;
        let plain = sample_ternary::<DEGREE, _>(params.hamming_weight, rng);
        let s = RnsPolyRing::from_integer_coeffs(&plain, params.basis.clone());
        Ok(SecretKey { s })
    }
}

/* //! Secret Key (sk): Sample a "small" polynomial s(X) from R.
//! "Small" means its coefficients are small (e.g., chosen from {-1, 0, 1}).
use crate::PolyRing;
use rand::Rng;

/// A secret key in the CKKS scheme
pub struct SecretKey {
    /// The secret polynomial with ternary coefficients
    pub s: PolyRing,
}

/// Parameters specific to secret key generation
pub struct SecretKeyParams {
    /// Ring degree (must be power of 2)
    pub ring_dim: usize,
    /// Modulus for the polynomial coefficients
    pub modulus: u64,
    /// Number of non-zero coefficients in the secret key (Hamming weight)
    pub hamming_weight: usize,
}

impl SecretKey {
    /// Generate a new secret key with ternary coefficients {-1, 0, 1}
    pub fn generate<R: Rng>(params: &SecretKeyParams, rng: &mut R) -> Self {
        // Validate parameters before starting
        if params.hamming_weight > params.ring_dim {
            panic!(
                "Invalid parameters: hamming_weight ({}) cannot be larger than ring_degree ({})",
                params.hamming_weight, params.ring_dim
            );
        }
        // Initialize coefficients to zero
        let mut coeffs = vec![0u64; params.ring_dim];
        let mut count = 0;

        // Add exactly hamming_weight non-zero coefficients
        while count < params.hamming_weight {
            let idx = rng.random_range(0..params.ring_dim);

            // Only set if position is still zero
            if coeffs[idx] == 0 {
                // Generate either 1 or -1 (represented as modulus-1) with equal probability
                coeffs[idx] = if rng.random_bool(0.5) {
                    1
                } else {
                    params.modulus - 1 // -1 mod q
                };
                count += 1;
            }
        }

        // Create the polynomial
        let s = PolyRing::from_coeffs(&coeffs, params.modulus, params.ring_dim);

        Self { s }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_secret_key_hamming_weight() {
        let params = SecretKeyParams {
            ring_dim: 128,
            modulus: 65537,
            hamming_weight: 40,
        };
        let mut rng = StdRng::seed_from_u64(0);
        let sk = SecretKey::generate(&params, &mut rng);
        let non_zero_count = sk.s.into_iter().filter(|&&x| x != 0).count();

        assert_eq!(non_zero_count, params.hamming_weight);
    }

    #[test]
    fn test_coefficient_values() {
        let params = SecretKeyParams {
            ring_dim: 64,
            modulus: 65537,
            hamming_weight: 20,
        };
        let mut rng = StdRng::seed_from_u64(0);
        let sk = SecretKey::generate(&params, &mut rng);

        for &coeff in &sk.s {
            assert!(
                coeff == 0 || coeff == 1 || coeff == params.modulus - 1,
                "Coefficient value should be 0, 1, or -1 (mod q)"
            );
        }
    }

    #[test]
    fn test_value_distribution() {
        // Generate multiple keys and check distribution of 1/-1 values
        let params = SecretKeyParams {
            ring_dim: 1024,
            modulus: 65537,
            hamming_weight: 500,
        };

        let mut ones_count = 0;
        let mut neg_ones_count = 0;

        for seed in 0..10 {
            // Generate 10 keys
            let mut rng = StdRng::seed_from_u64(seed);
            let sk = SecretKey::generate(&params, &mut rng);

            for &coeff in &sk.s {
                if coeff == 1 {
                    ones_count += 1;
                } else if coeff == params.modulus - 1 {
                    neg_ones_count += 1;
                }
            }
        }

        // We expect roughly equal numbers of 1 and -1
        let ratio = ones_count as f64 / (ones_count + neg_ones_count) as f64;
        assert!(
            (ratio - 0.5).abs() < 0.1,
            "1 and -1 should be equally distributed"
        );
    }
} */
