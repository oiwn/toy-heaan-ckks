//! Secret Key (sk): Sample a "small" polynomial s(X) from R.
//! "Small" means its coefficients are small (e.g., chosen from {-1, 0, 1}).
use crate::PolyRing;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

/// A secret key in the CKKS scheme
pub struct SecretKey {
    /// The secret polynomial with ternary coefficients
    pub s: PolyRing,
}

/// Parameters specific to secret key generation
pub struct SecretKeyParams {
    /// Ring degree (must be power of 2)
    pub ring_degree: usize,
    /// Modulus for the polynomial coefficients
    pub modulus: u64,
    /// Number of non-zero coefficients in the secret key (Hamming weight)
    pub hamming_weight: usize,
}

impl SecretKey {
    /// Generate a new secret key with ternary coefficients {-1, 0, 1}
    pub fn generate<R: Rng>(params: &SecretKeyParams, rng: &mut R) -> Self {
        // Initialize coefficients to zero
        let mut coeffs = vec![0u64; params.ring_degree];
        let mut count = 0;

        // Add exactly hamming_weight non-zero coefficients
        while count < params.hamming_weight {
            let idx = rng.random_range(0..params.ring_degree);

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
        let s = PolyRing::from_coeffs(&coeffs, params.modulus, params.ring_degree);

        Self { s }
    }
}
