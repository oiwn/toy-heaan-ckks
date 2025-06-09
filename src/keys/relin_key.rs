use crate::{
    PolyRing, PublicKeyParams, SecretKey, generate_error_poly,
    generate_random_poly, negate_poly,
};
use rand::Rng;

/// Parameters for relinearization key generation.
pub struct RelinearizationKeyParams {
    /// Ring dimension (degree of the polynomial ring, n).
    pub ring_dim: usize,
    /// Modulus q.
    pub modulus: u64,
    /// Error variance σ².
    pub error_variance: f64,
}

/// Relinearization key used to transform ciphertexts after multiplication
pub struct RelinearizationKey {
    /// b component: -a*s + e + s^2
    pub b: PolyRing,
    /// a component: random polynomial
    pub a: PolyRing,
}

impl RelinearizationKey {
    /// Generate a relinearization key from a secret key and parameters.
    pub fn from_secret_key<R: Rng>(
        secret_key: &SecretKey,
        params: &RelinearizationKeyParams,
        rng: &mut R,
    ) -> Self {
        let ring_dim = secret_key.s.ring_dim();
        assert_eq!(
            ring_dim, params.ring_dim,
            "Ring dimension mismatch: secret_key.ring_dim()={} but params.ring_dim={}",
            ring_dim, params.ring_dim
        );

        let s = &secret_key.s;
        let s_squared = s.clone() * s.clone();

        let a =
            generate_random_poly(params.ring_dim, params.modulus, ring_dim, rng);
        let e = generate_error_poly(
            params.ring_dim,
            params.modulus,
            params.error_variance,
            ring_dim,
            rng,
        );

        let b =
            negate_poly(&(a.clone() * s.clone()), params.modulus) + e + s_squared;

        Self { a, b }
    }
}

impl From<&PublicKeyParams> for RelinearizationKeyParams {
    fn from(pk: &PublicKeyParams) -> Self {
        RelinearizationKeyParams {
            ring_dim: pk.ring_dim,
            modulus: pk.modulus,
            error_variance: pk.error_variance,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SecretKey;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_relinearization_key_properties() {
        let modulus = 1231u64;
        let ring_dim = 4;
        let params = RelinearizationKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.0,
        };

        // Secret key with known coeffs [1, 0, −1, 0]
        let mut s_coeffs = vec![0u64; ring_dim];
        s_coeffs[0] = 1;
        s_coeffs[2] = modulus - 1;
        let s = PolyRing::from_coeffs(&s_coeffs, modulus, ring_dim);
        let secret_key = SecretKey { s };

        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let relin_key =
            RelinearizationKey::from_secret_key(&secret_key, &params, &mut rng);

        // Verify b + a*s ≈ s² within expected error bounds
        let s_squared = secret_key.s.clone() * secret_key.s.clone();
        let a_times_s = relin_key.a.clone() * secret_key.s.clone();
        let verification = relin_key.b.clone() + a_times_s;

        let expected_error_bound = (3.0 * 3.0_f64.sqrt()).ceil() as i64; // 3σ bound = 6
        let half_modulus = modulus / 2;

        for (i, &coeff) in verification.coefficients.iter().enumerate() {
            let s2 = s_squared.coefficients[i];
            // 1) subtract mod q
            let raw_diff = if coeff >= s2 {
                coeff - s2
            } else {
                modulus - (s2 - coeff)
            };
            // 2) recenter into (−q/2..q/2)
            let signed_diff = if raw_diff > half_modulus {
                (raw_diff as i64) - (modulus as i64)
            } else {
                raw_diff as i64
            };

            assert!(
                signed_diff.abs() <= expected_error_bound,
                "Coefficient {} error too large: {} > {}",
                i,
                signed_diff,
                expected_error_bound
            );
        }
    }

    #[test]
    fn test_relinearization_key_to_delete() {
        // Create a secret key with known coefficients
        let modulus = 1231u64;
        let ring_dim = 4;
        let params = RelinearizationKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.0,
        };

        // Manually create a secret key with coefficients [1, 0, -1, 0]
        let mut s_coeffs = vec![0u64; ring_dim];
        s_coeffs[0] = 1;
        s_coeffs[2] = modulus - 1; // -1 mod q
        let s = PolyRing::from_coeffs(&s_coeffs, modulus, ring_dim);
        let secret_key = SecretKey { s };

        // We'll use a deterministic RNG that will produce specific values
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Create relinearization key
        let relin_key =
            RelinearizationKey::from_secret_key(&secret_key, &params, &mut rng);

        // Calculate s^2
        let s_squared = secret_key.s.clone() * secret_key.s.clone();

        // Verify key property: b + a*s ≈ s^2
        let a_times_s = relin_key.a.clone() * secret_key.s.clone();
        let verification = relin_key.b.clone() + a_times_s;

        // The difference should be the error term
        // We can't check exact equality due to random error,
        // but we can check that difference is small

        // Extract the difference into signed values
        let mut diff_values = Vec::new();
        let half_modulus = modulus / 2;

        for (&v1, &v2) in verification.into_iter().zip(s_squared.into_iter()) {
            let diff = if v1 >= v2 {
                v1 - v2
            } else {
                modulus - (v2 - v1)
            };

            // Convert to signed representation
            let signed_diff = if diff > half_modulus {
                (diff as i64) - (modulus as i64)
            } else {
                diff as i64
            };

            diff_values.push(signed_diff);
        }

        // Check that differences are within expected error bounds
        let expected_error_bound = (3.0 * 3.0_f64.sqrt()).ceil() as i64;
        for diff in diff_values {
            assert!(
                diff.abs() <= expected_error_bound,
                "Error too large: {}, expected at most: {}",
                diff,
                expected_error_bound
            );
        }
    }
}
