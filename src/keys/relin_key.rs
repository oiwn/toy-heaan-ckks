use crate::{
    PolyRing, SecretKey, generate_error_poly, generate_random_poly, negate_poly,
};
use rand::Rng;

/// Relinearization key used to transform ciphertexts after multiplication
pub struct RelinearizationKey {
    /// b component: -a*s + e + s^2
    pub b: PolyRing,
    /// a component: random polynomial
    pub a: PolyRing,
}

impl RelinearizationKey {
    pub fn from_secret_key<R: Rng>(
        secret_key: &SecretKey,
        modulus: u64,
        error_variance: f64,
        rng: &mut R,
    ) -> Self {
        let n = secret_key.s.len();
        let ring_dim = secret_key.s.ring_dim();

        // Compute s^2
        let s_squared = secret_key.s.clone() * secret_key.s.clone();

        // Generate random polynomial a
        let a = generate_random_poly(n, modulus, ring_dim, rng);

        // Generate error polynomial e
        let e = generate_error_poly(n, modulus, error_variance, ring_dim, rng);

        // Compute b = -(a*s) + e + s^2
        let a_times_s = a.clone() * secret_key.s.clone();
        let neg_a_s = negate_poly(&a_times_s, modulus);
        let b = neg_a_s + e + s_squared;

        Self { a, b }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SecretKey;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_relinearization_key() {
        // Create a secret key with known coefficients
        let modulus = 1231u64;
        let ring_dim = 4;

        // Manually create a secret key with coefficients [1, 0, -1, 0]
        let mut s_coeffs = vec![0u64; ring_dim];
        s_coeffs[0] = 1;
        s_coeffs[2] = modulus - 1; // -1 mod q
        let s = PolyRing::from_coeffs(&s_coeffs, modulus, ring_dim);
        let secret_key = SecretKey { s };

        // We'll use a deterministic RNG that will produce specific values
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Create relinearization key
        let relin_key = RelinearizationKey::from_secret_key(
            &secret_key,
            modulus,
            3.0, // error variance
            &mut rng,
        );

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
