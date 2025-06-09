use crate::{
    PolyRing, SecretKey, generate_error_poly, generate_random_poly, negate_poly,
};
use rand::Rng;

/// Parameters for public key generation
pub struct PublicKeyParams {
    /// Number of coefficients in polynomial
    pub ring_dim: usize,
    /// Modulus for polynomial coefficients
    pub modulus: u64,
    /// Variance for error distribution
    pub error_variance: f64,
}

/// A public key in the CKKS scheme
pub struct PublicKey {
    /// The first component of the public key (b = -a*s + e)
    pub b: PolyRing,
    /// The second component of the public key (random polynomial a)
    pub a: PolyRing,
}

impl PublicKey {
    /// Generate a public key from a secret key using the provided parameters
    pub fn from_secret_key<R: Rng>(
        secret_key: &SecretKey,
        params: &PublicKeyParams,
        rng: &mut R,
    ) -> Self {
        // Get the ring dimension from the secret key
        let ring_dim = secret_key.s.ring_dim();

        // Generate random polynomial a
        let a =
            generate_random_poly(params.ring_dim, params.modulus, ring_dim, rng);

        // Generate small error polynomial e
        let e = generate_error_poly(
            params.ring_dim,
            params.modulus,
            params.error_variance,
            ring_dim,
            rng,
        );

        // Compute b = -(a * s) + e
        let a_times_s = a.clone() * secret_key.s.clone();
        let neg_a_s = negate_poly(&a_times_s, params.modulus);
        let b = neg_a_s + e;

        PublicKey { a, b }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SecretKey, SecretKeyParams};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_public_key_generation() {
        // Setup parameters
        let ring_dim = 8;
        let modulus = 65537u64;

        // Create secret key
        let sk_params = SecretKeyParams {
            ring_dim,
            modulus,
            hamming_weight: 4,
        };
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let secret_key = SecretKey::generate(&sk_params, &mut rng);

        // Create public key
        let pk_params = PublicKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.0,
        };
        let mut rng = ChaCha20Rng::seed_from_u64(43);
        let public_key =
            PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng);

        // Verify dimensions and modulus
        assert_eq!(public_key.a.len(), ring_dim);
        assert_eq!(public_key.b.len(), ring_dim);
        assert_eq!(public_key.a.modulus(), modulus);
    }

    #[test]
    fn test_public_key_relation() {
        // Setup parameters
        let ring_dim = 8;
        let modulus = 65537u64;

        // Create secret key and public key
        let sk_params = SecretKeyParams {
            ring_dim,
            modulus,
            hamming_weight: 4,
        };
        let pk_params = PublicKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.0,
        };

        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let secret_key = SecretKey::generate(&sk_params, &mut rng);
        let public_key =
            PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng);

        // Extract error term: e = b + (a * s)
        let a_times_s = public_key.a.clone() * secret_key.s.clone();
        let extracted_e = public_key.b.clone() + a_times_s;

        // Verify error has small coefficients
        let error_bound = (3.0 * (3.0f64.sqrt() * 3.0)).round() as u64;
        let half_modulus = modulus / 2;

        for &coeff in &extracted_e {
            let signed_coeff = if coeff > half_modulus {
                (coeff as i64) - (modulus as i64)
            } else {
                coeff as i64
            };

            assert!(
                signed_coeff.unsigned_abs() <= error_bound,
                "Error coefficient too large: {}, bound: {}",
                signed_coeff,
                error_bound
            );
        }
    }

    #[test]
    fn test_deterministic_generation() {
        // Setup
        let ring_dim = 8;
        let modulus = 65537u64;

        // Create secret key
        let sk_params = SecretKeyParams {
            ring_dim,
            modulus,
            hamming_weight: 4,
        };
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let secret_key = SecretKey::generate(&sk_params, &mut rng);

        // Create public keys with same seed
        let pk_params = PublicKeyParams {
            ring_dim,
            modulus,
            error_variance: 3.0,
        };

        let mut rng1 = ChaCha20Rng::seed_from_u64(123);
        let mut rng2 = ChaCha20Rng::seed_from_u64(123);

        let pk1 = PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng1);
        let pk2 = PublicKey::from_secret_key(&secret_key, &pk_params, &mut rng2);

        // Keys should be identical with same seed
        for i in 0..ring_dim {
            assert_eq!(pk1.a.into_iter().nth(i), pk2.a.into_iter().nth(i));
            assert_eq!(pk1.b.into_iter().nth(i), pk2.b.into_iter().nth(i));
        }
    }
}
