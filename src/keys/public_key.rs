use crate::{
    PolyRing, RelinearizationKey, SecretKey, generate_error_poly,
    generate_random_poly, negate_poly,
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

/// A public key in the CKKS scheme
pub struct PublicKey {
    /// The first component of the public key (b = -a·s + e)
    pub b: PolyRing,
    /// The second component of the public key (random polynomial a)
    pub a: PolyRing,

    /// Relinearization key for multiplication
    pub relin_key: RelinearizationKey,
}

/// Parameters for public key generation
pub struct PublicKeyParams {
    /// Polynomial degree (must be power of 2)
    pub n: usize,
    /// Modulus for the polynomial coefficients
    pub modulus: u64,
    /// Variance for error distribution used in public key generation
    pub error_variance: f64,
    /// Base for relinearization key decomposition (e.g., 2^16)
    pub relin_base: u64,
    /// Number of components for relinearization key
    pub relin_components: usize,
}

impl PublicKey {
    /// Generate a public key from a secret key using the provided parameters
    pub fn from_secret_key(
        secret_key: &SecretKey,
        params: &PublicKeyParams,
    ) -> Self {
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

        // Compute b = -(a * s) + e
        let modulus = params.modulus;
        let a_times_s = a.clone() * secret_key.s.clone();

        // Create negative of (a * s)
        let neg_a_s = negate_poly(&a_times_s, modulus);

        // b = -(a·s) + e
        let b = neg_a_s + e;

        // Generate relinearization key
        let relin_key = RelinearizationKey::from_secret_key(
            secret_key,
            params.modulus,
            params.relin_base,
            params.relin_components,
        );

        PublicKey { a, b, relin_key }
    }
}
