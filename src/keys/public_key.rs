use crate::{PolyRing, SecretKey};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

/// A public key in the CKKS scheme
pub struct PublicKey {
    /// The first component of the public key (b = -a·s + e)
    pub b: PolyRing,
    /// The second component of the public key (random polynomial a)
    pub a: PolyRing,
}

/// Parameters for public key generation
pub struct PublicKeyParams {
    /// Polynomial degree (must be power of 2)
    pub n: usize,
    /// Modulus for the polynomial coefficients
    pub modulus: u64,
    /// Variance for error distribution used in public key generation
    pub error_variance: f64,
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

        // Create negative of (a·s)
        let neg_a_s = negate_poly(&a_times_s, modulus);

        // b = -(a·s) + e
        let b = neg_a_s + e;

        PublicKey { a, b }
    }
}

// Helper functions needed for public key generation

/// Generate a polynomial with random coefficients in [0, modulus)
fn generate_random_poly<T: Rng>(n: usize, modulus: u64, rng: &mut T) -> PolyRing {
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Generate a random value modulo q
        let random_val = rng.random::<u64>() % modulus;
        coeffs.push(random_val);
    }

    PolyRing::from_coeffs(&coeffs, modulus)
}

/// Generate a polynomial with small random coefficients (error polynomial)
/// using discretized Gaussian distribution
fn generate_error_poly<T: Rng>(
    n: usize,
    modulus: u64,
    variance: f64,
    rng: &mut T,
) -> PolyRing {
    // Get standard deviation from variance
    let std_dev = variance.sqrt();
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Simplified approximation of Gaussian noise
        // Sum three random values and shift to get approximately normal distribution
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
