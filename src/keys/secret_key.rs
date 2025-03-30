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
    /// Polynomial degree (must be power of 2)
    pub n: usize,
    /// Modulus for the polynomial coefficients
    pub modulus: u64,
    /// Controls the sparsity of the ternary distribution (0 to 1)
    /// Higher values mean more zeros in the polynomial
    pub sparsity: f64,
}

impl SecretKey {
    /// Generate a new secret key with ternary coefficients {-1, 0, 1}
    pub fn generate(params: &SecretKeyParams) -> Self {
        // Create RNG with proper seeding
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]); // In production use entropy-based seed

        // Generate ternary polynomial for secret key
        let mut coeffs = Vec::with_capacity(params.n);

        for _ in 0..params.n {
            // First decide if coefficient is zero based on sparsity
            if rng.random::<f64>() < params.sparsity {
                coeffs.push(0);
            } else {
                // Otherwise, generate -1 or 1 with equal probability
                let coeff = if rng.random::<bool>() {
                    1u64
                } else {
                    // -1 is represented as (modulus - 1) in modular arithmetic
                    params.modulus - 1
                };
                coeffs.push(coeff);
            }
        }

        // Create the polynomial
        let s = PolyRing::from_coeffs(&coeffs, params.modulus);

        Self { s }
    }
}
