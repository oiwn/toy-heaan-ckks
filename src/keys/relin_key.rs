use crate::{PolyRing, SecretKey, generate_error_poly, generate_random_poly};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

/// Relinearization key used to transform ciphertexts after multiplication
#[derive(Clone)]
pub struct RelinearizationKey {
    /// The relinearization key components: a polynomial vector containing
    /// encryptions of s² under the public key
    pub components: Vec<(PolyRing, PolyRing)>,
    /// The decomposition base used during key generation
    pub base: u64,
    /// Number of decomposition components
    pub num_components: usize,
}

impl RelinearizationKey {
    /// Generate a relinearization key from a secret key
    pub fn from_secret_key(
        secret_key: &SecretKey,
        modulus: u64,
        base: u64,
        num_components: usize,
    ) -> Self {
        let mut rng = ChaCha20Rng::from_seed([2u8; 32]); // Different seed from other keys
        let n = secret_key.s.len();

        // Compute s²
        let s_squared = secret_key.s.clone() * secret_key.s.clone();

        // For each component of the decomposition base
        let mut components = Vec::with_capacity(num_components);

        // In a production system we'd decompose s² into digits with respect to the base
        // For simplicity in this toy implementation, we'll just encrypt s² directly

        for i in 0..num_components {
            // Generate random polynomial a_i
            let a_i = generate_random_poly(n, modulus, &mut rng);

            // Generate small error polynomial e_i
            let e_i = generate_error_poly(n, modulus, 3.0, &mut rng);

            // Scale s² by base^i (for the digit decomposition)
            let mut s_squared_scaled = s_squared.clone();
            for _ in 0..i {
                // This is simplified - we'd normally do digital decomposition
                s_squared_scaled = multiply_by_scalar(&s_squared_scaled, base);
            }

            // Compute b_i = -(a_i · s) + e_i + base^i · s²
            let a_i_times_s = a_i.clone() * secret_key.s.clone();
            let mut b_i = negate_poly(&a_i_times_s, modulus);
            b_i = b_i + e_i;
            b_i = b_i + s_squared_scaled;

            components.push((b_i, a_i));
        }

        Self {
            components,
            base,
            num_components,
        }
    }
}

// Helper function to multiply polynomial by scalar
fn multiply_by_scalar(poly: &PolyRing, scalar: u64) -> PolyRing {
    let modulus = poly.modulus();
    let mut result = Vec::with_capacity(poly.len());

    for &coeff in poly {
        let scaled = (coeff as u128 * scalar as u128) % modulus as u128;
        result.push(scaled as u64);
    }

    PolyRing::from_unsigned_coeffs(&result, modulus, 8)
}

// You'll need the negate_poly function from your common.rs
// Either move it back to public or reimplement it here
fn negate_poly(poly: &PolyRing, modulus: u64) -> PolyRing {
    let mut neg_coeffs = Vec::with_capacity(poly.len());

    for &coeff in poly {
        if coeff == 0 {
            neg_coeffs.push(0);
        } else {
            neg_coeffs.push(modulus - coeff);
        }
    }

    PolyRing::from_unsigned_coeffs(&neg_coeffs, modulus, 8)
}
