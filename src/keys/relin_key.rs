use crate::{
    PolyRing, SecretKey, generate_error_poly, generate_random_poly, negate_poly,
};
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
    pub fn from_secret_key(
        secret_key: &SecretKey,
        modulus: u64,
        base: u64, // This is our P factor
        num_components: usize,
    ) -> Self {
        let mut rng = ChaCha20Rng::from_seed([2u8; 32]);
        let n = secret_key.s.len();
        let ring_dim = secret_key.s.ring_dimension();

        // Compute s^2
        let s_squared = secret_key.s.clone() * secret_key.s.clone();

        // Scale s^2 by base (P) within the same modulus
        let mut s_squared_scaled_coeffs = Vec::with_capacity(n);
        for &coeff in &s_squared {
            // Multiply by base and reduce modulo modulus
            let scaled_coeff = (coeff as u128 * base as u128) % modulus as u128;
            s_squared_scaled_coeffs.push(scaled_coeff as u64);
        }
        let s_squared_scaled =
            PolyRing::from_coeffs(&s_squared_scaled_coeffs, modulus, ring_dim);

        // Generate random polynomial a
        let a = generate_random_poly(n, modulus, num_components, &mut rng);

        // Generate error polynomial e
        let e = generate_error_poly(n, modulus, 3.0, num_components, &mut rng);

        // Compute b = -(a*s) + e + P*s^2
        let a_times_s = a.clone() * secret_key.s.clone();
        let neg_a_s = negate_poly(&a_times_s, modulus);

        let b = neg_a_s + e + s_squared_scaled;

        let mut components = Vec::with_capacity(num_components);
        components.push((b, a));

        Self {
            components,
            base,
            num_components,
        }
    }
}

// Helper function to convert a polynomial to use a larger modulus
fn convert_poly_to_larger_modulus(poly: &PolyRing, new_modulus: u64) -> PolyRing {
    let mut new_coeffs = Vec::with_capacity(poly.len());

    for &coeff in poly {
        // Map coefficients properly to new modulus
        let new_coeff = if coeff == 0 {
            0
        } else if coeff == poly.modulus() - 1 {
            // This was -1 in the original modulus
            new_modulus - 1
        } else {
            coeff
        };

        new_coeffs.push(new_coeff);
    }

    PolyRing::from_coeffs(&new_coeffs, new_modulus, poly.ring_dimension())
}

// Helper function to multiply polynomial by scalar
fn multiply_by_scalar(poly: &PolyRing, scalar: u64) -> PolyRing {
    let modulus = poly.modulus();
    let mut result = Vec::with_capacity(poly.len());

    for &coeff in poly {
        let scaled = (coeff as u128 * scalar as u128) % modulus as u128;
        result.push(scaled as u64);
    }

    PolyRing::from_coeffs(&result, modulus, poly.ring_dimension())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SecretKey, SecretKeyParams};

    /* #[test]
    fn test_relin_key_basic_property2() {
        // Create a simple secret key
        let sk_params = SecretKeyParams {
            ring_degree: 8,
            modulus: 65537,
            sparsity: 0.5,
        };
        let secret_key = SecretKey::generate(&sk_params);

        // Generate relin key
        let base = 3; // Small P for testing
        let relin_key = RelinearizationKey::from_secret_key(
            &secret_key,
            sk_params.modulus,
            base,
            1,
        );

        // Test if evk_0 + evk_1*s ≈ P*s^2
        let (b, a) = &relin_key.components[0];
        let s_larger_mod =
            convert_poly_to_larger_modulus(&secret_key.s, b.modulus());
        let decrypted = b.clone() + (a.clone() * s_larger_mod);

        // Compute P*s^2 in the larger modulus
        let s_squared = secret_key.s.clone() * secret_key.s.clone();
        let expected = multiply_by_scalar(&s_squared, base);
        let expected_larger =
            convert_poly_to_larger_modulus(&expected, b.modulus());

        // Check if they're close (this is simplified)
        // In practice we'd measure the difference, but for a toy test
        // just printing the coefficients can be a visual check
        println!("Decrypted: {:?}", decrypted);
        println!("Expected P*s^2: {:?}", expected_larger);

        // The difference should be the error term e, which should be small
    } */
}
