use crate::{
    PolyRing, PublicKey, SecretKey, generate_error_poly, generate_ternary_poly,
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

/// Represents a ciphertext in the CKKS scheme
#[derive(Clone, Debug)]
pub struct Ciphertext {
    /// c0 component of ciphertext
    pub c0: PolyRing,
    /// c1 component of ciphertext
    pub c1: PolyRing,
    /// Tracks the scaling factor used
    pub scale: f64,
}

impl Ciphertext {
    /// Homomorphic addition of ciphertexts
    pub fn add(&self, other: &Self) -> Self {
        // Check scaling factors and moduli
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale"
        );

        // Add corresponding polynomials component-wise
        let c0 = self.c0.clone() + other.c0.clone();
        let c1 = self.c1.clone() + other.c1.clone();

        Self {
            c0,
            c1,
            scale: self.scale,
        }
    }

    /// Rescales the ciphertext after multiplication
    /// This divides by the scaling factor and updates the scale
    pub fn rescale(&self, target_scale: f64) -> Self {
        // let modulus = self.c0.modulus();
        let scale_factor = self.scale / target_scale;

        // For simplicity, we'll assume scale_factor is a power of 2
        let rounding_factor = scale_factor.round() as u64;

        // Create scaled versions of c0 and c1
        let scaled_c0 = rescale_poly(&self.c0, rounding_factor);
        let scaled_c1 = rescale_poly(&self.c1, rounding_factor);

        Self {
            c0: scaled_c0,
            c1: scaled_c1,
            scale: target_scale,
        }
    }
}

/// Encrypt a plaintext polynomial using the public key
pub fn encrypt(
    plaintext: &PolyRing,
    public_key: &PublicKey,
    scale: f64,
) -> Ciphertext {
    let modulus = plaintext.modulus();
    let mut rng = ChaCha20Rng::from_seed([2u8; 32]); // Different seed from keys

    // Generate small random polynomials for encryption
    let e1 = generate_error_poly(plaintext.len(), modulus, 3.0, &mut rng);
    let e2 = generate_error_poly(plaintext.len(), modulus, 3.0, &mut rng);

    // Generate random ephemeral value
    let u = generate_ternary_poly(plaintext.len(), modulus, 0.5, &mut rng);

    // c0 = b*u + e1 + plaintext
    let c0 = (public_key.b.clone() * u.clone()) + e1 + plaintext.clone();

    // c1 = a*u + e2
    let c1 = (public_key.a.clone() * u) + e2;

    Ciphertext { c0, c1, scale }
}

/// Decrypt a ciphertext using the secret key
pub fn decrypt(ciphertext: &Ciphertext, secret_key: &SecretKey) -> PolyRing {
    // Compute c0 + c1*s
    ciphertext.c0.clone() + (ciphertext.c1.clone() * secret_key.s.clone())
}

/// Helper function to rescale a polynomial
fn rescale_poly(poly: &PolyRing, scale_factor: u64) -> PolyRing {
    let modulus = poly.modulus();
    let mut new_coeffs = Vec::with_capacity(poly.len());

    for &coeff in poly {
        // Divide by scale factor and round to the nearest integer
        let scaled = (coeff as f64 / scale_factor as f64).round() as u64 % modulus;
        new_coeffs.push(scaled);
    }

    PolyRing::from_coeffs(&new_coeffs, modulus)
}
