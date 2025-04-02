use crate::{
    PolyRing, PublicKey, RelinearizationKey, SecretKey, generate_error_poly,
    generate_ternary_poly,
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
    /// c2 component of ciphertext (coefficient of s²),
    /// present only after multiplication before relinearization
    pub c2: Option<PolyRing>,
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
            c2: None,
            scale: self.scale,
        }
    }
    pub fn relinearize(&self, relin_key: &RelinearizationKey) -> Self {
        // For a ciphertext (c0, c1, c2) of the form c0 + c1·s + c2·s²
        // we use the relin_key to replace c2·s² with an encryption of the same value

        // In this toy implementation, we'll just use the first component of the key
        // In a full implementation, we'd decompose c2 into digits according to the base

        if self.c2.is_none() {
            // No s² term, nothing to do
            return self.clone();
        }

        let c2 = self.c2.as_ref().unwrap();

        // Extract the first component of the relinearization key
        let (b0, a0) = &relin_key.components[0];

        // Compute c2 · (b0 + a0·s)
        let c0_new = self.c0.clone() + (c2.clone() * b0.clone());
        let c1_new = self.c1.clone() + (c2.clone() * a0.clone());

        Self {
            c0: c0_new,
            c1: c1_new,
            c2: None,
            scale: self.scale,
        }
    }

    /// Homomorphic multiplication of ciphertexts with relinearization
    pub fn mul(&self, other: &Self, relin_key: &RelinearizationKey) -> Self {
        // Check scaling factors
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale for multiplication"
        );

        // For tensor product of ciphertexts: (c0, c1) * (d0, d1)
        // Resulting in c0·d0 + (c0·d1 + c1·d0)·s + c1·d1·s²

        let d0 = self.c0.clone() * other.c0.clone();
        let d1 = (self.c0.clone() * other.c1.clone())
            + (self.c1.clone() * other.c0.clone());
        let d2 = self.c1.clone() * other.c1.clone();

        // Create intermediate ciphertext with s² term
        let result = Self {
            c0: d0,
            c1: d1,
            c2: Some(d2), // We need to add the c2 field to Ciphertext
            scale: self.scale * other.scale,
        };

        // Apply relinearization to eliminate the s² term
        let relinearized = result.relinearize(relin_key);
        // Then rescale back down to the original scale
        relinearized.rescale(self.scale)
    }

    /// Homomorphic multiplication of ciphertexts
    pub fn mul_old(&self, other: &Self) -> Self {
        // Check scaling factors and moduli
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale for multiplication"
        );

        // For tensor product of ciphertexts: (c0, c1) * (d0, d1)
        // c0*d0 + c0*d1*s + c1*d0*s + c1*d1*s²
        // We need to relinearize to bring it back to the form c0' + c1'*s

        // d0 = c0*d0 (constant term)
        let d0 = self.c0.clone() * other.c0.clone();

        // d1 = c0*d1 + c1*d0 (coefficient of s)
        let d1 = (self.c0.clone() * other.c1.clone())
            + (self.c1.clone() * other.c0.clone());

        // d2 = c1*d1 (coefficient of s²) - this requires relinearization
        // let d2 = self.c1.clone() * other.c1.clone();

        // For this toy implementation, we'll do a naive relinearization
        // In a production implementation, this would use evaluation keys

        // The scale doubles after multiplication
        let new_scale = self.scale * other.scale;

        Self {
            c0: d0,
            c1: d1,
            c2: None,
            scale: new_scale,
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
            c2: None,
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

    Ciphertext {
        c0,
        c1,
        c2: None,
        scale,
    }
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

    PolyRing::from_coeffs(&new_coeffs, modulus, poly.ring_degree())
}
