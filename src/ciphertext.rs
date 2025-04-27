use crate::{PolyRing, RelinearizationKey};

/// Represents a ciphertext in the CKKS scheme
#[derive(Clone, Debug)]
pub struct Ciphertext {
    /// c0 component of ciphertext
    pub c0: PolyRing,
    /// c1 component of ciphertext
    pub c1: PolyRing,
    /// c2 component of ciphertext (coefficient of s^2),
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

    pub fn mul(&self, other: &Self, relin_key: &RelinearizationKey) -> Self {
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale for multiplication"
        );

        // Same steps for creating d0, d1, d2...
        let d0 = self.c0.clone() * other.c0.clone();
        let d1 = (self.c0.clone() * other.c1.clone())
            + (self.c1.clone() * other.c0.clone());
        let d2 = self.c1.clone() * other.c1.clone();

        let intermediate = Self {
            c0: d0,
            c1: d1,
            c2: Some(d2),
            scale: self.scale * other.scale,
        };

        // Apply relinearization to eliminate the s^2 term
        let relinearized = intermediate.relinearize(relin_key);

        // Rescale back down to the original scale
        relinearized.rescale(self.scale)
    }

    pub fn relinearize(&self, relin_key: &RelinearizationKey) -> Self {
        if self.c2.is_none() {
            // No s^2 term, nothing to do
            return self.clone();
        }

        let c2 = self.c2.as_ref().unwrap();

        // Use the b and a components directly from the relin_key
        let b0 = &relin_key.b;
        let a0 = &relin_key.a;

        // Compute c2 * (b0 + a0*s)
        let c0_new = self.c0.clone() + (c2.clone() * b0.clone());
        let c1_new = self.c1.clone() + (c2.clone() * a0.clone());

        Self {
            c0: c0_new,
            c1: c1_new,
            c2: None,
            scale: self.scale,
        }
    }

    pub fn rescale(&self, target_scale: f64) -> Self {
        // Calculate how much we need to scale down by
        let scale_factor = self.scale / target_scale;
        println!(
            "Rescaling from {} to {} (factor: {})",
            self.scale, target_scale, scale_factor
        );

        // Rescale the polynomials by dividing coefficients
        let c0_scaled = rescale_poly(&self.c0, scale_factor);
        let c1_scaled = rescale_poly(&self.c1, scale_factor);

        Self {
            c0: c0_scaled,
            c1: c1_scaled,
            c2: None,
            scale: target_scale,
        }
    }
}

fn rescale_poly(poly: &PolyRing, scale_factor: f64) -> PolyRing {
    let modulus = poly.modulus();
    let half_modulus = modulus / 2;
    let ring_dim = poly.ring_dim();

    // Convert scale_factor to power-of-2 shift amount
    let shift_bits = scale_factor.log2().round() as u32;

    let mut new_coeffs = Vec::with_capacity(poly.len());

    for &coeff in poly {
        // Convert to centered representation (-q/2, q/2)
        let centered_coeff = if coeff > half_modulus {
            (coeff as i128) - (modulus as i128)
        } else {
            coeff as i128
        };

        // Apply bit shifting with proper rounding
        let round_bit = 1i128 << (shift_bits - 1);
        let scaled_coeff = if centered_coeff >= 0 {
            (centered_coeff + round_bit) >> shift_bits
        } else {
            (centered_coeff - round_bit) >> shift_bits
        };

        // Convert back to modular representation
        let scaled = if scaled_coeff < 0 {
            (modulus as i128 + scaled_coeff) as u64
        } else {
            scaled_coeff as u64
        } % modulus;

        new_coeffs.push(scaled);
    }

    PolyRing::from_coeffs(&new_coeffs, modulus, ring_dim)
}
