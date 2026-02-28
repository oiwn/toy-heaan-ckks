//! CKKS encoder/decoder for `RnsPoly` plaintexts.
//!
//! Implements the canonical embedding:
//!   encode: scale → `special_idft` → round → `RnsPoly::from_coeffs`
//!   decode: `poly.to_coeffs()` → lift → `special_dft` → unscale
//!
//! The CRT reconstruction needed for decoding is fully handled by
//! `RnsBasis::reconstruct_centered_coeff` inside `to_coeffs()`, so no
//! `BigUint`/`BoxedUint` machinery is required here.

use std::sync::Arc;

use rustfft::num_complex::Complex64;

use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
use crate::rings::traits::PolyRing;

use super::special_fft::{
    VandermondeTables, build_slots_from_complex, special_dft, special_idft,
};

/// A plaintext polynomial produced by `CkksEncoder::encode`.
pub struct Plaintext<const DEGREE: usize> {
    /// The scaled, rounded polynomial in coefficient domain.
    pub poly: RnsPoly<DEGREE>,
    /// Number of bits in the scaling factor Δ = 2^scale_bits.
    pub scale_bits: u32,
    /// Number of slots encoded (determines how many values `decode` returns).
    pub slots: usize,
}

/// CKKS encoder/decoder targeting `RnsPoly<DEGREE>`.
///
/// # Example
/// ```ignore
/// let encoder = CkksEncoder::<8>::new(30);
/// let basis   = Arc::new(RnsBasis::<8>::new(primes)?);
/// let pt      = encoder.encode(&[1.5, -2.0], basis.clone());
/// let out     = encoder.decode(&pt);             // ≈ [1.5, -2.0]
/// ```
pub struct CkksEncoder<const DEGREE: usize> {
    scale_bits: u32,
    tables: VandermondeTables<DEGREE>,
}

impl<const DEGREE: usize> CkksEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> Self {
        assert!(
            DEGREE.is_power_of_two(),
            "CkksEncoder: DEGREE must be a power of two"
        );
        assert!(scale_bits > 0, "CkksEncoder: scale_bits must be positive");
        Self {
            scale_bits,
            tables: VandermondeTables::new(),
        }
    }

    pub fn scale_factor(&self) -> f64 {
        2f64.powi(self.scale_bits as i32)
    }

    pub fn max_slots(&self) -> usize {
        DEGREE / 2
    }

    // ── Encoding ─────────────────────────────────────────────────────────────

    /// Encodes real values into a `Plaintext`.
    ///
    /// Each value occupies one complex slot (imaginary = 0). Hermitian symmetry
    /// in the slot vector guarantees the IDFT yields real polynomial coefficients.
    /// At most `DEGREE/2` values may be encoded.
    pub fn encode(
        &self,
        values: &[f64],
        basis: Arc<RnsBasis<DEGREE>>,
    ) -> Plaintext<DEGREE> {
        assert!(
            values.len() <= DEGREE / 2,
            "encode: {} values exceed max slots {}",
            values.len(),
            DEGREE / 2
        );
        let delta = self.scale_factor();
        let complex: Vec<Complex64> = values
            .iter()
            .map(|&v| Complex64::new(v * delta, 0.0))
            .collect();
        self.encode_inner(&complex, values.len(), basis)
    }

    /// Encodes complex values into a `Plaintext`.
    pub fn encode_complex(
        &self,
        values: &[Complex64],
        basis: Arc<RnsBasis<DEGREE>>,
    ) -> Plaintext<DEGREE> {
        assert!(
            values.len() <= DEGREE / 2,
            "encode_complex: {} values exceed max slots {}",
            values.len(),
            DEGREE / 2
        );
        let delta = self.scale_factor();
        let scaled: Vec<Complex64> = values.iter().map(|v| v * delta).collect();
        self.encode_inner(&scaled, values.len(), basis)
    }

    fn encode_inner(
        &self,
        scaled: &[Complex64],
        slots: usize,
        basis: Arc<RnsBasis<DEGREE>>,
    ) -> Plaintext<DEGREE> {
        // Build conjugate-symmetric N-slot vector, then IDFT → coefficient vector.
        let slot_vec = build_slots_from_complex::<DEGREE>(scaled);
        let coeff_vec = special_idft::<DEGREE>(&slot_vec, &self.tables);

        // Coefficients are real by Hermitian symmetry; round to integers.
        let mut int_coeffs = [0i64; DEGREE];
        for (i, c) in coeff_vec.iter().enumerate() {
            int_coeffs[i] = c.re.round() as i64;
        }

        Plaintext {
            poly: RnsPoly::from_coeffs(&int_coeffs, basis),
            scale_bits: self.scale_bits,
            slots,
        }
    }

    // ── Decoding ─────────────────────────────────────────────────────────────

    /// Decodes a `Plaintext` back to real values.
    ///
    /// Error is bounded by the rounding noise from encoding (≈ 1/Δ per slot).
    pub fn decode(&self, pt: &Plaintext<DEGREE>) -> Vec<f64> {
        self.decode_complex(pt).into_iter().map(|s| s.re).collect()
    }

    /// Decodes a `Plaintext` back to complex values.
    pub fn decode_complex(&self, pt: &Plaintext<DEGREE>) -> Vec<Complex64> {
        let delta = 2f64.powi(pt.scale_bits as i32);

        // CRT-reconstruct centered integer coefficients.
        let int_coeffs = pt.poly.to_coeffs();

        // Lift to complex and apply the forward canonical embedding (DFT).
        let c_coeffs: Vec<Complex64> = int_coeffs
            .iter()
            .map(|&x| Complex64::new(x as f64, 0.0))
            .collect();
        let slot_vec = special_dft::<DEGREE>(&c_coeffs, &self.tables);

        // Return first `slots` values, unscaled by Δ.
        slot_vec
            .into_iter()
            .take(pt.slots)
            .map(|s| s / delta)
            .collect()
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // For degree 8, NTT-friendly primes must satisfy p ≡ 1 (mod 16).
    // 97 ≡ 1 (mod 16), 113 ≡ 1 (mod 16).
    // With scale_bits=5 (Δ=32), encoded coefficients fit comfortably mod 97.
    fn basis() -> Arc<RnsBasis<8>> {
        Arc::new(RnsBasis::new(vec![97, 113]).unwrap())
    }

    #[test]
    fn roundtrip_real_values() {
        let encoder = CkksEncoder::<8>::new(5);
        let values = [1.0f64, -1.0, 0.5, -0.5];

        let pt = encoder.encode(&values, basis());
        let decoded = encoder.decode(&pt);

        for (orig, got) in values.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig, got, epsilon = 0.1);
        }
    }

    #[test]
    fn roundtrip_complex_values() {
        let encoder = CkksEncoder::<8>::new(5);
        let values = [Complex64::new(1.0, 0.5), Complex64::new(-0.5, 0.25)];

        let pt = encoder.encode_complex(&values, basis());
        let decoded = encoder.decode_complex(&pt);

        for (orig, got) in values.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig.re, got.re, epsilon = 0.1);
            assert_relative_eq!(orig.im, got.im, epsilon = 0.1);
        }
    }

    #[test]
    fn single_value_roundtrip() {
        let encoder = CkksEncoder::<8>::new(5);
        let pt = encoder.encode(&[3.0], basis());
        let out = encoder.decode(&pt);
        assert_relative_eq!(out[0], 3.0, epsilon = 0.1);
    }

    #[test]
    fn slot_count_preserved() {
        let encoder = CkksEncoder::<8>::new(5);
        let pt = encoder.encode(&[1.0, 2.0, 3.0], basis());
        let out = encoder.decode(&pt);
        assert_eq!(out.len(), 3);
    }

    #[test]
    fn max_slots_is_half_degree() {
        let encoder = CkksEncoder::<8>::new(10);
        assert_eq!(encoder.max_slots(), 4);
    }

    #[test]
    #[should_panic(expected = "exceed max slots")]
    fn panics_on_too_many_values() {
        let encoder = CkksEncoder::<8>::new(5);
        encoder.encode(&[0.0; 5], basis()); // 5 > 8/2
    }
}
