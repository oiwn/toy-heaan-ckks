//! CKKS encoder/decoder over any [`PolyRing`] backend.
//!
//! Implements the canonical embedding:
//!   encode: scale → `special_idft` → round → [`PolyRing::from_rounded_f64`]
//!   decode: [`PolyRing::to_centered_f64`] → `special_dft` → unscale
//!
//! Coefficients are staged in f64, not i64: a rounded coefficient may be any
//! size below `Q/2`. The backend reduces it into the ring and reconstructs
//! it on decode (for RNS: per-prime residues and CRT), so neither direction
//! is limited to i64.

use rustfft::num_complex::Complex64;

use crate::PolyRing;
use crate::ckks::types::Plaintext;
use crate::rings::RingResult;

use super::special_fft::{
    VandermondeTables, build_slots_from_complex, special_dft, special_idft,
};

/// CKKS encoder/decoder for polynomials of degree `DEGREE`.
///
/// # Example
/// ```ignore
/// let encoder = CkksEncoder::<8>::new(30);
/// let basis   = Arc::new(RnsBasis::<8>::new(primes)?);
/// let pt: Plaintext<RnsPoly<8>, 8> = encoder.encode(&[1.5, -2.0], basis.clone());
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
    /// Panicking wrapper around [`CkksEncoder::try_encode`].
    ///
    /// # Panics
    /// Panics if more than `DEGREE/2` values are given, or if a scaled
    /// coefficient does not fit the modulus (see `try_encode`).
    pub fn encode<P: PolyRing<DEGREE>>(
        &self,
        values: &[f64],
        context: P::Context,
    ) -> Plaintext<P, DEGREE> {
        self.try_encode(values, context)
            .unwrap_or_else(|err| panic!("encode: {err}"))
    }

    /// Encodes real values into a `Plaintext`.
    ///
    /// Each value occupies one complex slot (imaginary = 0). Hermitian symmetry
    /// in the slot vector guarantees the IDFT yields real polynomial coefficients.
    /// At most `DEGREE/2` values may be encoded.
    ///
    /// Coefficients are staged in f64, so precision is ~2^-52 relative to the
    /// largest scaled coefficient: a slot of 1e6 next to a slot of 1e-4 costs
    /// the small slot about 1e-9 of absolute accuracy.
    ///
    /// # Errors
    /// Returns [`RingError::CoefficientTooLarge`](crate::rings::RingError) if
    /// a scaled coefficient is not finite or its magnitude reaches
    /// `2^(modulus_bits - 1)` (a conservative `Q/2` bound).
    ///
    /// # Panics
    /// Panics if more than `DEGREE/2` values are given.
    pub fn try_encode<P: PolyRing<DEGREE>>(
        &self,
        values: &[f64],
        context: P::Context,
    ) -> RingResult<Plaintext<P, DEGREE>> {
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
        self.encode_inner(&complex, values.len(), context)
    }

    /// Encodes complex values into a `Plaintext`.
    ///
    /// Panicking wrapper around [`CkksEncoder::try_encode_complex`].
    ///
    /// # Panics
    /// Panics if more than `DEGREE/2` values are given, or if a scaled
    /// coefficient does not fit the modulus.
    pub fn encode_complex<P: PolyRing<DEGREE>>(
        &self,
        values: &[Complex64],
        context: P::Context,
    ) -> Plaintext<P, DEGREE> {
        self.try_encode_complex(values, context)
            .unwrap_or_else(|err| panic!("encode_complex: {err}"))
    }

    /// Encodes complex values into a `Plaintext`.
    ///
    /// # Errors
    /// Returns [`RingError::CoefficientTooLarge`](crate::rings::RingError)
    /// under the same rule as [`CkksEncoder::try_encode`].
    ///
    /// # Panics
    /// Panics if more than `DEGREE/2` values are given.
    pub fn try_encode_complex<P: PolyRing<DEGREE>>(
        &self,
        values: &[Complex64],
        context: P::Context,
    ) -> RingResult<Plaintext<P, DEGREE>> {
        assert!(
            values.len() <= DEGREE / 2,
            "encode_complex: {} values exceed max slots {}",
            values.len(),
            DEGREE / 2
        );
        let delta = self.scale_factor();
        let scaled: Vec<Complex64> = values.iter().map(|v| v * delta).collect();
        self.encode_inner(&scaled, values.len(), context)
    }

    fn encode_inner<P: PolyRing<DEGREE>>(
        &self,
        scaled: &[Complex64],
        slots: usize,
        context: P::Context,
    ) -> RingResult<Plaintext<P, DEGREE>> {
        // Build conjugate-symmetric N-slot vector, then IDFT → coefficient vector.
        let slot_vec = build_slots_from_complex::<DEGREE>(scaled);
        let coeff_vec = special_idft::<DEGREE>(&slot_vec, &self.tables);

        // Coefficients are real by Hermitian symmetry; round to integers but
        // keep them in f64 so they may exceed the i64 range.
        let mut rounded = [0f64; DEGREE];
        for (slot, c) in rounded.iter_mut().zip(coeff_vec.iter()) {
            *slot = c.re.round();
        }
        self.plaintext_from_rounded(&rounded, slots, context)
    }

    // ── Scalar encoding ──────────────────────────────────────────────────────

    /// Encodes one real value into every slot (SEAL-style scalar encode).
    ///
    /// Panicking wrapper around [`CkksEncoder::try_encode_constant`].
    ///
    /// # Panics
    /// Panics if `value * 2^scale_bits` does not fit the modulus.
    pub fn encode_constant<P: PolyRing<DEGREE>>(
        &self,
        value: f64,
        context: P::Context,
    ) -> Plaintext<P, DEGREE> {
        self.try_encode_constant(value, context)
            .unwrap_or_else(|err| panic!("encode_constant: {err}"))
    }

    /// Encodes one real value into every slot (SEAL-style scalar encode).
    ///
    /// A value repeated in all `DEGREE/2` slots is the constant polynomial
    /// `round(value * 2^scale_bits)`, so this skips the DFT entirely: it is
    /// exact and O(DEGREE). Scalars must always fill every slot — spreading
    /// a constant over fewer slots gives a non-constant polynomial whose
    /// coefficients can round away (e.g. `0.5` at 1 bit).
    ///
    /// # Errors
    /// Returns [`RingError::CoefficientTooLarge`](crate::rings::RingError)
    /// under the same rule as [`CkksEncoder::try_encode`].
    pub fn try_encode_constant<P: PolyRing<DEGREE>>(
        &self,
        value: f64,
        context: P::Context,
    ) -> RingResult<Plaintext<P, DEGREE>> {
        let mut rounded = [0f64; DEGREE];
        rounded[0] = (value * self.scale_factor()).round();
        self.plaintext_from_rounded(&rounded, DEGREE / 2, context)
    }

    /// Hands the rounded coefficients to the backend (range check included).
    fn plaintext_from_rounded<P: PolyRing<DEGREE>>(
        &self,
        rounded: &[f64; DEGREE],
        slots: usize,
        context: P::Context,
    ) -> RingResult<Plaintext<P, DEGREE>> {
        Ok(Plaintext {
            poly: P::from_rounded_f64(rounded, &context)?,
            scale_bits: self.scale_bits,
            slots,
        })
    }

    // ── Decoding ─────────────────────────────────────────────────────────────

    /// Decodes a `Plaintext` back to real values.
    ///
    /// Error is bounded by the rounding noise from encoding (≈ 1/Δ per slot).
    pub fn decode<P: PolyRing<DEGREE>>(
        &self,
        pt: &Plaintext<P, DEGREE>,
    ) -> Vec<f64> {
        self.decode_complex(pt).into_iter().map(|s| s.re).collect()
    }

    /// Decodes a `Plaintext` back to complex values.
    pub fn decode_complex<P: PolyRing<DEGREE>>(
        &self,
        pt: &Plaintext<P, DEGREE>,
    ) -> Vec<Complex64> {
        let delta = 2f64.powi(pt.scale_bits as i32);

        // Centered coefficients straight into f64 (no i64 limit), then the
        // forward canonical embedding (DFT).
        let c_coeffs: Vec<Complex64> = pt
            .poly
            .to_centered_f64()
            .iter()
            .map(|&c| Complex64::new(c, 0.0))
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
    use std::sync::Arc;

    use super::*;
    use crate::rings::RingError;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
    use approx::assert_relative_eq;

    type Pt<const D: usize> = Plaintext<RnsPoly<D>, D>;

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

        let pt: Pt<8> = encoder.encode(&values, basis());
        let decoded = encoder.decode(&pt);

        for (orig, got) in values.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig, got, epsilon = 0.1);
        }
    }

    #[test]
    fn roundtrip_complex_values() {
        let encoder = CkksEncoder::<8>::new(5);
        let values = [Complex64::new(1.0, 0.5), Complex64::new(-0.5, 0.25)];

        let pt: Pt<8> = encoder.encode_complex(&values, basis());
        let decoded = encoder.decode_complex(&pt);

        for (orig, got) in values.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig.re, got.re, epsilon = 0.1);
            assert_relative_eq!(orig.im, got.im, epsilon = 0.1);
        }
    }

    #[test]
    fn single_value_roundtrip() {
        let encoder = CkksEncoder::<8>::new(5);
        let pt: Pt<8> = encoder.encode(&[3.0], basis());
        let out = encoder.decode(&pt);
        assert_relative_eq!(out[0], 3.0, epsilon = 0.1);
    }

    #[test]
    fn slot_count_preserved() {
        let encoder = CkksEncoder::<8>::new(5);
        let pt: Pt<8> = encoder.encode(&[1.0, 2.0, 3.0], basis());
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
        let _: Pt<8> = encoder.encode(&[0.0; 5], basis()); // 5 > 8/2
    }

    // ── Wide staging (coefficients beyond i64) ──────────────────────────────

    /// Four 61-bit NTT primes for DEGREE = 16: Q ~ 2^244.
    fn wide_basis() -> Arc<RnsBasis<16>> {
        let primes = crate::math::generate_primes(61, 4, 16);
        Arc::new(RnsBasis::new(primes).unwrap())
    }

    /// At scale 2^61, |x| >= 4 overflowed the old i64 staging and saturated.
    #[test]
    fn roundtrip_values_beyond_i64_staging() {
        let encoder = CkksEncoder::<16>::new(61);
        let values = [10.0f64, -250.0, 1.0e6, 6.49, -3.5, 0.0, 1.0e-4, 42.0];

        let pt: Pt<16> = encoder.try_encode(&values, wide_basis()).unwrap();
        let decoded = encoder.decode(&pt);

        // f64 staging keeps ~52 bits relative to the largest slot, so the
        // absolute error scales with max |value| (1e6 here -> ~1e-9).
        let tolerance = 1e-13 * values.iter().fold(1.0f64, |m, v| m.max(v.abs()));
        for (orig, got) in values.iter().zip(decoded.iter()) {
            assert!(
                (orig - got).abs() < tolerance,
                "expected {orig}, decoded {got}"
            );
        }
    }

    #[test]
    fn rejects_coefficients_at_half_modulus() {
        let encoder = CkksEncoder::<16>::new(1);
        let basis = wide_basis();
        let max_bits = basis.total_bits() - 1;

        // c in every slot encodes to the single coefficient c * 2^1, so
        // c = 2^(max_bits - 1) lands exactly on the 2^max_bits limit.
        let limit_value = 2f64.powi(max_bits as i32 - 1);
        let below = limit_value - 2f64.powi(max_bits as i32 - 53);
        let ok: RingResult<Pt<16>> = encoder.try_encode(&[below; 8], basis.clone());
        assert!(ok.is_ok());

        let err: RingError = encoder
            .try_encode::<RnsPoly<16>>(&[limit_value; 8], basis.clone())
            .unwrap_err();
        let RingError::CoefficientTooLarge {
            value,
            max_bits: bits,
        } = err
        else {
            panic!("expected CoefficientTooLarge, got {err}");
        };
        assert_eq!(bits, max_bits);
        assert!(value.abs() >= 2f64.powi(max_bits as i32));
        assert!(
            encoder
                .try_encode::<RnsPoly<16>>(&[-limit_value; 8], basis)
                .is_err()
        );
    }

    /// A scalar fills every slot and is exactly the constant polynomial.
    #[test]
    fn encode_constant_is_exact_constant_polynomial() {
        let basis = wide_basis();
        let encoder = CkksEncoder::<16>::new(1);

        // 0.5 at 1 bit -> the integer 1 in coefficient 0, zeros elsewhere.
        let pt: Pt<16> = encoder.encode_constant(0.5, basis.clone());
        let mut want = [0f64; 16];
        want[0] = 1.0;
        assert_eq!(pt.poly.to_centered_f64(), want);

        let encoder = CkksEncoder::<16>::new(61);
        let pt: Pt<16> = encoder.encode_constant(-635.25, basis);
        let decoded = encoder.decode(&pt);
        assert_eq!(decoded.len(), 8, "a scalar fills all DEGREE/2 slots");
        for got in decoded {
            assert!((got + 635.25).abs() < 1e-9, "decoded {got}");
        }
    }

    #[test]
    fn encode_constant_rejects_out_of_range() {
        let basis = wide_basis();
        let encoder = CkksEncoder::<16>::new(61);
        assert!(
            encoder
                .try_encode_constant::<RnsPoly<16>>(2f64.powi(200), basis.clone())
                .is_err()
        );
        assert!(
            encoder
                .try_encode_constant::<RnsPoly<16>>(f64::NAN, basis)
                .is_err()
        );
    }

    #[test]
    fn rejects_non_finite_values() {
        let encoder = CkksEncoder::<16>::new(61);
        for bad in [f64::NAN, f64::INFINITY] {
            let result: RingResult<Pt<16>> =
                encoder.try_encode(&[bad], wide_basis());
            assert!(result.is_err());
        }
    }

    #[test]
    #[should_panic(expected = "does not fit the modulus")]
    fn encode_panics_with_clear_message() {
        // Basis 97 * 113 holds |coeff| < 2^11; 1000 * 2^5 does not fit.
        let encoder = CkksEncoder::<8>::new(5);
        let _: Pt<8> = encoder.encode(&[1000.0], basis());
    }
}
