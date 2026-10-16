//! The ring backend trait: the only view the CKKS layer has of polynomials.
//!
//! [`PolyRing`] describes `R_Q = Z_Q[X]/(X^N + 1)` at one level of a modulus
//! chain. A backend decides how it stores elements (RNS channels, one big
//! integer, NTT or coefficient form); nothing in the trait depends on that
//! choice. Concrete backends live in the `backends` module.

use std::ops::{AddAssign, MulAssign, Neg, SubAssign};

use rand::Rng;

use super::errors::RingResult;

/// A polynomial in `Z_Q[X]/(X^DEGREE + 1)` at one level of a modulus chain.
///
/// `DEGREE` is a power of two. Arithmetic is modulo `Q` with `X^DEGREE = -1`:
/// `AddAssign`, `SubAssign` and `Neg` work coefficient-wise, `MulAssign` is
/// the negacyclic product. Both operands must sit at the same level.
pub trait PolyRing<const DEGREE: usize>:
    Clone
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + Neg<Output = Self>
{
    /// Ring parameters for one level (e.g. the modulus chain node).
    type Context: Clone;

    // ── Context ──────────────────────────────────────────────────────────

    /// The zero polynomial.
    fn zero(context: &Self::Context) -> Self;

    /// The context this polynomial lives in.
    fn context(&self) -> &Self::Context;

    // ── Coefficients ─────────────────────────────────────────────────────

    /// Builds a polynomial from signed coefficients, reduced modulo `Q`
    /// (`-1 -> Q - 1`). Uses the first `DEGREE` entries of `coeffs`.
    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self;

    /// Coefficients centered in `(-Q/2, Q/2]`.
    ///
    /// # Panics
    /// Panics if a centered coefficient does not fit in `i64`.
    fn to_coeffs(&self) -> [i64; DEGREE];

    // ── Encoding boundary ────────────────────────────────────────────────

    /// Builds a polynomial from integer-valued `f64` coefficients (the
    /// encoder's rounded output), which may exceed the `i64` range.
    ///
    /// # Errors
    /// Returns [`RingError::CoefficientTooLarge`](super::RingError) if a
    /// coefficient is not finite or `|c| >= 2^(modulus_bits - 1)`.
    fn from_rounded_f64(
        coeffs: &[f64; DEGREE],
        context: &Self::Context,
    ) -> RingResult<Self>;

    /// Coefficients centered in `(-Q/2, Q/2]` as `f64` (exact up to 53 bits,
    /// rounded beyond); the decoder's input.
    fn to_centered_f64(&self) -> [f64; DEGREE];

    // ── Sampling ─────────────────────────────────────────────────────────

    /// Coefficients uniform modulo `Q`.
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self;

    /// Coefficients from a rounded Gaussian with mean 0 and `std_dev`.
    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;

    /// Exactly `hamming_weight` coefficients set to `+1` or `-1`, the rest 0.
    fn sample_ternary<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;

    // ── Automorphism ─────────────────────────────────────────────────────

    /// The ring automorphism `X -> X^exponent` (`exponent` odd).
    fn automorphism(&self, exponent: u64) -> Self;

    /// Rotates the CKKS slots by `k`: `X -> X^(5^k mod 2N)`.
    ///
    /// A negative `k` applies `X -> X^(5^|k|)` followed by the conjugation
    /// `X -> X^(2N - 1)`.
    fn rotate_slots(&self, k: i32) -> Self {
        let modulus_2n = 2 * DEGREE as u64;
        let mut exponent = 1u64;
        for _ in 0..k.unsigned_abs() {
            exponent = exponent * 5 % modulus_2n;
        }
        let rotated = self.automorphism(exponent);
        if k >= 0 {
            rotated
        } else {
            rotated.automorphism(modulus_2n - 1)
        }
    }

    // ── Modulus chain ────────────────────────────────────────────────────

    /// Position in the modulus chain; level 1 is the bottom.
    fn level(&self) -> usize;

    /// A lower bound on `log2(Q)` at the current level.
    fn modulus_bits(&self) -> u32;

    /// Bit length of the divisor [`PolyRing::rescale`] divides by.
    fn rescale_bits(&self) -> u32;

    /// Divides by the top divisor of the chain and moves one level down.
    ///
    /// # Errors
    /// Returns [`RingError::InvalidLevel`](super::RingError) at level 1.
    fn rescale(&self) -> RingResult<Self>;

    /// Moves down to `level` without dividing (reduces modulo the smaller
    /// modulus).
    ///
    /// # Errors
    /// Returns [`RingError::InvalidLevel`](super::RingError) if `level` is 0
    /// or above the current level.
    fn mod_switch_to(&self, level: usize) -> RingResult<Self>;

    // ── Gadget ───────────────────────────────────────────────────────────

    /// Number of gadget terms at this level for base `2^digit_bits`.
    fn gadget_len(context: &Self::Context, digit_bits: u32) -> usize;

    /// Splits `self` into `gadget_len` small polynomials (coefficients below
    /// `2^digit_bits`) such that, for every `t`,
    ///
    /// ```text
    /// sum_k decompose(x)[k] * t.gadget_term(k) == x * t
    /// ```
    fn decompose(&self, digit_bits: u32) -> Vec<Self>;

    /// `self` scaled by the `index`-th gadget factor; see
    /// [`PolyRing::decompose`].
    fn gadget_term(&self, index: usize, digit_bits: u32) -> Self;
}
