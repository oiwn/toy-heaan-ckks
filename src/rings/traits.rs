//! Traits defining the minimum interface for polynomial rings in CKKS.
//!
//! This module defines three core traits required for CKKS homomorphic encryption:
//! - [`PolyRing`]: Basic polynomial ring operations (add, multiply, negate)
//! - [`PolySampler`]: Sampling operations for key generation and encryption
//! - [`PolyRescale`]: Rescaling operation for managing ciphertext scale after multiplication
//!
//! All CKKS operations work over polynomial rings `R_q = Z_q[X]/(X^N + 1)`.
//! Concrete implementations are provided in the `backends` module.

use rand::Rng;
use std::ops::{AddAssign, MulAssign, Neg};

/// Core polynomial ring trait for CKKS operations.
///
/// Represents polynomials in `R_q = Z_q[X]/(X^DEGREE + 1)` with:
/// - Arithmetic operations modulo q (coefficient modulus)
/// - Reduction via `X^DEGREE = -1`
///
/// # Type Parameters
/// - `DEGREE`: Polynomial degree (must be a power of 2)
///
/// # Required Operations
/// - `AddAssign`: Coefficient-wise addition modulo q
/// - `MulAssign`: Polynomial multiplication with `X^DEGREE = -1` reduction
/// - `Neg`: Coefficient-wise negation modulo q
/// - `Clone`: Required for copying polynomials
pub trait PolyRing<const DEGREE: usize>:
    Clone
    + for<'a> AddAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + Neg<Output = Self>
{
    /// Context type containing ring parameters (e.g., modulus q).
    ///
    /// Different implementations may store additional parameters like
    /// precomputed constants for optimization.
    type Context;

    /// Creates the zero polynomial (all coefficients = 0).
    fn zero(context: &Self::Context) -> Self;

    /// Creates a polynomial from signed integer coefficients.
    ///
    /// # Design Note: i64 vs u64
    /// Uses `i64` for the interface to allow natural representation of signed values
    /// (e.g., -1, -2, etc.), even though coefficients are stored internally as `u64`
    /// in the range [0, q). Negative values are automatically converted to their
    /// modular equivalents: -1 → q-1, -2 → q-2, etc.
    ///
    /// # Parameters
    /// - `coeffs`: Signed coefficient array (length ≤ DEGREE)
    /// - `context`: Ring context containing modulus
    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self;

    /// Extracts polynomial coefficients as signed integers.
    ///
    /// Converts internal representation to centered representation:
    /// values in `(q/2, q)` are mapped to negative integers `(-q/2, 0)`.
    fn to_coeffs(&self) -> [i64; DEGREE];

    /// Returns a reference to the ring context.
    fn context(&self) -> &Self::Context;
}

/// Sampling operations for polynomial rings.
///
/// Provides cryptographically secure sampling methods required for CKKS:
/// - Uniform sampling for public randomness
/// - Gaussian/discrete noise for encryption security
/// - Ternary sampling for efficient secret keys
///
/// See: <https://en.wikipedia.org/wiki/Ring_learning_with_errors>
pub trait PolySampler<const DEGREE: usize>: PolyRing<DEGREE> {
    /// Samples a polynomial with uniformly random coefficients in `[0, q)`.
    ///
    /// Used for generating public randomness in encryption and key generation.
    ///
    /// See: <https://en.wikipedia.org/wiki/Discrete_uniform_distribution>
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self;

    /// Samples a polynomial with coefficients from a discrete Gaussian distribution.
    ///
    /// Each coefficient is drawn from a Gaussian with mean 0 and standard deviation
    /// `std_dev`, then rounded to the nearest integer and reduced modulo `q`.
    ///
    /// # Parameters
    /// - `std_dev`: Standard deviation of the Gaussian distribution (typically 3.2)
    ///
    /// See: <https://en.wikipedia.org/wiki/Discrete_Gaussian_sampling>
    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;

    /// Samples a sparse ternary polynomial with coefficients in `{-1, 0, 1}`.
    ///
    /// Generates a polynomial with exactly `hamming_weight` non-zero coefficients,
    /// where each non-zero coefficient is randomly +1/-1. Used for secret keys.
    ///
    /// # Parameters
    /// - `hamming_weight`: Number of non-zero coefficients (typically DEGREE/2)
    ///
    /// See: <https://en.wikipedia.org/wiki/Hamming_weight>
    fn sample_tribits<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;

    /// Samples encryption noise from a discrete Gaussian distribution.
    ///
    /// Similar to `sample_gaussian`, but parameterized by variance instead of
    /// standard deviation for consistency with CKKS literature.
    ///
    /// # Parameters
    /// - `variance`: Variance `σ^2` of the noise distribution
    /// `std_dev = sqrt(variance)`
    /// TODO: need to figure out if we actually need it....
    fn sample_noise<R: Rng>(
        variance: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;
}

pub trait PolyRescale<const DEGREE: usize> {
    fn rescale_assign(&mut self, scale_factor: f64);
}
