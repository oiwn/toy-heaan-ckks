// mod arithmetic;
// mod basis;
// mod display;
// mod ntt;
mod backends;
// mod poly_ring;

// pub use basis::{RnsBasis, RnsBasisBuilder, RnsError, RnsResult};
// pub use ntt::NttTables;
pub use backends::NaivePolyRing;

use rand::Rng;
use std::ops::{AddAssign, MulAssign, Neg};

// Core polynomial ring trait - all CKKS operations work on this
pub trait PolyRing<const DEGREE: usize>:
    Clone
    + for<'a> AddAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + Neg<Output = Self>
{
    type Context;

    fn zero() -> Self;
    fn from_coeffs(coeffs: &[u64]) -> Self;

    // Convert to coeffs if possible
    fn to_coeffs(&self) -> [u64; DEGREE];
}

// Sampling trait - provides common sampling operations for polynomials
pub trait PolySampler<const DEGREE: usize>: PolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(&self, rng: &mut R, max_coeff: u64) -> Self;
    fn sample_gaussian<R: Rng>(&self, rng: &mut R, std_dev: f64) -> Self;
    fn sample_tribits<R: Rng>(&self, rng: &mut R, hamming_weight: usize) -> Self;
    fn sample_noise<R: Rng>(&self, rng: &mut R, variance: f64) -> Self;
}
