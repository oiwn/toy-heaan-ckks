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

    fn zero(context: &Self::Context) -> Self;
    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self;
    fn to_coeffs(&self) -> [i64; DEGREE];

    fn context(&self) -> &Self::Context;
}

// Sampling trait - provides common sampling operations for polynomials
pub trait PolySampler<const DEGREE: usize>: PolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self;
    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;
    fn sample_tribits<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;
    fn sample_noise<R: Rng>(
        variance: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self;
}

// How to rescale polynomail.
pub trait PolyRescale<const DEGREE: usize> {
    fn rescale_assign(&mut self, scale_factor: f64);
}
