use crate::{
    PolyRescale, PolyRing, PolySampler,
    math::{gaussian_coefficients, ternary_coefficients, uniform_coefficients},
};
use rand::Rng;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use std::ops::{AddAssign, MulAssign, Neg};

/// Extended context for Kim's HEAAN multiplication algorithm
/// Supports multiple moduli: q (base), Q (extended), QQ, qQ for key switching
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NaiveContext {
    pub logq: u32, // Base modulus bits (e.g., 20)
    pub logQ: u32, // Extended modulus bits (e.g., 20)
}

impl NaiveContext {
    /// Create new context with logq and logQ
    pub fn new(logq: u32, logQ: u32) -> Self {
        // Validate bit budget fits in u64
        assert!(logq <= 32, "logq too large: {} > 32 bits", logq);
        assert!(logQ <= 32, "logQ too large: {} > 32 bits", logQ);
        assert!(
            logq + logQ <= 40,
            "qQ = {} + {} = {} > 40 bits, exceeds u64 capacity",
            logq,
            logQ,
            logq + logQ
        );

        Self { logq, logQ }
    }

    /// Base ciphertext modulus: q = 2^logq
    pub fn q(&self) -> u64 {
        1u64 << self.logq
    }

    /// Extended modulus: Q = 2^logQ  
    pub fn Q(&self) -> u64 {
        1u64 << self.logQ
    }

    /// Key switching modulus: qQ = q * Q = 2^(logq + logQ)
    pub fn qQ(&self) -> u64 {
        1u64 << (self.logq + self.logQ)
    }

    /// Double extended modulus: QQ = Q^2 (uses u128 for intermediate calculations)
    pub fn QQ(&self) -> u128 {
        let Q = self.Q() as u128;
        Q * Q
    }
}

/// Legacy single-modulus constructor for backward compatibility
impl From<u64> for NaiveContext {
    fn from(modulus: u64) -> Self {
        // Extract log2 of modulus for logq, set logQ = 0
        let logq = 64 - modulus.leading_zeros();
        Self::new(logq, 0)
    }
}

#[serde_as]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NaivePolyRing<const DEGREE: usize> {
    #[serde_as(as = "[_; DEGREE]")]
    pub coeffs: [u64; DEGREE],
    pub context: u64,
}

impl<const DEGREE: usize> NaivePolyRing<DEGREE> {
    pub fn with_modulus(modulus: u64) -> Self {
        Self {
            coeffs: [0; DEGREE],
            context: modulus,
        }
    }
}

impl<const DEGREE: usize> PolyRing<DEGREE> for NaivePolyRing<DEGREE> {
    type Context = u64; // modulus

    fn zero(context: &Self::Context) -> Self {
        Self {
            coeffs: [0; DEGREE],
            context: *context,
        }
    }

    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        assert!(
            coeffs.len() >= DEGREE,
            "Insufficient number of coefficients"
        );

        let mut poly_coeffs = [0u64; DEGREE];
        for (i, coeff) in coeffs.iter().enumerate() {
            poly_coeffs[i] = if *coeff >= 0 {
                *coeff as u64
            } else {
                (context - (-coeff as u64)) % context
            };
        }

        Self {
            coeffs: poly_coeffs,
            context: *context,
        }
    }

    fn context(&self) -> &Self::Context {
        &self.context
    }

    fn to_coeffs(&self) -> [i64; DEGREE] {
        let mut signed_coeffs = [0i64; DEGREE];
        let half = self.context / 2;

        for (i, &c) in self.coeffs.iter().enumerate() {
            signed_coeffs[i] = if c < half {
                c as i64
            } else {
                (c as i64).wrapping_sub(self.context as i64)
            };
        }

        signed_coeffs
    }
}

impl<const DEGREE: usize> AddAssign<&Self> for NaivePolyRing<DEGREE> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.context, rhs.context,
            "Cannot add polynomials with different moduli: {} vs {}",
            self.context, rhs.context
        );
        for i in 0..DEGREE {
            // self.coeffs[i] = (self.coeffs[i] + rhs.coeffs[i]) % self.modulus;
            let sum = (self.coeffs[i] as u128 + rhs.coeffs[i] as u128)
                % self.context as u128;
            self.coeffs[i] = sum.try_into().unwrap();
        }
    }
}

impl<const DEGREE: usize> MulAssign<&Self> for NaivePolyRing<DEGREE> {
    fn mul_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.context, rhs.context,
            "Cannot add polynomials with different moduli: {} vs {}",
            self.context, rhs.context
        );
        // Schoolbook multiplication with X^DEGREE + 1 reduction
        let mut result = [0u64; DEGREE];

        for i in 0..DEGREE {
            for j in 0..DEGREE {
                let coeff_pos = i + j;
                let coeff_val = (self.coeffs[i] as u128 * rhs.coeffs[j] as u128)
                    % self.context as u128;

                if coeff_pos < DEGREE {
                    result[coeff_pos] = ((result[coeff_pos] as u128 + coeff_val)
                        % self.context as u128)
                        as u64;
                } else {
                    // X^DEGREE = -1, so X^(DEGREE+k) = -X^k
                    let wrapped_pos = coeff_pos - DEGREE;
                    result[wrapped_pos] =
                        ((result[wrapped_pos] as u128 + self.context as u128
                            - coeff_val)
                            % self.context as u128) as u64;
                }
            }
        }

        self.coeffs = result.map(|x| x);
    }
}

impl<const DEGREE: usize> Neg for NaivePolyRing<DEGREE> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for coeff in &mut self.coeffs {
            *coeff = (self.context - *coeff) % self.context;
        }
        self
    }
}

impl<const DEGREE: usize> PolySampler<DEGREE> for NaivePolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self {
        Self {
            coeffs: uniform_coefficients::<DEGREE, R>(*context, rng),
            context: *context,
        }
    }

    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        Self {
            coeffs: gaussian_coefficients::<DEGREE, R>(std_dev, *context, rng),
            context: *context,
        }
    }

    fn sample_tribits<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let ternary = ternary_coefficients::<DEGREE, R>(hamming_weight, rng);
        let coeffs = ternary.map(|x| if x == -1 { context - 1 } else { x as u64 });
        Self {
            coeffs,
            context: *context,
        }
    }

    fn sample_noise<R: Rng>(
        variance: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        Self::sample_gaussian(variance.sqrt(), context, rng)
    }
}

impl<const DEGREE: usize> PolyRescale<DEGREE> for NaivePolyRing<DEGREE> {
    fn rescale_assign(&mut self, scale_factor: f64) {
        // Convert scale_factor to bit shift amount (more accurate than floating point division)
        let shift_bits = (scale_factor.log2().round() as u32).min(63);

        for coeff in &mut self.coeffs {
            *coeff = modular_right_shift(*coeff, shift_bits, self.context);
        }
    }
}

/// Performs modular-aware right shift, handling negative values correctly
/// This implements HEAAN's rightShiftAndEqual behavior
fn modular_right_shift(value: u64, shift_bits: u32, modulus: u64) -> u64 {
    if shift_bits == 0 {
        return value;
    }

    let half_modulus = modulus / 2;

    // Check if this represents a negative value in modular arithmetic
    if value > half_modulus {
        // This is a negative value: convert to signed, shift, convert back
        let negative_value = -((modulus - value) as i64);
        let shifted_negative = negative_value >> shift_bits;

        if shifted_negative >= 0 {
            shifted_negative as u64
        } else {
            // Convert back to modular representation
            (modulus as i64 + shifted_negative) as u64
        }
    } else {
        // This is a positive value: simple right shift
        value >> shift_bits
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_MODULUS: u64 = (1u64 << 50) - 27; // Same as your example

    #[test]
    fn test_basic_polynomial_multiplication() {
        println!("\n🧪 Testing basic polynomial multiplication in NaivePolyRing");

        // Test simple constant multiplication: 2048 * 3072 = 6291456
        let poly1 = NaivePolyRing::<8>::from_coeffs(
            &[2048, 0, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );
        let poly2 = NaivePolyRing::<8>::from_coeffs(
            &[3072, 0, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );

        println!("   poly1 coeffs: {:?}", &poly1.coeffs);
        println!("   poly2 coeffs: {:?}", &poly2.coeffs);

        let mut result = poly1.clone();
        result *= &poly2;

        println!("   result coeffs: {:?}", &result.coeffs);
        println!("   result[0]: {}", result.coeffs[0]);

        let expected = 2048u64 * 3072u64; // = 6291456
        println!("   expected[0]: {}", expected);

        // Check if it's correct
        assert_eq!(
            result.coeffs[0], expected,
            "Polynomial multiplication failed: got {}, expected {}",
            result.coeffs[0], expected
        );

        // Other coefficients should be zero for this simple case
        for i in 1..8 {
            assert_eq!(
                result.coeffs[i], 0,
                "Coefficient {} should be 0, got {}",
                i, result.coeffs[i]
            );
        }

        println!("   ✅ Basic constant multiplication works!");
    }

    #[test]
    fn test_polynomial_multiplication_with_higher_terms() {
        println!("\n🧪 Testing polynomial multiplication with higher degree terms");

        // Test: (1 + 2x) * (3 + 4x) = 3 + 10x + 8x²
        let poly1 = NaivePolyRing::<8>::from_coeffs(
            &[1, 2, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );
        let poly2 = NaivePolyRing::<8>::from_coeffs(
            &[3, 4, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );

        let mut result = poly1.clone();
        result *= &poly2;

        println!("   (1 + 2x) * (3 + 4x) = {:?}", &result.coeffs[0..3]);

        // Expected: 3 + 10x + 8x²
        assert_eq!(result.coeffs[0], 3, "Constant term should be 3");
        assert_eq!(result.coeffs[1], 10, "x term should be 10");
        assert_eq!(result.coeffs[2], 8, "x² term should be 8");

        println!("   ✅ Higher degree multiplication works!");
    }

    #[test]
    fn test_polynomial_multiplication_modular_reduction() {
        println!("\n🧪 Testing polynomial multiplication with modular reduction");

        // Test large coefficients that might overflow
        let large_val1 = TEST_MODULUS / 2;
        let large_val2 = 3;

        let poly1 = NaivePolyRing::<8>::from_coeffs(
            &[large_val1 as i64, 0, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );
        let poly2 = NaivePolyRing::<8>::from_coeffs(
            &[large_val2 as i64, 0, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );

        let mut result = poly1.clone();
        result *= &poly2;

        let expected = (large_val1 * large_val2) % TEST_MODULUS;

        println!("   large_val1: {}", large_val1);
        println!("   large_val2: {}", large_val2);
        println!("   result[0]: {}", result.coeffs[0]);
        println!("   expected: {}", expected);

        assert_eq!(result.coeffs[0], expected, "Modular multiplication failed");

        println!("   ✅ Modular reduction works!");
    }

    #[test]
    fn test_cyclic_polynomial_multiplication() {
        println!(
            "\n🧪 Testing cyclic polynomial multiplication (X^n + 1 reduction)"
        );

        // Test if high-degree terms wrap around correctly
        // In ring R = Z[X]/(X^8 + 1), we have X^8 = -1
        // So X^7 * X = X^8 = -1

        let poly1 = NaivePolyRing::<8>::from_coeffs(
            &[0, 0, 0, 0, 0, 0, 0, 1], // X^7
            &TEST_MODULUS,
        );
        let poly2 = NaivePolyRing::<8>::from_coeffs(
            &[0, 1, 0, 0, 0, 0, 0, 0], // X
            &TEST_MODULUS,
        );

        let mut result = poly1.clone();
        result *= &poly2;

        println!("   X^7 * X = {:?}", &result.coeffs);

        // X^7 * X = X^8 = -1 in Z[X]/(X^8 + 1)
        // So result should be [-1, 0, 0, 0, 0, 0, 0, 0]
        // In modular arithmetic: -1 ≡ modulus - 1
        let expected_neg_one = TEST_MODULUS - 1;

        assert_eq!(
            result.coeffs[0], expected_neg_one,
            "X^8 reduction failed: got {}, expected {}",
            result.coeffs[0], expected_neg_one
        );

        for i in 1..8 {
            assert_eq!(
                result.coeffs[i], 0,
                "Higher terms should be 0 after X^8 reduction"
            );
        }

        println!("   ✅ Cyclic reduction X^8 = -1 works!");
    }

    #[test]
    fn test_multiplication_associativity() {
        println!("\n🧪 Testing multiplication associativity");

        let poly1 = NaivePolyRing::<8>::from_coeffs(
            &[1, 1, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );
        let poly2 = NaivePolyRing::<8>::from_coeffs(
            &[2, 1, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );
        let poly3 = NaivePolyRing::<8>::from_coeffs(
            &[1, 2, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );

        // Test (poly1 * poly2) * poly3 = poly1 * (poly2 * poly3)
        let mut left = poly1.clone();
        left *= &poly2;
        left *= &poly3;

        let mut temp = poly2.clone();
        temp *= &poly3;
        let mut right = poly1.clone();
        right *= &temp;

        println!("   (p1 * p2) * p3 = {:?}", &left.coeffs[0..4]);
        println!("   p1 * (p2 * p3) = {:?}", &right.coeffs[0..4]);

        for i in 0..8 {
            assert_eq!(
                left.coeffs[i], right.coeffs[i],
                "Associativity failed at coefficient {}",
                i
            );
        }

        println!("   ✅ Multiplication is associative!");
    }

    #[test]
    fn test_ckks_realistic_multiplication() {
        println!("\n🧪 Testing CKKS-realistic polynomial multiplication");

        // Simulate what happens in your CKKS example
        // Two polynomials representing scaled constants
        let scale = 1024.0;
        let val1 = 2.0;
        let val2 = 3.0;

        let coeff1 = (val1 * scale) as i64; // 2048
        let coeff2 = (val2 * scale) as i64; // 3072

        let poly1 = NaivePolyRing::<8>::from_coeffs(
            &[coeff1, 0, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );
        let poly2 = NaivePolyRing::<8>::from_coeffs(
            &[coeff2, 0, 0, 0, 0, 0, 0, 0],
            &TEST_MODULUS,
        );

        println!("   Input poly1[0]: {} (represents {:.1})", coeff1, val1);
        println!("   Input poly2[0]: {} (represents {:.1})", coeff2, val2);

        let mut result = poly1.clone();
        result *= &poly2;

        println!("   Raw result[0]: {}", result.coeffs[0]);

        let expected_raw = coeff1 as u64 * coeff2 as u64; // 2048 * 3072 = 6291456
        let expected_value = (result.coeffs[0] as f64) / (scale * scale); // Should be ~6.0

        println!("   Expected raw: {}", expected_raw);
        println!("   Computed value: {:.6}", expected_value);
        println!("   Expected value: {:.1}", val1 * val2);

        assert_eq!(
            result.coeffs[0], expected_raw,
            "CKKS polynomial multiplication failed"
        );

        let error = (expected_value - (val1 * val2)).abs();
        assert!(
            error < 0.001,
            "CKKS value computation failed: error {:.6}",
            error
        );

        println!("   ✅ CKKS-style polynomial multiplication works perfectly!");
    }
}
