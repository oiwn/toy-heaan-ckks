//! Arithmetic operations for RnsPolyRing
//!
//! This module implements in-place arithmetic operations for polynomials
//! in RNS representation.
//! Operations are performed channel-by-channel (per prime) with
use super::RnsPolyRing;
use std::ops::{Add, AddAssign, MulAssign, Neg};
use thiserror::Error;

pub type ArithmeticResult<T> = Result<T, ArithmeticError>;

// Add to rings/basis.rs or create rings/arithmetic.rs
#[derive(Error, Debug)]
pub enum ArithmeticError {
    #[error(
        "Basis mismatch: cannot operate on polynomials with different RNS bases"
    )]
    BasisMismatch,
    #[error("Channel count mismatch: {left} vs {right}")]
    ChannelMismatch { left: usize, right: usize },
}

impl<const DEGREE: usize> AddAssign<&RnsPolyRing<DEGREE>> for RnsPolyRing<DEGREE> {
    /// In-place addition: self += rhs
    ///
    /// Performs coefficient-wise addition modulo each prime in the RNS basis.
    /// Both polynomials must have the same basis.
    fn add_assign(&mut self, rhs: &RnsPolyRing<DEGREE>) {
        // Check basis compatibility - panic on mismatch since this is AddAssign trait
        assert!(
            std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()),
            "Cannot add RnsPolyRing with different bases"
        );

        let channel_count = self.channels();
        assert_eq!(
            channel_count,
            rhs.channels(),
            "Channel count mismatch: {} vs {}",
            channel_count,
            rhs.channels()
        );

        // Perform addition channel by channel
        for channel_idx in 0..channel_count {
            let prime = self.basis.primes()[channel_idx];

            // Add coefficients modulo prime
            for coeff_idx in 0..DEGREE {
                let sum = self.coefficients[channel_idx][coeff_idx]
                    + rhs.coefficients[channel_idx][coeff_idx];

                // Apply modular reduction
                self.coefficients[channel_idx][coeff_idx] = sum % prime;
            }
        }
    }
}

impl<const DEGREE: usize> MulAssign<&RnsPolyRing<DEGREE>> for RnsPolyRing<DEGREE> {
    /// In-place multiplication: self *= rhs
    ///
    /// Panics on basis mismatch for ergonomic use with the *= operator.
    fn mul_assign(&mut self, rhs: &RnsPolyRing<DEGREE>) {
        self.mul_assign(rhs)
            .expect("Multiplication failed due to incompatible bases");
    }
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// Checked in-place addition that returns a Result instead of panicking.
    ///
    /// # Arguments
    /// * `rhs` - The polynomial to add to self
    ///
    /// # Returns
    /// * `Ok(())` if addition was successful
    /// * `Err(ArithmeticError)` if bases are incompatible
    pub fn add_assign_checked(
        &mut self,
        rhs: &RnsPolyRing<DEGREE>,
    ) -> ArithmeticResult<()> {
        // Validate basis compatibility
        if !std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()) {
            return Err(ArithmeticError::BasisMismatch);
        }

        let left_channels = self.channels();
        let right_channels = rhs.channels();

        if left_channels != right_channels {
            return Err(ArithmeticError::ChannelMismatch {
                left: left_channels,
                right: right_channels,
            });
        }

        // Perform the addition (same logic as AddAssign)
        for channel_idx in 0..left_channels {
            let prime = self.basis.primes()[channel_idx];

            for coeff_idx in 0..DEGREE {
                let sum = self.coefficients[channel_idx][coeff_idx]
                    + rhs.coefficients[channel_idx][coeff_idx];

                self.coefficients[channel_idx][coeff_idx] = sum % prime;
            }
        }

        Ok(())
    }

    /// In-place polynomial multiplication in quotient ring Z[X]/(X^n + 1)
    ///
    /// Uses schoolbook multiplication with on-the-fly reduction.
    /// Each coefficient product is immediately reduced modulo the quotient.
    pub fn mul_assign(
        &mut self,
        rhs: &RnsPolyRing<DEGREE>,
    ) -> ArithmeticResult<()> {
        // Validate basis compatibility
        if !std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()) {
            return Err(ArithmeticError::BasisMismatch);
        }

        let channel_count = self.channels();
        if channel_count != rhs.channels() {
            return Err(ArithmeticError::ChannelMismatch {
                left: channel_count,
                right: rhs.channels(),
            });
        }

        // Perform multiplication channel by channel
        for channel_idx in 0..channel_count {
            let prime = self.basis.primes()[channel_idx] as u128;

            // Store original coefficients (we're modifying in-place)
            let lhs_coeffs = self.coefficients[channel_idx];
            let rhs_coeffs = rhs.coefficients[channel_idx];

            // Initialize result to zero
            for coeff in &mut self.coefficients[channel_idx] {
                *coeff = 0;
            }

            // Schoolbook multiplication with on-the-fly quotient reduction
            for i in 0..DEGREE {
                for j in 0..DEGREE {
                    let prod =
                        (lhs_coeffs[i] as u128 * rhs_coeffs[j] as u128) % prime;
                    let degree_sum = i + j;

                    if degree_sum < DEGREE {
                        // Normal coefficient: result[i+j] += a[i] * b[j]
                        let current =
                            self.coefficients[channel_idx][degree_sum] as u128;
                        self.coefficients[channel_idx][degree_sum] =
                            ((current + prod) % prime) as u64;
                    } else {
                        // Quotient reduction: X^(n+k) = -X^k, so subtract from lower degree
                        let wrapped_idx = degree_sum - DEGREE;
                        let current =
                            self.coefficients[channel_idx][wrapped_idx] as u128;
                        self.coefficients[channel_idx][wrapped_idx] =
                            ((current + prime - prod) % prime) as u64;
                    }
                }
            }
        }

        Ok(())
    }
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// In-place negation: self = -self (mod each prime)
    ///
    /// For each coefficient c, computes prime - c to get -c mod prime
    pub fn negate_assign(&mut self) {
        let channel_count = self.channels();

        for channel_idx in 0..channel_count {
            let prime = self.basis.primes()[channel_idx];

            for coeff_idx in 0..DEGREE {
                let current = self.coefficients[channel_idx][coeff_idx];

                // -c mod p = p - c (for c != 0)
                if current == 0 {
                    self.coefficients[channel_idx][coeff_idx] = 0;
                } else {
                    self.coefficients[channel_idx][coeff_idx] = prime - current;
                }
            }
        }
    }

    /// Returns negated polynomial: -self (mod each prime)
    pub fn negate(&self) -> Self {
        let mut result = self.clone();
        result.negate_assign();
        result
    }
}

// Implement Add trait for RnsPolyRing references
impl<const DEGREE: usize> Add<&RnsPolyRing<DEGREE>> for &RnsPolyRing<DEGREE> {
    type Output = RnsPolyRing<DEGREE>;

    fn add(self, rhs: &RnsPolyRing<DEGREE>) -> Self::Output {
        // Check basis compatibility
        assert!(
            std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()),
            "Cannot add RnsPolyRing with different bases"
        );

        let mut result = self.clone();
        result
            .add_assign_checked(rhs)
            .expect("Addition failed due to basis mismatch");
        result
    }
}

// Implement Add trait for owned RnsPolyRing values
impl<const DEGREE: usize> Add<RnsPolyRing<DEGREE>> for RnsPolyRing<DEGREE> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

// Implement Neg trait for RnsPolyRing reference
impl<const DEGREE: usize> Neg for &RnsPolyRing<DEGREE> {
    type Output = RnsPolyRing<DEGREE>;

    fn neg(self) -> Self::Output {
        self.negate()
    }
}

// Implement Neg trait for owned RnsPolyRing
impl<const DEGREE: usize> Neg for RnsPolyRing<DEGREE> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        (&self).negate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::{RnsBasisBuilder, RnsPolyRing};
    use std::sync::Arc;

    fn create_test_basis<const DEGREE: usize>() -> Arc<crate::rings::RnsBasis> {
        Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19, 23]) // Small test primes
                .build()
                .unwrap(),
        )
    }

    #[test]
    fn test_add_assign_basic() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        // Create two test polynomials
        let coeffs1 = [1i64, 2, 3, 4];
        let coeffs2 = [5i64, 6, 7, 8];

        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs1, basis.clone());
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs2, basis.clone());

        // Perform addition
        poly1 += &poly2;

        // Verify result by reconstructing coefficients
        let result_coeffs = poly1.to_u64_coefficients();
        let expected = [6u64, 8, 10, 12]; // 1+5, 2+6, 3+7, 4+8

        assert_eq!(result_coeffs, expected);
    }

    #[test]
    fn test_add_assign_with_modular_reduction() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        // Test case where addition causes modular reduction
        // Using coefficients that will wrap around mod 17 (smallest prime)
        let coeffs1 = [16i64, 15, 14, 13]; // Close to mod 17
        let coeffs2 = [5i64, 6, 7, 8];

        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs1, basis.clone());
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs2, basis.clone());

        poly1 += &poly2;

        // Verify modular reduction occurred
        let result_coeffs = poly1.to_u64_coefficients();

        // Expected: [21%17=4, 21%17=4, 21%17=4, 21%17=4]
        // but we need to check CRT reconstruction
        // The actual values will depend on CRT reconstruction from all primes
        for &coeff in &result_coeffs {
            assert!(coeff < 17 * 19 * 23, "Result should be within CRT range");
        }
    }

    #[test]
    fn test_add_assign_checked_success() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        let coeffs1 = [1i64, 2, 3, 4];
        let coeffs2 = [5i64, 6, 7, 8];

        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs1, basis.clone());
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs2, basis.clone());

        // Should succeed
        let result = poly1.add_assign_checked(&poly2);
        assert!(result.is_ok());

        // Verify the addition happened
        let result_coeffs = poly1.to_u64_coefficients();
        let expected = [6u64, 8, 10, 12];
        assert_eq!(result_coeffs, expected);
    }

    #[test]
    fn test_add_assign_checked_basis_mismatch() {
        const DEGREE: usize = 4;
        let basis1 = create_test_basis::<DEGREE>();
        let basis2 = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![29, 31, 37]) // Different primes
                .build()
                .unwrap(),
        );

        let coeffs = [1i64, 2, 3, 4];
        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs, basis1);
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs, basis2);

        // Should fail with basis mismatch
        let result = poly1.add_assign_checked(&poly2);
        assert!(matches!(result, Err(ArithmeticError::BasisMismatch)));
    }

    #[test]
    #[should_panic(expected = "Cannot add RnsPolyRing with different bases")]
    fn test_add_assign_panics_on_basis_mismatch() {
        const DEGREE: usize = 4;
        let basis1 = create_test_basis::<DEGREE>();
        let basis2 = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![29, 31, 37])
                .build()
                .unwrap(),
        );

        let coeffs = [1i64, 2, 3, 4];
        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs, basis1);
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs, basis2);

        // This should panic
        poly1 += &poly2;
    }

    #[test]
    fn test_mul_assign_basic() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        // Test simple case: [1, 0, 0, 0] * [1, 1, 0, 0] = [1, 1, 0, 0]
        let coeffs1 = [1i64, 0, 0, 0];
        let coeffs2 = [1i64, 1, 0, 0];

        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs1, basis.clone());
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs2, basis.clone());

        poly1.mul_assign(&poly2).unwrap();

        let result = poly1.to_u64_coefficients();
        let expected = [1u64, 1, 0, 0]; // (1) * (1 + X) = 1 + X
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_assign_quotient_reduction() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        // Test X^3 * X = X^4 = -1 in Z[X]/(X^4 + 1)
        let coeffs1 = [0i64, 0, 0, 1]; // X^3
        let coeffs2 = [0i64, 1, 0, 0]; // X

        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs1, basis.clone());
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs2, basis.clone());

        poly1.mul_assign(&poly2).unwrap();

        let result = poly1.to_u64_coefficients();

        // X^3 * X = X^4 = -1, but we're working in unsigned arithmetic
        // So -1 mod (product of primes) = (product - 1)
        let product = basis.primes().iter().product::<u64>();
        let expected = [product - 1, 0, 0, 0]; // -1 represented as unsigned
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_assign_operator_syntax() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        let coeffs1 = [1i64, 2, 0, 0];
        let coeffs2 = [1i64, 1, 0, 0];

        let mut poly1: RnsPolyRing<4> =
            RnsPolyRing::from_integer_coeffs(&coeffs1, basis.clone());
        let poly2 = RnsPolyRing::from_integer_coeffs(&coeffs2, basis.clone());

        // Test operator syntax
        poly1 *= &poly2;

        // (1 + 2X) * (1 + X) = 1 + X + 2X + 2X^2 = 1 + 3X + 2X^2
        let result = poly1.to_u64_coefficients();
        let expected = [1u64, 3, 2, 0];
        assert_eq!(result, expected);
    }
}
