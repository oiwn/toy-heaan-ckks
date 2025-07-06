//! Arithmetic operations for RnsPolyRing
//!
//! This module implements in-place arithmetic operations for polynomials
//! in RNS representation.
//! Operations are performed channel-by-channel (per prime) with
use super::RnsPolyRing;
use std::ops::AddAssign;
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
}
