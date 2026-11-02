//! Positive, finite CKKS encoding scales.
//!
//! Scales are real factors, not rounded bit counts. In particular, rescaling
//! divides by the actual ring divisor. The representation has the same
//! 53-bit significand precision as the encoder's `f64` coefficients.

use thiserror::Error;

/// A scale could not be represented as a positive, finite `f64`.
#[derive(Clone, Copy, Debug, Error, PartialEq)]
#[error("scale must be positive and finite, got {value}")]
pub struct ScaleError {
    /// Rejected value, including overflow or underflow from scale arithmetic.
    pub value: f64,
}

/// A positive, finite factor interpreting encoded integers as CKKS values.
///
/// The private representation ensures that checked arithmetic cannot create
/// zero, negative, infinite or NaN scales. Equality compares the represented
/// factors exactly; use [`Self::is_compatible`] for addition's roundoff check.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Scale(f64);

impl Scale {
    /// Construct a scale from its real factor.
    ///
    /// # Errors
    /// Rejects zero, negative and non-finite factors.
    pub fn new(value: f64) -> Result<Self, ScaleError> {
        if value.is_finite() && value > 0.0 {
            Ok(Self(value))
        } else {
            Err(ScaleError { value })
        }
    }

    /// Construct a power-of-two scale. Negative exponents are permitted.
    ///
    /// # Errors
    /// Rejects exponents whose result overflows or underflows `f64`.
    pub fn power_of_two(exponent: i32) -> Result<Self, ScaleError> {
        let value = match exponent {
            -1074..=-1023 => f64::from_bits(1u64 << (exponent + 1074)),
            -1022..=1023 => f64::from_bits(((exponent + 1023) as u64) << 52),
            ..=-1075 => 0.0,
            _ => f64::INFINITY,
        };
        Self::new(value)
    }

    /// The real encoding factor.
    #[must_use]
    pub fn value(self) -> f64 {
        self.0
    }

    /// Base-two logarithm for diagnostics; this is not rounded metadata.
    #[must_use]
    pub fn log2(self) -> f64 {
        self.0.log2()
    }

    /// The exponent `k` if this scale is exactly `2^k`, otherwise `None`.
    #[must_use]
    pub fn exact_log2(self) -> Option<i32> {
        // Positive finite, so log2 is finite and well inside i32.
        let exponent = self.0.log2().round() as i32;
        (Self::power_of_two(exponent).ok() == Some(self)).then_some(exponent)
    }

    /// Multiply encoding factors, as for ciphertext or plaintext products.
    ///
    /// # Errors
    /// Rejects an overflowing or underflowing result.
    pub fn checked_mul(self, other: Self) -> Result<Self, ScaleError> {
        Self::new(self.0 * other.0)
    }

    /// Divide encoding factors, as for rescaling by an actual modulus.
    ///
    /// # Errors
    /// Rejects an overflowing or underflowing result.
    pub fn checked_div(self, other: Self) -> Result<Self, ScaleError> {
        Self::new(self.0 / other.0)
    }

    /// Square root used when scheduling scales backwards through squaring.
    #[must_use]
    pub fn sqrt(self) -> Self {
        // Every positive finite f64 has a positive finite square root.
        Self(self.0.sqrt())
    }

    /// Whether scales differ only by floating-point arithmetic roundoff.
    ///
    /// Uses a symmetric relative tolerance of `64 * f64::EPSILON`. This
    /// does not modify either scale or align materially different scales.
    #[must_use]
    pub fn is_compatible(self, other: Self) -> bool {
        (self.0 - other.0).abs() <= 64.0 * f64::EPSILON * self.0.max(other.0)
    }
}

#[cfg(test)]
mod tests {
    use super::Scale;

    #[test]
    fn rejects_invalid_factors() {
        for value in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(Scale::new(value).is_err());
        }
        assert!(Scale::new(f64::from_bits(1)).is_ok());
        assert!(Scale::new(f64::MAX).is_ok());
    }

    #[test]
    fn preserves_real_factors_and_fractional_bits() {
        let scale = Scale::new(3.5).unwrap();
        assert_eq!(scale.value(), 3.5);
        assert_eq!(scale.log2(), 3.5f64.log2());
        assert_eq!(Scale::power_of_two(-3).unwrap().value(), 0.125);
        assert_eq!(
            Scale::power_of_two(-1074).unwrap().value(),
            f64::from_bits(1)
        );
        assert_eq!(
            Scale::power_of_two(-1022).unwrap().value(),
            f64::MIN_POSITIVE
        );
        assert!(Scale::power_of_two(1023).is_ok());
        assert!(Scale::power_of_two(i32::MIN).is_err());
        assert!(Scale::power_of_two(i32::MAX).is_err());
        assert!(Scale::power_of_two(1024).is_err());
        assert!(Scale::power_of_two(-1075).is_err());
    }

    #[test]
    fn rescaling_uses_actual_divisor() {
        let scale = Scale::power_of_two(30).unwrap();
        let product = scale.checked_mul(scale).unwrap();
        let divisor = Scale::new(1_073_750_017.0).unwrap();
        let rescaled = product.checked_div(divisor).unwrap();
        assert_eq!(rescaled.value(), 2.0f64.powi(60) / divisor.value());
        assert!(!rescaled.is_compatible(scale));
    }

    #[test]
    fn exact_log2_detects_powers_of_two() {
        for exponent in [-1074, -1022, -3, 0, 1, 44, 61, 1023] {
            let scale = Scale::power_of_two(exponent).unwrap();
            assert_eq!(scale.exact_log2(), Some(exponent));
        }
        for value in [3.0, 1.5 * 2f64.powi(44), 2f64.powi(44) + 1.0] {
            assert_eq!(Scale::new(value).unwrap().exact_log2(), None);
        }
    }

    #[test]
    fn checked_arithmetic_reports_overflow_and_underflow() {
        let large = Scale::new(f64::MAX).unwrap();
        let small = Scale::new(f64::from_bits(1)).unwrap();
        let two = Scale::new(2.0).unwrap();
        assert!(large.checked_mul(two).is_err());
        assert!(small.checked_div(two).is_err());
        assert!(large.checked_div(small).is_err());
        assert!(small.checked_mul(small).is_err());
    }

    #[test]
    fn compatibility_is_symmetric_and_relative() {
        for value in [1e-200, 1.0, 1e200] {
            let scale = Scale::new(value).unwrap();
            let near = Scale::new(value * (1.0 + 16.0 * f64::EPSILON)).unwrap();
            let far = Scale::new(value * (1.0 + 1e-10)).unwrap();
            assert!(scale.is_compatible(near));
            assert!(near.is_compatible(scale));
            assert!(!scale.is_compatible(far));
            assert!(!far.is_compatible(scale));
        }
    }

    #[test]
    fn backward_square_schedule_recovers_target_scale() {
        let target = Scale::new(2.0f64.powi(44)).unwrap();
        let divisor = Scale::new(1_073_750_017.0).unwrap();
        let input = target.checked_mul(divisor).unwrap().sqrt();
        let output = input
            .checked_mul(input)
            .unwrap()
            .checked_div(divisor)
            .unwrap();
        assert!(output.is_compatible(target));
    }
}
