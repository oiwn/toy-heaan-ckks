use crate::fixed::FixedPoint;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use std::ops::{Add, Div, Mul, Sub};

/// Configuration constants for CKKS fixed-point arithmetic
pub trait CKKSConfig: Clone {
    /// Number of bits used for fixed-point representation
    const PRECISION: i64;
    /// Number of bits used for display/output
    const OUTPUT_PRECISION: i64 = 10;
    /// Minimum exponent allowed before overflow error
    const MIN_EXPONENT: i64 = i64::MIN / 2;
    /// Maximum exponent allowed before overflow error
    const MAX_EXPONENT: i64 = i64::MAX / 2;
}

/// Default configuration suitable for most uses
#[derive(Debug, Clone)]
pub struct DefaultConfig;

impl CKKSConfig for DefaultConfig {
    const PRECISION: i64 = 150;
}

/// Fixed-point number implementation for CKKS
#[derive(Debug, Clone)]
pub struct CKKSFixed<C: CKKSConfig> {
    /// The mantissa stored as a big integer
    mantissa: BigInt,
    /// The binary exponent
    exponent: i64,
    /// Type parameter for configuration
    _config: std::marker::PhantomData<C>,
}

impl<C: CKKSConfig> CKKSFixed<C> {
    /// Checks if exponent is within valid bounds
    fn check_bounds(&self) -> bool {
        self.exponent > C::MIN_EXPONENT && self.exponent < C::MAX_EXPONENT
    }

    pub fn new(mantissa: BigInt, exponent: i64) -> Self {
        let mut num = Self {
            mantissa,
            exponent,
            _config: std::marker::PhantomData,
        };
        num.normalize(); // Normalize during creation
        assert!(num.check_bounds(), "Exponent out of bounds");
        num
    }
}

// Implement basic arithmetic operations
impl<C: CKKSConfig> Add for CKKSFixed<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        // Align exponents before adding
        if self.exponent > rhs.exponent {
            let shift = self.exponent - rhs.exponent;
            let shifted_mantissa = self.mantissa << shift;
            Self::new(shifted_mantissa + rhs.mantissa, rhs.exponent)
        } else {
            let shift = rhs.exponent - self.exponent;
            let shifted_mantissa = rhs.mantissa << shift;
            Self::new(self.mantissa + shifted_mantissa, self.exponent)
        }
    }
}

impl<C: CKKSConfig> Sub for CKKSFixed<C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        // Similar to add but with subtraction
        if self.exponent > rhs.exponent {
            let shift = self.exponent - rhs.exponent;
            let shifted_mantissa = self.mantissa << shift;
            Self::new(shifted_mantissa - rhs.mantissa, rhs.exponent)
        } else {
            let shift = rhs.exponent - self.exponent;
            let shifted_mantissa = rhs.mantissa << shift;
            Self::new(self.mantissa - shifted_mantissa, self.exponent)
        }
    }
}

impl<C: CKKSConfig> Mul for CKKSFixed<C> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        // Multiply mantissas and add exponents
        Self::new(self.mantissa * rhs.mantissa, self.exponent + rhs.exponent)
    }
}

impl<C: CKKSConfig> Div for CKKSFixed<C> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        assert!(!rhs.is_zero(), "Division by zero");

        // Scale up for precision
        let scaled_mantissa = self.mantissa << C::PRECISION;
        Self::new(
            scaled_mantissa / rhs.mantissa,
            self.exponent - rhs.exponent - C::PRECISION,
        )
    }
}

impl<C: CKKSConfig> FixedPoint for CKKSFixed<C> {
    type IntegerType = BigInt;

    fn precision_bits() -> i64 {
        C::PRECISION
    }

    fn from_mantissa_exp(mantissa: BigInt, exponent: i64) -> Self {
        Self::new(mantissa, exponent)
    }

    fn from_float(_value: f64) -> Self {
        // TODO: Implement proper float conversion
        todo!()
    }

    fn to_float(&self) -> f64 {
        // TODO: Implement proper float conversion
        todo!()
    }

    fn scale_up(&self, bits: i64) -> Self {
        let new_exp = self.exponent + bits;
        Self::new(self.mantissa.clone(), new_exp)
    }

    fn scale_down(&self, bits: i64) -> Self {
        let new_exp = self.exponent - bits;
        Self::new(self.mantissa.clone(), new_exp)
    }

    fn mantissa(&self) -> BigInt {
        self.mantissa.clone()
    }

    fn exponent(&self) -> i64 {
        self.exponent
    }

    fn abs(&self) -> Self {
        Self::new(self.mantissa.abs(), self.exponent)
    }

    fn normalize(&mut self) {
        if self.mantissa.is_zero() {
            self.exponent = 0;
            return;
        }

        // Remove trailing zeros from mantissa
        let trailing_zeros = self.mantissa.trailing_zeros().unwrap_or(0) as i64;

        if trailing_zeros > 0 {
            self.mantissa >>= trailing_zeros;
            self.exponent += trailing_zeros;
        }

        // Check precision and scale if needed
        let mantissa_bits = self.mantissa.bits() as i64;
        if mantissa_bits > C::PRECISION {
            let excess_bits = mantissa_bits - C::PRECISION;
            self.mantissa >>= excess_bits;
            self.exponent += excess_bits;
        }
    }
}

impl<C: CKKSConfig> Zero for CKKSFixed<C> {
    fn zero() -> Self {
        Self::new(BigInt::zero(), 0)
    }

    fn is_zero(&self) -> bool {
        self.mantissa.is_zero()
    }
}

impl<C: CKKSConfig> One for CKKSFixed<C> {
    fn one() -> Self {
        Self::new(BigInt::one(), 0)
    }

    fn is_one(&self) -> bool {
        self.mantissa == BigInt::one() && self.exponent == 0
    }
}

impl<C: CKKSConfig> PartialEq for CKKSFixed<C> {
    fn eq(&self, other: &Self) -> bool {
        // First normalize both numbers
        let mut a = self.clone();
        let mut b = other.clone();
        a.normalize();
        b.normalize();

        // After normalization, equal numbers should have identical representations
        a.mantissa == b.mantissa && a.exponent == b.exponent
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_traits::{One, Zero};
    use proptest::prelude::*;

    type Fixed = CKKSFixed<DefaultConfig>;

    #[test]
    fn test_zero() {
        let zero = Fixed::zero();
        assert!(zero.is_zero());
        assert_eq!(zero.mantissa(), BigInt::from(0));
        assert_eq!(zero.exponent(), 0);
    }

    #[test]
    fn test_one() {
        let one = Fixed::one();
        assert!(one.is_one());
        assert_eq!(one.mantissa(), BigInt::from(1));
        assert_eq!(one.exponent(), 0);
    }

    #[test]
    fn test_basic_arithmetic() {
        let a = Fixed::new(BigInt::from(5), 0);
        let b = Fixed::new(BigInt::from(3), 0);

        // Addition
        let sum = a.clone() + b.clone();
        assert_eq!(sum.mantissa(), BigInt::from(8));
        assert_eq!(sum.exponent(), 0);

        // Subtraction
        let diff = a.clone() - b.clone();
        assert_eq!(diff.mantissa(), BigInt::from(2));
        assert_eq!(diff.exponent(), 0);

        // Multiplication
        let prod = a.clone() * b.clone();
        assert_eq!(prod.mantissa(), BigInt::from(15));
        assert_eq!(prod.exponent(), 0);

        // Division
        let quot = a / b;
        // Should be roughly 1.666... after scaling
        assert!(quot.mantissa() > BigInt::from(1));
        assert!(quot.mantissa() < BigInt::from(2));
    }

    #[test]
    fn test_different_exponents() {
        let a = Fixed::new(BigInt::from(5), 1); // 5 * 2^1 = 10
        let b = Fixed::new(BigInt::from(3), 0); // 3 * 2^0 = 3

        let sum = a.clone() + b.clone();
        assert_eq!(sum.mantissa(), BigInt::from(13));
        assert_eq!(sum.exponent(), 0);

        let diff = a.clone() - b.clone();
        assert_eq!(diff.mantissa(), BigInt::from(7));
        assert_eq!(diff.exponent(), 0);
    }

    #[test]
    fn test_normalization() {
        let mut a = Fixed::new(BigInt::from(8), 0);
        a.normalize();
        assert_eq!(a.mantissa(), BigInt::from(1));
        assert_eq!(a.exponent(), 3); // 8 = 1 * 2^3
    }

    #[test]
    fn test_rounding() {
        let a = Fixed::new(BigInt::from(123456), 0);
        let rounded = a.round_to(3);
        assert_eq!(rounded.mantissa(), BigInt::from(15)); // 123456 ≈ 15 * 2^13
        assert_eq!(rounded.exponent(), 13);
    }

    #[test]
    #[should_panic(expected = "Division by zero")]
    fn test_division_by_zero() {
        let a = Fixed::new(BigInt::from(5), 0);
        let b = Fixed::zero();
        let _ = a / b;
    }

    #[test]
    fn test_abs() {
        let pos = Fixed::new(BigInt::from(5), 0);
        let neg = Fixed::new(BigInt::from(-5), 0);

        assert_eq!(pos.abs().mantissa(), BigInt::from(5));
        assert_eq!(neg.abs().mantissa(), BigInt::from(5));
    }

    proptest! {
        #[test]
        fn test_add_commutative(
            a in -1000i64..1000i64,
            b in -1000i64..1000i64
        ) {
            let fa = Fixed::new(BigInt::from(a), 0);
            let fb = Fixed::new(BigInt::from(b), 0);

            let sum1 = fa.clone() + fb.clone();
            let sum2 = fb + fa;

            prop_assert_eq!(sum1.mantissa(), sum2.mantissa());
            prop_assert_eq!(sum1.exponent(), sum2.exponent());
        }

        #[test]
        fn test_mul_commutative(
            a in -1000i64..1000i64,
            b in -1000i64..1000i64
        ) {
            let fa = Fixed::new(BigInt::from(a), 0);
            let fb = Fixed::new(BigInt::from(b), 0);

            let prod1 = fa.clone() * fb.clone();
            let prod2 = fb * fa;

            prop_assert_eq!(prod1.mantissa(), prod2.mantissa());
            prop_assert_eq!(prod1.exponent(), prod2.exponent());
        }
    }
}
