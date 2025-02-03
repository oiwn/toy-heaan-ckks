use crate::fixed::FixedPoint;

use num_bigint::{BigInt, Sign};
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::cell::Cell;
use std::ops::{Add, Div, Mul, Sub};

thread_local! {
    /// Number of bits used for fixed-point representation
    static PRECISION: Cell<i64> = Cell::new(150);
    static SCALE_FACTOR: Cell<i64> = Cell::new(30);
    static MIN_PRECISION: Cell<i64> = Cell::new(20);
    /// Number of bits used for display/output
    static OUTPUT_PRECISION: Cell<i64> = Cell::new(10);
}

/// Minimum exponent allowed before overflow error
const MIN_EXPONENT: i64 = i64::MIN / 2;
/// Maximum exponent allowed before overflow error
const MAX_EXPONENT: i64 = i64::MAX / 2;

/// Default configuration suitable for most uses
#[derive(Debug, Clone)]
pub struct DefaultConfig;

/// Fixed-point number implementation for CKKS
#[derive(Debug, Clone)]
pub struct CKKSFixed {
    mantissa: BigInt,
    exponent: i64,
}

impl CKKSFixed {
    /// Checks if exponent is within valid bounds
    fn check_bounds(&self) -> bool {
        self.exponent > MIN_EXPONENT && self.exponent < MAX_EXPONENT
    }

    pub fn new(mantissa: BigInt, exponent: i64) -> Self {
        let mut num = Self { mantissa, exponent };
        num.normalize(); // Normalize during creation
        assert!(num.check_bounds(), "Exponent out of bounds");
        num
    }

    pub fn rescale(&mut self, target_precision: i64) {
        let current_bits = self.mantissa.bits() as i64;
        if current_bits > target_precision {
            let shift = current_bits - target_precision;
            self.mantissa >>= shift;
            self.exponent += shift;
        }
    }

    /// Convert i64 to fixed-point
    pub fn from_i64(n: i64) -> Self {
        let scale = SCALE_FACTOR.with(|s| s.get());
        CKKSFixed {
            // left-shift is like multiply 2^scale
            mantissa: BigInt::from(n) << scale,
            // nagative exponent to divide when decode
            exponent: -scale,
        }
    }

    /// Convert fixed-point to f64
    pub fn to_f64(&self) -> f64 {
        // calculate as mantissa * 2^(exponent)
        let mantissa_f64 = self.mantissa.to_f64().unwrap();
        mantissa_f64 * 2f64.powi(self.exponent as i32)
    }
}

// Implement basic arithmetic operations
impl Add for CKKSFixed {
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

impl Sub for CKKSFixed {
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

impl Mul for CKKSFixed {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let precision = PRECISION.with(|p| p.get());
        // Multiply mantissas
        let mut mantissa = self.mantissa * rhs.mantissa;
        // Adjust precision
        let bits = mantissa.bits() as i64;
        if bits > precision {
            mantissa >>= bits - precision;
        }
        Self::new(mantissa, self.exponent + rhs.exponent)
    }
}

impl Div for CKKSFixed {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        assert!(!rhs.is_zero(), "Division by zero");

        // Scale up for precision
        let scaled_mantissa = self.mantissa << PRECISION.with(|p| p.get());
        Self::new(
            scaled_mantissa / rhs.mantissa,
            self.exponent - rhs.exponent - PRECISION.with(|p| p.get()),
        )
    }
}

impl FixedPoint for CKKSFixed {
    type IntegerType = BigInt;

    fn precision_bits() -> i64 {
        PRECISION.with(|p| p.get())
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

        let mantissa_bits = self.mantissa.bits() as i64;

        // Only scale down if we exceed precision
        if mantissa_bits > PRECISION.with(|p| p.get()) {
            let scale_down = mantissa_bits - PRECISION.with(|p| p.get());
            self.mantissa >>= scale_down;
            self.exponent += scale_down;
        }
    }
}

impl Zero for CKKSFixed {
    fn zero() -> Self {
        Self::new(BigInt::zero(), 0)
    }

    fn is_zero(&self) -> bool {
        self.mantissa.is_zero()
    }
}

impl One for CKKSFixed {
    fn one() -> Self {
        Self::new(BigInt::one(), 0)
    }

    fn is_one(&self) -> bool {
        self.mantissa == BigInt::one() && self.exponent == 0
    }
}

impl PartialEq for CKKSFixed {
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
    // use proptest::prelude::*;

    #[test]
    fn test_zero() {
        let zero = CKKSFixed::zero();
        assert!(zero.is_zero());
        assert_eq!(zero.mantissa(), BigInt::from(0));
        assert_eq!(zero.exponent(), 0);
    }

    #[test]
    fn test_one() {
        let one = CKKSFixed::one();
        assert!(one.is_one());
        assert_eq!(one.mantissa(), BigInt::from(1));
        assert_eq!(one.exponent(), 0);
    }

    #[test]
    fn test_special_values() {
        // Test zero
        let zero = CKKSFixed::from_i64(0);
        assert!(zero.is_zero());
        assert_eq!(zero.to_f64(), 0.0);

        // Test one
        let one = CKKSFixed::from_i64(1);
        assert_eq!(one.to_f64(), 1.0);

        // Test i64 max
        let large = CKKSFixed::from_i64(i64::MAX);
        assert_eq!(large.to_f64(), i64::MAX.to_f64().unwrap());

        // Test i64 min
        let small = CKKSFixed::from_i64(i64::MIN);
        assert_eq!(small.to_f64(), i64::MIN.to_f64().unwrap());
    }

    /* #[test]
    fn test_float_conversion() {
        // Test basic numbers
        let test_values = vec![
            1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 0.333333, -0.333333, 1234.5678, -1234.5678,
        ];

        for val in test_values {
            let fixed = CKKSFixed::from_float(val);
            let back_to_float = fixed.to_float();

            // Print debug info
            println!("Original: {}", val);
            println!(
                "As fixed: mantissa={}, exp={}",
                fixed.mantissa(),
                fixed.exponent()
            );
            println!("Back to float: {}", back_to_float);
            println!(
                "Relative error: {}",
                (val - back_to_float).abs() / val.abs()
            );

            // Check relative error is small
            // We use a relatively large epsilon because we expect some precision loss
            assert!(
                (val - back_to_float).abs() / val.abs() < 1e-10,
                "Failed conversion for {}",
                val
            );
        }
    } */

    /* #[test]
    fn test_precision_maintenance() {
        // Set a specific precision for testing
        PRECISION.with(|p| p.set(53)); // Double precision

        let original = 1.0 / 3.0; // A number that can't be exactly represented
        let fixed = CKKSFixed::from_float(original);
        let reconstructed = fixed.to_float();

        println!("Original: {}", original);
        println!("Mantissa: {}", fixed.mantissa());
        println!("Exponent: {}", fixed.exponent());
        println!("Reconstructed: {}", reconstructed);
        println!("Absolute error: {}", (original - reconstructed).abs());
        println!(
            "Relative error: {}",
            (original - reconstructed).abs() / original
        );

        // Should maintain precision at least to 1e-10
        assert!((original - reconstructed).abs() / original < 1e-10);
    }

    #[test]
    fn test_normalization() {
        // Create number 8 = 2³
        let mut a = CKKSFixed {
            mantissa: BigInt::from(8),
            exponent: 0,
        };

        // Before normalization: 8 * 2⁰
        assert_eq!(a.mantissa, BigInt::from(8));
        assert_eq!(a.exponent, 0);

        a.normalize();

        // After normalization should still be 8 * 2⁰ since it's under precision limit
        assert_eq!(a.mantissa, BigInt::from(8));
        assert_eq!(a.exponent, 0);

        // Now let's test with a number exceeding precision
        let big_num = BigInt::from(1) << (PRECISION.with(|p| p.get()) + 10); // Exceeds precision by 10 bits
        let mut b = CKKSFixed {
            mantissa: big_num,
            exponent: 0,
        };

        b.normalize();

        // Should be scaled down by 10 bits
        assert_eq!(b.mantissa.bits() as i64, PRECISION.with(|p| p.get()));
        assert_eq!(b.exponent, 10);
    } */

    /* #[test]
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
    } */
}
