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

    /// Creates a new fixed-point number with bounds checking
    pub fn new(mantissa: BigInt, exponent: i64) -> Self {
        let mut num = Self {
            mantissa,
            exponent,
            _config: std::marker::PhantomData,
        };
        num.normalize();
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
        assert!(!FixedPoint::is_zero(&rhs), "Division by zero");

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
}
