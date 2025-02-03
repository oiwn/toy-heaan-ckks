use num_traits::{One, Zero};
use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Sub};

/// Trait for handling fixed-point arithmetic operations
pub trait FixedPoint:
    Sized
    + Clone
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Zero
    + One
{
    /// The underlying integer type used for the mantissa
    type IntegerType: Clone + Zero + PartialEq;

    /// Returns the number of bits used for fixed-point representation
    fn precision_bits() -> i64;

    /// Creates a new fixed-point number from an integer mantissa and exponent
    fn from_mantissa_exp(mantissa: Self::IntegerType, exponent: i64) -> Self;

    /// Converts a floating-point value to fixed-point representation
    fn from_float(value: f64) -> Self;

    /// Converts the fixed-point value to a floating-point number
    fn to_float(&self) -> f64;

    /// Scales the value up by a power of 2
    fn scale_up(&self, bits: i64) -> Self;

    /// Scales the value down by a power of 2
    fn scale_down(&self, bits: i64) -> Self;

    /// Returns the mantissa of the fixed-point number
    fn mantissa(&self) -> Self::IntegerType;

    /// Returns the exponent of the fixed-point number
    fn exponent(&self) -> i64;

    /// Returns the absolute value
    fn abs(&self) -> Self;

    /// Normalizes the representation to use minimal exponent
    /// while preserving precision
    fn normalize(&mut self);

    /// Rounds the fixed-point number to a specified precision
    fn round_to(&self, precision: i64) -> Self {
        let mut rounded = self.clone();
        rounded.normalize();
        match precision.cmp(&self.exponent()) {
            Ordering::Less => rounded.scale_down(self.exponent() - precision),
            Ordering::Greater => rounded.scale_up(precision - self.exponent()),
            Ordering::Equal => rounded,
        }
    }
}
