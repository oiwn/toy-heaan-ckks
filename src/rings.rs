use core::ops::{Add, Mul, Rem};
use std::iter::IntoIterator;

/// Represents a polynomial in a ring Z_q[X]/(X^n + 1) where:
/// q - coefficients modulo
/// n - ring dimension, should be power of 2
///
/// In this quotient ring, polynomials have degree at most n-1,
/// and X^n = -1 (which is used during multiplication)
#[derive(Debug, Clone, PartialEq)]
pub struct PolyRing {
    coefficients: Vec<u64>,
    modulus: u64,    // q
    ring_dim: usize, // n
}

impl PolyRing {
    pub fn new_empty(modulus: u64, ring_dim: usize) -> Self {
        Self {
            coefficients: Vec::with_capacity(ring_dim),
            modulus,
            ring_dim,
        }
    }

    pub fn zero(modulus: u64, ring_dim: usize) -> Self {
        let coefficients = vec![0; ring_dim];

        Self {
            coefficients,
            modulus,
            ring_dim,
        }
    }

    pub fn from_coeffs<T>(coeffs: &[T], modulus: u64, ring_dim: usize) -> Self
    where
        T: Copy,
        i64: TryFrom<T>,
        <i64 as TryFrom<T>>::Error: std::fmt::Debug,
    {
        let mut poly = Self::zero(modulus, ring_dim);

        for (i, &coeff) in coeffs.iter().enumerate().take(ring_dim) {
            // Convert to i64, panicking if the value is out of range
            let value = i64::try_from(coeff).expect("Coefficient out of i64 range");
            // Proper modular reduction for any integer
            poly.coefficients[i] =
                ((value % modulus as i64 + modulus as i64) % modulus as i64) as u64;
        }

        poly
    }

    pub fn reduce_mod(&mut self) {
        for coeff in &mut self.coefficients {
            *coeff %= self.modulus;
        }
    }

    pub fn ring_dimension(&self) -> usize {
        self.ring_dim
    }

    pub fn polynomial_degree(&self) -> u64 {
        if self.coefficients.is_empty() {
            0
        } else {
            (self.coefficients.len() - 1) as u64
        }
    }

    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// Get the length (number of coefficients)
    pub fn len(&self) -> usize {
        self.coefficients.len()
    }

    /// Check if the polynomial has no coefficients
    pub fn is_empty(&self) -> bool {
        self.coefficients.is_empty()
    }

    pub fn reduce_coeffs(&mut self) {
        for coeff in &mut self.coefficients {
            *coeff = coeff.rem(&self.modulus);
        }
    }
}

impl Add for PolyRing {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus, "Incompatible moduli");

        let max_len = self.coefficients.len().max(rhs.coefficients.len());
        self.coefficients.resize(max_len, 0);

        for (i, rhs_coeff) in rhs.coefficients.into_iter().enumerate() {
            self.coefficients[i] = self.coefficients[i]
                .checked_add(rhs_coeff)
                .expect("Addition overflow")
                .rem(&self.modulus);
        }

        self
    }
}

impl Mul for PolyRing {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus, "Incompatible moduli");
        assert_eq!(self.ring_dim, rhs.ring_dim, "Incompatible ring degrees");

        // First do standard polynomial multiplication
        let mut result = vec![0u64; 2 * self.ring_dim];

        for (i, &a) in self.coefficients.iter().enumerate() {
            for (j, &b) in rhs.coefficients.iter().enumerate() {
                let prod = (a as u128 * b as u128) % self.modulus as u128;
                result[i + j] = ((result[i + j] as u128 + prod)
                    % (self.modulus as u128))
                    as u64;
            }
        }

        // Now reduce modulo X^n + 1
        let mut reduced = vec![0u64; self.ring_dim];

        for i in 0..result.len() {
            if i < self.ring_dim {
                reduced[i] = result[i];
            } else {
                // For X^j where j >= n, we use X^j = -X^(j-n) because X^n = -1
                let idx = i % self.ring_dim;
                // Subtract because X^n = -1, so X^(n+k) = -X^k
                if result[i] != 0 {
                    reduced[idx] =
                        (reduced[idx] + self.modulus - result[i]) % self.modulus;
                }
            }
        }

        // Create the reduced polynomial
        Self {
            coefficients: reduced,
            modulus: self.modulus,
            ring_dim: self.ring_dim,
        }
    }
}

impl<'a> IntoIterator for &'a PolyRing {
    type Item = &'a u64;
    type IntoIter = std::slice::Iter<'a, u64>;

    fn into_iter(self) -> Self::IntoIter {
        self.coefficients.iter()
    }
}

#[cfg(test)]
mod operation_tests {
    use super::*;

    #[test]
    fn test_degree() {
        let modulus = 17;

        let p = PolyRing::from_coeffs(&[1, 2, 3], modulus, 8);
        assert_eq!(p.polynomial_degree(), 7);
    }

    #[test]
    fn test_addition() {
        // Test p1 = 3 + 2x, p2 = 7 + 5x
        // coefficients in reverse order
        let p1 = PolyRing::from_coeffs(&[3, 2], 17, 4);
        let p2 = PolyRing::from_coeffs(&[7, 5], 17, 4);

        let sum = p1 + p2;
        assert_eq!(sum.coefficients[0], 10); // (3 + 7) mod 17
        assert_eq!(sum.coefficients[1], 7); // (2 + 5) mod 17
    }

    #[test]
    fn test_addition_overflow() {
        // test p1 = 5 + 4x, p2 = 4 + 3x
        let p1 = PolyRing::from_coeffs(&[5, 4], 6, 4);
        let p2 = PolyRing::from_coeffs(&[4, 3], 6, 4);

        // (5 + 4x) + (3 + 4x) = (9 + 7x) mod 6 = 3 + x
        let sum = p1 + p2;
        assert_eq!(sum.coefficients[0], 3); // 9 mod 6
        assert_eq!(sum.coefficients[1], 1); // 7 mod 6
    }

    #[test]
    fn test_addition_negative() {
        // Test p1 = 3 - 2x, p2 = -4 + 5x mod 7
        let p1 = PolyRing::from_coeffs(&[3, -2], 7, 4); // becomes [3, 5]
        let p2 = PolyRing::from_coeffs(&[-4, 5], 7, 4); // becomes [3, 5]

        let sum = p1 + p2;
        // Expected: (3 - 4) = -1 = 6, (-2 + 5) = 3
        assert_eq!(sum.coefficients[0], 6);
        assert_eq!(sum.coefficients[1], 3);
    }

    #[test]
    fn test_basic_multiplication() {
        // p1 = x + 2, p2 = x + 3
        let p1 = PolyRing::from_coeffs(&[2, 1], 17, 2);
        let p2 = PolyRing::from_coeffs(&[3, 1], 17, 2);

        let product = p1 * p2;
        // In regular polynomial multiplication: Result would be x^2 + 5x + 6
        // In Z[X]/(X^n + 1) with n=2: We use X^2 = -1, which gives
        // x^2 + 5x + 6 = -1 + 5x + 6 = 5 + 5x
        assert_eq!(product.coefficients[0], 5); // 6 - 1 = 5 mod 17
        assert_eq!(product.coefficients[1], 5); // 5 mod 17
    }

    #[test]
    fn test_quotient_ring_multiplication_with_sagemath() {
        let modulus = 100003; // Large prime to avoid overflow issues in test

        // p1 = 5 + 6x + 7x^2 + 8x^3 (coefficients in reverse order in our implementation)
        let p1 = PolyRing::from_coeffs(&[5, 6, 7, 8], modulus, 4);

        // p2 = 1 + 2x + 3x^2 + 4x^3
        let p2 = PolyRing::from_coeffs(&[1, 2, 3, 4], modulus, 4);

        // Multiply polynomials
        let result = p1 * p2;

        // Expected result from SageMath: 60*x^3 + 2*x^2 - 36*x - 56
        // In modular arithmetic, negative values are represented as (modulus - value)
        let expected_coeffs = [
            modulus - 56, // -56 mod modulus
            modulus - 36, // -36 mod modulus
            2,
            60,
        ];
        println!("{:?}", result);

        // Check the result
        assert_eq!(
            result.coefficients.len(),
            4,
            "Result should have degree = 4"
        );
        // for i in 0..4 {
        for (i, item) in expected_coeffs.iter().enumerate() {
            assert_eq!(
                result.coefficients[i], *item,
                "Coefficient at position {} doesn't match",
                i
            );
        }
    }

    #[test]
    fn test_addition_with_different_degrees() {
        // p1 = 2x^2 + 3x + 4, p2 = 5x + 7
        let p1 = PolyRing::from_coeffs(&[4, 3, 2], 17, 3);
        let p2 = PolyRing::from_coeffs(&[7, 5], 17, 2);

        let sum = p1 + p2;
        assert_eq!(sum.coefficients[0], 11); // (4 + 7) mod 17
        assert_eq!(sum.coefficients[1], 8); // (3 + 5) mod 17
        assert_eq!(sum.coefficients[2], 2); // 2 mod 17
    }

    /*

    #[test]
    fn test_multiplication_with_reduction() {
        let modulus = 7;

        // p1 = 3x + 4, p2 = 2x + 5
        let p1 = create_test_poly(&[4, 3], modulus, 8);
        let p2 = create_test_poly(&[5, 2], modulus, 8);

        let product = p1 * p2;
        // Raw result would be 6x^2 + 23x + 20, which needs reduction mod 7
        assert_eq!(product.coeffs[0], 6); // 20 mod 7
        assert_eq!(product.coeffs[1], 2); // 23 mod 7
        assert_eq!(product.coeffs[2], 6); // 6 mod 7
    }


    #[test]
    fn test_polynomial_ring_multiplication() {
        // For a ring Z[X]/(X^4 + 1)
        let modulus = 1231231237; // A large prime for testing

        // Create two polynomials
        // p1 = 1 + 2x + 3x^2 + 4x^3
        let p1 = PolyRing::from_coeffs(&[1, 2, 3, 4], modulus, 4);

        // p2 = 5 + 6x + 7x^2 + 8x^3
        let p2 = PolyRing::from_coeffs(&[5, 6, 7, 8], modulus, 4);

        // Multiply them
        let result = p1.clone() * p2.clone();

        // In the ring Z[X]/(X^4 + 1), the result should be:
        // (1 + 2x + 3x^2 + 4x^3) * (5 + 6x + 7x^2 + 8x^3)
        // = 5 + 16x + 34x^2 + 60x^3 + 61x^4 + 52x^5 + 32x^6 + 32x^7
        // After reduction with X^4 = -1:
        // = 5 + 16x + 34x^2 + 60x^3 + 61(-1) + 52(-x) + 32(-x^2) + 32(-x^3)
        // = 5 - 61 + 16x - 52x + 34x^2 - 32x^2 + 60x^3 - 32x^3
        // = (5 - 61) + (16 - 52)x + (34 - 32)x^2 + (60 - 32)x^3
        // = -56 - 36x + 2x^2 + 28x^3

        // Let's check each coefficient
        assert_eq!(result.coeffs[0], (modulus - 56) % modulus); // -56 mod q
        assert_eq!(result.coeffs[1], (modulus - 36) % modulus); // -36 mod q
        assert_eq!(result.coeffs[2], 2); // 2
        assert_eq!(result.coeffs[3], 28); // 28
    }

    #[test]
    fn test_multiplicative_identity() {
        let modulus = 17;
        let ring_degree = 8;
        let p1 = create_test_poly(&[1, 2, 3], modulus, ring_degree);

        // Create polynomial representing 1 with proper length
        let mut one_coeffs = vec![0u64; ring_degree as usize];
        one_coeffs[0] = 1;
        let one = PolyRing::from_coeffs(&one_coeffs, modulus, ring_degree);

        let prod1 = p1.clone() * one.clone();
        let prod2 = one * p1.clone();

        assert_eq!(
            prod1,
            p1.clone(),
            "Multiplying by one should not change polynomial"
        );
        assert_eq!(prod2, p1, "Multiplying by one should not change polynomial");
    } */
}

#[cfg(test)]
mod ring_property_tests {
    // use super::*;

    /* #[test]
    fn test_addition_associativity() {
        let modulus = 17;

        let p1 = create_test_poly(&[1, 2], modulus);
        let p2 = create_test_poly(&[3, 4], modulus);
        let p3 = create_test_poly(&[5, 6], modulus);

        // (p1 + p2) + p3
        let sum1 = (p1.clone() + p2.clone()) + p3.clone();

        // p1 + (p2 + p3)
        let sum2 = p1 + (p2 + p3);

        assert_eq!(sum1, sum2, "Addition should be associative");
    }

    #[test]
    fn test_addition_commutativity() {
        let modulus = 17;

        let p1 = create_test_poly(&[1, 2, 3], modulus);
        let p2 = create_test_poly(&[4, 5, 6], modulus);

        let sum1 = p1.clone() + p2.clone();
        let sum2 = p2 + p1;

        assert_eq!(sum1, sum2, "Addition should be commutative");
    }

    #[test]
    fn test_multiplication_associativity() {
        let modulus = 17;

        let p1 = create_test_poly(&[1, 2], modulus);
        let p2 = create_test_poly(&[3, 4], modulus);
        let p3 = create_test_poly(&[5, 6], modulus);

        // (p1 * p2) * p3
        let prod1 = (p1.clone() * p2.clone()) * p3.clone();

        // p1 * (p2 * p3)
        let prod2 = p1 * (p2 * p3);

        assert_eq!(prod1, prod2, "Multiplication should be associative");
    }

    #[test]
    fn test_multiplication_commutativity() {
        let modulus = 17;

        let p1 = create_test_poly(&[1, 2], modulus);
        let p2 = create_test_poly(&[3, 4], modulus);

        let prod1 = p1.clone() * p2.clone();
        let prod2 = p2 * p1;

        assert_eq!(prod1, prod2, "Multiplication should be commutative");
    }

    #[test]
    fn test_distributive_property() {
        let modulus = 17;

        let p1 = create_test_poly(&[1, 2], modulus);
        let p2 = create_test_poly(&[3, 4], modulus);
        let p3 = create_test_poly(&[5, 6], modulus);

        // p1 * (p2 + p3)
        let left = p1.clone() * (p2.clone() + p3.clone());

        // (p1 * p2) + (p1 * p3)
        let right = (p1.clone() * p2) + (p1 * p3);

        assert_eq!(
            left, right,
            "Multiplication should distribute over addition"
        );
    }

    #[test]
    fn test_additive_identity() {
        let modulus = 17;
        let p1 = create_test_poly(&[1, 2, 3], modulus);

        // Create zero polynomial
        let zero = PolyRing::from_coeffs(&[], modulus, 8);

        let sum1 = p1.clone() + zero.clone();
        let sum2 = zero + p1.clone();

        assert_eq!(sum1, p1.clone(), "Adding zero should not change polynomial");
        assert_eq!(sum2, p1, "Adding zero should not change polynomial");
    }

    #[test]
    fn test_multiplicative_identity() {
        let modulus = 17;
        let p1 = create_test_poly(&[1, 2, 3], modulus);

        // Create polynomial representing 1
        let one = create_test_poly(&[1], modulus);

        let prod1 = p1.clone() * one.clone();
        let prod2 = one * p1.clone();

        assert_eq!(
            prod1,
            p1.clone(),
            "Multiplying by one should not change polynomial"
        );
        assert_eq!(prod2, p1, "Multiplying by one should not change polynomial");
    } */
}
