use core::ops::{Add, Mul, Rem};
use std::iter::IntoIterator;

/// Represents a polynomial in a ring Z[X]/(X^N + 1) with coefficients modulo q
/// This kind of rings called "quotient rings"
#[derive(Debug, Clone, PartialEq)]
pub struct PolyRing {
    coeffs: Vec<u64>,
    modulus: u64,     // q
    ring_degree: u64, // Ring degree (where X^n + 1 is the modulus)
}

impl PolyRing {
    /// Create a new polynomial with given modulus
    /// Returns None if modulus is zero
    pub fn new(modulus: u64, ring_degree: u64) -> Self {
        Self {
            coeffs: Vec::new(),
            modulus,
            ring_degree,
        }
    }

    /// Returns an iterator over the coefficients
    pub fn iter(&self) -> std::slice::Iter<'_, u64> {
        self.coeffs.iter()
    }

    /// Create from coefficients with modular reduction
    pub fn from_coeffs(coeffs: &[u64], modulus: u64, degree: u64) -> Self {
        let n = degree as usize;
        let mut poly = Self::new(modulus, degree);

        // Copy the coefficients
        poly.coeffs = coeffs.to_vec();

        // Pad with zeros if needed to reach the ring degree
        poly.coeffs.resize(n, 0);

        // Reduce coefficients modulo q
        poly.reduce_coeffs();

        poly
    }

    pub fn ring_degree(&self) -> u64 {
        self.ring_degree
    }

    pub fn poly_degree(&self) -> u64 {
        if self.coeffs.is_empty() {
            0
        } else {
            (self.coeffs.len() - 1) as u64
        }
    }

    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// Get the length (number of coefficients)
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Check if the polynomial has no coefficients
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn reduce_coeffs(&mut self) {
        for coeff in &mut self.coeffs {
            *coeff = coeff.rem(&self.modulus);
        }
    }
}

impl Add for PolyRing {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus, "Incompatible moduli");

        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        self.coeffs.resize(max_len, 0);

        for (i, rhs_coeff) in rhs.coeffs.into_iter().enumerate() {
            self.coeffs[i] = self.coeffs[i]
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
        assert_eq!(
            self.ring_degree, rhs.ring_degree,
            "Incompatible ring degrees"
        );

        let n = self.ring_degree as usize;

        // First do standard polynomial multiplication
        let mut result = vec![0u64; 2 * n];

        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in rhs.coeffs.iter().enumerate() {
                let prod = (a as u128 * b as u128) % self.modulus as u128;
                result[i + j] = ((result[i + j] as u128 + prod)
                    % (self.modulus as u128))
                    as u64;
            }
        }

        // Now reduce modulo X^n + 1
        let mut reduced = vec![0u64; n];

        for i in 0..result.len() {
            if i < n {
                reduced[i] = result[i];
            } else {
                // For X^j where j >= n, we use X^j = -X^(j-n) because X^n = -1
                let idx = i % n;
                // Subtract because X^n = -1, so X^(n+k) = -X^k
                if result[i] != 0 {
                    reduced[idx] =
                        (reduced[idx] + self.modulus - result[i]) % self.modulus;
                }
            }
        }

        // Create the reduced polynomial
        Self {
            coeffs: reduced,
            modulus: self.modulus,
            ring_degree: self.ring_degree,
        }
    }
}

impl<'a> IntoIterator for &'a PolyRing {
    type Item = &'a u64;
    type IntoIter = std::slice::Iter<'a, u64>;

    fn into_iter(self) -> Self::IntoIter {
        self.coeffs.iter()
    }
}

// Helper functions to convert between coefficient vectors and polynomials
pub fn coeffs_to_poly(coeffs: &[i64], modulus: u64, degree: u64) -> PolyRing {
    let u_coeffs: Vec<u64> = coeffs
        .iter()
        .map(|&c| {
            if c < 0 {
                modulus - ((-c) as u64 % modulus)
            } else {
                c as u64 % modulus
            }
        })
        .collect();

    PolyRing::from_coeffs(&u_coeffs, modulus, degree)
}

#[cfg(test)]
mod operation_tests {
    use super::*;

    // Helper function to create test polynomials
    fn create_test_poly(coeffs: &[u64], modulus: u64, degree: u64) -> PolyRing {
        PolyRing::from_coeffs(coeffs, modulus, degree)
    }

    #[test]
    fn test_basic_addition() {
        let modulus = 17;

        // Test p1 = 2x + 3, p2 = 5x + 7
        let p1 = create_test_poly(&[3, 2], modulus, 8); // coefficients in reverse order
        let p2 = create_test_poly(&[7, 5], modulus, 8);

        let sum = p1 + p2;
        assert_eq!(sum.coeffs[0], 10); // (3 + 7) mod 17
        assert_eq!(sum.coeffs[1], 7); // (2 + 5) mod 17
    }

    #[test]
    fn test_addition_with_different_degrees() {
        let modulus = 17;

        // p1 = 2x^2 + 3x + 4, p2 = 5x + 7
        let p1 = create_test_poly(&[4, 3, 2], modulus, 8);
        let p2 = create_test_poly(&[7, 5], modulus, 8);

        let sum = p1 + p2;
        assert_eq!(sum.coeffs[0], 11); // (4 + 7) mod 17
        assert_eq!(sum.coeffs[1], 8); // (3 + 5) mod 17
        assert_eq!(sum.coeffs[2], 2); // 2 mod 17
    }

    #[test]
    fn test_quotient_ring_multiplication_with_sagemath() {
        let modulus = 100003; // Large prime to avoid overflow issues in test

        // p1 = 5 + 6x + 7x^2 + 8x^3 (coefficients in reverse order in our implementation)
        let p1 = PolyRing::from_coeffs(&[5, 6, 7, 8], modulus, 8);

        // p2 = 1 + 2x + 3x^2 + 4x^3
        let p2 = PolyRing::from_coeffs(&[1, 2, 3, 4], modulus, 8);

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

        // Check the result
        assert_eq!(result.coeffs.len(), 4, "Result should have degree < 4");
        for i in 0..4 {
            assert_eq!(
                result.coeffs[i], expected_coeffs[i],
                "Coefficient at position {} doesn't match",
                i
            );
        }
    }

    #[test]
    fn test_basic_multiplication() {
        let modulus = 17;

        // p1 = x + 2, p2 = x + 3
        let p1 = create_test_poly(&[2, 1], modulus, 8);
        let p2 = create_test_poly(&[3, 1], modulus, 8);

        let product = p1 * p2;
        // In regular polynomial multiplication: Result would be x^2 + 5x + 6
        // In Z[X]/(X^n + 1) with n=2: We use X^2 = -1, which gives
        // x^2 + 5x + 6 = -1 + 5x + 6 = 5x + 5
        assert_eq!(product.coeffs[0], 5); // 6 - 1 = 5 mod 17
        assert_eq!(product.coeffs[1], 5); // 5 mod 17
    }

    #[test]
    fn test_basic_multiplication_old() {
        let modulus = 17;

        // p1 = x + 2, p2 = x + 3
        let p1 = create_test_poly(&[2, 1], modulus, 8);
        let p2 = create_test_poly(&[3, 1], modulus, 8);

        let product = p1 * p2;
        // Result should be x^2 + 5x + 6
        assert_eq!(product.coeffs[0], 6); // (2 * 3) mod 17
        assert_eq!(product.coeffs[1], 5); // (2 * 1 + 1 * 3) mod 17
        assert_eq!(product.coeffs[2], 1); // (1 * 1) mod 17
    }

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
    fn test_degree() {
        let modulus = 17;

        let p1 = create_test_poly(&[1], modulus, 8);
        assert_eq!(p1.poly_degree(), 0);

        let p2 = create_test_poly(&[1, 2, 3], modulus, 8);
        assert_eq!(p2.poly_degree(), 2);

        let p3 = create_test_poly(&[], modulus, 8);
        assert_eq!(p3.poly_degree(), 0);
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
    }
}

#[cfg(test)]
mod ring_property_tests {
    use super::*;

    fn create_test_poly(coeffs: &[u64], modulus: u64) -> PolyRing {
        // let coeffs = coeffs.iter().map(|&x| x).collect::<Vec<_>>();
        PolyRing::from_coeffs(coeffs, modulus, 8)
    }

    #[test]
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
    }
}
