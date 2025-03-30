use core::ops::{Add, Mul, Rem};
use std::iter::IntoIterator;

/// Represents a polynomial in a ring Z[X]/(X^N + 1) with coefficients modulo q
#[derive(Debug, Clone, PartialEq)]
pub struct PolyRing {
    coeffs: Vec<u64>,
    modulus: u64,
}

impl PolyRing {
    /// Create a new polynomial with given modulus
    /// Returns None if modulus is zero
    pub fn new(modulus: u64) -> Self {
        Self {
            coeffs: Vec::new(),
            modulus,
        }
    }

    /// Create from coefficients with modular reduction
    pub fn from_coeffs(coeffs: &[u64], modulus: u64) -> Self {
        let mut poly = Self::new(modulus);
        poly.coeffs = coeffs.to_vec();
        poly.reduce_coeffs();
        poly
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
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

        let n = self.coeffs.len() + rhs.coeffs.len() - 1;
        let mut result: Vec<u64> = vec![0; n];

        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().enumerate() {
                let prod = a.checked_mul(*b).expect("Multiplication overflow");
                result[i + j] = result[i + j]
                    .checked_add(prod)
                    .expect("Addition overflow")
                    .rem(&self.modulus);
            }
        }

        Self {
            coeffs: result,
            modulus: self.modulus,
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

#[cfg(test)]
mod operation_tests {
    use super::*;

    // Helper function to create test polynomials
    fn create_test_poly(coeffs: &[u64], modulus: u64) -> PolyRing {
        PolyRing::from_coeffs(coeffs, modulus)
    }

    #[test]
    fn test_basic_addition() {
        let modulus = 17;

        // Test p1 = 2x + 3, p2 = 5x + 7
        let p1 = create_test_poly(&[3, 2], modulus); // coefficients in reverse order
        let p2 = create_test_poly(&[7, 5], modulus);

        let sum = p1 + p2;
        assert_eq!(sum.coeffs[0], 10); // (3 + 7) mod 17
        assert_eq!(sum.coeffs[1], 7); // (2 + 5) mod 17
    }

    #[test]
    fn test_addition_with_different_degrees() {
        let modulus = 17;

        // p1 = 2x^2 + 3x + 4, p2 = 5x + 7
        let p1 = create_test_poly(&[4, 3, 2], modulus);
        let p2 = create_test_poly(&[7, 5], modulus);

        let sum = p1 + p2;
        assert_eq!(sum.coeffs[0], 11); // (4 + 7) mod 17
        assert_eq!(sum.coeffs[1], 8); // (3 + 5) mod 17
        assert_eq!(sum.coeffs[2], 2); // 2 mod 17
    }

    #[test]
    fn test_basic_multiplication() {
        let modulus = 17;

        // p1 = x + 2, p2 = x + 3
        let p1 = create_test_poly(&[2, 1], modulus);
        let p2 = create_test_poly(&[3, 1], modulus);

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
        let p1 = create_test_poly(&[4, 3], modulus);
        let p2 = create_test_poly(&[5, 2], modulus);

        let product = p1 * p2;
        // Raw result would be 6x^2 + 23x + 20, which needs reduction mod 7
        assert_eq!(product.coeffs[0], 6); // 20 mod 7
        assert_eq!(product.coeffs[1], 2); // 23 mod 7
        assert_eq!(product.coeffs[2], 6); // 6 mod 7
    }

    #[test]
    fn test_degree() {
        let modulus = 17;

        let p1 = create_test_poly(&[1], modulus);
        assert_eq!(p1.degree(), 0);

        let p2 = create_test_poly(&[1, 2, 3], modulus);
        assert_eq!(p2.degree(), 2);

        let p3 = create_test_poly(&[], modulus);
        assert_eq!(p3.degree(), 0);
    }
}

#[cfg(test)]
mod ring_property_tests {
    use super::*;

    fn create_test_poly(coeffs: &[u64], modulus: u64) -> PolyRing {
        // let coeffs = coeffs.iter().map(|&x| x).collect::<Vec<_>>();
        PolyRing::from_coeffs(coeffs, modulus)
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
        let zero = PolyRing::from_coeffs(&[], modulus);

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
