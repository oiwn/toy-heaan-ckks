use core::ops::{Add, Mul};
use crypto_bigint::{CheckedAdd, CheckedMul, NonZero, Uint};
use std::iter::IntoIterator;

/// Represents a polynomial in a ring Z[X]/(X^N + 1) with coefficients modulo q
#[derive(Debug, Clone, PartialEq)]
pub struct PolyRing<const LIMBS: usize> {
    coeffs: Vec<Uint<LIMBS>>,
    // Store modulus as NonZero to avoid repeated conversions
    modulus: NonZero<Uint<LIMBS>>,
}

impl<const LIMBS: usize> PolyRing<LIMBS> {
    /// Create a new polynomial with given modulus
    /// Returns None if modulus is zero
    pub fn new(modulus: Uint<LIMBS>) -> Option<Self> {
        NonZero::new(modulus)
            .map(|nz_modulus| Self {
                coeffs: Vec::new(),
                modulus: nz_modulus,
            })
            .into()
    }

    /// Create from coefficients with modular reduction
    /// Returns None if modulus is zero
    pub fn from_coeffs(
        coeffs: Vec<Uint<LIMBS>>,
        modulus: Uint<LIMBS>,
    ) -> Option<Self> {
        let mut poly = Self::new(modulus)?;
        poly.coeffs = coeffs;
        poly.reduce_coeffs();
        Some(poly)
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

impl<const LIMBS: usize> Add for PolyRing<LIMBS> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus, "Incompatible moduli");

        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        self.coeffs.resize(max_len, Uint::<LIMBS>::ZERO);

        for (i, rhs_coeff) in rhs.coeffs.into_iter().enumerate() {
            self.coeffs[i] = self.coeffs[i]
                .checked_add(&rhs_coeff)
                .expect("Addition overflow")
                .rem(&self.modulus);
        }

        self
    }
}

impl<const LIMBS: usize> Mul for PolyRing<LIMBS> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus, "Incompatible moduli");

        let n = self.coeffs.len() + rhs.coeffs.len() - 1;
        let mut result = vec![Uint::<LIMBS>::ZERO; n];

        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().enumerate() {
                let prod = a.checked_mul(b).expect("Multiplication overflow");
                result[i + j] = result[i + j]
                    .checked_add(&prod)
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

impl<'a, const LIMBS: usize> IntoIterator for &'a PolyRing<LIMBS> {
    type Item = &'a Uint<LIMBS>;
    type IntoIter = std::slice::Iter<'a, Uint<LIMBS>>;

    fn into_iter(self) -> Self::IntoIter {
        self.coeffs.iter()
    }
}

#[cfg(test)]
mod operation_tests {
    use super::*;
    use crypto_bigint::{Uint, nlimbs};

    // Helper function to create test polynomials
    fn create_test_poly(
        coeffs: &[u64],
        modulus: Uint<{ nlimbs!(256) }>,
    ) -> PolyRing<{ nlimbs!(256) }> {
        let coeffs = coeffs
            .iter()
            .map(|&x| Uint::from_u64(x))
            .collect::<Vec<_>>();
        PolyRing::from_coeffs(coeffs, modulus).unwrap()
    }

    #[test]
    fn test_basic_addition() {
        let modulus = Uint::from_u64(17);

        // Test p1 = 2x + 3, p2 = 5x + 7
        let p1 = create_test_poly(&[3, 2], modulus); // coefficients in reverse order
        let p2 = create_test_poly(&[7, 5], modulus);

        let sum = p1 + p2;
        assert_eq!(sum.coeffs[0], Uint::from_u64(10)); // (3 + 7) mod 17
        assert_eq!(sum.coeffs[1], Uint::from_u64(7)); // (2 + 5) mod 17
    }

    #[test]
    fn test_addition_with_different_degrees() {
        let modulus = Uint::from_u64(17);

        // p1 = 2x^2 + 3x + 4, p2 = 5x + 7
        let p1 = create_test_poly(&[4, 3, 2], modulus);
        let p2 = create_test_poly(&[7, 5], modulus);

        let sum = p1 + p2;
        assert_eq!(sum.coeffs[0], Uint::from_u64(11)); // (4 + 7) mod 17
        assert_eq!(sum.coeffs[1], Uint::from_u64(8)); // (3 + 5) mod 17
        assert_eq!(sum.coeffs[2], Uint::from_u64(2)); // 2 mod 17
    }

    #[test]
    fn test_basic_multiplication() {
        let modulus = Uint::from_u64(17);

        // p1 = x + 2, p2 = x + 3
        let p1 = create_test_poly(&[2, 1], modulus);
        let p2 = create_test_poly(&[3, 1], modulus);

        let product = p1 * p2;
        // Result should be x^2 + 5x + 6
        assert_eq!(product.coeffs[0], Uint::from_u64(6)); // (2 * 3) mod 17
        assert_eq!(product.coeffs[1], Uint::from_u64(5)); // (2 * 1 + 1 * 3) mod 17
        assert_eq!(product.coeffs[2], Uint::from_u64(1)); // (1 * 1) mod 17
    }

    #[test]
    fn test_multiplication_with_reduction() {
        let modulus = Uint::from_u64(7);

        // p1 = 3x + 4, p2 = 2x + 5
        let p1 = create_test_poly(&[4, 3], modulus);
        let p2 = create_test_poly(&[5, 2], modulus);

        let product = p1 * p2;
        // Raw result would be 6x^2 + 23x + 20, which needs reduction mod 7
        assert_eq!(product.coeffs[0], Uint::from_u64(6)); // 20 mod 7
        assert_eq!(product.coeffs[1], Uint::from_u64(2)); // 23 mod 7
        assert_eq!(product.coeffs[2], Uint::from_u64(6)); // 6 mod 7
    }

    #[test]
    fn test_degree() {
        let modulus = Uint::from_u64(17);

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
    use crypto_bigint::{Uint, nlimbs};

    fn create_test_poly(
        coeffs: &[u64],
        modulus: Uint<{ nlimbs!(256) }>,
    ) -> PolyRing<{ nlimbs!(256) }> {
        let coeffs = coeffs
            .iter()
            .map(|&x| Uint::from_u64(x))
            .collect::<Vec<_>>();
        PolyRing::from_coeffs(coeffs, modulus).unwrap()
    }

    #[test]
    fn test_addition_associativity() {
        let modulus = Uint::from_u64(17);

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
        let modulus = Uint::from_u64(17);

        let p1 = create_test_poly(&[1, 2, 3], modulus);
        let p2 = create_test_poly(&[4, 5, 6], modulus);

        let sum1 = p1.clone() + p2.clone();
        let sum2 = p2 + p1;

        assert_eq!(sum1, sum2, "Addition should be commutative");
    }

    #[test]
    fn test_multiplication_associativity() {
        let modulus = Uint::from_u64(17);

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
        let modulus = Uint::from_u64(17);

        let p1 = create_test_poly(&[1, 2], modulus);
        let p2 = create_test_poly(&[3, 4], modulus);

        let prod1 = p1.clone() * p2.clone();
        let prod2 = p2 * p1;

        assert_eq!(prod1, prod2, "Multiplication should be commutative");
    }

    #[test]
    fn test_distributive_property() {
        let modulus = Uint::from_u64(17);

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
        let modulus = Uint::from_u64(17);
        let p1 = create_test_poly(&[1, 2, 3], modulus);

        // Create zero polynomial
        let zero = PolyRing::from_coeffs(vec![], modulus).unwrap();

        let sum1 = p1.clone() + zero.clone();
        let sum2 = zero + p1.clone();

        assert_eq!(sum1, p1.clone(), "Adding zero should not change polynomial");
        assert_eq!(sum2, p1, "Adding zero should not change polynomial");
    }

    #[test]
    fn test_multiplicative_identity() {
        let modulus = Uint::from_u64(17);
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
