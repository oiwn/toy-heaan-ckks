use crypto_bigint::{Uint, nlimbs};
use heaan_ring_utils::PolyRing;
use proptest::prelude::*;

const LIMBS: usize = 4; // Using 4 limbs as in the code

// Strategy to generate coefficients within our modulus
fn coeff_strategy() -> impl Strategy<Value = Uint<LIMBS>> {
    // Use a fixed range for generating coefficients
    (0..1000u64).prop_map(|x| Uint::<LIMBS>::from_u64(x))
}

// Strategy to generate polynomials
fn poly_strategy(max_degree: usize) -> impl Strategy<Value = PolyRing<LIMBS>> {
    let modulus = Uint::<LIMBS>::from_u64(17);
    prop::collection::vec(coeff_strategy(), 0..=max_degree)
        .prop_map(move |coeffs| PolyRing::from_coeffs(coeffs, modulus).unwrap())
}

proptest! {
    // Basic sanity checks
    #[test]
    fn test_degree_is_valid(
        coeffs in prop::collection::vec(0u64..1000u64, 0..10),
    ) {
        let modulus = Uint::<LIMBS>::from_u64(1231231237);
        let coeffs = coeffs.into_iter()
            .map(|x| Uint::<LIMBS>::from_u64(x))
            .collect::<Vec<_>>();

        let poly = PolyRing::from_coeffs(coeffs.clone(), modulus).unwrap();
        // Degree should be coeffs.len() - 1, unless all coefficients are 0
        if coeffs.is_empty() {
            prop_assert_eq!(poly.degree(), 0);
        } else {
            prop_assert_eq!(poly.degree(), coeffs.len() - 1);
        }
    }

    // Ring properties
    #[test]
    fn test_addition_associativity(
        p1 in poly_strategy(5),
        p2 in poly_strategy(5),
        p3 in poly_strategy(5),
    ) {
        let sum1 = (p1.clone() + p2.clone()) + p3.clone();
        let sum2 = p1 + (p2 + p3);
        prop_assert_eq!(sum1, sum2);
    }

    #[test]
    fn test_addition_commutativity(
        p1 in poly_strategy(5),
        p2 in poly_strategy(5),
    ) {
        let sum1 = p1.clone() + p2.clone();
        let sum2 = p2 + p1;
        prop_assert_eq!(sum1, sum2);
    }

    #[test]
    fn test_multiplication_associativity(
        p1 in poly_strategy(3),
        p2 in poly_strategy(3),
        p3 in poly_strategy(3),
    ) {
        let prod1 = (p1.clone() * p2.clone()) * p3.clone();
        let prod2 = p1 * (p2 * p3);
        prop_assert_eq!(prod1, prod2);
    }

    #[test]
    fn test_multiplication_commutativity(
        p1 in poly_strategy(3),
        p2 in poly_strategy(3),
    ) {
        let prod1 = p1.clone() * p2.clone();
        let prod2 = p2 * p1;
        prop_assert_eq!(prod1, prod2);
    }

    #[test]
    fn test_distributive_law(
        p1 in poly_strategy(3),
        p2 in poly_strategy(3),
        p3 in poly_strategy(3),
    ) {
        let left = p1.clone() * (p2.clone() + p3.clone());
        let right = (p1.clone() * p2) + (p1 * p3);
        prop_assert_eq!(left, right);
    }

    #[test]
    fn test_additive_identity(
        p in poly_strategy(5),
    ) {
        let modulus = Uint::<LIMBS>::from_u64(17);
        let zero = PolyRing::from_coeffs(vec![], modulus).unwrap();
        prop_assert_eq!(p.clone() + zero.clone(), p.clone());
        prop_assert_eq!(zero + p.clone(), p);
    }

    #[test]
    fn test_multiplicative_identity(
        p in poly_strategy(5),
    ) {
        let modulus = Uint::<LIMBS>::from_u64(17);
        let one = PolyRing::from_coeffs(
            vec![Uint::<LIMBS>::from_u64(1)],
            modulus
        ).unwrap();
        prop_assert_eq!(p.clone() * one.clone(), p.clone());
        prop_assert_eq!(one * p.clone(), p);
    }

    // Edge cases and special values
    #[test]
    fn test_zero_polynomial_properties(
        p in poly_strategy(5),
    ) {
        let modulus = Uint::<LIMBS>::from_u64(17);
        let zero = PolyRing::from_coeffs(vec![], modulus).unwrap();
        prop_assert_eq!(p.clone() * zero.clone(), zero.clone());
        prop_assert_eq!(zero.clone() * p, zero);
    }

    #[test]
    fn test_coefficient_reduction(
        coeffs in prop::collection::vec(0u64..10000u64, 0..5),
    ) {
        let modulus = Uint::<LIMBS>::from_u64(17);
        let coeffs = coeffs.into_iter()
            .map(|x| Uint::<LIMBS>::from_u64(x))
            .collect::<Vec<_>>();

        let poly = PolyRing::from_coeffs(coeffs, modulus).unwrap();

        // Check that all coefficients are properly reduced
        for coeff in &poly {
            prop_assert!(*coeff < modulus);
        }
    }
}
