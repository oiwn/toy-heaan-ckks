use crate::PolyRing;

#[derive(Debug)]
pub struct Plaintext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub poly: P,
    pub scale: f64,
}

pub struct Ciphertext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub c0: P,
    pub c1: P,
    pub scale: f64,
}

/* /// Represents a ciphertext in the CKKS scheme
#[derive(Clone, Debug)]
pub struct Ciphertext<const DEGREE: usize> {
    /// c0 component of ciphertext
    pub c0: RnsPolyRing<DEGREE>,
    /// c1 component of ciphertext
    pub c1: RnsPolyRing<DEGREE>,
    /// c2 component of ciphertext (coefficient of s^2),
    /// present only after multiplication before relinearization
    pub c2: Option<RnsPolyRing<DEGREE>>,
    /// Tracks the scaling factor used
    pub scale: f64,
}

impl<const DEGREE: usize> Ciphertext<DEGREE> {
    /// Create a new ciphertext from c0, c1 components
    pub fn new(
        c0: RnsPolyRing<DEGREE>,
        c1: RnsPolyRing<DEGREE>,
        scale: f64,
    ) -> Self {
        Self {
            c0,
            c1,
            c2: None,
            scale,
        }
    }

    /// Homomorphic addition of ciphertexts
    ///
    /// Adds two ciphertexts component-wise. Both ciphertexts must have the same scale.
    /// Returns a new ciphertext containing the sum.
    pub fn add(&self, other: &Self) -> Self {
        // Check scaling factors match
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale for addition: {} != {}",
            self.scale, other.scale
        );

        // Add corresponding polynomials component-wise
        let c0 = &self.c0 + &other.c0;
        let c1 = &self.c1 + &other.c1;

        Self {
            c0,
            c1,
            c2: None, // Addition preserves degree, so no c2 component
            scale: self.scale,
        }
    }

    /// In-place homomorphic addition
    ///
    /// Adds another ciphertext to this one in-place.
    pub fn add_assign(&mut self, other: &Self) {
        assert_eq!(
            self.scale, other.scale,
            "Ciphertexts must have the same scale for addition: {} != {}",
            self.scale, other.scale
        );

        self.c0
            .add_assign_checked(&other.c0)
            .expect("Failed to add c0 components");
        self.c1
            .add_assign_checked(&other.c1)
            .expect("Failed to add c1 components");

        // Clear any c2 component since addition preserves degree
        self.c2 = None;
    }

    /// Homomorphic negation
    ///
    /// Returns a new ciphertext containing the negation of this ciphertext.
    pub fn negate(&self) -> Self {
        Self {
            c0: -&self.c0,
            c1: -&self.c1,
            c2: self.c2.as_ref().map(|c2| -c2),
            scale: self.scale,
        }
    }

    /// In-place homomorphic negation
    ///
    /// Negates this ciphertext in-place.
    pub fn negate_assign(&mut self) {
        self.c0.negate_assign();
        self.c1.negate_assign();
        if let Some(ref mut c2) = self.c2 {
            c2.negate_assign();
        }
    }

    /// Check if this ciphertext has a c2 component (degree > 1)
    pub fn has_c2(&self) -> bool {
        self.c2.is_some()
    }

    /// Get the polynomial degree (number of components - 1)
    pub fn degree(&self) -> usize {
        if self.c2.is_some() { 2 } else { 1 }
    }
}

impl<const DEGREE: usize> Add for &Ciphertext<DEGREE> {
    type Output = Ciphertext<DEGREE>;

    fn add(self, other: Self) -> Self::Output {
        self.add(other)
    }
}

impl<const DEGREE: usize> Add for Ciphertext<DEGREE> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        (&self).add(&other)
    }
}

impl<const DEGREE: usize> Neg for &Ciphertext<DEGREE> {
    type Output = Ciphertext<DEGREE>;

    fn neg(self) -> Self::Output {
        self.negate()
    }
}

impl<const DEGREE: usize> Neg for Ciphertext<DEGREE> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        (&self).negate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{RnsBasisBuilder, RnsPolyRing};
    use std::sync::Arc;

    fn create_test_basis<const DEGREE: usize>() -> Arc<crate::RnsBasis> {
        Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19, 23])
                .build()
                .unwrap(),
        )
    }

    #[test]
    fn test_ciphertext_creation() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        let c0: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_i64_slice(&[1, 2, 3, 4], basis.clone());
        let c1 = RnsPolyRing::from_i64_slice(&[5, 6, 7, 8], basis.clone());
        let scale = 1024.0;

        let ct = Ciphertext::new(c0, c1, scale);

        assert_eq!(ct.scale, scale);
        assert!(!ct.has_c2());
        assert_eq!(ct.degree(), 1);
    }

    #[test]
    fn test_ciphertext_addition() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        let c0_1: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_i64_slice(&[1, 2, 3, 4], basis.clone());
        let c1_1 = RnsPolyRing::from_i64_slice(&[5, 6, 7, 8], basis.clone());
        let ct1 = Ciphertext::new(c0_1, c1_1, 1024.0);

        let c0_2 = RnsPolyRing::from_i64_slice(&[2, 3, 4, 5], basis.clone());
        let c1_2 = RnsPolyRing::from_i64_slice(&[6, 7, 8, 9], basis.clone());
        let ct2 = Ciphertext::new(c0_2, c1_2, 1024.0);

        let result = ct1.add(ct2);

        // Check that coefficients were added correctly
        let expected_c0 = [3u64, 5, 7, 9]; // 1+2, 2+3, 3+4, 4+5
        let expected_c1 = [11u64, 13, 15, 17]; // 5+6, 6+7, 7+8, 8+9

        assert_eq!(result.c0.to_u64_coefficients(), expected_c0);
        assert_eq!(result.c1.to_u64_coefficients(), expected_c1);
        assert_eq!(result.scale, 1024.0);
        assert!(!result.has_c2());
    }

    #[test]
    fn test_ciphertext_negation() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        let c0: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_i64_slice(&[1, 2, 3, 4], basis.clone());
        let c1 = RnsPolyRing::from_i64_slice(&[5, 6, 7, 8], basis.clone());
        let ct = Ciphertext::new(c0, c1, 1024.0);

        let neg_ct = ct.negate();

        // In modular arithmetic, -x = modulus - x for positive x
        let product = basis.primes().iter().product::<u64>();
        let expected_c0 = [product - 1, product - 2, product - 3, product - 4];
        let expected_c1 = [product - 5, product - 6, product - 7, product - 8];

        assert_eq!(neg_ct.c0.to_u64_coefficients(), expected_c0);
        assert_eq!(neg_ct.c1.to_u64_coefficients(), expected_c1);
        assert_eq!(neg_ct.scale, 1024.0);
    }

    #[test]
    fn test_operator_overloads() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        let c0_1: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_i64_slice(&[1, 2, 3, 4], basis.clone());
        let c1_1 = RnsPolyRing::from_i64_slice(&[5, 6, 7, 8], basis.clone());
        let ct1 = Ciphertext::new(c0_1, c1_1, 1024.0);

        let c0_2 = RnsPolyRing::from_i64_slice(&[2, 3, 4, 5], basis.clone());
        let c1_2 = RnsPolyRing::from_i64_slice(&[6, 7, 8, 9], basis.clone());
        let ct2 = Ciphertext::new(c0_2, c1_2, 1024.0);

        // Test addition operator
        let sum = &ct1 + &ct2;
        let expected_c0 = [3u64, 5, 7, 9];
        assert_eq!(sum.c0.to_u64_coefficients(), expected_c0);

        // Test negation operator
        let neg = -&ct1;
        let product = basis.primes().iter().product::<u64>();
        let expected_neg_c0 = [product - 1, product - 2, product - 3, product - 4];
        assert_eq!(neg.c0.to_u64_coefficients(), expected_neg_c0);
    }

    #[test]
    #[should_panic(expected = "Ciphertexts must have the same scale")]
    fn test_addition_different_scales_panics() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();

        let c0_1: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_i64_slice(&[1, 2, 3, 4], basis.clone());
        let c1_1 = RnsPolyRing::from_i64_slice(&[5, 6, 7, 8], basis.clone());
        let ct1 = Ciphertext::new(c0_1, c1_1, 1024.0);

        let c0_2 = RnsPolyRing::from_i64_slice(&[2, 3, 4, 5], basis.clone());
        let c1_2 = RnsPolyRing::from_i64_slice(&[6, 7, 8, 9], basis.clone());
        let ct2 = Ciphertext::new(c0_2, c1_2, 512.0); // Different scale

        let _result = ct1.add(ct2); // Should panic
    }
} */
