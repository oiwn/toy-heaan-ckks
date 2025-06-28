//! Core RnsElement type and basic operations

use super::{RnsBasis, RnsError, RnsResult};
use std::sync::Arc;

/// An integer represented in Residue Number System
#[derive(Debug, Clone)]
pub struct RnsElement {
    residues: Vec<u64>,
    basis: Arc<RnsBasis>,
}

impl RnsElement {
    pub fn new(residues: Vec<u64>, basis: Arc<RnsBasis>) -> RnsResult<Self> {
        if residues.len() != basis.prime_count() {
            return Err(RnsError::IncompatibleBases);
        }
        Ok(Self { residues, basis })
    }

    pub fn to_u64(&self) -> RnsResult<u64> {
        // Chinese Remainder Theorem reconstruction
        let sum: u128 = self
            .residues
            .iter()
            .zip(self.basis.primes())
            .map(|(&residue, &prime)| {
                let q_product = self.basis.modulus_product();
                let q_i = q_product / prime as u128;
                let inv = mod_inverse_u128(q_i, prime as u128)
                    .expect("Modular inverse should exist");
                (residue as u128) * q_i * inv
            })
            .sum();

        let result = sum % self.basis.modulus_product();

        if result > u64::MAX as u128 {
            return Err(RnsError::ModulusOverflow);
        }

        Ok(result as u64)
    }

    pub fn from_u64(value: u64, basis: Arc<RnsBasis>) -> Self {
        let residues = basis.primes().iter().map(|&prime| value % prime).collect();

        Self { residues, basis }
    }

    pub fn residues(&self) -> &[u64] {
        &self.residues
    }

    pub fn basis(&self) -> &Arc<RnsBasis> {
        &self.basis
    }

    pub fn prime_count(&self) -> usize {
        self.residues.len()
    }
}

fn mod_inverse_u128(a: u128, m: u128) -> Option<u128> {
    let (g, x, _) = extended_gcd_i128(a as i128, m as i128);
    if g != 1 {
        None
    } else {
        Some(((x % m as i128 + m as i128) % m as i128) as u128)
    }
}

fn extended_gcd_i128(a: i128, b: i128) -> (i128, i128, i128) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_gcd_i128(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rns::RnsBasisBuilder;

    #[test]
    fn test_from_u64_and_back() {
        let basis = test_basis();
        let original = 12345u64;

        let element = RnsElement::from_u64(original, basis);
        let reconstructed = element.to_u64().unwrap();

        assert_eq!(original, reconstructed);
    }

    #[test]
    fn test_zero_conversion() {
        let basis = test_basis();
        let element = RnsElement::from_u64(0, basis);

        assert_eq!(element.to_u64().unwrap(), 0);
        assert!(element.residues().iter().all(|&r| r == 0));
    }

    #[test]
    fn test_residues_are_correct() {
        let basis = test_basis();
        let value = 67890u64;
        let element = RnsElement::from_u64(value, basis.clone());

        for (&residue, &prime) in element.residues().iter().zip(basis.primes()) {
            assert_eq!(residue, value % prime);
            assert!(residue < prime);
        }
    }

    #[test]
    fn test_new_with_correct_residues() {
        let basis = test_basis();
        let residues = vec![1, 2, 3]; // Assuming 3 primes in test basis

        let element = RnsElement::new(residues.clone(), basis).unwrap();
        assert_eq!(element.residues(), &residues);
    }

    #[test]
    fn test_new_with_wrong_residue_count() {
        let basis = test_basis();
        let wrong_residues = vec![1, 2]; // Wrong count

        let result = RnsElement::new(wrong_residues, basis);
        assert!(matches!(result, Err(RnsError::IncompatibleBases)));
    }

    // Helper function for tests
    fn test_basis() -> Arc<RnsBasis> {
        RnsBasisBuilder::new(16)
            .with_prime_count(3)
            .build()
            .unwrap()
    }
}
