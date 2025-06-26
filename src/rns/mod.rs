//! Residue Number System (RNS) implementation for CKKS
//!
//! This module provides efficient arithmetic on large integers
//! by representing them as residues modulo small primes.

mod arithmetic;
mod basis;
mod element;

pub use arithmetic::RnsOps;
pub use basis::{RnsBasis, RnsBasisBuilder};
pub use element::RnsElement;

// Re-export common error types
pub type RnsResult<T> = Result<T, RnsError>;

#[derive(Debug, thiserror::Error)]
pub enum RnsError {
    #[error("Incompatible RNS bases")]
    IncompatibleBases,
    #[error("Invalid prime: {0}")]
    InvalidPrime(u64),
    #[error("Modulus too large for reconstruction")]
    ModulusOverflow,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_rns_operations() {
        // Create a simple RNS basis
        let basis = RnsBasisBuilder::new(60)
            .with_prime_count(5)
            .build()
            .unwrap();

        // Create some RNS elements
        let a = RnsElement::from_u64(12345, basis.clone());
        let b = RnsElement::from_u64(6789, basis.clone());

        // Test basic operations
        let sum = a.add(&b).unwrap();
        let diff = a.sub(&b).unwrap();
        let neg_a = a.neg();
        let scaled = a.mul_scalar(5);

        println!("a residues: {:?}", a.residues());
        println!("b residues: {:?}", b.residues());
        println!("sum residues: {:?}", sum.residues());

        // Basic sanity check - operations should produce valid residues
        for (&residue, &prime) in sum.residues().iter().zip(basis.primes()) {
            assert!(residue < prime);
        }
    }
}
