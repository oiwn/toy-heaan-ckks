//! Residue Number System (RNS) implementation for CKKS
//!
//! This module provides efficient arithmetic on large integers
//! by representing them as residues modulo small primes.

mod arithmetic;
mod basis;
mod element;

pub use basis::{RnsBasis, RnsBasisBuilder};

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
