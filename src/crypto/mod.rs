//! High-level CKKS cryptographic operations
//!
//! This module provides the main user-facing API for CKKS operations,
//! including key generation, encryption, decryption, and homomorphic operations.

pub mod builder;
pub mod engine;
pub mod errors;
pub mod operations;
pub mod types;

// Re-export the main types users need
pub use engine::CkksEngine;
pub use errors::{CkksError, CkksResult};
pub use operations::{decrypt, encrypt};
pub use types::{Ciphertext, Plaintext};
