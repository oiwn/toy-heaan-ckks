//! CKKS scheme: encoding, keys, encryption and homomorphic operations.

pub mod algorithms;
pub mod builder;
pub mod encoding;
pub mod engine;
pub mod errors;
pub mod evaluator;
pub mod keys;
pub mod ops;
pub mod types;

// Re-export the main types users need
pub use engine::CkksEngine;
pub use errors::{CkksError, CkksResult};
pub use types::{Ciphertext, Plaintext};
