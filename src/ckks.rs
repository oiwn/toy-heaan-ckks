//! CKKS scheme: encoding, keys, encryption and homomorphic operations.

pub mod builder;
pub mod encoding;
pub mod engine;
pub mod errors;
pub mod keys;
pub mod pipeline;
pub mod types;

// Re-export the main types users need
pub use engine::{CkksEngine, RnsGadgetRelinKey};
pub use errors::{CkksError, CkksResult};
pub use types::{Ciphertext, Plaintext};
