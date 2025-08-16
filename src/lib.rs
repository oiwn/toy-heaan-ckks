#![allow(clippy::needless_range_loop)]
#![allow(clippy::len_without_is_empty)]
#![allow(clippy::manual_memcpy)]
pub mod crypto;
pub mod encoding;
pub mod keys;
pub mod math;
pub mod rings;

// Re-export only the main types users need
pub use crypto::{Ciphertext, CkksEngine, CkksError, CkksResult, Plaintext};
pub use encoding::{Encoder, EncodingParams, RustFftEncoder, decode, encode};
pub use keys::{
    PublicKey, PublicKeyError, PublicKeyParams, RelinearizationKey,
    RelinearizationKeyError, RelinearizationKeyParams, SecretKey, SecretKeyError,
    SecretKeyParams,
};
pub use rings::backends::BigIntPolyRing;
pub use rings::backends::NaivePolyRing;
pub use rings::backends::RnsPolyRing;
pub use rings::{PolyRescale, PolyRing, PolySampler};
