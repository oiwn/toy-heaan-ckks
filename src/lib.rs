pub mod crypto;
pub mod encoding;
pub mod keys;
pub mod math;
pub mod rings;

// Re-export only the main types users need
pub use crypto::{CkksEngine, CkksError, CkksResult};
pub use encoding::{
    Encoder, EncoderType, EncodingParams, RustFftEncoder, decode, encode,
};
pub use keys::{
    PublicKey, PublicKeyError, PublicKeyParams, SecretKey, SecretKeyError,
    SecretKeyParams,
};
pub use rings::backends::NaivePolyRing;
pub use rings::{BackendType, PolyRing, PolySampler};
