// mod fft;
pub mod ckks_encoder;
pub mod special_fft;
// pub mod textbook;

use thiserror::Error;

// use crate::{Plaintext, PolyRing};
// pub use fft::{EncodingParams, RustFftEncoder, decode, encode};
// pub use textbook::{TextbookEncoder, TextbookEncodingParams};

/// Errors raised while encoding values into a CKKS plaintext.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum EncodingError {
    /// A scaled coefficient is not finite or does not fit in `(-Q/2, Q/2)`.
    #[error(
        "encoded coefficient {value:e} does not fit the modulus \
         (|coeff| must be < 2^{max_bits})"
    )]
    CoefficientTooLarge { value: f64, max_bits: u32 },
}

pub type EncodingResult<T> = Result<T, EncodingError>;

/* pub type EncodingResult<T> = Result<T, EncodingError>;

#[derive(Error, Debug)]
pub enum EncodingError {
    #[error("Invalid input: {message}")]
    InvalidInput { message: String },

    #[error("Ring degree {degree} not supported")]
    InvalidRingDegree { degree: usize },

    #[error("Input too long: got {got}, max {max}")]
    InputTooLong { got: usize, max: usize },

    #[error("Coefficient {value} out of range")]
    CoefficientOutOfRange { value: f64 },
}

pub trait Encoder<P: PolyRing<DEGREE>, const DEGREE: usize>: Send + Sync {
    fn encode(&self, values: &[f64], context: &P::Context) -> Plaintext<P, DEGREE>;
    fn decode(&self, plaintext: &Plaintext<P, DEGREE>) -> Vec<f64>;
} */
