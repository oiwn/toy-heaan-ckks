mod fft;
pub mod special_fft;
pub mod textbook;

use crate::{Plaintext, PolyRing};
pub use fft::{EncodingParams, RustFftEncoder, decode, encode};
pub use textbook::{TextbookEncoder, TextbookEncodingParams};
use thiserror::Error;

pub type EncodingResult<T> = Result<T, EncodingError>;

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
}
