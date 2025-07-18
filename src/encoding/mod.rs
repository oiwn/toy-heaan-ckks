mod bigint;
mod fft;

pub use fft::{EncodingParams, RustFftEncoder, decode, encode};
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

pub trait Encoder<const DEGREE: usize>: Send + Sync {
    fn encode(&self, values: &[f64]) -> EncodingResult<Vec<i64>>;
    fn decode(&self, coeffs: &[i64]) -> EncodingResult<Vec<f64>>;
    fn scale(&self) -> f64;
}

#[derive(Debug, Clone)]
pub enum EncoderType {
    RustFft,
    BigInt { scale_bits: u32 },
}
