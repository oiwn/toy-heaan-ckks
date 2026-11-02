use thiserror::Error;

use crate::rings::RingError;

#[derive(Error, Debug)]
pub enum CkksError {
    #[error(transparent)]
    Scale(#[from] crate::ckks::ScaleError),
    #[error("Invalid parameter: {message}")]
    InvalidParameter { message: String },

    /// Operand scales differ; both fields are `log2` of the scales.
    #[error("Scale mismatch: expected 2^{expected:.4}, got 2^{actual:.4}")]
    ScaleMismatch { expected: f64, actual: f64 },

    #[error("Level mismatch: expected {expected}, got {actual}")]
    LevelMismatch { expected: usize, actual: usize },

    #[error("Degree mismatch: expected {expected}, got {actual}")]
    DegreeMismatch { expected: usize, actual: usize },

    #[error("Ring error: {0}")]
    Ring(#[from] RingError),
}

pub type CkksResult<T> = Result<T, CkksError>;
