use thiserror::Error;

use crate::rings::backends::rns_ntt::RnsNttError;

#[derive(Error, Debug)]
pub enum CkksError {
    #[error("Invalid parameter: {message}")]
    InvalidParameter { message: String },

    #[error("Scale mismatch: expected {expected:.2}, got {actual:.2}")]
    ScaleMismatch { expected: f64, actual: f64 },

    #[error("Degree mismatch: expected {expected}, got {actual}")]
    DegreeMismatch { expected: usize, actual: usize },

    #[error("Ring error: {0}")]
    Ring(#[from] RnsNttError),
}

pub type CkksResult<T> = Result<T, CkksError>;
