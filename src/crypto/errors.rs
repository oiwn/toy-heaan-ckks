use thiserror::Error;

#[derive(Error, Debug)]
pub enum CkksError {
    #[error("Invalid parameter: {message}")]
    InvalidParameter { message: String },

    // #[error("Backend error: {source}")]
    // BackendError {
    //     #[from]
    //     source: crate::rings::RingError,
    // },

    // #[error("Key operation failed: {source}")]
    // KeyError {
    //     #[from]
    //     source: crate::keys::KeyError,
    // },
    #[error("Encoding failed: {source}")]
    EncodingError {
        #[from]
        source: crate::encoding::EncodingError,
    },

    #[error("Scale mismatch: expected {expected:.2}, got {actual:.2}")]
    ScaleMismatch { expected: f64, actual: f64 },

    #[error("Degree mismatch: expected {expected}, got {actual}")]
    DegreeMismatch { expected: usize, actual: usize },
}

pub type CkksResult<T> = Result<T, CkksError>;
