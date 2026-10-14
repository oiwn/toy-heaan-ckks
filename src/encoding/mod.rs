pub mod ckks_encoder;
pub mod special_fft;

use thiserror::Error;

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
