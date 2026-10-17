use thiserror::Error;

/// Errors a ring backend reports through [`PolyRing`](super::PolyRing).
///
/// Backend-agnostic: construction errors specific to one backend (bad
/// primes, wrong degree) stay in that backend's own error type.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum RingError {
    /// A level change asked for a level the chain does not have.
    #[error("level {level} is out of range (current level {max_level})")]
    InvalidLevel { level: usize, max_level: usize },
    /// A coefficient is not finite or does not fit in `(-Q/2, Q/2)`.
    #[error(
        "encoded coefficient {value:e} does not fit the modulus \
         (|coeff| must be < 2^{max_bits})"
    )]
    CoefficientTooLarge { value: f64, max_bits: u32 },
}

pub type RingResult<T> = Result<T, RingError>;
