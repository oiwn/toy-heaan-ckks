use thiserror::Error;

/// Errors raised while building an RNS-NTT basis.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RnsNttError {
    #[error("ring degree must be a power of two, got {degree}")]
    InvalidDegree { degree: usize },
    #[error("RNS basis must contain at least one modulus")]
    EmptyBasis,
    #[error("modulus {modulus} is not NTT-friendly for degree {degree}")]
    NonNttFriendlyModulus { modulus: u64, degree: usize },
    /// Shoup and Barrett arithmetic need `modulus < 2^62`.
    #[error("modulus {modulus} must be below 2^62")]
    ModulusTooLarge { modulus: u64 },
}

pub type RnsNttResult<T> = Result<T, RnsNttError>;
