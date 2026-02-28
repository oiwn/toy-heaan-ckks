use thiserror::Error;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RnsNttError {
    #[error("ring degree must be a power of two, got {degree}")]
    InvalidDegree { degree: usize },
    #[error("RNS basis must contain at least one modulus")]
    EmptyBasis,
    #[error("modulus {modulus} is not NTT-friendly for degree {degree}")]
    NonNttFriendlyModulus { modulus: u64, degree: usize },
    #[error("invalid mod-drop count {drop_count} for {channel_count} channels")]
    InvalidModDrop {
        drop_count: usize,
        channel_count: usize,
    },
    #[error("channel count mismatch: expected {expected}, got {actual}")]
    ChannelCountMismatch { expected: usize, actual: usize },
    #[error("coefficient {coefficient} is not reduced modulo {modulus}")]
    NonReducedCoefficient { coefficient: u64, modulus: u64 },
}

pub type RnsNttResult<T> = Result<T, RnsNttError>;
