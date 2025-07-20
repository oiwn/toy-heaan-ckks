pub mod backends;
pub mod traits;

pub use backends::{NaivePolyRing, PolyRingU256, RnsPolyRing};
pub use traits::{PolyRing, PolySampler};

#[derive(Debug, Clone)]
pub enum BackendType {
    Naive(u64),
    BigIntU256(crypto_bigint::NonZero<crypto_bigint::U256>),
    RNS(std::sync::Arc<backends::rns::RnsBasis>),
    // RNS-NTT
}

impl Default for BackendType {
    fn default() -> Self {
        Self::Naive(741507920154517877)
    }
}
