//! RNS basis management and prime generation
use super::RnsResult;
use std::sync::Arc;

/// Collection of primes forming an RNS basis
#[derive(Debug)]
pub struct RnsBasis {
    primes: Vec<u64>,
}

/// Builder for constructing RNS bases with specific properties
pub struct RnsBasisBuilder {
    target_bit_size: u32,
    ntt_friendly: bool,
    prime_count: Option<usize>,
}

impl RnsBasisBuilder {
    pub fn new(target_bit_size: u32) -> Self { /* */
    }
    pub fn build(self) -> RnsResult<Arc<RnsBasis>> { /* */
    }
}
