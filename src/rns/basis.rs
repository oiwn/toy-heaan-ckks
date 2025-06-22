//! RNS basis management and prime generation
use super::{RnsError, RnsResult};
use std::sync::Arc;

/// Collection of primes forming an RNS basis
#[derive(Debug)]
pub struct RnsBasis {
    primes: Vec<u64>,
}

impl RnsBasis {
    pub fn new(primes: Vec<u64>) -> RnsResult<Self> {
        if primes.is_empty() {
            return Err(RnsError::InvalidPrime(0));
        }
        Ok(Self { primes })
    }

    pub fn primes(&self) -> &[u64] {
        &self.primes
    }

    pub fn prime_count(&self) -> usize {
        self.primes.len()
    }
}

/// Builder for constructing RNS bases with specific properties
pub struct RnsBasisBuilder {
    target_bit_size: u32,
    prime_count: Option<usize>,
}

impl RnsBasisBuilder {
    pub fn new(target_bit_size: u32) -> Self {
        Self {
            target_bit_size,
            prime_count: None,
        }
    }

    pub fn with_prime_count(mut self, count: usize) -> Self {
        self.prime_count = Some(count);
        self
    }

    pub fn build(self) -> RnsResult<Arc<RnsBasis>> {
        // For now, just use some small primes for testing
        // TODO: Generate primes based on target_bit_size
        let primes = match self.prime_count.unwrap_or(3) {
            1 => vec![65537],
            2 => vec![65537, 65539],
            3 => vec![65537, 65539, 65551],
            n => {
                // Generate more primes as needed
                let mut primes = vec![65537, 65539, 65551];
                let mut candidate = 65557;
                while primes.len() < n {
                    if is_prime_simple(candidate) {
                        primes.push(candidate);
                    }
                    candidate += 2; // Only check odd numbers
                }
                primes
            }
        };

        let basis = RnsBasis::new(primes)?;
        Ok(Arc::new(basis))
    }
}

// Simple prime check for small numbers
fn is_prime_simple(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 {
        return true;
    }
    if n % 2 == 0 {
        return false;
    }

    let sqrt_n = (n as f64).sqrt() as u64 + 1;
    for i in (3..=sqrt_n).step_by(2) {
        if n % i == 0 {
            return false;
        }
    }
    true
}
