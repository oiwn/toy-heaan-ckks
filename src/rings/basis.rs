use super::{NttTables, generate_primes};
use thiserror::Error;

/// Errors for RNS basis construction and operations.
#[derive(Error, Debug)]
pub enum RnsError {
    /// Requested number of primes is zero or invalid.
    #[error("Invalid prime count: {0}")]
    InvalidPrimeCount(usize),

    /// Failed to find a suitable prime for the given constraints.
    #[error("Unable to find prime for bit size {0} and ring degree {1}")]
    PrimeGenerationFailed(usize, usize),

    /// Provided basis does not match expected channel count.
    #[error(
        "RNS basis mismatch: expected {expected_count} primes, got {actual_count}"
    )]
    BasisMismatch {
        expected_count: usize,
        actual_count: usize,
    },
    #[error("Invalid prime: {0}")]
    InvalidPrime(u64),
}

/// Result type for RNS operations.
pub type RnsResult<T> = Result<T, RnsError>;

/// An RNS basis containing prime moduli, NTT twiddles, and rescale factors.
#[derive(Debug)]
pub struct RnsBasis {
    /// Prime moduli for each residue channel.
    pub(crate) primes: Vec<u64>,
    /// NTT tables
    #[allow(dead_code)]
    pub(crate) ntt_tables: NttTables,
}

impl RnsBasis {
    // Primes
    pub fn primes(&self) -> &Vec<u64> {
        &self.primes
    }

    /// Number of RNS channels (primes).
    pub fn channel_count(&self) -> usize {
        self.primes.len()
    }

    /// Validate expected prime channel count.
    pub fn validate(&self, expected: usize) -> RnsResult<()> {
        let actual = self.primes.len();
        if actual != expected {
            Err(RnsError::BasisMismatch {
                expected_count: expected,
                actual_count: actual,
            })
        } else {
            Ok(())
        }
    }
}

/// Builder for constructing an RNS basis.
pub struct RnsBasisBuilder {
    ring_degree: usize,
    prime_bits: Vec<usize>,
    custom_primes: Option<Vec<u64>>,
}

impl RnsBasisBuilder {
    /// Create a new builder for given polynomial degree.
    pub fn new(ring_degree: usize) -> Self {
        Self {
            ring_degree,
            prime_bits: Vec::new(),
            custom_primes: None,
        }
    }

    /// Specify desired bit-size for each prime modulus.
    pub fn with_prime_bits(mut self, bits: Vec<usize>) -> Self {
        self.prime_bits = bits;
        self
    }

    /// Specify custom primes directly (for testing).
    pub fn with_custom_primes(mut self, primes: Vec<u64>) -> Self {
        self.custom_primes = Some(primes);
        self
    }

    pub fn build(self) -> RnsResult<RnsBasis> {
        let primes = self.get_or_generate_primes()?;
        let ntt_tables = NttTables::build_ntt_tables_for_primes(&primes)?;

        Ok(RnsBasis { primes, ntt_tables })
    }

    // Use custom primes or generated to required bit size
    fn get_or_generate_primes(&self) -> RnsResult<Vec<u64>> {
        if let Some(ref custom) = self.custom_primes {
            if custom.is_empty() {
                return Err(RnsError::InvalidPrimeCount(0));
            }
            return Ok(custom.clone());
        }

        if self.prime_bits.is_empty() {
            return Err(RnsError::InvalidPrimeCount(0));
        }

        let mut primes = Vec::new();
        for &bit_size in &self.prime_bits {
            let new_primes = generate_primes(bit_size, 1).map_err(|_| {
                RnsError::PrimeGenerationFailed(bit_size, self.ring_degree)
            })?;
            primes.extend(new_primes);
        }

        Ok(primes)
    }

    pub fn with_uniform_prime_bits(
        mut self,
        bit_size: usize,
        count: usize,
    ) -> Self {
        self.prime_bits = std::iter::repeat(bit_size).take(count).collect();
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::math::is_prime;

    #[test]
    fn test_builder_with_custom_primes() {
        let primes = vec![17, 257, 65537];
        let basis = RnsBasisBuilder::new(1024)
            .with_custom_primes(primes.clone())
            .build()
            .expect("should build RNS basis");

        assert_eq!(basis.primes(), &primes);
        assert_eq!(basis.channel_count(), 3);
    }

    #[test]
    fn test_builder_with_uniform_prime_bits() {
        let basis = RnsBasisBuilder::new(1024)
            .with_uniform_prime_bits(10, 3)
            .build()
            .expect("should build RNS basis");

        assert_eq!(basis.channel_count(), 3);
        for &p in basis.primes() {
            assert!(
                p >= (1 << 9) && p < (1 << 10),
                "prime bit size not in expected range"
            );
            assert!(is_prime(p), "not a prime");
        }
    }

    #[test]
    fn test_builder_fails_on_empty_custom_primes() {
        let err = RnsBasisBuilder::new(1024)
            .with_custom_primes(vec![])
            .build()
            .unwrap_err();

        match err {
            RnsError::InvalidPrimeCount(0) => (),
            _ => panic!("unexpected error: {:?}", err),
        }
    }

    #[test]
    fn test_builder_fails_on_empty_prime_bits() {
        let err = RnsBasisBuilder::new(1024).build().unwrap_err();

        match err {
            RnsError::InvalidPrimeCount(0) => (),
            _ => panic!("unexpected error: {:?}", err),
        }
    }

    #[test]
    fn test_validate_channel_count_ok() {
        let basis = RnsBasisBuilder::new(1024)
            .with_uniform_prime_bits(12, 2)
            .build()
            .unwrap();

        assert!(basis.validate(2).is_ok());
    }

    #[test]
    fn test_validate_channel_count_fail() {
        let basis = RnsBasisBuilder::new(1024)
            .with_uniform_prime_bits(12, 2)
            .build()
            .unwrap();

        let err = basis.validate(3).unwrap_err();
        match err {
            RnsError::BasisMismatch {
                expected_count: 3,
                actual_count: 2,
            } => (),
            _ => panic!("unexpected error: {:?}", err),
        }
    }
}
