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
        let count = self.prime_count.unwrap_or(3);
        let primes = generate_rns_primes(self.target_bit_size, count)?; // Add ? operator

        let basis = RnsBasis::new(primes)?;
        Ok(Arc::new(basis))
    }
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }

    // Check divisibility by numbers of form 6k ± 1
    let sqrt_n = ((n as f64).sqrt() as u64) + 1;
    let mut i = 5;
    while i <= sqrt_n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

fn generate_rns_primes(bit_size: u32, count: usize) -> RnsResult<Vec<u64>> {
    assert!(bit_size <= 63 && bit_size >= 4);

    let mut primes = Vec::with_capacity(count);
    let max_val = (1u64 << bit_size) - 1;
    let min_val = 1u64 << (bit_size - 1);
    let mut candidate = max_val - (max_val % 2 == 0) as u64; // Ensure odd

    while primes.len() < count && candidate >= min_val {
        if is_prime(candidate) {
            primes.push(candidate);
        }
        candidate = candidate.saturating_sub(2);
    }

    if primes.len() < count {
        return Err(RnsError::InvalidPrime(bit_size as u64));
    }

    Ok(primes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_prime_basic() {
        // Test known primes
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(is_prime(5));
        assert!(is_prime(7));
        assert!(is_prime(11));
        assert!(is_prime(13));
        assert!(is_prime(17));
        assert!(is_prime(19));

        // Test known composites
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(!is_prime(4));
        assert!(!is_prime(6));
        assert!(!is_prime(8));
        assert!(!is_prime(9));
        assert!(!is_prime(10));
        assert!(!is_prime(12));
        assert!(!is_prime(15));
        assert!(!is_prime(16));
    }

    #[test]
    fn test_is_prime_large() {
        // Test some larger known primes
        assert!(is_prime(65537)); // 2^16 + 1
        assert!(is_prime(982451653)); // Large prime

        // Test large composites
        assert!(!is_prime(65536)); // 2^16
        assert!(!is_prime(982451652)); // Even number
    }

    #[test]
    fn test_generate_rns_primes_small() {
        let primes = generate_rns_primes(8, 3).unwrap();

        assert_eq!(primes.len(), 3);

        // All should be prime
        for &p in &primes {
            assert!(is_prime(p), "Generated number {} is not prime", p);
        }

        // All should be within bit range
        for &p in &primes {
            assert!(p >= (1u64 << 7), "Prime {} is too small for 8-bit", p);
            assert!(p < (1u64 << 8), "Prime {} is too large for 8-bit", p);
        }

        // Should be in descending order (we generate backwards)
        for i in 1..primes.len() {
            assert!(
                primes[i - 1] > primes[i],
                "Primes should be in descending order"
            );
        }
    }

    #[test]
    fn test_generate_rns_primes_60bit() {
        let primes = generate_rns_primes(60, 2).unwrap();

        assert_eq!(primes.len(), 2);

        for &p in &primes {
            assert!(is_prime(p), "Generated number {} is not prime", p);
            assert!(p >= (1u64 << 59), "Prime {} is too small for 60-bit", p);
            assert!(p < (1u64 << 60), "Prime {} is too large for 60-bit", p);
        }
    }

    #[test]
    fn test_rns_basis_creation() {
        let primes = vec![17, 19, 23];
        let basis = RnsBasis::new(primes.clone()).unwrap();

        assert_eq!(basis.primes(), &[17, 19, 23]);
        assert_eq!(basis.prime_count(), 3);
    }

    #[test]
    fn test_rns_basis_empty_primes() {
        let result = RnsBasis::new(vec![]);
        assert!(result.is_err());

        if let Err(RnsError::InvalidPrime(0)) = result {
            // Expected error
        } else {
            panic!("Expected InvalidPrime(0) error");
        }
    }

    #[test]
    fn test_rns_basis_builder_default() {
        let basis = RnsBasisBuilder::new(16).build().unwrap();

        assert_eq!(basis.prime_count(), 3); // Default count

        // All primes should be valid
        for &p in basis.primes() {
            assert!(is_prime(p));
            assert!(p >= (1u64 << 15));
            assert!(p < (1u64 << 16));
        }
    }

    #[test]
    fn test_rns_basis_builder_custom_count() {
        let basis = RnsBasisBuilder::new(20)
            .with_prime_count(5)
            .build()
            .unwrap();

        assert_eq!(basis.prime_count(), 5);

        for &p in basis.primes() {
            assert!(is_prime(p));
            assert!(p >= (1u64 << 19));
            assert!(p < (1u64 << 20));
        }
    }

    #[test]
    fn test_rns_basis_builder_large_bit_size() {
        let basis = RnsBasisBuilder::new(32)
            .with_prime_count(2)
            .build()
            .unwrap();

        assert_eq!(basis.prime_count(), 2);

        for &p in basis.primes() {
            assert!(is_prime(p));
            assert!(p >= (1u64 << 31));
            assert!(p < (1u64 << 32));
        }

        println!("Generated 32-bit primes: {:?}", basis.primes());
    }

    #[test]
    fn test_prime_uniqueness() {
        let primes = generate_rns_primes(16, 10).unwrap();

        // Check all primes are unique
        for i in 0..primes.len() {
            for j in (i + 1)..primes.len() {
                assert_ne!(
                    primes[i], primes[j],
                    "Duplicate prime found: {}",
                    primes[i]
                );
            }
        }
    }

    #[test]
    fn test_insufficient_primes_for_bit_size() {
        // Try to generate more primes than exist in a small bit range
        let result = generate_rns_primes(4, 10); // 4-bit range has very few primes

        assert!(result.is_err());

        if let Err(RnsError::InvalidPrime(bit_size)) = result {
            assert_eq!(bit_size, 4);
        } else {
            panic!("Expected InvalidPrime error with bit_size 4");
        }
    }

    #[test]
    fn test_builder_insufficient_primes() {
        // Test that builder also returns error for impossible requests
        let result = RnsBasisBuilder::new(4).with_prime_count(10).build();

        assert!(result.is_err());
    }

    #[test]
    #[should_panic(expected = "assertion failed")]
    fn test_invalid_bit_size_too_large() {
        generate_rns_primes(64, 1).unwrap(); // Should panic - too large
    }

    #[test]
    #[should_panic(expected = "assertion failed")]
    fn test_invalid_bit_size_too_small() {
        generate_rns_primes(3, 1).unwrap(); // Should panic - too small
    }
}
