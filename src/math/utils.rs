pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_gcd(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

pub fn mod_inverse(a: u64, m: u64) -> u64 {
    let (g, x, _) = extended_gcd(a as i64, m as i64);
    assert_eq!(g, 1, "Numbers must be coprime");
    ((x % m as i64 + m as i64) % m as i64) as u64
}

pub fn crt_reconstruct(residues: &[u64], primes: &[u64]) -> u64 {
    let product: u64 = primes.iter().product();

    let mut result = 0u64;
    for (_i, (&residue, &prime)) in residues.iter().zip(primes.iter()).enumerate() {
        let partial_product = product / prime;
        let inverse = mod_inverse(partial_product, prime);
        result = (result + residue * partial_product * inverse) % product;
    }

    result
}

/// Checks whether a given number `n` is prime using trial division.
///
/// This function uses the 6k ± 1 optimization to skip unnecessary checks:
/// - Returns `false` for values less than 2.
/// - Returns `true` immediately for 2 and 3.
/// - Eliminates multiples of 2 and 3 early.
/// - Then checks divisibility by all numbers of the form 6k ± 1 up to √n.
///
/// # Arguments
///
/// * `n` - The number to check for primality.
///
/// # Returns
///
/// * `true` if `n` is a prime number.
/// * `false` otherwise.
///
/// # Examples
///
/// ```
/// use toy_heaan_ckks::rings::is_prime;
/// assert!(is_prime(7));
/// assert!(!is_prime(9));
/// ```
pub fn is_prime(n: u64) -> bool {
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

/// Generates a list of distinct prime numbers of a given bit size.
///
/// This function searches for odd prime numbers within the specified bit range,
/// starting from the maximum possible value and working downward. It uses
/// a brute-force method with `is_prime()` to check for primality.
///
/// # Arguments
///
/// * `bit_size` - The bit size of each prime (must be between 4 and 63 inclusive).
/// * `count` - The number of primes to generate.
///
/// # Returns
///
/// * `Ok(Vec<u64>)` containing the generated primes if successful.
/// * `Err(RnsError::InvalidPrime)` if not enough primes can be found within the range.
///
/// # Panics
///
/// * Panics if `bit_size` is not between 4 and 63 inclusive.
///
/// # Example
///
/// ```
/// use toy_heaan_ckks::math::{is_prime, generate_primes};
/// let primes = generate_primes(32, 3).unwrap();
/// assert_eq!(primes.len(), 3);
/// for p in primes {
///     assert!(is_prime(p));
///     assert!(p.leading_zeros() <= (64 - 32) as u32);
/// }
/// ```
pub fn generate_primes(bit_size: usize, count: usize) -> Vec<u64> {
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

    // if primes.len() < count {
    //     return Err(RnsError::InvalidPrime(bit_size as u64));
    // }

    primes
}

#[cfg(test)]
mod tests {
    use super::*;

    const KNOWN_PRIMES: [u64; 8] = [2, 3, 5, 7, 11, 13, 17, 19];
    const KNOWN_COMPOSITES: [u64; 10] = [0, 1, 4, 6, 8, 9, 10, 12, 15, 16];

    #[test]
    fn test_is_prime_basic() {
        // Test known primes
        for prime in KNOWN_PRIMES {
            assert!(is_prime(prime));
        }
        // Test known composites
        for prime in KNOWN_COMPOSITES {
            assert!(!is_prime(prime));
        }
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
        let primes = generate_primes(8, 3);

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
        let primes = generate_primes(60, 2);

        assert_eq!(primes.len(), 2);

        for &p in &primes {
            assert!(is_prime(p), "Generated number {} is not prime", p);
            assert!(p >= (1u64 << 59), "Prime {} is too small for 60-bit", p);
            assert!(p < (1u64 << 60), "Prime {} is too large for 60-bit", p);
        }
    }

    /* #[test]
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
    } */
}
