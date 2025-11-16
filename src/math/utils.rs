use super::primes::get_first_prime_down;

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
    for (&residue, &prime) in residues.iter().zip(primes.iter()) {
        let partial_product = product / prime;
        let inverse = mod_inverse(partial_product, prime);
        result = (result + residue * partial_product * inverse) % product;
    }

    result
}

/// Generate `count` distinct NTT-friendly primes with the requested bit width.
///
/// The search starts from the largest `bit_size`-wide integer and walks
/// downward in steps of `2 * degree`, ensuring every returned prime satisfies
/// the NTT constraint `p ≡ 1 (mod 2 * degree)`.
///
/// ```
/// use toy_heaan_ckks::math::{generate_primes, is_ntt_friendly_prime};
///
/// let degree = 1024;
/// let primes = generate_primes(32, 3, degree);
/// assert_eq!(primes.len(), 3);
/// for p in primes {
///     assert!(is_ntt_friendly_prime(p, degree as u64));
/// }
/// ```
pub fn generate_primes(bit_size: usize, count: usize, degree: u64) -> Vec<u64> {
    assert!((4..=63).contains(&bit_size));
    assert!(count > 0, "prime count must be positive");
    assert!(degree > 0, "degree must be positive");

    let upper_bound = (1u64 << bit_size).saturating_sub(1);
    let lower_bound = 1u64 << (bit_size - 1);

    let mut primes = Vec::with_capacity(count);
    let mut cursor = get_first_prime_down(upper_bound.saturating_add(1), degree)
        .unwrap_or_else(|| {
            panic!(
                "Failed to find NTT prime below {bit_size} bits for degree {degree}"
            )
        });

    while primes.len() < count {
        if cursor < lower_bound {
            break;
        }
        primes.push(cursor);
        match get_first_prime_down(cursor, degree) {
            Some(next) => cursor = next,
            None => break,
        }
    }

    assert!(
        primes.len() == count,
        "Unable to find {count} NTT primes with {bit_size}-bit ceiling for degree {degree}"
    );

    primes
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::primes::is_ntt_friendly_prime;

    #[test]
    fn generates_ntt_primes_in_range() {
        let degree = 1024;
        let primes = generate_primes(20, 3, degree);
        assert_eq!(primes.len(), 3);

        for &prime in &primes {
            assert!(prime >= (1u64 << 19));
            assert!(prime < (1u64 << 20));
            assert!(is_ntt_friendly_prime(prime, degree));
        }
    }

    #[test]
    #[should_panic(expected = "Unable to find")]
    fn panics_when_not_enough_primes() {
        // Too many primes requested for small bit range should panic
        let _ = generate_primes(4, 10, 2);
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
