//! Prime utilities for constructing NTT-friendly modulus sets.

const MILLER_RABIN_BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut acc = 1u64;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exp >>= 1;
    }
    acc
}

fn decompose(n: u64) -> (u64, u32) {
    let mut d = n;
    let mut r = 0;
    while d & 1 == 0 {
        d >>= 1;
        r += 1;
    }
    (d, r)
}

/// Deterministic Miller–Rabin for all `u64` inputs.
pub fn is_prime(n: u64) -> bool {
    match n {
        0 | 1 => return false,
        2 | 3 => return true,
        _ if n & 1 == 0 => return false,
        _ => {}
    }

    let (d, r) = decompose(n - 1);
    'bases: for &a in MILLER_RABIN_BASES.iter() {
        if a >= n {
            continue;
        }
        let mut x = mod_pow(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..r {
            x = mul_mod(x, x, n);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// Slow-but-clear reference test used for benches and debugging.
pub fn is_prime_reference(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n.is_multiple_of(2) || n.is_multiple_of(3) {
        return false;
    }
    let mut i = 5u64;
    while i.saturating_mul(i) <= n {
        if n.is_multiple_of(i) || n.is_multiple_of(i + 2) {
            return false;
        }
        i += 6;
    }
    true
}

/// Returns `true` when `p` is prime and congruent to `1 (mod 2n)`.
#[inline]
pub fn is_ntt_friendly_prime(p: u64, n: u64) -> bool {
    is_prime(p) && p % (2 * n) == 1
}

fn snap_up_to_congruence(value: u64, modulus: u64) -> u64 {
    if modulus == 0 {
        return value;
    }
    let remainder = value % modulus;
    if remainder == 1 {
        value
    } else {
        let delta = (modulus + 1 - remainder) % modulus;
        value + delta
    }
}

fn snap_down_to_congruence(value: u64, modulus: u64) -> Option<u64> {
    if modulus == 0 {
        return Some(value);
    }
    if value < 1 {
        return None;
    }
    let remainder = value % modulus;
    let delta = (remainder + modulus - 1) % modulus;
    value.checked_sub(delta)
}

/// Find the first NTT-friendly prime at or above `2^logq`.
pub fn get_first_prime_up(logq: u32, n: u64) -> u64 {
    if logq >= 64 {
        panic!("logq must be less than 64 to fit in u64");
    }
    assert!(n > 0, "n must be positive");

    let m = 2 * n; // Step size for NTT constraint
    let mut candidate = snap_up_to_congruence((1u64 << logq) + 1, m);

    loop {
        if is_prime(candidate) {
            return candidate;
        }
        candidate += m; // Step by 2N
    }
}

/// Walks downward in steps of `2n` to find the next NTT-friendly prime.
pub fn get_first_prime_down(bound: u64, n: u64) -> Option<u64> {
    if bound <= 2 || n == 0 {
        return None;
    }

    let step = 2 * n;
    let mut candidate = snap_down_to_congruence(bound.saturating_sub(1), step)?;

    loop {
        if candidate <= 2 {
            return None;
        }
        if is_prime(candidate) {
            return Some(candidate);
        }
        if candidate <= step + 1 {
            return None;
        }
        candidate = candidate.saturating_sub(step);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const KNOWN_SMALL_PRIMES: [u64; 8] = [2, 3, 5, 7, 11, 13, 17, 19];
    const KNOWN_SMALL_COMPOSITES: [u64; 10] = [0, 1, 4, 6, 8, 9, 10, 12, 15, 16];

    #[test]
    fn test_is_prime_basic() {
        for &prime in &KNOWN_SMALL_PRIMES {
            assert!(is_prime(prime));
            assert!(is_prime_reference(prime));
        }
        for &composite in &KNOWN_SMALL_COMPOSITES {
            assert!(!is_prime(composite));
            assert!(!is_prime_reference(composite));
        }
    }

    #[test]
    fn test_is_prime_large() {
        assert!(is_prime(65537));
        assert!(is_prime(982_451_653));
        assert!(is_prime(2_147_483_647));

        assert!(!is_prime(65536));
        assert!(!is_prime(982_451_654));
        assert!(!is_prime(2_147_483_648));
    }

    #[test]
    fn miller_rabin_matches_reference_on_range() {
        for n in 2..10_000 {
            assert_eq!(is_prime(n), is_prime_reference(n), "mismatch at {n}");
        }
    }

    #[test]
    fn test_ntt_friendly_condition() {
        // Test that 12289 is NTT-friendly for N=1024
        assert!(is_ntt_friendly_prime(12289, 1024));
        assert_eq!(12289 % (2 * 1024), 1);

        // Test that composite numbers are not NTT-friendly
        assert!(!is_ntt_friendly_prime(2049, 1024));
        assert!(!is_ntt_friendly_prime(4097, 1024));
    }

    #[test]
    fn first_prime_up_matches_reference() {
        let prime = get_first_prime_up(30, 1024);
        assert_eq!(prime, 1_073_750_017);
    }

    #[test]
    fn prime_down_descends() {
        let prime = get_first_prime_up(20, 1024);
        let below = get_first_prime_down(prime, 1024).unwrap();
        assert!(below < prime);
        assert!(is_ntt_friendly_prime(below, 1024));
    }

    #[test]
    fn test_prime_down_basic() {
        let n = 1024;
        let bound = 2u64.pow(31);
        if let Some(prime) = get_first_prime_down(bound, n) {
            assert!(prime < bound);
            assert!(is_ntt_friendly_prime(prime, n));
        }

        assert_eq!(get_first_prime_down(2, n), None);
        assert_eq!(get_first_prime_down(1, n), None);
    }
}
