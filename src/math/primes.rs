//! Prime utilities for constructing NTT-friendly modulus sets.
//!
//! This module uses the Miller-Rabin primality test as its fast prime check.
//! Miller-Rabin is a probabilistic test in general: it checks whether a number
//! behaves like a prime for a set of chosen bases.
//!
//! For `u64` values, using a fixed set of known bases makes the test
//! deterministic in practice for the full input range we care about. The
//! implementation below decomposes `n - 1` into `d * 2^r`, then verifies that
//! each base either:
//! - produces `1` or `n - 1` directly, or
//! - reaches `n - 1` after repeated squaring modulo `n`.
//!
//! If no base can witness compositeness, the number is treated as prime.
//! Reference:
//! https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test

// These bases are deterministic for all n < 318,665,857,834,031,151,167,461,
// which covers all u64 values.
// Source: https://miller-rabin.appspot.com/
const MILLER_RABIN_BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Computes `(a * b) mod modulus` using `u128` intermediate arithmetic.
fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    assert!(modulus > 0, "mul_mod: modulus must be positive");
    ((a as u128 * b as u128) % modulus as u128) as u64
}

/// Computes `base^exp mod modulus` via binary exponentiation.
fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    assert!(modulus > 0, "mod_pow: modulus must be positive");
    if modulus == 1 {
        return 0;
    }
    let mut acc = 1 % modulus;
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

/// Returns `(odd_part, power_of_two)` such that `n = odd_part * 2^power_of_two`.
fn decompose(n: u64) -> (u64, u32) {
    assert!(n > 0, "decompose: n must be positive");
    let mut d = n;
    let mut r = 0;
    while d & 1 == 0 {
        d >>= 1;
        r += 1;
    }
    (d, r)
}

/// Returns `true` if `n` is prime using deterministic Miller-Rabin on `u64`.
///
/// The function first handles small numbers directly, then writes `n - 1` as
/// `d * 2^r` with odd `d`. For each fixed base `a`, it checks whether:
/// - `a^d mod n` is `1` or `n - 1`, or
/// - repeated squaring reaches `n - 1` within `r - 1` steps.
///
/// If any base fails both checks, `n` is composite. Otherwise `n` is prime.
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

/// Slow-but-clear reference test using `6k +/- 1` trial division.
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

/// Returns `true` when `p` is an NTT-friendly prime for ring degree `n`.
///
/// Here "NTT-friendly" means:
/// - `p` is prime, and
/// - `p = 1 (mod 2n)` (equivalently, `2n` divides `p - 1`).
///
/// This condition ensures `Z_p` contains a primitive `2n`-th root of unity,
/// which is required for negacyclic NTT over `x^n + 1`.
#[inline]
pub fn is_ntt_friendly_prime(p: u64, n: u64) -> bool {
    assert!(n > 0, "is_ntt_friendly_prime: n must be positive");
    let modulus = n
        .checked_mul(2)
        .expect("is_ntt_friendly_prime: 2 * n must fit in u64");
    is_prime(p) && p % modulus == 1
}

/// Returns the smallest `x >= value` such that `x % modulus == 1`.
fn snap_up_to_congruence(value: u64, modulus: u64) -> u64 {
    assert!(modulus > 1, "snap_up_to_congruence: modulus must be greater than 1");
    let remainder = value % modulus;
    if remainder == 1 {
        value
    } else {
        let delta = (modulus + 1 - remainder) % modulus;
        value
            .checked_add(delta)
            .expect("snap_up_to_congruence: overflow while stepping upward")
    }
}

/// Returns the largest `x <= value` such that `x % modulus == 1`.
fn snap_down_to_congruence(value: u64, modulus: u64) -> u64 {
    assert!(
        modulus > 1,
        "snap_down_to_congruence: modulus must be greater than 1"
    );
    let remainder = value % modulus;
    let delta = (remainder + modulus - 1) % modulus;
    value
        .checked_sub(delta)
        .expect("snap_down_to_congruence: underflow while stepping downward")
}

/// Returns the first NTT-friendly prime `p` such that `p >= 2^logq`.
///
/// The search only visits candidates `p = 1 (mod 2n)`.
///
/// # Panics
///
/// Panics if `logq >= 64`, `n == 0`, `2 * n` overflows `u64`, or stepping
/// upward overflows before a prime is found.
pub fn get_first_prime_up(logq: u32, n: u64) -> u64 {
    assert!(logq < 64, "get_first_prime_up: logq must be less than 64");
    assert!(n > 0, "get_first_prime_up: n must be positive");

    let step = n
        .checked_mul(2)
        .expect("get_first_prime_up: 2 * n must fit in u64");
    let mut candidate = snap_up_to_congruence((1u64 << logq) + 1, step);

    loop {
        if is_prime(candidate) {
            return candidate;
        }
        candidate = candidate
            .checked_add(step)
            .expect("get_first_prime_up: overflow while stepping upward");
    }
}

/// Returns the largest NTT-friendly prime `p` such that `p < bound`.
///
/// The search only visits candidates `p = 1 (mod 2n)`.
///
/// Returns `None` if no such prime exists below `bound`.
///
/// # Panics
///
/// Panics if `n == 0` or `2 * n` overflows `u64`.
pub fn get_first_prime_down(bound: u64, n: u64) -> Option<u64> {
    assert!(n > 0, "get_first_prime_down: n must be positive");
    if bound <= 2 {
        return None;
    }

    let step = n
        .checked_mul(2)
        .expect("get_first_prime_down: 2 * n must fit in u64");
    let mut candidate = snap_down_to_congruence(bound.saturating_sub(1), step);

    loop {
        if candidate <= 2 {
            return None;
        }
        if is_prime(candidate) {
            return Some(candidate);
        }
        candidate = candidate.checked_sub(step)?;
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
    fn mul_mod_matches_widened_reference() {
        let a = u64::MAX - 11;
        let b = u64::MAX - 17;
        let modulus = 1_073_750_017u64;
        let expected = ((a as u128 * b as u128) % modulus as u128) as u64;
        assert_eq!(mul_mod(a, b, modulus), expected);
    }

    #[test]
    #[should_panic(expected = "mul_mod: modulus must be positive")]
    fn mul_mod_panics_on_zero_modulus() {
        let _ = mul_mod(5, 7, 0);
    }

    #[test]
    fn mod_pow_handles_edge_cases() {
        assert_eq!(mod_pow(2, 0, 17), 1);
        assert_eq!(mod_pow(5, 0, 1), 0);
        assert_eq!(mod_pow(0, 5, 17), 0);
        assert_eq!(mod_pow(7, 1, 19), 7);
    }

    #[test]
    #[should_panic(expected = "mod_pow: modulus must be positive")]
    fn mod_pow_panics_on_zero_modulus() {
        let _ = mod_pow(2, 10, 0);
    }

    #[test]
    fn decompose_splits_power_of_two_factor() {
        assert_eq!(decompose(24), (3, 3));
        assert_eq!(decompose(40), (5, 3));
        assert_eq!(decompose(1), (1, 0));
    }

    #[test]
    #[should_panic(expected = "decompose: n must be positive")]
    fn decompose_panics_on_zero() {
        let _ = decompose(0);
    }

    #[test]
    fn snap_up_to_congruence_aligns_upward() {
        assert_eq!(snap_up_to_congruence(1, 8), 1);
        assert_eq!(snap_up_to_congruence(10, 8), 17);
        assert_eq!(snap_up_to_congruence(17, 8), 17);
    }

    #[test]
    #[should_panic(expected = "snap_up_to_congruence: modulus must be greater than 1")]
    fn snap_up_to_congruence_panics_on_zero_modulus() {
        let _ = snap_up_to_congruence(10, 0);
    }

    #[test]
    #[should_panic(expected = "snap_up_to_congruence: modulus must be greater than 1")]
    fn snap_up_to_congruence_panics_on_modulus_one() {
        let _ = snap_up_to_congruence(10, 1);
    }

    #[test]
    #[should_panic(expected = "snap_up_to_congruence: overflow while stepping upward")]
    fn snap_up_to_congruence_panics_on_overflow() {
        let _ = snap_up_to_congruence(u64::MAX, 4);
    }

    #[test]
    fn snap_down_to_congruence_aligns_downward() {
        assert_eq!(snap_down_to_congruence(17, 8), 17);
        assert_eq!(snap_down_to_congruence(16, 8), 9);
        assert_eq!(snap_down_to_congruence(8, 8), 1);
    }

    #[test]
    #[should_panic(expected = "snap_down_to_congruence: modulus must be greater than 1")]
    fn snap_down_to_congruence_panics_on_zero_modulus() {
        let _ = snap_down_to_congruence(10, 0);
    }

    #[test]
    #[should_panic(expected = "snap_down_to_congruence: modulus must be greater than 1")]
    fn snap_down_to_congruence_panics_on_modulus_one() {
        let _ = snap_down_to_congruence(10, 1);
    }

    #[test]
    #[should_panic(expected = "snap_down_to_congruence: underflow while stepping downward")]
    fn snap_down_to_congruence_panics_on_underflow() {
        let _ = snap_down_to_congruence(0, 8);
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
    fn test_is_prime_tricky_composites() {
        // Carmichael numbers and strong pseudoprimes for small base sets.
        let tricky = [561u64, 1_105, 1_729, 3_215_031_751];
        for &n in &tricky {
            assert!(!is_prime(n), "expected composite: {n}");
        }
    }

    #[test]
    fn test_is_prime_near_u64_limit() {
        assert!(!is_prime(u64::MAX));
        assert!(is_prime(18_446_744_073_709_551_557));
    }

    #[test]
    fn miller_rabin_matches_reference_on_selected_ranges() {
        let ranges: [(u64, u64); 4] = [
            (2, 14),
            (90, 114),
            (10_000, 10_024),
            (1_000_000, 1_000_024),
        ];

        for (start, end) in ranges {
            for n in start..=end {
                assert_eq!(is_prime(n), is_prime_reference(n), "mismatch at {n}");
            }
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
    #[should_panic(expected = "is_ntt_friendly_prime: n must be positive")]
    fn ntt_friendly_panics_on_zero_degree() {
        let _ = is_ntt_friendly_prime(12289, 0);
    }

    #[test]
    #[should_panic(expected = "is_ntt_friendly_prime: 2 * n must fit in u64")]
    fn ntt_friendly_panics_on_degree_overflow() {
        let too_large_n = (u64::MAX / 2) + 1;
        let _ = is_ntt_friendly_prime(12289, too_large_n);
    }

    #[test]
    fn first_prime_up_matches_reference() {
        let prime = get_first_prime_up(30, 1024);
        assert_eq!(prime, 1_073_750_017);
    }

    #[test]
    #[should_panic(expected = "get_first_prime_up: logq must be less than 64")]
    fn first_prime_up_panics_on_large_logq() {
        let _ = get_first_prime_up(64, 1024);
    }

    #[test]
    #[should_panic(expected = "get_first_prime_up: n must be positive")]
    fn first_prime_up_panics_on_zero_degree() {
        let _ = get_first_prime_up(30, 0);
    }

    #[test]
    #[should_panic(expected = "get_first_prime_up: 2 * n must fit in u64")]
    fn first_prime_up_panics_on_degree_overflow() {
        let too_large_n = (u64::MAX / 2) + 1;
        let _ = get_first_prime_up(30, too_large_n);
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

    #[test]
    #[should_panic(expected = "get_first_prime_down: n must be positive")]
    fn first_prime_down_panics_on_zero_degree() {
        let _ = get_first_prime_down(100, 0);
    }

    #[test]
    #[should_panic(expected = "get_first_prime_down: 2 * n must fit in u64")]
    fn first_prime_down_panics_on_degree_overflow() {
        let too_large_n = (u64::MAX / 2) + 1;
        let _ = get_first_prime_down(100, too_large_n);
    }
}
