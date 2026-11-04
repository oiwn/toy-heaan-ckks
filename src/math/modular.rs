//! Modular arithmetic on `u64` residues: the reference path.
//!
//! Every product widens to `u128` and reduces with `%`. That is exact and
//! simple, but the 128-bit remainder is a software routine, so it is the
//! slow path. Faster variants (Shoup multiplication for a fixed operand,
//! Barrett reduction for general products, task 0052) live next to these
//! functions and are checked and benchmarked against them
//! (`benches/arith.rs`).

/// `(a * b) mod modulus`, via `u128` intermediate arithmetic.
///
/// # Panics
/// Panics if `modulus` is zero.
#[inline]
#[must_use]
pub fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    assert!(modulus > 0, "mul_mod: modulus must be positive");
    ((u128::from(a) * u128::from(b)) % u128::from(modulus)) as u64
}

/// `(a + b) mod modulus` for residues `a, b < modulus`.
#[inline]
#[must_use]
pub fn add_mod(a: u64, b: u64, modulus: u64) -> u64 {
    debug_assert!(a < modulus && b < modulus, "add_mod: inputs not reduced");
    let sum = a + b;
    if sum >= modulus { sum - modulus } else { sum }
}

/// `(a - b) mod modulus` for residues `a, b < modulus`.
#[inline]
#[must_use]
pub fn sub_mod(a: u64, b: u64, modulus: u64) -> u64 {
    debug_assert!(a < modulus && b < modulus, "sub_mod: inputs not reduced");
    if a >= b { a - b } else { a + modulus - b }
}

/// `base^exponent mod modulus`, by binary exponentiation.
///
/// # Panics
/// Panics if `modulus` is zero.
#[must_use]
pub fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    assert!(modulus > 0, "pow_mod: modulus must be positive");
    if modulus == 1 {
        return 0;
    }
    let mut acc = 1;
    base %= modulus;
    while exponent > 0 {
        if exponent & 1 == 1 {
            acc = mul_mod(acc, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn add_and_sub_wrap_at_the_modulus() {
        let modulus = 97;
        assert_eq!(add_mod(96, 5, modulus), 4);
        assert_eq!(add_mod(3, 5, modulus), 8);
        assert_eq!(sub_mod(3, 5, modulus), 95);
        assert_eq!(sub_mod(5, 3, modulus), 2);
    }

    #[test]
    fn pow_mod_handles_edge_cases() {
        assert_eq!(pow_mod(2, 0, 17), 1);
        assert_eq!(pow_mod(5, 0, 1), 0);
        assert_eq!(pow_mod(0, 5, 17), 0);
        assert_eq!(pow_mod(7, 1, 19), 7);
        assert_eq!(pow_mod(3, 16, 17), 1);
    }

    #[test]
    #[should_panic(expected = "pow_mod: modulus must be positive")]
    fn pow_mod_panics_on_zero_modulus() {
        let _ = pow_mod(2, 10, 0);
    }
}
