//! Modular arithmetic on `u64` residues.
//!
//! [`mul_mod`] is the reference path: the product widens to `u128` and
//! reduces with `%`. That is exact and simple, but the 128-bit remainder is
//! a software routine, so it is slow. Two division-free variants sit next to
//! it (as in OpenFHE and Lattigo):
//!
//! - **Shoup** ([`mul_shoup`]) multiplies by a *fixed* operand `w` with a
//!   precomputed companion `floor(w * 2^64 / q)`: NTT twiddles, twist
//!   factors, scalar constants.
//! - **Barrett** ([`Barrett`]) reduces a general product with a per-modulus
//!   constant `floor(2^128 / q)`.
//!
//! Both return the fully reduced residue in `[0, q)`, which is unique, so
//! they agree with [`mul_mod`] bit for bit. They need `q < 2^62`
//! ([`MAX_FAST_MODULUS`]). Tests compare them with the reference and
//! `benches/arith.rs` measures all three.

/// Fast paths require `modulus < MAX_FAST_MODULUS = 2^62`, so the
/// intermediate `r < 2q` of Shoup and Barrett fits in a `u64`.
pub const MAX_FAST_MODULUS: u64 = 1 << 62;

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
/// Uses [`Barrett`] below [`MAX_FAST_MODULUS`] and the `u128` reference
/// above it; both give the same result.
///
/// # Panics
/// Panics if `modulus` is zero.
#[must_use]
pub fn pow_mod(base: u64, exponent: u64, modulus: u64) -> u64 {
    assert!(modulus > 0, "pow_mod: modulus must be positive");
    if modulus == 1 {
        return 0;
    }
    if modulus < MAX_FAST_MODULUS {
        return Barrett::new(modulus).pow(base, exponent);
    }
    pow_mod_reference(base, exponent, modulus)
}

/// Binary exponentiation over the `u128` reference [`mul_mod`].
///
/// # Panics
/// Panics if `modulus` is zero.
#[must_use]
pub fn pow_mod_reference(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    assert!(modulus > 0, "pow_mod: modulus must be positive");
    let mut acc = 1 % modulus;
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

/// Shoup companion of a fixed operand: `floor(operand * 2^64 / modulus)`.
///
/// Computed once per constant (one `u128` division), then reused by every
/// [`mul_shoup`] with that operand.
///
/// # Panics
/// Panics unless `operand < modulus < 2^62`.
#[must_use]
pub fn shoup_companion(operand: u64, modulus: u64) -> u64 {
    assert!(
        modulus < MAX_FAST_MODULUS,
        "shoup_companion: modulus must be below 2^62"
    );
    assert!(operand < modulus, "shoup_companion: operand not reduced");
    ((u128::from(operand) << 64) / u128::from(modulus)) as u64
}

/// `(value * operand) mod modulus` for a fixed `operand` with its
/// [`shoup_companion`], without division.
///
/// The companion estimates the quotient `value * operand / modulus` from
/// below by at most one, so `value * operand - estimate * modulus` lies in
/// `[0, 2 * modulus)` and one conditional subtraction finishes. The
/// subtraction is computed modulo 2^64, which is exact because the true
/// difference is below `2 * modulus < 2^63`. `value` may be any `u64`.
#[inline]
#[must_use]
pub fn mul_shoup(value: u64, operand: u64, companion: u64, modulus: u64) -> u64 {
    debug_assert!(modulus < MAX_FAST_MODULUS && operand < modulus);
    let estimate = ((u128::from(value) * u128::from(companion)) >> 64) as u64;
    let remainder = value
        .wrapping_mul(operand)
        .wrapping_sub(estimate.wrapping_mul(modulus));
    if remainder >= modulus {
        remainder - modulus
    } else {
        remainder
    }
}

/// Barrett reduction for one modulus: general products without division.
///
/// Stores `ratio = floor(2^128 / modulus)` as two 64-bit words. For an input
/// `x < modulus * 2^64`, the quotient estimate `floor(x * ratio / 2^128)`
/// misses `floor(x / modulus)` by at most one, so one conditional
/// subtraction gives the residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Barrett {
    modulus: u64,
    ratio_hi: u64,
    ratio_lo: u64,
}

impl Barrett {
    /// Precomputes the Barrett constant for `modulus`.
    ///
    /// # Panics
    /// Panics unless `2 <= modulus < 2^62`.
    #[must_use]
    pub fn new(modulus: u64) -> Self {
        assert!(
            (2..MAX_FAST_MODULUS).contains(&modulus),
            "Barrett: modulus must be in [2, 2^62)"
        );
        // floor(2^128 / q) from floor((2^128 - 1) / q): they differ only
        // when q divides 2^128, i.e. when the remainder is q - 1.
        let modulus_wide = u128::from(modulus);
        let mut ratio = u128::MAX / modulus_wide;
        if u128::MAX % modulus_wide == modulus_wide - 1 {
            ratio += 1;
        }
        Self {
            modulus,
            ratio_hi: (ratio >> 64) as u64,
            ratio_lo: ratio as u64,
        }
    }

    /// The modulus this constant reduces by.
    #[must_use]
    pub const fn modulus(&self) -> u64 {
        self.modulus
    }

    /// `value mod modulus` for `value < modulus * 2^64`.
    #[inline]
    #[must_use]
    pub fn reduce_u128(&self, value: u128) -> u64 {
        let (value_hi, value_lo) = ((value >> 64) as u64, value as u64);
        debug_assert!(value_hi < self.modulus, "Barrett: input too large");
        // High 128 bits of value * ratio (a 256-bit product), word by word.
        let wide = |a: u64, b: u64| u128::from(a) * u128::from(b);
        let low_low = wide(value_lo, self.ratio_lo) >> 64;
        let low_high = wide(value_lo, self.ratio_hi);
        let high_low = wide(value_hi, self.ratio_lo);
        let middle =
            low_low + u128::from(low_high as u64) + u128::from(high_low as u64);
        let estimate = wide(value_hi, self.ratio_hi)
            + (low_high >> 64)
            + (high_low >> 64)
            + (middle >> 64);
        // The difference to the true residue is below 2q < 2^63, so the
        // low 64 bits are exact.
        let remainder =
            value_lo.wrapping_sub((estimate as u64).wrapping_mul(self.modulus));
        if remainder >= self.modulus {
            remainder - self.modulus
        } else {
            remainder
        }
    }

    /// `(a * b) mod modulus` for `a, b < modulus`.
    #[inline]
    #[must_use]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.modulus && b < self.modulus);
        self.reduce_u128(u128::from(a) * u128::from(b))
    }

    /// `base^exponent mod modulus`, by binary exponentiation.
    #[must_use]
    pub fn pow(&self, base: u64, mut exponent: u64) -> u64 {
        let mut acc = 1;
        let mut base = base % self.modulus;
        while exponent > 0 {
            if exponent & 1 == 1 {
                acc = self.mul(acc, base);
            }
            base = self.mul(base, base);
            exponent >>= 1;
        }
        acc
    }
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

    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use crate::math::generate_primes;

    /// NTT primes of each width the fast paths support, plus small and
    /// non-prime moduli (Barrett and Shoup work for any modulus).
    fn test_moduli() -> Vec<u64> {
        let mut moduli: Vec<u64> = [30, 50, 59, 60, 61, 62]
            .iter()
            .flat_map(|&bits| generate_primes(bits, 2, 1024))
            .collect();
        moduli.extend([2, 3, 17, 1 << 40, MAX_FAST_MODULUS - 1]);
        moduli
    }

    /// Seeded random residues plus the edge values 0, 1, q - 1.
    fn operands(modulus: u64, rng: &mut ChaCha20Rng) -> Vec<u64> {
        let mut values = vec![0, 1 % modulus, modulus - 1];
        values.extend((0..200).map(|_| rng.random::<u64>() % modulus));
        values
    }

    #[test]
    fn fast_paths_match_u128_reference() {
        let mut rng = ChaCha20Rng::seed_from_u64(52);
        for modulus in test_moduli() {
            let barrett = Barrett::new(modulus);
            let values = operands(modulus, &mut rng);
            for &w in &values {
                let companion = shoup_companion(w, modulus);
                for &a in &values {
                    let expected = mul_mod(a, w, modulus);
                    assert_eq!(barrett.mul(a, w), expected, "barrett q={modulus}");
                    assert_eq!(
                        mul_shoup(a, w, companion, modulus),
                        expected,
                        "shoup q={modulus} a={a} w={w}"
                    );
                }
                // Shoup accepts unreduced values.
                for a in [u64::MAX, u64::MAX - 1, modulus, 2 * modulus] {
                    assert_eq!(
                        mul_shoup(a, w, companion, modulus),
                        mul_mod(a, w, modulus)
                    );
                }
            }
        }
    }

    #[test]
    fn barrett_reduces_wide_inputs() {
        let mut rng = ChaCha20Rng::seed_from_u64(53);
        for modulus in test_moduli() {
            let barrett = Barrett::new(modulus);
            let wide = u128::from(modulus);
            let limit = wide << 64;
            let mut inputs =
                vec![0, (wide - 1) * (wide - 1), limit - 1, wide << 32, wide];
            inputs.extend((0..200).map(|_| rng.random::<u128>() % limit));
            for value in inputs {
                assert_eq!(barrett.reduce_u128(value), (value % wide) as u64);
            }
        }
    }

    #[test]
    fn pow_paths_match_reference() {
        let mut rng = ChaCha20Rng::seed_from_u64(54);
        for modulus in test_moduli().into_iter().chain([u64::MAX - 58]) {
            for _ in 0..20 {
                let base = rng.random::<u64>();
                let exponent = rng.random::<u64>();
                assert_eq!(
                    pow_mod(base, exponent, modulus),
                    pow_mod_reference(base, exponent, modulus)
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "Barrett: modulus must be in [2, 2^62)")]
    fn barrett_rejects_large_modulus() {
        let _ = Barrett::new(MAX_FAST_MODULUS + 1);
    }
}
