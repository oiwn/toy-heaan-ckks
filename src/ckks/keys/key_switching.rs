//! Gadget key switching: the one operation behind relinearization and
//! rotation (fhetextbook D-3.8; OpenFHE `EvalKey`, BV variant).
//!
//! A polynomial `d` whose decryption contribution is `d * s_old` is turned
//! into a pair `(c0', c1')` with `c0' + c1' * s ~ d * s_old`. The key holds
//! one encryption of `s_old` per gadget term, so `d` can be split into small
//! digits before it meets the key noise.

use rand::Rng;

use crate::PolyRing;
use crate::ckks::keys::SecretKey;

/// Digit width for base-2^`DIGIT_BITS` gadget decomposition.
///
/// Each key-switch digit is bounded by `2^DIGIT_BITS`, so the noise added per
/// accumulated term is ~`2^DIGIT_BITS · σ` instead of ~`q_i · σ` (the old
/// CRT-channel-as-digit scheme). With 61-bit primes this gives
/// `ceil(61 / 20) = 4` digits per channel (one fewer gadget key entry than
/// DIGIT_BITS = 15) at moderately higher noise per switch; demo params
/// (12 primes) have headroom for it.
pub const DIGIT_BITS: u32 = 20;

/// Switches `d * s_old` to a ciphertext under `s`.
///
/// Entry `k` of the key satisfies
///
/// ```text
/// b[k] + a[k] * s ~ s_old.gadget_term(k)   (mod Q)
/// ```
///
/// For RNS, term `k = (i, j)` is digit `j` of channel `i` and the gadget
/// factor is `e_i * 2^(j * DIGIT_BITS)` (`e_i`: the CRT basis element).
/// The key is bound to the level it was generated at.
#[derive(Clone, Debug)]
pub struct KeySwitchingKey<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    /// Uniformly random polynomials, one per gadget term.
    pub a: Vec<P>,
    /// `b[k] = -(a[k] * s) + e_k + s_old.gadget_term(k)`.
    pub b: Vec<P>,
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> KeySwitchingKey<P, DEGREE> {
    /// Encrypts every gadget term of `s_old` under `sk`, at `sk`'s level.
    pub fn generate<R: Rng>(
        sk: &SecretKey<P, DEGREE>,
        s_old: &P,
        error_std: f64,
        rng: &mut R,
    ) -> Self {
        let context = sk.poly.context();
        let terms = P::gadget_len(context, DIGIT_BITS);
        let mut a = Vec::with_capacity(terms);
        let mut b = Vec::with_capacity(terms);

        for k in 0..terms {
            let a_k = P::sample_uniform(context, rng);
            let e_k = P::sample_gaussian(error_std, context, rng);

            // b = -(a · s) + e + s_old scaled by gadget factor k
            let mut a_times_s = a_k.clone();
            a_times_s *= &sk.poly;
            let mut b_k = -a_times_s;
            b_k += &e_k;
            b_k += &s_old.gadget_term(k, DIGIT_BITS);

            a.push(a_k);
            b.push(b_k);
        }

        Self { a, b }
    }

    /// Returns `(c0', c1')` with `c0' + c1' * s ~ d * s_old`:
    ///
    /// ```text
    /// c0' = sum_k decompose(d)[k] * b[k]
    /// c1' = sum_k decompose(d)[k] * a[k]
    /// ```
    ///
    /// Every digit coefficient is below `2^DIGIT_BITS`, which keeps the noise
    /// picked up from the key errors small regardless of the modulus size.
    /// `d` must sit at the key's level.
    pub fn switch(&self, d: &P) -> (P, P) {
        let digits = d.decompose(DIGIT_BITS);
        debug_assert_eq!(
            digits.len(),
            self.a.len(),
            "key switch: key generated at another level"
        );

        let mut c0 = P::zero(d.context());
        let mut c1 = P::zero(d.context());
        for ((digit, a_k), b_k) in digits.iter().zip(&self.a).zip(&self.b) {
            let mut tb = digit.clone();
            tb *= b_k;
            c0 += &tb;

            let mut ta = digit.clone();
            ta *= a_k;
            c1 += &ta;
        }
        (c0, c1)
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use super::*;
    use crate::ckks::keys::SecretKeyParams;
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

    const N: usize = 16;

    /// Three 30-bit primes (Q ~ 2^90), a secret key and a seeded RNG.
    fn setup() -> (Arc<RnsBasis<N>>, SecretKey<RnsPoly<N>, N>, ChaCha20Rng) {
        let basis =
            Arc::new(RnsBasis::new(generate_primes(30, 3, N as u64)).unwrap());
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let params = SecretKeyParams::new(N / 2).unwrap();
        let sk = SecretKey::generate(&params, &basis, &mut rng).unwrap();
        (basis, sk, rng)
    }

    /// Max |coeff| of `c0 + c1 * s - d * s_old` for a uniform `d`.
    fn switch_residual(
        basis: &Arc<RnsBasis<N>>,
        sk: &SecretKey<RnsPoly<N>, N>,
        s_old: &RnsPoly<N>,
        rng: &mut ChaCha20Rng,
    ) -> i64 {
        let ksk = KeySwitchingKey::generate(sk, s_old, 3.2, rng);
        let d = RnsPoly::sample_uniform(basis, rng);
        let (c0, mut lhs) = ksk.switch(&d);
        lhs *= &sk.poly;
        lhs += &c0;
        let mut want = d;
        want *= s_old;
        lhs -= &want;
        lhs.to_coeffs().iter().map(|c| c.abs()).max().unwrap()
    }

    // Noise bound: 6 terms * N * 2^20 (digit) * ~20 (error) < 2^32, far
    // below Q ~ 2^90; a wrong switch leaves Q-sized coefficients.
    const NOISE_BOUND: i64 = 1 << 40;

    #[test]
    fn switch_from_s_squared_leaves_small_noise() {
        let (basis, sk, mut rng) = setup();
        let mut s_sq = sk.poly.clone();
        s_sq *= &sk.poly;
        let residual = switch_residual(&basis, &sk, &s_sq, &mut rng);
        assert!(residual < NOISE_BOUND, "residual {residual}");
    }

    #[test]
    fn switch_from_rotated_key_leaves_small_noise() {
        let (basis, sk, mut rng) = setup();
        let s_rot = sk.poly.rotate_slots(3);
        let residual = switch_residual(&basis, &sk, &s_rot, &mut rng);
        assert!(residual < NOISE_BOUND, "residual {residual}");
    }
}
