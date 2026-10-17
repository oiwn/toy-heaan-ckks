//! Chebyshev-basis function approximation.
//!
//! [`fit_chebyshev`] interpolates a function on `[lower, upper]`;
//! [`eval_chebyshev_clenshaw`] is the plaintext reference and
//! [`eval_chebyshev_encrypted`] the homomorphic evaluator (Clenshaw on
//! ciphertexts).

use rand::Rng;

use crate::PolyRing;
use crate::ckks::evaluator::LeveledEvaluator;
use crate::ckks::types::Ciphertext;

// ── Encrypted Chebyshev evaluation ───────────────────────────────────────────

/// Evaluates `sum_k coeffs[k] * T_k(t)` on every encrypted slot of `v_ct`,
/// where `t` maps `[lower, upper]` to `[-1, 1]`, via encrypted Clenshaw.
///
/// Recurrence (Chebyshev basis throughout — converting to power basis is
/// ill-conditioned for small `lower`):
///
/// ```text
/// two_t = 2 * (2*v - (hi+lo)) / (hi-lo)     -- plaintext affine map
/// b_n   = c_n,  b_{n+1} = 0
/// b_k   = two_t * b_{k+1} - b_{k+2} + c_k   -- one ct x ct mul per step
/// result = t * b_1 - b_2 + c_0              -- final ct x ct mul
/// ```
///
/// `c_0` is not a recurrence term; the loop stops at `b_1`. `two_t * b_{k+1}`
/// is rescaled one level down, so `b_{k+2}` is dropped to the new level with
/// [`Ciphertext::match_level`] (exact mod switch) before the subtraction.
///
/// Levels consumed: `degree + 1` (affine map + `degree - 1` recurrence muls +
/// final mul). The final mul needs >= 3 RNS channels of modulus so its raw
/// product (scale ~2^(2*logp)) does not wrap.
///
/// # Panics
/// Panics if `coeffs` is empty or the ciphertext runs out of RNS channels
/// before the evaluation completes.
pub fn eval_chebyshev_encrypted<P, R, const DEGREE: usize>(
    v_ct: &Ciphertext<P, DEGREE>,
    coeffs: &[f64],
    lower: f64,
    upper: f64,
    evaluator: &mut LeveledEvaluator<P, DEGREE>,
    rng: &mut R,
) -> Ciphertext<P, DEGREE>
where
    P: PolyRing<DEGREE>,
    R: Rng,
{
    assert!(!coeffs.is_empty(), "coefficients must not be empty");
    let degree = coeffs.len() - 1;

    // two_t = 4/(hi-lo) * v - 2*(hi+lo)/(hi-lo), then Clenshaw needs only
    // ciphertext multiplications by two_t.
    let span = upper - lower;
    let two_t = v_ct
        .mul_const(4.0 / span, v_ct.logp)
        .rescale()
        .expect("eval_chebyshev_encrypted: rescale")
        .add_const(-2.0 * (upper + lower) / span);

    // b_n = c_n (constant), b_{n+1} = 0 — both as trivial ciphertexts at
    // two_t's level.
    let mut b_next_next = two_t.zero_like();
    let mut b_next = b_next_next.add_const(coeffs[degree]);

    // b_k = two_t * b_{k+1} - b_{k+2} + c_k, for k = n-1 down to 1.
    // (c_0 is NOT a recurrence term; it enters only the final step.)
    for &c_k in coeffs[1..degree].iter().rev() {
        // two_t must meet b_{k+1} at its level (exact mod switch).
        let two_t_aligned = match_level(&two_t, &b_next);

        let mut b_k = evaluator.mul(&two_t_aligned, &b_next, rng);

        // b_{k+2}: drop to the level the multiplication just landed on.
        b_next_next = match_level(&b_next_next, &b_k);
        b_k = (&b_k - &b_next_next).add_const(c_k);

        b_next_next = b_next;
        b_next = b_k;
    }

    // result = t * b_1 - b_2 + c_0 (b_next = b_1, b_next_next = b_2), with
    // t = two_t / 2. The halving is level-free: product * (0.5 encoded at
    // 1 bit -> integer 1) halves the slots and lifts logp by 1, while b_2 *
    // (1.0 at 1 bit -> integer 2) keeps its slots unchanged at the lifted
    // scale — so all three terms meet at logp + 1 without consuming a prime.
    let two_t_aligned = match_level(&two_t, &b_next);
    let product = evaluator.mul(&two_t_aligned, &b_next, rng);

    let t_times_b0 = product.mul_const(0.5, 1);
    let b_one = match_level(&b_next_next, &t_times_b0).mul_const(1.0, 1);

    (&t_times_b0 - &b_one).add_const(coeffs[0])
}

/// [`Ciphertext::match_level`] for the Clenshaw terms, which only ever move
/// down the chain.
fn match_level<P: PolyRing<DEGREE>, const DEGREE: usize>(
    ct: &Ciphertext<P, DEGREE>,
    target: &Ciphertext<P, DEGREE>,
) -> Ciphertext<P, DEGREE> {
    ct.match_level(target)
        .expect("eval_chebyshev_encrypted: term above the target level")
}

/// Maps `x` in `[lower, upper]` to `t` in `[-1, 1]` (Chebyshev domain).
pub fn to_chebyshev_domain(x: f64, lower: f64, upper: f64) -> f64 {
    (2.0 * x - (upper + lower)) / (upper - lower)
}

/// Inverse of [`to_chebyshev_domain`].
pub fn from_chebyshev_domain(t: f64, lower: f64, upper: f64) -> f64 {
    0.5 * ((upper - lower) * t + upper + lower)
}

/// Chebyshev interpolation of `f` on `[lower, upper]`.
///
/// Samples `f` at the `degree + 1` Chebyshev nodes and returns coefficients
/// `c_k` with `f(x) ~ sum_k c_k * T_k(t(x))`, evaluated by
/// [`eval_chebyshev_clenshaw`] (plaintext) or [`eval_chebyshev_encrypted`].
/// Keeping the result in Chebyshev basis (never converting to power basis)
/// keeps the fit well conditioned, e.g. for `sqrt` near a small `lower`.
///
/// # Panics
/// Panics if `upper <= lower`.
pub fn fit_chebyshev(
    f: impl Fn(f64) -> f64,
    lower: f64,
    upper: f64,
    degree: usize,
) -> Vec<f64> {
    assert!(upper > lower, "fit_chebyshev: upper must exceed lower");

    let node_count = degree + 1;
    let mut values = Vec::with_capacity(node_count);

    for j in 0..node_count {
        let theta = std::f64::consts::PI * (j as f64 + 0.5) / node_count as f64;
        let t = theta.cos();
        let x = from_chebyshev_domain(t, lower, upper);
        values.push(f(x));
    }

    let mut coeffs = vec![0.0; node_count];
    for (k, coeff) in coeffs.iter_mut().enumerate() {
        let sum = values
            .iter()
            .enumerate()
            .map(|(j, value)| {
                let theta =
                    std::f64::consts::PI * (j as f64 + 0.5) / node_count as f64;
                value * (k as f64 * theta).cos()
            })
            .sum::<f64>();

        *coeff = if k == 0 {
            sum / node_count as f64
        } else {
            2.0 * sum / node_count as f64
        };
    }

    coeffs
}

/// Evaluates `sum_k c_k * T_k(t)` on `t` in `[-1, 1]` via the Clenshaw
/// recurrence.
///
/// This is the plaintext reference for the encrypted evaluator; it produces
/// bit-identical arithmetic (up to rounding) to the homomorphic version.
///
/// # Panics
/// Panics if `coeffs` is empty.
pub fn eval_chebyshev_clenshaw(t: f64, coeffs: &[f64]) -> f64 {
    assert!(!coeffs.is_empty(), "coefficients must not be empty");

    let mut b_next = 0.0;
    let mut b_next_next = 0.0;

    for &coeff in coeffs.iter().skip(1).rev() {
        let b = 2.0 * t * b_next - b_next_next + coeff;
        b_next_next = b_next;
        b_next = b;
    }

    t * b_next - b_next_next + coeffs[0]
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use super::*;
    use crate::ckks::encoding::ckks_encoder::CkksEncoder;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

    /// Encrypted Clenshaw must track plaintext Clenshaw (same coefficients,
    /// same domain map) to within the accumulated CKKS noise.
    ///
    /// Toy params: N = 16, degree 3 -> levels = affine(1) + steps(2) +
    /// final(1) = 4. The final mul must run at >= 3 channels: its raw
    /// product sits at scale 2^(2*61) = 2^122, which only fits the modulus
    /// (Q ~ 2^(61*chans)) without wrapping when chans >= 3. Budget
    /// degree + 4 primes.
    #[test]
    fn encrypted_clenshaw_matches_plaintext() {
        const N: usize = 16;
        const SLOTS: usize = 8;
        const SCALE: u32 = 61;

        let (lower, upper) = (0.01, 1.0);
        let degree = 3;
        let coeffs = fit_chebyshev(f64::sqrt, lower, upper, degree);

        let input: [f64; SLOTS] = [0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95];
        let expected: Vec<f64> = input
            .iter()
            .map(|&x| {
                let t = to_chebyshev_domain(x, lower, upper);
                eval_chebyshev_clenshaw(t, &coeffs)
            })
            .collect();

        let primes = generate_primes(SCALE as usize, degree + 4, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));

        let params = CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        };
        let engine =
            CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
        let mut rng = ChaCha20Rng::seed_from_u64(1234);

        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");

        let encoder = CkksEncoder::<N>::new(SCALE);
        let pt = encoder.encode(&input, basis);
        let ct = engine.encrypt(&pt, &pk, &mut rng);

        let mut evaluator = LeveledEvaluator::new(params, sk);
        let result = eval_chebyshev_encrypted(
            &ct,
            &coeffs,
            lower,
            upper,
            &mut evaluator,
            &mut rng,
        );

        assert_eq!(
            result.level(),
            3,
            "degree-3 eval from 7 primes: affine + 2 steps + final = 4 levels"
        );

        let decoded = encoder.decode(&evaluator.decrypt(&result));
        for (slot, (&got, &want)) in decoded.iter().zip(&expected).enumerate() {
            assert!(
                (got - want).abs() < 1e-3,
                "slot {slot}: encrypted {got} vs plaintext {want}"
            );
        }
    }

    /// Regression: the input fills only 2 of 8 slots (zero-padded). Scalars
    /// used to be broadcast to just the filled slots, which made the
    /// level-free `0.5 @ 1 bit` halving round to the zero polynomial.
    #[test]
    fn encrypted_clenshaw_with_partially_filled_slots() {
        const N: usize = 16;
        const SCALE: u32 = 61;

        let (lower, upper) = (1e-4, 6.4e-3);
        let degree = 3;
        let coeffs = fit_chebyshev(f64::sqrt, lower, upper, degree);
        let input = [8.2e-4, 3.0e-3];

        let primes = generate_primes(SCALE as usize, degree + 4, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
        let params = CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        };
        let engine =
            CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
        let mut rng = ChaCha20Rng::seed_from_u64(99);
        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");

        let encoder = CkksEncoder::<N>::new(SCALE);
        let pt = encoder.encode(&input, basis);
        let ct = engine.encrypt(&pt, &pk, &mut rng);

        let mut evaluator = LeveledEvaluator::new(params, sk);
        let result = eval_chebyshev_encrypted(
            &ct,
            &coeffs,
            lower,
            upper,
            &mut evaluator,
            &mut rng,
        );

        let decoded = encoder.decode(&evaluator.decrypt(&result));
        for (slot, &x) in input.iter().enumerate() {
            let want = eval_chebyshev_clenshaw(
                to_chebyshev_domain(x, lower, upper),
                &coeffs,
            );
            assert!(
                (decoded[slot] - want).abs() < 1e-9,
                "slot {slot}: encrypted {} vs plaintext {want}",
                decoded[slot]
            );
        }
    }

    #[test]
    fn chebyshev_domain_map_roundtrips() {
        for x in [0.01, 0.1, 0.5, 0.99, 1.0] {
            let t = to_chebyshev_domain(x, 0.01, 1.0);
            assert!(
                (-1.0..=1.0).contains(&t),
                "x={x} maps outside [-1, 1]: t={t}"
            );
            assert!((from_chebyshev_domain(t, 0.01, 1.0) - x).abs() < 1e-15);
        }
    }

    #[test]
    fn fit_sqrt_degree7_accuracy() {
        let coeffs = fit_chebyshev(f64::sqrt, 0.01, 1.0, 7);

        let mut max_err = 0.0f64;
        for i in 0..2001 {
            let x = 0.01 + (1.0 - 0.01) * i as f64 / 2000.0;
            let t = to_chebyshev_domain(x, 0.01, 1.0);
            let err = (eval_chebyshev_clenshaw(t, &coeffs) - x.sqrt()).abs();
            max_err = max_err.max(err);
        }

        // Degree-7 fit of sqrt on [0.01, 1]: worst error near the lower
        // bound.
        assert!(max_err < 0.02, "degree-7 fit max error {max_err:.4e}");
    }

    #[test]
    fn clenshaw_matches_direct_chebyshev_sum() {
        let coeffs = fit_chebyshev(f64::sqrt, 0.01, 1.0, 7);

        for i in -5..=5 {
            let t = i as f64 / 5.0;

            // Direct sum c_k * T_k(t) via the T recurrence.
            let (mut t_prev, mut t_cur) = (1.0, t);
            let mut direct = coeffs[0] + coeffs[1] * t;
            for &c in coeffs.iter().skip(2) {
                let t_next = 2.0 * t * t_cur - t_prev;
                direct += c * t_next;
                t_prev = t_cur;
                t_cur = t_next;
            }

            let clenshaw = eval_chebyshev_clenshaw(t, &coeffs);
            assert!(
                (direct - clenshaw).abs() < 1e-12,
                "t={t}: direct {direct} vs clenshaw {clenshaw}"
            );
        }
    }
}
