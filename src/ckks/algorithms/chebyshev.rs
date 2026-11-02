//! Chebyshev-basis function approximation.
//!
//! [`fit_chebyshev`] interpolates a function on `[lower, upper]`;
//! [`eval_chebyshev_clenshaw`] is the plaintext reference and
//! [`eval_chebyshev_encrypted`] the homomorphic evaluator (a compiled
//! [`ChebyshevPlan`] executed on ciphertexts).

use rand::Rng;

use crate::PolyRing;
use crate::ckks::CkksResult;
use crate::ckks::algorithms::bootstrap::polynomial::ChebyshevPlan;
use crate::ckks::evaluator::LeveledEvaluator;
use crate::ckks::types::Ciphertext;

// ── Encrypted Chebyshev evaluation ───────────────────────────────────────────

/// Evaluates `sum_k coeffs[k] * T_k(t)` on every encrypted slot of `v_ct`,
/// where `t` maps `[lower, upper]` to `[-1, 1]`.
///
/// 1. Affine map `t = (2*v - (hi+lo)) / (hi-lo)`: one constant
///    multiplication via [`Ciphertext::mul_const_to_scale`], so `t` keeps
///    `v_ct`'s scale.
/// 2. A [`ChebyshevPlan`] for `coeffs` (Chebyshev basis throughout; the
///    power basis is ill-conditioned for small `lower`) is compiled for `t`
///    and executed. The output lands at `v_ct`'s scale.
///
/// Levels consumed: `1 + depth`, with `depth = floor(log2(degree)) + 1` for
/// degree >= 1 (0 for a constant). Relinearization keys come from
/// `evaluator`, generated per level on first use.
///
/// Slots must lie in `[lower, upper]`; this is a caller precondition.
///
/// # Errors
/// Returns an error for empty or non-finite coefficients, insufficient
/// levels or modulus headroom (see [`ChebyshevPlan::compile`]), or a failed
/// ciphertext operation.
pub fn eval_chebyshev_encrypted<P, R, const DEGREE: usize>(
    v_ct: &Ciphertext<P, DEGREE>,
    coeffs: &[f64],
    lower: f64,
    upper: f64,
    evaluator: &mut LeveledEvaluator<P, DEGREE>,
    rng: &mut R,
) -> CkksResult<Ciphertext<P, DEGREE>>
where
    P: PolyRing<DEGREE>,
    R: Rng,
{
    let span = upper - lower;
    let t = v_ct
        .mul_const_to_scale(2.0 / span, v_ct.scale)?
        .add_const(-(upper + lower) / span);
    let plan = ChebyshevPlan::compile(&t.c0, t.scale, coeffs, v_ct.scale)?;
    plan.execute(&t, |a, b| Ok(a.mul(b, evaluator.relin_key(a.level(), rng))))
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
    use crate::ckks::Scale;
    use crate::ckks::encoding::ckks_encoder::CkksEncoder;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

    fn test_scale(bits: u32) -> Scale {
        Scale::power_of_two(bits as i32).expect("valid scale")
    }

    /// Encrypted evaluation must track plaintext Clenshaw (same coefficients,
    /// same domain map) to within the accumulated CKKS noise.
    ///
    /// Toy params: N = 16, degree 3 -> levels = affine(1) + plan depth(2) =
    /// 3. Products sit at scale 2^(2*61) = 2^122, which only fits the
    /// modulus (Q ~ 2^(61*chans)) without wrapping when chans >= 3; the plan
    /// checks that headroom.
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
            scale: test_scale(SCALE),
        };
        let engine =
            CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
        let mut rng = ChaCha20Rng::seed_from_u64(1234);

        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");

        let encoder = CkksEncoder::<N>::new(test_scale(SCALE));
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
        )
        .expect("encrypted Chebyshev");

        assert_eq!(
            result.level(),
            4,
            "degree-3 eval from 7 primes: affine + plan depth 2 = 3 levels"
        );
        assert!(result.scale.is_compatible(ct.scale));

        let decoded = encoder.decode(&evaluator.decrypt(&result));
        for (slot, (&got, &want)) in decoded.iter().zip(&expected).enumerate() {
            assert!(
                (got - want).abs() < 1e-3,
                "slot {slot}: encrypted {got} vs plaintext {want}"
            );
        }
    }

    /// Regression: the input fills only 2 of 8 slots (zero-padded). Scalars
    /// used to be broadcast to just the filled slots, which made small
    /// constant encodings round to the zero polynomial; plan constants must
    /// fill every slot.
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
            scale: test_scale(SCALE),
        };
        let engine =
            CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
        let mut rng = ChaCha20Rng::seed_from_u64(99);
        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");

        let encoder = CkksEncoder::<N>::new(test_scale(SCALE));
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
        )
        .expect("encrypted Chebyshev");

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
