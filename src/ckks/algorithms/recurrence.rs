//! A bounded nonlinear state update for illustrating repeated refresh.
//!
//! For real encrypted slots u and x in [-1,1], compute
//! x_next = x/2 + x*x/8 + sign*u/4, with alternating sign.
//! Evaluation uses public keys only. Slot bounds are caller preconditions.

use std::collections::HashMap;

use rand::Rng;

use crate::PolyRing;
use crate::ckks::errors::{CkksError, CkksResult};
use crate::ckks::keys::{RelinearizationKey, SecretKey};
use crate::ckks::types::Ciphertext;

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

/// Fixed-scale educational recurrence over a shared modulus chain.
pub struct RecurrencePlan<P: PolyRing<N>, const N: usize> {
    context: P::Context,
}

/// Level-specific public evaluation keys, generated before computation.
pub struct RecurrenceKeys<P: PolyRing<N>, const N: usize> {
    context: P::Context,
    relin: HashMap<usize, RelinearizationKey<P, N>>,
}

impl<P: PolyRing<N>, const N: usize> RecurrencePlan<P, N> {
    /// Construct the toy recurrence with scale2^44 and near-2^61 divisors.
    ///
    /// # Errors
    /// Rejects unsupported degrees, short chains or incompatible divisors.
    pub fn new(context: P::Context) -> CkksResult<Self> {
        let zero = P::zero(&context);
        if ![16, 64].contains(&N) || zero.level() < 2 {
            return Err(invalid(
                "recurrence requires degree16/64 and at least two primes",
            ));
        }
        for level in 1..=zero.level() {
            if zero
                .mod_switch_to(level)?
                .rescale_divisor()
                .abs_diff(1u64 << 61)
                > 1u64 << 20
            {
                return Err(invalid("recurrence requires near-2^61 divisors"));
            }
        }
        Ok(Self { context })
    }

    /// Whether another step needs refresh for its verified 91-bit product.
    /// Two near-2^61 primes provide a modulus near2^122; after rescale the
    /// scale2^30 result and exact lift to2^44 fit the bottom prime.
    /// This only checks depth; `step` validates contexts and scales.
    pub fn needs_refresh(&self, state: &Ciphertext<P, N>) -> bool {
        state.level() < 2
    }

    /// Generate all recurrence keys once during client setup.
    ///
    /// # Errors
    /// Rejects a mismatched secret context or invalid noise standard deviation.
    pub fn generate_keys<R: Rng>(
        &self,
        secret: &SecretKey<P, N>,
        error_std: f64,
        rng: &mut R,
    ) -> CkksResult<RecurrenceKeys<P, N>> {
        if !P::zero(&self.context).same_context(&secret.poly)
            || !error_std.is_finite()
            || error_std <= 0.0
        {
            return Err(invalid("invalid recurrence key setup"));
        }
        let mut relin = HashMap::new();
        for level in 2..=secret.poly.level() {
            relin.insert(
                level,
                RelinearizationKey::generate(
                    &secret.mod_switch_to(level)?,
                    error_std,
                    rng,
                ),
            );
        }
        Ok(RecurrenceKeys {
            context: self.context.clone(),
            relin,
        })
    }

    /// Compute one round, consuming one prime and returning scale2^44.
    /// Round zero uses a positive driver; odd rounds use a negative driver.
    /// Inputs are not modified. No secret key or RNG is used.
    ///
    /// # Errors
    /// Rejects unsupported scales, wrong contexts, missing keys or insufficient
    /// levels. Real slot magnitudes <=1 are a hidden-data precondition.
    pub fn step(
        &self,
        state: &Ciphertext<P, N>,
        driver: &Ciphertext<P, N>,
        round: usize,
        keys: &RecurrenceKeys<P, N>,
    ) -> CkksResult<Ciphertext<P, N>> {
        if self.needs_refresh(state)
            || state.logp != 44
            || driver.logp != 44
            || !state.c0.same_context(&state.c1)
            || !driver.c0.same_context(&driver.c1)
            || !P::zero(&self.context).same_context(&P::zero(&keys.context))
        {
            return Err(invalid(
                "invalid recurrence inputs, keys or insufficient levels",
            ));
        }
        let zero = P::zero(&self.context);
        if !state.c0.same_context(&zero.mod_switch_to(state.level())?)
            || !driver.c0.same_context(&zero.mod_switch_to(driver.level())?)
        {
            return Err(invalid("recurrence inputs belong to a different chain"));
        }
        let driver = driver.mod_switch_to(state.level())?;
        let key = keys
            .relin
            .get(&state.level())
            .ok_or_else(|| invalid("missing recurrence key"))?;
        // Encode 1/8 exactly with three fractional bits. All weighted terms
        // meet at the square's scale plus those three bits before rescaling.
        let square = state.mul(state, key).mul_const(0.125, 3);
        let linear_bits = square.logp - state.logp;
        let linear = state.mul_const(0.5, linear_bits);
        let signed_driver = driver.mul_const(
            if round.is_multiple_of(2) { 0.25 } else { -0.25 },
            linear_bits,
        );
        let result = (&(&square + &linear) + &signed_driver).rescale()?;
        let increase = 44u32
            .checked_sub(result.logp)
            .ok_or_else(|| invalid("unsupported recurrence output scale"))?;
        Ok(if increase == 0 {
            result
        } else {
            result.mul_const(1.0, increase)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::algorithms::bootstrap::{BootstrapParameters, BootstrapPlan};
    use crate::ckks::encoding::ckks_encoder::CkksEncoder;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha20Rng;
    use std::sync::Arc;

    fn scenario<const N: usize>(count: usize, refresh: bool, seed: u64) {
        let primes = generate_primes(61, count, N as u64);
        let parameters = BootstrapParameters::toy(primes[0]);
        let context = Arc::new(RnsBasis::<N>::new(primes.clone()).unwrap());
        let other = Arc::new(RnsBasis::<N>::new(primes).unwrap());
        let plan = RecurrencePlan::new(context.clone()).unwrap();
        let engine = CkksEngine::new(
            context.clone(),
            CkksParams {
                error_variance: 3.2,
                hamming_weight: N / 2,
                scale_bits: 44,
            },
        );
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let secret = engine.generate_secret_key(&mut rng).unwrap();
        let public = engine.generate_public_key(&secret, &mut rng).unwrap();
        let mut keys = plan
            .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
            .unwrap();
        let bootstrap = refresh
            .then(|| BootstrapPlan::new(context.clone(), parameters).unwrap());
        let boot_keys = bootstrap
            .as_ref()
            .map(|p| p.generate_keys(&secret, 3.2f64.sqrt(), &mut rng).unwrap());
        let encoder = CkksEncoder::new(44);
        let values: Vec<_> = (0..N / 2)
            .map(|j| match j {
                0 => -1.0,
                1 => 0.0,
                2 => 1.0,
                _ => rng.random_range(-1.0..1.0),
            })
            .collect();
        let driver = engine.encrypt(
            &encoder.encode(&values, context.clone()),
            &public,
            &mut rng,
        );
        let snapshot = driver.clone();
        // Exercise the lowest supported level at boundary slot values,
        // not only the smaller values reached after many recurrence rounds.
        if seed == 0 {
            let low = driver.mod_switch_to(2).unwrap();
            for round in [0, 1] {
                let next = plan.step(&low, &driver, round, &keys).unwrap();
                assert_eq!((next.level(), next.logp), (1, 44));
                let decoded = encoder.decode(&CkksEngine::decrypt(
                    &next,
                    &secret.mod_switch_to(1).unwrap(),
                ));
                for (a, u) in decoded.iter().zip(&values) {
                    let want = u / 2.0
                        + u * u / 8.0
                        + if round == 0 { u / 4.0 } else { -u / 4.0 };
                    assert!(a.is_finite() && (a - want).abs() < 1e-3);
                }
            }
        }
        let mut state = driver.clone();
        let mut expected = values.clone();
        let key = keys.relin.remove(&count).unwrap();
        assert!(plan.step(&state, &driver, 0, &keys).is_err());
        keys.relin.insert(count, key);
        let mut bad = state.clone();
        bad.logp = 43;
        assert!(plan.step(&bad, &driver, 0, &keys).is_err());
        bad = state.clone();
        bad.c0 = RnsPoly::zero(&other);
        assert!(plan.step(&bad, &driver, 0, &keys).is_err());
        bad.c1 = bad.c0.clone();
        assert!(plan.step(&bad, &driver, 0, &keys).is_err());
        assert!(
            plan.step(&state.mod_switch_to(1).unwrap(), &driver, 0, &keys)
                .is_err()
        );
        let mut rounds = 0;
        let mut refreshes = 0;
        let check = |state: &Ciphertext<RnsPoly<N>, N>, expected: &[f64]| {
            let decoded = encoder.decode(&CkksEngine::decrypt(
                state,
                &secret.mod_switch_to(state.level()).unwrap(),
            ));
            for (a, b) in decoded.iter().zip(expected) {
                assert!(
                    a.is_finite() && (a - b).abs() < 1e-3,
                    "N={N} seed={seed} got={a} want={b}"
                );
            }
        };
        for round in 0..32 {
            if plan.needs_refresh(&state) {
                let Some(bootstrap) = &bootstrap else { break };
                state = bootstrap
                    .bootstrap(&state, boot_keys.as_ref().unwrap())
                    .unwrap();
                refreshes += 1;
                check(&state, &expected);
            }
            let before = state.clone();
            let level = state.level();
            state = plan.step(&state, &driver, round, &keys).unwrap();
            assert_eq!((state.level(), state.logp), (level - 1, 44));
            // Check the saved input against a second deterministic evaluation.
            let repeated = plan.step(&before, &driver, round, &keys).unwrap();
            for (a, b) in [(&state.c0, &repeated.c0), (&state.c1, &repeated.c1)] {
                let mut delta = a.clone();
                delta -= b;
                assert_eq!(delta.to_coeffs(), [0; N]);
            }
            for (x, u) in expected.iter_mut().zip(&values) {
                *x = (*x) / 2.0
                    + (*x).powi(2) / 8.0
                    + if round % 2 == 0 { *u / 4.0 } else { -*u / 4.0 };
                assert!((-0.625..=0.875).contains(x));
            }
            check(&state, &expected);
            rounds += 1;
        }
        assert_eq!(rounds, if count == 18 && !refresh { 17 } else { 32 });
        assert_eq!(refreshes, if refresh { 3 } else { 0 });
        for (a, b) in [(&driver.c0, &snapshot.c0), (&driver.c1, &snapshot.c1)] {
            let mut delta = a.clone();
            delta -= b;
            assert_eq!(delta.to_coeffs(), [0; N]);
        }
    }

    #[test]
    fn recurrence_comparisons_16() {
        for seed in 0..2 {
            for (count, refresh) in [(18, true), (18, false), (34, false)] {
                scenario::<16>(count, refresh, seed);
            }
        }
    }
    #[test]
    fn recurrence_comparisons_64() {
        for seed in 0..2 {
            for (count, refresh) in [(18, true), (18, false), (34, false)] {
                scenario::<64>(count, refresh, seed);
            }
        }
    }
}
