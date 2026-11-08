//! A bounded nonlinear state update for illustrating repeated refresh.
//!
//! For real encrypted slots u and x in [-1,1], compute
//! x_next = x/2 + x*x/8 + sign*u/4, with alternating sign.
//! Evaluation uses public keys only. Slot bounds are caller preconditions.

use std::collections::HashMap;

use rand::Rng;

use crate::PolyRing;
use crate::ckks::Scale;
use crate::ckks::errors::{CkksError, CkksResult};
use crate::ckks::keys::{RelinearizationKey, SecretKey};
use crate::ckks::types::Ciphertext;

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

/// Accepted input scales lie within this many bits of the plan scale.
const INPUT_WINDOW_BITS: f64 = 1.0;
/// `1/8` is encoded exactly with three fractional bits.
const EIGHTH_BITS: i32 = 3;

/// Fixed-scale educational recurrence over a shared modulus chain.
pub struct RecurrencePlan<P: PolyRing<N>, const N: usize> {
    context: P::Context,
    scale: Scale,
    refresh_level: usize,
}

/// Level-specific public evaluation keys, generated before computation.
pub struct RecurrenceKeys<P: PolyRing<N>, const N: usize> {
    context: P::Context,
    relin: HashMap<usize, RelinearizationKey<P, N>>,
}

impl<P: PolyRing<N>, const N: usize> RecurrencePlan<P, N> {
    /// Construct the recurrence for states and drivers at `scale`.
    ///
    /// The chain decides how deep a step may go. A step at level L squares
    /// the state (scale up to `2 * scale_bits + 2`, plus 3 bits for `1/8`),
    /// rescales by the level's prime and lifts the result back to `scale`.
    /// The refresh level is the lowest L >= 2 whose modulus holds that
    /// product and whose next level holds the lifted result. Every prime
    /// from there up must be wide enough that the rescaled result does not
    /// exceed `scale`, so the lift is a non-negative power of two.
    /// Precision lost in that rescale is not checked (the caller's choice of
    /// scale and chain, as with bootstrap).
    ///
    /// # Errors
    /// Rejects a chain with no level that fits the product, or a prime too
    /// narrow for the lift; the message names the level.
    pub fn new(context: P::Context, scale: Scale) -> CkksResult<Self> {
        let zero = P::zero(&context);
        let scale_bits = scale.log2();
        let product_bits =
            2.0 * (scale_bits + INPUT_WINDOW_BITS) + f64::from(EIGHTH_BITS);
        let fits = |bits: f64, level: usize| -> CkksResult<bool> {
            Ok(bits + 1.0 < f64::from(zero.mod_switch_to(level)?.modulus_bits()))
        };
        let mut refresh_level = None;
        for level in 2..=zero.level() {
            if fits(product_bits, level)? && fits(scale_bits + 0.5, level - 1)? {
                refresh_level = Some(level);
                break;
            }
        }
        let refresh_level = refresh_level.ok_or_else(|| {
            invalid(&format!(
                "recurrence chain of {} levels has no level holding a \
                 2^{product_bits:.0} product",
                zero.level()
            ))
        })?;
        for level in refresh_level..=zero.level() {
            let divisor = zero.mod_switch_to(level)?.rescale_divisor() as f64;
            // The rescaled result must not exceed `scale` by half a bit.
            if product_bits - divisor.log2() > scale_bits + 0.5 {
                return Err(invalid(&format!(
                    "recurrence prime at level {level} is too narrow to \
                     rescale a 2^{product_bits:.0} product back to scale"
                )));
            }
        }
        Ok(Self {
            context,
            scale,
            refresh_level,
        })
    }

    /// The scale every step returns to, within half a bit.
    pub const fn scale(&self) -> Scale {
        self.scale
    }

    /// The lowest level a state may have before a step.
    pub const fn refresh_level(&self) -> usize {
        self.refresh_level
    }

    /// Whether the state is too low in the chain for another step.
    /// This only checks depth; `step` validates contexts and scales.
    pub fn needs_refresh(&self, state: &Ciphertext<P, N>) -> bool {
        state.level() < self.refresh_level
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
        for level in self.refresh_level..=secret.poly.level() {
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

    /// Compute one round, consuming one prime.
    /// Round zero uses a positive driver; odd rounds use a negative driver.
    /// Inputs are not modified. No secret key or RNG is used.
    ///
    /// State and driver scales may be any value less than one bit from the
    /// plan scale: rescaling divides by the actual prime, so tracked scales
    /// drift slightly from powers of two. The output is lifted by an exact
    /// `2^k` to within half a bit of the plan scale, which keeps the scale
    /// bounded across rounds without snapping it.
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
        let target = self.scale;
        let supported =
            |scale: Scale| (scale.log2() - target.log2()).abs() < INPUT_WINDOW_BITS;
        if self.needs_refresh(state)
            || !supported(state.scale)
            || !supported(driver.scale)
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
        // meet at the square's scale plus those three bits before rescaling:
        // each one is encoded at `square / its own scale`.
        let square = state
            .mul(state, key)
            .mul_const(0.125, Scale::power_of_two(EIGHTH_BITS)?);
        let linear = state.mul_const(0.5, square.scale.checked_div(state.scale)?);
        let signed_driver = driver.mul_const(
            if round.is_multiple_of(2) { 0.25 } else { -0.25 },
            square.scale.checked_div(driver.scale)?,
        );
        let result = square.add(&linear)?.add(&signed_driver)?.rescale()?;
        // Lift by the exact power of two that brings the scale within half
        // a bit of the target (rounding, not flooring: a result just above
        // 2^30 must lift to about 2^44, not to about 2^43, at target 2^44).
        let increase = target.checked_div(result.scale)?.log2().round();
        if increase < 0.0 {
            return Err(invalid("unsupported recurrence output scale"));
        }
        Ok(if increase == 0.0 {
            result
        } else {
            // A small non-negative integer: the cast is exact.
            result.mul_const(1.0, Scale::power_of_two(increase as i32)?)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
    use crate::ckks::encoding::ckks_encoder::CkksEncoder;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha20Rng;
    use std::sync::Arc;

    fn pow2(exponent: i32) -> Scale {
        Scale::power_of_two(exponent).unwrap()
    }

    /// `step` lifts its output to within half a bit of the plan scale.
    fn in_output_window(scale: Scale, target: Scale) -> bool {
        (scale.log2() - target.log2()).abs() <= 0.5
    }

    /// Rounds completed (of 32) and refreshes, derived from the levels: a
    /// step needs `refresh_level` and drops one level; a refresh returns
    /// to `output_level`.
    fn expected_rounds(
        top: usize,
        refresh_level: usize,
        output_level: Option<usize>,
    ) -> (usize, usize) {
        let first = top + 1 - refresh_level;
        match output_level {
            _ if first >= 32 => (32, 0),
            None => (first, 0),
            Some(output) => {
                let per_cycle = output + 1 - refresh_level;
                (32, (32 - first).div_ceil(per_cycle))
            }
        }
    }

    fn scenario<const N: usize>(count: usize, refresh: bool, seed: u64) {
        let primes = generate_primes(61, count, N as u64);
        let config = BootstrapConfig::reference(N).unwrap();
        let context = Arc::new(RnsBasis::<N>::new(primes.clone()).unwrap());
        let other = Arc::new(RnsBasis::<N>::new(primes).unwrap());
        let bootstrap = refresh
            .then(|| BootstrapPlan::new(context.clone(), config.clone()).unwrap());
        let scale = bootstrap
            .as_ref()
            .map_or(pow2(44), |p| p.schedule().message_scale());
        let plan =
            RecurrencePlan::<RnsPoly<N>, N>::new(context.clone(), scale).unwrap();
        let engine = CkksEngine::new(
            context.clone(),
            CkksParams {
                error_variance: 3.2,
                hamming_weight: config.max_secret_weight,
                scale,
            },
        );
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let secret = engine.generate_secret_key(&mut rng).unwrap();
        let public = engine.generate_public_key(&secret, &mut rng).unwrap();
        let mut keys = plan
            .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
            .unwrap();
        let boot_keys = bootstrap
            .as_ref()
            .map(|p| p.generate_keys(&secret, 3.2f64.sqrt(), &mut rng).unwrap());
        let encoder = CkksEncoder::new(scale);
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
            let lowest = plan.refresh_level();
            let low = driver.mod_switch_to(lowest).unwrap();
            for round in [0, 1] {
                let next = plan.step(&low, &driver, round, &keys).unwrap();
                assert_eq!(next.level(), lowest - 1);
                assert!(in_output_window(next.scale, scale));
                let decoded = encoder.decode(&CkksEngine::decrypt(
                    &next,
                    &secret.mod_switch_to(lowest - 1).unwrap(),
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
        bad.scale = scale.checked_div(pow2(2)).unwrap();
        assert!(plan.step(&bad, &driver, 0, &keys).is_err());
        bad = state.clone();
        bad.c0 = RnsPoly::zero(&other);
        assert!(plan.step(&bad, &driver, 0, &keys).is_err());
        bad.c1 = bad.c0.clone();
        assert!(plan.step(&bad, &driver, 0, &keys).is_err());
        let below = state.mod_switch_to(plan.refresh_level() - 1).unwrap();
        assert!(plan.step(&below, &driver, 0, &keys).is_err());
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
            assert_eq!(state.level(), level - 1);
            assert!(in_output_window(state.scale, scale));
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
        let output = bootstrap.as_ref().map(|p| p.schedule().output_level());
        assert_eq!(
            (rounds, refreshes),
            expected_rounds(count, plan.refresh_level(), output)
        );
        for (a, b) in [(&driver.c0, &snapshot.c0), (&driver.c1, &snapshot.c1)] {
            let mut delta = a.clone();
            delta -= b;
            assert_eq!(delta.to_coeffs(), [0; N]);
        }
    }

    #[test]
    fn derived_rounds_match_the_reference_profile() {
        // 61-bit primes at 2^44: steps need level 2, bootstrap returns
        // level 7, so 18 primes give 17 rounds unrefreshed, or 32 rounds
        // with three refreshes.
        let context =
            Arc::new(RnsBasis::<16>::new(generate_primes(61, 18, 16)).unwrap());
        let plan =
            RecurrencePlan::<RnsPoly<16>, 16>::new(context, pow2(44)).unwrap();
        assert_eq!(plan.refresh_level(), 2);
        assert_eq!(expected_rounds(18, 2, None), (17, 0));
        assert_eq!(expected_rounds(18, 2, Some(7)), (32, 3));
        assert_eq!(expected_rounds(34, 2, None), (32, 0));
    }

    /// A 50-bit chain at scale 2^40 runs without refresh to its last level.
    #[test]
    fn mixed_width_chain() {
        const N: usize = 32;
        let scale = pow2(40);
        let context =
            Arc::new(RnsBasis::<N>::new(generate_primes(50, 6, N as u64)).unwrap());
        let plan =
            RecurrencePlan::<RnsPoly<N>, N>::new(context.clone(), scale).unwrap();
        assert_eq!(plan.refresh_level(), 2);
        let engine = CkksEngine::new(
            context.clone(),
            CkksParams {
                error_variance: 3.2,
                hamming_weight: N / 2,
                scale,
            },
        );
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let secret = engine.generate_secret_key(&mut rng).unwrap();
        let public = engine.generate_public_key(&secret, &mut rng).unwrap();
        let keys = plan
            .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
            .unwrap();
        let values: Vec<_> =
            (0..N / 2).map(|_| rng.random_range(-1.0..1.0)).collect();
        let encoder = CkksEncoder::new(scale);
        let driver = engine.encrypt(
            &encoder.encode(&values, context.clone()),
            &public,
            &mut rng,
        );
        let mut state = driver.clone();
        let mut expected = values.clone();
        let mut round = 0;
        while !plan.needs_refresh(&state) {
            state = plan.step(&state, &driver, round, &keys).unwrap();
            assert!(in_output_window(state.scale, scale));
            for (x, u) in expected.iter_mut().zip(&values) {
                *x = *x / 2.0
                    + *x * *x / 8.0
                    + if round % 2 == 0 { *u / 4.0 } else { -*u / 4.0 };
            }
            let decoded = encoder.decode(&CkksEngine::decrypt(
                &state,
                &secret.mod_switch_to(state.level()).unwrap(),
            ));
            for (a, b) in decoded.iter().zip(&expected) {
                assert!((a - b).abs() < 1e-3, "round={round} got={a} want={b}");
            }
            round += 1;
        }
        assert_eq!(round, 5);
    }

    #[test]
    fn rejects_chains_that_cannot_hold_a_step() {
        let chain = |bits: usize, count: usize| {
            Arc::new(RnsBasis::<16>::new(generate_primes(bits, count, 16)).unwrap())
        };
        let error = |context, exponent| {
            RecurrencePlan::<RnsPoly<16>, 16>::new(context, pow2(exponent))
                .err()
                .expect("rejected")
                .to_string()
        };
        // One prime: no level to rescale into.
        assert!(error(chain(61, 1), 44).contains("no level holding"));
        // Two 61-bit primes cannot hold a 2^125 product.
        assert!(error(chain(61, 2), 60).contains("no level holding"));
        // Three primes hold it at level 3, but a 61-bit prime leaves the
        // rescaled result above 2^60.
        assert!(error(chain(61, 3), 60).contains("level 3 is too narrow"));
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
    fn recurrence_comparisons_32() {
        for (count, refresh) in [(18, true), (18, false)] {
            scenario::<32>(count, refresh, 0);
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
