//! Educational, fully packed complex CKKS bootstrap.
//!
//! Settings come from [`BootstrapConfig`]; the levels, scales and keys are
//! derived from the modulus chain, which must hold every stage with product
//! headroom and primes close enough to the working scale. Approximation
//! accuracy for a given secret weight and polynomial is the caller's
//! choice (see [`BootstrapConfig::reference`]). Evaluation uses public keys;
//! secret-key diagnostics belong to callers. Hidden slot magnitudes must be
//! at most one.

use std::collections::HashMap;

use rand::Rng;

use crate::PolyRing;
use crate::ckks::Scale;
use crate::ckks::errors::{CkksError, CkksResult};
use crate::ckks::keys::{KeySwitchingKey, RelinearizationKey, SecretKey};
use crate::ckks::types::Ciphertext;

pub mod executor;
pub mod linear;
pub mod polynomial;
pub mod schedule;
pub mod transforms;

use linear::EncodedTransform;
pub use schedule::{
    BootstrapConfig, BootstrapKeyRequirements, BootstrapSchedule, BootstrapStage,
    StageSchedule, reference_ratio_bits,
};
use transforms::FactoredTransform;

fn invalid_parameter(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

/// Public precomputation; contains no secret key.
pub struct BootstrapPlan<P: PolyRing<N>, const N: usize> {
    config: BootstrapConfig,
    schedule: BootstrapSchedule,
    context: P::Context,
    /// Bottom prime `q0`, read from the chain.
    bottom_modulus: u64,
    /// Coefficients-to-slots groups (normalization fused into the first).
    forward: EncodedTransform<P, N>,
    /// Slots-to-coefficients groups (restoration fused into the last).
    inverse: EncodedTransform<P, N>,
}

/// Evaluation keys bound to one plan context. Contains no secret key.
pub struct BootstrapKeys<P: PolyRing<N>, const N: usize> {
    context: P::Context,
    /// Ciphertext level -> relinearization key for multiplication at that level.
    relin: HashMap<usize, RelinearizationKey<P, N>>,
    /// (Ciphertext level, positive slot offset) -> rotation switching key.
    rotations: HashMap<(usize, usize), KeySwitchingKey<P, N>>,
    /// Ciphertext level -> conjugation switching key for splitting complex slots.
    conjugation: HashMap<usize, KeySwitchingKey<P, N>>,
}

impl<P: PolyRing<N>, const N: usize> BootstrapKeys<P, N> {
    /// Check that every key in `required` is present.
    ///
    /// # Errors
    /// Names the first missing key kind and level.
    pub fn check(&self, required: &BootstrapKeyRequirements) -> CkksResult<()> {
        if let Some(level) = required
            .relinearization
            .iter()
            .find(|l| !self.relin.contains_key(l))
        {
            return Err(invalid_parameter(&format!(
                "missing relinearization key at level {level}"
            )));
        }
        if let Some((level, offset)) = required
            .rotations
            .iter()
            .find(|k| !self.rotations.contains_key(k))
        {
            return Err(invalid_parameter(&format!(
                "missing rotation key ({offset}) at level {level}"
            )));
        }
        if let Some(level) = required
            .conjugation
            .iter()
            .find(|l| !self.conjugation.contains_key(l))
        {
            return Err(invalid_parameter(&format!(
                "missing conjugation key at level {level}"
            )));
        }
        Ok(())
    }
}

impl<P: PolyRing<N>, const N: usize> BootstrapPlan<P, N> {
    /// Validate `config` against the chain, derive the schedule and encode
    /// the transform diagonals at their scheduled levels.
    ///
    /// # Errors
    /// Rejects invalid settings (see [`Self::validate_config`]) and
    /// schedules that do not fit the chain (see [`BootstrapSchedule::derive`]).
    pub fn new(context: P::Context, config: BootstrapConfig) -> CkksResult<Self> {
        let zero = P::zero(&context);
        Self::validate_config(&config)?;
        let bottom_modulus = zero.mod_switch_to(1)?.rescale_divisor();
        let coeffs_to_slots = FactoredTransform::slots_to_coefficients(
            N,
            config.coeffs_to_slots_levels,
        )?;
        let slots_to_coeffs = FactoredTransform::coefficients_to_slots(
            N,
            config.slots_to_coeffs_levels,
        )?;
        let schedule = BootstrapSchedule::derive(
            &zero,
            &config,
            &coeffs_to_slots,
            &slots_to_coeffs,
        )?;
        let radius = config.radius();
        let bottom = bottom_modulus as f64;
        // The interval normalization shrinks the signal by 2*radius. Fusing
        // it into the last group keeps the signal large while earlier groups
        // add rotation and rescale noise.
        let forward = EncodedTransform::encode(
            &zero,
            &coeffs_to_slots,
            &schedule
                .stages_of(BootstrapStage::CoeffsToSlots)
                .collect::<Vec<_>>(),
            1.0,
            config.working_scale.value() / bottom / (2.0 * radius),
        )?;
        let inverse = EncodedTransform::encode(
            &zero,
            &slots_to_coeffs,
            &schedule
                .stages_of(BootstrapStage::SlotsToCoeffs)
                .collect::<Vec<_>>(),
            1.0,
            bottom / schedule.message_scale().value() / std::f64::consts::TAU,
        )?;
        Ok(Self {
            config,
            schedule,
            context,
            bottom_modulus,
            forward,
            inverse,
        })
    }

    /// The settings this plan was built from.
    pub fn config(&self) -> &BootstrapConfig {
        &self.config
    }

    /// Levels, scales and keys of every stage.
    pub fn schedule(&self) -> &BootstrapSchedule {
        &self.schedule
    }

    /// Bottom prime `q0` that a mod raise wraps around.
    pub fn bottom_modulus(&self) -> u64 {
        self.bottom_modulus
    }

    /// Checks the settings that do not depend on the chain: a power-of-two
    /// ring degree of at least 8, transform budgets `1..=log2(N/2)`, secret
    /// weight `1..=N` and polynomial degree at least 1. The schedule then
    /// checks the chain: the message ratio against `q0`, chain length,
    /// product headroom and scale drift.
    fn validate_config(config: &BootstrapConfig) -> CkksResult<()> {
        if N < 8 || !N.is_power_of_two() {
            return Err(invalid_parameter(
                "bootstrap ring degree must be a power of two >= 8",
            ));
        }
        let budgets = 1..=(N / 2).ilog2() as usize;
        if !budgets.contains(&config.coeffs_to_slots_levels)
            || !budgets.contains(&config.slots_to_coeffs_levels)
        {
            return Err(invalid_parameter(
                "bootstrap transform budgets must be in 1..=log2(N/2)",
            ));
        }
        if !(1..=N).contains(&config.max_secret_weight) {
            return Err(invalid_parameter(
                "bootstrap secret weight must be in 1..=N",
            ));
        }
        if config.polynomial_degree == 0 {
            return Err(invalid_parameter(
                "bootstrap polynomial degree must be at least 1",
            ));
        }
        Ok(())
    }

    /// Generate level-specific evaluation keys during client setup.
    ///
    /// # Errors
    /// Rejects wrong contexts and invalid noise parameters.
    pub fn generate_keys<R: Rng>(
        &self,
        secret: &SecretKey<P, N>,
        error_std: f64,
        rng: &mut R,
    ) -> CkksResult<BootstrapKeys<P, N>> {
        if !secret.poly.same_context(&P::zero(&self.context))
            || !error_std.is_finite()
            || error_std <= 0.0
        {
            return Err(invalid_parameter(
                "invalid bootstrap secret context or key noise",
            ));
        }
        let mut keys = BootstrapKeys {
            context: self.context.clone(),
            relin: HashMap::new(),
            rotations: HashMap::new(),
            conjugation: HashMap::new(),
        };
        let coefficients = secret.poly.to_coeffs();
        if coefficients.iter().any(|c| ![-1, 0, 1].contains(c))
            || coefficients.iter().filter(|&&c| c != 0).count()
                > self.config.max_secret_weight
        {
            return Err(invalid_parameter(
                "bootstrap requires a ternary secret within the configured weight",
            ));
        }
        let required = self.schedule.required_keys();
        for level in required.levels() {
            let sk = secret.mod_switch_to(level)?;
            if required.relinearization.contains(&level) {
                keys.relin.insert(
                    level,
                    RelinearizationKey::generate(&sk, error_std, rng),
                );
            }
            for &(_, offset) in required.rotations.range((level, 0)..(level + 1, 0))
            {
                let offset_i32 = i32::try_from(offset)
                    .map_err(|_| invalid_parameter("rotation offset too large"))?;
                let rotated = sk.poly.rotate_slots(offset_i32);
                keys.rotations.insert(
                    (level, offset),
                    KeySwitchingKey::generate(&sk, &rotated, error_std, rng),
                );
            }
            if required.conjugation.contains(&level) {
                let conjugated = sk.poly.automorphism((2 * N - 1) as u64);
                keys.conjugation.insert(
                    level,
                    KeySwitchingKey::generate(&sk, &conjugated, error_std, rng),
                );
            }
        }
        Ok(keys)
    }

    fn switch(
        &self,
        ct: &Ciphertext<P, N>,
        exponent: u64,
        key: &KeySwitchingKey<P, N>,
    ) -> Ciphertext<P, N> {
        let mut c0 = ct.c0.automorphism(exponent);
        let (extra, c1) = key.switch(&ct.c1.automorphism(exponent));
        c0 += &extra;
        Ciphertext {
            c0,
            c1,
            scale: ct.scale,
        }
    }

    fn transform(
        &self,
        ct: &Ciphertext<P, N>,
        transform: &EncodedTransform<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ciphertext<P, N>> {
        transform.evaluate(ct, |ct, offset| {
            let key =
                keys.rotations.get(&(ct.level(), offset)).ok_or_else(|| {
                    invalid_parameter(&format!(
                        "missing rotation key ({offset}) at level {}",
                        ct.level()
                    ))
                })?;
            let exponent =
                (0..offset).fold(1u64, |exponent, _| exponent * 5 % (2 * N) as u64);
            Ok(self.switch(ct, exponent, key))
        })
    }

    fn multiply(
        &self,
        a: &Ciphertext<P, N>,
        b: &Ciphertext<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ciphertext<P, N>> {
        let level = a.level().min(b.level());
        if level < 3 {
            return Err(invalid_parameter("insufficient multiplication headroom"));
        }
        let key = keys
            .relin
            .get(&level)
            .ok_or_else(|| invalid_parameter("missing relinearization key"))?;
        a.mod_switch_to(level)?
            .mul(&b.mod_switch_to(level)?, key)
            .rescale()
    }

    fn eval_mod(
        &self,
        ct: &Ciphertext<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ciphertext<P, N>> {
        // Forward diagonals already apply the interval normalization.
        let t = ct.add_const(-0.25 / self.config.radius());
        // The schedule compiled the cosine polynomial for this stage's level
        // and scale; the executor checks every node against it.
        let mut result = self.schedule.eval_mod_plan().execute(&t, |a, b| {
            let key = keys
                .relin
                .get(&a.level())
                .ok_or_else(|| invalid_parameter("missing relinearization key"))?;
            Ok(a.mul(b, key))
        })?;
        for _ in 0..self.config.double_angles {
            let square = self.multiply(&result, &result, keys)?;
            result = square.add(&square)?.add_const(-1.0);
        }
        Ok(result)
    }

    /// Refresh a fully packed complex ciphertext, preparing its scale and level.
    ///
    /// Accepts any level in this plan's chain and any scale `s` from 2 up to
    /// about `message_scale * 2^0.5`, where `message_scale` is
    /// [`BootstrapSchedule::message_scale`]. Preparation lifts the input by the
    /// exact power of two `2^k`, `k = round(log2(message_scale / s))`, so the
    /// lifted scale `s'` is within half a bit of the message scale. The
    /// output is at [`BootstrapSchedule::output_level`] with scale
    /// [`BootstrapSchedule::output_scale`]`(s')`: tracked, not snapped.
    /// Preparation preserves represented values but cannot recover lost
    /// precision. The input is not modified.
    ///
    /// # Errors
    /// Rejects incompatible contexts, unsupported scales and missing keys.
    /// Slot magnitudes <=1 are a hidden-data precondition.
    pub fn bootstrap(
        &self,
        input: &Ciphertext<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ciphertext<P, N>> {
        keys.check(self.schedule.required_keys())?;
        let prepared = self.prepare_input(input, keys)?;
        self.bootstrap_prepared(&prepared, keys)
    }

    fn prepare_input(
        &self,
        input: &Ciphertext<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ciphertext<P, N>> {
        let target = self.schedule.message_scale();
        let increase = target.checked_div(input.scale)?.log2().round();
        if input.scale.log2() < 1.0 || increase < 0.0 {
            return Err(invalid_parameter("invalid bootstrap ciphertext or keys"));
        }
        if !input.c0.same_context(&input.c1)
            || !P::zero(&keys.context).same_context(&P::zero(&self.context))
        {
            return Err(invalid_parameter("invalid bootstrap ciphertext or keys"));
        }
        let expected = P::zero(&self.context).mod_switch_to(input.level())?;
        if !input.c0.same_context(&expected) {
            return Err(invalid_parameter(
                "ciphertext belongs to a different context chain",
            ));
        }
        let normalized = if increase == 0.0 {
            input.clone()
        } else {
            // Scale both the integer representation and its interpretation
            // by the same exact power of two, without changing slot values.
            // `increase` is a small non-negative integer (log2 of a ratio of
            // finite scales), so the cast is exact.
            input.mul_const(1.0, Scale::power_of_two(increase as i32)?)
        };
        normalized.mod_switch_to(1)
    }

    fn bootstrap_prepared(
        &self,
        input: &Ciphertext<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ciphertext<P, N>> {
        let raised = Ciphertext {
            c0: input.c0.mod_raise(&self.context)?,
            c1: input.c1.mod_raise(&self.context)?,
            scale: self.config.working_scale,
        };
        let packed = self.transform(&raised, &self.forward, keys)?;
        let conjugation = keys
            .conjugation
            .get(&packed.level())
            .ok_or_else(|| invalid_parameter("missing conjugation key"))?;
        let conjugate = self.switch(&packed, (2 * N - 1) as u64, conjugation);
        // Half factors were fused into the forward diagonals.
        let real = packed.add(&conjugate)?;
        let imaginary = Self::multiply_by_i(&packed.sub(&conjugate)?, true);
        let real = self.eval_mod(&real, keys)?;
        let imaginary = self.eval_mod(&imaginary, keys)?;
        let combined = real.add(&Self::multiply_by_i(&imaginary, false))?;
        let mut output = self.transform(&combined, &self.inverse, keys)?;
        // The inverse diagonals assume the message scale, so the tracked
        // output scale decodes to `x * s' / message_scale`. Recording
        // `tracked * s' / message_scale` decodes to `x` exactly.
        output.scale = self.schedule.output_scale(input.scale)?;
        Ok(output)
    }

    fn multiply_by_i(input: &Ciphertext<P, N>, negative: bool) -> Ciphertext<P, N> {
        // X^(N/2) evaluates to i at every root indexed by 5^j.
        let mut coefficients = vec![0i64; N];
        coefficients[N / 2] = if negative { -1 } else { 1 };
        let unit = P::from_coeffs(&coefficients, input.c0.context());
        let mut c0 = input.c0.clone();
        c0 *= &unit;
        let mut c1 = input.c1.clone();
        c1 *= &unit;
        Ciphertext {
            c0,
            c1,
            scale: input.scale,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::encoding::ckks_encoder::CkksEncoder;
    use crate::ckks::types::Plaintext;
    use num_complex::Complex64;
    use std::collections::BTreeSet;

    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use std::sync::Arc;

    fn pow2(exponent: i32) -> Scale {
        Scale::power_of_two(exponent).unwrap()
    }

    fn assert_same<const N: usize>(
        left: &Ciphertext<RnsPoly<N>, N>,
        right: &Ciphertext<RnsPoly<N>, N>,
    ) {
        assert_eq!((left.level(), left.scale), (right.level(), right.scale));
        for (a, b) in [(&left.c0, &right.c0), (&left.c1, &right.c1)] {
            let mut difference = a.clone();
            difference -= b;
            assert_eq!(difference.to_coeffs(), [0; N]);
        }
    }

    fn check_stages<const N: usize>() {
        let primes = generate_primes(61, 18, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes.clone()).unwrap());
        let plan = BootstrapPlan::<RnsPoly<N>, N>::new(
            basis.clone(),
            BootstrapConfig::reference(N).unwrap(),
        )
        .unwrap();
        check_schedule(&plan);
        check_invalid_configs(&basis);
        check_multi_group_schedule(&basis);
        let encoder = CkksEncoder::new(pow2(44));
        let engine = CkksEngine::new(
            basis.clone(),
            CkksParams {
                error_variance: 3.2,
                hamming_weight: N / 2,
                scale: pow2(44),
            },
        );
        check_modulus_raise(&basis, &primes);
        let mut rng = ChaCha20Rng::seed_from_u64(0);
        let secret = engine.generate_secret_key(&mut rng).unwrap();
        let public = engine.generate_public_key(&secret, &mut rng).unwrap();
        let mut keys = plan
            .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
            .unwrap();
        let required = plan.schedule().required_keys();
        assert_eq!(
            keys.relin.keys().copied().collect::<BTreeSet<_>>(),
            required.relinearization
        );
        assert_eq!(
            keys.rotations.keys().copied().collect::<BTreeSet<_>>(),
            required.rotations
        );
        assert_eq!(
            keys.conjugation.keys().copied().collect::<BTreeSet<_>>(),
            required.conjugation
        );

        check_modular_reduction(&plan, &basis, &secret, &keys);
        let values: Vec<_> = (0..N / 2)
            .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
            .collect();
        let initial = engine.encrypt(
            &encoder.encode_complex(&values, basis.clone()),
            &public,
            &mut rng,
        );

        check_forward_transform(&plan, &basis, &secret, &keys, &values);
        check_preparation(&plan, &initial, &mut keys);
        let mut lower_scale = initial.clone();
        lower_scale.scale = pow2(27);
        check_preparation(&plan, &lower_scale, &mut keys);
    }

    fn check_modulus_raise<const N: usize>(
        basis: &Arc<RnsBasis<N>>,
        primes: &[u64],
    ) {
        let other = Arc::new(RnsBasis::<N>::new(primes.to_vec()).unwrap());
        let bottom = basis.at_level(1).unwrap();
        let mut coeffs = [0i64; N];
        coeffs[0] = (primes[0] / 2) as i64;
        coeffs[1] = -coeffs[0];
        coeffs[2] = -1;
        let poly = RnsPoly::from_coeffs(&coeffs, &bottom);
        assert_eq!(poly.mod_raise(basis).unwrap().to_coeffs(), coeffs);
        assert!(poly.mod_raise(&other).is_err());
        assert!(poly.mod_raise(&bottom).is_err());
        assert!(RnsPoly::zero(basis).mod_raise(basis).is_err());
    }

    fn check_schedule<const N: usize>(plan: &BootstrapPlan<RnsPoly<N>, N>) {
        use BootstrapStage::*;
        let schedule = plan.schedule();
        // CtS 1 + polynomial 5 + double angles 4 + StC 1.
        assert_eq!(schedule.depth(), 11);
        assert!(schedule.depth() <= 12);
        assert_eq!(schedule.output_level(), 18 - schedule.depth());
        let kinds: Vec<_> = schedule.stages().iter().map(|s| s.stage).collect();
        assert_eq!(
            kinds,
            [
                ModRaise,
                CoeffsToSlots,
                ConjugateSplit,
                EvalModPolynomial,
                DoubleAngle,
                DoubleAngle,
                DoubleAngle,
                DoubleAngle,
                SlotsToCoeffs
            ]
        );
        for pair in schedule.stages().windows(2) {
            assert_eq!(pair[1].input_level, pair[0].output_level);
            assert_eq!(
                pair[1].divisors.len(),
                pair[1].input_level - pair[1].output_level
            );
        }
        let output = schedule.output_scale(pow2(44)).unwrap();
        assert!((output.log2() - 44.0).abs() < 1e-9);
        let required = schedule.required_keys();
        assert_eq!(required.conjugation, BTreeSet::from([17]));
        assert_eq!(required.rotations.len(), 2 * (N / 2 - 1));
        assert!(
            required
                .relinearization
                .is_superset(&BTreeSet::from([12, 11, 10, 9]))
        );
    }

    fn check_invalid_configs<const N: usize>(basis: &Arc<RnsBasis<N>>) {
        let reference = BootstrapConfig::reference(N).unwrap();
        let invalid = [
            BootstrapConfig {
                slots_to_coeffs_levels: (N / 2).ilog2() as usize + 1,
                ..reference.clone()
            },
            BootstrapConfig {
                coeffs_to_slots_levels: 0,
                ..reference.clone()
            },
            BootstrapConfig {
                message_ratio_bits: 0,
                ..reference.clone()
            },
            BootstrapConfig {
                message_ratio_bits: 61,
                ..reference.clone()
            },
            BootstrapConfig {
                polynomial_degree: 0,
                ..reference.clone()
            },
            BootstrapConfig {
                max_secret_weight: 0,
                ..reference.clone()
            },
            BootstrapConfig {
                max_secret_weight: N + 1,
                ..reference.clone()
            },
        ];
        for config in invalid {
            assert!(
                BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), config).is_err()
            );
        }
        // A working scale one bit above the primes doubles its drift with
        // every squaring.
        let oversized = BootstrapConfig {
            working_scale: pow2(62),
            ..reference.clone()
        };
        let error = BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), oversized)
            .err()
            .expect("oversized working scale drifts")
            .to_string();
        assert!(error.contains("drift"), "{error}");
        // Alternating 60/61-bit primes drift the reference working scale.
        let alternating: Vec<_> = generate_primes(61, 9, N as u64)
            .into_iter()
            .zip(generate_primes(60, 9, N as u64))
            .flat_map(|(wide, narrow)| [wide, narrow])
            .collect();
        let alternating = Arc::new(RnsBasis::<N>::new(alternating).unwrap());
        let error =
            BootstrapPlan::<RnsPoly<N>, N>::new(alternating, reference.clone())
                .err()
                .expect("alternating widths drift")
                .to_string();
        assert!(error.contains("drift"), "{error}");
        // Full radix-two budgets need 2*log2(N/2) + 9 levels below the top,
        // and the output keeps one spare level for a squaring: at N=64 that
        // is 19 > 16, so the schedule rejects the 18-prime chain.
        let full = (N / 2).ilog2() as usize;
        let deepest = BootstrapConfig {
            coeffs_to_slots_levels: full,
            slots_to_coeffs_levels: full,
            ..reference.clone()
        };
        let fits = 2 * full + 9 <= 16;
        assert_eq!(
            BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), deepest).is_ok(),
            fits
        );
        // The reference consumes 11 levels and its output must be squared
        // once more: 13 primes are the shortest chain that holds it. With 12
        // the output lands on level 1 and the error names it.
        for (count, fits) in [(13, true), (12, false), (11, false)] {
            let chain = Arc::new(
                RnsBasis::<N>::new(generate_primes(61, count, N as u64)).unwrap(),
            );
            let result =
                BootstrapPlan::<RnsPoly<N>, N>::new(chain, reference.clone());
            assert_eq!(result.is_ok(), fits, "{count} primes");
            if count == 12 {
                let error = result.err().expect("no headroom").to_string();
                assert!(
                    error.contains("output at level 1 cannot be squared"),
                    "{error}"
                );
            }
        }
    }

    fn check_multi_group_schedule<const N: usize>(basis: &Arc<RnsBasis<N>>) {
        let config = BootstrapConfig {
            coeffs_to_slots_levels: 3,
            slots_to_coeffs_levels: 3,
            ..BootstrapConfig::reference(N).unwrap()
        };
        let plan =
            BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), config).unwrap();
        let schedule = plan.schedule();
        assert_eq!(schedule.depth(), 3 + 5 + 4 + 3);
        let required = &schedule.required_keys().rotations;
        for (encoded, stage) in [
            (&plan.forward, BootstrapStage::CoeffsToSlots),
            (&plan.inverse, BootstrapStage::SlotsToCoeffs),
        ] {
            let stages: Vec<_> = schedule.stages_of(stage).collect();
            assert_eq!(encoded.groups().len(), 3);
            for (group, stage) in encoded.groups().iter().zip(stages) {
                assert_eq!(group.level, stage.input_level);
                let offsets: BTreeSet<_> = group
                    .diagonals
                    .iter()
                    .map(|(offset, _)| *offset)
                    .filter(|offset| *offset != 0)
                    .collect();
                let needed: BTreeSet<_> = required
                    .range((group.level, 0)..(group.level + 1, 0))
                    .map(|(_, offset)| *offset)
                    .collect();
                assert_eq!(offsets, needed);
            }
        }
        // Only the last slots-to-coefficients group lands on the message scale.
        let scales: Vec<_> = schedule
            .stages_of(BootstrapStage::SlotsToCoeffs)
            .map(|s| s.plaintext_scale.unwrap())
            .collect();
        assert_eq!(scales, [pow2(61), pow2(61), pow2(44)]);
    }
    fn check_forward_transform<const N: usize>(
        plan: &BootstrapPlan<RnsPoly<N>, N>,
        basis: &Arc<RnsBasis<N>>,
        secret: &SecretKey<RnsPoly<N>, N>,
        keys: &BootstrapKeys<RnsPoly<N>, N>,
        values: &[Complex64],
    ) {
        let encoder = CkksEncoder::new(pow2(44));
        let pt: Plaintext<RnsPoly<N>, N> =
            encoder.encode_complex(values, basis.clone());
        let coefficients = pt.poly.to_centered_f64();
        let trivial = Ciphertext {
            c0: pt.poly,
            c1: RnsPoly::zero(basis),
            scale: pow2(44),
        };
        let transformed = plan.transform(&trivial, &plan.forward, keys).unwrap();
        let decoded = encoder.decode_complex(&CkksEngine::decrypt(
            &transformed,
            &secret.mod_switch_to(transformed.level()).unwrap(),
        ));
        for (j, actual) in decoded.iter().enumerate() {
            let index = j.reverse_bits() >> (usize::BITS - (N / 2).ilog2());
            let want =
                Complex64::new(coefficients[index], coefficients[index + N / 2])
                    / 2f64.powi(44)
                    * 2f64.powi(61)
                    / plan.bottom_modulus as f64;
            let want = want / (2.0 * plan.config.radius());
            assert!((*actual - want).norm() < 1e-9);
        }
    }
    fn check_preparation<const N: usize>(
        plan: &BootstrapPlan<RnsPoly<N>, N>,
        initial: &Ciphertext<RnsPoly<N>, N>,
        keys: &mut BootstrapKeys<RnsPoly<N>, N>,
    ) {
        for level in 1..=18 {
            let at_level = initial.mod_switch_to(level).unwrap();
            let prepared = plan.prepare_input(&at_level, keys).unwrap();
            let bits = initial.scale.exact_log2().unwrap();
            let normalized = if bits == 44 {
                initial.clone()
            } else {
                initial.mul_const(1.0, pow2(44 - bits))
            };
            let explicit = normalized.mod_switch_to(1).unwrap();
            assert_same(&prepared, &explicit);
        }

        let double_angle = plan
            .schedule()
            .stage(BootstrapStage::DoubleAngle)
            .unwrap()
            .input_level;
        let key = keys.relin.remove(&double_angle).unwrap();
        assert!(plan.bootstrap(initial, keys).is_err());
        keys.relin.insert(double_angle, key);

        let prepared = plan.prepare_input(initial, keys).unwrap();
        assert_same(&plan.prepare_input(&prepared, keys).unwrap(), &prepared);
    }

    fn check_modular_reduction<const N: usize>(
        plan: &BootstrapPlan<RnsPoly<N>, N>,
        basis: &Arc<RnsBasis<N>>,
        secret: &SecretKey<RnsPoly<N>, N>,
        keys: &BootstrapKeys<RnsPoly<N>, N>,
    ) {
        // Feed EvalMod at the level and scale its plan was compiled for.
        let split = plan
            .schedule()
            .stage(BootstrapStage::ConjugateSplit)
            .unwrap();
        let ctx = basis.at_level(split.output_level).unwrap();
        let enc = CkksEncoder::new(split.output_scale);
        let radius = plan.config.radius();
        let bound = plan.config.interval_bound();
        for base in (0..2 * bound + 1).step_by(N / 2) {
            let values: Vec<_> = (0..N / 2)
                .map(|j| {
                    let integer = (base + j).min(2 * bound) as f64 - bound as f64;
                    integer
                        + if j % 2 == 0 {
                            2f64.powi(-17)
                        } else {
                            -2f64.powi(-17)
                        }
                })
                .collect();
            let pt: Plaintext<RnsPoly<N>, N> = enc.encode(
                &values.iter().map(|x| x / radius).collect::<Vec<_>>(),
                ctx.clone(),
            );
            let ct = Ciphertext {
                c0: pt.poly,
                c1: RnsPoly::zero(&ctx),
                scale: split.output_scale,
            };
            let reduced = plan.eval_mod(&ct, keys).unwrap();
            let decoded =
                CkksEncoder::new(reduced.scale).decode(&CkksEngine::decrypt(
                    &reduced,
                    &secret.mod_switch_to(reduced.level()).unwrap(),
                ));
            for (j, value) in decoded.into_iter().enumerate() {
                let restored = value / (std::f64::consts::TAU * 2f64.powi(-17));
                assert!(
                    (restored - if j % 2 == 0 { 1.0 } else { -1.0 }).abs() < 1e-4
                );
            }
        }
    }
    /// Build a plan and check its derived depth and message scale.
    fn accepted<const N: usize>(
        bits: usize,
        config: BootstrapConfig,
        polynomial_depth: usize,
        message_bits: i32,
    ) {
        let basis = Arc::new(
            RnsBasis::<N>::new(generate_primes(bits, 18, N as u64)).unwrap(),
        );
        let plan = BootstrapPlan::<RnsPoly<N>, N>::new(basis, config.clone())
            .unwrap_or_else(|e| panic!("N={N} bits={bits} {config:?}: {e}"));
        let schedule = plan.schedule();
        assert_eq!(
            schedule.depth(),
            config.coeffs_to_slots_levels
                + polynomial_depth
                + config.double_angles as usize
                + config.slots_to_coeffs_levels
        );
        assert_eq!(schedule.output_level(), 18 - schedule.depth());
        assert_eq!(schedule.message_scale(), pow2(message_bits));
    }

    #[test]
    fn reference_follows_degree() {
        let ratios: Vec<_> = [16, 64, 128, 256, 1024, 8192, 1 << 14]
            .into_iter()
            .map(reference_ratio_bits)
            .collect();
        assert_eq!(ratios, [17, 17, 10, 9, 7, 4, 4]);
        for (degree, weight, angles) in [(64, 32, 4), (128, 64, 5), (8192, 64, 5)] {
            let config = BootstrapConfig::reference(degree).unwrap();
            assert_eq!(config.max_secret_weight, weight, "N={degree}");
            assert_eq!(config.double_angles, angles, "N={degree}");
        }
    }

    #[test]
    fn accepted_configs() {
        accepted::<32>(61, BootstrapConfig::reference(32).unwrap(), 5, 44);
        // Above N=64 the reference ratio is 17 - log2 N: 2^51 at N=128.
        accepted::<128>(61, BootstrapConfig::reference(128).unwrap(), 5, 51);
        accepted::<1024>(
            61,
            BootstrapConfig {
                coeffs_to_slots_levels: 2,
                slots_to_coeffs_levels: 2,
                ..BootstrapConfig::reference(1024).unwrap()
            },
            5,
            54,
        );
        accepted::<16>(
            60,
            BootstrapConfig {
                working_scale: pow2(60),
                ..BootstrapConfig::reference(16).unwrap()
            },
            5,
            43,
        );
        accepted::<16>(
            61,
            BootstrapConfig {
                polynomial_degree: 63,
                double_angles: 3,
                ..BootstrapConfig::reference(16).unwrap()
            },
            6,
            44,
        );
        accepted::<64>(
            61,
            BootstrapConfig {
                max_secret_weight: 16,
                message_ratio_bits: 20,
                ..BootstrapConfig::reference(64).unwrap()
            },
            5,
            41,
        );
    }

    #[test]
    fn private_stages_16() {
        check_stages::<16>();
    }
    #[test]
    fn private_stages_64() {
        check_stages::<64>();
    }
}
