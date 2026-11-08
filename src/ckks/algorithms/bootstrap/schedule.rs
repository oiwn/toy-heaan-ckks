//! Bootstrap configuration and the level/scale schedule derived from it.
//!
//! The schedule simulates every stage on the actual modulus chain with the
//! same `Scale` operations the encrypted evaluation performs, so planned and
//! actual scales agree up to `Scale::is_compatible`.

use std::collections::BTreeSet;
use std::f64::consts::TAU;

use super::polynomial::{ChebyshevPlan, Operation};
use super::transforms::{DiagonalFactor, FactoredTransform};
use crate::PolyRing;
use crate::ckks::algorithms::chebyshev::fit_chebyshev;
use crate::ckks::{CkksError, CkksResult, Scale};

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

/// Largest distance, in bits, between the working scale and any rescaled
/// EvalMod or double-angle scale. Rescale divides by the actual prime, so a
/// prime unlike the working scale shifts the scale, and every squaring
/// doubles the shift; past this bound the polynomial loses its precision.
pub const MAX_SCALE_DRIFT_BITS: f64 = 1.0;

/// Largest secret weight of the reference settings. The mod-raise overflow
/// bound, and with it the EvalMod interval, grows with the weight; 64
/// (radius 34) is what the large-degree measurements use.
pub const MAX_REFERENCE_WEIGHT: usize = 64;

/// Message ratio of the reference settings at small degrees (message scale
/// 2^44 on a 61-bit `q0`); recurrence and the small-degree demos rely on it.
pub const SMALL_DEGREE_RATIO_BITS: u32 = 17;

/// Largest degree that keeps [`SMALL_DEGREE_RATIO_BITS`].
pub const SMALL_DEGREE_LIMIT: usize = 64;

/// Smallest reference ratio; below it the sine linearization error grows
/// 16x per 2 bits (measured).
pub const MIN_REFERENCE_RATIO_BITS: u32 = 4;

/// Reference message ratio for ring degree `ring_degree`.
///
/// The bootstrap error is U-shaped in the ratio: above the optimum the
/// noise before EvalMod is amplified by `q0 / message_scale` (doubling per
/// bit), below it the cubic sine linearization error takes over. The
/// optimum drops about one bit per doubling of N, because large-N message
/// coefficients are small while key-switching noise grows. Fitted to the
/// sweep in `specs/research/large-degree-bootstrap.md` (20-bit digits,
/// weight min(N/2, 64)): `17 - log2 N` above N=64, at least 4, which gives
/// 10 at N=128 and 4 at N=8192. N <= 64 keeps 17 (worst measured 2.6e-4).
/// OpenFHE fits a similar line in log2 N for its correction factor. Hybrid
/// key switching (task 0048) lowers the noise, so the fit is redone there.
pub fn reference_ratio_bits(ring_degree: usize) -> u32 {
    if ring_degree <= SMALL_DEGREE_LIMIT {
        return SMALL_DEGREE_RATIO_BITS;
    }
    SMALL_DEGREE_RATIO_BITS
        .saturating_sub(ring_degree.ilog2())
        .max(MIN_REFERENCE_RATIO_BITS)
}

/// Numerical settings of a bootstrap.
#[derive(Clone, Debug, PartialEq)]
pub struct BootstrapConfig {
    /// `log2(q0 / message_scale)`, like Lattigo's `LogMessageRatio`. The
    /// message scale (what inputs are lifted to and what the inverse
    /// transform assumes) is derived from the chain's bottom prime `q0`, so
    /// the tradeoff between sine linearization error and EvalMod error
    /// amplification carries over to other prime widths.
    pub message_ratio_bits: u32,
    /// Scale of the raised ciphertext, the forward diagonals and EvalMod.
    pub working_scale: Scale,
    /// Degree of the Chebyshev cosine approximation.
    pub polynomial_degree: usize,
    /// Number of double-angle steps after the cosine polynomial.
    pub double_angles: u32,
    /// Levels used by the coefficients-to-slots transform.
    pub coeffs_to_slots_levels: usize,
    /// Levels used by the slots-to-coefficients transform.
    pub slots_to_coeffs_levels: usize,
    /// Largest Hamming weight of the ternary secret.
    pub max_secret_weight: usize,
}

impl BootstrapConfig {
    /// The reference settings, measured with near-2^61 primes: message
    /// ratio from [`reference_ratio_bits`] (2^17, message scale 2^44, up to
    /// N=64), working scale 2^61, degree 31, secret weight min(N/2, 64),
    /// four double angles (five above weight 32) and one-level transforms.
    ///
    /// # Errors
    /// Never fails for these constants; the signature follows `Scale`.
    pub fn reference(ring_degree: usize) -> CkksResult<Self> {
        let max_secret_weight = (ring_degree / 2).min(MAX_REFERENCE_WEIGHT);
        Ok(Self {
            message_ratio_bits: reference_ratio_bits(ring_degree),
            working_scale: Scale::power_of_two(61)?,
            polynomial_degree: 31,
            // A wider EvalMod interval needs one more halving of the period.
            double_angles: if max_secret_weight > 32 { 5 } else { 4 },
            coeffs_to_slots_levels: 1,
            slots_to_coeffs_levels: 1,
            max_secret_weight,
        })
    }

    /// The message scale for bottom prime `q0`: the power of two
    /// `2^(round(log2 q0) - message_ratio_bits)`.
    ///
    /// # Errors
    /// Rejects a zero ratio and ratios that leave a scale below 2.
    pub fn message_scale(&self, bottom_modulus: u64) -> CkksResult<Scale> {
        let bits = (bottom_modulus as f64).log2().round() as i64
            - i64::from(self.message_ratio_bits);
        if self.message_ratio_bits == 0 || bits < 1 {
            return Err(invalid(
                "message ratio must be at least 1 bit and leave a scale of 2 or more",
            ));
        }
        let bits = i32::try_from(bits)
            .map_err(|_| invalid("message scale is not representable"))?;
        Ok(Scale::power_of_two(bits)?)
    }

    /// Bound on the integer overflow `I` in `m + q0*I` after a mod raise,
    /// for a ternary secret of at most `max_secret_weight` nonzeros.
    pub fn interval_bound(&self) -> usize {
        (self.max_secret_weight + 2) / 2
    }

    /// Half-width of the EvalMod input interval, one above the bound.
    pub fn radius(&self) -> f64 {
        (self.interval_bound() + 1) as f64
    }

    /// Chebyshev coefficients of `cos(2*pi*radius*t / 2^double_angles)` on
    /// `[-1, 1]`; the double angles undo the period division.
    pub fn cosine_coefficients(&self) -> Vec<f64> {
        let radius = self.radius();
        let period = f64::from(self.double_angles).exp2();
        fit_chebyshev(
            |t| (TAU * radius * t / period).cos(),
            -1.0,
            1.0,
            self.polynomial_degree,
        )
    }
}

/// Named bootstrap stages in evaluation order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BootstrapStage {
    /// Reinterpret the bottom-level ciphertext in the full chain.
    ModRaise,
    /// One coefficients-to-slots group (diagonal products, one rescale).
    CoeffsToSlots,
    /// Split complex slots into real and imaginary parts (conjugation).
    ConjugateSplit,
    /// The compiled cosine polynomial.
    EvalModPolynomial,
    /// One double-angle step `2*x^2 - 1`.
    DoubleAngle,
    /// One slots-to-coefficients group (diagonal products, one rescale).
    SlotsToCoeffs,
}

/// Level and exact scale after one stage.
#[derive(Clone, Debug)]
pub struct StageSchedule {
    pub stage: BootstrapStage,
    pub input_level: usize,
    pub output_level: usize,
    pub output_scale: Scale,
    /// Scale of the stage's encoded diagonals (transform stages only).
    pub plaintext_scale: Option<Scale>,
    /// Rescale divisors consumed by the stage, top first.
    pub divisors: Vec<u64>,
}

/// Evaluation keys a schedule needs, by ciphertext level.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BootstrapKeyRequirements {
    pub relinearization: BTreeSet<usize>,
    /// (level, positive slot offset).
    pub rotations: BTreeSet<(usize, usize)>,
    pub conjugation: BTreeSet<usize>,
}

impl BootstrapKeyRequirements {
    /// Every level at which some key is needed.
    pub fn levels(&self) -> BTreeSet<usize> {
        let mut levels = self.relinearization.clone();
        levels.extend(self.rotations.iter().map(|(level, _)| *level));
        levels.extend(&self.conjugation);
        levels
    }
}

/// Levels, scales and keys of every bootstrap stage on one modulus chain.
#[derive(Clone, Debug)]
pub struct BootstrapSchedule {
    stages: Vec<StageSchedule>,
    keys: BootstrapKeyRequirements,
    eval_mod: ChebyshevPlan,
    message_scale: Scale,
}

impl BootstrapSchedule {
    /// Derive the schedule from the top-level ring element of a chain.
    ///
    /// # Errors
    /// Rejects chains too short for the stages, nodes without headroom,
    /// primes that drift the working scale by more than
    /// [`MAX_SCALE_DRIFT_BITS`] and unrepresentable scales.
    pub fn derive<P: PolyRing<N>, const N: usize>(
        top: &P,
        config: &BootstrapConfig,
        coeffs_to_slots: &FactoredTransform,
        slots_to_coeffs: &FactoredTransform,
    ) -> CkksResult<Self> {
        let message_scale =
            config.message_scale(top.mod_switch_to(1)?.rescale_divisor())?;
        let mut builder = Builder::<P, N> {
            top,
            stages: Vec::new(),
            keys: BootstrapKeyRequirements::default(),
            level: top.level(),
            scale: config.working_scale,
        };
        builder.stages.push(StageSchedule {
            stage: BootstrapStage::ModRaise,
            input_level: 1,
            output_level: builder.level,
            output_scale: builder.scale,
            plaintext_scale: None,
            divisors: Vec::new(),
        });
        for factor in coeffs_to_slots.factors() {
            builder.transform(
                BootstrapStage::CoeffsToSlots,
                factor,
                config.working_scale,
            )?;
        }

        builder.keys.conjugation.insert(builder.level);
        builder.push(BootstrapStage::ConjugateSplit, builder.level, Vec::new());

        let input = top.mod_switch_to(builder.level)?;
        let eval_mod = ChebyshevPlan::compile(
            &input,
            builder.scale,
            &config.cosine_coefficients(),
            builder.scale,
        )?;
        for node in eval_mod.values() {
            if matches!(node.operation, Operation::Input | Operation::Rescale(_)) {
                check_drift(config.working_scale, node.scale, node.level)?;
            }
        }
        builder
            .keys
            .relinearization
            .extend(eval_mod.required_relinearization_levels());
        let output = &eval_mod.values()[eval_mod.output()];
        let input_level = builder.level;
        let divisors = (output.level + 1..=input_level)
            .rev()
            .map(|level| Ok(top.mod_switch_to(level)?.rescale_divisor()))
            .collect::<CkksResult<_>>()?;
        builder.level = output.level;
        builder.scale = output.scale;
        builder.push(BootstrapStage::EvalModPolynomial, input_level, divisors);

        for _ in 0..config.double_angles {
            builder.keys.relinearization.insert(builder.level);
            let square = builder.scale.checked_mul(builder.scale)?;
            builder.headroom(square, 1.0)?;
            builder.rescale(BootstrapStage::DoubleAngle, square)?;
            check_drift(config.working_scale, builder.scale, builder.level)?;
        }
        // Earlier groups keep the working scale; the last one lands on the
        // message scale.
        let last = slots_to_coeffs.factors().len().saturating_sub(1);
        for (index, factor) in slots_to_coeffs.factors().iter().enumerate() {
            let diagonal_scale = if index == last {
                message_scale
            } else {
                config.working_scale
            };
            builder.transform(
                BootstrapStage::SlotsToCoeffs,
                factor,
                diagonal_scale,
            )?;
        }
        builder.output_headroom(message_scale)?;
        Ok(Self {
            stages: builder.stages,
            keys: builder.keys,
            eval_mod,
            message_scale,
        })
    }

    /// Scale inputs are lifted to; derived from `q0` and the message ratio.
    pub fn message_scale(&self) -> Scale {
        self.message_scale
    }

    pub fn stages(&self) -> &[StageSchedule] {
        &self.stages
    }

    /// The first stage of the given kind.
    pub fn stage(&self, stage: BootstrapStage) -> Option<&StageSchedule> {
        self.stages.iter().find(|s| s.stage == stage)
    }

    /// Every stage of the given kind, in evaluation order.
    pub fn stages_of(
        &self,
        stage: BootstrapStage,
    ) -> impl Iterator<Item = &StageSchedule> + '_ {
        self.stages.iter().filter(move |s| s.stage == stage)
    }

    pub fn required_keys(&self) -> &BootstrapKeyRequirements {
        &self.keys
    }

    pub fn eval_mod_plan(&self) -> &ChebyshevPlan {
        &self.eval_mod
    }

    fn last(&self) -> &StageSchedule {
        self.stages
            .last()
            .expect("a schedule has a mod-raise stage")
    }

    /// Levels consumed between the raised ciphertext and the output.
    pub fn depth(&self) -> usize {
        self.stages[0].output_level - self.output_level()
    }

    pub fn output_level(&self) -> usize {
        self.last().output_level
    }

    /// Output scale for an input prepared at `prepared` (within half a bit
    /// of the message scale): the inverse transform assumes the message
    /// scale, so the tracked scale is corrected by `prepared / message`.
    ///
    /// # Errors
    /// Rejects unrepresentable scales.
    pub fn output_scale(&self, prepared: Scale) -> CkksResult<Scale> {
        Ok(self
            .last()
            .output_scale
            .checked_mul(prepared)?
            .checked_div(self.message_scale)?)
    }
}

fn check_drift(working: Scale, scale: Scale, level: usize) -> CkksResult<()> {
    let drift = (scale.log2() - working.log2()).abs();
    if drift > MAX_SCALE_DRIFT_BITS {
        return Err(invalid(&format!(
            "bootstrap primes drift the working scale by {drift:.2} bits at level {level}"
        )));
    }
    Ok(())
}

struct Builder<'a, P: PolyRing<N>, const N: usize> {
    top: &'a P,
    stages: Vec<StageSchedule>,
    keys: BootstrapKeyRequirements,
    level: usize,
    scale: Scale,
}

impl<P: PolyRing<N>, const N: usize> Builder<'_, P, N> {
    fn push(
        &mut self,
        stage: BootstrapStage,
        input_level: usize,
        divisors: Vec<u64>,
    ) {
        self.stages.push(StageSchedule {
            stage,
            input_level,
            output_level: self.level,
            output_scale: self.scale,
            plaintext_scale: None,
            divisors,
        });
    }

    /// A product of `scale` and slot magnitude `bound` must fit the modulus
    /// at the current level with a sign bit to spare (the planner's rule).
    fn headroom(&self, scale: Scale, bound: f64) -> CkksResult<()> {
        let bits = self.top.mod_switch_to(self.level)?.modulus_bits();
        if bound.max(1.0).log2() + scale.log2() + 1.0 >= f64::from(bits) {
            return Err(invalid(&format!(
                "bootstrap product lacks headroom at level {}",
                self.level
            )));
        }
        Ok(())
    }

    /// The output must survive one squaring and its rescale: a level to
    /// drop and room for `message_scale^2` (a smaller message ratio raises
    /// the message scale, so this is checked, not assumed).
    fn output_headroom(&self, message_scale: Scale) -> CkksResult<()> {
        let square = message_scale.checked_mul(message_scale)?;
        if self.level < 2 || self.headroom(square, 1.0).is_err() {
            return Err(invalid(&format!(
                "bootstrap output at level {} cannot be squared: needs a \
                 level to rescale and room for a 2^{:.0} product",
                self.level,
                square.log2()
            )));
        }
        Ok(())
    }

    fn rescale(&mut self, stage: BootstrapStage, product: Scale) -> CkksResult<()> {
        if self.level < 2 {
            return Err(invalid("bootstrap chain is too short for its stages"));
        }
        let divisor = self.top.mod_switch_to(self.level)?.rescale_divisor();
        let input_level = self.level;
        // Same conversion as `Ciphertext::rescale`.
        self.scale = product.checked_div(Scale::new(divisor as f64)?)?;
        self.level -= 1;
        self.push(stage, input_level, vec![divisor]);
        Ok(())
    }

    fn transform(
        &mut self,
        stage: BootstrapStage,
        factor: &DiagonalFactor,
        diagonal_scale: Scale,
    ) -> CkksResult<()> {
        let level = self.level;
        self.keys.rotations.extend(
            factor
                .diagonals()
                .keys()
                .filter(|offset| **offset != 0)
                .map(|offset| (level, *offset)),
        );
        let product = self.scale.checked_mul(diagonal_scale)?;
        self.headroom(product, 1.0)?;
        self.rescale(stage, product)?;
        if let Some(last) = self.stages.last_mut() {
            last.plaintext_scale = Some(diagonal_scale);
        }
        Ok(())
    }
}
