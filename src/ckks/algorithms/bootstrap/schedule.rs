//! Bootstrap configuration and the level/scale schedule derived from it.
//!
//! The schedule simulates every stage on the actual modulus chain with the
//! same `Scale` operations the encrypted evaluation performs, so planned and
//! actual scales agree up to `Scale::is_compatible`.

use std::collections::BTreeSet;
use std::f64::consts::TAU;

use super::polynomial::ChebyshevPlan;
use super::transforms::{DiagonalFactor, FactoredTransform};
use crate::PolyRing;
use crate::ckks::algorithms::chebyshev::fit_chebyshev;
use crate::ckks::{CkksError, CkksResult, Scale};

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

/// Numerical settings of a bootstrap.
#[derive(Clone, Debug, PartialEq)]
pub struct BootstrapConfig {
    /// Scale inputs are lifted to, and the scale the inverse transform
    /// assumes for the restored message.
    pub message_scale: Scale,
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
    /// The settings validated at N=16/64 with 18 near-2^61 primes: message
    /// scale 2^44, working scale 2^61, degree 31, four double angles,
    /// one-level transforms and secret weight up to N/2.
    ///
    /// # Errors
    /// Never fails for these constants; the signature follows `Scale`.
    pub fn reference(ring_degree: usize) -> CkksResult<Self> {
        Ok(Self {
            message_scale: Scale::power_of_two(44)?,
            working_scale: Scale::power_of_two(61)?,
            polynomial_degree: 31,
            double_angles: 4,
            coeffs_to_slots_levels: 1,
            slots_to_coeffs_levels: 1,
            max_secret_weight: ring_degree / 2,
        })
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
    /// Rejects chains too short for the stages and unrepresentable scales.
    pub fn derive<P: PolyRing<N>, const N: usize>(
        top: &P,
        config: &BootstrapConfig,
        coeffs_to_slots: &FactoredTransform,
        slots_to_coeffs: &FactoredTransform,
    ) -> CkksResult<Self> {
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
            builder.rescale(BootstrapStage::DoubleAngle, square)?;
        }
        for factor in slots_to_coeffs.factors() {
            builder.transform(
                BootstrapStage::SlotsToCoeffs,
                factor,
                config.message_scale,
            )?;
        }
        Ok(Self {
            stages: builder.stages,
            keys: builder.keys,
            eval_mod,
            message_scale: config.message_scale,
        })
    }

    pub fn stages(&self) -> &[StageSchedule] {
        &self.stages
    }

    /// The first stage of the given kind.
    pub fn stage(&self, stage: BootstrapStage) -> Option<&StageSchedule> {
        self.stages.iter().find(|s| s.stage == stage)
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
            divisors,
        });
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
        self.rescale(stage, product)
    }
}
