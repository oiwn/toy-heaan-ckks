//! Target-scale Chebyshev operation graphs, compiled before evaluation.

use std::collections::{BTreeMap, BTreeSet};

use crate::PolyRing;
use crate::ckks::{CkksError, CkksResult, Scale};

/// A reference to an earlier operation in a polynomial graph.
pub type ValueId = usize;

/// An operation with explicit dependencies. Every dependency precedes its user.
#[derive(Clone, Debug)]
pub enum Operation {
    /// The normalized Chebyshev-domain input.
    Input,
    /// A trivial ciphertext at the node's level and scale.
    Constant(f64),
    /// Drop primes without dividing or changing scale.
    ModSwitch(ValueId),
    /// Multiply and relinearize using a key at the node's level.
    Multiply(ValueId, ValueId),
    /// Compute `2*input - 1` at the same scale.
    DoubleMinusOne(ValueId),
    /// Encode the coefficient at `plaintext_scale`, add the constant at the
    /// resulting product scale. The next operation performs the rescale.
    Linear {
        input: ValueId,
        coefficient: f64,
        constant: f64,
        plaintext_scale: Scale,
    },
    /// Divide by the prime at the source's level.
    Rescale(ValueId),
    /// Add two operands already at compatible levels and scales.
    Add(ValueId, ValueId),
}

/// Planned metadata and a conservative noiseless slot-magnitude bound.
#[derive(Clone, Debug)]
pub struct PlannedValue {
    pub operation: Operation,
    pub level: usize,
    pub scale: Scale,
    pub magnitude_bound: f64,
}

/// A polynomial evaluation graph. Levels use the ring convention: bottom=1.
#[derive(Clone, Debug)]
pub struct ChebyshevPlan {
    values: Vec<PlannedValue>,
    output: ValueId,
    relin_levels: BTreeSet<usize>,
}

impl ChebyshevPlan {
    /// Compile a polynomial with input restricted to [-1,1].
    ///
    /// # Errors
    /// Rejects empty/non-finite coefficients, insufficient levels, scale
    /// overflow, incompatible branches and insufficient noiseless headroom.
    /// Incoming noise remains a separate caller precondition.
    pub fn compile<P: PolyRing<N>, const N: usize>(
        input: &P,
        input_scale: Scale,
        coefficients: &[f64],
        target_scale: Scale,
    ) -> CkksResult<Self> {
        if coefficients.is_empty() || coefficients.iter().any(|c| !c.is_finite()) {
            return Err(invalid(
                "polynomial coefficients must be finite and nonempty",
            ));
        }
        let degree = coefficients.len() - 1;
        let depth = if degree == 0 {
            0
        } else {
            degree.ilog2() as usize + 1
        };
        let output_level = input
            .level()
            .checked_sub(depth)
            .filter(|level| *level > 0)
            .ok_or_else(|| invalid("polynomial has insufficient levels"))?;
        let mut compiler = Compiler {
            ring: input,
            plan: Self {
                values: Vec::new(),
                output: 0,
                relin_levels: BTreeSet::new(),
            },
            powers: BTreeMap::new(),
        };
        let source =
            compiler.push(Operation::Input, input.level(), input_scale, 1.0)?;
        compiler.powers.insert(1, source);
        let output =
            compiler.polynomial(coefficients, output_level, target_scale)?;
        compiler.plan.output = output;
        Ok(compiler.plan)
    }

    pub fn values(&self) -> &[PlannedValue] {
        &self.values
    }
    pub fn output(&self) -> ValueId {
        self.output
    }
    pub fn required_relinearization_levels(&self) -> &BTreeSet<usize> {
        &self.relin_levels
    }

    /// Evaluate graph semantics without encryption or encoding quantization.
    ///
    /// # Errors
    /// Rejects inputs outside the Chebyshev domain [-1,1].
    pub fn evaluate_reference(&self, input: f64) -> CkksResult<f64> {
        if !input.is_finite() || !(-1.0..=1.0).contains(&input) {
            return Err(invalid("polynomial reference input must be in [-1,1]"));
        }
        let mut values = Vec::with_capacity(self.values.len());
        for node in &self.values {
            let value = match node.operation {
                Operation::Input => input,
                Operation::Constant(value) => value,
                Operation::ModSwitch(a) | Operation::Rescale(a) => values[a],
                Operation::Multiply(a, b) => values[a] * values[b],
                Operation::DoubleMinusOne(a) => 2.0 * values[a] - 1.0,
                Operation::Linear {
                    input,
                    coefficient,
                    constant,
                    ..
                } => coefficient * values[input] + constant,
                Operation::Add(a, b) => values[a] + values[b],
            };
            values.push(value);
        }
        Ok(values[self.output])
    }
}

struct Compiler<'a, P> {
    ring: &'a P,
    plan: ChebyshevPlan,
    powers: BTreeMap<usize, ValueId>,
}

impl<P> Compiler<'_, P> {
    fn value(&self, id: ValueId) -> &PlannedValue {
        &self.plan.values[id]
    }
}

impl<P> Compiler<'_, P> {
    fn push<const N: usize>(
        &mut self,
        operation: Operation,
        level: usize,
        scale: Scale,
        bound: f64,
    ) -> CkksResult<ValueId>
    where
        P: PolyRing<N>,
    {
        self.check_headroom(level, scale, bound)?;
        let id = self.plan.values.len();
        self.plan.values.push(PlannedValue {
            operation,
            level,
            scale,
            magnitude_bound: bound,
        });
        Ok(id)
    }

    fn at_level<const N: usize>(
        &mut self,
        id: ValueId,
        level: usize,
    ) -> CkksResult<ValueId>
    where
        P: PolyRing<N>,
    {
        let node = self.value(id).clone();
        if node.level < level {
            return Err(invalid("polynomial branch cannot raise its level"));
        }
        if node.level == level {
            return Ok(id);
        }
        self.push(
            Operation::ModSwitch(id),
            level,
            node.scale,
            node.magnitude_bound,
        )
    }

    fn divisor<const N: usize>(&self, level: usize) -> CkksResult<Scale>
    where
        P: PolyRing<N>,
    {
        Ok(Scale::new(
            self.ring.mod_switch_to(level)?.rescale_divisor() as f64,
        )?)
    }

    fn rescale<const N: usize>(&mut self, source: ValueId) -> CkksResult<ValueId>
    where
        P: PolyRing<N>,
    {
        let node = self.value(source).clone();
        let level = node
            .level
            .checked_sub(1)
            .filter(|v| *v > 0)
            .ok_or_else(|| invalid("polynomial cannot rescale bottom level"))?;
        self.push(
            Operation::Rescale(source),
            level,
            node.scale.checked_div(self.divisor(node.level)?)?,
            node.magnitude_bound,
        )
    }

    fn power<const N: usize>(&mut self, index: usize) -> CkksResult<ValueId>
    where
        P: PolyRing<N>,
    {
        if let Some(id) = self.powers.get(&index) {
            return Ok(*id);
        }
        let half = self.power(index / 2)?;
        let node = self.value(half).clone();
        let square = self.push(
            Operation::Multiply(half, half),
            node.level,
            node.scale.checked_mul(node.scale)?,
            1.0,
        )?;
        self.plan.relin_levels.insert(node.level);
        // The exact function is T_(2k), whose magnitude is <=1 on [-1,1].
        let doubled = self.push(
            Operation::DoubleMinusOne(square),
            node.level,
            self.value(square).scale,
            1.0,
        )?;
        let result = self.rescale(doubled)?;
        self.powers.insert(index, result);
        Ok(result)
    }

    fn check_headroom<const N: usize>(
        &self,
        level: usize,
        scale: Scale,
        bound: f64,
    ) -> CkksResult<()>
    where
        P: PolyRing<N>,
    {
        let bits = self.ring.mod_switch_to(level)?.modulus_bits();
        if !bound.is_finite()
            || bound.max(1.0).log2() + scale.log2() + 1.0 >= f64::from(bits)
        {
            return Err(invalid("polynomial node lacks encoding/product headroom"));
        }
        Ok(())
    }

    fn linear<const N: usize>(
        &mut self,
        coefficients: &[f64],
        level: usize,
        target: Scale,
    ) -> CkksResult<ValueId>
    where
        P: PolyRing<N>,
    {
        let source = self.at_level(0, level + 1)?;
        let plaintext_scale = target
            .checked_mul(self.divisor(level + 1)?)?
            .checked_div(self.value(source).scale)?;
        self.check_headroom(level + 1, plaintext_scale, coefficients[1].abs())?;
        let scale = self.value(source).scale.checked_mul(plaintext_scale)?;
        let bound = coefficients.iter().map(|c| c.abs()).sum();
        let linear = self.push(
            Operation::Linear {
                input: source,
                coefficient: coefficients[1],
                constant: coefficients[0],
                plaintext_scale,
            },
            level + 1,
            scale,
            bound,
        )?;
        self.rescale(linear)
    }

    fn polynomial<const N: usize>(
        &mut self,
        coefficients: &[f64],
        level: usize,
        target: Scale,
    ) -> CkksResult<ValueId>
    where
        P: PolyRing<N>,
    {
        if coefficients.len() == 1 {
            return self.push(
                Operation::Constant(coefficients[0]),
                level,
                target,
                coefficients[0].abs(),
            );
        }
        if coefficients.len() == 2 {
            return self.linear(coefficients, level, target);
        }
        let degree = coefficients.len() - 1;
        let split = 1usize << degree.ilog2();
        let mut quotient = vec![0.0; coefficients.len() - split];
        let mut remainder = coefficients[..split].to_vec();
        quotient[0] = coefficients[split];
        for k in 1..quotient.len() {
            quotient[k] = 2.0 * coefficients[split + k];
            remainder[split - k] -= coefficients[split + k];
        }
        let power = self.power(split)?;
        let power = self.at_level(power, level + 1)?;
        let qscale = target
            .checked_mul(self.divisor(level + 1)?)?
            .checked_div(self.value(power).scale)?;
        let q = self.polynomial(&quotient, level + 1, qscale)?;
        let bound = self.value(q).magnitude_bound;
        let product = self.push(
            Operation::Multiply(power, q),
            level + 1,
            self.value(power).scale.checked_mul(self.value(q).scale)?,
            bound,
        )?;
        self.plan.relin_levels.insert(level + 1);
        let product = self.rescale(product)?;
        let r = self.polynomial(&remainder, level, self.value(product).scale)?;
        if !self.value(r).scale.is_compatible(self.value(product).scale) {
            return Err(invalid("polynomial branches have incompatible scales"));
        }
        self.push(
            Operation::Add(product, r),
            level,
            self.value(product).scale,
            self.value(product).magnitude_bound + self.value(r).magnitude_bound,
        )
    }
}

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::algorithms::chebyshev::eval_chebyshev_clenshaw;
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
    use std::sync::Arc;

    #[test]
    fn invalid_coefficients_levels_and_headroom_return_errors() {
        let context = Arc::new(RnsBasis::<8>::new(vec![97, 113, 193]).unwrap());
        let ring = RnsPoly::zero(&context);
        let scale = Scale::power_of_two(5).unwrap();
        for coefficients in [vec![], vec![f64::NAN], vec![f64::INFINITY]] {
            assert!(
                ChebyshevPlan::compile(&ring, scale, &coefficients, scale).is_err()
            );
        }
        assert!(ChebyshevPlan::compile(&ring, scale, &[1.0; 32], scale).is_err());
        let huge = Scale::power_of_two(100).unwrap();
        assert!(ChebyshevPlan::compile(&ring, huge, &[1.0], scale).is_err());
        let tiny = Scale::new(1e-100).unwrap();
        assert!(ChebyshevPlan::compile(&ring, tiny, &[0.0, 1.0], scale).is_err());
        let plan = ChebyshevPlan::compile(&ring, scale, &[1.0], scale).unwrap();
        assert!(plan.evaluate_reference(1.01).is_err());
        assert!(plan.evaluate_reference(f64::NAN).is_err());
    }

    #[test]
    fn graph_matches_clenshaw_across_degrees_and_mixed_divisors() {
        let mut primes = generate_primes(61, 20, 32);
        primes.extend(generate_primes(60, 10, 32));
        let context = Arc::new(RnsBasis::<32>::new(primes).unwrap());
        let ring = RnsPoly::zero(&context);
        let scale = Scale::power_of_two(61).unwrap();
        for degree in [0, 1, 2, 3, 15, 31, 63] {
            let coefficients: Vec<_> = (0..=degree)
                .map(|j| (j as f64 * 0.7).cos() / 16.0)
                .collect();
            let plan =
                ChebyshevPlan::compile(&ring, scale, &coefficients, scale).unwrap();
            assert!(plan.values()[plan.output()].scale.is_compatible(scale));
            for (id, node) in plan.values().iter().enumerate() {
                match node.operation {
                    Operation::Multiply(a, b) => {
                        assert!(a < id && b < id);
                        assert_eq!(plan.values()[a].level, node.level);
                        assert_eq!(plan.values()[b].level, node.level);
                        let expected = plan.values()[a]
                            .scale
                            .checked_mul(plan.values()[b].scale)
                            .unwrap();
                        assert!(expected.is_compatible(node.scale));
                        assert!(
                            plan.required_relinearization_levels()
                                .contains(&node.level)
                        );
                    }
                    Operation::Rescale(a) => {
                        assert!(a < id);
                        assert_eq!(plan.values()[a].level, node.level + 1);
                        let divisor = ring
                            .mod_switch_to(node.level + 1)
                            .unwrap()
                            .rescale_divisor();
                        let expected = plan.values()[a]
                            .scale
                            .checked_div(Scale::new(divisor as f64).unwrap())
                            .unwrap();
                        assert!(expected.is_compatible(node.scale));
                    }
                    Operation::Add(a, b) => {
                        assert!(a < id && b < id);
                        assert_eq!(plan.values()[a].level, node.level);
                        assert_eq!(plan.values()[b].level, node.level);
                        assert!(
                            plan.values()[a]
                                .scale
                                .is_compatible(plan.values()[b].scale)
                        );
                    }
                    _ => {}
                }
            }
            for j in 0..101 {
                let x = -1.0 + 2.0 * j as f64 / 100.0;
                let expected = eval_chebyshev_clenshaw(x, &coefficients);
                assert!(
                    (plan.evaluate_reference(x).unwrap() - expected).abs() < 1e-10
                );
            }
        }
    }
}
