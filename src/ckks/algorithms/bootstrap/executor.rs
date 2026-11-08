//! Encrypted execution of compiled Chebyshev plans.
//!
//! The executor interprets the same [`Operation`] graph that
//! [`ChebyshevPlan::evaluate_reference`] evaluates in plaintext. After every
//! node it checks the ciphertext's level and scale against the plan, so a
//! divergence is reported at the node where it happens.

use super::polynomial::{ChebyshevPlan, Operation, PlannedValue};
use crate::PolyRing;
use crate::ckks::encoding::ckks_encoder::CkksEncoder;
use crate::ckks::types::{Ciphertext, Plaintext};
use crate::ckks::{CkksError, CkksResult, Scale};

impl ChebyshevPlan {
    /// Evaluate the plan on an encrypted input.
    ///
    /// `input` must match the plan's input node (level and scale).
    /// `multiply(a, b)` returns the relinearized product of two operands at
    /// the same level, without rescaling; the plan rescales explicitly.
    ///
    /// # Errors
    /// Returns an error if a node's actual level or scale differs from the
    /// plan, if a constant does not fit the modulus, or if `multiply` or a
    /// ciphertext operation fails.
    ///
    /// # Panics
    /// Panics if a product scale overflows `f64`; compiled plans bound every
    /// scale well below that.
    pub fn execute<P, F, const N: usize>(
        &self,
        input: &Ciphertext<P, N>,
        mut multiply: F,
    ) -> CkksResult<Ciphertext<P, N>>
    where
        P: PolyRing<N>,
        F: FnMut(
            &Ciphertext<P, N>,
            &Ciphertext<P, N>,
        ) -> CkksResult<Ciphertext<P, N>>,
    {
        let mut values: Vec<Ciphertext<P, N>> =
            Vec::with_capacity(self.values().len());
        for (id, node) in self.values().iter().enumerate() {
            let value = match node.operation {
                Operation::Input => input.clone(),
                Operation::Constant(constant) => {
                    let mut zero = input.mod_switch_to(node.level)?.zero_like();
                    zero.scale = node.scale;
                    zero.add_plain(&constant_plaintext(
                        &zero, constant, node.scale,
                    )?)?
                }
                Operation::ModSwitch(a) => values[a].mod_switch_to(node.level)?,
                Operation::Multiply(a, b) => multiply(&values[a], &values[b])?,
                Operation::DoubleMinusOne(a) => {
                    let doubled = values[a].add(&values[a])?;
                    let one = constant_plaintext(&doubled, 1.0, doubled.scale)?;
                    doubled.sub_plain(&one)?
                }
                Operation::Linear {
                    input: a,
                    coefficient,
                    constant,
                    plaintext_scale,
                } => {
                    let factor = constant_plaintext(
                        &values[a],
                        coefficient,
                        plaintext_scale,
                    )?;
                    let product = values[a].mul_plain(&factor);
                    let offset =
                        constant_plaintext(&product, constant, product.scale)?;
                    product.add_plain(&offset)?
                }
                Operation::Rescale(a) => values[a].rescale()?,
                Operation::Add(a, b) => values[a].add(&values[b])?,
            };
            check_node(id, node, &value)?;
            values.push(value);
        }
        Ok(values.swap_remove(self.output()))
    }
}

/// `value` in every slot, encoded at `scale` in `ct`'s context.
fn constant_plaintext<P: PolyRing<N>, const N: usize>(
    ct: &Ciphertext<P, N>,
    value: f64,
    scale: Scale,
) -> CkksResult<Plaintext<P, N>> {
    Ok(CkksEncoder::<N>::new(scale)
        .try_encode_constant(value, ct.c0.context().clone())?)
}

/// Compare a computed node with its planned level and scale.
fn check_node<P: PolyRing<N>, const N: usize>(
    id: usize,
    planned: &PlannedValue,
    actual: &Ciphertext<P, N>,
) -> CkksResult<()> {
    if actual.level() == planned.level && actual.scale.is_compatible(planned.scale)
    {
        return Ok(());
    }
    Err(CkksError::InvalidParameter {
        message: format!(
            "polynomial node {id} ({:?}): planned level {} scale 2^{:.4}, \
             got level {} scale 2^{:.4}",
            planned.operation,
            planned.level,
            planned.scale.log2(),
            actual.level(),
            actual.scale.log2()
        ),
    })
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use super::*;
    use crate::ckks::algorithms::chebyshev::eval_chebyshev_clenshaw;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::ckks::keys::{RelinearizationKey, SecretKey};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

    const N: usize = 16;
    type Ct = Ciphertext<RnsPoly<N>, N>;

    struct Setup {
        encoder: CkksEncoder<N>,
        secret: SecretKey<RnsPoly<N>, N>,
        input: Ct,
        values: Vec<f64>,
        rng: ChaCha20Rng,
    }

    /// Twelve primes, bottom to top: three 61-bit, three 60-bit, six 61-bit.
    /// The power tree rescales by top primes close to the 2^61 scale (a
    /// mismatch doubles with every squaring); the degree-63 product rescale
    /// into level 5 divides by a 60-bit prime.
    fn setup() -> Setup {
        let wide = generate_primes(61, 9, N as u64);
        let mut primes = wide[..3].to_vec();
        primes.extend(generate_primes(60, 3, N as u64));
        primes.extend(&wide[3..]);
        let basis = Arc::new(RnsBasis::<N>::new(primes).unwrap());
        let scale = Scale::power_of_two(61).unwrap();
        let engine = CkksEngine::<RnsPoly<N>, N>::new(
            basis.clone(),
            CkksParams {
                error_variance: 3.2,
                hamming_weight: N / 2,
                scale,
            },
        );
        let mut rng = ChaCha20Rng::seed_from_u64(40);
        let secret = engine.generate_secret_key(&mut rng).unwrap();
        let public = engine.generate_public_key(&secret, &mut rng).unwrap();
        let encoder = CkksEncoder::new(scale);
        let values: Vec<f64> = (0..N / 2)
            .map(|j| -1.0 + 2.0 * j as f64 / (N / 2 - 1) as f64)
            .collect();
        let input =
            engine.encrypt(&encoder.encode(&values, basis), &public, &mut rng);
        Setup {
            encoder,
            secret,
            input,
            values,
            rng,
        }
    }

    impl Setup {
        fn run(&mut self, plan: &ChebyshevPlan, input: &Ct) -> CkksResult<Ct> {
            let mut keys: HashMap<usize, RelinearizationKey<RnsPoly<N>, N>> =
                HashMap::new();
            let (secret, rng) = (&self.secret, &mut self.rng);
            plan.execute(input, |a, b| {
                let level = a.level();
                let key = keys.entry(level).or_insert_with(|| {
                    let sk = secret.mod_switch_to(level).unwrap();
                    RelinearizationKey::generate(&sk, 3.2f64.sqrt(), rng)
                });
                Ok(a.mul(b, key))
            })
        }

        fn decrypt(&self, ct: &Ct) -> Vec<f64> {
            let sk = self.secret.mod_switch_to(ct.level()).unwrap();
            self.encoder.decode(&CkksEngine::decrypt(ct, &sk))
        }
    }

    #[test]
    fn encrypted_plan_matches_clenshaw_on_mixed_divisors() {
        let mut s = setup();
        let input = s.input.clone();
        for degree in [0usize, 1, 2, 3, 15, 31, 63] {
            let coefficients: Vec<_> = (0..=degree)
                .map(|j| (j as f64 * 0.7).cos() / 16.0)
                .collect();
            let plan = ChebyshevPlan::compile(
                &input.c0,
                input.scale,
                &coefficients,
                input.scale,
            )
            .unwrap();
            let output = s.run(&plan, &input).unwrap();
            let planned = &plan.values()[plan.output()];
            assert_eq!(output.level(), planned.level);
            assert!(output.scale.is_compatible(input.scale));
            let decoded = s.decrypt(&output);
            for (x, got) in s.values.iter().zip(&decoded) {
                let want = eval_chebyshev_clenshaw(*x, &coefficients);
                assert!(
                    (got - want).abs() < 1e-6,
                    "degree {degree}: x={x} got {got} want {want}"
                );
            }
        }
    }

    #[test]
    fn mismatched_input_metadata_is_a_node_error() {
        let mut s = setup();
        let input = s.input.clone();
        let coefficients = [0.1, 0.2, 0.3];
        let plan = ChebyshevPlan::compile(
            &input.c0,
            input.scale,
            &coefficients,
            input.scale,
        )
        .unwrap();

        let mut wrong_scale = input.clone();
        wrong_scale.scale = Scale::power_of_two(60).unwrap();
        let lower = input.mod_switch_to(input.level() - 1).unwrap();
        for bad in [wrong_scale, lower] {
            let err = s.run(&plan, &bad).unwrap_err();
            assert!(err.to_string().contains("polynomial node 0"), "{err}");
        }
    }
}
