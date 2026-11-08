//! Encrypted evaluation of factored linear transforms, one level per group.
//!
//! Each group of a [`FactoredTransform`] is a sparse set of cyclic
//! diagonals. A group is evaluated as `sum_k rotate(ct, k) * d_k` followed by
//! one rescale, rotating only by the offsets whose diagonal is nonzero.

use super::schedule::StageSchedule;
use super::transforms::FactoredTransform;
use crate::PolyRing;
use crate::ckks::encoding::ckks_encoder::CkksEncoder;
use crate::ckks::types::{Ciphertext, Plaintext};
use crate::ckks::{CkksError, CkksResult};

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

/// One encoded group: the input level and its nonzero diagonals.
pub struct EncodedGroup<P: PolyRing<N>, const N: usize> {
    pub level: usize,
    /// (positive slot offset, diagonal plaintext), offsets ascending.
    pub diagonals: Vec<(usize, Plaintext<P, N>)>,
}

/// A factored transform encoded at scheduled levels and scales.
pub struct EncodedTransform<P: PolyRing<N>, const N: usize> {
    groups: Vec<EncodedGroup<P, N>>,
}

impl<P: PolyRing<N>, const N: usize> EncodedTransform<P, N> {
    /// Encode group `i` at `stages[i].input_level` and
    /// `stages[i].plaintext_scale`. Group 0's values are multiplied by
    /// `first_factor` and the last group's by `last_factor` (normalization
    /// and restoration fused into the end factors).
    ///
    /// # Errors
    /// Rejects a stage count that differs from the group count, stages
    /// without a plaintext scale, invalid levels and unencodable values.
    pub fn encode(
        top: &P,
        transform: &FactoredTransform,
        stages: &[&StageSchedule],
        first_factor: f64,
        last_factor: f64,
    ) -> CkksResult<Self> {
        let factors = transform.factors();
        if factors.len() != stages.len() {
            return Err(invalid("transform groups do not match the schedule"));
        }
        let last = factors.len() - 1;
        let groups = factors
            .iter()
            .zip(stages)
            .enumerate()
            .map(|(index, (factor, stage))| {
                let mut multiplier = 1.0;
                if index == 0 {
                    multiplier *= first_factor;
                }
                if index == last {
                    multiplier *= last_factor;
                }
                let scale = stage.plaintext_scale.ok_or_else(|| {
                    invalid("transform stage has no plaintext scale")
                })?;
                let context =
                    top.mod_switch_to(stage.input_level)?.context().clone();
                let encoder = CkksEncoder::new(scale);
                let diagonals = factor
                    .diagonals()
                    .iter()
                    .map(|(offset, values)| {
                        let values: Vec<_> =
                            values.iter().map(|v| *v * multiplier).collect();
                        Ok((
                            *offset,
                            encoder.try_encode_complex(&values, context.clone())?,
                        ))
                    })
                    .collect::<CkksResult<_>>()?;
                Ok(EncodedGroup {
                    level: stage.input_level,
                    diagonals,
                })
            })
            .collect::<CkksResult<_>>()?;
        Ok(Self { groups })
    }

    pub fn groups(&self) -> &[EncodedGroup<P, N>] {
        &self.groups
    }

    /// Apply every group in order. `rotate(ct, offset)` must return `ct`
    /// with slots rotated left by `offset` (slot `j` takes slot `j+offset`).
    ///
    /// # Errors
    /// Rejects an input not at the first group's level, and propagates
    /// rotation, scale and rescale errors.
    pub fn evaluate<F>(
        &self,
        input: &Ciphertext<P, N>,
        mut rotate: F,
    ) -> CkksResult<Ciphertext<P, N>>
    where
        F: FnMut(&Ciphertext<P, N>, usize) -> CkksResult<Ciphertext<P, N>>,
    {
        let mut ct = input.clone();
        for group in &self.groups {
            if ct.level() != group.level {
                return Err(invalid(
                    "transform input is not at the scheduled group level",
                ));
            }
            let Some((_, first)) = group.diagonals.first() else {
                return Err(invalid("transform group has no diagonals"));
            };
            let mut sum = ct.zero_like();
            sum.scale = sum.scale.checked_mul(first.scale)?;
            for (offset, diagonal) in &group.diagonals {
                let rotated = if *offset == 0 {
                    ct.clone()
                } else {
                    rotate(&ct, *offset)?
                };
                sum = sum.add(&rotated.mul_plain(diagonal))?;
            }
            ct = sum.rescale()?;
        }
        Ok(ct)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, BTreeSet};
    use std::sync::Arc;

    use num_complex::Complex64;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use super::*;
    use crate::ckks::Scale;
    use crate::ckks::algorithms::bootstrap::schedule::BootstrapStage;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::ckks::keys::KeySwitchingKey;
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

    /// Near the bootstrap working scale, so key-switch noise is negligible;
    /// one spare prime keeps the last product (2^55 * 2^61) below Q.
    const SCALE_BITS: i32 = 55;

    /// Schedule stages for a standalone transform: one rescale per group,
    /// diagonals at 2^61 so the ciphertext scale stays near its input.
    fn stages(top: usize, groups: usize) -> Vec<StageSchedule> {
        let diagonal = Scale::power_of_two(61).unwrap();
        (0..groups)
            .map(|group| StageSchedule {
                stage: BootstrapStage::CoeffsToSlots,
                input_level: top - group,
                output_level: top - group - 1,
                output_scale: diagonal,
                plaintext_scale: Some(diagonal),
                divisors: Vec::new(),
            })
            .collect()
    }

    fn check<const N: usize>(budgets: &[usize]) {
        let slots = N / 2;
        let values: Vec<_> = (0..slots)
            .map(|j| Complex64::from_polar(0.9, j as f64 * 0.37))
            .collect();
        for &budget in budgets {
            for inverse in [false, true] {
                let transform = if inverse {
                    FactoredTransform::coefficients_to_slots(N, budget)
                } else {
                    FactoredTransform::slots_to_coefficients(N, budget)
                }
                .unwrap();
                // Keep intermediate magnitudes near 1 in both directions.
                let normalization = if inverse { 1.0 / slots as f64 } else { 1.0 };
                let mut expected = transform.apply(&values).unwrap();
                for value in &mut expected {
                    *value *= normalization;
                }
                check_one::<N>(&transform, &values, &expected, normalization);
            }
        }
    }

    fn check_one<const N: usize>(
        transform: &FactoredTransform,
        values: &[Complex64],
        expected: &[Complex64],
        normalization: f64,
    ) {
        let groups = transform.factors().len();
        let top = groups + 2;
        let primes = generate_primes(61, top, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).unwrap());
        let scale = Scale::power_of_two(SCALE_BITS).unwrap();
        let engine = CkksEngine::<RnsPoly<N>, N>::new(
            basis.clone(),
            CkksParams {
                error_variance: 3.2,
                hamming_weight: N / 2,
                scale,
            },
        );
        let mut rng = ChaCha20Rng::seed_from_u64(N as u64 + groups as u64);
        let secret = engine.generate_secret_key(&mut rng).unwrap();
        let public = engine.generate_public_key(&secret, &mut rng).unwrap();
        let zero = RnsPoly::zero(&basis);
        let plan_stages = stages(top, groups);
        let stage_refs: Vec<_> = plan_stages.iter().collect();
        let encoded = EncodedTransform::<RnsPoly<N>, N>::encode(
            &zero,
            transform,
            &stage_refs,
            normalization,
            1.0,
        )
        .unwrap();

        let mut keys = BTreeMap::new();
        for group in encoded.groups() {
            let sk = secret.mod_switch_to(group.level).unwrap();
            let offsets: BTreeSet<_> = group
                .diagonals
                .iter()
                .map(|(offset, _)| *offset)
                .filter(|offset| *offset != 0)
                .collect();
            if groups == slots_log2::<N>() {
                assert!(offsets.len() <= 2, "radix-two group has <= 2 rotations");
            }
            for offset in offsets {
                let rotated = sk.poly.rotate_slots(offset as i32);
                keys.insert(
                    (group.level, offset),
                    KeySwitchingKey::generate(
                        &sk,
                        &rotated,
                        3.2f64.sqrt(),
                        &mut rng,
                    ),
                );
            }
        }

        let encoder = CkksEncoder::new(scale);
        let ct = engine.encrypt(
            &encoder.encode_complex(values, basis.clone()),
            &public,
            &mut rng,
        );
        let output = encoded
            .evaluate(&ct, |ct, offset| {
                let key = keys.get(&(ct.level(), offset)).unwrap();
                let exponent =
                    (0..offset).fold(1u64, |e, _| e * 5 % (2 * N) as u64);
                let mut c0 = ct.c0.automorphism(exponent);
                let (extra, c1) = key.switch(&ct.c1.automorphism(exponent));
                c0 += &extra;
                Ok(Ciphertext {
                    c0,
                    c1,
                    scale: ct.scale,
                })
            })
            .unwrap();
        assert_eq!(output.level(), 2);
        let decoded = CkksEncoder::new(output.scale).decode_complex(
            &CkksEngine::decrypt(&output, &secret.mod_switch_to(2).unwrap()),
        );
        for (j, (actual, want)) in decoded.iter().zip(expected).enumerate() {
            assert!(
                (*actual - *want).norm() < 1e-6,
                "N={N} groups={groups} slot={j}: {actual} vs {want}"
            );
        }
    }

    fn slots_log2<const N: usize>() -> usize {
        (N / 2).ilog2() as usize
    }

    #[test]
    fn encrypted_groups_match_plaintext_16() {
        check::<16>(&[1, 2, 3]);
    }

    #[test]
    fn encrypted_groups_match_plaintext_32() {
        check::<32>(&[1, 2, 4]);
    }

    #[test]
    fn encrypted_groups_match_plaintext_64() {
        check::<64>(&[1, 2, 5]);
    }

    #[test]
    fn encrypted_groups_match_plaintext_128() {
        check::<128>(&[1, 2, 6]);
    }

    #[test]
    fn mismatched_level_and_stage_count_are_errors() {
        const N: usize = 16;
        let basis =
            Arc::new(RnsBasis::<N>::new(generate_primes(61, 3, N as u64)).unwrap());
        let zero = RnsPoly::zero(&basis);
        let transform = FactoredTransform::slots_to_coefficients(N, 2).unwrap();
        let scale = Scale::power_of_two(SCALE_BITS).unwrap();
        let one = stages(3, 1);
        assert!(
            EncodedTransform::<RnsPoly<N>, N>::encode(
                &zero,
                &transform,
                &one.iter().collect::<Vec<_>>(),
                1.0,
                1.0
            )
            .is_err()
        );
        let two = stages(3, 2);
        let encoded = EncodedTransform::<RnsPoly<N>, N>::encode(
            &zero,
            &transform,
            &two.iter().collect::<Vec<_>>(),
            1.0,
            1.0,
        )
        .unwrap();
        let low = Ciphertext {
            c0: zero.mod_switch_to(2).unwrap(),
            c1: zero.mod_switch_to(2).unwrap(),
            scale,
        };
        assert!(encoded.evaluate(&low, |ct, _| Ok(ct.clone())).is_err());
    }
}
