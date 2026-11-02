//! Sparse factors of the fully packed canonical CKKS transform.
//!
//! Coefficient halves use bit-reversed order internally. The same permutation
//! is applied to real and imaginary halves; slotwise modular reduction
//! commutes with it, so encoding and decoding cancel the permutation.

use std::collections::BTreeMap;

use num_complex::Complex64;

use crate::ckks::{CkksError, CkksResult};

/// Nonzero cyclic diagonals: output[j] = sum d[j] * input[j + offset].
#[derive(Clone, Debug)]
pub struct DiagonalFactor {
    slots: usize,
    diagonals: BTreeMap<usize, Vec<Complex64>>,
}

impl DiagonalFactor {
    /// Nonzero diagonals indexed by positive slot rotation.
    pub fn diagonals(&self) -> &BTreeMap<usize, Vec<Complex64>> {
        &self.diagonals
    }

    /// Apply a factor to a plaintext reference vector.
    ///
    /// # Errors
    /// Rejects a vector with a different slot count.
    pub fn apply(&self, input: &[Complex64]) -> CkksResult<Vec<Complex64>> {
        if input.len() != self.slots {
            return Err(invalid("transform input has wrong slot count"));
        }
        let mut output = vec![Complex64::default(); self.slots];
        for (offset, diagonal) in &self.diagonals {
            for (j, value) in output.iter_mut().enumerate() {
                *value += diagonal[j] * input[(j + offset) % self.slots];
            }
        }
        Ok(output)
    }

    fn identity(slots: usize) -> Self {
        Self {
            slots,
            diagonals: BTreeMap::from([(0, vec![Complex64::new(1.0, 0.0); slots])]),
        }
    }

    // Compose next(self(input)) directly in diagonal form, never as a matrix.
    fn followed_by(&self, next: &Self) -> Self {
        let mut diagonals: BTreeMap<usize, Vec<Complex64>> = BTreeMap::new();
        for (left_offset, left) in &next.diagonals {
            for (right_offset, right) in &self.diagonals {
                let offset = (left_offset + right_offset) % self.slots;
                let diagonal = diagonals
                    .entry(offset)
                    .or_insert_with(|| vec![Complex64::default(); self.slots]);
                for j in 0..self.slots {
                    diagonal[j] += left[j] * right[(j + left_offset) % self.slots];
                }
            }
        }
        diagonals
            .retain(|_, values| values.iter().any(|v| *v != Complex64::default()));
        Self {
            slots: self.slots,
            diagonals,
        }
    }

    fn adjoint_halved(&self) -> Self {
        let diagonals = self
            .diagonals
            .iter()
            .map(|(offset, values)| {
                let reverse = (self.slots - offset) % self.slots;
                let values = (0..self.slots)
                    .map(|j| values[(j + reverse) % self.slots].conj() * 0.5)
                    .collect();
                (reverse, values)
            })
            .collect();
        Self {
            slots: self.slots,
            diagonals,
        }
    }
}

/// Canonical transform compiled into a configurable number of sparse groups.
#[derive(Clone, Debug)]
pub struct FactoredTransform {
    slots: usize,
    factors: Vec<DiagonalFactor>,
}

impl FactoredTransform {
    /// Map slots to bit-reversed coefficient halves.
    ///
    /// # Errors
    /// Requires a power-of-two degree >=4 and a budget in 1..=log2(N/2).
    pub fn slots_to_coefficients(
        degree: usize,
        level_budget: usize,
    ) -> CkksResult<Self> {
        let stages = butterflies(degree)?;
        let stages: Vec<_> = stages
            .iter()
            .rev()
            .map(DiagonalFactor::adjoint_halved)
            .collect();
        Self::group(stages, degree / 2, level_budget)
    }

    /// Map bit-reversed coefficient halves back to slots.
    ///
    /// # Errors
    /// Requires a power-of-two degree >=4 and a budget in 1..=log2(N/2).
    pub fn coefficients_to_slots(
        degree: usize,
        level_budget: usize,
    ) -> CkksResult<Self> {
        Self::group(butterflies(degree)?, degree / 2, level_budget)
    }

    fn group(
        stages: Vec<DiagonalFactor>,
        slots: usize,
        budget: usize,
    ) -> CkksResult<Self> {
        if budget == 0 || budget > stages.len() {
            return Err(invalid(
                "transform level budget is outside 1..=log2(slots)",
            ));
        }
        let mut factors = Vec::with_capacity(budget);
        let mut stages = stages.into_iter();
        let count = slots.ilog2() as usize;
        for group in 0..budget {
            let size = count / budget + usize::from(group < count % budget);
            let mut factor = DiagonalFactor::identity(slots);
            for stage in stages.by_ref().take(size) {
                factor = factor.followed_by(&stage);
            }
            factors.push(factor);
        }
        Ok(Self { slots, factors })
    }

    /// Ordered factors, each consuming one rescale level when encrypted.
    pub fn factors(&self) -> &[DiagonalFactor] {
        &self.factors
    }

    /// Total number of stored complex entries across all diagonals.
    pub fn entry_count(&self) -> usize {
        self.factors
            .iter()
            .map(|f| f.diagonals.len() * self.slots)
            .sum()
    }

    /// Apply the factors to a plaintext reference vector.
    ///
    /// # Errors
    /// Rejects a vector with a different slot count.
    pub fn apply(&self, input: &[Complex64]) -> CkksResult<Vec<Complex64>> {
        if input.len() != self.slots {
            return Err(invalid("transform input has wrong slot count"));
        }
        let mut output = input.to_vec();
        for factor in &self.factors {
            output = factor.apply(&output)?;
        }
        Ok(output)
    }
}

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

fn butterflies(degree: usize) -> CkksResult<Vec<DiagonalFactor>> {
    if degree < 4 || !degree.is_power_of_two() || degree.checked_mul(2).is_none() {
        return Err(invalid("transform degree must be a power of two >=4"));
    }
    let slots = degree / 2;
    let mut factors = Vec::with_capacity(slots.ilog2() as usize);
    for stage in 1..=slots.ilog2() {
        let width = 1usize << stage;
        let half = width / 2;
        let mut main = vec![Complex64::default(); slots];
        let mut positive = main.clone();
        let mut negative = main.clone();
        let mut exponent = 1usize;
        for j in 0..half {
            let angle =
                std::f64::consts::TAU * exponent as f64 / (4 * width) as f64;
            let root = Complex64::from_polar(1.0, angle);
            for base in (0..slots).step_by(width) {
                main[base + j] = Complex64::new(1.0, 0.0);
                main[base + j + half] = -root;
                positive[base + j] = root;
                negative[base + j + half] = Complex64::new(1.0, 0.0);
            }
            exponent = ((exponent as u128 * 5) % (4 * width) as u128) as usize;
        }
        let mut diagonals = BTreeMap::from([(0, main), (half, positive)]);
        let reverse = (slots - half) % slots;
        let diagonal = diagonals
            .entry(reverse)
            .or_insert_with(|| vec![Complex64::default(); slots]);
        for (left, right) in diagonal.iter_mut().zip(negative) {
            *left += right;
        }
        factors.push(DiagonalFactor { slots, diagonals });
    }
    Ok(factors)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_dense_reference_and_roundtrip() {
        for degree in [16usize, 32, 64, 128] {
            let slots = degree / 2;
            let input: Vec<_> = (0..slots)
                .map(|j| {
                    Complex64::new(j as f64 / slots as f64, (j as f64 * 0.7).sin())
                })
                .collect();
            let mut expected = vec![Complex64::default(); slots];
            let mut exponent = 1;
            for value in &mut expected {
                for k in 0..slots {
                    let reverse = k.reverse_bits() >> (usize::BITS - slots.ilog2());
                    *value += input[reverse]
                        * Complex64::from_polar(
                            1.0,
                            std::f64::consts::PI * (exponent * k) as f64
                                / degree as f64,
                        );
                }
                exponent = exponent * 5 % (2 * degree);
            }
            for budget in 1..=slots.ilog2() as usize {
                let forward =
                    FactoredTransform::coefficients_to_slots(degree, budget)
                        .unwrap();
                let inverse =
                    FactoredTransform::slots_to_coefficients(degree, budget)
                        .unwrap();
                let actual = forward.apply(&input).unwrap();
                for (a, b) in actual.iter().zip(&expected) {
                    assert!((*a - *b).norm() < 1e-10);
                }
                let restored = inverse.apply(&actual).unwrap();
                for (a, b) in restored.iter().zip(&input) {
                    assert!((*a - *b).norm() < 1e-10);
                }
            }
        }
    }

    #[test]
    fn large_transform_has_linear_storage_per_stage() {
        let degree = 1024;
        let depth = (degree / 2usize).ilog2() as usize;
        let forward =
            FactoredTransform::coefficients_to_slots(degree, depth).unwrap();
        let inverse =
            FactoredTransform::slots_to_coefficients(degree, depth).unwrap();
        assert!(forward.entry_count() <= 3 * degree / 2 * depth);
        let input: Vec<_> = (0..degree / 2)
            .map(|j| Complex64::new((j as f64).sin(), (j as f64 * 0.3).cos()))
            .collect();
        let actual = inverse.apply(&forward.apply(&input).unwrap()).unwrap();
        for (a, b) in actual.iter().zip(input) {
            assert!((*a - b).norm() < 1e-10);
        }
    }

    #[test]
    fn invalid_configuration_returns_errors() {
        for degree in [0, 1, 2, 3, 15, usize::MAX] {
            assert!(FactoredTransform::coefficients_to_slots(degree, 1).is_err());
        }
        assert!(FactoredTransform::coefficients_to_slots(16, 0).is_err());
        assert!(FactoredTransform::coefficients_to_slots(16, 4).is_err());
        assert!(
            FactoredTransform::coefficients_to_slots(16, 1)
                .unwrap()
                .apply(&[])
                .is_err()
        );
    }
}
