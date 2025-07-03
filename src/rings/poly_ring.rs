use super::basis::RnsBasis;
use core::ops::{Add, Mul, Rem};
use std::{iter::IntoIterator, sync::Arc};

/// An RNS-encoded polynomial in ℤ[X]/(X^DEGREE + 1) using const generics.
///
/// Coefficients are stored modulus-by-modulus:
/// - Outer `Vec`: one entry per prime (channel) in the RNS basis.
/// - Inner `[u64; DEGREE]`: the residues of all `DEGREE` slots for that prime.
///
/// Arithmetic (NTT, add, mul, rescale) is performed per channel. The polynomial
/// degree (number of slots) is the compile-time constant `DEGREE`.
#[derive(Debug, Clone)]
pub struct RnsPolyRing<const DEGREE: usize> {
    /// `coefficients[c][i]` = i-th slot residue mod basis.primes[c]
    pub coefficients: Vec<[u64; DEGREE]>,
    /// Shared RNS basis with primes and NTT tables.
    pub basis: Arc<RnsBasis>,
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// Create the zero polynomial: all residues = 0.
    pub fn zero(basis: Arc<RnsBasis>) -> Self {
        let channels = basis.channel_count();
        Self {
            coefficients: vec![[0; DEGREE]; channels],
            basis,
        }
    }

    /// Number of slots (polynomial degree = DEGREE).
    pub fn len(&self) -> usize {
        DEGREE
    }

    /// Number of RNS channels (primes).
    pub fn channels(&self) -> usize {
        self.coefficients.len()
    }

    pub fn from_integer_coeffs(ints: &[i64], basis: Arc<RnsBasis>) -> Self {
        assert_eq!(ints.len(), DEGREE, "Input slice length must match DEGREE");
        let mut channels: Vec<[u64; DEGREE]> =
            Vec::with_capacity(basis.primes().len());
        for &q in basis.primes() {
            let mut arr = [0u64; DEGREE];
            let q_i64 = q as i64;
            for i in 0..DEGREE {
                // ensure non-negative residue
                let v = ((ints[i] % q_i64 + q_i64) % q_i64) as u64;
                arr[i] = v;
            }
            channels.push(arr);
        }
        Self {
            coefficients: channels,
            basis,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::rings::NttTables;

    use super::*;
    use std::sync::Arc;

    /// Build a minimal RNS basis with `channels` identical primes.
    fn dummy_basis<const D: usize>(channels: usize) -> Arc<RnsBasis> {
        let primes = vec![17u64; channels];
        let ntt_tables = NttTables::build_ntt_tables_for_primes(&primes).unwrap();
        Arc::new(RnsBasis { primes, ntt_tables })
    }

    #[test]
    fn zero_has_correct_dimensions() {
        const D: usize = 4;
        let basis = dummy_basis::<D>(3);
        let poly = RnsPolyRing::<D>::zero(basis.clone());
        assert_eq!(poly.len(), D);
        assert_eq!(poly.channels(), 3);
        // All residues are zero
        for channel in &poly.coefficients {
            for &res in channel.iter() {
                assert_eq!(res, 0, "Expected zero residue");
            }
        }
        // Check shared basis Arc count
        assert_eq!(Arc::strong_count(&basis), 2);
    }
}
