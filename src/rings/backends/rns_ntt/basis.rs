use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::math::is_ntt_friendly_prime;

use super::errors::{RnsNttError, RnsNttResult};

#[derive(Debug, Clone)]
pub struct NttTable<const DEGREE: usize> {
    pub forward_roots: [u64; DEGREE],
    pub inverse_roots: [u64; DEGREE],
    pub n_inv: u64,
    pub modulus: u64,
    /// `twist_factors[j] = ψ^j mod q` where ψ is the primitive 2N-th root of unity.
    /// Used to pre-twist before the forward NTT so that pointwise multiplication
    /// in NTT domain corresponds to negacyclic convolution in Z[X]/(X^N + 1).
    pub twist_factors: [u64; DEGREE],
    /// `untwist_factors[j] = ψ^{-j} mod q`. Applied after the inverse NTT to
    /// undo the pre-twist.
    pub untwist_factors: [u64; DEGREE],
}

impl<const DEGREE: usize> NttTable<DEGREE> {
    pub fn new(modulus: u64) -> RnsNttResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(RnsNttError::InvalidDegree { degree: DEGREE });
        }
        if !is_ntt_friendly_prime(modulus, DEGREE as u64) {
            return Err(RnsNttError::NonNttFriendlyModulus {
                modulus,
                degree: DEGREE,
            });
        }

        // ψ = primitive 2N-th root of unity (order 2N).
        // ω = ψ² = primitive N-th root of unity (order N).
        //
        // The butterfly twiddle factors use ω so that the NTT computes the
        // standard length-N circular convolution. The pre/post-twist factors
        // use ψ so that `to_ntt_domain`/`to_coeff_domain` implement the
        // negacyclic transform for Z[X]/(X^N + 1).
        let order = 2 * DEGREE;
        let psi = find_primitive_root(modulus, order); // order 2N
        let omega = mod_pow(psi, 2, modulus); // order N
        let _bit_count = DEGREE.trailing_zeros() as usize;

        let omega_inv = mod_inverse(omega, modulus);

        // Natural-power roots: forward_roots[i] = ω^i, inverse_roots[i] = ω^{-i}.
        // The cooley_tukey_ntt accesses roots[offset * step], which equals
        // ω^{offset * (N/len)} — exactly the twiddle needed for the DIT NTT
        // to compute the standard circular DFT at evaluation points ω^0,...,ω^{N-1}.
        let mut forward_roots = [1u64; DEGREE];
        for (i, item) in forward_roots.iter_mut().enumerate().skip(1) {
            *item = mod_pow(omega, i as u64, modulus);
        }

        let mut inverse_roots = [1u64; DEGREE];
        for (i, item) in inverse_roots.iter_mut().enumerate().skip(1) {
            *item = mod_pow(omega_inv, i as u64, modulus);
        }

        let n_inv = mod_inverse(DEGREE as u64, modulus);

        // ψ^j for j = 0..N: pre-twist before forward NTT.
        let psi_inv = mod_inverse(psi, modulus);
        let mut twist_factors = [1u64; DEGREE];
        for (j, item) in twist_factors.iter_mut().enumerate().skip(1) {
            *item = mod_pow(psi, j as u64, modulus);
        }

        // ψ^{-j} for j = 0..N: post-untwist after inverse NTT.
        let mut untwist_factors = [1u64; DEGREE];
        for (j, item) in untwist_factors.iter_mut().enumerate().skip(1) {
            *item = mod_pow(psi_inv, j as u64, modulus);
        }

        Ok(Self {
            forward_roots,
            inverse_roots,
            n_inv,
            modulus,
            twist_factors,
            untwist_factors,
        })
    }
}

/// RNS basis: a set of NTT-friendly prime moduli with precomputed NTT tables.
///
/// Invariant: `moduli.len() == ntt_tables.len()` and
/// `ntt_tables[i].modulus == moduli[i]` for all `i`.
#[derive(Debug, Clone)]
pub struct RnsBasis<const DEGREE: usize> {
    moduli: Vec<u64>,
    ntt_tables: Vec<NttTable<DEGREE>>,
}

impl<const DEGREE: usize> RnsBasis<DEGREE> {
    pub fn new(moduli: Vec<u64>) -> RnsNttResult<Self> {
        if moduli.is_empty() {
            return Err(RnsNttError::EmptyBasis);
        }
        let mut ntt_tables = Vec::with_capacity(moduli.len());
        for &modulus in &moduli {
            ntt_tables.push(NttTable::new(modulus)?);
        }
        Ok(Self { moduli, ntt_tables })
    }

    pub fn moduli(&self) -> &[u64] {
        &self.moduli
    }

    pub fn ntt_table(&self, channel: usize) -> &NttTable<DEGREE> {
        &self.ntt_tables[channel]
    }

    pub fn channel_count(&self) -> usize {
        self.moduli.len()
    }

    /// Returns a new basis with the last `drop_count` channels removed.
    pub fn drop_last(&self, drop_count: usize) -> RnsNttResult<Self> {
        let channel_count = self.channel_count();
        if drop_count >= channel_count {
            return Err(RnsNttError::InvalidModDrop {
                drop_count,
                channel_count,
            });
        }
        let keep = channel_count - drop_count;
        Ok(Self {
            moduli: self.moduli[..keep].to_vec(),
            ntt_tables: self.ntt_tables[..keep].to_vec(),
        })
    }

    /// Returns the total bit-size of the composite modulus Q = q_0 · … · q_{L-1}.
    ///
    /// Computed as `Σ floor(log2(q_i))`, which gives the number of bits needed
    /// to represent Q. Used to set `logq` when calling `CkksEngine::encrypt`.
    pub fn total_bits(&self) -> u32 {
        self.moduli
            .iter()
            .map(|&q| u64::BITS - q.leading_zeros() - 1)
            .sum()
    }

    /// CRT-reconstructs a single coefficient and centers it in `(-Q/2, Q/2]`.
    ///
    /// Given per-channel residues `r[i] = a mod q_i`, recovers `a mod Q`
    /// where `Q = q_0 * … * q_{L-1}`, then maps values above `Q/2` to negative.
    ///
    /// Uses `num_bigint::BigUint` internally, so large RNS chains are not
    /// limited by `u128` during reconstruction. The final centered coefficient
    /// must still fit in `i64`, because the current encoder/decoder interface
    /// returns integer coefficients as `[i64; DEGREE]`.
    ///
    /// # Panics
    /// Panics if the centered reconstructed coefficient does not fit in `i64`.
    pub fn reconstruct_centered_coeff(&self, residues: &[u64]) -> i64 {
        let (negative, magnitude) = self.reconstruct_centered_magnitude(residues);
        if negative {
            let magnitude = biguint_to_i64(&magnitude);
            -magnitude
        } else {
            biguint_to_i64(&magnitude)
        }
    }

    /// Like [`RnsBasis::reconstruct_centered_coeff`] but returns an `f64`,
    /// so it also works when the full modulus `Q` exceeds the `i64` range
    /// (high-level ciphertexts in deep pipelines).
    ///
    /// The returned value is exact for every coefficient that fits in the
    /// 53-bit mantissa of `f64`; larger ones are rounded.
    ///
    /// # Panics
    /// Panics if the centered magnitude does not fit in `f64`.
    pub fn reconstruct_centered_coeff_f64(&self, residues: &[u64]) -> f64 {
        let (negative, magnitude) = self.reconstruct_centered_magnitude(residues);
        let magnitude = magnitude
            .to_f64()
            .expect("reconstructed coefficient exceeds f64 range");
        if negative { -magnitude } else { magnitude }
    }

    /// Mixed-radix CRT reconstruction shared by the `i64` / `f64` variants.
    ///
    /// Returns `(is_negative, |value|)` where `value` is the centered
    /// representative in `(-Q/2, Q/2]` of the residue vector.
    fn reconstruct_centered_magnitude(&self, residues: &[u64]) -> (bool, BigUint) {
        debug_assert_eq!(residues.len(), self.moduli.len());

        let mut value = BigUint::from(residues[0]);
        let mut prefix_product = BigUint::from(self.moduli[0]);

        for (idx, &prime) in self.moduli.iter().enumerate().skip(1) {
            let current_mod = (&value % prime)
                .to_u64()
                .expect("value mod u64 prime always fits in u64");
            let delta = sub_mod(residues[idx], current_mod, prime);
            if delta != 0 {
                let prefix_mod = (&prefix_product % prime)
                    .to_u64()
                    .expect("prefix_product mod u64 prime always fits in u64");
                let prefix_inv = mod_inverse(prefix_mod, prime);
                let adjustment =
                    ((delta as u128 * prefix_inv as u128) % prime as u128) as u64;
                if adjustment != 0 {
                    value += &prefix_product * BigUint::from(adjustment);
                }
            }
            prefix_product *= prime;
        }

        let half_modulus = &prefix_product >> 1usize;

        if value > half_modulus {
            (true, prefix_product - value)
        } else {
            (false, value)
        }
    }
}

fn biguint_to_i64(value: &BigUint) -> i64 {
    value
        .to_i64()
        .expect("reconstructed coefficient exceeds i64 range")
}

// ─── Private number-theory helpers ───────────────────────────────────────────

fn mod_pow(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    base %= modulus;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = ((result as u128 * base as u128) % modulus as u128) as u64;
        }
        base = ((base as u128 * base as u128) % modulus as u128) as u64;
        exponent >>= 1;
    }
    result
}

fn sub_mod(value: u64, subtrahend: u64, modulus: u64) -> u64 {
    if value >= subtrahend {
        (value - subtrahend) % modulus
    } else {
        (value + modulus - subtrahend) % modulus
    }
}

pub(super) fn mod_inverse(value: u64, modulus: u64) -> u64 {
    fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
        if a == 0 {
            (b, 0, 1)
        } else {
            let (gcd, x1, y1) = extended_gcd(b % a, a);
            (gcd, y1 - (b / a) * x1, x1)
        }
    }
    let (gcd, x, _) = extended_gcd(value as i128, modulus as i128);
    assert_eq!(gcd, 1, "mod_inverse: values must be coprime");
    ((x % modulus as i128 + modulus as i128) % modulus as i128) as u64
}

/// Finds a primitive `order`-th root of unity in Z_modulus.
///
/// # Panics
/// Cannot panic when `modulus` is an NTT-friendly prime for the corresponding
/// degree, since such primes are guaranteed to have the required root.
fn find_primitive_root(modulus: u64, order: usize) -> u64 {
    let exponent = (modulus - 1) / order as u64;
    let factors = distinct_prime_factors(order);

    'candidate: for candidate in 2..modulus {
        let root = mod_pow(candidate, exponent, modulus);
        if root == 1 {
            continue;
        }
        for &factor in &factors {
            if mod_pow(root, (order / factor) as u64, modulus) == 1 {
                continue 'candidate;
            }
        }
        return root;
    }

    panic!(
        "find_primitive_root: no root found for modulus {modulus}, order {order}"
    );
}

fn distinct_prime_factors(mut value: usize) -> Vec<usize> {
    let mut factors = Vec::new();
    let mut d = 2usize;
    while d * d <= value {
        if value.is_multiple_of(d) {
            factors.push(d);
            while value.is_multiple_of(d) {
                value /= d;
            }
        }
        d += 1;
    }
    if value > 1 {
        factors.push(value);
    }
    factors
}

pub(super) fn reverse_bits(value: usize, bit_count: usize) -> usize {
    value.reverse_bits() >> (usize::BITS as usize - bit_count)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_ntt_table_for_friendly_prime() {
        let table = NttTable::<8>::new(17).unwrap();
        assert_eq!(table.modulus, 17);
        assert_eq!(table.forward_roots[0], 1);
    }

    #[test]
    fn rejects_non_friendly_modulus() {
        let result = NttTable::<8>::new(19);
        assert!(matches!(
            result,
            Err(RnsNttError::NonNttFriendlyModulus {
                modulus: 19,
                degree: 8
            })
        ));
    }

    #[test]
    fn rejects_empty_basis() {
        assert!(matches!(
            RnsBasis::<8>::new(vec![]),
            Err(RnsNttError::EmptyBasis)
        ));
    }

    #[test]
    fn drop_last_reduces_channel_count() {
        let basis = RnsBasis::<8>::new(vec![17, 97, 113]).unwrap();
        let reduced = basis.drop_last(1).unwrap();
        assert_eq!(reduced.channel_count(), 2);
        assert_eq!(reduced.moduli(), &[17, 97]);
    }

    #[test]
    fn drop_last_all_channels_fails() {
        let basis = RnsBasis::<8>::new(vec![17, 97]).unwrap();
        assert!(matches!(
            basis.drop_last(2),
            Err(RnsNttError::InvalidModDrop { .. })
        ));
    }

    #[test]
    fn reconstruct_centered_single_channel() {
        let basis = RnsBasis::<8>::new(vec![97]).unwrap();
        assert_eq!(basis.reconstruct_centered_coeff(&[5]), 5);
        // 96 ≡ -1 (mod 97)
        assert_eq!(basis.reconstruct_centered_coeff(&[96]), -1);
    }

    #[test]
    fn reconstruct_centered_two_channels() {
        let basis = RnsBasis::<8>::new(vec![17, 97]).unwrap();
        // value = 3: residues [3, 3]
        assert_eq!(basis.reconstruct_centered_coeff(&[3, 3]), 3);
        // value = -7: 17-7=10, 97-7=90
        assert_eq!(basis.reconstruct_centered_coeff(&[10, 90]), -7);
    }

    #[test]
    fn reconstruct_centered_many_channels_exceeds_u128_modulus() {
        let primes = crate::math::generate_primes(30, 5, 8);
        let basis = RnsBasis::<8>::new(primes.clone()).unwrap();

        let positive = 1_234_567i64;
        let positive_residues: Vec<u64> = primes
            .iter()
            .map(|&prime| positive.rem_euclid(prime as i64) as u64)
            .collect();
        assert_eq!(
            basis.reconstruct_centered_coeff(&positive_residues),
            positive
        );

        let negative = -987_654i64;
        let negative_residues: Vec<u64> = primes
            .iter()
            .map(|&prime| negative.rem_euclid(prime as i64) as u64)
            .collect();
        assert_eq!(
            basis.reconstruct_centered_coeff(&negative_residues),
            negative
        );
    }
}
