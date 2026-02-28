use crate::math::is_ntt_friendly_prime;

use super::errors::{RnsNttError, RnsNttResult};

#[derive(Debug, Clone)]
pub struct NttTable<const DEGREE: usize> {
    pub forward_roots: [u64; DEGREE],
    pub inverse_roots: [u64; DEGREE],
    pub n_inv: u64,
    pub modulus: u64,
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

        let order = 2 * DEGREE;
        let primitive_root = find_primitive_root(modulus, order);
        let bit_count = DEGREE.trailing_zeros() as usize;

        let mut forward_roots = [1u64; DEGREE];
        for (index, root) in forward_roots.iter_mut().enumerate().skip(1) {
            let bit_reversed = reverse_bits(index, bit_count);
            *root = mod_pow(primitive_root, bit_reversed as u64, modulus);
        }

        let primitive_root_inverse = mod_inverse(primitive_root, modulus);
        let mut inverse_roots = [1u64; DEGREE];
        for (index, root) in inverse_roots.iter_mut().enumerate().skip(1) {
            let bit_reversed = reverse_bits(index, bit_count);
            *root = mod_pow(primitive_root_inverse, bit_reversed as u64, modulus);
        }

        let n_inv = mod_inverse(DEGREE as u64, modulus);
        Ok(Self {
            forward_roots,
            inverse_roots,
            n_inv,
            modulus,
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

    /// CRT-reconstructs a single coefficient and centers it in `(-Q/2, Q/2]`.
    ///
    /// Given per-channel residues `r[i] = a mod q_i`, recovers `a mod Q`
    /// where `Q = q_0 * … * q_{L-1}`, then maps values above `Q/2` to negative.
    ///
    /// # Overflow note
    /// Uses u128 arithmetic. The intermediate product `r * (Q/q_i)` is bounded
    /// by Q (since r < q_i and Q/q_i = Q/q_i), so this is exact for Q < 2^128.
    pub fn reconstruct_centered_coeff(&self, residues: &[u64]) -> i64 {
        debug_assert_eq!(residues.len(), self.moduli.len());
        let q: u128 = self.moduli.iter().map(|&m| m as u128).product();
        let mut acc = 0u128;

        for (&r, &m) in residues.iter().zip(&self.moduli) {
            // Q_i = Q / q_i  (exact: q_i divides Q)
            let qi = q / m as u128;
            // r * Q_i < Q  (since r < q_i and Q_i = Q / q_i)
            let r_qi = r as u128 * qi;
            let qi_inv = mod_inverse((qi % m as u128) as u64, m);
            let term = (r_qi % q * qi_inv as u128) % q;
            acc = (acc + term) % q;
        }

        if acc > q / 2 {
            (acc as i128 - q as i128) as i64
        } else {
            acc as i64
        }
    }
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

fn mod_inverse(value: u64, modulus: u64) -> u64 {
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
}
