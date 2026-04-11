use super::{
    basis::{NttTable, RnsBasis, mod_inverse, reverse_bits},
    errors::{RnsNttError, RnsNttResult},
};
use crate::{
    math::sampling::{ternary_coefficients, uniform_coefficients},
    rings::traits::{PolyRing, PolySampler},
};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::{
    ops::{AddAssign, MulAssign, Neg},
    sync::Arc,
};

/// A polynomial in `Z_{q_0} x … x Z_{q_{L-1}}[X] / (X^N + 1)`.
///
/// Stores one `[u64; DEGREE]` array per RNS channel. The `in_ntt_domain` flag
/// tracks whether the arrays hold coefficient-domain or NTT-domain values.
///
/// # Invariants
/// - `channels.len() == basis.channel_count()`
/// - Every `channels[i][j] < basis.moduli()[i]`
/// - `DEGREE` is a power of two
#[derive(Clone, Debug)]
pub struct RnsPoly<const DEGREE: usize> {
    channels: Vec<[u64; DEGREE]>,
    basis: Arc<RnsBasis<DEGREE>>,
    in_ntt_domain: bool,
}

// ─── Constructors ─────────────────────────────────────────────────────────────

impl<const DEGREE: usize> RnsPoly<DEGREE> {
    /// Creates the zero polynomial in coefficient domain.
    pub fn zero(basis: Arc<RnsBasis<DEGREE>>) -> Self {
        let channels = vec![[0u64; DEGREE]; basis.channel_count()];
        Self {
            channels,
            basis,
            in_ntt_domain: false,
        }
    }

    /// Creates a polynomial from signed integer coefficients.
    ///
    /// Each coefficient is reduced into `[0, q_i)` per channel via `rem_euclid`.
    /// Accepts slices of length ≥ DEGREE; only the first DEGREE elements are used.
    pub fn from_coeffs(coeffs: &[i64], basis: Arc<RnsBasis<DEGREE>>) -> Self {
        assert!(
            coeffs.len() >= DEGREE,
            "from_coeffs: need at least {DEGREE} coefficients, got {}",
            coeffs.len()
        );
        let mut channels = vec![[0u64; DEGREE]; basis.channel_count()];
        for (ch, channel) in channels.iter_mut().enumerate() {
            let q = basis.moduli()[ch] as i128;
            for (i, coeff) in coeffs.iter().take(DEGREE).enumerate() {
                channel[i] = (*coeff as i128).rem_euclid(q) as u64;
            }
        }
        Self {
            channels,
            basis,
            in_ntt_domain: false,
        }
    }

    /// Creates a polynomial from pre-built channel arrays.
    ///
    /// Returns an error if the channel count doesn't match the basis, or if any
    /// coefficient is not reduced (i.e., ≥ the corresponding modulus).
    pub fn from_channels(
        channels: Vec<[u64; DEGREE]>,
        basis: Arc<RnsBasis<DEGREE>>,
        in_ntt_domain: bool,
    ) -> RnsNttResult<Self> {
        let expected = basis.channel_count();
        let actual = channels.len();
        if actual != expected {
            return Err(RnsNttError::ChannelCountMismatch { expected, actual });
        }
        for (ch, channel) in channels.iter().enumerate() {
            let q = basis.moduli()[ch];
            for &c in channel {
                if c >= q {
                    return Err(RnsNttError::NonReducedCoefficient {
                        coefficient: c,
                        modulus: q,
                    });
                }
            }
        }
        Ok(Self {
            channels,
            basis,
            in_ntt_domain,
        })
    }

    // Internal: skips the O(N·L) reducedness check. Only for use where
    // coefficients are constructed correctly by this module.
    fn new_unchecked(
        channels: Vec<[u64; DEGREE]>,
        basis: Arc<RnsBasis<DEGREE>>,
        in_ntt_domain: bool,
    ) -> Self {
        Self {
            channels,
            basis,
            in_ntt_domain,
        }
    }
}

// ─── Accessors & domain conversion ───────────────────────────────────────────

impl<const DEGREE: usize> RnsPoly<DEGREE> {
    pub fn channels(&self) -> &[[u64; DEGREE]] {
        &self.channels
    }

    pub fn basis(&self) -> &Arc<RnsBasis<DEGREE>> {
        &self.basis
    }

    pub fn is_ntt_domain(&self) -> bool {
        self.in_ntt_domain
    }

    /// Converts to NTT domain in-place (no-op if already there).
    ///
    /// Applies a pre-twist (multiply coefficient `j` by ψ^j) before the
    /// forward NTT, so that pointwise multiplication in this domain equals
    /// negacyclic convolution in Z[X]/(X^N + 1).
    pub fn to_ntt_domain(&mut self) {
        if self.in_ntt_domain {
            return;
        }
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            let table = self.basis.ntt_table(ch);
            for (j, v) in channel.iter_mut().enumerate() {
                *v = mul_mod(*v, table.twist_factors[j], table.modulus);
            }
            forward_ntt(channel, table);
        }
        self.in_ntt_domain = true;
    }

    /// Converts to coefficient domain in-place (no-op if already there).
    ///
    /// Applies the inverse NTT followed by a post-untwist (multiply coefficient
    /// `j` by ψ^{-j}), undoing the pre-twist applied in `to_ntt_domain`.
    pub fn to_coeff_domain(&mut self) {
        if !self.in_ntt_domain {
            return;
        }
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            let table = self.basis.ntt_table(ch);
            inverse_ntt(channel, table);
            for (j, v) in channel.iter_mut().enumerate() {
                *v = mul_mod(*v, table.untwist_factors[j], table.modulus);
            }
        }
        self.in_ntt_domain = false;
    }

    /// Returns a new polynomial with the last `drop_count` RNS channels removed.
    pub fn mod_drop_last(&self, drop_count: usize) -> RnsNttResult<Self> {
        let reduced_basis = Arc::new(self.basis.drop_last(drop_count)?);
        let keep = reduced_basis.channel_count();
        Ok(Self::new_unchecked(
            self.channels[..keep].to_vec(),
            reduced_basis,
            self.in_ntt_domain,
        ))
    }

    /// Divides every coefficient by the last RNS prime `q_L` and drops that channel.
    ///
    /// Uses a caller-supplied `new_basis` (the reduced basis with `q_L` dropped) so that
    /// multiple polynomials (e.g. both ciphertext components) can share the same
    /// `Arc<RnsBasis>` instance — required for operations that later `Arc::ptr_eq`-check
    /// their bases.
    ///
    /// See [`rescale`] for a convenience wrapper that creates its own `Arc`.
    pub fn rescale_into(&self, new_basis: Arc<RnsBasis<DEGREE>>) -> RnsNttResult<Self> {
        let num_channels = self.basis.channel_count();
        if num_channels < 2 {
            return Err(RnsNttError::InvalidModDrop {
                drop_count: 1,
                channel_count: num_channels,
            });
        }

        // Work in coefficient domain without mutating self.
        let tmp;
        let channels: &[[u64; DEGREE]] = if self.in_ntt_domain {
            tmp = {
                let mut clone = self.clone();
                clone.to_coeff_domain();
                clone
            };
            &tmp.channels
        } else {
            &self.channels
        };

        let last = num_channels - 1;
        let q_last = self.basis.moduli()[last];

        let mut new_channels = vec![[0u64; DEGREE]; last];
        for i in 0..last {
            let q_i = self.basis.moduli()[i];
            let q_last_inv = mod_inverse(q_last % q_i, q_i);
            for j in 0..DEGREE {
                let c_i = channels[i][j];
                let c_l = channels[last][j] % q_i;
                let diff = sub_mod(c_i, c_l, q_i);
                new_channels[i][j] = mul_mod(diff, q_last_inv, q_i);
            }
        }

        Ok(Self::new_unchecked(new_channels, new_basis, false))
    }

    /// Divides every coefficient by the last RNS prime `q_L` and drops that channel.
    ///
    /// Used after homomorphic multiplication to bring `logp` back to its pre-multiplication
    /// value. The operation is *exact* (no rounding error): subtracting `c mod q_L` from `c`
    /// leaves a value exactly divisible by `q_L`.
    ///
    /// # Formula (per remaining channel `i`, per coefficient `j`)
    /// ```text
    /// new[i][j] = (c[i][j] − c[L][j] mod q_i) · q_L⁻¹  mod q_i
    /// ```
    ///
    /// The polynomial is converted to coefficient domain before rescaling; the result is
    /// always returned in coefficient domain.
    ///
    /// When rescaling multiple polynomials that will later be used together (e.g. the two
    /// ciphertext components), prefer [`rescale_into`] with a shared `Arc<RnsBasis>`.
    pub fn rescale(&self) -> RnsNttResult<Self> {
        let new_basis = Arc::new(self.basis.drop_last(1)?);
        self.rescale_into(new_basis)
    }
}

// ─── Arithmetic ───────────────────────────────────────────────────────────────

impl<const DEGREE: usize> AddAssign<&RnsPoly<DEGREE>> for RnsPoly<DEGREE> {
    /// Coefficient-wise addition modulo each `q_i`.
    ///
    /// Works in both coefficient and NTT domain. Both operands must share the
    /// same basis and be in the same domain.
    fn add_assign(&mut self, rhs: &RnsPoly<DEGREE>) {
        debug_assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "add_assign: basis mismatch"
        );
        debug_assert_eq!(
            self.in_ntt_domain, rhs.in_ntt_domain,
            "add_assign: domain mismatch"
        );
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            let q = self.basis.moduli()[ch];
            for (a, &b) in channel.iter_mut().zip(rhs.channels[ch].iter()) {
                *a = add_mod(*a, b, q);
            }
        }
    }
}

impl<const DEGREE: usize> MulAssign<&RnsPoly<DEGREE>> for RnsPoly<DEGREE> {
    /// Polynomial multiplication in `Z[X]/(X^N + 1)` via NTT, per channel.
    ///
    /// Both operands must share the same basis and be in the **same domain**:
    ///
    /// - **NTT domain**: pointwise multiply per channel; result stays in NTT
    ///   domain. O(N·L). Use this path when chaining multiple multiplications.
    /// - **Coefficient domain**: forward-NTT `self` in place, forward-NTT a
    ///   temporary copy of `rhs`, pointwise multiply, then inverse-NTT `self`
    ///   back to coefficient domain. O(N·log(N)·L).
    fn mul_assign(&mut self, rhs: &RnsPoly<DEGREE>) {
        debug_assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "mul_assign: basis mismatch"
        );
        debug_assert_eq!(
            self.in_ntt_domain, rhs.in_ntt_domain,
            "mul_assign: domain mismatch (both must be in the same domain)"
        );

        if self.in_ntt_domain {
            // Both already in NTT domain: O(N·L) pointwise multiply.
            for ch in 0..self.basis.channel_count() {
                let q = self.basis.moduli()[ch];
                for (a, &b) in
                    self.channels[ch].iter_mut().zip(rhs.channels[ch].iter())
                {
                    *a = mul_mod(*a, b, q);
                }
            }
        } else {
            // Coefficient domain: fwd-NTT self, fwd-NTT rhs clone, pointwise
            // multiply, inv-NTT self back to coefficient domain.
            self.to_ntt_domain();

            let mut rhs_ntt = rhs.channels.clone();
            for (ch, channel) in rhs_ntt.iter_mut().enumerate() {
                let table = self.basis.ntt_table(ch);
                for (j, v) in channel.iter_mut().enumerate() {
                    *v = mul_mod(*v, table.twist_factors[j], table.modulus);
                }
                forward_ntt(channel, table);
            }

            for ch in 0..self.basis.channel_count() {
                let q = self.basis.moduli()[ch];
                for (a, &b) in
                    self.channels[ch].iter_mut().zip(rhs_ntt[ch].iter())
                {
                    *a = mul_mod(*a, b, q);
                }
            }

            self.to_coeff_domain();
        }
    }
}

impl<const DEGREE: usize> RnsPoly<DEGREE> {
    /// Schoolbook O(N²) polynomial multiplication in `Z[X]/(X^N + 1)`.
    ///
    /// Both operands must be in coefficient domain and share the same basis.
    /// Kept as a reference implementation for correctness testing against the
    /// NTT-based path.
    pub fn mul_assign_naive(&mut self, rhs: &RnsPoly<DEGREE>) {
        debug_assert!(
            !self.in_ntt_domain && !rhs.in_ntt_domain,
            "mul_assign_naive: requires coefficient domain"
        );
        debug_assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "mul_assign_naive: basis mismatch"
        );
        for ch in 0..self.basis.channel_count() {
            let q = self.basis.moduli()[ch];
            let lhs = self.channels[ch];
            let rhs_ch = rhs.channels[ch];
            let mut result = [0u64; DEGREE];
            for i in 0..DEGREE {
                for j in 0..DEGREE {
                    let prod = mul_mod(lhs[i], rhs_ch[j], q);
                    if i + j < DEGREE {
                        result[i + j] = add_mod(result[i + j], prod, q);
                    } else {
                        // X^N = -1 in Z[X]/(X^N + 1)
                        result[i + j - DEGREE] =
                            sub_mod(result[i + j - DEGREE], prod, q);
                    }
                }
            }
            self.channels[ch] = result;
        }
    }
}

impl<const DEGREE: usize> Neg for RnsPoly<DEGREE> {
    type Output = Self;

    /// Coefficient-wise negation modulo each `q_i`. Works in both domains.
    fn neg(mut self) -> Self {
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            let q = self.basis.moduli()[ch];
            for c in channel.iter_mut() {
                if *c != 0 {
                    *c = q - *c;
                }
            }
        }
        self
    }
}

// ─── PolyRing trait ───────────────────────────────────────────────────────────

impl<const DEGREE: usize> PolyRing<DEGREE> for RnsPoly<DEGREE> {
    type Context = Arc<RnsBasis<DEGREE>>;

    fn zero(context: &Self::Context) -> Self {
        Self::zero(context.clone())
    }

    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        Self::from_coeffs(coeffs, context.clone())
    }

    /// CRT-reconstructs each coefficient and returns them centered in `(-Q/2, Q/2]`.
    ///
    /// If the polynomial is in NTT domain, a temporary clone is converted to
    /// coefficient domain first.
    fn to_coeffs(&self) -> [i64; DEGREE] {
        // Ensure we operate on coefficient-domain data without mutating self.
        let tmp;
        let channels: &[[u64; DEGREE]] = if self.in_ntt_domain {
            tmp = {
                let mut clone = self.clone();
                clone.to_coeff_domain();
                clone
            };
            &tmp.channels
        } else {
            &self.channels
        };

        let mut result = [0i64; DEGREE];
        let mut residues = vec![0u64; self.basis.channel_count()];
        for i in 0..DEGREE {
            for (ch, channel) in channels.iter().enumerate() {
                residues[ch] = channel[i];
            }
            result[i] = self.basis.reconstruct_centered_coeff(&residues);
        }
        result
    }

    fn context(&self) -> &Self::Context {
        &self.basis
    }
}

// ─── PolySampler trait ────────────────────────────────────────────────────────

impl<const DEGREE: usize> PolySampler<DEGREE> for RnsPoly<DEGREE> {
    /// Samples with coefficients uniform in `[0, q_i)` per channel.
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self {
        let mut channels = vec![[0u64; DEGREE]; context.channel_count()];
        for (ch, channel) in channels.iter_mut().enumerate() {
            *channel = uniform_coefficients::<DEGREE, _>(context.moduli()[ch], rng);
        }
        Self::new_unchecked(channels, context.clone(), false)
    }

    /// Samples noise from N(0, std_dev), rounded and CRT-encoded per channel.
    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let normal = Normal::new(0.0, std_dev)
            .expect("sample_gaussian: std_dev must be finite and positive");
        let mut noise = [0i64; DEGREE];
        for n in noise.iter_mut() {
            *n = normal.sample(rng).round() as i64;
        }
        Self::from_coeffs(&noise, context.clone())
    }

    /// Samples a ternary polynomial with exactly `hamming_weight` non-zero coefficients.
    fn sample_tribits<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let ternary = ternary_coefficients::<DEGREE, _>(hamming_weight, rng);
        Self::from_coeffs(&ternary, context.clone())
    }

    fn sample_noise<R: Rng>(
        variance: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        Self::sample_gaussian(variance.sqrt(), context, rng)
    }
}

// ─── NTT kernels (private) ────────────────────────────────────────────────────

fn forward_ntt<const DEGREE: usize>(
    coeffs: &mut [u64; DEGREE],
    table: &NttTable<DEGREE>,
) {
    bit_reverse_permute(coeffs);
    cooley_tukey_ntt(coeffs, &table.forward_roots, table.modulus);
}

fn inverse_ntt<const DEGREE: usize>(
    evals: &mut [u64; DEGREE],
    table: &NttTable<DEGREE>,
) {
    bit_reverse_permute(evals);
    cooley_tukey_ntt(evals, &table.inverse_roots, table.modulus);
    for v in evals.iter_mut() {
        *v = mul_mod(*v, table.n_inv, table.modulus);
    }
}

fn cooley_tukey_ntt<const DEGREE: usize>(
    values: &mut [u64; DEGREE],
    roots: &[u64; DEGREE],
    modulus: u64,
) {
    let mut len = 2;
    while len <= DEGREE {
        let half = len / 2;
        let step = DEGREE / len;
        for start in (0..DEGREE).step_by(len) {
            for offset in 0..half {
                let left = start + offset;
                let right = left + half;
                let twiddle = roots[offset * step];
                let t = mul_mod(values[right], twiddle, modulus);
                let u = values[left];
                values[left] = add_mod(u, t, modulus);
                values[right] = sub_mod(u, t, modulus);
            }
        }
        len *= 2;
    }
}

fn bit_reverse_permute<const DEGREE: usize>(values: &mut [u64; DEGREE]) {
    let bits = DEGREE.trailing_zeros() as usize;
    for i in 0..DEGREE {
        let j = reverse_bits(i, bits);
        if i < j {
            values.swap(i, j);
        }
    }
}

// ─── Modular arithmetic helpers (private) ─────────────────────────────────────

fn add_mod(a: u64, b: u64, q: u64) -> u64 {
    let s = a + b;
    if s >= q { s - q } else { s }
}

fn sub_mod(a: u64, b: u64, q: u64) -> u64 {
    if a >= b { a - b } else { a + q - b }
}

fn mul_mod(a: u64, b: u64, q: u64) -> u64 {
    ((a as u128 * b as u128) % q as u128) as u64
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    fn basis_17_97() -> Arc<RnsBasis<8>> {
        Arc::new(RnsBasis::new(vec![17, 97]).unwrap())
    }

    fn basis_three() -> Arc<RnsBasis<8>> {
        Arc::new(RnsBasis::new(vec![17, 97, 113]).unwrap())
    }

    // ── Construction ──────────────────────────────────────────────────────────

    #[test]
    fn zero_poly_is_all_zeros() {
        let poly = RnsPoly::<8>::zero(basis_17_97());
        for ch in poly.channels() {
            assert!(ch.iter().all(|&c| c == 0));
        }
        assert!(!poly.is_ntt_domain());
    }

    #[test]
    fn from_coeffs_reduces_correctly() {
        let basis = basis_17_97();
        let poly = RnsPoly::<8>::from_coeffs(&[-1, 2, -3, 4, 0, 0, 0, 0], basis);
        // channel 0: mod 17 → [-1→16, 2→2, -3→14, 4→4, …]
        assert_eq!(poly.channels()[0][0], 16);
        assert_eq!(poly.channels()[0][1], 2);
        assert_eq!(poly.channels()[0][2], 14);
        // channel 1: mod 97 → [-1→96, …]
        assert_eq!(poly.channels()[1][0], 96);
    }

    #[test]
    fn from_channels_rejects_unreduced_coefficient() {
        let basis = basis_17_97();
        let bad = vec![[17u64; 8], [0u64; 8]]; // 17 ≥ q_0 = 17
        assert!(matches!(
            RnsPoly::from_channels(bad, basis, false),
            Err(RnsNttError::NonReducedCoefficient { .. })
        ));
    }

    #[test]
    fn from_channels_rejects_wrong_channel_count() {
        let basis = basis_17_97();
        let too_few = vec![[0u64; 8]]; // need 2
        assert!(matches!(
            RnsPoly::from_channels(too_few, basis, false),
            Err(RnsNttError::ChannelCountMismatch { .. })
        ));
    }

    // ── NTT roundtrip ─────────────────────────────────────────────────────────

    #[test]
    fn ntt_roundtrip_preserves_coefficients() {
        let basis = basis_17_97();
        let mut poly =
            RnsPoly::<8>::from_coeffs(&[1, -2, 3, 4, -5, 6, 7, -8], basis);
        let original = poly.channels().to_vec();

        poly.to_ntt_domain();
        assert!(poly.is_ntt_domain());
        poly.to_coeff_domain();
        assert!(!poly.is_ntt_domain());

        assert_eq!(poly.channels(), original.as_slice());
    }

    #[test]
    fn to_ntt_is_idempotent() {
        let basis = basis_17_97();
        let mut poly = RnsPoly::<8>::from_coeffs(&[1, 2, 3, 4, 5, 6, 7, 8], basis);
        poly.to_ntt_domain();
        let after_first = poly.channels().to_vec();
        poly.to_ntt_domain(); // no-op
        assert_eq!(poly.channels(), after_first.as_slice());
    }

    // ── Mod-drop ──────────────────────────────────────────────────────────────

    #[test]
    fn mod_drop_removes_last_channels() {
        let basis = basis_three();
        let poly = RnsPoly::<8>::from_coeffs(&[1, 2, 3, 4, 5, 6, 7, 8], basis);
        let dropped = poly.mod_drop_last(1).unwrap();
        assert_eq!(dropped.basis().channel_count(), 2);
        assert_eq!(dropped.channels().len(), 2);
    }

    // ── Arithmetic ────────────────────────────────────────────────────────────

    #[test]
    fn add_assign_computes_correct_sum() {
        let basis = basis_17_97();
        let mut a =
            RnsPoly::<8>::from_coeffs(&[1, 2, 3, 4, 5, 6, 7, 8], basis.clone());
        let b = RnsPoly::<8>::from_coeffs(&[8, 7, 6, 5, 4, 3, 2, 1], basis);
        a += &b;
        // All sums are 9; check channel 0 (mod 17) and channel 1 (mod 97)
        assert!(a.channels()[0].iter().all(|&c| c == 9));
        assert!(a.channels()[1].iter().all(|&c| c == 9));
    }

    #[test]
    fn add_assign_wraps_at_modulus() {
        let basis = basis_17_97();
        let mut a =
            RnsPoly::<8>::from_coeffs(&[16, 0, 0, 0, 0, 0, 0, 0], basis.clone());
        let b = RnsPoly::<8>::from_coeffs(&[2, 0, 0, 0, 0, 0, 0, 0], basis);
        a += &b;
        // (16 + 2) mod 17 = 1
        assert_eq!(a.channels()[0][0], 1);
    }

    #[test]
    fn neg_negates_coefficients() {
        let basis = basis_17_97();
        let poly = RnsPoly::<8>::from_coeffs(&[3, 0, 0, 0, 0, 0, 0, 0], basis);
        let neg = -poly;
        // -3 mod 17 = 14
        assert_eq!(neg.channels()[0][0], 14);
        // 0 stays 0
        assert_eq!(neg.channels()[0][1], 0);
    }

    #[test]
    fn mul_assign_schoolbook_small() {
        // (1 + x) * (1 + x) = 1 + 2x + x^2  in Z[x]/(x^8+1)
        let basis = basis_17_97();
        let mut a =
            RnsPoly::<8>::from_coeffs(&[1, 1, 0, 0, 0, 0, 0, 0], basis.clone());
        let b = RnsPoly::<8>::from_coeffs(&[1, 1, 0, 0, 0, 0, 0, 0], basis);
        a *= &b;
        // coefficients: [1, 2, 1, 0, …]
        let coeffs = a.to_coeffs();
        assert_eq!(coeffs[0], 1);
        assert_eq!(coeffs[1], 2);
        assert_eq!(coeffs[2], 1);
        assert!(coeffs[3..].iter().all(|&c| c == 0));
    }

    #[test]
    fn mul_assign_wraps_around_quotient() {
        // x^7 * x = x^8 = -1  in Z[x]/(x^8+1)
        let basis = basis_17_97();
        let mut a =
            RnsPoly::<8>::from_coeffs(&[0, 0, 0, 0, 0, 0, 0, 1], basis.clone());
        let b = RnsPoly::<8>::from_coeffs(&[0, 1, 0, 0, 0, 0, 0, 0], basis);
        a *= &b;
        let coeffs = a.to_coeffs();
        assert_eq!(coeffs[0], -1);
        assert!(coeffs[1..].iter().all(|&c| c == 0));
    }

    // ── PolyRing trait ────────────────────────────────────────────────────────

    #[test]
    fn to_coeffs_roundtrips_from_coeffs() {
        let basis = basis_17_97();
        let input = [1i64, -2, 3, -4, 5, -6, 7, -8];
        let poly = RnsPoly::<8>::from_coeffs(&input, basis);
        assert_eq!(poly.to_coeffs(), input);
    }

    #[test]
    fn to_coeffs_works_from_ntt_domain() {
        let basis = basis_17_97();
        let input = [1i64, -2, 3, -4, 5, -6, 7, -8];
        let mut poly = RnsPoly::<8>::from_coeffs(&input, basis);
        poly.to_ntt_domain();
        // to_coeffs should auto-convert without mutating poly
        assert_eq!(poly.to_coeffs(), input);
        assert!(poly.is_ntt_domain(), "poly should still be in NTT domain");
    }

    // ── PolySampler trait ─────────────────────────────────────────────────────

    #[test]
    fn sample_uniform_stays_in_range() {
        let basis = basis_17_97();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let poly = RnsPoly::<8>::sample_uniform(&basis, &mut rng);
        for (ch, channel) in poly.channels().iter().enumerate() {
            let q = basis.moduli()[ch];
            for &c in channel {
                assert!(c < q);
            }
        }
    }

    #[test]
    fn mul_assign_ntt_domain_matches_coeff_domain() {
        // Multiply in coefficient domain, then compare to multiplying after
        // converting both operands to NTT domain. Results must be identical.
        let basis = basis_17_97();
        let a_coeffs = [1i64, -2, 3, -4, 5, -6, 7, -8];
        let b_coeffs = [2i64, 1, -1, 3, 0, -2, 4, 1];

        // Coefficient-domain path (goes through NTT internally).
        let mut a_coeff = RnsPoly::<8>::from_coeffs(&a_coeffs, basis.clone());
        let b_coeff = RnsPoly::<8>::from_coeffs(&b_coeffs, basis.clone());
        a_coeff *= &b_coeff;
        let expected = a_coeff.to_coeffs();

        // NTT-domain path: convert both first, then pointwise multiply.
        let mut a_ntt = RnsPoly::<8>::from_coeffs(&a_coeffs, basis.clone());
        let mut b_ntt = RnsPoly::<8>::from_coeffs(&b_coeffs, basis.clone());
        a_ntt.to_ntt_domain();
        b_ntt.to_ntt_domain();
        a_ntt *= &b_ntt;
        assert!(a_ntt.is_ntt_domain(), "result should stay in NTT domain");
        let got = a_ntt.to_coeffs();

        assert_eq!(expected, got);
    }

    #[test]
    fn mul_assign_matches_naive() {
        // NTT-based mul_assign must produce the same result as mul_assign_naive.
        let basis = basis_17_97();
        let a_coeffs = [3i64, 1, -2, 0, 4, -1, 2, -3];
        let b_coeffs = [1i64, -1, 0, 2, -2, 3, 1, -1];

        let mut a_ntt = RnsPoly::<8>::from_coeffs(&a_coeffs, basis.clone());
        let b = RnsPoly::<8>::from_coeffs(&b_coeffs, basis.clone());
        a_ntt *= &b;

        let mut a_naive = RnsPoly::<8>::from_coeffs(&a_coeffs, basis.clone());
        let b_naive = RnsPoly::<8>::from_coeffs(&b_coeffs, basis.clone());
        a_naive.mul_assign_naive(&b_naive);

        assert_eq!(a_ntt.to_coeffs(), a_naive.to_coeffs());
    }

    #[test]
    fn sample_tribits_has_correct_hamming_weight() {
        let basis = basis_17_97();
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let hw = 3;
        let poly = RnsPoly::<8>::sample_tribits(hw, &basis, &mut rng);
        // Reconstruct and count non-zero values in channel 0 (mod 17)
        // ternary {-1,0,1} → {16,0,1} mod 17; 0 stays 0
        let non_zero = poly.channels()[0].iter().filter(|&&c| c != 0).count();
        assert_eq!(non_zero, hw);
    }

    // ── Rescale ───────────────────────────────────────────────────────────────

    #[test]
    fn rescale_drops_channel_count() {
        let basis = basis_three();
        let poly = RnsPoly::<8>::from_coeffs(&[1, 2, 3, 4, 5, 6, 7, 8], basis);
        let rescaled = poly.rescale().unwrap();
        assert_eq!(rescaled.basis().channel_count(), 2);
        assert_eq!(rescaled.channels().len(), 2);
        assert!(!rescaled.is_ntt_domain());
    }

    #[test]
    fn rescale_on_single_channel_errors() {
        let basis = Arc::new(RnsBasis::<8>::new(vec![17]).unwrap());
        let poly = RnsPoly::<8>::from_coeffs(&[1, 2, 3, 4, 5, 6, 7, 8], basis);
        assert!(matches!(
            poly.rescale(),
            Err(RnsNttError::InvalidModDrop { .. })
        ));
    }

    #[test]
    fn rescale_is_exact_division_by_last_prime() {
        // Construct a polynomial whose integer coefficients are multiples of q_last=113.
        // After rescaling by 113, the result should equal the original / 113.
        //
        // basis: q0=17, q1=97, q2=113.  We want c = 113 * v for small v.
        // Choose c = 113 * 2 = 226 (fits in signed i64).
        let basis = basis_three(); // moduli [17, 97, 113]
        let v: i64 = 2;
        let c: i64 = 113 * v;

        // All-constant polynomial: every coefficient is c.
        let coeffs: Vec<i64> = vec![c; 8];
        let poly = RnsPoly::<8>::from_coeffs(&coeffs, basis);
        let rescaled = poly.rescale().unwrap();

        // After rescale by q_last=113, every coefficient should equal v=2.
        let out = rescaled.to_coeffs();
        assert!(
            out.iter().all(|&x| x == v),
            "expected all coefficients to be {v}, got {out:?}"
        );
    }

    #[test]
    fn rescale_from_ntt_domain_matches_coeff_domain() {
        // Rescaling from NTT domain should give the same result as from coeff domain.
        let basis = basis_three();
        let coeffs = [3i64, -1, 2, 0, -4, 5, -2, 1];

        let poly_coeff = RnsPoly::<8>::from_coeffs(&coeffs, basis.clone());
        let expected = poly_coeff.rescale().unwrap().to_coeffs();

        let mut poly_ntt = RnsPoly::<8>::from_coeffs(&coeffs, basis);
        poly_ntt.to_ntt_domain();
        let got = poly_ntt.rescale().unwrap().to_coeffs();

        assert_eq!(expected, got);
    }
}
