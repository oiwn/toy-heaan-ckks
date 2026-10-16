use super::basis::{NttTable, RnsBasis, mod_inverse, reverse_bits};
use crate::{
    math::sampling::{ternary_coefficients, uniform_coefficients},
    rings::{
        errors::{RingError, RingResult},
        traits::PolyRing,
    },
};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::{
    borrow::Cow,
    ops::{AddAssign, MulAssign, Neg, SubAssign},
    sync::Arc,
};

/// A polynomial in `Z_{q_0} x … x Z_{q_{L-1}}[X] / (X^N + 1)`.
///
/// Stores one `[u64; DEGREE]` array per RNS channel. Polynomials rest in NTT
/// form (OpenFHE/SEAL style): add, sub, neg and mul are pointwise per
/// channel. Coefficient form is only used inside the operations that need
/// it: construction and sampling, decoding, gadget decomposition, rescale
/// (last channel only) and automorphisms. The private `in_ntt_form` flag
/// records which form the arrays hold; every polynomial this module hands
/// out is in NTT form.
///
/// # Invariants
/// - `channels.len() == basis.channel_count()`
/// - Every `channels[i][j] < basis.moduli()[i]`
/// - `DEGREE` is a power of two
#[derive(Clone, Debug)]
pub struct RnsPoly<const DEGREE: usize> {
    channels: Vec<[u64; DEGREE]>,
    basis: Arc<RnsBasis<DEGREE>>,
    in_ntt_form: bool,
}

// ─── Construction & form conversion (private) ────────────────────────────────

impl<const DEGREE: usize> RnsPoly<DEGREE> {
    /// The basis (modulus chain node) this polynomial lives on.
    pub fn basis(&self) -> &Arc<RnsBasis<DEGREE>> {
        &self.basis
    }

    /// Wraps coefficient-form channels and moves them to NTT form.
    fn from_coeff_channels(
        channels: Vec<[u64; DEGREE]>,
        basis: Arc<RnsBasis<DEGREE>>,
    ) -> Self {
        debug_assert_eq!(channels.len(), basis.channel_count());
        let mut poly = Self {
            channels,
            basis,
            in_ntt_form: false,
        };
        poly.make_ntt_form();
        poly
    }

    /// Moves to NTT form in place (no-op if already there).
    ///
    /// Applies a pre-twist (multiply coefficient `j` by psi^j) before the
    /// forward NTT, so that pointwise multiplication in this form equals
    /// negacyclic convolution in Z[X]/(X^N + 1).
    fn make_ntt_form(&mut self) {
        if self.in_ntt_form {
            return;
        }
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            forward_twisted(channel, self.basis.ntt_table(ch));
        }
        self.in_ntt_form = true;
    }

    /// Moves to coefficient form in place (no-op if already there).
    fn make_coeff_form(&mut self) {
        if !self.in_ntt_form {
            return;
        }
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            inverse_untwisted(channel, self.basis.ntt_table(ch));
        }
        self.in_ntt_form = false;
    }

    /// `self` in NTT form, cloning only if a conversion is needed.
    fn ntt_form(&self) -> Cow<'_, Self> {
        if self.in_ntt_form {
            Cow::Borrowed(self)
        } else {
            let mut poly = self.clone();
            poly.make_ntt_form();
            Cow::Owned(poly)
        }
    }

    /// `self` in coefficient form, cloning only if a conversion is needed.
    fn coeff_form(&self) -> Cow<'_, Self> {
        if self.in_ntt_form {
            let mut poly = self.clone();
            poly.make_coeff_form();
            Cow::Owned(poly)
        } else {
            Cow::Borrowed(self)
        }
    }

    /// Per-coefficient residue vectors, for CRT reconstruction.
    fn for_each_residue_vector(&self, mut visit: impl FnMut(usize, &[u64])) {
        let poly = self.coeff_form();
        let mut residues = vec![0u64; poly.basis.channel_count()];
        for idx in 0..DEGREE {
            for (residue, channel) in residues.iter_mut().zip(&poly.channels) {
                *residue = channel[idx];
            }
            visit(idx, &residues);
        }
    }
}

// ─── Arithmetic ───────────────────────────────────────────────────────────────

impl<const DEGREE: usize> AddAssign<&RnsPoly<DEGREE>> for RnsPoly<DEGREE> {
    /// Pointwise addition modulo each `q_i`; both operands share the basis.
    fn add_assign(&mut self, rhs: &RnsPoly<DEGREE>) {
        debug_assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "add_assign: basis mismatch"
        );
        self.make_ntt_form();
        let rhs = rhs.ntt_form();
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            let q = self.basis.moduli()[ch];
            for (a, &b) in channel.iter_mut().zip(rhs.channels[ch].iter()) {
                *a = add_mod(*a, b, q);
            }
        }
    }
}

impl<const DEGREE: usize> SubAssign<&RnsPoly<DEGREE>> for RnsPoly<DEGREE> {
    /// Pointwise subtraction modulo each `q_i`; both operands share the basis.
    fn sub_assign(&mut self, rhs: &RnsPoly<DEGREE>) {
        debug_assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "sub_assign: basis mismatch"
        );
        self.make_ntt_form();
        let rhs = rhs.ntt_form();
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            let q = self.basis.moduli()[ch];
            for (a, &b) in channel.iter_mut().zip(rhs.channels[ch].iter()) {
                *a = sub_mod(*a, b, q);
            }
        }
    }
}

impl<const DEGREE: usize> MulAssign<&RnsPoly<DEGREE>> for RnsPoly<DEGREE> {
    /// Negacyclic product in `Z[X]/(X^N + 1)`: pointwise per channel in NTT
    /// form, O(N·L).
    fn mul_assign(&mut self, rhs: &RnsPoly<DEGREE>) {
        debug_assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "mul_assign: basis mismatch"
        );
        self.make_ntt_form();
        let rhs = rhs.ntt_form();
        for (ch, channel) in self.channels.iter_mut().enumerate() {
            let q = self.basis.moduli()[ch];
            for (a, &b) in channel.iter_mut().zip(rhs.channels[ch].iter()) {
                *a = mul_mod(*a, b, q);
            }
        }
    }
}

impl<const DEGREE: usize> Neg for RnsPoly<DEGREE> {
    type Output = Self;

    /// Negation modulo each `q_i`; the same in both forms.
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
        // Zero is zero in both forms.
        Self {
            channels: vec![[0u64; DEGREE]; context.channel_count()],
            basis: context.clone(),
            in_ntt_form: true,
        }
    }

    fn context(&self) -> &Self::Context {
        &self.basis
    }

    /// Each coefficient is reduced into `[0, q_i)` per channel.
    ///
    /// # Panics
    /// Panics if `coeffs` has fewer than `DEGREE` entries.
    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        assert!(
            coeffs.len() >= DEGREE,
            "from_coeffs: need at least {DEGREE} coefficients, got {}",
            coeffs.len()
        );
        let mut channels = vec![[0u64; DEGREE]; context.channel_count()];
        for (ch, channel) in channels.iter_mut().enumerate() {
            let q = context.moduli()[ch] as i128;
            for (i, coeff) in coeffs.iter().take(DEGREE).enumerate() {
                channel[i] = (*coeff as i128).rem_euclid(q) as u64;
            }
        }
        Self::from_coeff_channels(channels, context.clone())
    }

    /// CRT-reconstructs each coefficient and centers it in `(-Q/2, Q/2]`.
    fn to_coeffs(&self) -> [i64; DEGREE] {
        let mut result = [0i64; DEGREE];
        self.for_each_residue_vector(|idx, residues| {
            result[idx] = self.basis.reconstruct_centered_coeff(residues);
        });
        result
    }

    /// Reduces each value into every prime straight from its f64 mantissa
    /// and exponent (the approach SEAL takes for > 64-bit coefficients).
    /// The limit `2^(total_bits - 1)` is a conservative `Q/2` bound.
    fn from_rounded_f64(
        coeffs: &[f64; DEGREE],
        context: &Self::Context,
    ) -> RingResult<Self> {
        let max_bits = context.total_bits().saturating_sub(1);
        let limit = 2f64.powi(max_bits as i32);
        if let Some(&value) = coeffs
            .iter()
            .find(|value| !value.is_finite() || value.abs() >= limit)
        {
            return Err(RingError::CoefficientTooLarge { value, max_bits });
        }

        let channels = context
            .moduli()
            .iter()
            .map(|&prime| {
                let mut channel = [0u64; DEGREE];
                for (residue, &value) in channel.iter_mut().zip(coeffs.iter()) {
                    *residue = rounded_to_residue(value, prime);
                }
                channel
            })
            .collect();
        Ok(Self::from_coeff_channels(channels, context.clone()))
    }

    /// CRT-reconstructs straight into f64, so it works when `Q` exceeds the
    /// `i64` range.
    fn to_centered_f64(&self) -> [f64; DEGREE] {
        let mut result = [0f64; DEGREE];
        self.for_each_residue_vector(|idx, residues| {
            result[idx] = self.basis.reconstruct_centered_coeff_f64(residues);
        });
        result
    }

    /// Coefficients uniform in `[0, q_i)` per channel.
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self {
        let channels = context
            .moduli()
            .iter()
            .map(|&q| uniform_coefficients::<DEGREE, _>(q, rng))
            .collect();
        Self::from_coeff_channels(channels, context.clone())
    }

    /// Noise from N(0, std_dev), rounded and reduced per channel.
    ///
    /// # Panics
    /// Panics if `std_dev` is not finite and positive.
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
        Self::from_coeffs(&noise, context)
    }

    fn sample_ternary<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let ternary = ternary_coefficients::<DEGREE, _>(hamming_weight, rng);
        Self::from_coeffs(&ternary, context)
    }

    /// Permutes coefficients: `a_i` moves to `j = i * exponent mod 2N`,
    /// negated when `j >= N` (since `X^N = -1`).
    fn automorphism(&self, exponent: u64) -> Self {
        let modulus_2n = (2 * DEGREE) as u64;
        let exponent = exponent % modulus_2n;
        if exponent == 0 {
            // X^0 = 1: identity automorphism
            return self.clone();
        }

        let source = self.coeff_form();
        let mut channels = vec![[0u64; DEGREE]; self.basis.channel_count()];
        for (ch, new_channel) in channels.iter_mut().enumerate() {
            let q = self.basis.moduli()[ch];
            for (i, &coeff) in source.channels[ch].iter().enumerate() {
                if coeff == 0 {
                    continue;
                }
                let j_full = (i as u64 * exponent) % modulus_2n;
                let j = (j_full % DEGREE as u64) as usize;
                new_channel[j] = if j_full >= DEGREE as u64 {
                    q - coeff
                } else {
                    coeff
                };
            }
        }
        Self::from_coeff_channels(channels, self.basis.clone())
    }

    /// Number of RNS channels (primes) still carried.
    fn level(&self) -> usize {
        self.basis.channel_count()
    }

    /// `sum floor(log2(q_i))`, see [`RnsBasis::total_bits`].
    fn modulus_bits(&self) -> u32 {
        self.basis.total_bits()
    }

    /// Bit length of the last prime, `floor(log2(q_last)) + 1`.
    fn rescale_bits(&self) -> u32 {
        let q_last = *self.basis.moduli().last().expect("non-empty basis");
        u64::BITS - q_last.leading_zeros()
    }

    /// Divides every coefficient by the last prime `q_L` and drops it.
    ///
    /// Exact (no rounding): `c - (c mod q_L)` is divisible by `q_L`. Per
    /// remaining channel `i`:
    ///
    /// ```text
    /// new_i = (c_i - (c_L mod q_i)) * q_L^-1  mod q_i
    /// ```
    ///
    /// Only the last channel leaves NTT form: it goes back to coefficients,
    /// is reduced modulo each `q_i` and re-enters NTT form under `q_i`. The
    /// formula is linear, so the rest stays pointwise. The result lives on
    /// the next basis of the shared chain.
    fn rescale(&self) -> RingResult<Self> {
        let num_channels = self.basis.channel_count();
        let new_basis =
            self.basis.next().cloned().ok_or(RingError::InvalidLevel {
                level: num_channels - 1,
                max_level: num_channels,
            })?;

        let source = self.ntt_form();
        let last = num_channels - 1;
        let q_last = self.basis.moduli()[last];
        let mut last_coeffs = source.channels[last];
        inverse_untwisted(&mut last_coeffs, self.basis.ntt_table(last));

        let mut channels = vec![[0u64; DEGREE]; last];
        for (i, new_channel) in channels.iter_mut().enumerate() {
            let q_i = self.basis.moduli()[i];
            let q_last_inv = mod_inverse(q_last % q_i, q_i);
            let mut last_mod_qi = last_coeffs.map(|c| c % q_i);
            forward_twisted(&mut last_mod_qi, self.basis.ntt_table(i));
            for (j, value) in new_channel.iter_mut().enumerate() {
                let diff = sub_mod(source.channels[i][j], last_mod_qi[j], q_i);
                *value = mul_mod(diff, q_last_inv, q_i);
            }
        }

        Ok(Self {
            channels,
            basis: new_basis,
            in_ntt_form: true,
        })
    }

    /// Keeps the first `level` channels, without dividing.
    ///
    /// Channels are independent, so this is exact in either form. The result
    /// uses the shared chain basis for `level`.
    fn mod_switch_to(&self, level: usize) -> RingResult<Self> {
        let basis = self.basis.at_level(level).ok_or(RingError::InvalidLevel {
            level,
            max_level: self.basis.channel_count(),
        })?;
        Ok(Self {
            channels: self.channels[..level].to_vec(),
            basis,
            in_ntt_form: self.in_ntt_form,
        })
    }

    /// `sum_i ceil(bits(q_i) / digit_bits)`: one term per (channel, digit).
    fn gadget_len(context: &Self::Context, digit_bits: u32) -> usize {
        context
            .moduli()
            .iter()
            .map(|&q| digit_count(q, digit_bits))
            .sum()
    }

    /// RNS gadget decomposition, terms ordered by (channel `i`, digit `j`).
    ///
    /// Term `(i, j)` is the `j`-th base-2^`digit_bits` digit of the channel-`i`
    /// residues, `digit_j(alpha_i)`, copied into every channel. Since
    /// `alpha_i = sum_j digit_j(alpha_i) * 2^(j·digit_bits)` and `x` is the
    /// CRT combination of its residues `alpha_i`, pairing with
    /// [`PolyRing::gadget_term`] rebuilds `x * t`.
    fn decompose(&self, digit_bits: u32) -> Vec<Self> {
        let source = self.coeff_form();
        let moduli = self.basis.moduli();
        let mask = (1u64 << digit_bits) - 1;
        let mut terms =
            Vec::with_capacity(Self::gadget_len(&self.basis, digit_bits));

        for (i, &q_i) in moduli.iter().enumerate() {
            for j in 0..digit_count(q_i, digit_bits) {
                let shift = j as u32 * digit_bits;
                let digits = source.channels[i].map(|c| (c >> shift) & mask);
                let channels =
                    moduli.iter().map(|&q| digits.map(|d| d % q)).collect();
                terms.push(Self::from_coeff_channels(channels, self.basis.clone()));
            }
        }
        terms
    }

    /// Term `(i, j)`: channel `i` times `2^(j·digit_bits)`, every other
    /// channel zero (`e_i · 2^(j·digit_bits) · self`, where `e_i` is the CRT
    /// basis element). Linear per channel, so it works in NTT form directly.
    ///
    /// # Panics
    /// Panics if `index >= gadget_len`.
    fn gadget_term(&self, index: usize, digit_bits: u32) -> Self {
        let mut remaining = index;
        for (i, &q_i) in self.basis.moduli().iter().enumerate() {
            let digits = digit_count(q_i, digit_bits);
            if remaining < digits {
                let factor = pow2_mod(remaining as u64 * digit_bits as u64, q_i);
                let mut channels = vec![[0u64; DEGREE]; self.basis.channel_count()];
                channels[i] = self.channels[i].map(|c| mul_mod(c, factor, q_i));
                return Self {
                    channels,
                    basis: self.basis.clone(),
                    in_ntt_form: self.in_ntt_form,
                };
            }
            remaining -= digits;
        }
        panic!("gadget_term: index {index} out of range");
    }
}

/// Number of base-2^`digit_bits` digits needed for residues in `[0, modulus)`.
fn digit_count(modulus: u64, digit_bits: u32) -> usize {
    let modulus_bits = u64::BITS - modulus.leading_zeros();
    modulus_bits.div_ceil(digit_bits) as usize
}

/// Reduces an integer-valued `f64` into `[0, prime)`.
///
/// Values below 2^63 go through `i128`. Larger ones are exact
/// `mantissa * 2^exponent` products (53-bit mantissa), reduced as
/// `(mantissa mod p) * (2^exponent mod p) mod p`.
fn rounded_to_residue(value: f64, prime: u64) -> u64 {
    debug_assert!(value.is_finite() && value == value.trunc());

    let magnitude = value.abs();
    let residue = if magnitude < 2f64.powi(63) {
        ((magnitude as i128) % prime as i128) as u64
    } else {
        let bits = magnitude.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i64 - 1075;
        let mantissa = (bits & ((1u64 << 52) - 1)) | (1u64 << 52);
        debug_assert!(exponent > 0, "values >= 2^63 have a positive exponent");
        mul_mod(mantissa % prime, pow2_mod(exponent as u64, prime), prime)
    };

    if value < 0.0 && residue != 0 {
        prime - residue
    } else {
        residue
    }
}

// ─── NTT kernels (private) ────────────────────────────────────────────────────

/// Pre-twist by psi^j, then forward NTT: coefficients -> NTT form.
fn forward_twisted<const DEGREE: usize>(
    channel: &mut [u64; DEGREE],
    table: &NttTable<DEGREE>,
) {
    for (j, v) in channel.iter_mut().enumerate() {
        *v = mul_mod(*v, table.twist_factors[j], table.modulus);
    }
    forward_ntt(channel, table);
}

/// Inverse NTT, then post-untwist by psi^-j: NTT form -> coefficients.
fn inverse_untwisted<const DEGREE: usize>(
    channel: &mut [u64; DEGREE],
    table: &NttTable<DEGREE>,
) {
    inverse_ntt(channel, table);
    for (j, v) in channel.iter_mut().enumerate() {
        *v = mul_mod(*v, table.untwist_factors[j], table.modulus);
    }
}

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

fn pow2_mod(mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1 % modulus;
    let mut base = 2 % modulus;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mul_mod(result, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    result
}

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

    fn poly(coeffs: &[i64], basis: &Arc<RnsBasis<8>>) -> RnsPoly<8> {
        RnsPoly::from_coeffs(coeffs, basis)
    }

    /// Coefficient-form residues, channel by channel.
    fn residues(poly: &RnsPoly<8>) -> Vec<[u64; 8]> {
        poly.coeff_form().channels.clone()
    }

    /// Schoolbook O(N^2) negacyclic product, the reference for `mul_assign`.
    fn mul_naive(a: &RnsPoly<8>, b: &RnsPoly<8>) -> Vec<[u64; 8]> {
        let (a, b) = (residues(a), residues(b));
        let moduli = [17u64, 97];
        let mut out = vec![[0u64; 8]; 2];
        for ch in 0..2 {
            let q = moduli[ch];
            for i in 0..8 {
                for j in 0..8 {
                    let prod = mul_mod(a[ch][i], b[ch][j], q);
                    if i + j < 8 {
                        out[ch][i + j] = add_mod(out[ch][i + j], prod, q);
                    } else {
                        // X^N = -1 in Z[X]/(X^N + 1)
                        out[ch][i + j - 8] = sub_mod(out[ch][i + j - 8], prod, q);
                    }
                }
            }
        }
        out
    }

    // ── Construction & form ───────────────────────────────────────────────────

    #[test]
    fn zero_poly_is_all_zeros() {
        let zero = RnsPoly::<8>::zero(&basis_17_97());
        assert_eq!(zero.to_coeffs(), [0; 8]);
    }

    #[test]
    fn from_coeffs_reduces_correctly() {
        let p = poly(&[-1, 2, -3, 4, 0, 0, 0, 0], &basis_17_97());
        let r = residues(&p);
        // channel 0: mod 17 -> [-1->16, 2->2, -3->14, 4->4, …]
        assert_eq!(&r[0][..3], &[16, 2, 14]);
        // channel 1: mod 97 -> [-1->96, …]
        assert_eq!(r[1][0], 96);
    }

    #[test]
    fn polys_rest_in_ntt_form() {
        let basis = basis_three();
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let a = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis);
        assert!(a.in_ntt_form);
        assert!(RnsPoly::<8>::sample_uniform(&basis, &mut rng).in_ntt_form);
        assert!(a.rescale().unwrap().in_ntt_form);
        assert!(a.automorphism(3).in_ntt_form);
        assert!(a.decompose(4).iter().all(|t| t.in_ntt_form));
    }

    #[test]
    fn ntt_roundtrip_preserves_coefficients() {
        let mut p = poly(&[1, -2, 3, 4, -5, 6, 7, -8], &basis_17_97());
        let ntt = p.channels.clone();
        p.make_coeff_form();
        let coeffs = p.channels.clone();
        p.make_coeff_form(); // no-op
        assert_eq!(p.channels, coeffs);
        p.make_ntt_form();
        assert_eq!(p.channels, ntt);
        assert_eq!(p.to_coeffs(), [1, -2, 3, 4, -5, 6, 7, -8]);
    }

    #[test]
    fn mixed_forms_add_like_ntt_forms() {
        let basis = basis_17_97();
        let a = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis);
        let b = poly(&[8, 7, 6, 5, 4, 3, 2, 1], &basis);
        let mut b_coeff = b.clone();
        b_coeff.make_coeff_form();

        let mut sum = a.clone();
        sum += &b;
        let mut mixed = a;
        mixed += &b_coeff;
        assert!(mixed.in_ntt_form);
        assert_eq!(mixed.channels, sum.channels);
    }

    // ── Modulus chain ─────────────────────────────────────────────────────────

    #[test]
    fn chain_reports_level_and_bits() {
        let basis = basis_three(); // 17, 97, 113: floor bits 4 + 6 + 6
        let p = poly(&[1; 8], &basis);
        assert_eq!(p.level(), 3);
        assert_eq!(p.modulus_bits(), 16);
        assert_eq!(p.rescale_bits(), 7); // 113 < 2^7
        let down = p.rescale().unwrap();
        assert_eq!(down.level(), 2);
        assert_eq!(down.rescale_bits(), 7); // 97 < 2^7
    }

    #[test]
    fn mod_switch_keeps_first_channels_on_the_chain_basis() {
        let basis = basis_three();
        let p = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis);
        let switched = p.mod_switch_to(2).unwrap();
        assert_eq!(residues(&switched), residues(&p)[..2]);
        assert!(Arc::ptr_eq(switched.basis(), basis.next().unwrap()));
    }

    #[test]
    fn mod_switch_rejects_out_of_range_levels() {
        let p = poly(&[1; 8], &basis_three());
        for level in [0, 4] {
            assert!(matches!(
                p.mod_switch_to(level),
                Err(RingError::InvalidLevel { .. })
            ));
        }
    }

    // ── Arithmetic ────────────────────────────────────────────────────────────

    #[test]
    fn add_assign_computes_correct_sum() {
        let basis = basis_17_97();
        let mut a = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis);
        a += &poly(&[8, 7, 6, 5, 4, 3, 2, 1], &basis);
        assert_eq!(a.to_coeffs(), [9; 8]);
    }

    #[test]
    fn add_assign_wraps_at_modulus() {
        let basis = basis_17_97();
        let mut a = poly(&[16, 0, 0, 0, 0, 0, 0, 0], &basis);
        a += &poly(&[2, 0, 0, 0, 0, 0, 0, 0], &basis);
        // (16 + 2) mod 17 = 1
        assert_eq!(residues(&a)[0][0], 1);
    }

    #[test]
    fn sub_assign_computes_correct_difference() {
        let basis = basis_17_97();
        let mut a = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis);
        a -= &poly(&[8, 7, 6, 5, 4, 3, 2, 1], &basis);
        assert_eq!(a.to_coeffs(), [-7, -5, -3, -1, 1, 3, 5, 7]);
    }

    #[test]
    fn neg_negates_coefficients() {
        let neg = -poly(&[3, 0, 0, 0, 0, 0, 0, 0], &basis_17_97());
        // -3 mod 17 = 14; 0 stays 0
        assert_eq!(&residues(&neg)[0][..2], &[14, 0]);
    }

    #[test]
    fn mul_assign_schoolbook_small() {
        // (1 + x) * (1 + x) = 1 + 2x + x^2  in Z[x]/(x^8+1)
        let basis = basis_17_97();
        let mut a = poly(&[1, 1, 0, 0, 0, 0, 0, 0], &basis);
        a *= &poly(&[1, 1, 0, 0, 0, 0, 0, 0], &basis);
        assert_eq!(a.to_coeffs(), [1, 2, 1, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn mul_assign_wraps_around_quotient() {
        // x^7 * x = x^8 = -1  in Z[x]/(x^8+1)
        let basis = basis_17_97();
        let mut a = poly(&[0, 0, 0, 0, 0, 0, 0, 1], &basis);
        a *= &poly(&[0, 1, 0, 0, 0, 0, 0, 0], &basis);
        assert_eq!(a.to_coeffs(), [-1, 0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn mul_assign_matches_naive() {
        let basis = basis_17_97();
        let mut a = poly(&[3, 1, -2, 0, 4, -1, 2, -3], &basis);
        let b = poly(&[1, -1, 0, 2, -2, 3, 1, -1], &basis);
        let expected = mul_naive(&a, &b);
        a *= &b;
        assert_eq!(residues(&a), expected);
    }

    // ── Coefficients & encoding boundary ──────────────────────────────────────

    #[test]
    fn to_coeffs_roundtrips_from_coeffs() {
        let input = [1i64, -2, 3, -4, 5, -6, 7, -8];
        assert_eq!(poly(&input, &basis_17_97()).to_coeffs(), input);
    }

    /// Four 61-bit NTT primes for DEGREE = 16: Q ~ 2^244.
    fn wide_basis() -> Arc<RnsBasis<16>> {
        let primes = crate::math::generate_primes(61, 4, 16);
        Arc::new(RnsBasis::new(primes).unwrap())
    }

    #[test]
    fn rounded_f64_roundtrips_beyond_i64() {
        let basis = wide_basis();
        let mut coeffs = [0f64; 16];
        coeffs[0] = 2f64.powi(63);
        coeffs[1] = -2f64.powi(150);
        coeffs[2] = -987_654_321.0;
        coeffs[3] = (2f64.powi(53) - 1.0) * 2f64.powi(70);
        let p = RnsPoly::<16>::from_rounded_f64(&coeffs, &basis).unwrap();
        assert_eq!(p.to_centered_f64(), coeffs);
    }

    #[test]
    fn rounded_f64_rejects_half_modulus_and_non_finite() {
        let basis = wide_basis();
        let max_bits = basis.total_bits() - 1;
        for bad in [2f64.powi(max_bits as i32), f64::NAN, f64::INFINITY] {
            let mut coeffs = [0f64; 16];
            coeffs[5] = bad;
            assert!(matches!(
                RnsPoly::<16>::from_rounded_f64(&coeffs, &basis),
                Err(RingError::CoefficientTooLarge { max_bits: bits, .. })
                    if bits == max_bits
            ));
        }
    }

    /// The f64 mantissa/exponent path must match exact BigInt reduction.
    #[test]
    fn residues_match_bigint_reduction() {
        use num_bigint::BigInt;

        // Exact integer value of an integer-valued f64: mantissa * 2^exponent.
        fn exact(value: f64) -> BigInt {
            if value == 0.0 {
                return BigInt::from(0);
            }
            let bits = value.abs().to_bits();
            let exponent = ((bits >> 52) & 0x7ff) as i64 - 1075;
            let mantissa = BigInt::from((bits & ((1u64 << 52) - 1)) | (1u64 << 52));
            let magnitude = if exponent >= 0 {
                mantissa << exponent as usize
            } else {
                mantissa >> (-exponent) as usize
            };
            if value < 0.0 { -magnitude } else { magnitude }
        }

        let primes = crate::math::generate_primes(61, 3, 16);
        let samples = [
            0.0f64,
            1.0,
            -1.0,
            2f64.powi(62),
            2f64.powi(63),
            -2f64.powi(63),
            2f64.powi(100),
            -2f64.powi(150),
            (2f64.powi(53) - 1.0) * 2f64.powi(70),
            -(2f64.powi(52) + 3.0) * 2f64.powi(90),
            123_456_789_012_345.0,
            -987_654_321.0,
        ];

        for &value in &samples {
            let exact = exact(value);
            for &prime in &primes {
                let p = BigInt::from(prime);
                let want = ((&exact % &p) + &p) % &p;
                let got = rounded_to_residue(value, prime);
                assert_eq!(BigInt::from(got), want, "value {value:e} mod {prime}");
            }
        }
    }

    // ── Sampling ──────────────────────────────────────────────────────────────

    #[test]
    fn sample_uniform_stays_in_range() {
        let basis = basis_17_97();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let p = RnsPoly::<8>::sample_uniform(&basis, &mut rng);
        for (ch, channel) in residues(&p).iter().enumerate() {
            assert!(channel.iter().all(|&c| c < basis.moduli()[ch]));
        }
    }

    #[test]
    fn sample_ternary_has_correct_hamming_weight() {
        let basis = basis_17_97();
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let p = RnsPoly::<8>::sample_ternary(3, &basis, &mut rng);
        let coeffs = p.to_coeffs();
        assert_eq!(coeffs.iter().filter(|&&c| c != 0).count(), 3);
        assert!(coeffs.iter().all(|c| c.abs() <= 1));
    }

    // ── Automorphism ──────────────────────────────────────────────────────────

    #[test]
    fn automorphism_identity_preserves_coefficients() {
        let coeffs = [1i64, 2, 3, 4, 5, 6, 7, 8];
        let p = poly(&coeffs, &basis_17_97());
        // X -> X^1, and X -> X^(2N) = X^16 (mod 2N = 16)
        assert_eq!(p.automorphism(1).to_coeffs(), coeffs);
        assert_eq!(p.automorphism(16).to_coeffs(), coeffs);
    }

    #[test]
    fn automorphism_applies_sign_change_correctly() {
        // X -> X^9 with N = 8: x^1 -> x^9 = -x, so 1 + x -> 1 - x
        let p = poly(&[1, 1, 0, 0, 0, 0, 0, 0], &basis_17_97());
        assert_eq!(p.automorphism(9).to_coeffs(), [1, -1, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn rotate_slots_matches_automorphism_by_powers_of_five() {
        let p = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis_17_97());
        // 5^2 = 25 = 9 (mod 16); negative k adds the conjugation X -> X^15
        assert_eq!(p.rotate_slots(2).to_coeffs(), p.automorphism(9).to_coeffs());
        let conj = p.automorphism(5).automorphism(15);
        assert_eq!(p.rotate_slots(-1).to_coeffs(), conj.to_coeffs());
    }

    // ── Gadget ────────────────────────────────────────────────────────────────

    #[test]
    fn gadget_terms_rebuild_the_product() {
        let basis = basis_three(); // 17, 97, 113 -> 2 + 2 + 2 digits of 4 bits
        assert_eq!(RnsPoly::<8>::gadget_len(&basis, 4), 6);

        let x = poly(&[5, -3, 40, 7, -60, 2, 0, 11], &basis);
        let t = poly(&[1, 0, -2, 3, 0, 1, -1, 4], &basis);
        let terms = x.decompose(4);
        assert_eq!(terms.len(), 6);
        assert!(
            terms
                .iter()
                .all(|d| d.to_coeffs().iter().all(|&c| (0..16).contains(&c)))
        );

        let mut sum = RnsPoly::<8>::zero(&basis);
        for (k, digit) in terms.iter().enumerate() {
            let mut term = digit.clone();
            term *= &t.gadget_term(k, 4);
            sum += &term;
        }
        let mut product = x;
        product *= &t;
        assert_eq!(sum.to_coeffs(), product.to_coeffs());
    }

    // ── Rescale ───────────────────────────────────────────────────────────────

    #[test]
    fn rescale_drops_channel_count() {
        let p = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis_three());
        let rescaled = p.rescale().unwrap();
        assert_eq!(rescaled.basis().channel_count(), 2);
        assert_eq!(rescaled.channels.len(), 2);
    }

    #[test]
    fn rescale_on_single_channel_errors() {
        let basis = Arc::new(RnsBasis::<8>::new(vec![17]).unwrap());
        let p = poly(&[1, 2, 3, 4, 5, 6, 7, 8], &basis);
        assert!(matches!(p.rescale(), Err(RingError::InvalidLevel { .. })));
    }

    #[test]
    fn rescales_of_different_polys_share_one_basis() {
        let basis = basis_three();
        let a = poly(&[1; 8], &basis).rescale().unwrap();
        let b = poly(&[2; 8], &basis).rescale().unwrap();
        assert!(Arc::ptr_eq(a.basis(), b.basis()));
    }

    #[test]
    fn rescale_is_exact_division_by_last_prime() {
        // basis 17, 97, 113: c = 113 * 2 rescales to exactly 2.
        let p = poly(&[113 * 2; 8], &basis_three());
        assert_eq!(p.rescale().unwrap().to_coeffs(), [2; 8]);
    }

    #[test]
    fn rescale_from_coefficient_form_matches_ntt_form() {
        let p = poly(&[3, -1, 2, 0, -4, 5, -2, 1], &basis_three());
        let mut p_coeff = p.clone();
        p_coeff.make_coeff_form();
        assert_eq!(
            p_coeff.rescale().unwrap().channels,
            p.rescale().unwrap().channels
        );
    }
}
