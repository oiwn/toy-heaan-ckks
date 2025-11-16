use super::poly::{NttTable, RnsBasis, crt_reconstruct};
use crate::math::sampling::{
    gaussian_coefficients, ternary_coefficients, uniform_coefficients,
};
use crate::{PolyModSwitch, PolyRing, PolySampler};
use rand::Rng;
use std::{
    fmt,
    ops::{Add, AddAssign, MulAssign, Neg},
    sync::Arc,
};

#[derive(Clone, Debug)]
pub struct RnsNttPoly<const DEGREE: usize> {
    pub(crate) coefficients: Vec<[u64; DEGREE]>,
    pub(crate) in_ntt_domain: bool,
    pub(crate) basis: Arc<RnsBasis>,
}

impl<const DEGREE: usize> RnsNttPoly<DEGREE> {
    pub fn zero(basis: Arc<RnsBasis>) -> Self {
        let mut channels = Vec::with_capacity(basis.channel_count());
        for _ in 0..basis.channel_count() {
            channels.push([0u64; DEGREE]);
        }
        Self {
            coefficients: channels,
            in_ntt_domain: true,
            basis,
        }
    }

    pub fn from_i64_slice(coeffs: &[i64], basis: Arc<RnsBasis>) -> Self {
        assert!(
            coeffs.len() >= DEGREE,
            "expected at least {DEGREE} coefficients, got {}",
            coeffs.len()
        );
        let mut fixed = [0i64; DEGREE];
        fixed.copy_from_slice(&coeffs[..DEGREE]);
        let channels = reduce_signed_coeffs(&fixed, &basis);
        Self {
            coefficients: channels,
            in_ntt_domain: false,
            basis,
        }
    }

    /// Builds a polynomial directly from coefficient-domain residues laid out
    /// per prime channel.
    pub fn from_channels(
        channels: Vec<[u64; DEGREE]>,
        basis: Arc<RnsBasis>,
    ) -> Self {
        assert_eq!(
            channels.len(),
            basis.channel_count(),
            "channel count mismatch for RNS polynomial construction"
        );
        Self {
            coefficients: channels,
            in_ntt_domain: false,
            basis,
        }
    }

    pub fn from_u64_slice(coeffs: &[u64; DEGREE], basis: Arc<RnsBasis>) -> Self {
        let channels = reduce_unsigned_coeffs(coeffs, &basis);
        Self {
            coefficients: channels,
            in_ntt_domain: false,
            basis,
        }
    }

    pub fn channels(&self) -> usize {
        self.coefficients.len()
    }

    pub fn is_ntt_domain(&self) -> bool {
        self.in_ntt_domain
    }

    pub fn to_ntt_domain(&mut self) {
        if self.in_ntt_domain {
            return;
        }
        for (idx, channel) in self.coefficients.iter_mut().enumerate() {
            let table = self.basis.ntt_table(idx);
            forward_ntt(channel, table);
        }
        self.in_ntt_domain = true;
    }

    pub fn to_coeff_domain(&mut self) {
        if !self.in_ntt_domain {
            return;
        }
        for (idx, channel) in self.coefficients.iter_mut().enumerate() {
            let table = self.basis.ntt_table(idx);
            inverse_ntt(channel, table);
        }
        self.in_ntt_domain = false;
    }

    fn coeffs_as_i64(&self) -> [i64; DEGREE] {
        debug_assert!(
            !self.in_ntt_domain,
            "coefficients must be in coefficient domain"
        );
        let primes = self.basis.primes();
        let product: u128 = primes.iter().fold(1u128, |acc, &p| acc * p as u128);
        let half = product / 2;
        let mut result = [0i64; DEGREE];

        for coeff_idx in 0..DEGREE {
            let residues: Vec<u64> = self
                .coefficients
                .iter()
                .map(|channel| channel[coeff_idx])
                .collect();
            let combined = crt_reconstruct(&residues, primes) as u128;
            if combined >= half {
                result[coeff_idx] = (combined as i128 - product as i128) as i64;
            } else {
                result[coeff_idx] = combined as i64;
            }
        }

        result
    }

    pub fn coefficient_to_u64(&self, position: usize) -> u64 {
        assert!(position < DEGREE, "coefficient index out of bounds");
        let mut coeff_domain = self.clone();
        coeff_domain.to_coeff_domain();
        let residues: Vec<u64> = coeff_domain
            .coefficients
            .iter()
            .map(|channel| channel[position])
            .collect();
        crt_reconstruct(&residues, self.basis.primes())
    }

    pub fn to_u64_coefficients(&self) -> Vec<u64> {
        let mut coeff_domain = self.clone();
        coeff_domain.to_coeff_domain();
        let primes = coeff_domain.basis.primes().clone();
        (0..DEGREE)
            .map(|idx| {
                let residues: Vec<u64> = coeff_domain
                    .coefficients
                    .iter()
                    .map(|channel| channel[idx])
                    .collect();
                crt_reconstruct(&residues, &primes)
            })
            .collect()
    }

    pub fn to_i64_coefficients(&self) -> Vec<i64> {
        let mut coeff_domain = self.clone();
        coeff_domain.to_coeff_domain();
        coeff_domain.coeffs_as_i64().to_vec()
    }

    pub fn residue_at(&self, channel_idx: usize, coeff_idx: usize) -> u64 {
        assert!(channel_idx < self.channels(), "channel index out of bounds");
        assert!(coeff_idx < DEGREE, "coefficient index out of bounds");
        assert!(
            !self.in_ntt_domain,
            "residue access requires coefficient domain"
        );
        self.coefficients[channel_idx][coeff_idx]
    }
}

impl<const DEGREE: usize> PolyRing<DEGREE> for RnsNttPoly<DEGREE> {
    type Context = Arc<RnsBasis>;

    fn zero(context: &Self::Context) -> Self {
        Self::zero(context.clone())
    }

    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        Self::from_i64_slice(coeffs, context.clone())
    }

    fn to_coeffs(&self) -> [i64; DEGREE] {
        let mut copy = self.clone();
        copy.to_coeff_domain();
        copy.coeffs_as_i64()
    }

    fn context(&self) -> &Self::Context {
        &self.basis
    }
}

impl<const DEGREE: usize> AddAssign<&RnsNttPoly<DEGREE>> for RnsNttPoly<DEGREE> {
    fn add_assign(&mut self, rhs: &RnsNttPoly<DEGREE>) {
        assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "basis mismatch for NTT addition"
        );
        assert!(
            self.in_ntt_domain && rhs.in_ntt_domain,
            "addition requires NTT domain"
        );
        assert_eq!(self.channels(), rhs.channels(), "channel count mismatch");

        for (channel_idx, channel) in self.coefficients.iter_mut().enumerate() {
            let prime = self.basis.primes()[channel_idx];
            for (lhs_coeff, &rhs_coeff) in
                channel.iter_mut().zip(rhs.coefficients[channel_idx].iter())
            {
                let sum = *lhs_coeff + rhs_coeff;
                *lhs_coeff = if sum >= prime { sum - prime } else { sum };
            }
        }
    }
}

impl<const DEGREE: usize> Add<&RnsNttPoly<DEGREE>> for &RnsNttPoly<DEGREE> {
    type Output = RnsNttPoly<DEGREE>;

    fn add(self, rhs: &RnsNttPoly<DEGREE>) -> Self::Output {
        let mut result = self.clone();
        result += rhs;
        result
    }
}

impl<const DEGREE: usize> MulAssign<&RnsNttPoly<DEGREE>> for RnsNttPoly<DEGREE> {
    fn mul_assign(&mut self, rhs: &RnsNttPoly<DEGREE>) {
        assert!(
            Arc::ptr_eq(&self.basis, &rhs.basis),
            "basis mismatch for NTT multiplication"
        );
        assert_eq!(self.channels(), rhs.channels(), "channel count mismatch");
        assert!(
            self.in_ntt_domain && rhs.in_ntt_domain,
            "multiplication requires NTT domain"
        );

        for (channel_idx, channel) in self.coefficients.iter_mut().enumerate() {
            let prime = self.basis.primes()[channel_idx];
            let rhs_channel = &rhs.coefficients[channel_idx];
            for (lhs_coeff, &rhs_coeff) in channel.iter_mut().zip(rhs_channel) {
                *lhs_coeff = mul_mod(*lhs_coeff, rhs_coeff, prime);
            }
        }
        // Result stays in the transform domain.
        self.in_ntt_domain = true;
    }
}

impl<const DEGREE: usize> Neg for RnsNttPoly<DEGREE> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.negate_assign();
        self
    }
}

impl<const DEGREE: usize> Neg for &RnsNttPoly<DEGREE> {
    type Output = RnsNttPoly<DEGREE>;

    fn neg(self) -> Self::Output {
        let mut clone = self.clone();
        clone.negate_assign();
        clone
    }
}

impl<const DEGREE: usize> RnsNttPoly<DEGREE> {
    pub fn negate_assign(&mut self) {
        for (channel_idx, channel) in self.coefficients.iter_mut().enumerate() {
            let prime = self.basis.primes()[channel_idx];
            for coeff in channel.iter_mut() {
                if *coeff != 0 {
                    *coeff = prime - *coeff;
                }
            }
        }
    }
}

impl<const DEGREE: usize> PolySampler<DEGREE> for RnsNttPoly<DEGREE> {
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self {
        let max_prime = *context.primes().iter().max().unwrap_or(&1);
        let coeffs = uniform_coefficients::<DEGREE, _>(max_prime, rng);
        let mut poly = Self::from_u64_slice(&coeffs, context.clone());
        poly.to_ntt_domain();
        poly
    }

    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let modulus_product: u64 = context.primes().iter().product();
        let coeffs =
            gaussian_coefficients::<DEGREE, _>(std_dev, modulus_product, rng);
        let mut poly = Self::from_u64_slice(&coeffs, context.clone());
        poly.to_ntt_domain();
        poly
    }

    fn sample_tribits<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let ternary = ternary_coefficients::<DEGREE, _>(hamming_weight, rng);
        let mut poly = Self::from_i64_slice(&ternary, context.clone());
        poly.to_ntt_domain();
        poly
    }

    fn sample_noise<R: Rng>(
        variance: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        Self::sample_gaussian(variance.sqrt(), context, rng)
    }
}

impl<const DEGREE: usize> fmt::Display for RnsNttPoly<DEGREE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let domain = if self.in_ntt_domain { "NTT" } else { "COEFF" };
        write!(f, "RnsNttPoly<{DEGREE}>[{domain}](")?;
        let preview = DEGREE.min(4);
        for idx in 0..preview {
            if idx > 0 {
                write!(f, ", ")?;
            }
            write!(f, "[")?;
            for (channel_idx, channel) in self.coefficients.iter().enumerate() {
                if channel_idx > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", channel[idx])?;
            }
            write!(f, "]")?;
        }
        if DEGREE > preview {
            write!(f, ", …")?;
        }
        write!(f, ")")
    }
}

impl<const DEGREE: usize> PolyModSwitch<DEGREE> for RnsNttPoly<DEGREE> {
    fn mod_switch(&self, new_context: &Self::Context) -> Self {
        let current_count = self.channels();
        let new_count = new_context.channel_count();

        assert!(
            new_count <= current_count,
            "new basis cannot have more primes than current basis \
             (current {current_count}, new {new_count})"
        );

        for (idx, &new_prime) in new_context.primes().iter().enumerate() {
            let current_prime = self.basis.primes()[idx];
            assert_eq!(
                new_prime, current_prime,
                "new basis primes must match prefix; mismatch at index {idx}"
            );
        }

        let truncated = self.coefficients[..new_count].to_vec();

        Self {
            coefficients: truncated,
            in_ntt_domain: self.in_ntt_domain,
            basis: new_context.clone(),
        }
    }
}

fn reduce_signed_coeffs<const DEGREE: usize>(
    coeffs: &[i64; DEGREE],
    basis: &Arc<RnsBasis>,
) -> Vec<[u64; DEGREE]> {
    let mut channels = Vec::with_capacity(basis.channel_count());
    for &prime in basis.primes() {
        let modulus = prime as i64;
        let mut channel = [0u64; DEGREE];
        for (idx, &coeff) in coeffs.iter().enumerate() {
            channel[idx] = ((coeff % modulus + modulus) % modulus) as u64;
        }
        channels.push(channel);
    }
    channels
}

fn reduce_unsigned_coeffs<const DEGREE: usize>(
    coeffs: &[u64; DEGREE],
    basis: &Arc<RnsBasis>,
) -> Vec<[u64; DEGREE]> {
    let mut channels = Vec::with_capacity(basis.channel_count());
    for &prime in basis.primes() {
        let mut channel = [0u64; DEGREE];
        for (idx, &coeff) in coeffs.iter().enumerate() {
            channel[idx] = coeff % prime;
        }
        channels.push(channel);
    }
    channels
}

fn forward_ntt<const DEGREE: usize>(coeffs: &mut [u64; DEGREE], table: &NttTable) {
    let n = DEGREE;
    assert!(n.is_power_of_two(), "degree must be power of two");
    let mut t = n >> 1;
    let mut m = 1;
    while m < n {
        for i in 0..m {
            let j1 = 2 * i * t;
            let j2 = j1 + t - 1;
            let twiddle = table.forward_roots[m + i];
            for j in j1..=j2 {
                let u = coeffs[j];
                let v = mul_mod(coeffs[j + t], twiddle, table.prime);
                let sum = u + v;
                coeffs[j] = if sum >= table.prime {
                    sum - table.prime
                } else {
                    sum
                };
                coeffs[j + t] = if u >= v { u - v } else { table.prime - (v - u) };
            }
        }
        m <<= 1;
        t >>= 1;
    }
}

fn inverse_ntt<const DEGREE: usize>(coeffs: &mut [u64; DEGREE], table: &NttTable) {
    let n = DEGREE;
    assert!(n.is_power_of_two(), "degree must be power of two");
    let mut t = 1;
    let mut h = n >> 1;
    while h > 0 {
        let mut j1 = 0;
        for i in 0..h {
            let j2 = j1 + t - 1;
            let twiddle = table.inverse_roots[h + i];
            for j in j1..=j2 {
                let u = coeffs[j];
                let v = coeffs[j + t];
                let sum = u + v;
                coeffs[j] = if sum >= table.prime {
                    sum - table.prime
                } else {
                    sum
                };
                let diff = if u >= v { u - v } else { table.prime - (v - u) };
                coeffs[j + t] = mul_mod(diff, twiddle, table.prime);
            }
            j1 += 2 * t;
        }
        t <<= 1;
        h >>= 1;
    }
    for value in coeffs.iter_mut() {
        *value = mul_mod(*value, table.n_inv, table.prime);
    }
}

fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::backends::rns::{
        RnsPolyRing, params::toy_basis, params::toy_basis_with_channels,
    };
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn ntt_roundtrip_preserves_coefficients() {
        const DEGREE: usize = 16;
        let basis = toy_basis::<DEGREE>().expect("toy basis");
        let mut coeffs = [0i64; DEGREE];
        for i in 0..DEGREE {
            coeffs[i] = (i as i64 % 5) - 2;
        }

        let mut poly = RnsNttPoly::from_i64_slice(&coeffs, basis);
        let original: [i64; DEGREE] = poly.to_coeffs();
        poly.to_ntt_domain();
        poly.to_coeff_domain();
        assert_eq!(original, poly.to_coeffs());
    }

    #[test]
    fn ntt_multiplication_matches_schoolbook() {
        const DEGREE: usize = 16;
        let basis = toy_basis::<DEGREE>().expect("toy basis");
        let coeffs_a = [1i64; DEGREE];
        let coeffs_b = [2i64; DEGREE];

        let mut lhs_ntt: RnsNttPoly<DEGREE> =
            RnsNttPoly::from_i64_slice(&coeffs_a, basis.clone());
        let mut rhs_ntt: RnsNttPoly<DEGREE> =
            RnsNttPoly::from_i64_slice(&coeffs_b, basis.clone());
        lhs_ntt.to_ntt_domain();
        rhs_ntt.to_ntt_domain();
        let mut product_ntt = lhs_ntt.clone();
        product_ntt *= &rhs_ntt;
        let product_coeffs = product_ntt.to_coeffs();

        let mut lhs = RnsPolyRing::from_i64_slice(&coeffs_a, basis.clone());
        let rhs = RnsPolyRing::from_i64_slice(&coeffs_b, basis);
        lhs *= &rhs;
        assert_eq!(product_coeffs, lhs.to_coeffs());
    }

    #[test]
    fn ntt_addition_requires_domain() {
        const DEGREE: usize = 16;
        let basis = toy_basis_with_channels::<DEGREE>(2).expect("basis");
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let coeffs_a = uniform_coefficients::<DEGREE, _>(50, &mut rng);
        let coeffs_b = uniform_coefficients::<DEGREE, _>(50, &mut rng);

        let mut lhs: RnsNttPoly<DEGREE> =
            RnsNttPoly::from_u64_slice(&coeffs_a, basis.clone());
        let mut rhs: RnsNttPoly<DEGREE> =
            RnsNttPoly::from_u64_slice(&coeffs_b, basis);
        lhs.to_ntt_domain();
        rhs.to_ntt_domain();
        let mut sum = lhs.clone();
        sum += &rhs;
        let sum_coeffs = sum.to_coeffs();

        lhs.to_coeff_domain();
        rhs.to_coeff_domain();
        let mut expected = lhs.to_coeffs();
        let rhs_coeffs = rhs.to_coeffs();
        for (idx, value) in expected.iter_mut().enumerate() {
            *value += rhs_coeffs[idx];
        }
        assert_eq!(sum_coeffs, expected);
    }
}
