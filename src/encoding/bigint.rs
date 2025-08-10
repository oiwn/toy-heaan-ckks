#![allow(dead_code)]
use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{Plaintext, PolyRingU256};
use crypto_bigint::{NonZero, U256};
use num_complex::Complex64;
use rand::Rng;
use std::f64::consts::PI;

pub struct BigIntEncodingParams<const DEGREE: usize> {
    scale_bits: u32,
}

pub struct BigIntEncoder<const DEGREE: usize> {
    params: BigIntEncodingParams<DEGREE>,
}

#[derive(Clone)]
struct QuantizedPolyU256 {
    re: Vec<U256>,
    im: Vec<U256>,
}

impl<const DEGREE: usize> BigIntEncodingParams<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }
        assert!(scale_bits <= 52, "scale_bits must be ≤ 52 with f64 backend");
        Ok(Self { scale_bits })
    }

    #[inline]
    pub fn delta(&self) -> f64 {
        (1u64 << self.scale_bits) as f64
    }
}

impl<const DEGREE: usize> BigIntEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        Ok(Self {
            params: BigIntEncodingParams::new(scale_bits)?,
        })
    }

    /// π^{-1}: Expand C^{N/2} to Hermitian space H
    fn pi_inverse(&self, z: &[f64]) -> Vec<Complex64> {
        let n = z.len();
        let mut result = Vec::with_capacity(2 * n);

        // First N/2 elements: convert real to complex
        for &val in z {
            result.push(Complex64::new(val, 0.0));
        }

        // Second N/2 elements: complex conjugates in reverse order
        for &val in z.iter().rev() {
            result.push(Complex64::new(val, 0.0)); // Since input is real, conjugate is same
        }

        result
    }

    /// Create the σ(R) lattice basis
    fn create_sigma_r_basis(&self) -> Vec<Vec<Complex64>> {
        let n = DEGREE / 2;
        let xi = primitive_root(n as u32);

        // Vandermonde matrix: sigma(1), sigma(X), ..., sigma(X^{N-1})
        (0..DEGREE)
            .map(|j| {
                (0..n)
                    .map(|i| xi.powu(((2 * i + 1) * j) as u32))
                    .chain((0..n).map(|i| xi.powu(((2 * i + 1) * j) as u32).conj()))
                    .collect()
            })
            .collect()
    }

    /// Reconstruct vector from lattice basis coordinates
    fn reconstruct_from_basis(
        &self,
        coords: &[i64],
        basis: &[Vec<Complex64>],
    ) -> Vec<Complex64> {
        let mut result = vec![Complex64::new(0.0, 0.0); basis[0].len()];

        for (coord, basis_vec) in coords.iter().zip(basis.iter()) {
            let weight = *coord as f64;
            for (i, &basis_val) in basis_vec.iter().enumerate() {
                result[i] += basis_val * weight;
            }
        }

        result
    }
}

impl<const DEGREE: usize> Encoder<PolyRingU256<DEGREE>, DEGREE>
    for BigIntEncoder<DEGREE>
{
    fn encode(
        &self,
        values: &[f64],
        context: &NonZero<U256>,
    ) -> Plaintext<PolyRingU256<DEGREE>, DEGREE> {
        let n = DEGREE / 2;
        let mut padded = values.to_vec();
        padded.resize(n, 0.0);

        // For now, use the original vanilla implementation
        // TODO: Add full lattice projection in separate commit
        let complex_coeffs = sigma_inv(&padded, n as u32);
        let qp = quantize_u256(&complex_coeffs, context, self.params.delta());

        // Pack into polynomial coefficients
        let mut coeffs = [U256::ZERO; DEGREE];
        for i in 0..n {
            coeffs[i] = qp.re[i];
            coeffs[i + n] = qp.im[i];
        }

        let poly = PolyRingU256::<DEGREE>::from_u256_coeffs(&coeffs, context);
        Plaintext {
            poly,
            scale: self.params.delta(),
        }
    }

    fn decode(
        &self,
        plaintext: &Plaintext<PolyRingU256<DEGREE>, DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.coefficients();
        let modulus = plaintext.poly.modulus(); // NonZero<U256>

        let n = DEGREE / 2;

        let qp = QuantizedPolyU256 {
            re: coeffs[0..n].to_vec(),
            im: coeffs[n..2 * n].to_vec(),
        };

        let deq = dequantize_u256(&qp, modulus, n as u32, plaintext.scale);

        // let deq = dequantize_u256(&qp, modulus, n as u32, self.params.delta());
        let restored = sigma(&deq, n as u32);
        restored
    }
}

// xi = e^{2pi*i / (2*N)}
fn primitive_root(n: u32) -> Complex64 {
    let m = 2 * n;
    Complex64::from_polar(1.0, 2.0 * PI / (m as f64))
}

// r_i = xi^{2i+1}, i = 0..N-1
fn odd_roots(n: u32) -> Vec<Complex64> {
    let xi = primitive_root(n);
    (0..n).map(|i| xi.powu((2 * i + 1) as u32)).collect()
}

/// sigma^{-1}: b (len N, real) -> coeffs (len N, complex)
/// c_j = (1/N) * sum_i b_i * conj(r_i)^j
fn sigma_inv(b: &[f64], n: u32) -> Vec<Complex64> {
    assert_eq!(b.len(), n as usize);
    let roots = odd_roots(n);
    let mut coeffs = vec![Complex64::new(0.0, 0.0); n as usize];

    for j in 0..n as usize {
        let mut acc = Complex64::new(0.0, 0.0);
        for (i, &bi) in b.iter().enumerate() {
            let mut pow = Complex64::new(1.0, 0.0);
            let r_conj = roots[i].conj();
            for _ in 0..j {
                pow *= r_conj;
            }
            acc += pow * bi;
        }
        coeffs[j] = acc / n as f64;
    }
    coeffs
}

/// sigma: coeffs -> b (len N, real part)
fn sigma(coeffs: &[Complex64], n: u32) -> Vec<f64> {
    assert_eq!(coeffs.len(), n as usize);
    let roots = odd_roots(n);
    let mut out = vec![0.0f64; n as usize];

    for i in 0..n as usize {
        let mut acc = Complex64::new(0.0, 0.0);
        let mut pow = Complex64::new(1.0, 0.0);
        for (j, &c) in coeffs.iter().enumerate() {
            if j > 0 {
                pow *= roots[i];
            }
            acc += c * pow;
        }
        out[i] = acc.re;
    }
    out
}

fn quantize_u256(
    coeffs: &[Complex64],
    q: &NonZero<U256>,
    scale: f64,
) -> QuantizedPolyU256 {
    let re = coeffs
        .iter()
        .map(|c| to_mod_u256((c.re * scale).round() as i128, q))
        .collect();
    let im = coeffs
        .iter()
        .map(|c| to_mod_u256((c.im * scale).round() as i128, q))
        .collect();
    QuantizedPolyU256 { re, im }
}

fn dequantize_u256(
    qp: &QuantizedPolyU256,
    q: NonZero<U256>,
    n: u32,
    scale: f64,
) -> Vec<Complex64> {
    assert_eq!(qp.re.len(), n as usize);
    assert_eq!(qp.im.len(), n as usize);

    (0..n as usize)
        .map(|i| {
            let re = from_mod_u256_centered(qp.re[i], q.get()) as f64 / scale;
            let im = from_mod_u256_centered(qp.im[i], q.get()) as f64 / scale;
            Complex64::new(re, im)
        })
        .collect()
}

fn compute_basis_coordinates(
    z: &[Complex64],
    basis: &[Vec<Complex64>],
) -> Vec<f64> {
    basis
        .iter()
        .map(|b| {
            let dot_product: Complex64 =
                z.iter().zip(b.iter()).map(|(zi, bi)| zi * bi.conj()).sum();
            let norm_squared: f64 = b.iter().map(|bi| bi.norm_sqr()).sum();
            dot_product.re / norm_squared
        })
        .collect()
}

fn coordinate_wise_random_rounding(coords: &[f64], rng: &mut impl Rng) -> Vec<i64> {
    coords
        .iter()
        .map(|&c| {
            let floor_c = c.floor();
            let frac = c - floor_c;
            if rng.random::<f64>() < frac {
                floor_c as i64 + 1
            } else {
                floor_c as i64
            }
        })
        .collect()
}

#[inline]
fn to_mod_u256(v: i128, q: &NonZero<U256>) -> U256 {
    if v >= 0 {
        U256::from_u128(v as u128).rem_vartime(q)
    } else {
        let mag = U256::from_u128((-v) as u128).rem_vartime(q);
        q.get().wrapping_sub(&mag)
    }
}

#[inline]
fn from_mod_u256_centered(x: U256, q: U256) -> i128 {
    let half = q >> 1;
    if x > half {
        let mag = q.wrapping_sub(&x);
        -(u256_to_u128(mag) as i128)
    } else {
        u256_to_u128(x) as i128
    }
}

#[inline]
fn u256_to_u128(x: U256) -> u128 {
    let words = x.to_words();
    ((words[1] as u128) << 64) | (words[0] as u128)
}
