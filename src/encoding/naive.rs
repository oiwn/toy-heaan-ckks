#![allow(dead_code)]
use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{NaivePolyRing, Plaintext, PolyRing};
use num_complex::Complex64;
use std::f64::consts::PI;

pub struct NaiveEncodingParams<const DEGREE: usize> {
    scale_bits: u32,
}

pub struct NaiveEncoder<const DEGREE: usize> {
    params: NaiveEncodingParams<DEGREE>,
}

#[derive(Clone)]
struct QuantizedPolyU64 {
    re: Vec<u64>,
    im: Vec<u64>,
}

impl<const DEGREE: usize> NaiveEncodingParams<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }
        // For u64, we can go higher than 52 bits since we're not limited by f64 mantissa
        assert!(scale_bits <= 60, "scale_bits must be ≤ 60 with u64 backend");
        Ok(Self { scale_bits })
    }

    #[inline]
    pub fn delta(&self) -> f64 {
        (1u64 << self.scale_bits) as f64
    }
}

impl<const DEGREE: usize> NaiveEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        Ok(Self {
            params: NaiveEncodingParams::new(scale_bits)?,
        })
    }

    /// π^{-1}: Expand C^{N/2} to Hermitian space H (same as BigInt version)
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

    /// Create the σ(R) lattice basis (same structure as BigInt)
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

impl<const DEGREE: usize> Encoder<NaivePolyRing<DEGREE>, DEGREE>
    for NaiveEncoder<DEGREE>
{
    fn encode(
        &self,
        values: &[f64],
        context: &u64, // NaivePolyRing uses u64 modulus
    ) -> Plaintext<NaivePolyRing<DEGREE>, DEGREE> {
        let n = DEGREE / 2;
        let mut padded = values.to_vec();
        padded.resize(n, 0.0);

        // For now, use the vanilla implementation (same as BigInt approach)
        let complex_coeffs = sigma_inv(&padded, n as u32);
        let qp = quantize_u64(&complex_coeffs, *context, self.params.delta());

        // Pack into polynomial coefficients: [re0, re1, ..., im0, im1, ...]
        let mut coeffs = [0u64; DEGREE];
        for i in 0..n {
            coeffs[i] = qp.re[i];
            coeffs[i + n] = qp.im[i];
        }

        // Convert u64 to i64 for NaivePolyRing
        let i64_coeffs = coeffs.map(|x| x as i64);
        let poly = NaivePolyRing::<DEGREE>::from_coeffs(&i64_coeffs, context);

        Plaintext {
            poly,
            scale: self.params.delta(),
        }
    }

    fn decode(
        &self,
        plaintext: &Plaintext<NaivePolyRing<DEGREE>, DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.coeffs; // [u64; DEGREE] - already u64!
        let modulus = plaintext.poly.context();

        let n = DEGREE / 2;

        // Coefficients are already u64, no conversion needed
        let qp = QuantizedPolyU64 {
            re: coeffs[0..n].to_vec(),
            im: coeffs[n..2 * n].to_vec(),
        };

        let deq = dequantize_u64(&qp, *modulus, n as u32, self.params.delta());

        sigma(&deq, n as u32)
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
    (0..n).map(|i| xi.powu(2 * i + 1)).collect()
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

fn quantize_u64(coeffs: &[Complex64], q: u64, scale: f64) -> QuantizedPolyU64 {
    let re = coeffs
        .iter()
        .map(|c| to_mod_u64((c.re * scale).round() as i64, q))
        .collect();
    let im = coeffs
        .iter()
        .map(|c| to_mod_u64((c.im * scale).round() as i64, q))
        .collect();
    QuantizedPolyU64 { re, im }
}

fn dequantize_u64(
    qp: &QuantizedPolyU64,
    q: u64,
    n: u32,
    scale: f64,
) -> Vec<Complex64> {
    assert_eq!(qp.re.len(), n as usize);
    assert_eq!(qp.im.len(), n as usize);

    (0..n as usize)
        .map(|i| {
            let re = from_mod_u64_centered(qp.re[i], q) as f64 / scale;
            let im = from_mod_u64_centered(qp.im[i], q) as f64 / scale;
            Complex64::new(re, im)
        })
        .collect()
}

#[inline]
fn to_mod_u64(v: i64, q: u64) -> u64 {
    let q_i64 = q as i64;
    if v >= 0 {
        (v % q_i64) as u64
    } else {
        let mag = (-v) % q_i64;
        q - (mag as u64)
    }
}

#[inline]
fn from_mod_u64_centered(x: u64, q: u64) -> i64 {
    let half = q / 2;
    if x > half {
        let mag = q - x;
        -(mag as i64)
    } else {
        x as i64
    }
}

// Helper functions for lattice projection (unused for now, but ready for Part 2 CKKS)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naive_encoder_roundtrip() {
        let scale_bits = 30;
        let modulus = (1u64 << 50) - 27; // Same as your working example

        let encoder = NaiveEncoder::<8>::new(scale_bits).unwrap();

        let values = vec![1.0, 2.0, 3.0, 4.0];

        // Encode f64 → u64 polynomial
        let plaintext = encoder.encode(&values, &modulus);

        // Decode u64 polynomial → f64
        let decoded = encoder.decode(&plaintext);

        println!("Original: {:?}", values);
        println!("Decoded:  {:?}", &decoded[..values.len()]);

        // Check accuracy
        for (orig, dec) in values.iter().zip(decoded.iter()) {
            let error = (orig - dec).abs();
            assert!(
                error < 1e-10,
                "Error too large: {} vs {}, diff: {}",
                orig,
                dec,
                error
            );
        }
    }

    #[test]
    fn test_quantization_roundtrip() {
        let q = (1u64 << 50) - 27;
        let scale = (1u64 << 30) as f64;

        let complex_vals = vec![
            Complex64::new(1.5, -2.3),
            Complex64::new(0.0, 4.7),
            Complex64::new(-3.1, 0.0),
        ];

        let quantized = quantize_u64(&complex_vals, q, scale);
        let dequantized = dequantize_u64(&quantized, q, 3, scale);

        for (orig, restored) in complex_vals.iter().zip(dequantized.iter()) {
            let re_error = (orig.re - restored.re).abs();
            let im_error = (orig.im - restored.im).abs();

            assert!(re_error < 1e-9, "Real part error too large: {}", re_error);
            assert!(im_error < 1e-9, "Imag part error too large: {}", im_error);
        }
    }
}
/* use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{NaivePolyRing, Plaintext, PolyRing};
use num_complex::Complex64;
use std::f64::consts::PI;

pub struct EncodingParams<const DEGREE: usize> {
    scale_bits: u32,
}

pub struct NaiveEncoder<const DEGREE: usize> {
    params: EncodingParams<DEGREE>,
}

#[derive(Clone)]
struct QuantizedPoly {
    re: Vec<u64>,
    im: Vec<u64>,
}

impl<const DEGREE: usize> EncodingParams<DEGREE> {
    /// Creates new encoding parameters.
    ///
    /// # Arguments
    /// * `scale_bits` - Number of bits for scaling factor (2^scale_bits)
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        // Ring degree must be power of 2 for FFT and CKKS
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }
        Ok(Self { scale_bits })
    }

    /// Gets the scaling factor as a float (2^scale_bits)
    pub fn delta(&self) -> f64 {
        (1u64 << self.scale_bits) as f64
    }
}

impl<const DEGREE: usize> NaiveEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        // Return Result<Self>
        let params = EncodingParams::new(scale_bits)?;
        Ok(Self { params })
    }

    fn primitive_root(&self) -> Complex64 {
        let m = DEGREE;
        Complex64::from_polar(1.0, 2.0 * PI / (m as f64))
    }
}

impl<const DEGREE: usize> Encoder<NaivePolyRing<DEGREE>, DEGREE>
    for NaiveEncoder<DEGREE>
{
    fn encode(
        &self,
        values: &[f64],
        context: &<NaivePolyRing<DEGREE> as PolyRing<DEGREE>>::Context,
    ) -> Plaintext<NaivePolyRing<DEGREE>, DEGREE> {
        let n = DEGREE / 2; // Number of slots

        // Pad values to n length if needed
        let mut padded_values = values.to_vec();
        padded_values.resize(n, 0.0);

        // Apply sigma_inv to get complex coefficients
        let complex_coeffs = sigma_inv(&padded_values, n as u32);

        // Quantize to get integer coefficients
        let qp = quantize(&complex_coeffs, *context, self.params.delta());

        // Convert QuantizedPoly to polynomial coefficients
        let mut poly_coeffs = [0u64; DEGREE];
        for i in 0..n {
            poly_coeffs[i] = qp.re[i];
            poly_coeffs[i + n] = qp.im[i];
        }

        let poly =
            NaivePolyRing::from_coeffs(&poly_coeffs.map(|x| x as i64), context);

        Plaintext {
            poly,
            scale: self.params.delta(),
        }
    }

    fn decode(
        &self,
        plaintext: &Plaintext<NaivePolyRing<DEGREE>, DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.coeffs;
        let modulus = plaintext.poly.context();

        let n = DEGREE / 2; // Number of slots (N/2 for CKKS)
        let qp = QuantizedPoly {
            re: coeffs[0..n].to_vec(),
            im: coeffs[n..2 * n].to_vec(), // Assuming second half is imaginary
        };

        let dequantized = dequantize(&qp, *modulus, n as u32, self.params.delta());
        let restored = sigma(&dequantized, n as u32);
        restored
    }
}

// xi = e^{2pi*i / (2*N)}
fn primitive_root(n: u32) -> Complex64 {
    let m = 2 * n;
    Complex64::from_polar(1.0, 2.0 * std::f64::consts::PI / (m as f64))
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
            // conj(r_i)^j
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

fn quantize(coeffs: &[Complex64], q: u64, scale: f64) -> QuantizedPoly {
    let q = q as i64;
    let re = coeffs.iter().map(|c| to_mod_q(c.re * scale, q)).collect();
    let im = coeffs.iter().map(|c| to_mod_q(c.im * scale, q)).collect();
    QuantizedPoly { re, im }
}

fn dequantize(qp: &QuantizedPoly, q: u64, n: u32, scale: f64) -> Vec<Complex64> {
    assert_eq!(qp.re.len(), n as usize);
    assert_eq!(qp.im.len(), n as usize);
    (0..n as usize)
        .map(|i| {
            let re = from_mod_q(qp.re[i], q) as f64 / scale;
            let im = from_mod_q(qp.im[i], q) as f64 / scale;
            Complex64::new(re, im)
        })
        .collect()
}

#[inline]
fn to_mod_q(x: f64, q: i64) -> u64 {
    let mut r = x.round() as i64 % q;
    if r < 0 {
        r += q;
    }
    r as u64
}

#[inline]
fn from_mod_q(x: u64, q: u64) -> i64 {
    let q_i = q as i64;
    let half = (q >> 1) as i64;
    let mut v = x as i64;
    if v > half {
        v -= q_i;
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naive_encoder_roundtrip() {
        // Test parameters similar to your working example
        let scale_bits = 30;
        let modulus = (1u64 << 50) - 27; // Same as your working example

        let encoder = NaiveEncoder::<8>::new(scale_bits).unwrap();
        let context = modulus;

        // Test values
        let values = vec![1.0, 2.0, 3.0, 4.0];

        // Encode
        let plaintext = encoder.encode(&values, &context);

        // Decode
        let decoded = encoder.decode(&plaintext);

        println!("Original: {:?}", values);
        println!("Decoded:  {:?}", &decoded[..values.len()]);

        // Check accuracy (same tolerance as your working example)
        for (orig, dec) in values.iter().zip(decoded.iter()) {
            let error = (orig - dec).abs();
            assert!(
                error < 1e-10,
                "Error too large: {} vs {}, diff: {}",
                orig,
                dec,
                error
            );
        }
    }
} */
