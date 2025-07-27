use crate::encoding::{Encoder, EncodingError, EncodingResult};
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
}
