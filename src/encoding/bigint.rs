//! BigInt encoder for U256 backend
//!
//! This encoder implements CKKS encoding using direct mathematical DFT computation
//! optimized for large integer arithmetic, avoiding rustfft dependency.

use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{Plaintext, PolyRing};
use std::f64::consts::PI;

/// BigInt encoder parameters optimized for large moduli
#[derive(Debug, Clone)]
pub struct BigIntEncodingParams<const DEGREE: usize> {
    /// Number of bits used for scaling factor
    scale_bits: u32,
    /// Precomputed roots of unity for DFT
    roots: Vec<(f64, f64)>, // (real, imag) pairs
}

/// BigInt encoder implementation for U256 and other large integer backends
pub struct BigIntEncoder<const DEGREE: usize> {
    params: BigIntEncodingParams<DEGREE>,
}

impl<const DEGREE: usize> BigIntEncodingParams<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }

        // For BigInt, we can handle larger scale values
        if scale_bits > 120 {
            return Err(EncodingError::InvalidInput {
                message: format!(
                    "Scale bits {} too large for BigInt encoder",
                    scale_bits
                ),
            });
        }

        // Precompute roots of unity: ω = e^(-2πi(2k+1)/(2N)) for k = 0..N/2-1
        let mut roots = Vec::with_capacity(DEGREE);
        for i in 0..DEGREE / 2 {
            let angle = -2.0 * PI * (2 * i + 1) as f64 / (2.0 * DEGREE as f64);
            let real = angle.cos();
            let imag = angle.sin();
            roots.push((real, imag));
        }
        // Add conjugate roots for second half
        for i in 0..DEGREE / 2 {
            let (real, imag) = roots[i];
            roots.push((real, -imag)); // conjugate
        }

        Ok(Self { scale_bits, roots })
    }

    /// Get scaling factor as f64 (for calculations)
    pub fn delta(&self) -> f64 {
        2.0_f64.powi(self.scale_bits as i32)
    }

    /// Maximum number of values that can be encoded
    pub fn max_slots(&self) -> usize {
        DEGREE / 2
    }
}

impl<const DEGREE: usize> BigIntEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        let params = BigIntEncodingParams::new(scale_bits)?;
        Ok(Self { params })
    }

    /// Custom DFT implementation for encoding
    /// Computes polynomial evaluation at roots of unity
    fn custom_dft(&self, values: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let n = values.len();
        let mut result = vec![(0.0, 0.0); n];

        for k in 0..n {
            let mut sum_real = 0.0;
            let mut sum_imag = 0.0;

            for j in 0..n {
                let (val_real, val_imag) = values[j];
                let (root_real, root_imag) = self.params.roots[k];

                // Compute root^j
                let mut power_real = 1.0;
                let mut power_imag = 0.0;
                for _ in 0..j {
                    let new_real = power_real * root_real - power_imag * root_imag;
                    let new_imag = power_real * root_imag + power_imag * root_real;
                    power_real = new_real;
                    power_imag = new_imag;
                }

                // Multiply value by root^j and accumulate
                sum_real += val_real * power_real - val_imag * power_imag;
                sum_imag += val_real * power_imag + val_imag * power_real;
            }

            result[k] = (sum_real, sum_imag);
        }

        result
    }

    /// Custom inverse DFT implementation for decoding
    /// Performs polynomial interpolation from evaluations
    fn custom_idft(&self, values: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let n = values.len();
        let mut result = vec![(0.0, 0.0); n];

        for k in 0..n {
            let mut sum_real = 0.0;
            let mut sum_imag = 0.0;

            for j in 0..n {
                let (val_real, val_imag) = values[j];
                let (root_real, root_imag) = self.params.roots[j];

                // Compute root^(-k) = conjugate(root)^k
                let conj_root_real = root_real;
                let conj_root_imag = -root_imag;

                let mut power_real = 1.0;
                let mut power_imag = 0.0;
                for _ in 0..k {
                    let new_real =
                        power_real * conj_root_real - power_imag * conj_root_imag;
                    let new_imag =
                        power_real * conj_root_imag + power_imag * conj_root_real;
                    power_real = new_real;
                    power_imag = new_imag;
                }

                // Multiply value by conjugate(root)^k and accumulate
                sum_real += val_real * power_real - val_imag * power_imag;
                sum_imag += val_real * power_imag + val_imag * power_real;
            }

            // Normalize by 1/N
            result[k] = (sum_real / n as f64, sum_imag / n as f64);
        }

        result
    }

    /// Encode values using mathematical DFT approach
    fn encode_values(&self, values: &[f64]) -> EncodingResult<[i64; DEGREE]> {
        let max_slots = self.params.max_slots();

        if values.len() > max_slots {
            return Err(EncodingError::InputTooLong {
                got: values.len(),
                max: max_slots,
            });
        }

        // Create complex input vector with conjugate symmetry
        let mut complex_values = vec![(0.0, 0.0); DEGREE];

        // Fill first half with input values (as real numbers)
        for (i, &val) in values.iter().enumerate() {
            complex_values[i] = (val, 0.0);
        }

        // Apply conjugate symmetry for second half to ensure real coefficients
        for i in 1..DEGREE / 2 {
            if i < values.len() {
                let (real, imag) = complex_values[i];
                complex_values[DEGREE - i] = (real, -imag); // conjugate
            }
        }

        // Apply inverse DFT to get polynomial coefficients
        let coeffs_complex = self.custom_idft(&complex_values);

        // Scale and convert to integers, taking only real parts
        let delta = self.params.delta();
        let mut result = [0i64; DEGREE];

        for (i, (real, _imag)) in coeffs_complex.iter().enumerate() {
            let scaled_val = real * delta;

            // Check for overflow in i64 range
            if scaled_val.abs() > i64::MAX as f64 / 2.0 {
                return Err(EncodingError::CoefficientOutOfRange {
                    value: scaled_val,
                });
            }

            result[i] = scaled_val.round() as i64;
        }

        Ok(result)
    }

    /// Decode coefficients back to values using mathematical DFT
    fn decode_coeffs(&self, coeffs: &[i64]) -> Vec<f64> {
        // Convert coefficients to complex numbers
        let mut complex_coeffs = vec![(0.0, 0.0); DEGREE];
        for (i, &coeff) in coeffs.iter().enumerate() {
            complex_coeffs[i] = (coeff as f64, 0.0);
        }

        // Apply DFT to evaluate polynomial at roots of unity
        let evaluations = self.custom_dft(&complex_coeffs);

        // Extract real parts and scale down
        let delta = self.params.delta();
        let max_slots = self.params.max_slots();

        evaluations
            .iter()
            .take(max_slots)
            .map(|(real, _imag)| real / delta)
            .collect()
    }
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Encoder<P, DEGREE>
    for BigIntEncoder<DEGREE>
{
    fn encode(&self, values: &[f64], context: &P::Context) -> Plaintext<P, DEGREE> {
        let coeffs = self.encode_values(values).expect("BigInt encoding failed");

        let poly = P::from_coeffs(&coeffs, context);

        Plaintext {
            poly,
            scale: self.params.delta(),
        }
    }

    fn decode(&self, plaintext: &Plaintext<P, DEGREE>) -> Vec<f64> {
        let coeffs = plaintext.poly.to_coeffs();
        self.decode_coeffs(&coeffs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bigint_encoding_roundtrip() {
        let encoder = BigIntEncoder::<8>::new(30).unwrap(); // Smaller scale for testing
        let values = vec![1.5, 2.5, 3.5];

        let coeffs = encoder.encode_values(&values).unwrap();
        let decoded = encoder.decode_coeffs(&coeffs);

        for (orig, decoded) in values.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig, decoded, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_roots_computation() {
        let params = BigIntEncodingParams::<8>::new(40).unwrap();
        assert_eq!(params.roots.len(), 8);

        // Verify roots are on unit circle
        for (real, imag) in &params.roots {
            let magnitude = (real * real + imag * imag).sqrt();
            assert_relative_eq!(magnitude, 1.0, epsilon = 1e-10);
        }
    }
}
