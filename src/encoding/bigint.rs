use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crypto_bigint::{BigInt, ToBigInt};
use num_complex::Complex64;
use std::f64::consts::PI;

pub struct BigIntEncoder<const DEGREE: usize> {
    scale_bits: u32,
    scale: BigInt,
}

impl<const DEGREE: usize> BigIntEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }

        let scale = BigInt::from(2u64).pow(scale_bits);
        Ok(Self { scale_bits, scale })
    }

    /// Prepare roots: ζⱼ = e^(-2πi(2j+1)/(2N)) for j = 0..N/2
    fn prepare_roots(&self) -> Vec<Complex64> {
        let mut roots = vec![Complex64::new(0.0, 0.0); DEGREE];

        for i in 0..(DEGREE / 2) {
            let angle = -2.0 * PI * (2 * i + 1) as f64 / (2.0 * DEGREE as f64);
            roots[i] = Complex64::from_polar(1.0, angle);
            roots[DEGREE - i - 1] = roots[i].conj();
        }
        roots
    }

    /// Custom inverse DFT using Vandermonde matrix approach
    fn custom_idft(
        &self,
        values: &[Complex64],
        roots: &[Complex64],
    ) -> Vec<Complex64> {
        let n = values.len();
        let mut result = vec![Complex64::new(0.0, 0.0); n];

        // W^T @ values / N where W[n][k] = roots[n]^(-k)
        for i in 0..n {
            for k in 0..n {
                let power = roots[i].powc(Complex64::new(-(k as f64), 0.0));
                result[i] += power * values[k];
            }
            result[i] /= n as f64;
        }
        result
    }

    /// Custom DFT using Vandermonde matrix approach  
    fn custom_dft(&self, coeffs: &[i64], roots: &[Complex64]) -> Vec<Complex64> {
        let n = coeffs.len();
        let mut result = vec![Complex64::new(0.0, 0.0); n];

        // W @ coeffs where W[n][k] = roots[n]^k
        for i in 0..n {
            for k in 0..n {
                let power = roots[i].powc(Complex64::new(k as f64, 0.0));
                result[i] += power * (coeffs[k] as f64);
            }
        }
        result
    }
}

impl<const DEGREE: usize> Encoder<DEGREE> for BigIntEncoder<DEGREE> {
    fn encode(&self, values: &[f64]) -> EncodingResult<Vec<i64>> {
        if values.len() > DEGREE / 2 {
            return Err(EncodingError::InputTooLong {
                got: values.len(),
                max: DEGREE / 2,
            });
        }

        let roots = self.prepare_roots();

        // Create complex vector with conjugate symmetry
        let mut complex_values = vec![Complex64::new(0.0, 0.0); DEGREE];
        for (i, &val) in values.iter().enumerate() {
            complex_values[i] = Complex64::new(val, 0.0);
        }

        // Add conjugate values for real polynomial coefficients
        for i in 0..values.len() {
            complex_values[DEGREE - i - 1] = complex_values[i].conj();
        }

        // Apply inverse DFT to get polynomial coefficients
        let coeffs_complex = self.custom_idft(&complex_values, &roots);

        // Scale and convert to integers
        let mut result = Vec::with_capacity(DEGREE);
        for coeff in coeffs_complex {
            let scaled = coeff.re * self.scale.to_f64().unwrap();
            let rounded = scaled.round() as i64;

            // Handle negative numbers properly
            if rounded < 0 {
                result.push(rounded);
            } else {
                result.push(rounded);
            }
        }

        Ok(result)
    }

    fn decode(&self, coeffs: &[i64]) -> EncodingResult<Vec<f64>> {
        if coeffs.len() != DEGREE {
            return Err(EncodingError::InvalidInput {
                message: format!(
                    "Expected {} coefficients, got {}",
                    DEGREE,
                    coeffs.len()
                ),
            });
        }

        let roots = self.prepare_roots();
        let scale_f64 = self.scale.to_f64().unwrap();

        // Apply DFT to evaluate polynomial at roots
        let values_complex = self.custom_dft(coeffs, &roots);

        // Extract first N/2 values and scale down
        let mut result = Vec::with_capacity(DEGREE / 2);
        for i in 0..(DEGREE / 2) {
            let decoded = values_complex[i].re / scale_f64;
            result.push(decoded);
        }

        Ok(result)
    }

    fn scale(&self) -> f64 {
        self.scale.to_f64().unwrap()
    }
}
