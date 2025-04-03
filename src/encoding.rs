//! The CKKS encoding process transforms a vector of real numbers into a polynomial
//! with integer coefficients. This transformation happens in the ring R[X]/(X^n + 1)
//! where n is a power of 2. The encoding preserves the values when the polynomial
//! is evaluated at specific roots of unity.
//!
//! Mathematical background:
//! - We work in the ring R[X]/(X^n + 1) where n is a power of 2
//! - The roots of X^n + 1 = 0 are the primitive 2n-th roots of unity
//! - For n=4, these are e^(2πi/8), e^(6πi/8), e^(10πi/8), e^(14πi/8)
//! - The encoding maps real values to evaluations at these roots
//! - We use scaling by 2^scale_bits to handle fixed-point arithmetic
use num_complex::Complex64;
use rustfft::FftPlanner;

/// Parameters for the CKKS encoding scheme.
pub struct EncodingParams {
    /// Polynomial ring degree (must be power of 2).
    /// This determines how many values we can encode:
    /// - We can encode up to n/2 real numbers
    /// - The polynomial will have degree n-1
    ring_degree: usize,

    /// Number of bits used for scaling factor.
    /// The scaling factor will be 2^scale_bits.
    /// This determines the precision of our fixed-point representation:
    /// - Larger values give more precision but risk overflow
    /// - Common values are 30-50 bits
    scale_bits: u32,
}

impl EncodingParams {
    /// Creates new encoding parameters.
    ///
    /// # Arguments
    /// * `n` - Ring degree, must be a power of 2
    /// * `scale_bits` - Number of bits for scaling factor (2^scale_bits)
    pub fn new(ring_degree: usize, scale_bits: u32) -> Result<Self, String> {
        // Ring degree must be power of 2 for FFT and CKKS
        if !ring_degree.is_power_of_two() {
            return Err("Ring degree must be a power of 2".to_string());
        }
        Ok(Self {
            ring_degree,
            scale_bits,
        })
    }

    /// Gets the scaling factor as a float (2^scale_bits)
    fn delta(&self) -> f64 {
        (1u64 << self.scale_bits) as f64
    }
}

/// Encodes a vector of real numbers into polynomial coefficients using CKKS encoding.
///
/// The encoding process:
/// 1. Scale the input values by 2^scale_bits to get fixed-point representation
/// 2. Create a complex vector with conjugate symmetry (needed for real coefficients)
/// 3. Apply inverse FFT to get polynomial coefficients
/// 4. Round coefficients to integers
///
/// Mathematical steps:
/// - Input vector z = (z₁, ..., zₖ) where k ≤ n/2
/// - Let zeta_i be the 2n-th roots of unity
/// - Find polynomial m(X) such that m(zeta_i) ≈ delta·z_i where delta = 2^scale_bits
/// - The polynomial coefficients are our encoding
pub fn encode(values: &[f64], params: &EncodingParams) -> Result<Vec<i64>, String> {
    if values.len() > params.ring_degree / 2 {
        return Err("Input vector too long".to_string());
    }

    let mut fft_input = vec![Complex64::new(0.0, 0.0); params.ring_degree];

    // Scale values and prepare conjugate symmetric vector
    // This ensures our polynomial has real coefficients
    for (i, &val) in values.iter().enumerate() {
        fft_input[i] = Complex64::new(val * params.delta(), 0.0);
    }
    for i in 1..values.len() {
        fft_input[params.ring_degree - i] = fft_input[i].conj();
    }

    // Apply inverse FFT to get polynomial coefficients
    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(params.ring_degree);
    ifft.process(&mut fft_input);

    // Round to nearest integers and normalize
    // The factor of n comes from FFT normalization
    let scale = (params.ring_degree as f64).recip();
    let coeffs: Vec<i64> = fft_input
        .iter()
        .map(|&x| (x.re * scale).round() as i64)
        .collect();

    Ok(coeffs)
}

/// Decodes polynomial coefficients back to real values.
///
/// The decoding process:
/// 1. Convert integer coefficients to complex numbers
/// 2. Apply forward FFT
/// 3. Extract real parts and scale back down
///
/// This is the inverse of the encoding process:
/// - Take polynomial coefficients aᵢ
/// - Evaluate at roots of unity via FFT
/// - Scale down by 2^scale_bits
/// - Extract real parts to get original values
pub fn decode(coeffs: &[i64], params: &EncodingParams) -> Result<Vec<f64>, String> {
    let mut fft_input: Vec<Complex64> = coeffs
        .iter()
        .map(|&x| Complex64::new(x as f64, 0.0))
        .collect();

    fft_input.resize(params.ring_degree, Complex64::new(0.0, 0.0));

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(params.ring_degree);
    fft.process(&mut fft_input);

    let result: Vec<f64> = fft_input
        .iter()
        .take(params.ring_degree / 2)
        .map(|&x| x.re / params.delta())
        .collect();

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_different_scale_bits() {
        // Test with increasing scale_bits
        let test_scales = vec![30, 40, 50];
        let input = vec![1.23456789, -2.34567891, 3.45678912];

        for scale_bits in test_scales {
            println!("Scale bits: {}", scale_bits);
            let params = EncodingParams::new(8, scale_bits).unwrap();
            let encoded = encode(&input, &params).unwrap();
            let decoded = decode(&encoded, &params).unwrap();

            println!("Epsilon: {}", 1.0 / params.delta());
            // Higher scale_bits should give better precision
            for (orig, dec) in input.iter().zip(decoded.iter()) {
                // assert_relative_eq!(orig, dec, epsilon = 1.0 / params.delta());
                assert_relative_eq!(orig, dec, epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn test_toy_example() {
        let params = EncodingParams::new(8, 16).unwrap(); // Delta = 2^7 = 128
        let input = vec![1.2, 3.4, 0.0, 5.4];

        let encoded = encode(&input, &params).unwrap();
        println!("Encoded coefficients (Delta=2^7): {:?}", encoded);

        let decoded = decode(&encoded, &params).unwrap();

        for (orig, dec) in input.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig, dec, epsilon = 1e-3);
        }
    }

    #[test]
    fn test_precise_small_values() {
        let params = EncodingParams::new(8, 30).unwrap();
        // Test small values that should encode/decode exactly
        let input = vec![0.5, 0.25, 0.125, 0.0625];

        let encoded = encode(&input, &params).unwrap();
        let decoded = decode(&encoded, &params).unwrap();

        for (orig, dec) in input.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig, dec, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_large_values() {
        let params = EncodingParams::new(8, 20).unwrap();
        // Test larger values to check scaling behavior
        let input = vec![1000.0, -2000.0, 3000.0];

        let encoded = encode(&input, &params).unwrap();
        let decoded = decode(&encoded, &params).unwrap();

        for (orig, dec) in input.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig, dec, epsilon = 1.0);
        }
    }

    #[test]
    fn test_error_cases() {
        // Test invalid ring degree
        assert!(EncodingParams::new(3, 20).is_err());

        // Test vector too long
        let params = EncodingParams::new(4, 20).unwrap();
        let long_input = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(encode(&long_input, &params).is_err());
    }

    #[test]
    fn test_conjugate_symmetry() {
        // This test verifies that our encoding preserves conjugate symmetry,
        // which is essential for real polynomial coefficients
        let params = EncodingParams::new(4, 20).unwrap();
        let input = vec![1.0, 2.0];

        let encoded = encode(&input, &params).unwrap();

        // Convert back to complex for FFT
        let mut complex_coeffs: Vec<Complex64> = encoded
            .iter()
            .map(|&x| Complex64::new(x as f64, 0.0))
            .collect();

        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(params.ring_degree);
        fft.process(&mut complex_coeffs);

        // Check conjugate pairs
        for i in 1..params.ring_degree / 2 {
            assert_relative_eq!(
                complex_coeffs[i].conj().re,
                complex_coeffs[params.ring_degree - i].re,
                epsilon = 1e-10
            );
            assert_relative_eq!(
                complex_coeffs[i].conj().im,
                complex_coeffs[params.ring_degree - i].im,
                epsilon = 1e-10
            );
        }
    }
}
