//! Fixed BigInt encoder for U256 backend
//!
//! This implementation follows the Vandermonde matrix approach like the working Python version,
//! ensuring proper f64 → U256 → f64 conversion flow with 50-bit scaling.

use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{Plaintext, PolyRing, PolyRingU256};
use crypto_bigint::U256;
use num_complex::Complex64;
use std::f64::consts::PI;

/// BigInt encoder parameters using Vandermonde matrix approach
#[derive(Debug, Clone)]
pub struct BigIntEncodingParams<const DEGREE: usize> {
    /// Number of bits used for scaling factor (should be 50 for your use case)
    scale_bits: u32,
    /// Precomputed primitive roots ξ^(2k+1) for k = 0..N/2-1
    xi_powers: Vec<Complex64>,
    /// The primitive M-th root of unity (ξ = e^(2πi/M) where M = 2*DEGREE)
    xi: Complex64,
}

/// BigInt encoder using Vandermonde matrix approach
pub struct BigIntEncoder<const DEGREE: usize> {
    params: BigIntEncodingParams<DEGREE>,
}

impl<const DEGREE: usize> BigIntEncodingParams<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }

        if scale_bits > 60 {
            return Err(EncodingError::InvalidInput {
                message: format!(
                    "Scale bits {} too large for BigInt encoder",
                    scale_bits
                ),
            });
        }

        let m = 2 * DEGREE; // M = 2N in CKKS

        // Compute primitive M-th root of unity: ξ = e^(2πi/M)
        let xi = Complex64::new(0.0, 2.0 * PI / (m as f64)).exp();

        // Precompute ξ^(2k+1) for k = 0..N/2-1 (the roots we evaluate at)
        let mut xi_powers = Vec::with_capacity(DEGREE / 2);
        for k in 0..(DEGREE / 2) {
            let power = 2 * k + 1;
            xi_powers.push(xi.powi(power as i32));
        }

        Ok(Self {
            scale_bits,
            xi_powers,
            xi,
        })
    }

    /// Get scaling factor as f64
    pub fn delta(&self) -> f64 {
        2.0_f64.powi(self.scale_bits as i32)
    }

    /// Maximum number of values that can be encoded (N/2 for CKKS)
    pub fn max_slots(&self) -> usize {
        DEGREE / 2
    }
}

impl<const DEGREE: usize> BigIntEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        let params = BigIntEncodingParams::new(scale_bits)?;
        Ok(Self { params })
    }

    /// Build Vandermonde matrix from xi^(2k+1) roots
    /// Each row i contains powers: `[1, root^1, root^2, ..., root^(N/2-1)]`
    fn build_vandermonde_matrix(&self) -> Vec<Vec<Complex64>> {
        let n = DEGREE / 2;
        let mut matrix = Vec::with_capacity(n);

        for root in &self.params.xi_powers {
            let mut row = Vec::with_capacity(n);
            let mut power = Complex64::new(1.0, 0.0); // root^0 = 1

            for _ in 0..n {
                row.push(power);
                power *= root; // root^j
            }
            matrix.push(row);
        }

        matrix
    }

    /// Solve linear system A * coeffs = values using Gaussian elimination
    /// This is the sigma_inverse operation from the Python code
    fn solve_vandermonde_system(
        &self,
        values: &[Complex64],
    ) -> EncodingResult<Vec<Complex64>> {
        let n = values.len();
        if n > DEGREE / 2 {
            return Err(EncodingError::InputTooLong {
                got: n,
                max: DEGREE / 2,
            });
        }

        // Build augmented matrix [A | b] where A is Vandermonde matrix
        let vandermonde = self.build_vandermonde_matrix();
        let mut augmented = Vec::with_capacity(n);

        for i in 0..n {
            let mut row = vandermonde[i].clone();
            row.push(values[i]); // Append the value as the last column
            augmented.push(row);
        }

        // Gaussian elimination with partial pivoting
        for i in 0..n {
            // Find pivot
            let mut max_row = i;
            for k in (i + 1)..n {
                if augmented[k][i].norm() > augmented[max_row][i].norm() {
                    max_row = k;
                }
            }

            // Swap rows if needed
            if max_row != i {
                augmented.swap(i, max_row);
            }

            // Check for singular matrix
            if augmented[i][i].norm() < 1e-12 {
                return Err(EncodingError::InvalidInput {
                    message: "Singular Vandermonde matrix".to_string(),
                });
            }

            // Make diagonal element 1
            let pivot = augmented[i][i];
            for j in 0..=n {
                augmented[i][j] /= pivot;
            }

            for k in 0..n {
                if k != i {
                    let factor = augmented[k][i];
                    for j in 0..=n {
                        let temp = factor * augmented[i][j]; // Store the value first
                        augmented[k][j] -= temp;
                    }
                }
            }
        }

        // Extract solution (last column)
        let mut coeffs = Vec::with_capacity(DEGREE / 2);
        for i in 0..n {
            coeffs.push(augmented[i][n]);
        }

        // Pad with zeros if needed
        while coeffs.len() < DEGREE / 2 {
            coeffs.push(Complex64::new(0.0, 0.0));
        }

        Ok(coeffs)
    }

    /// Encode f64 values to U256 coefficients using CKKS canonical embedding
    pub fn encode_values(&self, values: &[f64]) -> EncodingResult<Vec<U256>> {
        if values.len() > self.params.max_slots() {
            return Err(EncodingError::InputTooLong {
                got: values.len(),
                max: self.params.max_slots(),
            });
        }

        // Convert f64 to Complex64 (real values only)
        let mut complex_values = Vec::with_capacity(values.len());
        for &val in values {
            complex_values.push(Complex64::new(val, 0.0));
        }

        // Solve Vandermonde system to get polynomial coefficients
        let coeffs_complex = self.solve_vandermonde_system(&complex_values)?;

        // Extend to full polynomial with conjugate symmetry for real coefficients
        let mut full_coeffs = vec![Complex64::new(0.0, 0.0); DEGREE];

        // Set first N/2 coefficients
        for (i, &coeff) in coeffs_complex.iter().enumerate() {
            full_coeffs[i] = coeff;
        }

        // Apply conjugate symmetry: coeffs[N-i] = conjugate(coeffs[i]) for i > 0
        for i in 1..(DEGREE / 2) {
            if i < coeffs_complex.len() {
                full_coeffs[DEGREE - i] = coeffs_complex[i].conj();
            }
        }

        // Scale by 2^scale_bits and convert to U256
        let delta = self.params.delta();
        let mut result = Vec::with_capacity(DEGREE);

        for coeff in &full_coeffs {
            // Should be real due to conjugate symmetry, but take real part anyway
            let scaled_val = coeff.re * delta;

            // Check for overflow
            if scaled_val.abs() > (u128::MAX as f64) {
                return Err(EncodingError::CoefficientOutOfRange {
                    value: scaled_val,
                });
            }

            // Convert to U256 (handling negative values using two's complement)
            let rounded = scaled_val.round();
            let u256_val = if rounded >= 0.0 {
                U256::from(rounded as u128)
            } else {
                // For negative values, use modular arithmetic
                // This will be properly handled when we create the polynomial with modulus
                U256::from((-rounded) as u128)
            };

            result.push(u256_val);
        }

        Ok(result)
    }

    /// Decode U256 coefficients back to f64 values
    pub fn decode_coeffs(&self, coeffs: &[U256], modulus: U256) -> Vec<f64> {
        let half_modulus = modulus >> 1;
        let scale_bits = self.params.scale_bits as usize;

        // Convert U256 coefficients back to Complex64 with proper sign handling
        let mut complex_coeffs = Vec::with_capacity(DEGREE);

        for &coeff in coeffs.iter().take(DEGREE) {
            // Handle two's complement representation
            let signed_val = if coeff > half_modulus {
                // Negative number
                let positive_part = modulus - coeff;
                -(u256_to_f64(&positive_part))
            } else {
                u256_to_f64(&coeff)
            };

            // Descale by dividing by 2^scale_bits
            let descaled = signed_val / (2.0_f64.powi(scale_bits as i32));
            complex_coeffs.push(Complex64::new(descaled, 0.0));
        }

        // Evaluate polynomial at the primitive roots ξ^(2k+1)
        // This is the sigma operation from the Python code
        let mut result = Vec::with_capacity(self.params.max_slots());

        for root in &self.params.xi_powers {
            let mut eval = Complex64::new(0.0, 0.0);
            let mut power = Complex64::new(1.0, 0.0); // root^0

            // Evaluate polynomial: sum(coeff[j] * root^j) for j = 0..N-1
            for &coeff in &complex_coeffs {
                eval += coeff * power;
                power *= root;
            }

            // Take real part (should be real due to CKKS construction)
            result.push(eval.re);
        }

        result
    }
}

/// Convert U256 to f64 (helper function) WTF????
fn u256_to_f64(val: &U256) -> f64 {
    let words = val.as_words(); // Get u64 words
    let mut result = 0.0f64;
    let mut base = 1.0f64;

    for &word in words {
        result += (word as f64) * base;
        base *= 18446744073709551616.0; // 2^64
    }

    result
}

/// Convert f64 to U256 with proper overflow checking
fn f64_to_u256(val: f64) -> Result<U256, EncodingError> {
    if val < 0.0 || val > (u128::MAX as f64) {
        return Err(EncodingError::CoefficientOutOfRange { value: val });
    }

    Ok(U256::from(val as u128))
}

// Implement the Encoder trait
impl<const DEGREE: usize> Encoder<PolyRingU256<DEGREE>, DEGREE>
    for BigIntEncoder<DEGREE>
{
    fn encode(
        &self,
        values: &[f64],
        context: &<PolyRingU256<DEGREE> as PolyRing<DEGREE>>::Context,
    ) -> Plaintext<PolyRingU256<DEGREE>, DEGREE> {
        let coeffs = self.encode_values(values).expect("BigInt encoding failed");

        // coeffs are already U256, just handle modular arithmetic
        let u256_coeffs: Vec<U256> =
            coeffs.iter().map(|&c| c % context.get()).collect();

        let poly = PolyRingU256::<DEGREE>::from_u256_coeffs(&u256_coeffs, context);
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
        let modulus = plaintext.poly.modulus().get(); // .get() for NonZero<U256>
        self.decode_coeffs(&coeffs, modulus)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vandermonde_encoding_roundtrip() {
        let encoder = BigIntEncoder::<8>::new(50).unwrap(); // Use 50-bit scaling
        let values = vec![1.5, 2.5, 3.5, 4.5];

        let coeffs = encoder.encode_values(&values).unwrap();

        // Create a dummy modulus for testing (should be large enough)
        let modulus = U256::from(1u128 << 100); // Large test modulus
        let decoded = encoder.decode_coeffs(&coeffs, modulus);

        println!("Original: {:?}", values);
        println!("Decoded:  {:?}", &decoded[..values.len()]);

        for (orig, decoded) in values.iter().zip(decoded.iter()) {
            assert_relative_eq!(orig, decoded, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_vandermonde_matrix() {
        let encoder = BigIntEncoder::<8>::new(50).unwrap();
        let matrix = encoder.build_vandermonde_matrix();

        assert_eq!(matrix.len(), 4); // N/2 = 4 for DEGREE=8
        assert_eq!(matrix[0].len(), 4); // Each row has N/2 elements

        // First column should be all 1s (root^0)
        for row in &matrix {
            assert_relative_eq!(row[0].re, 1.0, epsilon = 1e-10);
            assert_relative_eq!(row[0].im, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_simple_encoding() {
        let encoder = BigIntEncoder::<8>::new(50).unwrap();

        // Test with simple values
        let values = vec![1.0];
        let coeffs = encoder.encode_values(&values).unwrap();

        println!("Coefficients for [1.0]: {:?}", coeffs);

        // Should not panic and should produce reasonable coefficients
        assert_eq!(coeffs.len(), 8);
    }

    #[test]
    fn debug_vandermonde() {
        let encoder = BigIntEncoder::<8>::new(50).unwrap();
        let matrix = encoder.build_vandermonde_matrix();

        // Print the matrix to see if it looks right
        for (i, row) in matrix.iter().enumerate() {
            println!("Row {}: {:?}", i, row);
        }
    }

    #[test]
    fn debug_roots() {
        let encoder = BigIntEncoder::<8>::new(50).unwrap();
        println!("Xi: {:?}", encoder.params.xi);
        for (i, root) in encoder.params.xi_powers.iter().enumerate() {
            println!("Root {}: {:?}", i, root);
        }
    }

    #[test]
    fn test_simple_case() {
        let encoder = BigIntEncoder::<8>::new(50).unwrap();
        let values = vec![1.0]; // Just one value
        let coeffs = encoder.encode_values(&values).unwrap();
        let decoded = encoder.decode_coeffs(&coeffs, U256::from_u8(17));

        println!("Input: {:?}", values);
        println!("Coeffs: {:?}", coeffs);
        println!("Decoded: {:?}", decoded);
    }
}

/* //! BigInt encoder for U256 backend
//!
//! This encoder implements CKKS encoding using direct mathematical DFT computation
//! optimized for large integer arithmetic, avoiding rustfft dependency.

use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{Plaintext, PolyRing, PolyRingU256};
use crypto_bigint::U256;
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

impl<const DEGREE: usize> Encoder<PolyRingU256<DEGREE>, DEGREE>
    for BigIntEncoder<DEGREE>
{
    fn encode(
        &self,
        values: &[f64],
        context: &<PolyRingU256<DEGREE> as PolyRing<DEGREE>>::Context,
    ) -> Plaintext<PolyRingU256<DEGREE>, DEGREE> {
        let coeffs = self.encode_values(values).expect("BigInt encoding failed");
        let poly = PolyRingU256::<DEGREE>::from_coeffs(&coeffs, context);
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
        let modulus = plaintext.poly.modulus().get(); // .get() for NonZero<U256>
        let scale_bits = self.params.scale_bits as usize;
        let half_modulus = modulus >> 1;

        coeffs
            .iter()
            .map(|c| {
                let (unsigned, sign) = if *c > half_modulus {
                    (modulus - *c, -1.0)
                } else {
                    (*c, 1.0)
                };
                let scaled = unsigned >> scale_bits;
                sign * u256_to_f64(&scaled)
            })
            .collect()
    }
}

fn u256_to_f64(val: &U256) -> f64 {
    // If you use the crypto_bigint crate:
    let words = val.as_words();
    let mut result = 0.0f64;
    let mut base = 1.0f64;
    for &w in words {
        result += (w as f64) * base;
        base *= 18446744073709551616.0; // 2^64
    }
    result
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

    #[test]
    fn test_final_encoding_comparison() {
        use crate::encoding::fft::{
            EncodingParams as FftParams, encode as fft_encode,
        };

        let values = vec![1.0];
        println!("=== Final Encoding Comparison ===");

        // Working FFT approach
        let fft_params = FftParams::<8>::new(30).unwrap();
        let fft_coeffs = fft_encode(&values, &fft_params).unwrap();

        // Your BigInt approach
        let bigint_encoder = BigIntEncoder::<8>::new(30).unwrap();
        let bigint_coeffs = bigint_encoder.encode_values(&values).unwrap();

        println!("FFT coeffs:    {:?}", fft_coeffs);
        println!("BigInt coeffs: {:?}", bigint_coeffs);

        // These SHOULD be very similar (both encoding [1.0])
        for (i, (&fft_c, &big_c)) in
            fft_coeffs.iter().zip(bigint_coeffs.iter()).enumerate()
        {
            let diff = (fft_c - big_c).abs();
            println!(
                "coeff[{}]: FFT={}, BigInt={}, diff={}",
                i, fft_c, big_c, diff
            );
        }
    }
} */
