use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{BigIntPolyRing, Plaintext};
use crypto_bigint::{NonZero, U256};
use num_complex::Complex64;
use std::f64::consts::PI;

pub struct BigIntEncodingParams<const DEGREE: usize> {
    scale_bits: u32,
    rot_group: Vec<u64>,
    ksi_pows: Vec<Complex64>,
}

pub struct BigIntEncoder<const DEGREE: usize> {
    params: BigIntEncodingParams<DEGREE>,
}

impl<const DEGREE: usize> BigIntEncodingParams<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }
        // assert!(scale_bits <= 52, "scale_bits must be ≤ 52 with f64 backend");

        // Initialize Kim's rotation group and roots of unity
        let (rot_group, ksi_pows) = Self::init_kim_structures();

        Ok(Self {
            scale_bits,
            rot_group,
            ksi_pows,
        })
    }

    /// Kim's initialization: rotGroup and ksiPows
    fn init_kim_structures() -> (Vec<u64>, Vec<Complex64>) {
        let n = DEGREE; // N in Kim's code
        let nh = n / 2; // Nh in Kim's code
        let m = n * 2; // M in Kim's code

        // rotGroup: powers of 5 mod M, Kim's rotation group
        let mut rot_group = vec![0u64; nh];
        let mut five_pows = 1u64;
        for i in 0..nh {
            rot_group[i] = five_pows;
            five_pows = (five_pows * 5) % (m as u64);
        }

        // ksiPows: e^(2πij/M) for j = 0..M
        let mut ksi_pows = vec![Complex64::new(0.0, 0.0); m + 1];
        for j in 0..m {
            let angle = 2.0 * PI * (j as f64) / (m as f64);
            ksi_pows[j] = Complex64::new(angle.cos(), angle.sin());
        }
        ksi_pows[m] = ksi_pows[0]; // Wraparound

        (rot_group, ksi_pows)
    }
}

impl<const DEGREE: usize> BigIntEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        Ok(Self {
            params: BigIntEncodingParams::new(scale_bits)?,
        })
    }

    /// Kim's bit reversal algorithm for FFT
    fn bit_reverse(&self, vals: &mut [Complex64]) {
        let size = vals.len();
        let mut j = 0usize;

        for i in 1..size {
            let mut bit = size >> 1;
            while j >= bit {
                j -= bit;
                bit >>= 1;
            }
            j += bit;
            if i < j {
                vals.swap(i, j);
            }
        }
    }

    /// Kim's fftSpecialInvLazy implementation
    fn fft_special_inv_lazy(&self, vals: &mut [Complex64]) {
        let size = vals.len();
        let m = DEGREE * 2; // M in Kim's code

        let mut len = size;
        while len >= 1 {
            for i in (0..size).step_by(len) {
                let lenh = len >> 1;
                let lenq = len << 2;

                for j in 0..lenh {
                    let rot_idx = self.params.rot_group[j] % (lenq as u64);
                    let idx = ((lenq as u64 - rot_idx) * (m as u64) / (lenq as u64))
                        as usize;

                    let u = vals[i + j] + vals[i + j + lenh];
                    let mut v = vals[i + j] - vals[i + j + lenh];
                    v *= self.params.ksi_pows[idx];

                    vals[i + j] = u;
                    vals[i + j + lenh] = v;
                }
            }
            len >>= 1;
        }

        self.bit_reverse(vals);
    }

    /// Kim's fftSpecialInv with normalization
    fn fft_special_inv(&self, vals: &mut [Complex64]) {
        println!("fft_special_inv: input size = {}", vals.len());
        self.fft_special_inv_lazy(vals);
        println!("fft_special_inv: after lazy, normalizing...");
        let size = vals.len() as f64;
        for val in vals.iter_mut() {
            *val /= size;
        }
        println!("fft_special_inv: done");
    }

    /// Kim's fftSpecial for decoding
    fn fft_special(&self, vals: &mut [Complex64]) {
        let size = vals.len();
        let m = DEGREE * 2; // M in Kim's code

        self.bit_reverse(vals);

        let mut len = 2;
        while len <= size {
            for i in (0..size).step_by(len) {
                let lenh = len >> 1;
                let lenq = len << 2;

                for j in 0..lenh {
                    let rot_idx = self.params.rot_group[j] % (lenq as u64);
                    let idx = (rot_idx * (m as u64) / (lenq as u64)) as usize;

                    let u = vals[i + j];
                    let mut v = vals[i + j + lenh];
                    v *= self.params.ksi_pows[idx];

                    vals[i + j] = u + v;
                    vals[i + j + lenh] = u - v;
                }
            }
            len <<= 1;
        }
    }

    fn scale_down_to_real(
        &self,
        coeff: U256,
        logp: u32,
        modulus: &NonZero<U256>,
    ) -> f64 {
        let q = modulus.get();
        let half_q = q >> 1;

        // Convert to centered representation: map [0, q) -> [-(q/2), q/2)
        let centered_val = if coeff > half_q {
            // Negative: coeff - q (wrapped subtraction)
            let diff = q.wrapping_sub(&coeff);
            -(diff.to_words()[0] as f64) // Use only lower word for now
        } else {
            // Positive
            coeff.to_words()[0] as f64
        };

        // Apply Kim's approach: divide by 2^logp (scale factor)
        let scale_factor = (1u64 << logp) as f64;
        centered_val / scale_factor
    }

    /// Improved scale up that handles edge cases better
    fn scale_up_to_u256(&self, val: f64, logp: u32) -> U256 {
        let scale_factor = (1u64 << logp) as f64;
        let scaled = (val * scale_factor).round();

        if scaled >= 0.0 {
            // Positive value
            let uint_val = scaled as u64;
            U256::from_u64(uint_val)
        } else {
            // Negative value: compute modulus - |scaled|
            let abs_val = (-scaled) as u64;
            let abs_u256 = U256::from_u64(abs_val);
            // Note: This is a simplified version. In practice you'd need the actual modulus here
            U256::ZERO.wrapping_sub(&abs_u256)
        }
    }
}

impl<const DEGREE: usize> Encoder<BigIntPolyRing<DEGREE>, DEGREE>
    for BigIntEncoder<DEGREE>
{
    fn encode(
        &self,
        values: &[f64],
        context: &NonZero<U256>,
    ) -> Plaintext<BigIntPolyRing<DEGREE>, DEGREE> {
        let n = DEGREE / 2; // Nh in Kim's code
        let slots = values.len();

        // Assert no padding as requested
        assert!(slots <= n, "Too many values: got {}, max {}", slots, n);

        // Step 1: Convert real values to complex (imag = 0.0)
        let mut uvals = vec![Complex64::new(0.0, 0.0); slots];
        for (i, &val) in values.iter().enumerate() {
            uvals[i] = Complex64::new(val, 0.0); // Real part only
        }

        // Step 2: Kim's fftSpecialInv
        self.fft_special_inv(&mut uvals);

        // Step 3: Kim's gap-based coefficient distribution
        let mut coeffs = [U256::ZERO; DEGREE];
        let gap = n / slots;

        for i in 0..slots {
            let idx = i * gap; // Real parts: 0, gap, 2*gap, ...
            let jdx = n + i * gap; // Imag parts: Nh, Nh+gap, Nh+2*gap, ...

            coeffs[idx] =
                self.scale_up_to_u256(uvals[i].re, self.params.scale_bits);
            coeffs[jdx] =
                self.scale_up_to_u256(uvals[i].im, self.params.scale_bits);
        }

        let poly = BigIntPolyRing::<DEGREE>::from_u256_coeffs(&coeffs, context);
        Plaintext {
            poly,
            scale_bits: self.params.scale_bits,
        }
    }

    fn decode(
        &self,
        plaintext: &Plaintext<BigIntPolyRing<DEGREE>, DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.coefficients();
        let modulus = plaintext.poly.modulus();
        let n = DEGREE / 2;

        let logp = plaintext.scale_bits;

        println!("Decode debug: scale_bits={}", logp);

        let slots = n;
        let gap = n / slots;
        let mut uvals = vec![Complex64::new(0.0, 0.0); slots];

        for i in 0..slots {
            let idx = i * gap;
            let jdx = n + i * gap;

            let real_val = self.scale_down_to_real(coeffs[idx], logp, &modulus);
            let imag_val = self.scale_down_to_real(coeffs[jdx], logp, &modulus);

            println!(
                "  uvals[{}]: coeffs[{}]={:x} -> {:.6}, coeffs[{}]={:x} -> {:.6}",
                i,
                idx,
                coeffs[idx].to_words()[0],
                real_val,
                jdx,
                coeffs[jdx].to_words()[0],
                imag_val
            );

            uvals[i] = Complex64::new(real_val, imag_val);
        }

        self.fft_special(&mut uvals);
        let result: Vec<f64> = uvals.iter().map(|c| c.re).collect();
        println!("After fft_special: {:?}", result);
        result
    }

    /* fn decode(
        &self,
        plaintext: &Plaintext<BigIntPolyRing<DEGREE>, DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.coefficients();
        let modulus = plaintext.poly.modulus();
        let n = DEGREE / 2;

        // Calculate logp from the actual plaintext scale
        let actual_scale = plaintext.scale;
        let logp = (actual_scale.log2().round() as u32); // Convert scale back to logp

        let slots = n;
        let gap = n / slots;
        let mut uvals = vec![Complex64::new(0.0, 0.0); slots];

        for i in 0..slots {
            let idx = i * gap;
            let jdx = n + i * gap;

            let real_val = self.scale_down_to_real(coeffs[idx], logp, &modulus); // ✅ Use actual logp
            let imag_val = self.scale_down_to_real(coeffs[jdx], logp, &modulus); // ✅ Use actual logp

            uvals[i] = Complex64::new(real_val, imag_val);
        }

        self.fft_special(&mut uvals);
        uvals.iter().map(|c| c.re).collect()
    } */

    /* fn decode(
        &self,
        plaintext: &Plaintext<BigIntPolyRing<DEGREE>, DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.coefficients();
        let modulus = plaintext.poly.modulus();
        let n = DEGREE / 2;

        // Extract coefficients back with Kim's gap structure
        // For now, assume we encoded with gap=1 (full slots)
        // In a complete implementation, we'd need to track the original gap
        let slots = n; // Simplified assumption
        let gap = n / slots;

        let mut uvals = vec![Complex64::new(0.0, 0.0); slots];

        for i in 0..slots {
            let idx = i * gap;
            let jdx = n + i * gap;

            let real_val = self.scale_down_to_real(
                coeffs[idx],
                self.params.scale_bits,
                &modulus,
            );
            let imag_val = self.scale_down_to_real(
                coeffs[jdx],
                self.params.scale_bits,
                &modulus,
            );

            uvals[i] = Complex64::new(real_val, imag_val);
        }

        // Apply Kim's fftSpecial for decoding
        self.fft_special(&mut uvals);

        // Extract real parts (since we encoded real values)
        uvals.iter().map(|c| c.re).collect()
    } */
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::U256;

    #[test]
    fn test_bigint_encoder_roundtrip() {
        const DEGREE: usize = 8;
        const SCALE_BITS: u32 = 20;

        let encoder = BigIntEncoder::<DEGREE>::new(SCALE_BITS).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(1152921504606846883u128)).unwrap();

        // Start with just 1 value to debug
        let values = vec![1.5];
        println!("Starting encode with values: {:?}", values);

        let plaintext = encoder.encode(&values, &modulus);
        println!("Encode completed, starting decode...");

        let decoded = encoder.decode(&plaintext);
        println!("Decoded: {:?}", decoded);
    }

    #[test]
    fn test_max_slots_assertion() {
        const DEGREE: usize = 8;
        let encoder = BigIntEncoder::<DEGREE>::new(20).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(1152921504606846883u128)).unwrap();

        // Should work: 4 slots for DEGREE=8
        let values_ok = vec![1.0, 2.0, 3.0, 4.0];
        encoder.encode(&values_ok, &modulus);

        // Should panic: 5 slots > max 4
        let values_too_many = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        std::panic::catch_unwind(|| {
            encoder.encode(&values_too_many, &modulus);
        })
        .expect_err("Should panic with too many values");
    }
}
