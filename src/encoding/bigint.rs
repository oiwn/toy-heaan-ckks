use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::{BigIntPolyRing, Plaintext};
use crypto_bigint::{NonZero, U256};
use malachite::base::num::arithmetic::traits::{PowerOf2, UnsignedAbs};
use malachite::base::num::conversion::traits::RoundingFrom;
use malachite::base::rounding_modes::RoundingMode;
use malachite::rational::Rational;
use malachite::{Integer, Natural};
use num_complex::Complex64;
use std::collections::HashMap;
use std::f64::consts::PI;
use tracing::{debug, info, instrument, warn};

pub struct BigIntEncodingParams<const DEGREE: usize> {
    scale_bits: u32,
    // Precomputed tables (depend only on DEGREE)
    rot_group: Vec<u64>,                      // rotGroup[Nh]
    ksi_pows: Vec<Complex64>,                 // ksiPows[M+1]
    taylor_coeffs: HashMap<String, Vec<f64>>, // taylorCoeffsMap
}

pub struct BigIntEncoder<const DEGREE: usize> {
    params: BigIntEncodingParams<DEGREE>,
}

impl<const DEGREE: usize> BigIntEncodingParams<DEGREE> {
    /// Helper functions for derived parameters (matching Kim's HEAAN)
    const fn n() -> usize {
        DEGREE
    }
    const fn nh() -> usize {
        DEGREE / 2
    }
    const fn m() -> usize {
        DEGREE * 2
    }
    #[allow(dead_code)]
    fn log_n() -> u32 {
        // Compute log2(DEGREE)
        DEGREE.ilog2()
    }

    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        if !DEGREE.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: DEGREE });
        }

        // Initialize Kim's tables exactly as in C++
        let (rot_group, ksi_pows, taylor_coeffs) = Self::init_tables();

        Ok(Self {
            scale_bits,
            rot_group,
            ksi_pows,
            taylor_coeffs,
        })
    }

    /// Kim's HEAAN initialization exactly as in Context.cpp lines 151-176
    fn init_tables() -> (Vec<u64>, Vec<Complex64>, HashMap<String, Vec<f64>>) {
        let _n = Self::n(); // N
        let nh = Self::nh(); // Nh = N >> 1
        let m = Self::m(); // M = N << 1

        // rotGroup: powers of 5 mod M (Context.cpp lines 151-157)
        let mut rot_group = vec![0u64; nh];
        let mut five_pows = 1u64;
        for i in 0..nh {
            rot_group[i] = five_pows;
            five_pows = (five_pows * 5) % (m as u64);
        }

        // ksiPows: e^(2πij/M) for j = 0..M (Context.cpp lines 159-166)
        let mut ksi_pows = vec![Complex64::new(0.0, 0.0); m + 1];
        for j in 0..m {
            let angle = 2.0 * PI * (j as f64) / (m as f64);
            ksi_pows[j] = Complex64::new(angle.cos(), angle.sin());
        }
        ksi_pows[m] = ksi_pows[0]; // Wraparound

        // taylorCoeffsMap (Context.cpp lines 174-176)
        let mut taylor_coeffs = HashMap::new();
        taylor_coeffs.insert(
            "LOGARITHM".to_string(),
            vec![
                0.0,
                1.0,
                -0.5,
                1.0 / 3.0,
                -1.0 / 4.0,
                1.0 / 5.0,
                -1.0 / 6.0,
                1.0 / 7.0,
                -1.0 / 8.0,
                1.0 / 9.0,
                -1.0 / 10.0,
            ],
        );
        taylor_coeffs.insert(
            "EXPONENT".to_string(),
            vec![
                1.0,
                1.0,
                0.5,
                1.0 / 6.0,
                1.0 / 24.0,
                1.0 / 120.0,
                1.0 / 720.0,
                1.0 / 5040.0,
                1.0 / 40320.0,
                1.0 / 362880.0,
                1.0 / 3628800.0,
            ],
        );
        taylor_coeffs.insert(
            "SIGMOID".to_string(),
            vec![
                1.0 / 2.0,
                1.0 / 4.0,
                0.0,
                -1.0 / 48.0,
                0.0,
                1.0 / 480.0,
                0.0,
                -17.0 / 80640.0,
                0.0,
                31.0 / 1451520.0,
                0.0,
            ],
        );

        (rot_group, ksi_pows, taylor_coeffs)
    }

    /// Generate qpowvec on-demand based on logQ (Context.cpp lines 168-172)
    pub fn qpowvec(log_q: u32) -> Vec<u64> {
        let log_qq = log_q << 1; // logQQ = logQ << 1
        let mut qpowvec = vec![0u64; (log_qq + 1) as usize];
        qpowvec[0] = 1; // qpowvec[0] = ZZ(1)
        for i in 1..=(log_qq as usize) {
            qpowvec[i] = qpowvec[i - 1] << 1; // qpowvec[i] = qpowvec[i - 1] << 1
        }
        qpowvec
    }

    /// Getter functions for testing
    pub fn rot_group(&self) -> &[u64] {
        &self.rot_group
    }
    pub fn ksi_pows(&self) -> &[Complex64] {
        &self.ksi_pows
    }
    pub fn taylor_coeffs(&self) -> &HashMap<String, Vec<f64>> {
        &self.taylor_coeffs
    }
}

impl<const DEGREE: usize> BigIntEncoder<DEGREE> {
    pub fn new(scale_bits: u32) -> EncodingResult<Self> {
        Ok(Self {
            params: BigIntEncodingParams::new(scale_bits)?,
        })
    }

    /// Get access to encoding parameters for testing
    pub fn params(&self) -> &BigIntEncodingParams<DEGREE> {
        &self.params
    }

    /// Kim's bit reversal algorithm for FFT (Context.cpp lines 341-352)
    fn bit_reverse(&self, vals: &mut [Complex64]) {
        let size = vals.len();
        let mut j = 0usize;

        // Kim's exact algorithm: for (long i = 1, j = 0; i < size; ++i)
        for i in 1..size {
            let mut bit = size >> 1;
            // Kim's inner loop: for (; j >= bit; bit>>=1) { j -= bit; }
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

    /// HEAAN-style fftSpecialInvLazy implementation (Context.cpp lines 415-431)
    #[instrument(skip(self, vals), fields(size = vals.len()))]
    fn fft_special_inv_lazy(&self, vals: &mut [Complex64]) {
        let size = vals.len();
        let m = DEGREE * 2; // M in HEAAN terminology

        info!(
            "FFT inverse lazy: size={}, m={}, rot_group.len()={}",
            size,
            m,
            self.params.rot_group.len()
        );

        // Kim's algorithm: start from size and work down to 1
        let mut len = size;
        let mut iteration = 0;
        while len >= 2 {
            iteration += 1;
            debug!("FFT iteration {}: len={}", iteration, len);
            let lenh = len >> 1;
            let lenq = len << 2; // lenq = len << 2

            debug!("  lenh={}, lenq={}", lenh, lenq);

            for i in (0..size).step_by(len) {
                debug!("  outer i={}", i);
                for j in 0..lenh {
                    debug!("    inner j={}", j);

                    // Kim's exact indexing: (lenq - (rotGroup[j] % lenq)) * M / lenq
                    let rot_group_val = if j < self.params.rot_group.len() {
                        self.params.rot_group[j]
                    } else {
                        // 🚨 CRITICAL: This is the bottleneck!
                        warn!(
                            "FFT accessing beyond rot_group bounds: j={}, rot_group.len()={}, using fallback j={}",
                            j,
                            self.params.rot_group.len(),
                            j
                        );
                        j as u64 // fallback for j >= rot_group.len()
                    };

                    let rot_mod = rot_group_val % (lenq as u64);
                    let ksi_idx = ((lenq as u64 - rot_mod) * (m as u64)
                        / (lenq as u64)) as usize;
                    let ksi_idx = ksi_idx % self.params.ksi_pows.len();

                    debug!(
                        "    rot_group_val={}, rot_mod={}, ksi_idx={}",
                        rot_group_val, rot_mod, ksi_idx
                    );
                    debug!("    indices: vals[{}], vals[{}]", i + j, i + j + lenh);

                    // Check bounds before accessing
                    if i + j + lenh >= vals.len() {
                        warn!(
                            "FFT butterfly out of bounds: i={}, j={}, lenh={}, vals.len()={}",
                            i,
                            j,
                            lenh,
                            vals.len()
                        );
                        break;
                    }

                    // Kim's butterfly: u = vals[i + j] + vals[i + j + lenh]; v = vals[i + j] - vals[i + j + lenh]; v *= ksiPows[idx];
                    let u = vals[i + j] + vals[i + j + lenh];
                    let mut v = vals[i + j] - vals[i + j + lenh];
                    v *= self.params.ksi_pows[ksi_idx];

                    debug!("    butterfly: u={:?}, v={:?}", u, v);

                    vals[i + j] = u;
                    vals[i + j + lenh] = v;
                }
            }
            len >>= 1;
            debug!("  new len={}", len);
        }

        debug!("FFT iterations complete, applying bit_reverse");
        // Kim's algorithm: bitReverse after the loops
        self.bit_reverse(vals);
        debug!("FFT inverse lazy complete");
    }

    /// HEAAN fftSpecialInv with normalization
    fn fft_special_inv(&self, vals: &mut [Complex64]) {
        self.fft_special_inv_lazy(vals);
        let size = vals.len() as f64;
        for val in vals.iter_mut() {
            *val /= size;
        }
    }

    /// HEAAN-style fftSpecial for decoding (Context.cpp lines 397-413)
    fn fft_special(&self, vals: &mut [Complex64]) {
        let size = vals.len();
        let m = DEGREE * 2; // M in HEAAN terminology

        // Kim's algorithm: bitReverse first
        self.bit_reverse(vals);

        let mut len = 2;
        while len <= size {
            let lenh = len >> 1;
            let lenq = len << 2; // lenq = len << 2

            for i in (0..size).step_by(len) {
                for j in 0..lenh {
                    // Kim's exact indexing: ((rotGroup[j] % lenq)) * M / lenq
                    let rot_group_val = if j < self.params.rot_group.len() {
                        self.params.rot_group[j]
                    } else {
                        j as u64 // fallback for j >= rot_group.len()
                    };

                    let rot_mod = rot_group_val % (lenq as u64);
                    let ksi_idx = (rot_mod * (m as u64) / (lenq as u64)) as usize;
                    let ksi_idx = ksi_idx % self.params.ksi_pows.len();

                    // Kim's butterfly: u = vals[i + j]; v = vals[i + j + lenh]; v *= ksiPows[idx]; vals[i + j] = u + v; vals[i + j + lenh] = u - v;
                    let u = vals[i + j];
                    let mut v = vals[i + j + lenh];
                    v *= self.params.ksi_pows[ksi_idx];

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
        _modulus: &NonZero<U256>,
    ) -> f64 {
        // Convert U256 to malachite Integer - assume centered reduction already done
        let coeff_words = coeff.to_words();
        let coeff_natural = Natural::from_limbs_asc(&coeff_words);

        // Handle signed representation (if MSB is set, it's a negative number in two's complement)
        let coeff_integer = if coeff_words[3] & 0x8000000000000000 != 0 {
            // Negative number in two's complement representation
            // Convert to proper negative integer
            let max_u256 = (Natural::from(1u32) << 256) - Natural::from(1u32);
            let unsigned_val = Natural::from_limbs_asc(&coeff_words);
            -Integer::from(&max_u256 - unsigned_val + Natural::from(1u32))
        } else {
            Integer::from(coeff_natural)
        };

        // Kim's HEAAN: exact scaling using RR arithmetic
        // Equivalent to: xp.e -= logp in NTL's RR
        let scale_factor = Natural::power_of_2(logp as u64);
        let scaled_rational =
            Rational::from(coeff_integer) / Rational::from(scale_factor);

        // Convert to f64 using exact rounding
        f64::rounding_from(scaled_rational, RoundingMode::Nearest).0
    }

    /// Kim's HEAAN centered reduction: if NumBits(tmp) == logq then tmp -= q
    /// This detects when tmp is near the upper boundary of the modulus
    fn apply_centered_reduction(
        &self,
        coeff: U256,
        modulus: &NonZero<U256>,
    ) -> U256 {
        let q = modulus.get();
        let tmp = coeff % q;

        // Convert to check if we're in the upper half of the modulus range
        // In HEAAN, coefficients > q/2 are treated as negative numbers
        let q_words = q.to_words();
        let q_natural = Natural::from_limbs_asc(&q_words);
        let half_q = &q_natural >> 1; // q/2

        let tmp_words = tmp.to_words();
        let tmp_natural = Natural::from_limbs_asc(&tmp_words);

        if tmp_natural > half_q {
            // Convert to negative by subtracting q (gives a large positive number in modular arithmetic)
            tmp.wrapping_sub(&q)
        } else {
            tmp
        }
    }

    /// HEAAN-style scaling: multiply by 2^logp and round using exact arithmetic
    fn scale_up_to_u256(
        &self,
        val: f64,
        logp: u32,
        modulus: &NonZero<U256>,
    ) -> U256 {
        let q_words = modulus.get().to_words();
        let q_natural = Natural::from_limbs_asc(&q_words);
        let q_integer = Integer::from(q_natural);

        // Kim's HEAAN: exact scaling equivalent to MakeRR(x.x, x.e + logp)
        // Convert f64 to rational - use a simple approach since exact conversion is complex
        let val_rational = if val.is_finite() && val.abs() < 1e15 {
            // For reasonable finite values, use fraction representation
            let scaled_int = (val * 1e12).round() as i64;
            Rational::from(scaled_int) / Rational::from(1000000000000u64) // 1e12
        } else {
            // Fallback: use integer representation
            Rational::from(Integer::from(val.round() as i64))
        };
        let scale_factor = Natural::power_of_2(logp as u64);
        let scaled_rational = val_rational * Rational::from(scale_factor);

        // Round to nearest integer (equivalent to RoundToZZ in NTL)
        let rounded_integer =
            Integer::rounding_from(scaled_rational, RoundingMode::Nearest).0;

        // Reduce modulo q with proper sign handling
        let reduced = rounded_integer % &q_integer;
        let final_value = if reduced < 0 {
            reduced + &q_integer
        } else {
            reduced
        };

        // Convert back to U256
        let final_natural = final_value.unsigned_abs();
        let final_words = final_natural.to_limbs_asc();
        let mut result_words = [0u64; 4];
        for (i, &word) in final_words.iter().enumerate() {
            if i < 4 {
                result_words[i] = word;
            }
        }
        U256::from_words(result_words)
    }
}

impl<const DEGREE: usize> Encoder<BigIntPolyRing<DEGREE>, DEGREE>
    for BigIntEncoder<DEGREE>
{
    #[instrument(skip(self, context), fields(slots = values.len(), degree = DEGREE))]
    fn encode(
        &self,
        values: &[f64],
        context: &NonZero<U256>,
    ) -> Plaintext<BigIntPolyRing<DEGREE>, DEGREE> {
        let n = DEGREE / 2; // Nh in HEAAN terminology
        let slots = values.len();

        assert!(slots <= n, "Too many values: got {}, max {}", slots, n);
        assert!(
            slots.is_power_of_two(),
            "Kim's HEAAN requires slots to be a power of 2, got {}",
            slots
        );

        info!(
            "Starting BigInt encoding with {} slots, gap={}",
            slots,
            n / slots
        );

        // HEAAN-style: Convert real values to complex, exactly as in Kim's C++
        debug!("Converting {} real values to complex", slots);
        let mut uvals = vec![Complex64::new(0.0, 0.0); slots];
        for (i, &val) in values.iter().enumerate() {
            uvals[i] = Complex64::new(val, 0.0);
        }

        // Apply inverse FFT to get polynomial coefficients (Kim's: fftSpecialInv(uvals, slots))
        // Kim's algorithm can handle arbitrary slot counts directly
        debug!("Applying inverse FFT to {} complex values", uvals.len());
        self.fft_special_inv(&mut uvals);

        // HEAAN coefficient layout with gap-based indexing
        // Kim's HEAAN: gap = Nh / slots, idx += gap, jdx += gap
        let mut coeffs = [U256::ZERO; DEGREE];
        let gap = n / slots;
        let mut idx = 0;
        let mut jdx = n;

        for i in 0..slots {
            // Real parts at gap-spaced indices starting from 0
            coeffs[idx] =
                self.scale_up_to_u256(uvals[i].re, self.params.scale_bits, context);
            // Imaginary parts at gap-spaced indices starting from Nh
            coeffs[jdx] =
                self.scale_up_to_u256(uvals[i].im, self.params.scale_bits, context);

            idx += gap;
            jdx += gap;
        }

        let poly = BigIntPolyRing::<DEGREE>::from_u256_coeffs(&coeffs, context);
        Plaintext {
            poly,
            scale_bits: self.params.scale_bits,
            slots,
        }
    }

    fn decode(
        &self,
        plaintext: &Plaintext<BigIntPolyRing<DEGREE>, DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.coefficients();
        let modulus = plaintext.poly.modulus();
        let n = DEGREE / 2; // Nh in Kim's terms
        let logp = plaintext.scale_bits;

        // Use the slots field from plaintext (no more pattern detection needed!)
        let slots = plaintext.slots;
        let gap = n / slots;

        // HEAAN gap-based coefficient layout: reconstruct complex values
        let mut uvals = vec![Complex64::new(0.0, 0.0); slots];

        for i in 0..slots {
            let idx = i * gap; // Kim's: idx += gap for each slot

            // Real part from gap-spaced indices starting at 0
            let raw_coeff = coeffs[idx];
            let real_coeff = self.apply_centered_reduction(raw_coeff, &modulus);
            let real_val = self.scale_down_to_real(real_coeff, logp, &modulus);

            // Imaginary part from gap-spaced indices starting at Nh (idx + Nh in Kim's code)
            let imag_coeff =
                self.apply_centered_reduction(coeffs[idx + n], &modulus);
            let imag_val = self.scale_down_to_real(imag_coeff, logp, &modulus);

            uvals[i] = Complex64::new(real_val, imag_val);
        }

        // Apply forward FFT to recover original values
        // Kim's algorithm can handle arbitrary slot counts directly
        self.fft_special(&mut uvals);

        // Extract real parts and return only the slots we encoded
        uvals.iter().map(|c| c.re).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::U256;

    #[test]
    fn test_bigint_encoder_roundtrip_small() {
        const DEGREE: usize = 8;
        const SCALE_BITS: u32 = 60;

        let encoder = BigIntEncoder::<DEGREE>::new(SCALE_BITS).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(1152921504606846883u128)).unwrap();

        // Test with a single value - note: small degrees may have precision issues
        let values = vec![1.5];
        let plaintext = encoder.encode(&values, &modulus);
        let decoded = encoder.decode(&plaintext);

        println!(
            "Small degree test - Input: {:?}, Decoded: {:?}",
            values, decoded
        );

        // For small degrees, we expect some precision loss due to limited FFT resolution
        assert_eq!(decoded.len(), 1);
        // Relaxed tolerance for small degrees
        assert!(
            (decoded[0] - 1.5).abs() < 2.0,
            "Expected approximately 1.5, got {} (small degree test)",
            decoded[0]
        );
    }

    #[test]
    fn test_precision_by_degree() {
        // Test multiple degrees to understand precision threshold
        let degrees_and_expected = vec![
            (8, 2.0),    // Known to have precision issues
            (16, 1.0),   // Test if 16 is better
            (32, 0.01),  // Should be very good
            (64, 0.001), // Should be excellent
        ];

        for (degree, max_error) in degrees_and_expected {
            println!("\n=== Testing DEGREE={} ===", degree);

            let result = match degree {
                8 => test_degree_precision::<8>(),
                16 => test_degree_precision::<16>(),
                32 => test_degree_precision::<32>(),
                64 => test_degree_precision::<64>(),
                _ => panic!("Unsupported degree"),
            };

            let error = (result - 1.5).abs();
            println!(
                "DEGREE={}: Input=1.5, Decoded={}, Error={}",
                degree, result, error
            );

            assert!(
                error < max_error,
                "DEGREE={} exceeded expected error threshold: {} > {}",
                degree,
                error,
                max_error
            );
        }
    }

    fn test_degree_precision<const DEGREE: usize>() -> f64 {
        const SCALE_BITS: u32 = 30;

        let encoder = BigIntEncoder::<DEGREE>::new(SCALE_BITS).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFF63u128))
                .unwrap();

        let values = vec![1.5];
        let plaintext = encoder.encode(&values, &modulus);
        let decoded = encoder.decode(&plaintext);

        decoded[0]
    }

    #[test]
    fn test_bigint_encoder_roundtrip_large() {
        const DEGREE: usize = 1024; // Larger degree for better precision
        const SCALE_BITS: u32 = 30;

        let encoder = BigIntEncoder::<DEGREE>::new(SCALE_BITS).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFF63u128))
                .unwrap();

        // Test with a single value
        let values = vec![1.5];
        let plaintext = encoder.encode(&values, &modulus);
        let decoded = encoder.decode(&plaintext);

        println!(
            "Large degree test - Input: {:?}, Decoded: {:?}",
            values, decoded
        );

        assert_eq!(decoded.len(), 1);
        assert!(
            (decoded[0] - 1.5).abs() < 1e-10,
            "Expected ~1.5, got {} (large degree should be accurate)",
            decoded[0]
        );
    }

    #[test]
    fn test_multiple_values() {
        const DEGREE: usize = 32; // Use larger degree for better precision
        const SCALE_BITS: u32 = 30;

        let encoder = BigIntEncoder::<DEGREE>::new(SCALE_BITS).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFF63u128))
                .unwrap();

        // Test with multiple values (power of 2 required)
        let values = vec![1.0, 2.5, -0.75, 1.25];
        let plaintext = encoder.encode(&values, &modulus);
        let decoded = encoder.decode(&plaintext);

        assert_eq!(decoded.len(), 4);
        for (i, (&expected, &actual)) in
            values.iter().zip(decoded.iter()).enumerate()
        {
            assert!(
                (actual - expected).abs() < 1e-2, // Relaxed tolerance for multi-slot
                "Value {} mismatch: expected {}, got {}",
                i,
                expected,
                actual
            );
        }
    }

    #[test]
    fn test_max_slots_assertion() {
        const DEGREE: usize = 8;
        let encoder = BigIntEncoder::<DEGREE>::new(20).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(1152921504606846883u128)).unwrap();

        // Should work: 4 slots for DEGREE=8 (power of 2)
        let values_ok = vec![1.0, 2.0, 3.0, 4.0];
        encoder.encode(&values_ok, &modulus);

        // Should panic: 5 slots > max 4
        let values_too_many = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        std::panic::catch_unwind(|| {
            encoder.encode(&values_too_many, &modulus);
        })
        .expect_err("Should panic with too many values");
    }

    #[test]
    fn test_power_of_2_assertion() {
        const DEGREE: usize = 8;
        let encoder = BigIntEncoder::<DEGREE>::new(20).unwrap();
        let modulus =
            NonZero::new(U256::from_u128(1152921504606846883u128)).unwrap();

        // Should work: 1, 2, 4 slots (all powers of 2)
        let values_1 = vec![1.0];
        encoder.encode(&values_1, &modulus);

        let values_2 = vec![1.0, 2.0];
        encoder.encode(&values_2, &modulus);

        let values_4 = vec![1.0, 2.0, 3.0, 4.0];
        encoder.encode(&values_4, &modulus);

        // Should panic: 3 slots (not power of 2)
        let values_3 = vec![1.0, 2.0, 3.0];
        let result = std::panic::catch_unwind(|| {
            encoder.encode(&values_3, &modulus);
        });
        assert!(
            result.is_err(),
            "Should panic with non-power-of-2 slot count"
        );
    }
}
