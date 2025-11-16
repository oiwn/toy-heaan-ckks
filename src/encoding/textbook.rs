//! CKKS encoder/decoder implementation based on the FHE textbook
//! Vandermonde (special FFT) embedding. This module currently serves the RNS
//! backend exclusively.

use std::{cmp::Ordering, sync::Arc};

use crypto_bigint::{BoxedUint, Limb, NonZero};
use rustfft::num_complex::Complex64;

use crate::encoding::{Encoder, EncodingError, EncodingResult};
use crate::rings::backends::rns::{RnsBasis, RnsNttPoly, mod_inverse};
use crate::{Plaintext, PolyRing};

use super::special_fft::{
    VandermondeTables, build_slots_from_complex, special_dft, special_idft,
};

#[derive(Debug, Clone)]
pub struct CrtPrecomputation {
    modulus: BoxedUint,
    half_modulus: BoxedUint,
    prefix_products: Vec<BoxedUint>,
    prefix_inverses: Vec<u64>,
}

impl CrtPrecomputation {
    fn new(primes: &[u64]) -> EncodingResult<Self> {
        if primes.is_empty() {
            return Err(EncodingError::InvalidInput {
                message: "prime list cannot be empty".into(),
            });
        }

        let mut prefix_products = Vec::with_capacity(primes.len());
        let mut prefix_inverses = vec![0u64; primes.len()];

        let mut current = BoxedUint::from(primes[0]);
        prefix_products.push(current.clone());

        for (idx, &prime) in primes.iter().enumerate().skip(1) {
            let prime_nz =
                NonZero::new(Limb::from(prime)).expect("prime must be non-zero");
            let residue = current.rem_limb(prime_nz);
            let inv = mod_inverse(residue.0, prime);
            prefix_inverses[idx] = inv;
            current = current.mul(&BoxedUint::from(prime));
            prefix_products.push(current.clone());
        }

        let modulus = prefix_products
            .last()
            .cloned()
            .unwrap_or_else(|| BoxedUint::from(1u64));
        let half_modulus = modulus.shr_vartime(1).unwrap_or_else(|| {
            BoxedUint::zero_with_precision(modulus.bits_precision())
        });

        Ok(Self {
            modulus,
            half_modulus,
            prefix_products,
            prefix_inverses,
        })
    }
}

/// Shared parameter object for the textbook encoder.
#[derive(Debug, Clone)]
pub struct TextbookEncodingParams<const N: usize> {
    pub degree: usize,
    pub slots: usize,
    pub scale_bits: u32,
    pub primes: Arc<Vec<u64>>,
    pub vandermonde: Arc<VandermondeTables<N>>,
    pub crt: Arc<CrtPrecomputation>,
}

impl<const N: usize> TextbookEncodingParams<N> {
    pub fn new(scale_bits: u32, primes: Arc<Vec<u64>>) -> EncodingResult<Self> {
        if !N.is_power_of_two() {
            return Err(EncodingError::InvalidRingDegree { degree: N });
        }

        if scale_bits == 0 {
            return Err(EncodingError::InvalidInput {
                message: "scale_bits must be positive".into(),
            });
        }

        if primes.is_empty() {
            return Err(EncodingError::InvalidInput {
                message: "prime list cannot be empty".into(),
            });
        }

        let order = (2 * N) as u64;
        for &prime in primes.iter() {
            let modulus_minus_one = prime.checked_sub(1).ok_or_else(|| {
                EncodingError::InvalidInput {
                    message: format!("prime {prime} too small"),
                }
            })?;
            if modulus_minus_one % order != 0 {
                return Err(EncodingError::InvalidInput {
                    message: format!(
                        "prime {prime} is not NTT-friendly for degree {N}"
                    ),
                });
            }
        }

        let crt = Arc::new(CrtPrecomputation::new(primes.as_ref())?);

        Ok(Self {
            degree: N,
            slots: N / 2,
            scale_bits,
            primes,
            vandermonde: Arc::new(VandermondeTables::new()),
            crt,
        })
    }

    pub fn scaling_factor(&self) -> f64 {
        2f64.powi(self.scale_bits as i32)
    }

    pub fn max_slots(&self) -> usize {
        self.slots
    }

    pub fn validate_slot_count(&self, requested: usize) -> EncodingResult<()> {
        if requested > self.slots {
            Err(EncodingError::InputTooLong {
                got: requested,
                max: self.slots,
            })
        } else {
            Ok(())
        }
    }
}

/// Textbook encoder that targets `RnsNttPoly` plaintexts.
pub struct TextbookEncoder<const N: usize> {
    pub params: TextbookEncodingParams<N>,
}

impl<const N: usize> TextbookEncoder<N> {
    pub fn new(params: TextbookEncodingParams<N>) -> Self {
        Self { params }
    }

    pub fn encode_complex(
        &self,
        values: &[Complex64],
        basis: &Arc<RnsBasis>,
    ) -> EncodingResult<Plaintext<RnsNttPoly<N>, N>> {
        self.params.validate_slot_count(values.len())?;
        self.ensure_basis(basis)?;

        let scale = self.params.scaling_factor();
        let scaled: Vec<Complex64> = values
            .iter()
            .map(|value| Complex64::new(value.re * scale, value.im * scale))
            .collect();
        let slots = build_slots_from_complex::<N>(&scaled);
        let coeffs = special_idft::<N>(&slots, &self.params.vandermonde);
        let rounded = self.round_coefficients(&coeffs)?;
        let mut poly = RnsNttPoly::from_i64_slice(&rounded, basis.clone());
        poly.to_ntt_domain();
        Ok(Plaintext {
            poly,
            scale_bits: self.params.scale_bits,
            slots: values.len(),
        })
    }

    pub fn decode_complex(
        &self,
        plaintext: &Plaintext<RnsNttPoly<N>, N>,
    ) -> EncodingResult<Vec<Complex64>> {
        let mut coeff_poly = plaintext.poly.clone();
        coeff_poly.to_coeff_domain();

        let channel_count = coeff_poly.channels();
        let mut residues = vec![0u64; channel_count];
        let mut coeffs = vec![Complex64::new(0.0, 0.0); N];

        for coeff_idx in 0..N {
            for (channel_idx, channel) in coeff_poly.coefficients.iter().enumerate()
            {
                residues[channel_idx] = channel[coeff_idx];
            }
            let real = self.reconstruct_real_value(&residues)?;
            coeffs[coeff_idx] = Complex64::new(real, 0.0);
        }

        let slots = special_dft::<N>(&coeffs, &self.params.vandermonde);
        let scale = self.params.scaling_factor();
        Ok(slots
            .into_iter()
            .take(plaintext.slots)
            .map(|slot| slot / scale)
            .collect())
    }
}

impl<const N: usize> TextbookEncoder<N> {
    fn ensure_basis(&self, basis: &Arc<RnsBasis>) -> EncodingResult<()> {
        if basis.primes() != self.params.primes.as_ref() {
            return Err(EncodingError::InvalidInput {
                message: "basis primes do not match encoding params".into(),
            });
        }
        Ok(())
    }

    fn round_coefficients(&self, coeffs: &[Complex64]) -> EncodingResult<[i64; N]> {
        const IMAG_EPS: f64 = 1e-6;
        let mut rounded = [0i64; N];
        for (idx, value) in coeffs.iter().enumerate() {
            if value.im.abs() > IMAG_EPS {
                println!(
                    "round_coefficients: coeff {idx} = {} + {}i",
                    value.re, value.im
                );
                return Err(EncodingError::InvalidInput {
                    message: format!(
                        "coefficient {idx} has imaginary part {} exceeding tolerance",
                        value.im
                    ),
                });
            }
            let real = value.re;
            if real.abs() > i64::MAX as f64 {
                return Err(EncodingError::CoefficientOutOfRange { value: real });
            }
            rounded[idx] = real.round() as i64;
        }
        Ok(rounded)
    }

    fn reconstruct_real_value(&self, residues: &[u64]) -> EncodingResult<f64> {
        let primes = self.params.primes.as_ref();
        if residues.len() != primes.len() {
            return Err(EncodingError::InvalidInput {
                message: format!(
                    "residue count {} does not match prime count {}",
                    residues.len(),
                    primes.len()
                ),
            });
        }

        let crt = self.params.crt.as_ref();
        let mut value = BoxedUint::from(residues[0]);

        for (idx, &prime) in primes.iter().enumerate().skip(1) {
            let prime_nz =
                NonZero::new(Limb::from(prime)).expect("prime must be non-zero");
            let current_mod = value.rem_limb(prime_nz);
            let delta = mod_sub(residues[idx], current_mod.0, prime);
            if delta == 0 {
                continue;
            }
            let inv = crt.prefix_inverses[idx];
            let adjustment = mod_mul(delta, inv, prime);
            if adjustment == 0 {
                continue;
            }
            let term =
                crt.prefix_products[idx - 1].mul(&BoxedUint::from(adjustment));
            value += term;
        }

        match value.cmp(&crt.half_modulus) {
            Ordering::Greater => {
                let mut diff = crt.modulus.clone();
                diff -= &value;
                Ok(-boxed_to_f64(&diff))
            }
            _ => Ok(boxed_to_f64(&value)),
        }
    }
}

impl<const N: usize> Encoder<RnsNttPoly<N>, N> for TextbookEncoder<N> {
    fn encode(
        &self,
        values: &[f64],
        context: &<RnsNttPoly<N> as PolyRing<N>>::Context,
    ) -> Plaintext<RnsNttPoly<N>, N> {
        let complex: Vec<Complex64> =
            values.iter().map(|&v| Complex64::new(v, 0.0)).collect();
        self.encode_complex(&complex, context)
            .expect("textbook encode not ready")
    }

    fn decode(&self, plaintext: &Plaintext<RnsNttPoly<N>, N>) -> Vec<f64> {
        let decoded = self
            .decode_complex(plaintext)
            .expect("textbook decode not ready");
        decoded.iter().map(|c| c.re).collect()
    }
}

fn mod_sub(value: u64, subtrahend: u64, modulus: u64) -> u64 {
    if value >= subtrahend {
        (value - subtrahend) % modulus
    } else {
        (value + modulus - subtrahend) % modulus
    }
}

fn mod_mul(lhs: u64, rhs: u64, modulus: u64) -> u64 {
    ((lhs as u128 * rhs as u128) % modulus as u128) as u64
}

fn boxed_to_f64(value: &BoxedUint) -> f64 {
    const LIMB_SCALE: f64 = 18446744073709551616.0; // 2^64
    let mut acc = 0.0;
    let mut factor = 1.0;
    for word in value.as_words() {
        acc += (*word as f64) * factor;
        factor *= LIMB_SCALE;
    }
    acc
}
