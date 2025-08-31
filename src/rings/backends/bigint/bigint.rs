use crate::{
    Ciphertext, PolyRescale, PolyRing, PolySampler, math::ternary_coefficients,
};
use crypto_bigint::{NonZero, U256, Zero};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::ops::{AddAssign, MulAssign, Neg};

pub use super::display::CompactPolyDisplay;

#[derive(Debug, Clone, PartialEq)]
pub struct BigIntPolyRing<const DEGREE: usize> {
    pub coeffs: [U256; DEGREE],
    modulus: NonZero<U256>,
}

impl<const DEGREE: usize> BigIntPolyRing<DEGREE> {
    pub fn with_modulus(modulus: NonZero<U256>) -> Self {
        Self {
            coeffs: [U256::ZERO; DEGREE],
            modulus,
        }
    }

    pub fn coefficients(&self) -> [U256; DEGREE] {
        self.coeffs
    }

    pub fn modulus(&self) -> NonZero<U256> {
        self.modulus
    }

    /// Create polynomial from U256 coefficients directly
    pub fn from_u256_coeffs(coeffs: &[U256], context: &NonZero<U256>) -> Self {
        let mut result_coeffs = [U256::ZERO; DEGREE];
        for (i, &coeff) in coeffs.iter().enumerate().take(DEGREE) {
            result_coeffs[i] = coeff % context.get(); // Ensure coefficient is in range
        }

        Self {
            coeffs: result_coeffs,
            modulus: *context,
        }
    }
}

impl<const DEGREE: usize> PolyRing<DEGREE> for BigIntPolyRing<DEGREE> {
    type Context = crypto_bigint::NonZero<U256>;

    fn zero(context: &Self::Context) -> Self {
        Self {
            coeffs: [U256::ZERO; DEGREE],
            modulus: *context,
        }
    }

    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        assert!(
            coeffs.len() >= DEGREE,
            "Insufficient number of coefficients"
        );

        let mut poly_coeffs = [U256::ZERO; DEGREE];
        for (i, &coeff) in coeffs.iter().take(DEGREE).enumerate() {
            poly_coeffs[i] = if coeff >= 0 {
                U256::from(coeff as u64)
            } else {
                // For negative coefficients: modulus - |coeff|
                let abs_coeff = U256::from((-coeff) as u64);
                context.wrapping_sub(&abs_coeff)
            };
        }

        Self {
            coeffs: poly_coeffs,
            modulus: *context,
        }
    }

    fn context(&self) -> &Self::Context {
        &self.modulus
    }

    fn to_coeffs(&self) -> [i64; DEGREE] {
        let mut signed_coeffs = [0i64; DEGREE];
        let half = self.modulus.wrapping_shr(1); // modulus / 2

        for (i, &c) in self.coeffs.iter().enumerate() {
            signed_coeffs[i] = if c < half {
                // Convert to u64 first, then to i64 (assuming it fits)
                // For toy implementation, we assume coefficients fit in i64
                let as_u64 = c.as_words()[0]; // Get lowest 64 bits
                as_u64 as i64
            } else {
                // Negative representation: -(modulus - c)
                let diff = self.modulus.wrapping_sub(&c);
                let as_u64 = diff.as_words()[0];
                -(as_u64 as i64)
            };
        }

        signed_coeffs
    }
}

impl<const DEGREE: usize> AddAssign<&Self> for BigIntPolyRing<DEGREE> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.modulus, rhs.modulus,
            "Cannot add polynomials with different moduli"
        );

        for i in 0..DEGREE {
            // Use crypto-bigint modular addition
            self.coeffs[i] = self.coeffs[i].add_mod(&rhs.coeffs[i], &self.modulus);
        }
    }
}

impl<const DEGREE: usize> MulAssign<&Self> for BigIntPolyRing<DEGREE> {
    fn mul_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.modulus, rhs.modulus,
            "Cannot multiply polynomials with different moduli"
        );

        // Schoolbook polynomial multiplication in Z[X]/(X^DEGREE + 1)
        let mut result = [U256::ZERO; DEGREE];

        for i in 0..DEGREE {
            for j in 0..DEGREE {
                let coeff_product =
                    self.coeffs[i].mul_mod(&rhs.coeffs[j], &self.modulus);

                if i + j < DEGREE {
                    // Normal coefficient
                    result[i + j] =
                        result[i + j].add_mod(&coeff_product, &self.modulus);
                } else {
                    // Wrap around with negation due to X^DEGREE = -1
                    let wrapped_idx = (i + j) - DEGREE;
                    result[wrapped_idx] =
                        result[wrapped_idx].sub_mod(&coeff_product, &self.modulus);
                }
            }
        }

        self.coeffs = result;
    }
}

impl<const DEGREE: usize> Neg for BigIntPolyRing<DEGREE> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for coeff in &mut self.coeffs {
            *coeff = self.modulus.wrapping_sub(coeff);
        }
        self
    }
}

impl<const DEGREE: usize> PolySampler<DEGREE> for BigIntPolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self {
        // Generate uniform coefficients using full U256 range
        let mut coeffs = [U256::ZERO; DEGREE];

        for coeff in &mut coeffs {
            // Build a full U256 from four u64 values for proper uniform distribution
            let words = [
                rng.random::<u64>(),
                rng.random::<u64>(),
                rng.random::<u64>(),
                rng.random::<u64>(),
            ];
            let random_u256 = U256::from_words(words);
            *coeff = random_u256.rem(context);
        }

        Self {
            coeffs,
            modulus: *context,
        }
    }

    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        // Sample centered Gaussian integers, map into Z_q using full U256 arithmetic
        let q = context.get();
        let normal = Normal::new(0.0, std_dev).expect("invalid Gaussian std_dev");
        let mut coeffs = [U256::ZERO; DEGREE];

        for c in &mut coeffs {
            let sample = normal.sample(rng).round() as i128; // centered integer
            if sample == 0 {
                *c = U256::ZERO;
                continue;
            }
            if sample > 0 {
                // Reduce positive sample modulo q
                let v = U256::from_u128(sample as u128);
                *c = v.rem(context);
            } else {
                // Map negative sample to q - (|v| mod q)
                let v_abs = U256::from_u128((-sample) as u128).rem(context);
                *c = if v_abs.is_zero().into() {
                    U256::ZERO
                } else {
                    q.wrapping_sub(&v_abs)
                };
            }
        }

        Self {
            coeffs,
            modulus: *context,
        }
    }

    fn sample_tribits<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let ternary = ternary_coefficients::<DEGREE, R>(hamming_weight, rng);
        let mut coeffs = [U256::ZERO; DEGREE];

        for (i, &val) in ternary.iter().enumerate() {
            coeffs[i] = if val == -1 {
                context.wrapping_sub(&U256::ONE) // modulus - 1
            } else {
                U256::from(val as u64)
            };
        }

        Self {
            coeffs,
            modulus: *context,
        }
    }

    fn sample_noise<R: Rng>(
        variance: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        Self::sample_gaussian(variance.sqrt(), context, rng)
    }
}

impl<const DEGREE: usize> PolyRescale<DEGREE> for BigIntPolyRing<DEGREE> {
    fn rescale_assign(&mut self, scale_factor: f64) {
        // Extract scale_bits from scale_factor (assuming scale_factor = 2^k)
        let k = (scale_factor.log2().round()) as u32;
        let modulus = self.modulus.get();

        for coeff in &mut self.coeffs {
            *coeff = rescale_coeff_pow2_u256(*coeff, modulus, k);
        }
    }
}

/// Round-to-nearest division by 2^k in Z_q, sign-aware, mapping back to [0, q).
fn rescale_coeff_pow2_u256(c: U256, q: U256, k: u32) -> U256 {
    let half = q >> 1;
    // centered lift to Z: x in (-(q/2), q/2]
    let (neg, x_abs) = if c > half {
        (true, q.wrapping_sub(&c))
    } else {
        (false, c)
    };

    // y = round(x_abs / 2^k) = (x_abs + 2^(k-1)) >> k
    let round = U256::ONE << (k - 1);
    let y = x_abs.saturating_add(&round) >> k;

    // map back to [0, q)
    if neg {
        if y.is_zero().into() {
            U256::ZERO
        } else {
            q.wrapping_sub(&y)
        }
    } else {
        y
    }
}

pub fn rescale_ciphertext_u256_inplace<const DEGREE: usize>(
    ct: &mut Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>,
    k: u32,
) {
    let q = ct.c0.modulus().get();

    for i in 0..DEGREE {
        let c0 = ct.c0.coeffs[i];
        let c1 = ct.c1.coeffs[i];

        let r0 = rescale_coeff_pow2_u256(c0, q, k);
        let r1 = rescale_coeff_pow2_u256(c1, q, k);

        ct.c0.coeffs[i] = r0;
        ct.c1.coeffs[i] = r1;
    }

    ct.scale_bits -= k;
}
