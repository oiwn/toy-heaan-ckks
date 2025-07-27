use crate::{
    PolyRing, PolySampler,
    math::{gaussian_coefficients, ternary_coefficients},
};
use crypto_bigint::{NonZero, U256};
use rand::Rng;
use std::ops::{AddAssign, MulAssign, Neg};

#[derive(Debug, Clone, PartialEq)]
pub struct PolyRingU256<const DEGREE: usize> {
    pub coeffs: [U256; DEGREE],
    modulus: NonZero<U256>,
}

impl<const DEGREE: usize> PolyRingU256<DEGREE> {
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

impl<const DEGREE: usize> PolyRing<DEGREE> for PolyRingU256<DEGREE> {
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

impl<const DEGREE: usize> AddAssign<&Self> for PolyRingU256<DEGREE> {
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

impl<const DEGREE: usize> MulAssign<&Self> for PolyRingU256<DEGREE> {
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

impl<const DEGREE: usize> Neg for PolyRingU256<DEGREE> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for coeff in &mut self.coeffs {
            *coeff = self.modulus.wrapping_sub(coeff);
        }
        self
    }
}

impl<const DEGREE: usize> PolySampler<DEGREE> for PolyRingU256<DEGREE> {
    fn sample_uniform<R: Rng>(rng: &mut R, context: &Self::Context) -> Self {
        // Generate uniform coefficients as u64, then convert to U256
        let mut coeffs = [U256::ZERO; DEGREE];

        for coeff in &mut coeffs {
            // Sample a random u64 and reduce modulo context
            let random_u64 = rng.random::<u64>();
            let random_u256 = U256::from(random_u64);
            *coeff = random_u256.rem(context);
        }

        Self {
            coeffs,
            modulus: *context,
        }
    }

    fn sample_gaussian<R: Rng>(
        rng: &mut R,
        std_dev: f64,
        context: &Self::Context,
    ) -> Self {
        // Use existing gaussian sampling for u64, then convert
        let gaussian_u64 =
            gaussian_coefficients::<DEGREE, R>(std_dev, context.as_words()[0], rng);
        let mut coeffs = [U256::ZERO; DEGREE];

        for (i, &val) in gaussian_u64.iter().enumerate() {
            coeffs[i] = U256::from(val);
        }

        Self {
            coeffs,
            modulus: *context,
        }
    }

    fn sample_tribits<R: Rng>(
        rng: &mut R,
        hamming_weight: usize,
        context: &Self::Context,
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
        rng: &mut R,
        variance: f64,
        context: &Self::Context,
    ) -> Self {
        Self::sample_gaussian(rng, variance.sqrt(), context)
    }
}
