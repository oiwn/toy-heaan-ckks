use crate::{
    PolyRing, PolySampler,
    math::{gaussian_coefficients, ternary_coefficients, uniform_coefficients},
};
use rand::Rng;
use std::ops::{AddAssign, MulAssign, Neg};

#[derive(Debug, Clone, PartialEq)]
pub struct NaivePolyRing<const DEGREE: usize> {
    pub coeffs: [u64; DEGREE],
    pub context: u64,
}

impl<const DEGREE: usize> NaivePolyRing<DEGREE> {
    pub fn default_modulus() -> u64 {
        return 741507920154517877u64;
    }

    pub fn with_modulus(modulus: u64) -> Self {
        Self {
            coeffs: [0; DEGREE],
            context: modulus,
        }
    }
}

impl<const DEGREE: usize> PolyRing<DEGREE> for NaivePolyRing<DEGREE> {
    type Context = u64; // Just the modulus for naive implementation

    fn zero(context: &Self::Context) -> Self {
        Self {
            coeffs: [0; DEGREE],
            context: *context,
        }
    }

    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        assert!(
            coeffs.len() >= DEGREE,
            "Insufficient number of coefficients"
        );

        let mut poly_coeffs = [0u64; DEGREE];
        for (i, coeff) in coeffs.iter().enumerate() {
            poly_coeffs[i] = if *coeff >= 0 {
                *coeff as u64
            } else {
                (context - (-coeff as u64)) % context
            };
        }

        Self {
            coeffs: poly_coeffs,
            context: *context,
        }
    }

    fn context(&self) -> &Self::Context {
        &self.context
    }

    fn to_coeffs(&self) -> [i64; DEGREE] {
        let mut signed_coeffs = [0i64; DEGREE];
        let half = self.context / 2;

        for (i, &c) in self.coeffs.iter().enumerate() {
            signed_coeffs[i] = if c < half {
                c as i64
            } else {
                (c as i64).wrapping_sub(self.context as i64)
            };
        }

        signed_coeffs
    }
}

impl<const DEGREE: usize> AddAssign<&Self> for NaivePolyRing<DEGREE> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.context, rhs.context,
            "Cannot add polynomials with different moduli: {} vs {}",
            self.context, rhs.context
        );
        for i in 0..DEGREE {
            // self.coeffs[i] = (self.coeffs[i] + rhs.coeffs[i]) % self.modulus;
            let sum = (self.coeffs[i] as u128 + rhs.coeffs[i] as u128)
                % self.context as u128;
            self.coeffs[i] = sum.try_into().unwrap();
        }
    }
}

impl<const DEGREE: usize> MulAssign<&Self> for NaivePolyRing<DEGREE> {
    fn mul_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.context, rhs.context,
            "Cannot add polynomials with different moduli: {} vs {}",
            self.context, rhs.context
        );
        // Schoolbook multiplication with X^DEGREE + 1 reduction
        let mut result = [0u64; DEGREE];

        for i in 0..DEGREE {
            for j in 0..DEGREE {
                let coeff_pos = i + j;
                let coeff_val = (self.coeffs[i] as u128 * rhs.coeffs[j] as u128)
                    % self.context as u128;

                if coeff_pos < DEGREE {
                    result[coeff_pos] = ((result[coeff_pos] as u128 + coeff_val)
                        % self.context as u128)
                        as u64;
                } else {
                    // X^DEGREE = -1, so X^(DEGREE+k) = -X^k
                    let wrapped_pos = coeff_pos - DEGREE;
                    result[wrapped_pos] =
                        ((result[wrapped_pos] as u128 + self.context as u128
                            - coeff_val)
                            % self.context as u128) as u64;
                }
            }
        }

        self.coeffs = result.map(|x| x as u64);
    }
}

impl<const DEGREE: usize> Neg for NaivePolyRing<DEGREE> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for coeff in &mut self.coeffs {
            *coeff = (self.context - *coeff) % self.context;
        }
        self
    }
}

impl<const DEGREE: usize> PolySampler<DEGREE> for NaivePolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(rng: &mut R, context: &Self::Context) -> Self {
        Self {
            coeffs: uniform_coefficients::<DEGREE, R>(*context, rng), // Use context instead of max_coeff
            context: *context,
        }
    }

    fn sample_gaussian<R: Rng>(rng: &mut R, std_dev: f64) -> Self {
        // Need a default context - this is a problem we need to solve
        let context = u64::MAX; // Or get from somewhere else
        Self {
            coeffs: gaussian_coefficients::<DEGREE, R>(std_dev, context, rng),
            context,
        }
    }

    fn sample_tribits<R: Rng>(rng: &mut R, hamming_weight: usize) -> Self {
        let context = u64::MAX; // Or get from somewhere else
        let ternary = ternary_coefficients::<DEGREE, R>(hamming_weight, rng);
        let coeffs = ternary.map(|x| if x == -1 { context - 1 } else { x as u64 });
        Self { coeffs, context }
    }

    fn sample_noise<R: Rng>(rng: &mut R, variance: f64) -> Self {
        Self::sample_gaussian(rng, variance.sqrt())
    }
}
