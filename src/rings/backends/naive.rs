use crate::PolyRing;
use std::ops::{AddAssign, MulAssign, Neg};

#[derive(Debug, Clone, PartialEq)]
pub struct NaivePolyRing<const DEGREE: usize> {
    pub coeffs: [u64; DEGREE],
    pub modulus: u64,
}

impl<const DEGREE: usize> NaivePolyRing<DEGREE> {
    pub fn default_modulus() -> u64 {
        return 741507920154517877u64;
    }
}

impl<const DEGREE: usize> PolyRing<DEGREE> for NaivePolyRing<DEGREE> {
    type Context = u64; // Just the modulus for naive implementation

    fn zero() -> Self {
        Self {
            coeffs: [0; DEGREE],
            modulus: Self::default_modulus(), // You'll need to define this
        }
    }

    fn from_coeffs(coeffs: &[u64]) -> Self {
        let mut result = [0u64; DEGREE];
        let len = coeffs.len().min(DEGREE);
        result[..len].copy_from_slice(&coeffs[..len]);

        Self {
            coeffs: result,
            modulus: Self::default_modulus(),
        }
    }

    fn to_coeffs(&self) -> [u64; DEGREE] {
        self.coeffs
    }
}

impl<const DEGREE: usize> AddAssign<&Self> for NaivePolyRing<DEGREE> {
    fn add_assign(&mut self, rhs: &Self) {
        for i in 0..DEGREE {
            self.coeffs[i] = (self.coeffs[i] + rhs.coeffs[i]) % self.modulus;
        }
    }
}

impl<const DEGREE: usize> MulAssign<&Self> for NaivePolyRing<DEGREE> {
    fn mul_assign(&mut self, rhs: &Self) {
        // Schoolbook multiplication with X^DEGREE + 1 reduction
        let mut result = [0u64; DEGREE];

        for i in 0..DEGREE {
            for j in 0..DEGREE {
                let coeff_pos = i + j;
                let coeff_val = (self.coeffs[i] as u128 * rhs.coeffs[j] as u128)
                    % self.modulus as u128;

                if coeff_pos < DEGREE {
                    result[coeff_pos] = ((result[coeff_pos] as u128 + coeff_val)
                        % self.modulus as u128)
                        as u64;
                } else {
                    // X^DEGREE = -1, so X^(DEGREE+k) = -X^k
                    let wrapped_pos = coeff_pos - DEGREE;
                    result[wrapped_pos] =
                        ((result[wrapped_pos] as u128 + self.modulus as u128
                            - coeff_val)
                            % self.modulus as u128) as u64;
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
            *coeff = (self.modulus - *coeff) % self.modulus;
        }
        self
    }
}
