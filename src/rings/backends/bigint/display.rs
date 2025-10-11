use super::BigIntPolyRing;
use std::fmt;

pub struct CompactPolyDisplay<'a, const DEGREE: usize> {
    poly: &'a BigIntPolyRing<DEGREE>,
    precision: usize,
}

impl<'a, const DEGREE: usize> CompactPolyDisplay<'a, DEGREE> {
    pub fn new(poly: &'a BigIntPolyRing<DEGREE>, precision: usize) -> Self {
        Self { poly, precision }
    }
}

impl<const DEGREE: usize> fmt::Display for CompactPolyDisplay<'_, DEGREE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let coeffs = self.poly.coefficients();

        write!(f, "BigIntPoly[{}](", DEGREE)?;

        let n = self.precision;

        if DEGREE <= 2 * n {
            // Small polynomial - show all coefficients
            for (i, &coeff) in coeffs.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", coeff)?;
            }
        } else {
            // Large polynomial - show first n, ..., last n
            for i in 0..n {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", coeffs[i])?;
            }

            write!(f, ", ...")?;

            for i in 0..n {
                let idx = DEGREE - n + i;
                write!(f, ", {}", coeffs[idx])?;
            }
        }

        write!(f, ")")
    }
}

impl<const DEGREE: usize> fmt::Display for BigIntPolyRing<DEGREE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        CompactPolyDisplay::new(self, 3).fmt(f)
    }
}
