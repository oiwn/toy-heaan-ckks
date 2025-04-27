use crate::PolyRing;
use std::fmt;

impl fmt::Display for PolyRing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Check if alternate format is requested (full polynomial)
        if f.alternate() {
            return self.fmt_full(f);
        }

        // Get number of coefficients to display from precision, default to 3
        let num_display = f.precision().unwrap_or(3);

        // Display truncated format
        self.fmt_truncated(f, num_display)
    }
}

impl PolyRing {
    // Format with truncated coefficient display
    fn fmt_truncated(
        &self,
        f: &mut fmt::Formatter<'_>,
        num_display: usize,
    ) -> fmt::Result {
        write!(
            f,
            "Poly(len={}, n={}, coeffs=[",
            self.len(),
            self.ring_dim()
        )?;

        let half_modulus = self.modulus() / 2;
        let len = self.len();

        if len <= num_display * 2 {
            // Display all coefficients if there are few enough
            let mut first = true;
            for &coeff in &self.coefficients {
                if !first {
                    write!(f, ", ")?;
                }
                first = false;

                // Convert to centered representation
                let value = if coeff > half_modulus {
                    -(self.modulus() as i64 - coeff as i64)
                } else {
                    coeff as i64
                };

                write!(f, "{}", value)?;
            }
        } else {
            // Display front and back coefficients with ellipsis
            for i in 0..num_display {
                if i > 0 {
                    write!(f, ", ")?;
                }

                let coeff = self.coefficients[i];
                let value = if coeff > half_modulus {
                    -(self.modulus() as i64 - coeff as i64)
                } else {
                    coeff as i64
                };

                write!(f, "{}", value)?;
            }

            // Ellipsis
            write!(f, ", …")?;

            // Show ending coefficients
            for i in (len - num_display)..len {
                write!(f, ", ")?;

                let coeff = self.coefficients[i];
                let value = if coeff > half_modulus {
                    -(self.modulus() as i64 - coeff as i64)
                } else {
                    coeff as i64
                };

                write!(f, "{}", value)?;
            }
        }

        write!(f, "])")
    }

    // Display the full polynomial in mathematical form
    fn fmt_full(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Keep your existing implementation here for the full polynomial display
        if self.coefficients.is_empty() {
            return write!(f, "0");
        }

        // Your existing code for full polynomial display...
        let half_modulus = self.modulus() / 2;
        let mut first = true;

        for (i, &coeff) in self.coefficients.iter().enumerate() {
            if coeff == 0 {
                continue;
            }

            // Convert to centered representation (-q/2, q/2)
            let value = if coeff > half_modulus {
                -(self.modulus() as i64 - coeff as i64)
            } else {
                coeff as i64
            };

            // Handle the sign
            if first {
                if value < 0 {
                    write!(f, "-")?;
                }
                first = false;
            } else if value < 0 {
                write!(f, " - ")?;
            } else {
                write!(f, " + ")?;
            };

            // Write coefficient and term
            let abs_value = value.abs();
            if i == 0 {
                // Constant term
                write!(f, "{}", abs_value)?;
            } else if i == 1 {
                // Linear term
                if abs_value == 1 {
                    write!(f, "x")?;
                } else {
                    write!(f, "{}*x", abs_value)?;
                }
            } else {
                // Higher degree terms
                if abs_value == 1 {
                    write!(f, "x^{}", i)?;
                } else {
                    write!(f, "{}x^{}", abs_value, i)?;
                }
            }
        }

        if first {
            // All coefficients were zero
            write!(f, "0")?;
        }

        Ok(())
    }
}
