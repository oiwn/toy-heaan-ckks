use super::{RnsBasis, RnsPolyRing};
use std::fmt;
use std::sync::Arc;

impl<const DEGREE: usize> fmt::Display for RnsPolyRing<DEGREE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Alternate (`{:#}`) triggers full expansion
        if f.alternate() {
            return self.fmt_full(f);
        }
        // Default: truncated with precision or 3
        let num = f.precision().unwrap_or(3);
        self.fmt_truncated(f, num)
    }
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// Truncated display: first `num` and last `num` CRT-reconstructed coefficients
    fn fmt_truncated(&self, f: &mut fmt::Formatter<'_>, num: usize) -> fmt::Result {
        let coeffs = self.crt();
        let len = coeffs.len();
        write!(f, "Poly<{}>[", DEGREE)?;

        if len <= num * 2 {
            // show all
            for (i, &c) in coeffs.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", c)?;
            }
        } else {
            // front
            for i in 0..num {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", coeffs[i])?;
            }
            write!(f, ", …")?;
            // back
            for i in (len - num)..len {
                write!(f, ", {}", coeffs[i])?;
            }
        }
        write!(f, "]")
    }

    /// Full display: CRT-reconstructed polynomial in standard form
    fn fmt_full(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let coeffs = self.crt();
        // All zero?
        if coeffs.iter().all(|&c| c == 0) {
            return write!(f, "0");
        }
        let mut first = true;
        for (i, &c) in coeffs.iter().enumerate() {
            if c == 0 {
                continue;
            }
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            match i {
                0 => write!(f, "{}", c)?,
                1 => {
                    if c == 1 {
                        write!(f, "x")?;
                    } else {
                        write!(f, "{}*x", c)?;
                    }
                }
                _ => {
                    if c == 1 {
                        write!(f, "x^{}", i)?;
                    } else {
                        write!(f, "{}*x^{}", c, i)?;
                    }
                }
            }
        }
        Ok(())
    }
}
