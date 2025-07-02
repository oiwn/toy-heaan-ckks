use super::RnsPolyRing;
use std::fmt;

impl<const DEGREE: usize> fmt::Display for RnsPolyRing<DEGREE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            return self.fmt_full(f);
        }
        let num = f.precision().unwrap_or(3);
        self.fmt_truncated(f, num)
    }
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// Truncated display of RNS-polynomial: shows first `num` and last `num` residues*slots
    fn fmt_truncated(&self, f: &mut fmt::Formatter<'_>, num: usize) -> fmt::Result {
        let deg = DEGREE;
        let channels = self.basis.channel_count();
        write!(f, "RnsPoly<{}>[", deg)?;

        if deg <= num * 2 {
            // show all slots
            for i in 0..deg {
                if i > 0 {
                    write!(f, ", ")?;
                }
                // gather residues for slot i
                write!(f, "[")?;
                for c in 0..channels {
                    if c > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", self.coefficients[c][i])?;
                }
                write!(f, "]")?;
                if i > 0 {
                    write!(f, "*x^{}", i)?;
                }
            }
        } else {
            // front
            for i in 0..num {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "[")?;
                for c in 0..channels {
                    if c > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", self.coefficients[c][i])?;
                }
                write!(f, "]*x^{}", i)?;
            }
            write!(f, ", …")?;
            // back
            for i in (deg - num)..deg {
                write!(f, ", [")?;
                for c in 0..channels {
                    if c > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", self.coefficients[c][i])?;
                }
                write!(f, "]*x^{}", i)?;
            }
        }

        write!(f, "]")
    }

    /// Full display of RNS-polynomial: all residues*slots in polynomial form
    fn fmt_full(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let deg = DEGREE;
        let channels = self.basis.channel_count();
        let mut first = true;
        for i in 0..deg {
            // gather residues
            let mut slot = String::new();
            slot.push('[');
            for c in 0..channels {
                if c > 0 {
                    slot.push_str(", ");
                }
                slot.push_str(&self.coefficients[c][i].to_string());
            }
            slot.push(']');

            if !first {
                write!(f, " + ")?;
            }
            first = false;
            write!(f, "{}", slot)?;
            if i > 0 {
                write!(f, "*x^{}", i)?;
            }
        }
        Ok(())
    }
}
