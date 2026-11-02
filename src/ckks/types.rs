use crate::PolyRing;
use crate::ckks::Scale;
use crate::ckks::errors::CkksResult;

#[derive(Debug)]
pub struct Plaintext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub poly: P,
    /// Encoding factor: the polynomial holds `round(value * scale)`.
    pub scale: Scale,
    pub slots: usize, // Number of encoded slots (determines decode output length)
}

/// CKKS Ciphertext with HEAAN-style scale tracking
///
/// `scale` is the encoding factor the decrypted polynomial is divided by.
/// The modulus level is not stored; it is read off the basis (`level()`,
/// `logq()`).
///
/// - **Mod switch** drops primes without dividing: `scale` stays.
/// - **Rescale** divides by the last prime: `scale` and `logq` both shrink.
///
/// # Reference
/// - HEAAN implementation: <https://github.com/kimandrik/HEAAN>
/// - FHE Textbook: c04-glwe-modulus-switch.tex
#[derive(Clone, Debug)]
pub struct Ciphertext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub c0: P,
    pub c1: P,
    /// Encoding factor of the decrypted message.
    pub scale: Scale,
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Ciphertext<P, DEGREE> {
    /// Position in the modulus chain (for RNS: primes still carried).
    pub fn level(&self) -> usize {
        self.c0.level()
    }

    /// Modulus bits at the current level (a lower bound on `log2(Q)`).
    pub fn logq(&self) -> u32 {
        self.c0.modulus_bits()
    }

    /// Divides by the last prime and drops it.
    ///
    /// Call this after a multiplication to bring the scale from `scale^2`
    /// back to about `scale`. The result lives on the shared chain basis one
    /// level down.
    ///
    /// The recorded scale is divided by the actual prime. The prime is
    /// converted to `f64` (relative rounding about 2^-53), the same way
    /// polynomial plans compute it, so planned and actual scales agree.
    ///
    /// # Errors
    /// Returns an error at level 1 (no prime left to drop), or if the scale
    /// would underflow.
    pub fn rescale(&self) -> CkksResult<Self> {
        let divisor = Scale::new(self.c0.rescale_divisor() as f64)?;
        Ok(Self {
            c0: self.c0.rescale()?,
            c1: self.c1.rescale()?,
            scale: self.scale.checked_div(divisor)?,
        })
    }

    /// Drops primes down to `level` without dividing; `scale` is unchanged
    /// and the same secret key still decrypts.
    ///
    /// # Errors
    /// Returns an error if `level` is 0 or above the current level.
    pub fn mod_switch_to(&self, level: usize) -> CkksResult<Self> {
        Ok(Self {
            c0: self.c0.mod_switch_to(level)?,
            c1: self.c1.mod_switch_to(level)?,
            scale: self.scale,
        })
    }

    /// Mod-switches down to `other`'s level so the two can be combined.
    ///
    /// # Errors
    /// Returns an error if `other` sits at a higher level than `self`.
    pub fn match_level(&self, other: &Self) -> CkksResult<Self> {
        self.mod_switch_to(other.level())
    }
}
