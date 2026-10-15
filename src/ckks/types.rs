use crate::PolyRing;
use crate::ckks::errors::CkksResult;
use crate::rings::backends::rns_ntt::RnsPoly;

#[derive(Debug)]
pub struct Plaintext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub poly: P,
    pub scale_bits: u32,
    pub slots: usize, // Number of encoded slots (determines decode output length)
}

/// CKKS Ciphertext with HEAAN-style scale tracking
///
/// `logp` is the precision: the encoding scale is `2^logp`. The modulus
/// level is not stored; it is read off the basis (`level()`, `logq()`).
///
/// - **Mod switch** drops primes without dividing: `logp` stays.
/// - **Rescale** divides by the last prime: `logp` and `logq` both shrink.
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
    pub logp: u32, // Precision bits (encoding scale = 2^logp)
}

impl<const DEGREE: usize> Ciphertext<RnsPoly<DEGREE>, DEGREE> {
    /// Number of RNS channels (primes) the ciphertext still carries.
    pub fn level(&self) -> usize {
        self.c0.basis().channel_count()
    }

    /// Modulus bits at the current level (floor-based, see
    /// `RnsBasis::total_bits`).
    pub fn logq(&self) -> u32 {
        self.c0.basis().total_bits()
    }

    /// Divides by the last prime and drops it, lowering `logp` by its bit
    /// length.
    ///
    /// Call this after a multiplication to bring the scale from `2 * logp`
    /// back to about `logp`. The result lives on the shared chain basis one
    /// level down.
    ///
    /// # Errors
    /// Returns an error at level 1 (no prime left to drop).
    pub fn rescale(&self) -> CkksResult<Self> {
        let q_last = *self.c0.basis().moduli().last().expect("non-empty basis");
        // Bit length, floor(log2(q_last)) + 1: division by a ~2^61 prime
        // lowers the scale by about 61 bits.
        let bits_dropped = u64::BITS - q_last.leading_zeros();
        Ok(Self {
            c0: self.c0.rescale()?,
            c1: self.c1.rescale()?,
            logp: self.logp - bits_dropped,
        })
    }

    /// Drops primes down to `level` without dividing; `logp` is unchanged and
    /// the same secret key still decrypts.
    ///
    /// # Errors
    /// Returns an error if `level` is 0 or above the current level.
    pub fn mod_switch_to(&self, level: usize) -> CkksResult<Self> {
        Ok(Self {
            c0: self.c0.mod_switch_to(level)?,
            c1: self.c1.mod_switch_to(level)?,
            logp: self.logp,
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
