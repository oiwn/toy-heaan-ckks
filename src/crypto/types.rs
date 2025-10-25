use crate::{PolyModSwitch, PolyRing};

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
/// Following HEAAN's design, we separate precision and modulus parameters:
/// - `logp`: Precision parameter (stays constant during modulus switching)
/// - `logq`: Current modulus level in bits (changes during modulus switching)
///
/// This approach cleanly separates:
/// - **Modulus switching** (noise management): only `logq` changes
/// - **Rescaling** (after multiplication): both `logp` and `logq` change
///
/// # Reference
/// - HEAAN implementation: https://github.com/kimandrik/HEAAN
/// - FHE Textbook: c04-glwe-modulus-switch.tex
#[derive(Debug)]
pub struct Ciphertext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub c0: P,
    pub c1: P,
    pub logp: u32, // Precision bits (encoding scale = 2^logp)
    pub logq: u32, // Modulus bits (current modulus level)
}

impl<P, const DEGREE: usize> Ciphertext<P, DEGREE>
where
    P: PolyModSwitch<DEGREE>,
{
    /// Modulus switching for ciphertexts (HEAAN's modDown operation)
    ///
    /// Applies modulus switching to both ciphertext components (c0, c1) for
    /// noise management. Following HEAAN's design:
    /// - Precision parameter `logp` stays UNCHANGED
    /// - Modulus level `logq` is updated to reflect new modulus
    /// - Same secret key decrypts at the new modulus
    ///
    /// This differs from rescaling (which changes both logp and logq).
    ///
    /// # Parameters
    /// - `new_context`: Target modulus context
    /// - `new_logq`: New modulus level in bits (log2 of new modulus)
    ///
    /// # Returns
    /// New ciphertext at the target modulus with unchanged precision
    ///
    /// # Reference
    /// - HEAAN: modDownBy/modDownTo operations
    /// - FHE Textbook: c04-glwe-modulus-switch.tex
    pub fn mod_switch(&self, new_context: &P::Context, new_logq: u32) -> Self {
        Self {
            c0: self.c0.mod_switch(new_context),
            c1: self.c1.mod_switch(new_context),
            logp: self.logp, // Precision parameter unchanged
            logq: new_logq,  // Modulus level updated
        }
    }
}
