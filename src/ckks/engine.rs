use super::builder::CkksEngineBuilder;
use super::types::{Ciphertext, Plaintext};
use crate::rings::backends::rns_ntt::RnsPoly;
use crate::rings::traits::PolyAutomorphism;
use crate::{
    PolyRing, PolySampler, PublicKey, PublicKeyError, PublicKeyParams,
    RelinearizationKey, RelinearizationKeyError, RelinearizationKeyParams,
    SecretKey, SecretKeyError, SecretKeyParams,
};
use rand::Rng;

pub struct CkksEngine<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    context: P::Context,
    pub params: CkksParams<DEGREE>,
}

#[derive(Debug, Clone)]
pub struct CkksParams<const DEGREE: usize> {
    pub error_variance: f64,
    pub hamming_weight: usize,
    pub scale_bits: u32,
}

impl<P, const DEGREE: usize> CkksEngine<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    pub fn builder() -> CkksEngineBuilder<DEGREE> {
        CkksEngineBuilder::new()
    }

    pub fn new(context: P::Context, params: CkksParams<DEGREE>) -> Self {
        Self { context, params }
    }

    pub fn context(&self) -> &P::Context {
        &self.context
    }

    pub fn generate_secret_key<R: Rng>(
        &self,
        rng: &mut R,
    ) -> Result<SecretKey<P, DEGREE>, SecretKeyError> {
        let sk_params = SecretKeyParams::new(self.params.hamming_weight)?;
        SecretKey::generate(&sk_params, &self.context, rng)
    }

    pub fn generate_public_key<R: Rng>(
        &self,
        secret_key: &SecretKey<P, DEGREE>,
        rng: &mut R,
    ) -> Result<PublicKey<P, DEGREE>, PublicKeyError> {
        let pk_params = PublicKeyParams::new(3.2)?;
        PublicKey::generate(secret_key, &pk_params, &self.context, rng)
    }

    /// Generate a relinearization key for ciphertext multiplication
    ///
    /// The relinearization key enables converting degree-2 ciphertexts (produced by multiplication)
    /// back to degree-1 ciphertexts while preserving the encrypted plaintext.
    ///
    /// # Arguments
    /// * `secret_key` - The secret key used to generate the relinearization key
    /// * `rng` - Random number generator for sampling
    ///
    /// # Returns
    /// * `Result<RelinearizationKey<P, DEGREE>, RelinearizationKeyError>` - The generated key or error
    pub fn generate_relinearization_key<R: Rng>(
        &self,
        secret_key: &SecretKey<P, DEGREE>,
        rng: &mut R,
    ) -> Result<RelinearizationKey<P, DEGREE>, RelinearizationKeyError> {
        // Use same error variance as for other key generation
        let relin_params =
            RelinearizationKeyParams::new(self.params.error_variance.sqrt())?;
        RelinearizationKey::generate(secret_key, &relin_params, &self.context, rng)
    }

    // Encryption/Decryption
    pub fn encrypt<R: Rng>(
        &self,
        plaintext: &Plaintext<P, DEGREE>,
        public_key: &PublicKey<P, DEGREE>,
        rng: &mut R,
    ) -> Ciphertext<P, DEGREE> {
        let u = P::sample_tribits(self.params.hamming_weight, &self.context, rng);
        let e0 = P::sample_gaussian(self.params.error_variance, &self.context, rng);
        let e1 = P::sample_gaussian(self.params.error_variance, &self.context, rng);

        // c0 = b * u + e0 + m
        let mut c0 = public_key.b.clone();
        c0 *= &u;
        c0 += &e0;
        c0 += &plaintext.poly;

        // c1 = a * u + e1
        let mut c1 = public_key.a.clone();
        c1 *= &u;
        c1 += &e1;

        Ciphertext {
            c0,
            c1,
            logp: plaintext.scale_bits, // Precision from encoding
        }
    }

    pub fn decrypt(
        ciphertext: &Ciphertext<P, DEGREE>,
        secret_key: &SecretKey<P, DEGREE>,
    ) -> Plaintext<P, DEGREE> {
        // m = c0 + c1 * s
        let mut result = ciphertext.c1.clone();
        result *= &secret_key.poly;
        result += &ciphertext.c0;

        Plaintext {
            poly: result,
            scale_bits: ciphertext.logp, // Use precision parameter
            slots: DEGREE / 2, // After decryption, assume max slots (decoder will handle actual slots)
        }
    }

    pub fn mul_ciphertexts(
        ct1: &Ciphertext<P, DEGREE>,
        ct2: &Ciphertext<P, DEGREE>,
        relin_key: &RelinearizationKey<P, DEGREE>,
    ) -> Ciphertext<P, DEGREE> {
        // Step 1: Raw multiplication
        // (c0 + c1 * s) * (c0' + c1's) = d0 + d1*s + d2*s^2

        // d0 = c0 * c0'
        let mut d0 = ct1.c0.clone();
        d0 *= &ct2.c0;

        // d1 = c0*c1' + c1*c0'
        let mut d1_part1 = ct1.c0.clone();
        d1_part1 *= &ct2.c1;
        let mut d1_part2 = ct1.c1.clone();
        d1_part2 *= &ct2.c0;
        let mut d1 = d1_part1;
        d1 += &d1_part2;

        // d2 = c1 * c1'
        let mut d2 = ct1.c1.clone();
        d2 *= &ct2.c1;

        // Step 2: Relinearization (internal)
        // Transform (d0, d1, d2) → (d0', d1') using relin_key
        // c0_new = d0 + rk.b * d2
        let mut rk_b_times_d2 = relin_key.b.clone();
        rk_b_times_d2 *= &d2;
        let mut c0_new = d0;
        c0_new += &rk_b_times_d2;

        // c1_new = d1 + rk.a * d2
        let mut rk_a_times_d2 = relin_key.a.clone();
        rk_a_times_d2 *= &d2;
        let mut c1_new = d1;
        c1_new += &rk_a_times_d2;

        // Step 3: Scale calculation
        // After multiplication, scale becomes 2^(logp1 + logp2); the modulus
        // level is unchanged (rescaling would reduce it)
        let new_logp = ct1.logp + ct2.logp;

        Ciphertext {
            c0: c0_new,
            c1: c1_new,
            logp: new_logp,
        }
    }
}

/// Digit width for base-2^`DIGIT_BITS` gadget decomposition.
///
/// Each key-switch digit is bounded by `2^DIGIT_BITS`, so the noise added per
/// accumulated term is ~`2^DIGIT_BITS · σ` instead of ~`q_i · σ` (the old
/// CRT-channel-as-digit scheme). With 61-bit primes this gives
/// `ceil(61 / 20) = 4` digits per channel (one fewer gadget key entry than
/// DIGIT_BITS = 15) at moderately higher noise per switch; demo params
/// (12 primes) have headroom for it.
pub const DIGIT_BITS: u32 = 20;

/// Number of base-2^`DIGIT_BITS` digits needed for residues in `[0, modulus)`.
fn digit_count(modulus: u64) -> usize {
    let modulus_bits = u64::BITS - modulus.leading_zeros();
    modulus_bits.div_ceil(DIGIT_BITS) as usize
}

/// `(a, b)` polynomial matrices indexed `[channel][digit]`.
type GadgetKeyMatrix<const DEGREE: usize> =
    (Vec<Vec<RnsPoly<DEGREE>>>, Vec<Vec<RnsPoly<DEGREE>>>);

/// Relinearization key using base-2^`DIGIT_BITS` RNS gadget decomposition.
///
/// `a[i][j]` / `b[i][j]` hold the key for digit `j` of channel `i`:
///
/// ```text
/// b[i][j] + a[i][j] · s ≈ e_i · 2^(j·DIGIT_BITS) · s²   (mod Q)
/// ```
///
/// where `e_i` is the CRT basis polynomial that equals 1 in channel `i` and 0
/// in all other channels. Writing the channel-`i` residues of `d2` as
/// `α_i = Σ_j digit_j(α_i) · 2^(j·DIGIT_BITS)` gives
///
/// ```text
/// Σ_i Σ_j digit_j(α_i(d2)) · key[i][j] ≈ d2 · s²
/// ```
///
/// with every digit coefficient below `2^DIGIT_BITS`, so the relin noise is
/// ~`num_channels · num_digits · 2^DIGIT_BITS · σ` — independent of the prime
/// size, unlike the previous per-channel decomposition where each "digit" was
/// a full ~61-bit residue.
#[derive(Clone, Debug)]
pub struct RnsGadgetRelinKey<const DEGREE: usize> {
    pub a: Vec<Vec<RnsPoly<DEGREE>>>,
    pub b: Vec<Vec<RnsPoly<DEGREE>>>,
}

/// Rotation key using base-2^`DIGIT_BITS` RNS gadget decomposition.
///
/// Structurally identical to `RnsGadgetRelinKey` but encodes `s_k = s(X^{5^k})`
/// (the rotated secret key) instead of `s²`. Entry `[i][j]` satisfies:
///
/// ```text
/// b[i][j] + a[i][j] · s ≈ e_i · 2^(j·DIGIT_BITS) · s_k   (mod Q)
/// ```
///
/// After applying the automorphism σ_k to a ciphertext `(c0, c1)`, the result
/// is encrypted under `s_k` rather than `s`. Gadget key-switching converts it
/// back to encryption under `s`:
///
/// ```text
/// c0' = σ_k(c0) + Σ_i Σ_j digit_j(α_i(σ_k(c1))) · b[i][j]
/// c1' =           Σ_i Σ_j digit_j(α_i(σ_k(c1))) · a[i][j]
/// ```
#[derive(Clone, Debug)]
pub struct RnsGadgetRotationKey<const DEGREE: usize> {
    pub a: Vec<Vec<RnsPoly<DEGREE>>>,
    pub b: Vec<Vec<RnsPoly<DEGREE>>>,
    /// The slot rotation offset this key was generated for.
    pub rotation: i32,
}

impl<const DEGREE: usize> CkksEngine<RnsPoly<DEGREE>, DEGREE> {
    /// Generates a gadget relinearization key from a secret key.
    ///
    /// Produces one `(a[i][j], b[i][j])` pair per RNS channel `i` and digit `j`
    /// where `b[i][j] + a[i][j]·s ≈ e_i·2^(j·DIGIT_BITS)·s²`.
    pub fn generate_gadget_relin_key<R: Rng>(
        &self,
        sk: &SecretKey<RnsPoly<DEGREE>, DEGREE>,
        rng: &mut R,
    ) -> RnsGadgetRelinKey<DEGREE> {
        // s² in coefficient domain.
        let mut s_sq = sk.poly.clone();
        s_sq *= &sk.poly;
        s_sq.to_coeff_domain();

        let (a, b) = self.generate_gadget_key_matrix(sk, &s_sq, rng);
        RnsGadgetRelinKey { a, b }
    }

    /// Generates a gadget rotation key for the given slot offset.
    ///
    /// Computes `s_k = s(X^{5^k mod 2N})` via the ring automorphism, then
    /// encodes it in the same digit × channel gadget structure used by the
    /// relinearization key:
    ///
    /// ```text
    /// b[i][j] + a[i][j] · s ≈ e_i · 2^(j·DIGIT_BITS) · s_k
    /// ```
    ///
    /// Negative `rotation` rotates slots right (toward lower indices).
    pub fn generate_gadget_rotation_key<R: Rng>(
        &self,
        sk: &SecretKey<RnsPoly<DEGREE>, DEGREE>,
        rotation: i32,
        rng: &mut R,
    ) -> RnsGadgetRotationKey<DEGREE> {
        // s_k = s(X^{5^rotation mod 2N}) in coefficient domain.
        // rotate_slots always returns coefficient domain.
        let s_k = sk.poly.rotate_slots(rotation);

        let (a, b) = self.generate_gadget_key_matrix(sk, &s_k, rng);
        RnsGadgetRotationKey { a, b, rotation }
    }

    /// Generates the digit × channel key matrix for gadget key switching.
    ///
    /// `target` (e.g. `s²` or a rotated secret key) must be in coefficient
    /// domain. Entry `[i][j]` satisfies:
    ///
    /// ```text
    /// b[i][j] + a[i][j] · s ≈ e_i · 2^(j·DIGIT_BITS) · target   (mod Q)
    /// ```
    ///
    /// Channel `i` gets `digit_count(q_i)` pairs, enough to cover any residue
    /// in `[0, q_i)` with base-2^`DIGIT_BITS` digits.
    fn generate_gadget_key_matrix<R: Rng>(
        &self,
        sk: &SecretKey<RnsPoly<DEGREE>, DEGREE>,
        target: &RnsPoly<DEGREE>,
        rng: &mut R,
    ) -> GadgetKeyMatrix<DEGREE> {
        debug_assert!(
            !target.is_ntt_domain(),
            "gadget key target must be in coefficient domain"
        );

        let basis = sk.poly.basis().clone();
        let l = basis.channel_count();
        let moduli: Vec<u64> = basis.moduli().to_vec();

        let mut a_matrix = Vec::with_capacity(l);
        let mut b_matrix = Vec::with_capacity(l);

        for i in 0..l {
            let q_i = moduli[i];
            let digits = digit_count(q_i);
            let mut a_row = Vec::with_capacity(digits);
            let mut b_row = Vec::with_capacity(digits);
            let mut pow2 = 1u64; // 2^(j·DIGIT_BITS) mod q_i

            for _ in 0..digits {
                // Plaintext for this digit: channel i = target_i · pow2,
                // all other channels = 0 (that is e_i · 2^(j·DIGIT_BITS) · target).
                let mut plain_channels: Vec<[u64; DEGREE]> =
                    vec![[0u64; DEGREE]; l];
                for (coeff, &t) in
                    plain_channels[i].iter_mut().zip(&target.channels()[i])
                {
                    *coeff = ((t as u128 * pow2 as u128) % q_i as u128) as u64;
                }
                let plain =
                    RnsPoly::from_channels(plain_channels, basis.clone(), false)
                        .expect("gadget plaintext: channel i is reduced mod q_i");

                let a_ij = RnsPoly::sample_uniform(&basis, rng);
                let e_ij = RnsPoly::sample_gaussian(
                    self.params.error_variance.sqrt(),
                    &basis,
                    rng,
                );

                // b = -(a · s) + e + plain
                let mut a_times_s = a_ij.clone();
                a_times_s *= &sk.poly;
                let mut b_ij = -a_times_s;
                b_ij += &e_ij;
                b_ij += &plain;

                a_row.push(a_ij);
                b_row.push(b_ij);

                pow2 = (((pow2 as u128) << DIGIT_BITS) % q_i as u128) as u64;
            }

            a_matrix.push(a_row);
            b_matrix.push(b_row);
        }

        (a_matrix, b_matrix)
    }

    /// Gadget key switch: computes the accumulated key contributions for `src`.
    ///
    /// `src` must be in coefficient domain. For each channel `i`, splits the
    /// residues into base-2^`DIGIT_BITS` digit polynomials and accumulates:
    ///
    /// ```text
    /// ks_c0 = Σ_i Σ_j digit_j(α_i(src)) · key_b[i][j]
    /// ks_c1 = Σ_i Σ_j digit_j(α_i(src)) · key_a[i][j]
    /// ```
    ///
    /// Every digit coefficient is below `2^DIGIT_BITS`, which keeps the noise
    /// picked up from the key errors small regardless of prime size.
    pub(crate) fn gadget_key_switch(
        src: &RnsPoly<DEGREE>,
        key_a: &[Vec<RnsPoly<DEGREE>>],
        key_b: &[Vec<RnsPoly<DEGREE>>],
    ) -> (RnsPoly<DEGREE>, RnsPoly<DEGREE>) {
        debug_assert!(
            !src.is_ntt_domain(),
            "gadget key switch source must be in coefficient domain"
        );

        let basis = src.basis().clone();
        let l = basis.channel_count();
        let moduli: Vec<u64> = basis.moduli().to_vec();
        let mask = (1u64 << DIGIT_BITS) - 1;

        let mut ks_c0 = RnsPoly::zero(basis.clone());
        let mut ks_c1 = RnsPoly::zero(basis.clone());

        for i in 0..l {
            debug_assert_eq!(
                key_a[i].len(),
                digit_count(moduli[i]),
                "gadget key digit count mismatch for channel {i}"
            );

            for (j, (a_ij, b_ij)) in key_a[i].iter().zip(&key_b[i]).enumerate() {
                // digit_j(α_i): the j-th base-2^DIGIT_BITS digit of each
                // channel-i residue.
                let shift = j as u32 * DIGIT_BITS;
                let mut digit_coeffs = [0u64; DEGREE];
                for (d, &c) in digit_coeffs.iter_mut().zip(&src.channels()[i]) {
                    *d = (c >> shift) & mask;
                }

                // Broadcast the digit polynomial to every channel.
                let alpha_channels: Vec<[u64; DEGREE]> = moduli
                    .iter()
                    .map(|&q| {
                        let mut ch = [0u64; DEGREE];
                        for (v, &d) in ch.iter_mut().zip(&digit_coeffs) {
                            *v = d % q;
                        }
                        ch
                    })
                    .collect();
                let alpha =
                    RnsPoly::from_channels(alpha_channels, basis.clone(), false)
                        .expect("gadget digits are reduced by construction");

                let mut tb = alpha.clone();
                tb *= b_ij;
                ks_c0 += &tb;

                let mut ta = alpha;
                ta *= a_ij;
                ks_c1 += &ta;
            }
        }

        (ks_c0, ks_c1)
    }
}
