use super::builder::CkksEngineBuilder;
use super::types::{Ciphertext, Plaintext};
use crate::rings::backends::rns_ntt::{RnsNttError, RnsPoly};
use crate::rings::traits::PolyAutomorphism;
use crate::{
    PolyRing, PolySampler, PublicKey, PublicKeyError, PublicKeyParams,
    RelinearizationKey, RelinearizationKeyError, RelinearizationKeyParams,
    SecretKey, SecretKeyError, SecretKeyParams,
};
use rand::Rng;
use std::sync::Arc;

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
        logq: u32,
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
            logq,                       // Modulus level
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

    // Homomorphic operations
    pub fn add_ciphertexts(
        ct1: &Ciphertext<P, DEGREE>,
        ct2: &Ciphertext<P, DEGREE>,
    ) -> Ciphertext<P, DEGREE> {
        let mut c0 = ct1.c0.clone();
        c0 += &ct2.c0;

        let mut c1 = ct1.c1.clone();
        c1 += &ct2.c1;

        // Both logp and logq should match for addition
        assert_eq!(ct1.logp, ct2.logp, "logp mismatch in addition");
        assert_eq!(ct1.logq, ct2.logq, "logq mismatch in addition");

        Ciphertext {
            c0,
            c1,
            logp: ct1.logp,
            logq: ct1.logq,
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

        // Step 3: Scale calculation - both logp and logq double
        // After multiplication, scale becomes 2^(logp1 + logp2)
        // Modulus stays at logq1 + logq2 (or should match for valid multiplication)
        assert_eq!(ct1.logq, ct2.logq, "logq mismatch in multiplication");
        let new_logp = ct1.logp + ct2.logp;

        Ciphertext {
            c0: c0_new,
            c1: c1_new,
            logp: new_logp,
            logq: ct1.logq, // Modulus unchanged (rescaling would reduce this)
        }
    }
}

/// Relinearization key using RNS gadget decomposition.
///
/// Holds one `(a_i, b_i)` pair per RNS channel. Each pair satisfies
///
/// ```text
/// b_i + a_i · s ≈ e_i · s²   (mod Q)
/// ```
///
/// where `e_i` is the CRT basis polynomial that equals 1 in channel `i` and 0 in all
/// other channels. The identity `d2 = Σ_i α_i(d2) · e_i` then gives
///
/// ```text
/// Σ_i α_i(d2) · key_i ≈ d2 · s²
/// ```
///
/// with relin noise ≈ `Σ_i α_i(d2) · e_i_rk` where `|α_i| < q_i`. For typical toy
/// parameters (N=16, q_i ≈ 2^31, σ=3.2) the noise after rescaling is ~5 × 10⁻⁸ — comparable
/// to the encryption noise floor.
#[derive(Clone, Debug)]
pub struct RnsGadgetRelinKey<const DEGREE: usize> {
    pub a: Vec<RnsPoly<DEGREE>>,
    pub b: Vec<RnsPoly<DEGREE>>,
}

/// Rotation key using RNS gadget decomposition.
///
/// Structurally identical to `RnsGadgetRelinKey` but encodes `s_k = s(X^{5^k})`
/// (the rotated secret key) instead of `s²`. Each pair satisfies:
///
/// ```text
/// b_i + a_i · s ≈ e_i · s_k   (mod Q)
/// ```
///
/// After applying the automorphism σ_k to a ciphertext `(c0, c1)`, the result
/// is encrypted under `s_k` rather than `s`. Gadget key-switching converts it
/// back to encryption under `s`:
///
/// ```text
/// c0' = σ_k(c0) + Σ_i α_i(σ_k(c1)) · b_i
/// c1' =           Σ_i α_i(σ_k(c1)) · a_i
/// ```
#[derive(Clone, Debug)]
pub struct RnsGadgetRotationKey<const DEGREE: usize> {
    pub a: Vec<RnsPoly<DEGREE>>,
    pub b: Vec<RnsPoly<DEGREE>>,
    /// The slot rotation offset this key was generated for.
    pub rotation: i32,
}

impl<const DEGREE: usize> CkksEngine<RnsPoly<DEGREE>, DEGREE> {
    /// Drops the last RNS prime and divides coefficients by it, reducing `logp` and `logq`.
    ///
    /// Both ciphertext components receive the **same** `Arc<RnsBasis>` so that subsequent
    /// polynomial operations (which `Arc::ptr_eq`-check their bases) succeed without issue.
    ///
    /// Call this immediately after `mul_ciphertexts_gadget` to bring the ciphertext scale
    /// back from `2·logp` to `logp`.
    pub fn rescale_ciphertext(
        ct: &Ciphertext<RnsPoly<DEGREE>, DEGREE>,
    ) -> Result<Ciphertext<RnsPoly<DEGREE>, DEGREE>, RnsNttError> {
        let q_last = *ct.c0.basis().moduli().last().expect("at least one channel");
        // bit_length(q_last) = floor(log2(q_last)) + 1.
        // Using floor alone would give a factor-of-2 error: for q_last ≈ 2^31 the
        // floor is 30 but actual division by q_last shifts the scale by 31 bits.
        let bits_dropped = u64::BITS - q_last.leading_zeros();

        // One shared Arc so that c0 and c1 (and any later-created sk_reduced) satisfy
        // Arc::ptr_eq checks in mul_assign / add_assign.
        let new_basis = Arc::new(ct.c0.basis().drop_last(1)?);

        Ok(Ciphertext {
            c0: ct.c0.rescale_into(new_basis.clone())?,
            c1: ct.c1.rescale_into(new_basis)?,
            logp: ct.logp - bits_dropped,
            logq: ct.logq - bits_dropped,
        })
    }

    /// Generates a gadget relinearization key from a secret key.
    ///
    /// Produces one `(a_i, b_i)` pair per RNS channel where
    /// `b_i + a_i·s ≈ e_i·s²` and `e_i` is the CRT basis polynomial for channel `i`.
    pub fn generate_gadget_relin_key<R: Rng>(
        &self,
        sk: &SecretKey<RnsPoly<DEGREE>, DEGREE>,
        rng: &mut R,
    ) -> RnsGadgetRelinKey<DEGREE> {
        let basis = sk.poly.basis().clone();
        let l = basis.channel_count();

        // s² in coefficient domain.
        let mut s_sq = sk.poly.clone();
        s_sq *= &sk.poly;
        s_sq.to_coeff_domain();

        let mut a_vec = Vec::with_capacity(l);
        let mut b_vec = Vec::with_capacity(l);

        for i in 0..l {
            // Build the plaintext polynomial e_i · s²:
            // channel i = s²_i (coeff-domain slice of s²), all other channels = 0.
            let s2_ch = s_sq.channels();
            let mut plain_channels: Vec<[u64; DEGREE]> = vec![[0u64; DEGREE]; l];
            plain_channels[i] = s2_ch[i]; // [u64; DEGREE] is Copy
            let plain_s2_i =
                RnsPoly::from_channels(plain_channels, basis.clone(), false)
                    .expect(
                        "gadget plaintext: channel i is already reduced mod q_i",
                    );

            let a_i = RnsPoly::sample_uniform(&basis, rng);
            let e_i = RnsPoly::sample_gaussian(
                self.params.error_variance.sqrt(),
                &basis,
                rng,
            );

            // b_i = -(a_i · s) + e_i + plain_s2_i
            let mut a_times_s = a_i.clone();
            a_times_s *= &sk.poly;
            let mut b_i = -a_times_s;
            b_i += &e_i;
            b_i += &plain_s2_i;

            a_vec.push(a_i);
            b_vec.push(b_i);
        }

        RnsGadgetRelinKey { a: a_vec, b: b_vec }
    }

    /// Generates a gadget rotation key for the given slot offset.
    ///
    /// Computes `s_k = s(X^{5^k mod 2N})` via the ring automorphism, then
    /// encodes it in the same RNS gadget structure used by the relinearization
    /// key. One `(a_i, b_i)` pair is produced per RNS channel, where:
    ///
    /// ```text
    /// b_i + a_i · s ≈ e_i · s_k
    /// ```
    ///
    /// Negative `rotation` rotates slots right (toward lower indices).
    pub fn generate_gadget_rotation_key<R: Rng>(
        &self,
        sk: &SecretKey<RnsPoly<DEGREE>, DEGREE>,
        rotation: i32,
        rng: &mut R,
    ) -> RnsGadgetRotationKey<DEGREE> {
        let basis = sk.poly.basis().clone();
        let l = basis.channel_count();

        // s_k = s(X^{5^rotation mod 2N}) in coefficient domain.
        // rotate_slots always returns coefficient domain.
        let s_k = sk.poly.rotate_slots(rotation);

        let mut a_vec = Vec::with_capacity(l);
        let mut b_vec = Vec::with_capacity(l);

        for i in 0..l {
            // Build e_i · s_k: channel i gets s_k's channel, all others zero.
            let s_k_channels = s_k.channels();
            let mut plain_channels: Vec<[u64; DEGREE]> = vec![[0u64; DEGREE]; l];
            plain_channels[i] = s_k_channels[i];
            let plain_s_k_i =
                RnsPoly::from_channels(plain_channels, basis.clone(), false)
                    .expect("rotation key plaintext: channel i is already reduced mod q_i");

            let a_i = RnsPoly::sample_uniform(&basis, rng);
            let e_i = RnsPoly::sample_gaussian(
                self.params.error_variance.sqrt(),
                &basis,
                rng,
            );

            // b_i = -(a_i · s) + e_i + plain_s_k_i
            let mut a_times_s = a_i.clone();
            a_times_s *= &sk.poly;
            let mut b_i = -a_times_s;
            b_i += &e_i;
            b_i += &plain_s_k_i;

            a_vec.push(a_i);
            b_vec.push(b_i);
        }

        RnsGadgetRotationKey { a: a_vec, b: b_vec, rotation }
    }

    /// Rotate ciphertext slots by `rotk.rotation` positions using gadget key switching.
    ///
    /// # Steps
    /// 1. Apply automorphism σ_k: `X → X^{5^k}` to both ciphertext components.
    ///    This rotates the plaintext slots but changes the decryption key from `s` to `s_k`.
    /// 2. Key-switch back to `s` using the gadget rotation key:
    ///    - Decompose `σ_k(c1)` channel-by-channel (same gadget decomposition as relinearization)
    ///    - Accumulate `Σ_i α_i(σ_k(c1)) · (rotk.b[i], rotk.a[i])` and add to `(σ_k(c0), 0)`
    ///
    /// The output is a fresh degree-1 ciphertext encrypted under `s` with the same
    /// `logp` and `logq` as the input — no level is consumed by rotation.
    pub fn rotate_ciphertext(
        ct: &Ciphertext<RnsPoly<DEGREE>, DEGREE>,
        rotk: &RnsGadgetRotationKey<DEGREE>,
    ) -> Ciphertext<RnsPoly<DEGREE>, DEGREE> {
        // Step 1: apply automorphism — result is in coefficient domain.
        let c0_rot = ct.c0.rotate_slots(rotk.rotation);
        let mut c1_rot = ct.c1.rotate_slots(rotk.rotation);
        c1_rot.to_coeff_domain(); // already coeff; explicit for clarity

        // Step 2: gadget key switch on c1_rot.
        let basis = ct.c0.basis().clone();
        let l = basis.channel_count();
        let moduli: Vec<u64> = basis.moduli().to_vec();

        let mut ks_c0 = RnsPoly::zero(basis.clone());
        let mut ks_c1 = RnsPoly::zero(basis.clone());

        for i in 0..l {
            // α_i(c1_rot): broadcast channel i across all channels, reducing mod q_j.
            let alpha_channels: Vec<[u64; DEGREE]> = (0..l)
                .map(|j| {
                    let qj = moduli[j];
                    let mut ch = [0u64; DEGREE];
                    for (k, v) in ch.iter_mut().enumerate() {
                        *v = c1_rot.channels()[i][k] % qj;
                    }
                    ch
                })
                .collect();
            let alpha_i =
                RnsPoly::from_channels(alpha_channels, basis.clone(), false)
                    .expect("alpha_i: residues are reduced by construction");

            let mut tb = alpha_i.clone();
            tb *= &rotk.b[i];
            ks_c0 += &tb;

            let mut ta = alpha_i;
            ta *= &rotk.a[i];
            ks_c1 += &ta;
        }

        let mut c0_new = c0_rot;
        c0_new += &ks_c0;

        Ciphertext {
            c0: c0_new,
            c1: ks_c1,
            logp: ct.logp,
            logq: ct.logq,
        }
    }

    /// Homomorphic multiplication using the RNS gadget relinearization key.
    ///
    /// Computes the tensor product `(c0 + c1·s)(c0' + c1'·s) = d0 + d1·s + d2·s²`,
    /// then replaces the `d2·s²` term with `Σ_i α_i(d2)·rlk_i`, yielding a valid
    /// degree-1 ciphertext.
    ///
    /// The output has `logp = ct1.logp + ct2.logp`; call `rescale_ciphertext` afterwards
    /// to bring it back to `ct1.logp`.
    pub fn mul_ciphertexts_gadget(
        ct1: &Ciphertext<RnsPoly<DEGREE>, DEGREE>,
        ct2: &Ciphertext<RnsPoly<DEGREE>, DEGREE>,
        rlk: &RnsGadgetRelinKey<DEGREE>,
    ) -> Ciphertext<RnsPoly<DEGREE>, DEGREE> {
        assert_eq!(ct1.logq, ct2.logq, "logq mismatch in gadget multiplication");

        // ── Tensor product ────────────────────────────────────────────────────
        let mut d0 = ct1.c0.clone();
        d0 *= &ct2.c0; // c0 · c0'

        let mut d1a = ct1.c0.clone();
        d1a *= &ct2.c1; // c0 · c1'
        let mut d1b = ct1.c1.clone();
        d1b *= &ct2.c0; // c1 · c0'
        let mut d1 = d1a;
        d1 += &d1b;

        let mut d2 = ct1.c1.clone();
        d2 *= &ct2.c1; // c1 · c1'
        d2.to_coeff_domain(); // must read channels[i] below

        // ── Gadget relinearization: replace d2·s² ────────────────────────────
        // c0 += Σ_i α_i(d2) · rlk.b[i]
        // c1 += Σ_i α_i(d2) · rlk.a[i]
        let basis = d2.basis().clone();
        let l = basis.channel_count();
        let moduli: Vec<u64> = basis.moduli().to_vec();

        let mut relin_c0 = RnsPoly::zero(basis.clone());
        let mut relin_c1 = RnsPoly::zero(basis.clone());

        for i in 0..l {
            // α_i(d2): broadcast channel i of d2 to all channels (reduce mod q_j).
            let alpha_channels: Vec<[u64; DEGREE]> = (0..l)
                .map(|j| {
                    let qj = moduli[j];
                    let mut ch = [0u64; DEGREE];
                    for (k, slot) in ch.iter_mut().enumerate() {
                        *slot = d2.channels()[i][k] % qj;
                    }
                    ch
                })
                .collect();
            let alpha_i =
                RnsPoly::from_channels(alpha_channels, basis.clone(), false)
                    .expect("alpha_i: residues are reduced by construction");

            let mut tb = alpha_i.clone();
            tb *= &rlk.b[i];
            relin_c0 += &tb;

            let mut ta = alpha_i;
            ta *= &rlk.a[i];
            relin_c1 += &ta;
        }

        d0 += &relin_c0;
        d1 += &relin_c1;

        Ciphertext {
            c0: d0,
            c1: d1,
            logp: ct1.logp + ct2.logp,
            logq: ct1.logq,
        }
    }
}
