//! Composable CKKS evaluation building blocks (slot sums, plaintext muls,
//! level alignment, Chebyshev evaluation).
//!
//! Everything here is built from engine primitives (`rotate`, `add`,
//! `mul_ciphertexts_gadget`, `mul_plaintext`, `rescale`) but encapsulates
//! the level/scale bookkeeping that the examples previously carried by hand:
//!
//! - [`sum_slots`] — replicate-and-add rotation tree over packed slots
//! - [`mul_plain_scalar`] — multiply every slot by a known constant
//! - [`align_level`] / [`reduce_secret_key`] — drop RNS channels so two
//!   ciphertexts (or a ciphertext and a key) live in the same basis
//! - [`decode_wide`] — decode a plaintext to real slot values
//! - [`fit_chebyshev`] + [`eval_chebyshev_clenshaw`] — Chebyshev-basis
//!   function approximation (plaintext reference for the encrypted
//!   evaluator [`eval_chebyshev_encrypted`])
//!
//! Application pipelines built from these ops (e.g. the sensor std-dev demo)
//! live in `examples/`.

use std::collections::HashMap;
use std::sync::Arc;

use rand::Rng;

use crate::crypto::engine::{
    CkksEngine, CkksParams, RnsGadgetRelinKey, RnsGadgetRotationKey,
};
use crate::crypto::types::{Ciphertext, Plaintext};
use crate::encoding::ckks_encoder::CkksEncoder;
use crate::keys::SecretKey;
use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

/// Rotation keys indexed by (positive) slot offset.
pub type RotationKeyMap<const DEGREE: usize> =
    HashMap<i32, RnsGadgetRotationKey<DEGREE>>;

/// Ciphertext type alias matching the RNS-NTT engine backend.
pub type Ct<const DEGREE: usize> = Ciphertext<RnsPoly<DEGREE>, DEGREE>;

// ── Slot aggregation ─────────────────────────────────────────────────────────

/// Rotation offsets needed by [`sum_slots`] for `slots` packed values:
/// `1, 2, 4, ..., slots/2`.
///
/// # Panics
/// Panics if `slots` is not a power-of-two >= 2.
#[must_use]
pub fn required_rotation_offsets(slots: usize) -> Vec<i32> {
    assert!(slots.is_power_of_two(), "slots must be a power of two");
    assert!(slots >= 2, "summing fewer than two slots is a no-op");

    (0..slots.ilog2()).map(|shift| 1i32 << shift).collect()
}

/// Generates the rotation key set needed to sum `slots` packed slots.
///
/// Convenience wrapper so callers do not hand-roll the offset list.
#[must_use]
pub fn generate_rotation_key_set<R: rand::Rng, const DEGREE: usize>(
    engine: &CkksEngine<RnsPoly<DEGREE>, DEGREE>,
    sk: &SecretKey<RnsPoly<DEGREE>, DEGREE>,
    slots: usize,
    rng: &mut R,
) -> RotationKeyMap<DEGREE> {
    required_rotation_offsets(slots)
        .into_iter()
        .map(|offset| {
            (offset, engine.generate_gadget_rotation_key(sk, offset, rng))
        })
        .collect()
}

/// Sums windows of `slots` consecutive slots via a binary rotation tree.
///
/// Slot `i` of the result holds `x_i + ... + x_{i+slots-1}` (indices mod
/// `DEGREE / 2`, since rotations wrap over all CKKS slots). So
/// `Enc([x0, ..., x_{slots-1}])` becomes `Enc([sum, ..., sum])` in every
/// slot only if `slots == DEGREE / 2` or the input is tiled (slot `i` holds
/// `x_{i mod slots}`). With unfilled slots left at zero, only slot 0 holds
/// the full sum. Consumes no levels; the level and basis are unchanged.
///
/// # Panics
/// Panics if `slots` is not a power of two, exceeds `DEGREE / 2`, or a
/// required rotation key is missing from `rotkeys`.
pub fn sum_slots<const DEGREE: usize>(
    ct: &Ct<DEGREE>,
    slots: usize,
    rotkeys: &RotationKeyMap<DEGREE>,
) -> Ct<DEGREE> {
    assert!(slots.is_power_of_two(), "slots must be a power of two");
    assert!(
        slots <= DEGREE / 2,
        "slots cannot exceed available CKKS slots"
    );

    let mut acc = ct.clone();

    for offset in required_rotation_offsets(slots) {
        let rotkey = rotkeys
            .get(&offset)
            .unwrap_or_else(|| panic!("missing rotation key for offset {offset}"));
        let rotated =
            CkksEngine::<RnsPoly<DEGREE>, DEGREE>::rotate_ciphertext(&acc, rotkey);
        acc =
            CkksEngine::<RnsPoly<DEGREE>, DEGREE>::add_ciphertexts(&acc, &rotated);
    }

    acc
}

// ── Scalar multiplication ────────────────────────────────────────────────────

/// Multiplies every slot by a known constant and rescales back to one scale.
///
/// The scalar is encoded into all slots ([`CkksEncoder::encode_constant`],
/// an exact constant polynomial) at `scalar_scale_bits` precision in the
/// ciphertext's current basis, multiplied via [`CkksEngine::mul_plaintext`]
/// (no relinearization needed), then rescaled.
///
/// # Panics
/// Panics if `scalar * 2^scalar_scale_bits` does not fit the basis modulus
/// or the rescale fails.
pub fn mul_plain_scalar<const DEGREE: usize>(
    ct: &Ct<DEGREE>,
    scalar: f64,
    scalar_scale_bits: u32,
) -> Ct<DEGREE> {
    let encoder = CkksEncoder::<DEGREE>::new(scalar_scale_bits);
    let pt = encoder.encode_constant(scalar, ct.c0.basis().clone());
    mul_plaintext_and_rescale(ct, &pt, "mul_plain_scalar")
}

/// [`mul_plain_scalar`] for a full slot vector of plaintext constants
/// (inner product with a known plaintext, slot-wise).
///
/// Constants are encoded at `scalar_scale_bits` directly; large scalars are
/// fine as long as `|constant| * 2^scalar_scale_bits` fits the basis modulus.
///
/// # Panics
/// Panics if `constants.len()` exceeds `DEGREE / 2`, a scaled constant does
/// not fit the basis modulus, or the rescale fails.
pub fn mul_plain_slots<const DEGREE: usize>(
    ct: &Ct<DEGREE>,
    constants: &[f64],
    scalar_scale_bits: u32,
) -> Ct<DEGREE> {
    assert!(
        constants.len() <= DEGREE / 2,
        "mul_plain_slots: {} constants exceed max slots {}",
        constants.len(),
        DEGREE / 2
    );

    let encoder = CkksEncoder::<DEGREE>::new(scalar_scale_bits);
    let pt = encoder.encode(constants, ct.c0.basis().clone());
    mul_plaintext_and_rescale(ct, &pt, "mul_plain_slots")
}

/// `ct * pt`, rescaled, with logp lifted back to `ct.logp` if needed.
fn mul_plaintext_and_rescale<const DEGREE: usize>(
    ct: &Ct<DEGREE>,
    pt: &Plaintext<RnsPoly<DEGREE>, DEGREE>,
    caller: &str,
) -> Ct<DEGREE> {
    let product = CkksEngine::<RnsPoly<DEGREE>, DEGREE>::mul_plaintext(ct, pt);
    let rescaled =
        CkksEngine::<RnsPoly<DEGREE>, DEGREE>::rescale_ciphertext(&product)
            .unwrap_or_else(|err| panic!("{caller}: rescale: {err}"));

    // A plaintext encoded with fewer bits than the dropped prime leaves logp
    // short after the rescale; lift it back so callers see an unchanged
    // scale. Level-free.
    if rescaled.logp < ct.logp {
        mul_plain_exact(&rescaled, 1.0, ct.logp - rescaled.logp)
    } else {
        rescaled
    }
}

/// Plaintext multiplication **without** the rescale: `logp` grows by
/// `scalar_scale_bits` and one prime is NOT consumed.
///
/// Encoding scalar `1.0` at `k` bits gives the constant polynomial `2^k`, so
/// the slot values are unchanged while `logp += k` — a level-free scale
/// lift. Encoding `0.5` at 1 bit gives the constant `1`, halving the slot
/// values — the level-free way to fold the `t = two_t / 2` of the Clenshaw
/// final step. Both rely on the scalar filling every slot
/// ([`CkksEncoder::encode_constant`]).
///
/// # Panics
/// Panics if `scalar * 2^scalar_scale_bits` does not fit the basis modulus.
pub fn mul_plain_exact<const DEGREE: usize>(
    ct: &Ct<DEGREE>,
    scalar: f64,
    scalar_scale_bits: u32,
) -> Ct<DEGREE> {
    let encoder = CkksEncoder::<DEGREE>::new(scalar_scale_bits);
    let pt = encoder.encode_constant(scalar, ct.c0.basis().clone());
    CkksEngine::<RnsPoly<DEGREE>, DEGREE>::mul_plaintext(ct, &pt)
}

// ── Level alignment ──────────────────────────────────────────────────────────

/// Drops RNS channels of `ct` so it lives in `target_basis` (which must be a
/// prefix of the ciphertext's current basis).
///
/// This is exact modulus switching per channel — the surviving residues are
/// already reduced — so it matches what `rescale_ciphertext` leaves behind
/// and lets a fresh (higher-level) ciphertext meet a computed one.
pub fn align_level<const DEGREE: usize>(
    ct: &Ct<DEGREE>,
    target_basis: Arc<RnsBasis<DEGREE>>,
    target_logq: u32,
) -> Ct<DEGREE> {
    let level = target_basis.channel_count();
    debug_assert!(
        ct.c0.basis().channel_count() >= level,
        "align_level: target basis has more channels than the ciphertext"
    );

    Ciphertext {
        c0: reduce_poly_to_basis(&ct.c0, &target_basis),
        c1: reduce_poly_to_basis(&ct.c1, &target_basis),
        logp: ct.logp,
        logq: target_logq,
    }
}

/// [`align_level`] for secret keys — truncate `sk` to `basis` channels.
pub fn reduce_secret_key<const DEGREE: usize>(
    sk: &SecretKey<RnsPoly<DEGREE>, DEGREE>,
    basis: Arc<RnsBasis<DEGREE>>,
) -> SecretKey<RnsPoly<DEGREE>, DEGREE> {
    let level = basis.channel_count();
    let channels: Vec<[u64; DEGREE]> = sk.poly.channels()[..level].to_vec();
    let poly = RnsPoly::from_channels(channels, basis, sk.poly.is_ntt_domain())
        .expect("reduce_secret_key: channels are already correctly reduced");
    SecretKey { poly }
}

fn reduce_poly_to_basis<const DEGREE: usize>(
    poly: &RnsPoly<DEGREE>,
    basis: &Arc<RnsBasis<DEGREE>>,
) -> RnsPoly<DEGREE> {
    let level = basis.channel_count();
    let channels: Vec<[u64; DEGREE]> = poly.channels()[..level].to_vec();
    RnsPoly::from_channels(channels, basis.clone(), poly.is_ntt_domain())
        .expect("reduce_poly_to_basis: channels are already correctly reduced")
}

// ── Wide decoding ────────────────────────────────────────────────────────────

/// Decodes a plaintext to real slot values.
///
/// Thin wrapper over [`CkksEncoder::decode`], which reconstructs
/// coefficients in `f64` and so handles moduli far beyond 64 bits.
pub fn decode_wide<const DEGREE: usize>(
    pt: &Plaintext<RnsPoly<DEGREE>, DEGREE>,
) -> Vec<f64> {
    CkksEncoder::<DEGREE>::new(pt.scale_bits).decode(pt)
}

// ── Leveled evaluation ───────────────────────────────────────────────────────

/// Ciphertext evaluator that hides the per-level key juggling.
///
/// A gadget key (relin or rotation) only decrypts correctly at the modulus
/// level it was generated for, so every ct x ct multiplication after a
/// rescale needs a fresh key for the reduced basis. This struct lazily
/// generates and caches those keys per channel count, driven by the original
/// secret key (which reduces exactly to any prefix basis).
pub struct LeveledEvaluator<const DEGREE: usize> {
    params: CkksParams<DEGREE>,
    sk: SecretKey<RnsPoly<DEGREE>, DEGREE>,
    /// Gadget relin keys indexed by channel count.
    rlks: HashMap<usize, RnsGadgetRelinKey<DEGREE>>,
    /// Gadget rotation keys indexed by (channel count, rotation offset).
    rotkeys: HashMap<(usize, i32), RnsGadgetRotationKey<DEGREE>>,
}

impl<const DEGREE: usize> LeveledEvaluator<DEGREE> {
    /// # Panics
    /// Panics later (lazily, at key generation) if `sk` does not decrypt the
    /// evaluated ciphertexts.
    pub fn new(
        params: CkksParams<DEGREE>,
        sk: SecretKey<RnsPoly<DEGREE>, DEGREE>,
    ) -> Self {
        Self {
            params,
            sk,
            rlks: HashMap::new(),
            rotkeys: HashMap::new(),
        }
    }

    /// ct x ct multiplication + rescale at the operands' current level.
    ///
    /// # Panics
    /// Panics if the operands' bases differ or the rescale fails.
    pub fn mul<R: Rng>(
        &mut self,
        a: &Ct<DEGREE>,
        b: &Ct<DEGREE>,
        rng: &mut R,
    ) -> Ct<DEGREE> {
        debug_assert!(
            a.c0.basis().channel_count() == b.c0.basis().channel_count(),
            "LeveledEvaluator::mul: operands at different levels"
        );

        let basis = a.c0.basis().clone();
        let level = basis.channel_count();
        let rlk = self.rlks.remove(&level).unwrap_or_else(|| {
            let engine = CkksEngine::new(basis.clone(), self.params.clone());
            let sk_l = reduce_secret_key(&self.sk, basis.clone());
            engine.generate_gadget_relin_key(&sk_l, rng)
        });
        let product = CkksEngine::<RnsPoly<DEGREE>, DEGREE>::mul_ciphertexts_gadget(
            a, b, &rlk,
        );
        self.rlks.insert(level, rlk);
        CkksEngine::<RnsPoly<DEGREE>, DEGREE>::rescale_ciphertext(&product)
            .expect("LeveledEvaluator::mul: rescale")
    }

    /// [`sum_slots`] with rotation keys generated (and cached) per level.
    ///
    /// # Panics
    /// Panics under the same conditions as [`sum_slots`].
    pub fn sum_slots<R: Rng>(
        &mut self,
        ct: &Ct<DEGREE>,
        slots: usize,
        rng: &mut R,
    ) -> Ct<DEGREE> {
        debug_assert!(
            !required_rotation_offsets(slots).is_empty(),
            "summing fewer than two slots is a no-op"
        );

        let basis = ct.c0.basis().clone();
        let level = basis.channel_count();
        let mut rotkeys = RotationKeyMap::new();
        for offset in required_rotation_offsets(slots) {
            let key = self.rotkeys.remove(&(level, offset)).unwrap_or_else(|| {
                let engine = CkksEngine::new(basis.clone(), self.params.clone());
                let sk_l = reduce_secret_key(&self.sk, basis.clone());
                engine.generate_gadget_rotation_key(&sk_l, offset, rng)
            });
            rotkeys.insert(offset, key);
        }

        let summed = sum_slots(ct, slots, &rotkeys);

        for (offset, key) in rotkeys {
            self.rotkeys.insert((level, offset), key);
        }
        summed
    }

    /// Decrypts `ct` by reducing the master secret key to its level.
    pub fn decrypt(&self, ct: &Ct<DEGREE>) -> Plaintext<RnsPoly<DEGREE>, DEGREE> {
        let sk_l = reduce_secret_key(&self.sk, ct.c0.basis().clone());
        CkksEngine::<RnsPoly<DEGREE>, DEGREE>::decrypt(ct, &sk_l)
    }
}

/// Builds a trivial ciphertext holding the constant `value` in every slot:
/// `(c0, c1) = (Enc(value), 0)`.
///
/// Decrypts to `value` under any key (c1 is zero), so adding/subtracting it
/// shifts every slot by `value` without a level or fresh encryption.
///
/// # Panics
/// Panics if `value` is not finite or `value * 2^scale_bits` does not fit
/// the basis modulus.
pub fn trivial_constant_ct<const DEGREE: usize>(
    value: f64,
    scale_bits: u32,
    basis: Arc<RnsBasis<DEGREE>>,
    logq: u32,
) -> Ct<DEGREE> {
    assert!(
        value.is_finite(),
        "trivial_constant_ct: value must be finite"
    );

    let encoder = CkksEncoder::<DEGREE>::new(scale_bits);
    let pt = encoder.encode_constant(value, basis.clone());

    let zeros = vec![[0u64; DEGREE]; basis.channel_count()];
    let c1 = RnsPoly::from_channels(zeros, basis, false)
        .expect("trivial_constant_ct: zeros are reduced");

    Ciphertext {
        c0: pt.poly,
        c1,
        logp: scale_bits,
        logq,
    }
}

// ── Encrypted Chebyshev evaluation ───────────────────────────────────────────

/// Evaluates `sum_k coeffs[k] * T_k(t)` on every encrypted slot of `v_ct`,
/// where `t` maps `[lower, upper]` to `[-1, 1]`, via encrypted Clenshaw.
///
/// Recurrence (Chebyshev basis throughout — converting to power basis is
/// ill-conditioned for small `lower`):
///
/// ```text
/// two_t = 2 * (2*v - (hi+lo)) / (hi-lo)     -- plaintext affine map
/// b_n   = c_n,  b_{n+1} = 0
/// b_k   = two_t * b_{k+1} - b_{k+2} + c_k   -- one ct x ct mul per step
/// result = t * b_1 - b_2 + c_0              -- final ct x ct mul
/// ```
///
/// `c_0` is not a recurrence term; the loop stops at `b_1`. `two_t * b_{k+1}`
/// is rescaled one level down, so `b_{k+2}` is dropped to the new level with
/// [`align_level`] (exact channel truncation) before the subtraction.
///
/// Levels consumed: `degree + 1` (affine map + `degree - 1` recurrence muls +
/// final mul). The final mul needs >= 3 RNS channels of modulus so its raw
/// product (scale ~2^(2*logp)) does not wrap.
///
/// # Panics
/// Panics if `coeffs` is empty or the ciphertext runs out of RNS channels
/// before the evaluation completes.
pub fn eval_chebyshev_encrypted<R: Rng, const DEGREE: usize>(
    v_ct: &Ct<DEGREE>,
    coeffs: &[f64],
    lower: f64,
    upper: f64,
    evaluator: &mut LeveledEvaluator<DEGREE>,
    rng: &mut R,
) -> Ct<DEGREE> {
    assert!(!coeffs.is_empty(), "coefficients must not be empty");
    let degree = coeffs.len() - 1;

    // two_t = 4/(hi-lo) * v - 2*(hi+lo)/(hi-lo), then Clenshaw needs only
    // ciphertext multiplications by two_t.
    let span = upper - lower;
    let two_t = mul_plain_scalar(v_ct, 4.0 / span, v_ct.logp);
    let two_t = sub_constant(&two_t, 2.0 * (upper + lower) / span);

    // b_n = c_n (constant), b_{n+1} = 0 — both as trivial ciphertexts at
    // two_t's level.
    let basis = two_t.c0.basis().clone();
    let logq = two_t.logq;
    let scale_bits = two_t.logp;
    let mut b_next =
        trivial_constant_ct(coeffs[degree], scale_bits, basis.clone(), logq);
    let mut b_next_next = trivial_constant_ct(0.0, scale_bits, basis, logq);

    // b_k = two_t * b_{k+1} - b_{k+2} + c_k, for k = n-1 down to 1.
    // (c_0 is NOT a recurrence term; it enters only the final step.)
    for &c_k in coeffs[1..degree].iter().rev() {
        // two_t must meet b_{k+1} at its level (exact truncation).
        let target = b_next.c0.basis().clone();
        let two_t_aligned = align_level(&two_t, target, b_next.logq);

        let mut b_k = evaluator.mul(&two_t_aligned, &b_next, rng);

        // b_{k+2}: drop to the level the multiplication just landed on.
        let target = b_k.c0.basis().clone();
        b_next_next = align_level(&b_next_next, target, b_k.logq);
        b_k = CkksEngine::<RnsPoly<DEGREE>, DEGREE>::sub_ciphertexts(
            &b_k,
            &b_next_next,
        );
        b_k = add_constant(&b_k, c_k);

        b_next_next = b_next;
        b_next = b_k;
    }

    // result = t * b_1 - b_2 + c_0 (b_next = b_1, b_next_next = b_2), with
    // t = two_t / 2. The halving is level-free: product * (0.5 encoded at
    // 1 bit -> integer 1) halves the slots and lifts logp by 1, while b_2 *
    // (1.0 at 1 bit -> integer 2) keeps its slots unchanged at the lifted
    // scale — so all three terms meet at logp + 1 without consuming a prime.
    let target = b_next.c0.basis().clone();
    let two_t_aligned = align_level(&two_t, target, b_next.logq);
    let product = evaluator.mul(&two_t_aligned, &b_next, rng);

    let t_times_b0 = mul_plain_exact(&product, 0.5, 1);

    let b_one =
        align_level(&b_next_next, t_times_b0.c0.basis().clone(), t_times_b0.logq);
    let b_one = mul_plain_exact(&b_one, 1.0, 1);

    let result =
        CkksEngine::<RnsPoly<DEGREE>, DEGREE>::sub_ciphertexts(&t_times_b0, &b_one);
    add_constant(&result, coeffs[0])
}

/// Shifts every slot down by `value` using a trivial ciphertext; consumes no
/// level.
fn sub_constant<const DEGREE: usize>(ct: &Ct<DEGREE>, value: f64) -> Ct<DEGREE> {
    add_constant(ct, -value)
}

/// Shifts every slot up by `value` using a trivial ciphertext; consumes no
/// level.
fn add_constant<const DEGREE: usize>(ct: &Ct<DEGREE>, value: f64) -> Ct<DEGREE> {
    let pos = trivial_constant_ct(value, ct.logp, ct.c0.basis().clone(), ct.logq);
    CkksEngine::<RnsPoly<DEGREE>, DEGREE>::add_ciphertexts(ct, &pos)
}

/// Maps `x` in `[lower, upper]` to `t` in `[-1, 1]` (Chebyshev domain).
pub fn to_chebyshev_domain(x: f64, lower: f64, upper: f64) -> f64 {
    (2.0 * x - (upper + lower)) / (upper - lower)
}

/// Inverse of [`to_chebyshev_domain`].
pub fn from_chebyshev_domain(t: f64, lower: f64, upper: f64) -> f64 {
    0.5 * ((upper - lower) * t + upper + lower)
}

/// Chebyshev interpolation of `f` on `[lower, upper]`.
///
/// Samples `f` at the `degree + 1` Chebyshev nodes and returns coefficients
/// `c_k` with `f(x) ~ sum_k c_k * T_k(t(x))`, evaluated by
/// [`eval_chebyshev_clenshaw`] (plaintext) or [`eval_chebyshev_encrypted`].
/// Keeping the result in Chebyshev basis (never converting to power basis)
/// keeps the fit well conditioned, e.g. for `sqrt` near a small `lower`.
///
/// # Panics
/// Panics if `upper <= lower`.
pub fn fit_chebyshev(
    f: impl Fn(f64) -> f64,
    lower: f64,
    upper: f64,
    degree: usize,
) -> Vec<f64> {
    assert!(upper > lower, "fit_chebyshev: upper must exceed lower");

    let node_count = degree + 1;
    let mut values = Vec::with_capacity(node_count);

    for j in 0..node_count {
        let theta = std::f64::consts::PI * (j as f64 + 0.5) / node_count as f64;
        let t = theta.cos();
        let x = from_chebyshev_domain(t, lower, upper);
        values.push(f(x));
    }

    let mut coeffs = vec![0.0; node_count];
    for (k, coeff) in coeffs.iter_mut().enumerate() {
        let sum = values
            .iter()
            .enumerate()
            .map(|(j, value)| {
                let theta =
                    std::f64::consts::PI * (j as f64 + 0.5) / node_count as f64;
                value * (k as f64 * theta).cos()
            })
            .sum::<f64>();

        *coeff = if k == 0 {
            sum / node_count as f64
        } else {
            2.0 * sum / node_count as f64
        };
    }

    coeffs
}

/// Evaluates `sum_k c_k * T_k(t)` on `t` in `[-1, 1]` via the Clenshaw
/// recurrence.
///
/// This is the plaintext reference for the encrypted evaluator; it produces
/// bit-identical arithmetic (up to rounding) to the homomorphic version.
///
/// # Panics
/// Panics if `coeffs` is empty.
pub fn eval_chebyshev_clenshaw(t: f64, coeffs: &[f64]) -> f64 {
    assert!(!coeffs.is_empty(), "coefficients must not be empty");

    let mut b_next = 0.0;
    let mut b_next_next = 0.0;

    for &coeff in coeffs.iter().skip(1).rev() {
        let b = 2.0 * t * b_next - b_next_next + coeff;
        b_next_next = b_next;
        b_next = b;
    }

    t * b_next - b_next_next + coeffs[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::generate_primes;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    /// Encrypted Clenshaw must track plaintext Clenshaw (same coefficients,
    /// same domain map) to within the accumulated CKKS noise.
    ///
    /// Toy params: N = 16, degree 3 -> levels = affine(1) + steps(2) +
    /// final(1) = 4. The final mul must run at >= 3 channels: its raw
    /// product sits at scale 2^(2*61) = 2^122, which only fits the modulus
    /// (Q ~ 2^(61*chans)) without wrapping when chans >= 3. Budget
    /// degree + 4 primes.
    #[test]
    fn encrypted_clenshaw_matches_plaintext() {
        const N: usize = 16;
        const SLOTS: usize = 8;
        const SCALE: u32 = 61;

        let (lower, upper) = (0.01, 1.0);
        let degree = 3;
        let coeffs = fit_chebyshev(f64::sqrt, lower, upper, degree);

        let input: [f64; SLOTS] = [0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95];
        let expected: Vec<f64> = input
            .iter()
            .map(|&x| {
                let t = to_chebyshev_domain(x, lower, upper);
                eval_chebyshev_clenshaw(t, &coeffs)
            })
            .collect();

        let primes = generate_primes(SCALE as usize, degree + 4, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
        let logq = basis.total_bits();

        let params = CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        };
        let engine =
            CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
        let mut rng = ChaCha20Rng::seed_from_u64(1234);

        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");

        let encoder = CkksEncoder::<N>::new(SCALE);
        let pt = encoder.encode(&input, basis);
        let ct = engine.encrypt(&pt, &pk, logq, &mut rng);

        let mut evaluator = LeveledEvaluator::new(params, sk);
        let result = eval_chebyshev_encrypted(
            &ct,
            &coeffs,
            lower,
            upper,
            &mut evaluator,
            &mut rng,
        );

        assert_eq!(
            result.c0.basis().channel_count(),
            3,
            "degree-3 eval from 7 primes: affine + 2 steps + final = 4 levels"
        );

        let decoded = decode_wide(&evaluator.decrypt(&result));
        for (slot, (&got, &want)) in decoded.iter().zip(&expected).enumerate() {
            assert!(
                (got - want).abs() < 1e-3,
                "slot {slot}: encrypted {got} vs plaintext {want}"
            );
        }
    }

    /// Regression: the input fills only 2 of 8 slots (zero-padded). Scalars
    /// used to be broadcast to just the filled slots, which made the
    /// level-free `0.5 @ 1 bit` halving round to the zero polynomial.
    #[test]
    fn encrypted_clenshaw_with_partially_filled_slots() {
        const N: usize = 16;
        const SCALE: u32 = 61;

        let (lower, upper) = (1e-4, 6.4e-3);
        let degree = 3;
        let coeffs = fit_chebyshev(f64::sqrt, lower, upper, degree);
        let input = [8.2e-4, 3.0e-3];

        let primes = generate_primes(SCALE as usize, degree + 4, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
        let logq = basis.total_bits();
        let params = CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        };
        let engine =
            CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
        let mut rng = ChaCha20Rng::seed_from_u64(99);
        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");

        let pt = CkksEncoder::<N>::new(SCALE).encode(&input, basis);
        let ct = engine.encrypt(&pt, &pk, logq, &mut rng);

        let mut evaluator = LeveledEvaluator::new(params, sk);
        let result = eval_chebyshev_encrypted(
            &ct,
            &coeffs,
            lower,
            upper,
            &mut evaluator,
            &mut rng,
        );

        let decoded = decode_wide(&evaluator.decrypt(&result));
        for (slot, &x) in input.iter().enumerate() {
            let want = eval_chebyshev_clenshaw(
                to_chebyshev_domain(x, lower, upper),
                &coeffs,
            );
            assert!(
                (decoded[slot] - want).abs() < 1e-9,
                "slot {slot}: encrypted {} vs plaintext {want}",
                decoded[slot]
            );
        }
    }

    #[test]
    fn chebyshev_domain_map_roundtrips() {
        for x in [0.01, 0.1, 0.5, 0.99, 1.0] {
            let t = to_chebyshev_domain(x, 0.01, 1.0);
            assert!(
                (-1.0..=1.0).contains(&t),
                "x={x} maps outside [-1, 1]: t={t}"
            );
            assert!((from_chebyshev_domain(t, 0.01, 1.0) - x).abs() < 1e-15);
        }
    }

    #[test]
    fn fit_sqrt_degree7_accuracy() {
        let coeffs = fit_chebyshev(f64::sqrt, 0.01, 1.0, 7);

        let mut max_err = 0.0f64;
        for i in 0..2001 {
            let x = 0.01 + (1.0 - 0.01) * i as f64 / 2000.0;
            let t = to_chebyshev_domain(x, 0.01, 1.0);
            let err = (eval_chebyshev_clenshaw(t, &coeffs) - x.sqrt()).abs();
            max_err = max_err.max(err);
        }

        // Degree-7 fit of sqrt on [0.01, 1]: worst error near the lower
        // bound.
        assert!(max_err < 0.02, "degree-7 fit max error {max_err:.4e}");
    }

    #[test]
    fn clenshaw_matches_direct_chebyshev_sum() {
        let coeffs = fit_chebyshev(f64::sqrt, 0.01, 1.0, 7);

        for i in -5..=5 {
            let t = i as f64 / 5.0;

            // Direct sum c_k * T_k(t) via the T recurrence.
            let (mut t_prev, mut t_cur) = (1.0, t);
            let mut direct = coeffs[0] + coeffs[1] * t;
            for &c in coeffs.iter().skip(2) {
                let t_next = 2.0 * t * t_cur - t_prev;
                direct += c * t_next;
                t_prev = t_cur;
                t_cur = t_next;
            }

            let clenshaw = eval_chebyshev_clenshaw(t, &coeffs);
            assert!(
                (direct - clenshaw).abs() < 1e-12,
                "t={t}: direct {direct} vs clenshaw {clenshaw}"
            );
        }
    }
}
