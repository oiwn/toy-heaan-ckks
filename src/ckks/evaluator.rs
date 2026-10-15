//! Per-level key management: gadget relinearization and rotation keys are
//! level-bound, so evaluation past a rescale needs keys for the lower basis.

use std::collections::HashMap;

use rand::Rng;

use crate::ckks::algorithms::sum_slots::{required_rotation_offsets, sum_slots};
use crate::ckks::engine::{
    CkksEngine, CkksParams, RnsGadgetRelinKey, RnsGadgetRotationKey,
};
use crate::ckks::keys::SecretKey;
use crate::ckks::types::{Ciphertext, Plaintext};
use crate::rings::backends::rns_ntt::RnsPoly;

type Ct<const DEGREE: usize> = Ciphertext<RnsPoly<DEGREE>, DEGREE>;

/// Rotation keys indexed by (positive) slot offset.
pub type RotationKeyMap<const DEGREE: usize> =
    HashMap<i32, RnsGadgetRotationKey<DEGREE>>;

/// Generates the rotation key set needed to sum `slots` packed slots.
///
/// Convenience wrapper so callers do not hand-roll the offset list.
#[must_use]
pub fn generate_rotation_key_set<R: Rng, const DEGREE: usize>(
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

/// Ciphertext evaluator that hides the per-level key juggling.
///
/// A gadget key (relin or rotation) only decrypts correctly at the modulus
/// level it was generated for, so every ct x ct multiplication after a
/// rescale needs a fresh key for the reduced basis. This struct lazily
/// generates and caches those keys per level, driven by the master secret
/// key (which mod-switches exactly to any lower level).
pub struct LeveledEvaluator<const DEGREE: usize> {
    params: CkksParams<DEGREE>,
    sk: SecretKey<RnsPoly<DEGREE>, DEGREE>,
    /// Gadget relin keys indexed by level.
    rlks: HashMap<usize, RnsGadgetRelinKey<DEGREE>>,
    /// Gadget rotation keys indexed by (level, rotation offset).
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
    /// Panics if the operands sit at different levels or the rescale fails.
    pub fn mul<R: Rng>(
        &mut self,
        a: &Ct<DEGREE>,
        b: &Ct<DEGREE>,
        rng: &mut R,
    ) -> Ct<DEGREE> {
        let level = a.level();
        let rlk = self.rlks.remove(&level).unwrap_or_else(|| {
            let engine = CkksEngine::new(a.c0.basis().clone(), self.params.clone());
            engine.generate_gadget_relin_key(&self.secret_key_at(level), rng)
        });
        let product = a.mul(b, &rlk);
        self.rlks.insert(level, rlk);
        product.rescale().expect("LeveledEvaluator::mul: rescale")
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
        let level = ct.level();
        let mut rotkeys = RotationKeyMap::new();
        for offset in required_rotation_offsets(slots) {
            let key = self.rotkeys.remove(&(level, offset)).unwrap_or_else(|| {
                let engine =
                    CkksEngine::new(ct.c0.basis().clone(), self.params.clone());
                engine.generate_gadget_rotation_key(
                    &self.secret_key_at(level),
                    offset,
                    rng,
                )
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
        CkksEngine::<RnsPoly<DEGREE>, DEGREE>::decrypt(
            ct,
            &self.secret_key_at(ct.level()),
        )
    }

    /// The master secret key mod-switched down to `level`.
    fn secret_key_at(&self, level: usize) -> SecretKey<RnsPoly<DEGREE>, DEGREE> {
        self.sk
            .mod_switch_to(level)
            .expect("LeveledEvaluator: ciphertext above the secret key level")
    }
}
