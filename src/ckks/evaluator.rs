//! Per-level key management: gadget relinearization and rotation keys are
//! level-bound, so evaluation past a rescale needs keys for the lower basis.

use std::collections::HashMap;

use rand::Rng;

use crate::PolyRing;
use crate::ckks::algorithms::sum_slots::{required_rotation_offsets, sum_slots};
use crate::ckks::engine::{CkksEngine, CkksParams};
use crate::ckks::keys::{RelinearizationKey, RotationKey, SecretKey};
use crate::ckks::types::{Ciphertext, Plaintext};

/// Rotation keys indexed by (positive) slot offset.
pub type RotationKeyMap<P, const DEGREE: usize> =
    HashMap<i32, RotationKey<P, DEGREE>>;

/// Generates the rotation key set needed to sum `slots` packed slots.
///
/// Convenience wrapper so callers do not hand-roll the offset list.
#[must_use]
pub fn generate_rotation_key_set<P, R, const DEGREE: usize>(
    engine: &CkksEngine<P, DEGREE>,
    sk: &SecretKey<P, DEGREE>,
    slots: usize,
    rng: &mut R,
) -> RotationKeyMap<P, DEGREE>
where
    P: PolyRing<DEGREE>,
    R: Rng,
{
    required_rotation_offsets(slots)
        .into_iter()
        .map(|offset| (offset, engine.generate_rotation_key(sk, offset, rng)))
        .collect()
}

/// Ciphertext evaluator that hides the per-level key juggling.
///
/// A gadget key (relin or rotation) only decrypts correctly at the modulus
/// level it was generated for, so every ct x ct multiplication after a
/// rescale needs a fresh key for the reduced basis. This struct lazily
/// generates and caches those keys per level, driven by the master secret
/// key (which mod-switches exactly to any lower level).
pub struct LeveledEvaluator<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    params: CkksParams<DEGREE>,
    sk: SecretKey<P, DEGREE>,
    /// Gadget relin keys indexed by level.
    rlks: HashMap<usize, RelinearizationKey<P, DEGREE>>,
    /// Gadget rotation keys indexed by (level, rotation offset).
    rotkeys: HashMap<(usize, i32), RotationKey<P, DEGREE>>,
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> LeveledEvaluator<P, DEGREE> {
    /// # Panics
    /// Panics later (lazily, at key generation) if `sk` does not decrypt the
    /// evaluated ciphertexts.
    pub fn new(params: CkksParams<DEGREE>, sk: SecretKey<P, DEGREE>) -> Self {
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
        a: &Ciphertext<P, DEGREE>,
        b: &Ciphertext<P, DEGREE>,
        rng: &mut R,
    ) -> Ciphertext<P, DEGREE> {
        let level = a.level();
        let rlk = self.rlks.remove(&level).unwrap_or_else(|| {
            RelinearizationKey::generate(
                &self.secret_key_at(level),
                self.key_error_std(),
                rng,
            )
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
        ct: &Ciphertext<P, DEGREE>,
        slots: usize,
        rng: &mut R,
    ) -> Ciphertext<P, DEGREE> {
        let level = ct.level();
        let mut rotkeys = RotationKeyMap::new();
        for offset in required_rotation_offsets(slots) {
            let key = self.rotkeys.remove(&(level, offset)).unwrap_or_else(|| {
                RotationKey::generate(
                    &self.secret_key_at(level),
                    offset,
                    self.key_error_std(),
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
    pub fn decrypt(&self, ct: &Ciphertext<P, DEGREE>) -> Plaintext<P, DEGREE> {
        CkksEngine::<P, DEGREE>::decrypt(ct, &self.secret_key_at(ct.level()))
    }

    /// Key error std-dev, as [`CkksEngine::generate_relin_key`] uses.
    fn key_error_std(&self) -> f64 {
        self.params.error_variance.sqrt()
    }

    /// The master secret key mod-switched down to `level`.
    fn secret_key_at(&self, level: usize) -> SecretKey<P, DEGREE> {
        self.sk
            .mod_switch_to(level)
            .expect("LeveledEvaluator: ciphertext above the secret key level")
    }
}
