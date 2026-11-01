//! Key switching for complex conjugation.
use rand::Rng;

use crate::PolyRing;
use crate::ckks::keys::{KeySwitchingKey, SecretKey};

/// Evaluation key for the automorphism X -> X^(2*N-1).
pub struct ConjugationKey<P: PolyRing<N>, const N: usize> {
    pub(crate) ksk: KeySwitchingKey<P, N>,
}

impl<P: PolyRing<N>, const N: usize> ConjugationKey<P, N> {
    /// Generate a level-specific conjugation key during client setup.
    pub fn generate<R: Rng>(
        secret: &SecretKey<P, N>,
        error_std: f64,
        rng: &mut R,
    ) -> Self {
        Self {
            ksk: KeySwitchingKey::generate(
                secret,
                &secret.poly.automorphism((2 * N - 1) as u64),
                error_std,
                rng,
            ),
        }
    }
}
