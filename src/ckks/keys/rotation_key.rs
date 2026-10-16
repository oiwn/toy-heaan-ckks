use rand::Rng;

use crate::PolyRing;
use crate::ckks::keys::{KeySwitchingKey, SecretKey};

/// Rotation key for one slot offset.
///
/// After the automorphism `X -> X^(5^k)` a ciphertext decrypts under the
/// rotated key `s_k = s(X^(5^k))`; key switching brings it back to `s`.
/// Bound to the level it was generated at.
#[derive(Clone, Debug)]
pub struct RotationKey<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    /// Key switching key from `s_k` to `s`.
    pub ksk: KeySwitchingKey<P, DEGREE>,
    /// The slot rotation offset this key was generated for.
    pub rotation: i32,
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> RotationKey<P, DEGREE> {
    /// Generates the key for `rotation` at `sk`'s level with key errors of
    /// `error_std`. A negative `rotation` rotates slots right.
    pub fn generate<R: Rng>(
        sk: &SecretKey<P, DEGREE>,
        rotation: i32,
        error_std: f64,
        rng: &mut R,
    ) -> Self {
        let s_k = sk.poly.rotate_slots(rotation);
        Self {
            ksk: KeySwitchingKey::generate(sk, &s_k, error_std, rng),
            rotation,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use super::*;
    use crate::ckks::keys::{DIGIT_BITS, SecretKeyParams};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

    #[test]
    fn first_term_encrypts_the_rotated_key() {
        const N: usize = 16;
        let basis =
            Arc::new(RnsBasis::new(generate_primes(30, 2, N as u64)).unwrap());
        let mut rng = ChaCha20Rng::seed_from_u64(13);
        let params = SecretKeyParams::new(N / 2).unwrap();
        let sk = SecretKey::<RnsPoly<N>, N>::generate(&params, &basis, &mut rng)
            .unwrap();

        for k in [1, 2, -1] {
            let rotk = RotationKey::generate(&sk, k, 3.2, &mut rng);
            assert_eq!(rotk.rotation, k);
            // b[0] + a[0]*s - s_k.gadget_term(0) is the key error e_0.
            let mut residual = rotk.ksk.a[0].clone();
            residual *= &sk.poly;
            residual += &rotk.ksk.b[0];
            residual -= &sk.poly.rotate_slots(k).gadget_term(0, DIGIT_BITS);
            let max = residual.to_coeffs().iter().map(|c| c.abs()).max();
            assert!(max.unwrap() < 30, "rotation {k}: error {max:?}");
        }
    }
}
