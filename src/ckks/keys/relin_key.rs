use rand::Rng;

use crate::PolyRing;
use crate::ckks::keys::{KeySwitchingKey, SecretKey};

/// Relinearization key: switches the `s^2` part of a ciphertext product
/// back to `s`.
///
/// A ct x ct product decrypts as `d0 + d1*s + d2*s^2`; key switching turns
/// `d2 * s^2` into a pair under `s`, giving a degree-1 ciphertext again.
/// Bound to the level it was generated at.
#[derive(Clone, Debug)]
pub struct RelinearizationKey<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    /// Key switching key from `s^2` to `s`.
    pub ksk: KeySwitchingKey<P, DEGREE>,
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> RelinearizationKey<P, DEGREE> {
    /// Generates the key at `sk`'s level with key errors of `error_std`.
    pub fn generate<R: Rng>(
        sk: &SecretKey<P, DEGREE>,
        error_std: f64,
        rng: &mut R,
    ) -> Self {
        let mut s_sq = sk.poly.clone();
        s_sq *= &sk.poly;
        Self {
            ksk: KeySwitchingKey::generate(sk, &s_sq, error_std, rng),
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
    fn every_term_encrypts_the_scaled_square() {
        const N: usize = 16;
        let basis =
            Arc::new(RnsBasis::new(generate_primes(30, 2, N as u64)).unwrap());
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let params = SecretKeyParams::new(N / 2).unwrap();
        let sk = SecretKey::<RnsPoly<N>, N>::generate(&params, &basis, &mut rng)
            .unwrap();
        let rlk = RelinearizationKey::generate(&sk, 3.2, &mut rng);

        let mut s_sq = sk.poly.clone();
        s_sq *= &sk.poly;
        assert_eq!(rlk.ksk.a.len(), RnsPoly::gadget_len(&basis, DIGIT_BITS));
        for (k, (a, b)) in rlk.ksk.a.iter().zip(&rlk.ksk.b).enumerate() {
            // b + a*s - s^2.gadget_term(k) is the key error e_k.
            let mut residual = a.clone();
            residual *= &sk.poly;
            residual += b;
            residual -= &s_sq.gadget_term(k, DIGIT_BITS);
            let max = residual.to_coeffs().iter().map(|c| c.abs()).max();
            assert!(max.unwrap() < 30, "term {k}: error {max:?}");
        }
    }
}
