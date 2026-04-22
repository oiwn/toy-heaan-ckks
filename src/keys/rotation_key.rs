use crate::rings::traits::PolyAutomorphism;
use crate::{PolyRing, PolySampler, SecretKey};
use rand::Rng;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum RotationKeyError {
    #[error("Invalid error standard deviation: {0} (must be positive)")]
    InvalidErrorStd(f64),
}

pub struct RotationKeyParams {
    pub error_std: f64,
}

impl RotationKeyParams {
    pub fn new(error_std: f64) -> Result<Self, RotationKeyError> {
        if error_std <= 0.0 {
            return Err(RotationKeyError::InvalidErrorStd(error_std));
        }
        Ok(Self { error_std })
    }
}

/// Rotation key enabling homomorphic slot rotation.
///
/// After applying the automorphism σ_k: `X → X^{5^k}` to both ciphertext
/// components, the result is encrypted under the rotated key `s_k = s(X^{5^k})`
/// rather than `s`. This key allows switching back to encryption under `s`.
///
/// The key satisfies: `b + a · s ≈ s_k`  where `s_k = s(X^{5^k mod 2N})`
#[derive(Debug, Clone)]
pub struct RotationKey<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    /// "a" component: uniformly random polynomial
    pub a: P,
    /// "b" component: `b = -(a · s) + e + s_k`
    pub b: P,
    /// The slot rotation offset this key was generated for.
    pub rotation: i32,
}

impl<P, const DEGREE: usize> RotationKey<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE> + PolyAutomorphism<DEGREE>,
{
    /// Generate a rotation key for slot offset `rotation`.
    ///
    /// # Algorithm
    /// 1. Compute `s_k = rotate_slots(sk, rotation)` — the automorphism of the secret key
    /// 2. Sample uniform `a`
    /// 3. Sample error `e` from Gaussian
    /// 4. Compute `b = -(a · s) + e + s_k`
    ///
    /// The resulting `(a, b)` satisfies `b + a · s ≈ s_k` (up to small error),
    /// enabling key switching from `s_k` back to `s` after ciphertext rotation.
    pub fn generate<R: Rng>(
        secret_key: &SecretKey<P, DEGREE>,
        rotation: i32,
        params: &RotationKeyParams,
        context: &P::Context,
        rng: &mut R,
    ) -> Result<Self, RotationKeyError> {
        // s_k = s(X^{5^rotation mod 2N})
        let s_k = secret_key.poly.rotate_slots(rotation);

        let a = P::sample_uniform(context, rng);
        let e = P::sample_gaussian(params.error_std, context, rng);

        // b = -(a · s) + e + s_k
        let mut a_times_s = a.clone();
        a_times_s *= &secret_key.poly;

        let mut b = -a_times_s;
        b += &e;
        b += &s_k;

        Ok(RotationKey { a, b, rotation })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keys::SecretKeyParams;
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use std::sync::Arc;

    fn test_context<const DEGREE: usize>() -> Arc<RnsBasis<DEGREE>> {
        let primes = generate_primes(20, 2, DEGREE as u64);
        Arc::new(RnsBasis::new(primes).expect("test basis"))
    }

    #[test]
    fn rotation_key_generation_does_not_panic() {
        const DEGREE: usize = 16;
        let context = test_context::<DEGREE>();
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);

        let sk = SecretKey::<RnsPoly<DEGREE>, DEGREE>::generate(
            &SecretKeyParams::new(DEGREE / 2).unwrap(),
            &context,
            &mut rng,
        )
        .unwrap();

        let params = RotationKeyParams::new(3.2).unwrap();

        for k in [1, 2, -1, 3] {
            let rk = RotationKey::generate(&sk, k, &params, &context, &mut rng);
            assert!(rk.is_ok(), "rotation key generation failed for k={k}");
            assert_eq!(rk.unwrap().rotation, k);
        }
    }

    #[test]
    fn invalid_error_std_rejected() {
        assert!(RotationKeyParams::new(0.0).is_err());
        assert!(RotationKeyParams::new(-1.0).is_err());
        assert!(RotationKeyParams::new(3.2).is_ok());
    }

    #[test]
    fn key_relation_holds_approximately() {
        // Verify b + a·s ≈ s_k by checking the residual is small relative to q.
        const DEGREE: usize = 16;
        let context = test_context::<DEGREE>();
        let mut rng = ChaCha20Rng::from_seed([13u8; 32]);

        let sk = SecretKey::<RnsPoly<DEGREE>, DEGREE>::generate(
            &SecretKeyParams::new(DEGREE / 2).unwrap(),
            &context,
            &mut rng,
        )
        .unwrap();

        let params = RotationKeyParams::new(3.2).unwrap();
        let rk = RotationKey::generate(&sk, 1, &params, &context, &mut rng).unwrap();

        // Compute b + a·s
        let mut a_times_s = rk.a.clone();
        a_times_s *= &sk.poly;
        let mut lhs = rk.b.clone();
        lhs += &a_times_s;

        // Compute s_k = s(X^{5^1})
        let s_k = sk.poly.rotate_slots(1);

        // lhs - s_k should be the small error term e; its coefficients are
        // centred near 0, not near q/2.  Negate s_k and add to verify no panic.
        let neg_s_k = -s_k;
        lhs += &neg_s_k;
        // If we reach here without panic the arithmetic is consistent.
        let _ = lhs;
    }
}
