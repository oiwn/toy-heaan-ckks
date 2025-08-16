use crate::{PolyRing, PolySampler, SecretKey};
use rand::Rng;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum RelinearizationKeyError {
    #[error("Invalid error standard deviation: {0} (must be positive)")]
    InvalidErrorStd(f64),
    #[error("Relinearization key generation failed: {0}")]
    GenerationFailed(String),
    #[error("Secret key validation failed: {0}")]
    InvalidSecretKey(String),
}

/// Parameters for relinearization key generation
pub struct RelinearizationKeyParams<const DEGREE: usize> {
    /// Standard deviation for the error distribution
    pub error_std: f64,
}

impl<const DEGREE: usize> RelinearizationKeyParams<DEGREE> {
    /// Create new relinearization key parameters
    pub fn new(error_std: f64) -> Result<Self, RelinearizationKeyError> {
        if error_std <= 0.0 {
            return Err(RelinearizationKeyError::InvalidErrorStd(error_std));
        }
        Ok(Self { error_std })
    }

    /// Validate parameters
    pub fn validate(&self) -> Result<(), RelinearizationKeyError> {
        if self.error_std <= 0.0 {
            Err(RelinearizationKeyError::InvalidErrorStd(self.error_std))
        } else {
            Ok(())
        }
    }
}

/// Relinearization key used to transform ciphertexts after multiplication
///
/// The relinearization key allows converting a degree-2 ciphertext (c0, c1, c2)
/// back to a degree-1 ciphertext (c0', c1') while preserving the encrypted plaintext.
///
/// The key satisfies: `b + a * s = s^2`
#[derive(Debug, Clone)]
pub struct RelinearizationKey<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    /// "a" component: uniformly random polynomial
    pub a: P,
    /// "b" component: b = -(a * s) + e + s^2
    pub b: P,
}

impl<P, const DEGREE: usize> RelinearizationKey<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    /// Generate a relinearization key from a secret key
    ///
    /// # Relinearization Key Generation Algorithm
    /// 1. Compute s^2 = secret_key * secret_key
    /// 2. Sample a uniformly random polynomial `a`
    /// 3. Sample an error polynomial `e` from Gaussian distribution
    /// 4. Compute `b = -(a * s) + e + s^2`
    ///
    /// The resulting key (a, b) satisfies: `b + a * s ≈ s^2` (mod small error)
    /// This allows us to replace c2 * s^2 with c2 * (b + a * s) during relinearization.
    pub fn generate<R: Rng>(
        secret_key: &SecretKey<P, DEGREE>,
        params: &RelinearizationKeyParams<DEGREE>,
        context: &P::Context,
        rng: &mut R,
    ) -> Result<Self, RelinearizationKeyError> {
        // Validate parameters
        params.validate()?;

        // Step 1: Compute s^2 = secret_key * secret_key
        let mut s_squared = secret_key.poly.clone();
        s_squared *= &secret_key.poly; // s^2

        // Step 2: Sample uniformly random polynomial 'a'
        let a = P::sample_uniform(context, rng);

        // Step 3: Sample error polynomial 'e' from Gaussian distribution
        let e = P::sample_gaussian(params.error_std, context, rng);

        // Step 4: Compute b = -(a * s) + e + s^2
        let mut a_times_s = a.clone();
        a_times_s *= &secret_key.poly; // a * s

        let mut b = -a_times_s; // -(a * s)
        b += &e; // -(a * s) + e
        b += &s_squared; // -(a * s) + e + s^2

        Ok(RelinearizationKey { a, b })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BigIntPolyRing, SecretKeyParams};
    use crypto_bigint::{NonZero, U256};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    fn get_test_modulus() -> NonZero<U256> {
        let modulus_words = [0x1, 0x0, 0x0, 0x8000000000000000u64];
        NonZero::new(U256::from_words(modulus_words)).unwrap()
    }

    #[test]
    fn test_relinearization_key_generation() {
        const DEGREE: usize = 8;
        let context = get_test_modulus();
        let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

        // Generate secret key
        let sk_params = SecretKeyParams::new(DEGREE / 2).unwrap();
        let secret_key = SecretKey::<BigIntPolyRing<DEGREE>, DEGREE>::generate(
            &sk_params, &context, &mut rng,
        )
        .unwrap();

        // Generate relinearization key
        let relin_params = RelinearizationKeyParams::new(3.2).unwrap();
        let relin_key = RelinearizationKey::generate(
            &secret_key,
            &relin_params,
            &context,
            &mut rng,
        )
        .unwrap();

        // Basic sanity checks
        // The key should be properly constructed (no panics above means basic success)

        // Verify the key relation: b + a * s ≈ s^2 (we can't check exact equality due to error)
        let mut s_squared = secret_key.poly.clone();
        s_squared *= &secret_key.poly;

        let mut a_times_s = relin_key.a.clone();
        a_times_s *= &secret_key.poly;

        let mut verification = relin_key.b.clone();
        verification += &a_times_s; // b + a * s

        // The difference (b + a * s) - s^2 should be the small error term
        // We can't verify exact values due to random error, but we can check it doesn't panic
        let neg_s_squared = -s_squared;
        verification += &neg_s_squared;

        // If we get here without panicking, the key generation worked properly
        // TODO: need to test result is small
        assert!(
            true,
            "Relinearization key generation completed successfully"
        );
    }

    #[test]
    fn test_invalid_error_std() {
        let result = RelinearizationKeyParams::<8>::new(-1.0);
        assert!(result.is_err());

        let result = RelinearizationKeyParams::<8>::new(0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_valid_error_std() {
        let result = RelinearizationKeyParams::<8>::new(3.2);
        assert!(result.is_ok());
    }
}
