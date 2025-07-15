use crate::{PolyRing, PolySampler, SecretKey};
use rand::Rng;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PublicKeyError {
    #[error("Invalid error standard deviation: {0} (must be positive)")]
    InvalidErrorStd(f64),
    #[error("Public key generation failed: {0}")]
    GenerationFailed(String),
    #[error("Secret key validation failed: {0}")]
    InvalidSecretKey(String),
    #[error("Public key parameter error: {0}")]
    InvalidParams(String),
}

/// Parameters for public key generation
pub struct PublicKeyParams<const DEGREE: usize> {
    /// Standard deviation for the error distribution
    pub error_std: f64,
}

impl<const DEGREE: usize> PublicKeyParams<DEGREE> {
    /// Create new public key parameters
    pub fn new(error_std: f64) -> Result<Self, PublicKeyError> {
        if error_std <= 0.0 {
            return Err(PublicKeyError::InvalidErrorStd(error_std));
        }
        Ok(Self { error_std })
    }

    /// Validate parameters
    pub fn validate(&self) -> Result<(), PublicKeyError> {
        if self.error_std <= 0.0 {
            Err(PublicKeyError::InvalidErrorStd(self.error_std))
        } else {
            Ok(())
        }
    }
}

/// Generic public key for CKKS scheme
#[derive(Debug, Clone)]
pub struct PublicKey<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    /// "a" component: uniformly random polynomial
    pub a: P,
    /// "b" component: b = -(a * s) + e
    pub b: P,
}

impl<P, const DEGREE: usize> PublicKey<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    /// Generate a new RLWE public key from a secret key with context
    ///
    /// # RLWE Public Key Generation
    /// 1. Sample a uniformly random polynomial `a`
    /// 2. Sample an error polynomial `e` from Gaussian distribution  
    /// 3. Compute `b = -(a * s) + e` where `s` is the secret key
    ///
    /// The public key (a, b) satisfies: `b + a * s ≈ 0` (mod small error)
    pub fn generate<R: Rng>(
        secret_key: &SecretKey<P, DEGREE>,
        params: &PublicKeyParams<DEGREE>,
        context: &P::Context,
        rng: &mut R,
    ) -> Result<Self, PublicKeyError> {
        // Validate parameters
        params.validate()?;

        // Sample uniformly random polynomial 'a'
        let a = P::sample_uniform(rng, context);

        // Sample error polynomial 'e' from Gaussian distribution
        let e = P::sample_gaussian(rng, params.error_std, context);

        // Compute b = -(a * s) + e
        let mut a_times_s = a.clone();
        a_times_s *= &secret_key.poly; // a * s

        let mut b = -a_times_s; // -(a * s)
        b += &e; // b = -(a * s) + e

        Ok(PublicKey { a, b })
    }
}
