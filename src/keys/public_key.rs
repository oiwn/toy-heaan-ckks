use crate::common::{sample_gaussian_poly, sample_uniform_poly};
use crate::keys::SecretKey;
use crate::rings::{RnsBasis, RnsPolyRing};
use rand::Rng;
use std::sync::Arc;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PublicKeyError {
    #[error("Secret key validation failed: {0}")]
    InvalidSecretKey(String),
    #[error("Public key parameter error: {0}")]
    InvalidParams(String),
}

/// Parameters for generating an RNS-encoded public key.
pub struct PublicKeyParams<const DEGREE: usize> {
    pub basis: Arc<RnsBasis>,
    /// Standard deviation for the error distribution
    pub error_std: f64,
}

impl<const DEGREE: usize> PublicKeyParams<DEGREE> {
    fn validate(&self) -> Result<(), PublicKeyError> {
        if !(self.error_std > 0.0) {
            Err(PublicKeyError::InvalidParams(
                "error_std must be positive".into(),
            ))
        } else {
            Ok(())
        }
    }
}

/// RNS-encoded public key (RLWE sample).
pub struct PublicKey<const DEGREE: usize> {
    /// "b" component: b = -(a * s) + e
    pub b: RnsPolyRing<DEGREE>,
    /// "a" component: uniformly random
    pub a: RnsPolyRing<DEGREE>,
}

impl<const DEGREE: usize> PublicKey<DEGREE> {
    /// Generate a new RLWE public key under RNS.
    pub fn generate<R: Rng>(
        secret_key: &SecretKey<DEGREE>,
        params: &PublicKeyParams<DEGREE>,
        rng: &mut R,
    ) -> Result<Self, PublicKeyError> {
        // Validate parameters
        params.validate()?;

        // Draw a uniformly random polynomial 'a'
        let a: RnsPolyRing<DEGREE> = sample_uniform_poly(params.basis.clone(), rng);

        // Draw error polynomial 'e' from Gaussian
        let e: RnsPolyRing<DEGREE> =
            sample_gaussian_poly(params.basis.clone(), params.error_std, rng);

        // Compute b = -(a * s) + e = -a * s + e
        let mut a_times_s = a.clone();
        a_times_s *= &secret_key.s; // a * s

        let neg_a_times_s = -&a_times_s; // -(a * s)

        let mut b = neg_a_times_s;
        b += &e; // b = -(a * s) + e

        Ok(PublicKey { b, a })
    }
}
