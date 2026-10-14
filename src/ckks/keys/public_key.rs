//! Public Key Generation for CKKS Scheme
//!
//! This module implements public key generation for the CKKS homomorphic encryption scheme.
//! Public keys in CKKS are generated using the Ring Learning With Errors (RLWE) problem,
//! which provides security for the encryption scheme.
//!
//! # RLWE Public Key Structure
//!
//! A public key consists of two polynomials `(a, b)` where:
//! - `a` is uniformly random over the polynomial ring
//! - `b = -(a × s) + e` where `s` is the secret key and `e` is small error
//!
//! The security relies on the difficulty of finding `s` given `(a, b)`.
//!
//! # Security Properties
//!
//! The public key satisfies the RLWE relation:
//! ```text
//! b + a × s ≈ e (mod q)
//! ```
//! where `e` is a small error polynomial sampled from a Gaussian distribution.
//!
//! # Example
//!
//! ```rust,ignore
//! # use rand::SeedableRng;
//! # use rand_chacha::ChaCha20Rng;
//! # use std::sync::Arc;
//! use toy_heaan_ckks::{
//!     PublicKey, PublicKeyParams, RnsNttPoly, SecretKey, SecretKeyParams,
//! };
//! use toy_heaan_ckks::rings::backends::rns::RnsBasisBuilder;
//!
//! const DEGREE: usize = 16;
//!
//! # fn example() -> Result<(), Box<dyn std::error::Error>> {
//! let basis = Arc::new(
//!     RnsBasisBuilder::new(DEGREE)
//!         .with_custom_primes(vec![97])
//!         .build()
//!         .unwrap(),
//! );
//! let mut rng = ChaCha20Rng::seed_from_u64(123);
//!
//! let sk_params = SecretKeyParams::<DEGREE>::new(DEGREE / 2)?;
//! let secret_key: SecretKey<RnsNttPoly<DEGREE>, DEGREE> =
//!     SecretKey::generate(&sk_params, &basis, &mut rng)?;
//!
//! let pk_params = PublicKeyParams::<DEGREE>::new(3.2)?;
//! let public_key: PublicKey<RnsNttPoly<DEGREE>, DEGREE> =
//!     PublicKey::generate(&secret_key, &pk_params, &basis, &mut rng)?;
//! # Ok(())
//! # }
//! ```
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
        // Sample uniformly random polynomial 'a'
        let a = P::sample_uniform(context, rng);

        // Sample error polynomial 'e' from Gaussian distribution
        let e = P::sample_gaussian(params.error_std, context, rng);

        // Compute b = -(a * s) + e
        let mut a_times_s = a.clone();
        a_times_s *= &secret_key.poly; // a * s

        let mut b = -a_times_s; // -(a * s)
        b += &e; // b = -(a * s) + e

        Ok(PublicKey { a, b })
    }
}
