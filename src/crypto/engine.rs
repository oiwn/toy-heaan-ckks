use super::builder::CkksEngineBuilder;
use super::types::{Ciphertext, Plaintext};
use crate::{
    PolyRing, PolySampler, PublicKey, PublicKeyError, PublicKeyParams,
    RelinearizationKey, RelinearizationKeyError, RelinearizationKeyParams,
    SecretKey, SecretKeyError, SecretKeyParams,
};
use rand::Rng;

pub struct CkksEngine<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    context: P::Context,
    pub params: CkksParams<DEGREE>,
}

#[derive(Debug, Clone)]
pub struct CkksParams<const DEGREE: usize> {
    pub error_variance: f64,
    pub hamming_weight: usize,
    pub scale_bits: u32,
}

impl<P, const DEGREE: usize> CkksEngine<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    pub fn builder() -> CkksEngineBuilder<DEGREE> {
        CkksEngineBuilder::new()
    }

    pub fn new(context: P::Context, params: CkksParams<DEGREE>) -> Self {
        Self { context, params }
    }

    pub fn context(&self) -> &P::Context {
        &self.context
    }

    pub fn generate_secret_key<R: Rng>(
        &self,
        rng: &mut R,
    ) -> Result<SecretKey<P, DEGREE>, SecretKeyError> {
        let sk_params = SecretKeyParams::new(self.params.hamming_weight)?;
        SecretKey::generate(&sk_params, &self.context, rng)
    }

    pub fn generate_public_key<R: Rng>(
        &self,
        secret_key: &SecretKey<P, DEGREE>,
        rng: &mut R,
    ) -> Result<PublicKey<P, DEGREE>, PublicKeyError> {
        let pk_params = PublicKeyParams::new(3.2)?;
        PublicKey::generate(secret_key, &pk_params, &self.context, rng)
    }

    /// Generate a relinearization key for ciphertext multiplication
    ///
    /// The relinearization key enables converting degree-2 ciphertexts (produced by multiplication)
    /// back to degree-1 ciphertexts while preserving the encrypted plaintext.
    ///
    /// # Arguments
    /// * `secret_key` - The secret key used to generate the relinearization key
    /// * `rng` - Random number generator for sampling
    ///
    /// # Returns
    /// * `Result<RelinearizationKey<P, DEGREE>, RelinearizationKeyError>` - The generated key or error
    pub fn generate_relinearization_key<R: Rng>(
        &self,
        secret_key: &SecretKey<P, DEGREE>,
        rng: &mut R,
    ) -> Result<RelinearizationKey<P, DEGREE>, RelinearizationKeyError> {
        // Use same error variance as for other key generation
        let relin_params =
            RelinearizationKeyParams::new(self.params.error_variance.sqrt())?;
        RelinearizationKey::generate(secret_key, &relin_params, &self.context, rng)
    }

    // Encryption/Decryption
    pub fn encrypt<R: Rng>(
        &self,
        plaintext: &Plaintext<P, DEGREE>,
        public_key: &PublicKey<P, DEGREE>,
        rng: &mut R,
    ) -> Ciphertext<P, DEGREE> {
        let u = P::sample_tribits(rng, self.params.hamming_weight, &self.context);
        let e0 = P::sample_gaussian(rng, self.params.error_variance, &self.context);
        let e1 = P::sample_gaussian(rng, self.params.error_variance, &self.context);

        // c0 = b * u + e0 + m
        let mut c0 = public_key.b.clone();
        c0 *= &u;
        c0 += &e0;
        c0 += &plaintext.poly;

        // c1 = a * u + e1
        let mut c1 = public_key.a.clone();
        c1 *= &u;
        c1 += &e1;

        Ciphertext {
            c0,
            c1,
            scale: plaintext.scale,
        }
    }

    pub fn decrypt(
        ciphertext: &Ciphertext<P, DEGREE>,
        secret_key: &SecretKey<P, DEGREE>,
    ) -> Plaintext<P, DEGREE> {
        // m = c0 + c1 * s
        let mut result = ciphertext.c1.clone();
        result *= &secret_key.poly;
        result += &ciphertext.c0;

        Plaintext {
            poly: result,
            scale: ciphertext.scale,
        }
    }

    // Homomorphic operations
    pub fn add_ciphertexts(
        ct1: &Ciphertext<P, DEGREE>,
        ct2: &Ciphertext<P, DEGREE>,
    ) -> Ciphertext<P, DEGREE> {
        let mut c0 = ct1.c0.clone();
        c0 += &ct2.c0;

        let mut c1 = ct1.c1.clone();
        c1 += &ct2.c1;

        // Scale should be the same for both
        assert_eq!(ct1.scale, ct2.scale);

        Ciphertext {
            c0,
            c1,
            scale: ct1.scale,
        }
    }
}
