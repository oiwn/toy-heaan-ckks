use super::builder::CkksEngineBuilder;
use super::types::{Ciphertext, Plaintext};
use crate::PolyRing;
use crate::ckks::keys::{
    PublicKey, PublicKeyError, PublicKeyParams, RelinearizationKey, RotationKey,
    SecretKey, SecretKeyError, SecretKeyParams,
};
use rand::Rng;

pub struct CkksEngine<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
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
    P: PolyRing<DEGREE>,
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

    /// Generates a relinearization key at `sk`'s level (see
    /// [`RelinearizationKey`]); key errors use `std_dev = sqrt(error_variance)`.
    pub fn generate_relin_key<R: Rng>(
        &self,
        sk: &SecretKey<P, DEGREE>,
        rng: &mut R,
    ) -> RelinearizationKey<P, DEGREE> {
        RelinearizationKey::generate(sk, self.params.error_variance.sqrt(), rng)
    }

    /// Generates a rotation key for slot offset `rotation` at `sk`'s level
    /// (see [`RotationKey`]). Negative `rotation` rotates slots right.
    pub fn generate_rotation_key<R: Rng>(
        &self,
        sk: &SecretKey<P, DEGREE>,
        rotation: i32,
        rng: &mut R,
    ) -> RotationKey<P, DEGREE> {
        RotationKey::generate(sk, rotation, self.params.error_variance.sqrt(), rng)
    }

    // Encryption/Decryption
    pub fn encrypt<R: Rng>(
        &self,
        plaintext: &Plaintext<P, DEGREE>,
        public_key: &PublicKey<P, DEGREE>,
        rng: &mut R,
    ) -> Ciphertext<P, DEGREE> {
        let u = P::sample_ternary(self.params.hamming_weight, &self.context, rng);
        let e0 = P::sample_gaussian(self.params.error_variance, &self.context, rng);
        let e1 = P::sample_gaussian(self.params.error_variance, &self.context, rng);

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
            logp: plaintext.scale_bits, // Precision from encoding
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
            scale_bits: ciphertext.logp, // Use precision parameter
            slots: DEGREE / 2, // After decryption, assume max slots (decoder will handle actual slots)
        }
    }
}
