use super::types::{Ciphertext, Plaintext};
use crate::{
    EncodingParams, PolyRing, PolySampler, PublicKey, PublicKeyError,
    PublicKeyParams, SecretKey, SecretKeyError, SecretKeyParams, decode, encode,
};
use rand::Rng;

pub struct CkksEngine<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    pub context: P::Context,
}

impl<P, const DEGREE: usize> CkksEngine<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    pub fn new(context: P::Context) -> Self {
        Self { context }
    }

    pub fn generate_secret_key<R: Rng>(
        params: &SecretKeyParams<DEGREE>,
        rng: &mut R,
    ) -> Result<SecretKey<P, DEGREE>, SecretKeyError> {
        SecretKey::generate(params, rng)
    }

    pub fn generate_public_key<R: Rng>(
        &self,
        secret_key: &SecretKey<P, DEGREE>,
        params: &PublicKeyParams<DEGREE>,
        rng: &mut R,
    ) -> Result<PublicKey<P, DEGREE>, PublicKeyError> {
        PublicKey::generate(secret_key, params, &self.context, rng)
    }

    // Encoding/Decoding
    pub fn encode(
        &self,
        values: &[f64],
        params: &EncodingParams<DEGREE>,
    ) -> Plaintext<P, DEGREE> {
        // Abstract encoding logic
        let scale = params.delta();
        let coeffs = Self::fft_encode(values, scale);
        let poly = P::from_coeffs(&coeffs, &self.context);
        Plaintext { poly, scale }
    }

    pub fn decode(
        &self,
        plaintext: &Plaintext<P, DEGREE>,
        _params: &EncodingParams<DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.to_coeffs();
        Self::fft_decode(&coeffs, plaintext.scale)
    }

    // Encryption/Decryption
    pub fn encrypt<R: Rng>(
        &self,
        plaintext: &Plaintext<P, DEGREE>,
        public_key: &PublicKey<P, DEGREE>,
        rng: &mut R,
    ) -> Ciphertext<P, DEGREE> {
        let u = P::sample_tribits(rng, DEGREE / 2);
        let e0 = P::sample_gaussian(rng, 3.0);
        let e1 = P::sample_gaussian(rng, 3.0);

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

    pub fn add_plaintext(
        ciphertext: &Ciphertext<P, DEGREE>,
        plaintext: &Plaintext<P, DEGREE>,
    ) -> Ciphertext<P, DEGREE> {
        let mut c0 = ciphertext.c0.clone();
        c0.add_assign(&plaintext.poly);

        Ciphertext {
            c0,
            c1: ciphertext.c1.clone(),
            scale: ciphertext.scale,
        }
    }

    fn fft_encode(values: &[f64], scale: f64) -> Vec<i64> {
        // recover scale_bits from scale = 2^scale_bits
        let bits = scale.log2().round() as u32;
        let params =
            EncodingParams::<DEGREE>::new(bits).expect("invalid encoding params");
        encode::<DEGREE>(values, &params)
            .expect("fft encode failed")
            .into()
    }

    fn fft_decode(coeffs: &[i64], scale: f64) -> Vec<f64> {
        let bits = scale.log2().round() as u32;
        let params =
            EncodingParams::<DEGREE>::new(bits).expect("invalid encoding params");
        decode::<DEGREE>(coeffs, &params).expect("fft decode failed")
    }
}
