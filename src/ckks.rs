use crate::{
    Ciphertext, EncodingParams, Plaintext, PolyRing, PolySampler, PublicKey,
    PublicKeyParams, SecretKey, SecretKeyParams, decode, encode,
};
use rand::Rng;
use std::marker::PhantomData;

pub struct CkksEngine<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    _phantom: PhantomData<P>,
}

impl<P, const DEGREE: usize> CkksEngine<P, DEGREE>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    // Key generation
    pub fn generate_secret_key<R: Rng>(
        params: &SecretKeyParams<DEGREE>,
        rng: &mut R,
    ) -> SecretKey<P, DEGREE> {
        let sampler = P::zero();
        let poly = sampler.sample_tribits(rng, params.hamming_weight);
        SecretKey { poly }
    }

    pub fn generate_public_key<R: Rng>(
        secret_key: &SecretKey<P, DEGREE>,
        params: &PublicKeyParams<DEGREE>,
        rng: &mut R,
    ) -> PublicKey<P, DEGREE> {
        let sampler = P::zero();
        let a = sampler.sample_uniform(rng, u64::MAX);
        let e = sampler.sample_gaussian(rng, params.error_std);
        // let a = P::sample_uniform(rng, u64::MAX);
        // let e = P::sample_gaussian(rng, params.error_std);

        // b = -(a * s + e)
        let mut b = a.clone();
        b *= &secret_key.poly;
        b += &e;
        b = -b;

        PublicKey { a, b }
    }

    // Encoding/Decoding
    pub fn encode(
        values: &[f64],
        params: &EncodingParams<DEGREE>,
    ) -> Plaintext<P, DEGREE> {
        // Abstract encoding logic
        let scale = params.delta();
        let coeffs = Self::fft_encode(values, scale);
        let poly = P::from_coeffs(&coeffs);
        Plaintext { poly, scale }
    }

    pub fn decode(
        plaintext: &Plaintext<P, DEGREE>,
        params: &EncodingParams<DEGREE>,
    ) -> Vec<f64> {
        let coeffs = plaintext.poly.to_coeffs();
        Self::fft_decode(&coeffs, plaintext.scale)
    }

    // Encryption/Decryption
    pub fn encrypt<R: Rng>(
        plaintext: &Plaintext<P, DEGREE>,
        public_key: &PublicKey<P, DEGREE>,
        rng: &mut R,
    ) -> Ciphertext<P, DEGREE> {
        let sampler = P::zero();
        let u = sampler.sample_tribits(rng, 2); // Small hamming weight
        let e0 = sampler.sample_gaussian(rng, 3.0);
        let e1 = sampler.sample_gaussian(rng, 3.0);

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

    fn fft_encode(values: &[f64], scale: f64) -> Vec<u64> {
        // recover scale_bits from scale = 2^scale_bits
        let bits = scale.log2().round() as u32;
        let params =
            EncodingParams::<DEGREE>::new(bits).expect("invalid encoding params");
        let int_coeffs =
            encode::<DEGREE>(values, &params).expect("fft encode failed");
        int_coeffs.iter().map(|&c| c as u64).collect()
    }

    fn fft_decode(coeffs: &[u64], scale: f64) -> Vec<f64> {
        let bits = scale.log2().round() as u32;
        let params =
            EncodingParams::<DEGREE>::new(bits).expect("invalid encoding params");
        let signed: Vec<i64> = coeffs.iter().map(|&c| c as i64).collect();
        decode::<DEGREE>(&signed, &params).expect("fft decode failed")
    }
}
