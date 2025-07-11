//! CKKS encryption and decryption operations.
//!
//! This module contains the core cryptographic operations for the CKKS homomorphic
//! encryption scheme, including functions to encrypt plaintexts and decrypt ciphertexts.
use crate::{
    Ciphertext, EncodingParams, Plaintext, PolyRing, PolySampler, PublicKey,
    SecretKey, decode, encode,
};
use rand::Rng;
use thiserror::Error;

/// Errors that can occur during encryption/decryption operations
#[derive(Error, Debug)]
pub enum EncryptionError {
    #[error("Incompatible polynomial contexts")]
    IncompatibleContexts,

    #[error("Scale mismatch: expected {expected}, got {actual}")]
    ScaleMismatch { expected: f64, actual: f64 },

    #[error("Polynomial sampling failed: {message}")]
    SamplingError { message: String },
}

pub type EncryptionResult<T> = Result<T, EncryptionError>;

/// Encrypts a plaintext using the CKKS encryption scheme.
///
/// # CKKS Encryption Process
/// 1. Sample a small random polynomial `u` (ternary with small hamming weight)
/// 2. Sample error polynomials `e0` and `e1` from Gaussian distribution
/// 3. Compute:
///    - `c0 = pk.b * u + e0 + plaintext.poly`  
///    - `c1 = pk.a * u + e1`
///
/// The resulting ciphertext satisfies: `c0 + c1 * sk ≈ plaintext` (mod noise)
///
/// # Arguments
/// * `plaintext` - The plaintext to encrypt
/// * `public_key` - Public key for encryption
/// * `rng` - Random number generator for sampling randomness
///
/// # Returns
/// * `Ok(Ciphertext)` - The encrypted ciphertext
/// * `Err(EncryptionError)` - If encryption fails
pub fn encrypt<P, const DEGREE: usize, R: Rng>(
    plaintext: &Plaintext<P, DEGREE>,
    public_key: &PublicKey<P, DEGREE>,
    rng: &mut R,
) -> EncryptionResult<Ciphertext<P, DEGREE>>
where
    P: PolyRing<DEGREE> + PolySampler<DEGREE>,
{
    let sampler = P::zero();
    // Sample ephemeral ternary polynomial with small hamming weight
    let u = sampler.sample_tribits(rng, 2); // Small hamming weight for security

    // Sample error polynomials from Gaussian distribution
    let e0 = sampler.sample_gaussian(rng, 3.0); // Standard deviation = 3.0
    let e1 = sampler.sample_gaussian(rng, 3.0);

    // Compute c0 = pk.b * u + e0 + plaintext
    let mut c0 = public_key.b.clone();
    // c0.mul_assign(&u); // pk.b * u
    c0 *= &u;
    c0.add_assign(&e0); // + e0
    c0.add_assign(&plaintext.poly); // + plaintext

    // Compute c1 = pk.a * u + e1
    let mut c1 = public_key.a.clone();
    c1.mul_assign(&u); // pk.a * u
    c1.add_assign(&e1); // + e1

    Ok(Ciphertext {
        c0,
        c1,
        scale: plaintext.scale,
    })
}

/// Decrypts a ciphertext using the CKKS decryption scheme.
///
/// # CKKS Decryption Process
/// Computes: `plaintext ≈ c0 + c1 * secret_key` (mod noise)
///
/// This works because during encryption we have:
/// - `c0 = pk.b * u + e0 + plaintext`
/// - `c1 = pk.a * u + e1`  
/// - `pk.b = -(pk.a * sk + e_pk)`
///
/// So: `c0 + c1 * sk = plaintext + error_terms`
///
/// # Arguments
/// * `ciphertext` - The ciphertext to decrypt
/// * `secret_key` - Secret key for decryption
///
/// # Returns
/// * `Plaintext` - The decrypted plaintext (with noise)
pub fn decrypt<P, const DEGREE: usize>(
    ciphertext: &Ciphertext<P, DEGREE>,
    secret_key: &SecretKey<P, DEGREE>,
) -> Plaintext<P, DEGREE>
where
    P: PolyRing<DEGREE>,
{
    // Compute c1 * secret_key
    let mut c1_times_sk = ciphertext.c1.clone();
    c1_times_sk.mul_assign(&secret_key.poly);

    // Compute result = c0 + c1 * secret_key
    let mut result = ciphertext.c0.clone();
    result.add_assign(&c1_times_sk);

    Plaintext {
        poly: result,
        scale: ciphertext.scale,
    }
}

/* #[cfg(test)]
mod tests {
    use super::*;
    use crate::{PublicKeyParams, RnsBasisBuilder, SecretKeyParams};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use std::sync::Arc;

    fn create_test_basis<const DEGREE: usize>() -> Arc<crate::RnsBasis> {
        Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19, 23])
                .build()
                .unwrap(),
        )
    }

    #[test]
    fn test_encrypt_decrypt_roundtrip() {
        const DEGREE: usize = 8;
        let scale_bits = 10; // 1024
        let scale = (1u64 << scale_bits) as f64;

        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![7829, 6761, 5693])
                .build()
                .unwrap(),
        );

        let mut rng = ChaCha20Rng::seed_from_u64(123);

        // Generate keys
        let sk_params: SecretKeyParams<DEGREE> = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 4,
        };
        let secret_key = SecretKey::generate(&sk_params, &mut rng).unwrap();

        let pk_params = PublicKeyParams {
            basis: basis.clone(),
            error_std: 3.0,
        };
        let public_key =
            PublicKey::generate(&secret_key, &pk_params, &mut rng).unwrap();

        // Encode input values [1, 2, 3, 4]
        let input = vec![1.0, 2.0, 3.0, 4.0];
        let enc_params = EncodingParams::<DEGREE>::new(scale_bits).unwrap();
        let encoded_coeffs = encode(&input, &enc_params).unwrap();
        let plaintext = RnsPolyRing::from_i64_slice(&encoded_coeffs, basis.clone());

        // Encrypt
        let ciphertext = encrypt(&plaintext, &public_key, scale, &mut rng);

        // Decrypt
        let decrypted = decrypt(&ciphertext, &secret_key);
        let coeffs = decrypted.to_i64_coefficients();

        // Decode
        let decoded = decode(&coeffs, &enc_params).unwrap();

        println!("Original input: {:?}", input);
        println!("Decoded output: {:?}", decoded);

        for (i, (&expected, &actual)) in input.iter().zip(&decoded).enumerate() {
            assert!(
                (expected - actual).abs() < 0.1,
                "Mismatch at {}: expected {}, got {}",
                i,
                expected,
                actual
            );
        }
    }
} */
