//! CKKS encryption and decryption operations.
//!
//! This module contains the core cryptographic operations for the CKKS homomorphic
//! encryption scheme, including functions to encrypt plaintexts and decrypt ciphertexts.
use crate::{
    Ciphertext, EncodingParams, PublicKey, RnsPolyRing, SecretKey, decode, encode,
    generate_ternary_poly, sample_gaussian_poly,
};
use rand::Rng;

/// Encrypts a plaintext polynomial using a public key.
///
/// # Encryption Process
/// 1. Generate small random error polynomials e1 and e2 (Gaussian)
/// 2. Generate a random ephemeral value u (small ternary polynomial)
/// 3. Compute c0 = b*u + e1 + plaintext
/// 4. Compute c1 = a*u + e2
///
/// # Arguments
/// * `plaintext` - The plaintext polynomial to encrypt (RnsPolyRing)
/// * `public_key` - The public key used for encryption
/// * `scale` - Scaling factor used in the CKKS scheme
/// * `rng` - Random number generator for encryption randomness
///
/// # Returns
/// A new ciphertext containing the encrypted plaintext
pub fn encrypt<const DEGREE: usize, R: Rng>(
    plaintext: &RnsPolyRing<DEGREE>,
    public_key: &PublicKey<DEGREE>,
    scale: f64,
    rng: &mut R,
) -> Ciphertext<DEGREE> {
    let basis = plaintext.basis.clone();

    // Generate small random error polynomials (Gaussian noise)
    let e1: RnsPolyRing<DEGREE> = sample_gaussian_poly(basis.clone(), 3.0, rng);
    let e2: RnsPolyRing<DEGREE> = sample_gaussian_poly(basis.clone(), 3.0, rng);

    // Generate random ephemeral ternary value u
    let u: RnsPolyRing<DEGREE> = generate_ternary_poly(basis.clone(), 0.5, rng);

    // c0 = b*u + e1 + plaintext
    let mut c0 = public_key.b.clone();
    c0 *= &u; // b * u
    c0 += &e1; // + e1
    c0 += plaintext; // + plaintext

    // c1 = a*u + e2
    let mut c1 = public_key.a.clone();
    c1 *= &u; // a * u
    c1 += &e2; // + e2

    Ciphertext::new(c0, c1, scale)
}

/// Decrypts a ciphertext using a secret key.
///
/// # Decryption Process
/// Computes c0 + c1*s, where:
/// * c0, c1 are components of the ciphertext
/// * s is the secret polynomial
///
/// # Arguments
/// * `ciphertext` - The ciphertext to decrypt
/// * `secret_key` - The secret key used for decryption
///
/// # Returns
/// The decrypted plaintext polynomial
pub fn decrypt<const DEGREE: usize>(
    ciphertext: &Ciphertext<DEGREE>,
    secret_key: &SecretKey<DEGREE>,
) -> RnsPolyRing<DEGREE> {
    // Compute c1*s
    let mut c1_times_s = ciphertext.c1.clone();
    c1_times_s *= &secret_key.s;

    // Compute c0 + c1*s
    let result = &ciphertext.c0 + &c1_times_s;

    result
}

#[cfg(test)]
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
        let enc_params = EncodingParams::new(DEGREE, scale_bits).unwrap();
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
}
