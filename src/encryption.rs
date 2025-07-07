//! CKKS encryption and decryption operations.
//!
//! This module contains the core cryptographic operations for the CKKS homomorphic
//! encryption scheme, including functions to encrypt plaintexts and decrypt ciphertexts.
use crate::{
    Ciphertext, PublicKey, RnsPolyRing, SecretKey, generate_ternary_poly,
    sample_gaussian_poly,
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
        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![7829, 6761, 5693])
                .build()
                .unwrap(),
        );

        let mut rng = ChaCha20Rng::seed_from_u64(123);

        // Generate keys
        let sk_params = SecretKeyParams {
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

        // Create test plaintext
        let plaintext_coeffs = [1i64, 2, 3, 4, 5, 6, 7, 8];
        let plaintext: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_integer_coeffs(&plaintext_coeffs, basis.clone());
        let plaintext_coeffs_from_rns = plaintext.to_u64_coefficients();

        println!(
            "Plaintext coeffs (CRT reconstructed): {:?}",
            plaintext_coeffs_from_rns
        );

        // Encrypt
        let ciphertext = encrypt(&plaintext, &public_key, 1024.0, &mut rng);

        // Decrypt
        let decrypted = decrypt(&ciphertext, &secret_key);

        // Check if decryption is close to original (allowing for noise)
        let original_coeffs = plaintext.to_u64_coefficients();
        let decrypted_coeffs = decrypted.to_u64_coefficients();

        println!("Original:  {:?}", original_coeffs);
        println!("Decrypted: {:?}", decrypted_coeffs);

        // The coefficients should be close (exact match in this small example)
        // In practice, there might be small noise, but for our test parameters
        // and small values, they should match exactly
        for (orig, dec) in original_coeffs.iter().zip(decrypted_coeffs.iter()) {
            let diff = if orig > dec { orig - dec } else { dec - orig };
            assert!(
                diff <= 5, // Allow small error due to noise
                "Decryption error too large: {} vs {}",
                orig,
                dec
            );
        }
    }

    #[test]
    fn test_ciphertext_scale_preservation() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();
        let mut rng = ChaCha20Rng::seed_from_u64(456);

        // Generate keys
        let sk_params = SecretKeyParams {
            basis: basis.clone(),
            hamming_weight: 2,
        };
        let secret_key = SecretKey::generate(&sk_params, &mut rng).unwrap();

        let pk_params = PublicKeyParams {
            basis: basis.clone(),
            error_std: 3.0,
        };
        let public_key =
            PublicKey::generate(&secret_key, &pk_params, &mut rng).unwrap();

        // Create test plaintext
        let plaintext_coeffs = [10i64, 20, 30, 40];
        let plaintext: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_integer_coeffs(&plaintext_coeffs, basis.clone());

        let test_scale = 2048.0;

        // Encrypt with specific scale
        let ciphertext = encrypt(&plaintext, &public_key, test_scale, &mut rng);

        // Verify scale is preserved
        assert_eq!(ciphertext.scale, test_scale);
    }
}
