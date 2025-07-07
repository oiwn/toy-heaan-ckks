//! CKKS encryption and decryption operations.
//!
//! This module contains the core cryptographic operations for the CKKS homomorphic
//! encryption scheme, including functions to encrypt plaintexts and decrypt ciphertexts.
use crate::{
    Ciphertext, PublicKey, RnsPolyRing, SecretKey, i64_to_rns, sample_gaussian_poly,
};
use rand::Rng;
use std::sync::Arc;

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
    let u: RnsPolyRing<DEGREE> = generate_ternary_rns(basis.clone(), 0.5, rng);

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

/// Generate a ternary polynomial in RNS representation.
///
/// Creates a polynomial with coefficients in {-1, 0, 1} with given probability
/// for non-zero coefficients, then converts to RNS representation.
///
/// # Arguments
/// * `basis` - The RNS basis to use
/// * `prob` - Probability that a coefficient is non-zero (typically 0.5)
/// * `rng` - Random number generator
///
/// # Returns
/// A ternary polynomial in RNS representation
pub fn generate_ternary_rns<const DEGREE: usize, R: Rng>(
    basis: Arc<crate::RnsBasis>,
    prob: f64,
    rng: &mut R,
) -> RnsPolyRing<DEGREE> {
    // Generate ternary coefficients as signed integers
    let mut ternary_coeffs = [0i64; DEGREE];

    for coeff in &mut ternary_coeffs {
        if rng.random::<f64>() < prob {
            *coeff = if rng.random::<bool>() { 1 } else { -1 };
        }
        // else remains 0
    }

    // Convert to RNS representation
    i64_to_rns(&ternary_coeffs, basis)
}

/// Generate an error polynomial with Gaussian distribution in RNS representation.
///
/// This is just an alias for the existing sample_gaussian_simple function
/// to maintain compatibility with the old API naming.
pub fn generate_error_rns<const DEGREE: usize, R: Rng>(
    basis: Arc<crate::RnsBasis>,
    std_dev: f64,
    rng: &mut R,
) -> RnsPolyRing<DEGREE> {
    sample_gaussian_poly(basis, std_dev, rng)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{PublicKeyParams, RnsBasisBuilder, SecretKeyParams};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    fn create_test_basis<const DEGREE: usize>() -> Arc<crate::RnsBasis> {
        Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19, 23])
                .build()
                .unwrap(),
        )
    }

    #[test]
    fn test_generate_ternary_rns() {
        const DEGREE: usize = 8;
        let basis = create_test_basis::<DEGREE>();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let ternary: RnsPolyRing<DEGREE> =
            generate_ternary_rns(basis.clone(), 0.5, &mut rng);

        // Verify it's a valid RNS polynomial
        assert_eq!(ternary.channels(), 3);
        assert_eq!(ternary.len(), DEGREE);

        // Verify coefficients are ternary when reconstructed
        let reconstructed = ternary.to_u64_coefficients();
        let product = basis.primes().iter().product::<u64>();

        for &coeff in &reconstructed {
            // Should be 0, 1, or product-1 (which represents -1)
            assert!(
                coeff == 0 || coeff == 1 || coeff == product - 1,
                "Invalid ternary coefficient: {}",
                coeff
            );
        }
    }

    #[test]
    fn test_encrypt_decrypt_roundtrip() {
        const DEGREE: usize = 4;
        let basis = create_test_basis::<DEGREE>();
        let mut rng = ChaCha20Rng::seed_from_u64(123);

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
        let plaintext_coeffs = [1i64, 2, 3, 4];
        let plaintext: RnsPolyRing<DEGREE> =
            RnsPolyRing::from_integer_coeffs(&plaintext_coeffs, basis.clone());

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

/* //! CKKS encryption and decryption operations.
//!
//! This module contains the core cryptographic operations for the CKKS homomorphic
//! encryption scheme, including functions to encrypt plaintexts and decrypt ciphertexts.
use crate::{
    Ciphertext, PolyRing, PublicKey, SecretKey, generate_error_poly,
    generate_ternary_poly,
};
use rand::Rng;

/// Encrypts a plaintext polynomial using a public key.
///
/// # Encryption Process
/// 1. Generate small random error polynomials e1 and e2
/// 2. Generate a random ephemeral value u (small ternary polynomial)
/// 3. Compute c0 = b*u + e1 + plaintext
/// 4. Compute c1 = a*u + e2
///
/// # Arguments
/// * `plaintext` - The plaintext polynomial to encrypt
/// * `public_key` - The public key used for encryption
/// * `scale` - Scaling factor used in the CKKS scheme
/// * `rng` - Random number generator for encryption randomness
///
/// # Returns
/// A new ciphertext containing the encrypted plaintext
pub fn encrypt<R: Rng>(
    plaintext: &PolyRing,
    public_key: &PublicKey,
    scale: f64,
    rng: &mut R,
) -> Ciphertext {
    let modulus = plaintext.modulus();

    // Generate small random polynomials for encryption
    let ring_dim = plaintext.ring_dim();
    let e1 = generate_error_poly(plaintext.len(), modulus, 3.0, ring_dim, rng);
    let e2 = generate_error_poly(plaintext.len(), modulus, 3.0, ring_dim, rng);

    println!("ring_dim: {}", ring_dim);
    println!("e1: {}", e1);
    println!("e2: {}", e2);

    // Generate random ephemeral value
    let u = generate_ternary_poly(plaintext.len(), modulus, 0.5, rng);
    println!("u: {}", u);

    // c0 = b*u + e1 + plaintext
    let c0 = (public_key.b.clone() * u.clone()) + e1 + plaintext.clone();
    println!("c0: {}", c0);

    // c1 = a*u + e2
    let c1 = (public_key.a.clone() * u) + e2;
    println!("c1: {}", c1);

    Ciphertext {
        c0,
        c1,
        c2: None,
        scale,
    }
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
pub fn decrypt(ciphertext: &Ciphertext, secret_key: &SecretKey) -> PolyRing {
    // Compute c0 + c1*s
    ciphertext.c0.clone() + (ciphertext.c1.clone() * secret_key.s.clone())
}

#[cfg(test)]
mod tests {
    // use super::*;

    // NOTE: This is wrong, Hadamard Product should be used
    // [2.0, 3.5, 5.0, 6.5] * [1.0, 2.0, 3.0, 4.0]
    // 26 x^6 + 39.5 x^5 + 42 x^4 + 35 x^3 + 18 x^2 + 7.5 x + 2
    // 2 + 7.5*x + 18*x^2 + 35*x^3 + 42*(-1) + 39.5*(-x) + 26*(-x^2)
    // (2-42) + (7.5 - 39.5)*x + (18-26)*x^2 + (35-26)*x^3
    // -40 - 32*x - 8*x^2 + 35*x^3

    /* #[test]
    fn test_relinearize_no_c2() {
        let c0 = PolyRing::from_coeffs(&[1, 2], 7, 2);
        let c1 = PolyRing::from_coeffs(&[3, 4], 7, 2);
        let ct = Ciphertext {
            c0: c0.clone(),
            c1: c1.clone(),
            c2: None,
            scale: 1.0,
        };

        // With no c2, relinearize should return the same ciphertext.
        let dummy_relin_key = RelinearizationKey {
            components: vec![],
            base: 2,
            num_components: 1,
        };
        let result = ct.relinearize(&dummy_relin_key);
        assert_eq!(result.c0, ct.c0);
        assert_eq!(result.c1, ct.c1);
        assert!(result.c2.is_none());
    }

    #[test]
    fn test_relinearize_with_c2() {
        let c0 = PolyRing::from_coeffs(&[1, 2], 7, 2);
        let c1 = PolyRing::from_coeffs(&[3, 4], 7, 2);
        let c2 = PolyRing::from_coeffs(&[2, 3], 7, 2);
        let ct = Ciphertext {
            c0: c0.clone(),
            c1: c1.clone(),
            c2: Some(c2.clone()),
            scale: 1.0,
        };

        // Create a relinearization key with one component (b0, a0)
        let b0 = PolyRing::from_coeffs(&[4, 5], 7, 2);
        let a0 = PolyRing::from_coeffs(&[6, 1], 7, 2);
        let relin_key = RelinearizationKey {
            components: vec![(b0.clone(), a0.clone())],
            base: 2,
            num_components: 1,
        };

        // Expected:
        // c2 * b0 = x = [0, 1]
        // c2 * a0 = 2 + 6x = [2, 6]
        // c0_new = [1, 2] + [0, 1] = [1, 3]
        // c1_new = [3, 4] + [2, 6] = [5, (4+6 mod7)= [5, 3]]
        let expected_c0 = PolyRing::from_coeffs(&[1, 3], 7, 2);
        let expected_c1 = PolyRing::from_coeffs(&[5, 3], 7, 2);

        let result = ct.relinearize(&relin_key);
        assert_eq!(result.c0, expected_c0);
        assert_eq!(result.c1, expected_c1);
        assert!(result.c2.is_none());
    } */
} */
