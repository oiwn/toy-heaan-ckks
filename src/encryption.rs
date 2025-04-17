use crate::{
    Ciphertext, PolyRing, PublicKey, SecretKey, generate_error_poly,
    generate_ternary_poly,
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

/// Encrypt a plaintext polynomial using the public key
pub fn encrypt(
    plaintext: &PolyRing,
    public_key: &PublicKey,
    scale: f64,
) -> Ciphertext {
    let modulus = plaintext.modulus();
    let mut rng = ChaCha20Rng::from_seed([2u8; 32]); // Different seed from keys

    // Generate small random polynomials for encryption
    let e1 = generate_error_poly(
        plaintext.len(),
        modulus,
        3.0,
        plaintext.len(),
        &mut rng,
    );
    let e2 = generate_error_poly(
        plaintext.len(),
        modulus,
        3.0,
        plaintext.len(),
        &mut rng,
    );

    // Generate random ephemeral value
    let u = generate_ternary_poly(plaintext.len(), modulus, 0.5, &mut rng);

    // c0 = b*u + e1 + plaintext
    let c0 = (public_key.b.clone() * u.clone()) + e1 + plaintext.clone();

    // c1 = a*u + e2
    let c1 = (public_key.a.clone() * u) + e2;

    Ciphertext {
        c0,
        c1,
        c2: None,
        scale,
    }
}

/// Decrypt a ciphertext using the secret key
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
}
