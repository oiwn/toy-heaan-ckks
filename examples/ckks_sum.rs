fn main() {}
/* use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::sync::Arc;
use toy_heaan_ckks::{
    Ciphertext, EncodingParams, PublicKey, PublicKeyParams, RnsBasisBuilder,
    RnsPolyRing, SecretKey, SecretKeyParams, decode, decrypt, encode, encrypt,
};

// Parameters setup
const DEGREE: usize = 8;
const SCALE_BITS: u32 = 30;
const SCALE: f64 = (1u64 << SCALE_BITS) as f64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    // Create RNS basis for polynomial operations
    let basis = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_custom_primes(vec![7829, 6761, 5693]) // Good primes for our degree
            .build()?,
    );

    // Generate Secret Key
    let sk_params: SecretKeyParams<DEGREE> = SecretKeyParams {
        basis: basis.clone(),
        hamming_weight: 4,
    };
    let secret_key = SecretKey::generate(&sk_params, &mut rng)?;

    // Generate Public Key from secret key
    let pk_params = PublicKeyParams {
        basis: basis.clone(),
        error_std: 3.0,
    };
    let public_key = PublicKey::generate(&secret_key, &pk_params, &mut rng)?;

    // Original floating-point values to encrypt
    let values1 = vec![1.5, 2.5, 3.5, 4.5];
    let values2 = vec![0.5, 1.0, 1.5, 2.0];

    println!("Values 1: {:?}", values1);
    println!("Values 2: {:?}", values2);

    // Step 1: ENCODE - Convert floating-point values to polynomial coefficients
    let encoding_params = EncodingParams::new(DEGREE, SCALE_BITS)?;

    let coeffs1 = encode(&values1, &encoding_params)?;
    let coeffs2 = encode(&values2, &encoding_params)?;

    println!("Encoded coeffs 1: {:?}", coeffs1);
    println!("Encoded coeffs 2: {:?}", coeffs2);

    // Convert coefficient vectors to RNS polynomial rings
    let poly1 = RnsPolyRing::from_i64_slice(&coeffs1, basis.clone());
    let poly2 = RnsPolyRing::from_i64_slice(&coeffs2, basis.clone());

    println!("Poly 1: {}", poly1);
    println!("Poly 2: {}", poly2);

    println!("Encoded to polynomials successfully");

    // Step 2: ENCRYPT - Encrypt the polynomial plaintexts
    let ct1 = encrypt(&poly1, &public_key, SCALE, &mut rng);
    let ct2 = encrypt(&poly2, &public_key, SCALE, &mut rng);

    println!("Encrypted both values successfully");

    // Step 3: HOMOMORPHIC OPERATION - Add the two ciphertexts
    let ct_sum = add_ciphertexts(&ct1, &ct2);

    println!("Performed homomorphic addition");

    // Step 4: DECRYPT - Decrypt the sum ciphertext
    let decrypted_poly = decrypt(&ct_sum, &secret_key);

    println!("Decrypted the sum successfully");

    // Step 5: DECODE - Convert decrypted polynomial back to floating-point values
    let decrypted_coeffs = decrypted_poly.to_i64_coefficients();
    let result = decode(&decrypted_coeffs, &encoding_params)?;

    // Calculate expected result for verification
    let expected: Vec<f64> = values1
        .iter()
        .zip(values2.iter())
        .map(|(a, b)| a + b)
        .collect();

    // Print results
    println!("\n=== CKKS Sum Results ===");
    println!("Expected sum: {:?}", expected);
    println!("Computed sum: {:?}", result);

    // Verify accuracy
    let max_error = expected
        .iter()
        .zip(result.iter())
        .map(|(exp, act)| (exp - act).abs())
        .fold(0.0, f64::max);

    println!("Maximum error: {:.2e}", max_error);

    if max_error < 1e-6 {
        println!("✅ Success! Error within acceptable bounds");
    } else {
        println!("❌ Warning: Error is higher than expected");
    }

    Ok(())
}

/// Homomorphic addition of two ciphertexts
/// In CKKS, addition is component-wise:
/// (c0_1, c1_1) + (c0_2, c1_2) = (c0_1+c0_2, c1_1+c1_2)
fn add_ciphertexts<const DEGREE: usize>(
    ct1: &Ciphertext<DEGREE>,
    ct2: &Ciphertext<DEGREE>,
) -> Ciphertext<DEGREE> {
    // Add corresponding components
    let c0_sum = &ct1.c0 + &ct2.c0;
    let c1_sum = &ct1.c1 + &ct2.c1;

    // Scale should be the same for both ciphertexts in addition
    let scale = ct1.scale;

    Ciphertext::new(c0_sum, c1_sum, scale)
} */
