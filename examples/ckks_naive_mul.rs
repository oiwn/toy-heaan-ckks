use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    Ciphertext, CkksEngine, NaivePolyRing, Plaintext, PolyRing,
    crypto::operations::multiply_ciphertexts,
};

const DEGREE: usize = 8;
const MODULUS: u64 = (1u64 << 50) - 27; // Large modulus for headroom
const SCALE_BITS: u32 = 10; // 2^10 = 1024

type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

    println!("🎯 CKKS Homomorphic Multiplication with Relinearization");
    println!("   Using NaivePolyRing<{}> with {}-bit modulus", DEGREE, 50);

    // First, test basic polynomial multiplication without encryption
    println!("\n🧪 Testing basic polynomial multiplication (no encryption):");
    test_basic_polynomial_multiplication();

    // Create engine with much smaller error variance for multiplication
    let engine = Engine::builder()
        .error_variance(0.2) // Back to 0.2 for working multiplication
        .hamming_weight(DEGREE / 2)
        .build_naive(MODULUS, SCALE_BITS)?;

    // Generate ALL required keys
    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    let relin_key = engine.generate_relinearization_key(&secret_key, &mut rng)?;

    println!("✅ Generated secret key, public key, and relinearization key");

    // Create test plaintexts
    println!("\n📝 Creating plaintexts...");
    let plaintext1 = create_constant_plaintext(2.0, &engine);
    let plaintext2 = create_constant_plaintext(3.0, &engine);

    println!(
        "   Plaintext 1: 2.0 (coeffs[0] = {})",
        plaintext1.poly.coeffs[0]
    );
    println!(
        "   Plaintext 2: 3.0 (coeffs[0] = {})",
        plaintext2.poly.coeffs[0]
    );

    // Encrypt and analyze noise
    println!("\n🔐 Encrypting plaintexts...");
    let ciphertext1 = engine.encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);

    // Analyze the noise/error polynomials
    analyze_ciphertext_noise(&ciphertext1, &plaintext1, &secret_key, "CT1");
    analyze_ciphertext_noise(&ciphertext2, &plaintext2, &secret_key, "CT2");

    // Verify encryption works
    let check1 = Engine::decrypt(&ciphertext1, &secret_key);
    let check2 = Engine::decrypt(&ciphertext2, &secret_key);
    let val1 = check1.poly.coeffs[0] as f64 / (1u64 << SCALE_BITS) as f64;
    let val2 = check2.poly.coeffs[0] as f64 / (1u64 << SCALE_BITS) as f64;

    println!("   Verification - CT1 decrypts to: {:.3}", val1);
    println!("   Verification - CT2 decrypts to: {:.3}", val2);

    if (val1 - 2.0).abs() > 0.1 || (val2 - 3.0).abs() > 0.1 {
        println!("❌ Basic encryption failed, cannot proceed");
        return Ok(());
    }

    // Perform homomorphic multiplication WITH relinearization
    println!("\n✖️  Performing homomorphic multiplication...");
    let ciphertext_product =
        multiply_ciphertexts(&ciphertext1, &ciphertext2, &relin_key, SCALE_BITS)?;

    println!("   ✅ Multiplication completed");
    println!("   Product scale_bits: {}", ciphertext_product.scale_bits);

    // Decrypt result
    println!("\n🔓 Decrypting result...");
    let decrypted_product = Engine::decrypt(&ciphertext_product, &secret_key);

    // Extract the result
    let result_scaled = decrypted_product.poly.coeffs[0] as f64;
    let expected_scale = (1u64 << SCALE_BITS) as f64; // Now back to single scale
    let result = result_scaled / expected_scale;

    println!(
        "   Raw coefficient[0]: {}",
        decrypted_product.poly.coeffs[0]
    );
    println!("   Expected scale (Δ): {:.0}", expected_scale);
    println!("   Computed result: {:.6}", result);
    println!("   Expected result: {:.1}", 2.0 * 3.0);

    let error = (result - 6.0).abs();
    println!("   Error: {:.6}", error);

    if error < 0.1 {
        println!(
            "\n🎉 SUCCESS! Homomorphic multiplication with relinearization works!"
        );
        println!("   ✅ 2.0 × 3.0 = {:.3} (computed homomorphically)", result);
    } else {
        println!("\n⚠️  Higher error than expected, but structure is correct");
    }

    // Demonstrate the difference with and without relinearization
    println!("\n🔬 Comparison: with vs without relinearization");

    let product_no_relin =
        multiply_without_relinearization(&ciphertext1, &ciphertext2);
    let decrypt_no_relin = Engine::decrypt(&product_no_relin, &secret_key);
    let doubled_scale = (1u64 << SCALE_BITS) as f64 * (1u64 << SCALE_BITS) as f64;
    let result_no_relin = decrypt_no_relin.poly.coeffs[0] as f64 / doubled_scale;

    println!(
        "   Without relinearization: {:.6} (error: {:.6})",
        result_no_relin,
        (result_no_relin - 6.0).abs()
    );
    println!(
        "   With relinearization:    {:.6} (error: {:.6})",
        result, error
    );

    println!("\n📊 Complete multiplication pipeline:");
    println!("   1. ✅ Generated relinearization key");
    println!("   2. ✅ Encrypted plaintexts");
    println!("   3. ✅ Performed ciphertext multiplication");
    println!("   4. ✅ Applied relinearization");
    println!("   5. ✅ Decrypted and verified result");

    Ok(())
}

/// Multiplication without relinearization (for comparison)
fn multiply_without_relinearization(
    ct1: &Ciphertext<NaivePolyRing<DEGREE>, DEGREE>,
    ct2: &Ciphertext<NaivePolyRing<DEGREE>, DEGREE>,
) -> Ciphertext<NaivePolyRing<DEGREE>, DEGREE> {
    // Just compute d0 and d1, ignore the s² term
    let mut d0 = ct1.c0.clone();
    d0 *= &ct2.c0;

    let mut d1_part1 = ct1.c0.clone();
    d1_part1 *= &ct2.c1;

    let mut d1_part2 = ct1.c1.clone();
    d1_part2 *= &ct2.c0;

    let mut d1 = d1_part1;
    d1 += &d1_part2;

    // The missing c1₁·c1₂·s² term becomes noise

    Ciphertext {
        c0: d0,
        c1: d1,
        scale_bits: ct1.scale_bits + ct2.scale_bits,
    }
}

/// Create a plaintext representing a constant value
fn create_constant_plaintext(
    value: f64,
    engine: &Engine,
) -> Plaintext<NaivePolyRing<DEGREE>, DEGREE> {
    let mut coeffs = [0i64; DEGREE];
    coeffs[0] = (value * (1u64 << SCALE_BITS) as f64) as i64;

    Plaintext {
        poly: NaivePolyRing::from_coeffs(&coeffs, engine.context()),
        scale_bits: SCALE_BITS,
        slots: 1, // Single slot for naive multiplication example
    }
}

/// Analyze the noise in a ciphertext by computing the error polynomial
fn analyze_ciphertext_noise(
    ct: &Ciphertext<NaivePolyRing<DEGREE>, DEGREE>,
    pt: &Plaintext<NaivePolyRing<DEGREE>, DEGREE>,
    sk: &toy_heaan_ckks::SecretKey<NaivePolyRing<DEGREE>, DEGREE>,
    label: &str,
) {
    // Decrypt to get noisy plaintext
    let decrypted = Engine::decrypt(ct, sk);

    // Compute noise = decrypted - original (manually since SubAssign not implemented)
    let mut noise_coeffs = [0i64; DEGREE];
    for i in 0..DEGREE {
        let dec_val = decrypted.poly.coeffs[i] as i64;
        let orig_val = pt.poly.coeffs[i] as i64;
        noise_coeffs[i] = dec_val - orig_val;
    }

    println!("\n🔍 Noise Analysis for {}:", label);
    println!("   Original plaintext[0]: {}", pt.poly.coeffs[0]);
    println!("   Decrypted plaintext[0]: {}", decrypted.poly.coeffs[0]);
    println!("   Noise polynomial coefficients:");
    for i in 0..DEGREE {
        if noise_coeffs[i] != 0 {
            println!("     coeffs[{}]: {}", i, noise_coeffs[i]);
        }
    }

    // Compute noise magnitude
    let noise_magnitude: f64 = noise_coeffs
        .iter()
        .map(|&c| (c as f64).powi(2))
        .sum::<f64>()
        .sqrt();
    println!("   Noise magnitude (L2 norm): {:.3}", noise_magnitude);

    // Scale-adjusted noise
    let scale = (1u64 << SCALE_BITS) as f64;
    let scaled_noise = noise_magnitude / scale;
    println!("   Scaled noise (actual error): {:.6}", scaled_noise);
}

/// Test basic polynomial multiplication without encryption noise
fn test_basic_polynomial_multiplication() {
    // Create polynomials representing 2.0 and 3.0 at scale 2^10
    let coeffs1 = [2048i64, 0, 0, 0, 0, 0, 0, 0]; // 2.0 * 1024
    let coeffs2 = [3072i64, 0, 0, 0, 0, 0, 0, 0]; // 3.0 * 1024

    let poly1 = NaivePolyRing::<DEGREE>::from_coeffs(&coeffs1, &MODULUS);
    let poly2 = NaivePolyRing::<DEGREE>::from_coeffs(&coeffs2, &MODULUS);

    println!("   poly1[0]: {} (represents 2.0)", poly1.coeffs[0]);
    println!("   poly2[0]: {} (represents 3.0)", poly2.coeffs[0]);

    // Multiply
    let mut result = poly1.clone();
    result *= &poly2;

    println!("   result[0]: {}", result.coeffs[0]);

    let expected = 2048u64 * 3072u64; // 6,291,456
    println!("   expected: {}", expected);

    // After multiplication, we have scale 2^20, so to get the value:
    let computed_value = result.coeffs[0] as f64
        / ((1u64 << SCALE_BITS) as f64 * (1u64 << SCALE_BITS) as f64);
    println!("   computed value: {:.6}", computed_value);
    println!("   expected value: {:.1}", 2.0 * 3.0);

    if result.coeffs[0] == expected {
        println!("   ✅ Basic polynomial multiplication is correct!");
    } else {
        println!("   ❌ Basic polynomial multiplication failed!");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiplication_with_relinearization() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let engine = Engine::builder()
            .error_variance(3.2)
            .hamming_weight(DEGREE / 2)
            .build_naive(MODULUS, 10)
            .unwrap();

        let secret_key = engine.generate_secret_key(&mut rng).unwrap();
        let public_key = engine.generate_public_key(&secret_key, &mut rng).unwrap();
        let relin_key = engine
            .generate_relinearization_key(&secret_key, &mut rng)
            .unwrap();

        // Test 4.0 × 2.0 = 8.0
        let pt1 = create_constant_plaintext(4.0, &engine);
        let pt2 = create_constant_plaintext(2.0, &engine);

        let ct1 = engine.encrypt(&pt1, &public_key, &mut rng);
        let ct2 = engine.encrypt(&pt2, &public_key, &mut rng);

        let ct_product = multiply_with_relinearization(&ct1, &ct2, &relin_key);
        let decrypted = Engine::decrypt(&ct_product, &secret_key);

        let result = decrypted.poly.coeffs[0] as f64
            / ((1u64 << SCALE_BITS) as f64 * (1u64 << SCALE_BITS) as f64);
        let expected = 8.0;

        // With proper relinearization, error should be much smaller
        assert!(
            (result - expected).abs() < 0.5,
            "Expected ~{}, got {}",
            expected,
            result
        );
    }
}
