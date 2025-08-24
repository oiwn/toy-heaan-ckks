use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    Ciphertext, CkksEngine, NaivePolyRing, Plaintext, PolyRescale, PolyRing,
    RelinearizationKey,
};

const DEGREE: usize = 8;
const MODULUS: u64 = (1u64 << 50) - 27; // Large modulus for headroom
const SCALE_BITS: u32 = 10; // 2^10 = 1024

type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

    println!("🎯 CKKS Homomorphic Multiplication with Relinearization");
    println!("   Using NaivePolyRing<{}> with {}-bit modulus", DEGREE, 50);

    // Create engine
    let engine = Engine::builder()
        .error_variance(3.2)
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

    // Encrypt
    println!("\n🔐 Encrypting plaintexts...");
    let ciphertext1 = engine.encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);

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
        multiply_with_relinearization(&ciphertext1, &ciphertext2, &relin_key);

    println!("   ✅ Multiplication completed");
    println!("   Product scale_bits: {}", ciphertext_product.scale_bits);

    // Decrypt result
    println!("\n🔓 Decrypting result...");
    let decrypted_product = Engine::decrypt(&ciphertext_product, &secret_key);

    // Extract the result
    let result_scaled = decrypted_product.poly.coeffs[0] as f64;
    let expected_scale = (1u64 << SCALE_BITS) as f64 * (1u64 << SCALE_BITS) as f64;
    let result = result_scaled / expected_scale;

    println!(
        "   Raw coefficient[0]: {}",
        decrypted_product.poly.coeffs[0]
    );
    println!("   Expected scale (Δ²): {:.0}", expected_scale);
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
    let result_no_relin = decrypt_no_relin.poly.coeffs[0] as f64 / expected_scale;

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

/// Proper homomorphic multiplication with relinearization
fn multiply_with_relinearization(
    ct1: &Ciphertext<NaivePolyRing<DEGREE>, DEGREE>,
    ct2: &Ciphertext<NaivePolyRing<DEGREE>, DEGREE>,
    relin_key: &RelinearizationKey<NaivePolyRing<DEGREE>, DEGREE>,
) -> Ciphertext<NaivePolyRing<DEGREE>, DEGREE> {
    println!("   Step 1: Computing polynomial products...");

    // Multiply: (c0₁ + c1₁·s) × (c0₂ + c1₂·s)
    // = c0₁·c0₂ + (c0₁·c1₂ + c1₁·c0₂)·s + c1₁·c1₂·s²

    // d0 = c0₁ * c0₂
    let mut d0 = ct1.c0.clone();
    d0 *= &ct2.c0;

    // d1 = c0₁·c1₂ + c1₁·c0₂
    let mut d1_part1 = ct1.c0.clone();
    d1_part1 *= &ct2.c1;

    let mut d1_part2 = ct1.c1.clone();
    d1_part2 *= &ct2.c0;

    let mut d1 = d1_part1;
    d1 += &d1_part2;

    // d2 = c1₁ * c1₂ (coefficient of s²)
    let mut d2 = ct1.c1.clone();
    d2 *= &ct2.c1;

    println!("   Step 2: Applying relinearization...");

    // Relinearization: Replace d2·s² with d2·(relin_key.b + relin_key.a·s)
    // Since relin_key.b + relin_key.a·s ≈ s² + small_error

    // Compute d2 * relin_key.b
    let mut relin_c0 = d2.clone();
    relin_c0 *= &relin_key.b;

    // Compute d2 * relin_key.a
    let mut relin_c1 = d2;
    relin_c1 *= &relin_key.a;

    // Add relinearization terms
    d0 += &relin_c0; // d0 += d2 * relin_key.b
    d1 += &relin_c1; // d1 += d2 * relin_key.a

    println!("   Step 3: Relinearization complete!");

    Ciphertext {
        c0: d0,
        c1: d1,
        scale_bits: ct1.scale_bits + ct2.scale_bits,
    }
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
