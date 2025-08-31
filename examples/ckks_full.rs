use crypto_bigint::{NonZero, U256};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    BigIntPolyRing, CkksEngine, Encoder, crypto::operations::multiply_ciphertexts,
    encoding,
};

const DEGREE: usize = 16;
const SCALE_BITS: u32 = 30; // Reduce scale for debugging
const Q50: u64 = (1u64 << 50) - 27;

type Engine = CkksEngine<BigIntPolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS BigInt U256 Backend Multiplication Demo");

    let modulus_u256 = U256::from_u128(Q50 as u128);
    let modulus = NonZero::new(modulus_u256).expect("Modulus should not be zero");
    let encoder = encoding::BigIntEncoder::new(SCALE_BITS)?;

    println!("✅ Using modulus: {}", modulus_u256);

    // Create CKKS Engine with BigInt U256 backend
    let engine = Engine::builder()
        .error_variance(0.2)
        .hamming_weight(DEGREE / 2)
        .build_bigint_u256(modulus, SCALE_BITS)?;
    println!("✅ Engine configured with BigInt U256 backend");

    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    let relin_key = engine.generate_relinearization_key(&secret_key, &mut rng)?;
    println!("✅ Secret and public keys generated with U256 arithmetic");

    // Input data
    let values1 = vec![1.0]; // Use only 1 slot for power-of-2 requirement
    let values2 = vec![2.0]; // Should give [2.0] after multiplication
    println!("📊 Input data1: {values1:?}");
    println!("📊 Input data2: {values2:?}");

    // Step 1: Encode: Vec<f64> -> Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext1 = encoder.encode(&values1, engine.context());
    println!(
        "Plaintext1 coefficient[0]: {:x}",
        plaintext1.poly.coeffs[0].to_words()[0]
    );
    let plaintext2 = encoder.encode(&values2, engine.context());
    println!(
        "Plaintext2 coefficient[0]: {:x}",
        plaintext2.poly.coeffs[0].to_words()[0]
    );
    println!("✅ Values encoded to plaintext with BigInt backend");

    // sanity: encode -> decode without encryption
    let decoded_plain = encoder.decode(&plaintext1);
    println!("Plain roundtrip: {:?}", &decoded_plain);
    println!(
        "  First few coeffs for plaintext1: {:?}",
        &plaintext1.poly.coeffs[0..4]
    );

    // Step 2: Encrypt: Plaintext → Ciphertext
    println!("\n🔒 Encrypting plaintext...");
    let ciphertext1 = engine.encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext using U256 operations");

    // Let's test with different target scales to see what's happening
    println!("\\nTesting multiplication with different rescaling approaches...");

    // Test 1: Use doubled scale bits (no rescaling)
    println!("Test 1: No rescaling (doubled scale bits)");
    let ciphertext_mul_no_rescale = multiply_ciphertexts(
        &ciphertext1,
        &ciphertext2,
        &relin_key,
        2 * SCALE_BITS,
    )
    .expect("Multiplication should succeed");
    let decrypted_no_rescale =
        Engine::decrypt(&ciphertext_mul_no_rescale, &secret_key);
    let decoded_no_rescale = encoder.decode(&decrypted_no_rescale);
    println!("  No rescaling result: {:?}", decoded_no_rescale);

    // Test 2: Use original scale bits (standard rescaling)
    println!("Test 2: Standard rescaling to original scale bits");
    let ciphertext_mul =
        multiply_ciphertexts(&ciphertext1, &ciphertext2, &relin_key, SCALE_BITS)
            .expect("Multiplication should succeed");

    println!("✅ Multiplication with relinearization and rescaling complete");
    println!("  Result scale_bits: {}", ciphertext_mul.scale_bits);
    println!("  c0[0]: {:x}", ciphertext_mul.c0.coeffs[0].to_words()[0]);
    println!("  c1[0]: {:x}", ciphertext_mul.c1.coeffs[0].to_words()[0]);

    // Step 3: Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext_mul, &secret_key);
    println!("✅ Ciphertext decrypted back to plaintext");

    println!("Debug scales:");
    println!(
        "  decrypted_plaintext.scale_bits: {}",
        decrypted_plaintext.scale_bits
    );

    let coeffs = decrypted_plaintext.poly.coefficients();
    println!("  First few coeffs: {:?}", &coeffs[0..4]);

    // Step 4: Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding back to floating-point...");
    let decoded_values = encoder.decode(&decrypted_plaintext);
    println!("✅ Plaintext1 decoded: {:?}", decoded_values);

    Ok(())
}
