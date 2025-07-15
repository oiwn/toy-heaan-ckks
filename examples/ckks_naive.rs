use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    CkksEngine, NaivePolyRing, encoding::EncoderType, rings::BackendType,
};

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 40;
const MODULUS: u64 = 741507920154517877u64;
type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS Abstract API Demo with Encryption");
    // Create CKKS Engine with context
    let engine = Engine::builder()
        .encoder(EncoderType::RustFft)
        .backend(BackendType::Naive(MODULUS))
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .scale_bits(SCALE_BITS)
        .build()?;
    println!("✅ Engine configured with builder pattern");

    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    println!(
        "✅ Secret and public keys generated: \n{secret_key:?}\n{public_key:?}"
    );

    // Input data (small vector to test)
    let values = vec![1.5, 2.5, 3.5];
    println!("\n📊 Input data: {:?}", values);

    // Step 1: Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext = engine.encode(&values);
    println!("✅ Values encoded to plaintext: {:?}", plaintext);

    // Step 2: Encrypt: Plaintext → Ciphertext
    println!("\n🔒 Encrypting plaintext...");
    let ciphertext = engine.encrypt(&plaintext, &public_key, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext: {:?}", ciphertext);

    // Step 3: Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext, &secret_key);
    println!(
        "✅ Ciphertext decrypted back to plaintext: {:?}",
        decrypted_plaintext
    );

    // Step 4: Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding back to floating-point...");
    let decoded_values = Engine::decode(
        EncoderType::RustFft,
        &decrypted_plaintext,
        engine.params.scale_bits,
    );
    println!("✅ Plaintext decoded: {:?}", decoded_values);

    // Display results
    println!("\n📊 Results:");
    println!("  Original: {:?}", values);
    println!("  Decoded:  {:?}", &decoded_values[..values.len()]);

    // Verify accuracy
    let max_error = values
        .iter()
        .zip(&decoded_values)
        .map(|(orig, decoded)| (orig - decoded).abs())
        .fold(0.0, f64::max);

    println!("  Max error: {:.2e}", max_error);

    if max_error < 1e-3 {
        println!("🎉 Success! Full CKKS encrypt/decrypt pipeline works!");
    } else {
        println!("⚠️  Warning: Error is higher than expected");
    }

    // Demonstrate homomorphic addition
    println!("\n➕ Testing homomorphic addition...");

    let values2 = vec![0.5, 1.0, 1.5];
    println!("Second input: {:?}", values2);

    let plaintext2 = engine.encode(&values2);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);

    // Homomorphic addition
    let ciphertext_sum = Engine::add_ciphertexts(&ciphertext, &ciphertext2);
    let decrypted_sum = Engine::decrypt(&ciphertext_sum, &secret_key);
    let decoded_sum = Engine::decode(
        EncoderType::RustFft,
        &decrypted_sum,
        engine.params.scale_bits,
    );

    let expected_sum: Vec<f64> =
        values.iter().zip(&values2).map(|(a, b)| a + b).collect();

    println!("Expected sum: {:?}", expected_sum);
    println!("Computed sum: {:?}", &decoded_sum[..values.len()]);

    let sum_error = expected_sum
        .iter()
        .zip(&decoded_sum)
        .map(|(exp, act)| (exp - act).abs())
        .fold(0.0, f64::max);

    println!("Sum error: {:.2e}", sum_error);

    if sum_error < 1e-3 {
        println!("🎉 Homomorphic addition works perfectly!");
    }

    println!("\n💡 Full CKKS pipeline completed:");
    println!("  1. ✅ Key generation (secret + public keys)");
    println!("  2. ✅ Encoding (float → plaintext)");
    println!("  3. ✅ Encryption (plaintext → ciphertext)");
    println!("  4. ✅ Homomorphic operations (encrypted addition)");
    println!("  5. ✅ Decryption (ciphertext → plaintext)");
    println!("  6. ✅ Decoding (plaintext → float)");

    Ok(())
}
