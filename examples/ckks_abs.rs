use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    CkksEngine, EncodingParams, NaivePolyRing, PublicKeyParams, SecretKeyParams,
    encoding::EncoderType, rings::BackendType,
};

const DEGREE: usize = 8;
type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔐 CKKS Abstract API Demo with Encryption");
    // Create CKKS Engine with context
    let engine = Engine::builder()
        .encoder(EncoderType::RustFft)
        .backend(BackendType::Naive)
        .error_variance(3.2)
        .hamming_weight(4)
        .scale_bits(30)
        .build()?;
    println!("✅ Engine configured with builder pattern");

    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    let enc_params = EncodingParams::<DEGREE>::new(30)?; // scale_bits = 30

    println!("\n🔑 Generating keys...");
    let sk_params = SecretKeyParams::new(4)?; // hamming weight = 4  
    let pk_params = PublicKeyParams::new(3.2)?; // error variance = 3.2
    let secret_key = Engine::generate_secret_key(&sk_params, &mut rng)?;
    let public_key =
        engine.generate_public_key(&secret_key, &pk_params, &mut rng)?;
    println!("✅ Secret and public keys generated");

    // Input data (small vector to test)
    let values = vec![1.5, 2.5, 3.5];
    println!("\n📊 Input data: {:?}", values);

    // Step 1: Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext = engine.encode(&values, &enc_params);
    println!("✅ Values encoded to plaintext");

    // Step 2: Encrypt: Plaintext → Ciphertext
    println!("\n🔒 Encrypting plaintext...");
    let ciphertext = engine.encrypt(&plaintext, &public_key, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext");

    // Step 3: Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext, &secret_key);
    println!("✅ Ciphertext decrypted back to plaintext");

    // Step 4: Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding back to floating-point...");
    let decoded_values = engine.decode(&decrypted_plaintext, &enc_params);
    println!("✅ Plaintext decoded");

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

    let plaintext2 = engine.encode(&values2, &enc_params);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);

    // Homomorphic addition
    let ciphertext_sum = Engine::add_ciphertexts(&ciphertext, &ciphertext2);
    let decrypted_sum = Engine::decrypt(&ciphertext_sum, &secret_key);
    let decoded_sum = engine.decode(&decrypted_sum, &enc_params);

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
