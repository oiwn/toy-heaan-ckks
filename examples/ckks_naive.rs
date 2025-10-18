use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{CkksEngine, Encoder, NaivePolyRing, RustFftEncoder};

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 20;
const MODULUS: u64 = 2u64.pow(61) - 1;
type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS Abstract API Demo with Encryption");
    // Create CKKS Engine with context
    let engine = Engine::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_naive(MODULUS, SCALE_BITS)?;
    let encoder = RustFftEncoder::new(engine.params.scale_bits)?;
    println!("✅ Engine configured with builder pattern");

    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    println!(
        "✅ Secret and public keys generated: \n\t{secret_key:?}\n\t{public_key:?}"
    );

    // Input data (small vector to test)
    let values = vec![1.5, 2.5, 3.5];
    println!("\n📊 Input data: {:?}", values);

    // Step 1: Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext = encoder.encode(&values, engine.context());
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
    let decoded_values = encoder.decode(&decrypted_plaintext);
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

    let plaintext2 = encoder.encode(&values2, engine.context());
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);

    // Homomorphic addition
    let ciphertext_sum = Engine::add_ciphertexts(&ciphertext, &ciphertext2);
    let decrypted_sum = Engine::decrypt(&ciphertext_sum, &secret_key);
    let decoded_sum = encoder.decode(&decrypted_sum);

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
