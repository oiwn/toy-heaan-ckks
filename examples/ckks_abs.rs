use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    CkksEngine, EncodingParams, NaivePolyRing, PublicKeyParams, SecretKeyParams,
};

const DEGREE: usize = 8;
type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Setup context for NaivePolyRing
    let modulus = (1u64 << 60) - 1; // Large prime-like modulus
    let context = modulus; // For NaivePolyRing, context is just the modulus

    println!("🔐 CKKS Abstract API Demo");
    // Create CKKS Engine with context
    let engine: CkksEngine<NaivePolyRing<DEGREE>, DEGREE> =
        CkksEngine::new(context);

    // Setup encoding parameters
    let scale_bits = 30;
    let enc_params = EncodingParams::<DEGREE>::new(scale_bits)?;

    // Initialize RNG with fixed seed for reproducible results
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

    // Setup encoding parameters
    let scale_bits = 30;
    let enc_params = EncodingParams::<DEGREE>::new(scale_bits)?;

    println!("✅ Engine and parameters configured");

    // Input data (small vector to test)
    let values = vec![1.5, 2.5, 3.5];
    println!("\n📊 Input data: {:?}", values);

    // Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext = engine.encode(&values, &enc_params);
    println!("✅ Values encoded to plaintext");

    // Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding back to floating-point...");
    let decoded_values = engine.decode(&plaintext, &enc_params);
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

    if max_error < 1e-6 {
        println!("🎉 Success! Encode/decode round trip works!");
    } else {
        println!("⚠️  Warning: Error is higher than expected");
    }

    // 🔑 Key generation
    // println!("\n🔑 Generating keys...");
    // let secret_key = Engine::generate_secret_key(&sk_params, &mut rng);
    // let public_key = Engine::generate_public_key(&secret_key, &pk_params, &mut rng);
    // println!("✅ Secret and public keys generated");

    /* // Input data
    let values1 = vec![1.5, 2.5, 3.5, 4.5];
    let values2 = vec![0.5, 1.0, 1.5, 2.0];
    let expected_sum: Vec<f64> =
        values1.iter().zip(&values2).map(|(a, b)| a + b).collect();

    println!("\n📊 Input data:");
    println!("  Values 1: {:?}", values1);
    println!("  Values 2: {:?}", values2);
    println!("  Expected sum: {:?}", expected_sum);

    // Encode floating-point values to plaintexts
    println!("\n🔢 Encoding values...");
    let plaintext1 = Engine::encode(&values1, &enc_params);
    let plaintext2 = Engine::encode(&values2, &enc_params);
    println!("✅ Values encoded to plaintexts");
    println!("Plaintext 1: {plaintext1:?}");
    println!("Plaintext 1: {plaintext2:?}");

    // Encrypt plaintexts to ciphertexts
    println!("\n🔒 Encrypting plaintexts...");
    let ciphertext1 = Engine::encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = Engine::encrypt(&plaintext2, &public_key, &mut rng);
    println!("✅ Plaintexts encrypted to ciphertexts");

    // Homomorphic addition
    // println!("\n➕ Performing homomorphic addition...");
    // let ciphertext_sum = Engine::add_ciphertexts(&ciphertext1, &ciphertext2);
    // println!("✅ Ciphertexts added homomorphically");

    // Decrypt the result
    println!("\n🔓 Decrypting result...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext1, &secret_key);
    println!("✅ Result decrypted");

    // Decode back to floating-point values
    println!("\n🔢 Decoding to floating-point...");
    let result_values = Engine::decode(&decrypted_plaintext, &enc_params);
    println!("✅ Result decoded");

    // Display results
    println!("\n📊 Results:");
    println!("  Expected: {:?}", expected_sum);
    println!("  Computed: {:?}", result_values);

    // Verify accuracy
    let max_error = expected_sum
        .iter()
        .zip(&result_values)
        .map(|(exp, act)| (exp - act).abs())
        .fold(0.0, f64::max);

    println!("  Max error: {:.2e}", max_error);

    if max_error < 1e-3 {
        println!("🎉 Success! CKKS homomorphic addition works correctly!");
    } else {
        println!("⚠️  Warning: Error is higher than expected");
    }

    // Show what makes this API great
    println!("\n💡 What makes this API powerful:");
    println!("  • Compile-time polynomial degree (const DEGREE)");
    println!(
        "  • Pluggable polynomial backends (RnsPoly, CoefficientPoly, NttPoly)"
    );
    println!("  • Type-safe operations (keys match ciphertext types)");
    println!("  • Zero-cost abstractions (no runtime overhead)");
    */

    Ok(())
}

fn small_params() -> (
    SecretKeyParams<DEGREE>,
    PublicKeyParams<DEGREE>,
    EncodingParams<DEGREE>,
) {
    let hamming_weight = 3;
    let error_std = 3.0;
    let scale_bits = 30;

    let sk_params: SecretKeyParams<DEGREE> = SecretKeyParams { hamming_weight };

    let pk_params: PublicKeyParams<DEGREE> = PublicKeyParams { error_std };

    let enc_params = EncodingParams::new(scale_bits).unwrap();

    (sk_params, pk_params, enc_params)
}
