use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{CkksEngine, Encoder, NaivePolyRing, encoding::NaiveEncoder};

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 30;
const MODULUS: u64 = (1u64 << 50) - 27;

type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS Naive Encoder Demo");

    let encoder = NaiveEncoder::new(SCALE_BITS)?;
    println!("✅ Using modulus: {} (50-bit)", MODULUS);

    // Create CKKS Engine with Naive backend
    let engine = Engine::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_naive(MODULUS, SCALE_BITS)?;
    println!("✅ Engine configured with Naive backend");

    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    println!("✅ Secret and public keys generated");

    // Input data (same as your working example)
    let values = vec![1.0, 2.0, 3.0, 4.0];
    println!("\n📊 Input data: {:?}", values);

    // Step 1: Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext = encoder.encode(&values, engine.context());
    println!("✅ Values encoded to plaintext with Naive encoder");

    // Step 2: Encrypt: Plaintext → Ciphertext
    println!("\n🔒 Encrypting plaintext...");
    let ciphertext = engine.encrypt(&plaintext, &public_key, SCALE_BITS, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext");

    // Step 3: Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext, &secret_key);
    println!("✅ Ciphertext decrypted back to plaintext");

    // Step 4: Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding back to floating-point...");
    let decoded_values = encoder.decode(&decrypted_plaintext);
    println!(
        "✅ Plaintext decoded: {:?}",
        &decoded_values[..values.len()]
    );

    // Display results
    println!("\n📊 Results:");
    println!("  Original: {:?}", values);
    println!("  Decoded:  {:?}", &decoded_values[..values.len()]);

    // Verify accuracy (similar to your working example)
    let max_error = values
        .iter()
        .zip(&decoded_values)
        .map(|(orig, decoded)| (orig - decoded).abs())
        .fold(0.0, f64::max);

    println!("  Max error: {:.2e}", max_error);

    if max_error < 1e-9 {
        println!(
            "🎉 Success! Full CKKS encrypt/decrypt pipeline works with Naive encoder!"
        );
    } else if max_error < 1e-6 {
        println!("✅ Good! CKKS pipeline works with acceptable precision");
    } else {
        println!("⚠️  Warning: Error is higher than expected");
    }

    // Demonstrate homomorphic addition
    println!("\n➕ Testing homomorphic addition...");

    let values2 = vec![-1.0, -2.0, -3.0, -4.0];
    println!("Second input: {:?}", values2);

    let plaintext2 = encoder.encode(&values2, engine.context());
    let ciphertext2 =
        engine.encrypt(&plaintext2, &public_key, SCALE_BITS, &mut rng);

    // Homomorphic addition
    let ciphertext_sum = Engine::add_ciphertexts(&ciphertext, &ciphertext2);
    let decrypted_sum = Engine::decrypt(&ciphertext_sum, &secret_key);
    let decoded_sum = encoder.decode(&decrypted_sum);

    // Expected result: [1,2,3,4] + [-1,-2,-3,-4] = [0,0,0,0]
    let expected_sum: Vec<f64> =
        values.iter().zip(&values2).map(|(a, b)| a + b).collect();

    println!("\n📊 Homomorphic Addition Results:");
    println!("  Input 1:    {:?}", values);
    println!("  Input 2:    {:?}", values2);
    println!("  Expected:   {:?}", expected_sum);
    println!("  Computed:   {:?}", &decoded_sum[..expected_sum.len()]);

    // Verify homomorphic addition accuracy
    let add_max_error = expected_sum
        .iter()
        .zip(&decoded_sum)
        .map(|(exp, comp)| (exp - comp).abs())
        .fold(0.0, f64::max);

    println!("  Add error:  {:.2e}", add_max_error);

    if add_max_error < 1e-9 {
        println!("🎉 Homomorphic addition successful with Naive encoder!");
    } else if add_max_error < 1e-6 {
        println!("✅ Homomorphic addition works with acceptable precision");
    } else {
        println!("⚠️  Warning: Addition error higher than expected");
    }

    Ok(())
}
