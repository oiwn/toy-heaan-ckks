use crypto_bigint::{NonZero, U256};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{CkksEngine, Encoder, PolyRingU256, encoding::BigIntEncoder};

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 60; // should be larger than usual
// Large 256-bit modulus (example: a 256-bit prime)
const MODULUS_HEX: &str =
    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";

type Engine = CkksEngine<PolyRingU256<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS BigInt U256 Backend Demo");

    // Create U256 modulus from hex string
    let modulus_u256 = U256::from_be_hex(MODULUS_HEX);
    let modulus = NonZero::new(modulus_u256).expect("Modulus should not be zero");
    let encoder = BigIntEncoder::new(SCALE_BITS)?;

    println!("✅ Using 256-bit modulus: 0x{}", MODULUS_HEX);

    // Create CKKS Engine with BigInt U256 backend
    let engine = Engine::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_bigint_u256(modulus, SCALE_BITS)?;
    println!("✅ Engine configured with BigInt U256 backend");

    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    println!("✅ Secret and public keys generated with U256 arithmetic");

    // Input data (small vector to test)
    let values = vec![1.5, 2.5, 3.5, 4.5];
    println!("\n📊 Input data: {:?}", values);

    // Step 1: Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext = encoder.encode(&values, engine.context());
    println!("✅ Values encoded to plaintext with BigInt backend");

    // Step 2: Encrypt: Plaintext → Ciphertext
    println!("\n🔒 Encrypting plaintext...");
    let ciphertext = engine.encrypt(&plaintext, &public_key, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext using U256 operations");

    // Step 3: Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext, &secret_key);
    println!("✅ Ciphertext decrypted back to plaintext");

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
        println!(
            "🎉 Success! Full CKKS encrypt/decrypt pipeline works with BigInt U256!"
        );
    } else {
        println!("⚠️  Warning: Error is higher than expected");
    }

    // Demonstrate homomorphic addition with BigInt backend
    println!("\n➕ Testing homomorphic addition with U256 arithmetic...");

    let values2 = vec![0.5, 1.0, 1.5, 0.0];
    println!("Second input: {:?}", values2);

    let plaintext2 = encoder.encode(&values2, engine.context());
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);

    // Homomorphic addition
    let ciphertext_sum = Engine::add_ciphertexts(&ciphertext, &ciphertext2);
    let decrypted_sum = Engine::decrypt(&ciphertext_sum, &secret_key);
    let decoded_sum = encoder.decode(&decrypted_sum);

    // Expected result
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

    if add_max_error < 1e-3 {
        println!("🎉 Homomorphic addition successful with BigInt U256 backend!");
    } else {
        println!("⚠️  Warning: Addition error higher than expected");
    }

    println!("\n🔬 BigInt Backend Comparison:");
    println!("  • Naive backend uses single u64 modulus");
    println!("  • BigInt backend uses 256-bit modulus");
    println!("  • Both backends produce equivalent CKKS functionality");
    println!(
        "  • BigInt backend supports much larger moduli for enhanced security"
    );

    Ok(())
}
