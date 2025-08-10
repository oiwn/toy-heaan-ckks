use crypto_bigint::{NonZero, U256};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::rings::backends::bigint::rescale_ciphertext_u256_inplace;
use toy_heaan_ckks::{CkksEngine, Encoder, PolyRingU256, encoding};

const DEGREE: usize = 8;
// const SCALE_BITS: u32 = 60; // should be larger than usual
const SCALE_BITS: u32 = 40;
const Q50: u64 = (1u64 << 50) - 27;
// const Q50_U256: U256 = U256::from_u64((1u64 << 50) - 27);
// const Q128: u128 = (1u128 << 100) - 3;

type Engine = CkksEngine<PolyRingU256<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS BigInt U256 Backend Multiplication Demo");

    let modulus_u256 = U256::from_u128(Q50 as u128);
    let modulus = NonZero::new(modulus_u256).expect("Modulus should not be zero");
    let encoder = encoding::BigIntEncoder::new(SCALE_BITS)?;

    println!("✅ Using 256-bit modulus: {}", modulus_u256);

    // Create CKKS Engine with BigInt U256 backend
    let engine = Engine::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_bigint_u256(modulus, SCALE_BITS)?;
    println!("✅ Engine configured with BigInt U256 backend");

    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    let relin_key = engine.generate_relinearization_key(&secret_key, &mut rng)?;
    println!("✅ Secret and public keys generated with U256 arithmetic");

    // Input data
    let values1 = vec![1.0, 2.0, 3.0, 4.0];
    let values2 = vec![5.0, 6.0, 7.0, 0.0];
    println!("\n📊 Input data1: {values1:?}");
    println!("\n📊 Input data2: {values2:?}");

    // Step 1: Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext1 = encoder.encode(&values1, engine.context());
    let plaintext2 = encoder.encode(&values2, engine.context());
    println!("✅ Values encoded to plaintext with BigInt backend");

    // sanity: encode -> decode without encryption
    let decoded_plain = encoder.decode(&plaintext1);
    println!("Plain roundtrip: {:?}", &decoded_plain[..values1.len()]);

    // Step 2: Encrypt: Plaintext → Ciphertext
    println!("\n🔒 Encrypting plaintext...");
    let ciphertext1 = engine.encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext using U256 operations");

    let mut ciphertext_mul =
        Engine::mul_ciphertexts(&ciphertext1, &ciphertext2, &relin_key);

    println!("Before rescale:");
    println!("  ct1.scale: {}", ciphertext1.scale);
    println!("  ct2.scale: {}", ciphertext2.scale);
    println!("  mul.scale: {}", ciphertext_mul.scale);
    println!(
        "  Expected Δ²: {}",
        (1u64 << SCALE_BITS) as f64 * (1u64 << SCALE_BITS) as f64
    );

    rescale_ciphertext_u256_inplace(&mut ciphertext_mul, SCALE_BITS);

    println!("After rescale:");
    println!("  mul.scale: {}", ciphertext_mul.scale);
    println!("  Expected Δ: {}", (1u64 << SCALE_BITS) as f64);

    // Step 3: Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext_mul, &secret_key);
    println!("✅ Ciphertext decrypted back to plaintext");

    println!("Debug scales:");
    println!("  encoder.delta(): {}", (1u64 << SCALE_BITS) as f64);
    println!("  decrypted_plaintext.scale: {}", decrypted_plaintext.scale);

    let coeffs = decrypted_plaintext.poly.coefficients();
    println!("  First few coeffs: {:?}", &coeffs[0..4]);

    // Step 4: Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding back to floating-point...");
    let decoded_values = encoder.decode(&decrypted_plaintext);
    println!("✅ Plaintext1 decoded: {:?}", decoded_values);

    /* // Display results
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
    } */

    Ok(())
}
