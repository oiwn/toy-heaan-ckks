use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    CkksEngine, Encoder, NaivePolyRing, crypto::operations::multiply_ciphertexts,
    encoding::RustFftEncoder,
};

const DEGREE: usize = 8;
const MODULUS: u64 = (1u64 << 50) - 27; // Large modulus for headroom
const SCALE_BITS: u32 = 10; // 2^10 = 1024

type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

    println!("🎯 CKKS Multiplication Test with 4-element vectors");
    println!("   Using NaivePolyRing<{}> with {}-bit modulus", DEGREE, 50);
    println!("   Scale factor: 2^{} = {}", SCALE_BITS, 1u64 << SCALE_BITS);

    // Create engine
    let engine = Engine::builder()
        .error_variance(0.2) // Same as working example
        .hamming_weight(DEGREE / 2)
        .build_naive(MODULUS, SCALE_BITS)?;

    // Generate keys
    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    let relin_key = engine.generate_relinearization_key(&secret_key, &mut rng)?;

    // Create encoder for 4 slots (DEGREE/2 = 4 complex numbers)
    let encoder = RustFftEncoder::<DEGREE>::new(SCALE_BITS)?;
    let num_slots = DEGREE / 2; // 4 slots for complex numbers

    // Create simple test vectors that are easy to verify
    println!("\n📝 Creating test vectors...");
    let vector1 = vec![1.0, 2.0, 3.0, 4.0]; // [1, 2, 3, 4]
    let vector2 = vec![2.0, 1.0, 1.0, 2.0]; // [2, 1, 1, 2]

    println!("   Vector 1: {:?}", vector1);
    println!("   Vector 2: {:?}", vector2);
    println!(
        "   Expected product: {:?}",
        vector1
            .iter()
            .zip(&vector2)
            .map(|(a, b)| a * b)
            .collect::<Vec<_>>()
    );

    // Encode to plaintexts
    println!("\n🔧 Encoding to plaintexts...");
    let plaintext1: toy_heaan_ckks::Plaintext<NaivePolyRing<DEGREE>, DEGREE> =
        encoder.encode(&vector1, engine.context());
    let plaintext2: toy_heaan_ckks::Plaintext<NaivePolyRing<DEGREE>, DEGREE> =
        encoder.encode(&vector2, engine.context());

    println!(
        "   Plaintext 1 coefficients: {:?}",
        &plaintext1.poly.coeffs[0..4]
    );
    println!(
        "   Plaintext 2 coefficients: {:?}",
        &plaintext2.poly.coeffs[0..4]
    );

    // Verify encoding by decoding
    let decoded1 = encoder.decode(&plaintext1);
    let decoded2 = encoder.decode(&plaintext2);
    println!("   Decoded vector 1: {:?}", &decoded1[0..num_slots]);
    println!("   Decoded vector 2: {:?}", &decoded2[0..num_slots]);

    // Encrypt
    println!("\n🔐 Encrypting plaintexts...");
    let ciphertext1 = engine.encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);

    // Verify encryption by decryption
    println!("\n🔍 Verifying encryption...");
    let decrypt1 = Engine::decrypt(&ciphertext1, &secret_key);
    let decrypt2 = Engine::decrypt(&ciphertext2, &secret_key);
    let verify1 = encoder.decode(&decrypt1);
    let verify2 = encoder.decode(&decrypt2);
    println!("   Decrypted vector 1: {:?}", &verify1[0..num_slots]);
    println!("   Decrypted vector 2: {:?}", &verify2[0..num_slots]);

    // Check encryption accuracy
    let error1: f64 = vector1
        .iter()
        .zip(&verify1[0..num_slots])
        .map(|(a, b)| (a - b).abs())
        .sum::<f64>()
        / vector1.len() as f64;
    let error2: f64 = vector2
        .iter()
        .zip(&verify2[0..num_slots])
        .map(|(a, b)| (a - b).abs())
        .sum::<f64>()
        / vector2.len() as f64;
    println!("   Average encryption error 1: {:.6}", error1);
    println!("   Average encryption error 2: {:.6}", error2);

    if error1 > 0.1 || error2 > 0.1 {
        println!("❌ Encryption error too high, stopping");
        return Ok(());
    }

    // Perform homomorphic multiplication
    println!("\n✖️  Performing homomorphic multiplication...");
    let ciphertext_product =
        multiply_ciphertexts(&ciphertext1, &ciphertext2, &relin_key, SCALE_BITS)?;

    println!("   Original scale_bits: {}", ciphertext1.scale_bits);
    println!("   Product scale_bits: {}", ciphertext_product.scale_bits);

    // Decrypt the result
    println!("\n🔓 Decrypting multiplication result...");
    let decrypted_product = Engine::decrypt(&ciphertext_product, &secret_key);

    // Decode the result
    let result_full = encoder.decode(&decrypted_product);
    let result_vector = &result_full[0..num_slots];

    // Calculate expected result
    let expected_product: Vec<f64> =
        vector1.iter().zip(&vector2).map(|(a, b)| a * b).collect();

    println!("\n📊 Results:");
    println!("   Input vector 1:    {:?}", vector1);
    println!("   Input vector 2:    {:?}", vector2);
    println!("   Expected product:  {:?}", expected_product);
    println!("   Computed product:  {:?}", result_vector);

    // Calculate errors
    let errors: Vec<f64> = expected_product
        .iter()
        .zip(result_vector)
        .map(|(exp, got)| (exp - got).abs())
        .collect();
    let avg_error = errors.iter().sum::<f64>() / errors.len() as f64;
    let max_error = errors.iter().fold(0.0f64, |a, &b| a.max(b));

    println!("   Element-wise errors: {:?}", errors);
    println!("   Average error: {:.6}", avg_error);
    println!("   Maximum error: {:.6}", max_error);

    if max_error < 0.1 {
        println!("\n🎉 SUCCESS! Homomorphic multiplication works correctly!");
        println!("   ✅ All elements computed within acceptable error bounds");
    } else {
        println!("\n⚠️  Some elements have higher error than expected");
        for (i, ((exp, got), err)) in expected_product
            .iter()
            .zip(result_vector)
            .zip(&errors)
            .enumerate()
        {
            if *err > 0.1 {
                println!(
                    "   Element {}: expected {:.3}, got {:.3}, error {:.6}",
                    i, exp, got, err
                );
            }
        }
    }

    // Show raw polynomial coefficients for debugging
    println!("\n🔬 Debug information:");
    println!(
        "   Raw decrypted coefficients: {:?}",
        &decrypted_product.poly.coeffs[0..4]
    );
    println!(
        "   Scale factor used: 2^{} = {}",
        SCALE_BITS,
        1u64 << SCALE_BITS
    );

    Ok(())
}
