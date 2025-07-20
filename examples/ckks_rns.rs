use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::sync::Arc;
use toy_heaan_ckks::{CkksEngine, Encoder, Plaintext, RustFftEncoder};

// Import the RNS backend components
use toy_heaan_ckks::rings::backends::rns::{RnsBasisBuilder, RnsPolyRing};

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 40;

type Engine = CkksEngine<RnsPolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS RNS Backend Demo");

    // Step 1: Create RNS basis with multiple primes
    println!("\n⚙️  Creating RNS basis...");
    let rns_basis = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_prime_bits(vec![17, 19, 23]) // Three primes of different bit sizes
            .build()?,
    );

    println!(
        "✅ RNS basis created with {} primes:",
        rns_basis.channel_count()
    );
    for (i, &prime) in rns_basis.primes().iter().enumerate() {
        println!("   Prime {}: {} ({} bits)", i, prime, prime.ilog2() + 1);
    }

    // Step 2: Create CKKS Engine with RNS context using builder
    println!("\n🏗️  Building CKKS engine with RNS backend...");
    let engine = Engine::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_rns(rns_basis.clone(), SCALE_BITS)?;

    let encoder = RustFftEncoder::new(engine.params.scale_bits)?;
    println!("✅ Engine configured with RNS backend");

    // Step 3: Generate keys
    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    println!("✅ Secret and public keys generated");

    // Print some debug info about the keys
    println!("   Secret key channels: {}", secret_key.poly.channels());
    println!("   Public key a channels: {}", public_key.a.channels());
    println!("   Public key b channels: {}", public_key.b.channels());

    // Step 4: Prepare test data
    let values1 = vec![1.5, 2.5, 3.5];
    let values2 = vec![0.5, 1.0, 1.5];
    println!("\n📊 Input data:");
    println!("   Values 1: {:?}", values1);
    println!("   Values 2: {:?}", values2);

    // Calculate expected results for verification
    let expected_sum: Vec<f64> =
        values1.iter().zip(&values2).map(|(a, b)| a + b).collect();
    println!("   Expected sum: {:?}", expected_sum);

    // Step 5: Encode and encrypt
    println!("\n🔢 Encoding and encrypting...");
    let plaintext1: Plaintext<RnsPolyRing<DEGREE>, DEGREE> =
        encoder.encode(&values1, engine.context());
    let plaintext2: Plaintext<RnsPolyRing<DEGREE>, DEGREE> =
        encoder.encode(&values2, engine.context());

    println!("✅ Values encoded to plaintexts");
    println!("   Plaintext 1 channels: {}", plaintext1.poly.channels());
    println!("   Plaintext 2 channels: {}", plaintext2.poly.channels());

    let ciphertext1 = engine.encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);
    println!("✅ Plaintexts encrypted to ciphertexts");

    // Step 6: Homomorphic addition
    println!("\n➕ Performing homomorphic addition...");
    let ciphertext_sum = Engine::add_ciphertexts(&ciphertext1, &ciphertext2);
    println!("✅ Homomorphic addition completed");
    println!(
        "   Result ciphertext c0 channels: {}",
        ciphertext_sum.c0.channels()
    );
    println!(
        "   Result ciphertext c1 channels: {}",
        ciphertext_sum.c1.channels()
    );

    // Step 7: Decrypt and decode
    println!("\n🔓 Decrypting and decoding result...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext_sum, &secret_key);
    let result = encoder.decode(&decrypted_plaintext);

    println!("✅ Decryption and decoding completed");
    println!("   Decrypted values: {:?}", result);

    // Step 8: Verify results
    println!("\n📋 Verification:");
    let max_error = expected_sum
        .iter()
        .zip(result.iter())
        .map(|(exp, act)| (exp - act).abs())
        .fold(0.0, f64::max);

    println!("   Expected: {:?}", expected_sum);
    println!("   Computed: {:?}", result);
    println!("   Maximum error: {:.2e}", max_error);

    if max_error < 1e-6 {
        println!("✅ Success! Error within acceptable bounds");
    } else {
        println!("❌ Warning: Error is higher than expected");
    }

    // Step 9: Demonstrate RNS-specific features
    println!("\n🔬 RNS-specific demonstrations:");

    // Show coefficient reconstruction via CRT
    println!("   Coefficient reconstruction example (first 4 coefficients):");
    for i in 0..4.min(DEGREE) {
        let reconstructed = ciphertext_sum.c0.coefficient_to_u64(i);
        println!("     c0[{}] = {} (reconstructed via CRT)", i, reconstructed);
    }

    // Show residues for each prime
    println!("   Residue representation for c0[0]:");
    for (channel_idx, &prime) in rns_basis.primes().iter().enumerate() {
        let residue = ciphertext_sum.c0.coefficients[channel_idx][0];
        println!("     c0[0] mod {} = {}", prime, residue);
    }

    // Test RNS polynomial arithmetic directly
    println!("\n🧮 Testing RNS polynomial arithmetic:");
    let poly_a: RnsPolyRing<DEGREE> =
        RnsPolyRing::from_i64_slice(&[1, 2, 3, 4, 0, 0, 0, 0], rns_basis.clone());
    let poly_b: RnsPolyRing<DEGREE> =
        RnsPolyRing::from_i64_slice(&[2, 3, 4, 5, 0, 0, 0, 0], rns_basis.clone());

    println!("   poly_a coefficients: {:?}", poly_a.to_i64_coefficients());
    println!("   poly_b coefficients: {:?}", poly_b.to_i64_coefficients());

    let poly_sum = &poly_a + &poly_b;
    println!("   poly_a + poly_b = {:?}", poly_sum.to_i64_coefficients());

    let mut poly_product = poly_a.clone();
    poly_product *= &poly_b;
    println!(
        "   poly_a * poly_b = {:?}",
        poly_product.to_i64_coefficients()
    );

    println!("\n🎉 RNS backend demonstration completed successfully!");

    Ok(())
}
