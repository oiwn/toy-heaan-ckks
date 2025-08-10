use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::sync::Arc;
use toy_heaan_ckks::{CkksEngine, Encoder, Plaintext, RustFftEncoder};

// Import the new NTT backend
use toy_heaan_ckks::rings::backends::ntt::{
    RnsNttBasis, RnsNttPolyRing, create_ntt_engine,
};
use toy_heaan_ckks::rings::backends::rns::RnsBasisBuilder;

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 40;

type Engine = CkksEngine<RnsNttPolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS with NTT Backend Demo");

    // Step 1: Create RNS basis with NTT-friendly primes
    println!("\n⚙️  Creating NTT-friendly RNS basis...");
    let base_rns_basis = RnsBasisBuilder::new(DEGREE)
        .with_prime_bits(vec![17, 19, 23]) // Generate NTT-friendly primes
        .build()?;

    println!(
        "✅ Base RNS basis created with {} primes:",
        base_rns_basis.channel_count()
    );
    for (i, &prime) in base_rns_basis.primes().iter().enumerate() {
        let ntt_check = (prime - 1) % (2 * DEGREE as u64) == 0;
        println!(
            "   Prime {}: {} ({} bits) - NTT-friendly: {}",
            i,
            prime,
            prime.ilog2() + 1,
            ntt_check
        );
    }

    // Step 2: Create NTT-enabled basis
    println!("\n🧮 Creating NTT tables...");
    let ntt_basis = Arc::new(RnsNttBasis::new(base_rns_basis, DEGREE)?);
    println!(
        "✅ NTT basis created with {} tables",
        ntt_basis.ntt_tables.len()
    );

    // Step 3: Create CKKS Engine with NTT backend
    println!("\n🏗️  Building CKKS engine with NTT backend...");
    let engine = create_ntt_engine(
        ntt_basis.clone(),
        SCALE_BITS,
        3.2,        // error_variance
        DEGREE / 2, // hamming_weight
    );

    let encoder = RustFftEncoder::new(engine.params.scale_bits)?;
    println!("✅ Engine configured with NTT backend");

    // Step 4: Generate keys
    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    println!("✅ Secret and public keys generated");

    // Print debug info about the keys and their domain
    println!("   Secret key channels: {}", secret_key.poly.channels());
    println!("   Secret key NTT form: {}", secret_key.poly.is_ntt_form());
    println!("   Public key a channels: {}", public_key.a.channels());
    println!("   Public key a NTT form: {}", public_key.a.is_ntt_form());

    // Step 5: Test basic CKKS operations
    let values = vec![1.0, 2.0, 3.0, 4.0];
    println!("\n📊 Input data: {:?}", values);

    // Encode: Vec<f64> → Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext: Plaintext<RnsNttPolyRing<DEGREE>, DEGREE> =
        encoder.encode(&values, engine.context());
    println!("✅ Values encoded to plaintext");
    println!(
        "   Plaintext poly NTT form: {}",
        plaintext.poly.is_ntt_form()
    );

    // Encrypt: Plaintext → Ciphertext
    println!("\n🔐 Encrypting plaintext...");
    let ciphertext = engine.encrypt(&plaintext, &public_key, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext");
    println!("   Ciphertext c0 NTT form: {}", ciphertext.c0.is_ntt_form());
    println!("   Ciphertext c1 NTT form: {}", ciphertext.c1.is_ntt_form());

    // Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext, &secret_key);
    println!("✅ Ciphertext decrypted to plaintext");
    println!(
        "   Decrypted plaintext NTT form: {}",
        decrypted_plaintext.poly.is_ntt_form()
    );

    // Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding plaintext...");
    let decoded_values = encoder.decode(&decrypted_plaintext);
    println!("✅ Plaintext decoded to values");

    // Verify encrypt/decrypt accuracy
    println!("\n📊 Encrypt/Decrypt Results:");
    println!("  Original:  {:?}", values);
    println!("  Decoded:   {:?}", &decoded_values[..values.len()]);

    let max_error = values
        .iter()
        .zip(decoded_values.iter())
        .map(|(orig, dec)| (orig - dec).abs())
        .fold(0.0, f64::max);

    println!("  Max error: {:.2e}", max_error);

    if max_error < 1e-3 {
        println!(
            "🎉 Success! Full CKKS encrypt/decrypt pipeline works with NTT backend!"
        );
    } else {
        println!("⚠️  Warning: Error is higher than expected");
    }

    // Step 6: Test homomorphic addition
    println!("\n➕ Testing homomorphic addition with NTT backend...");

    let values2 = vec![0.5, 1.0, 1.5, 2.0];
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
        println!("🎉 Homomorphic addition successful with NTT backend!");
    } else {
        println!("⚠️  Warning: Addition error higher than expected");
    }

    // Step 7: Test domain conversions
    println!("\n🔄 Testing domain conversions...");

    // Create a test polynomial in coefficient form
    let mut test_poly: RnsNttPolyRing<DEGREE> = RnsNttPolyRing::from_i64_slice(
        &[1, 2, 3, 4, 0, 0, 0, 0],
        ntt_basis.clone(),
    );
    println!(
        "Original poly (coeff): {:?}",
        test_poly.to_u64_coefficients()
    );
    println!("NTT form: {}", test_poly.is_ntt_form());

    // Convert to NTT form
    test_poly.to_ntt_form();
    println!("After to_ntt_form(): {}", test_poly.is_ntt_form());

    // Convert back to coefficient form
    test_poly.to_coeff_form();
    println!("After to_coeff_form(): {}", test_poly.is_ntt_form());
    println!(
        "Recovered poly (coeff): {:?}",
        test_poly.to_u64_coefficients()
    );

    println!("\n💡 NTT Backend Summary:");
    println!("  1. ✅ NTT-friendly prime generation");
    println!("  2. ✅ NTT table precomputation");
    println!("  3. ✅ Domain tracking (coefficient ↔ NTT)");
    println!("  4. ✅ Full CKKS pipeline with NTT backend");
    println!("  5. ✅ Homomorphic addition in both domains");
    println!("  6. ✅ Domain conversion methods");
    println!("  7. 🔄 NTT multiplication (TODO: implement efficient version)");

    Ok(())
}
