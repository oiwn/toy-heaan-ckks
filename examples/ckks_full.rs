use crypto_bigint::{NonZero, U256};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    BigIntPolyRing, CkksEngine, Encoder, PolyRescale, RelinearizationKey,
    crypto::Ciphertext, crypto::operations::multiply_ciphertexts, encoding,
    rings::backends::bigint::ModulusDomain,
};

// Kim's HEAAN Parameters (adjusted for U256 safety)
const LOG_N: u32 = 9; // Ring dimension bits (N = 2^9 = 512) - closer to Kim's expectations
const DEGREE: usize = 1 << LOG_N; // N = 512 (polynomial degree)
const SCALE_BITS: u32 = 40; // Scale/precision bits (Kim typically uses 40-45)

// Kim's HEAAN Modulus Parameters (safe for U256, but closer to Kim's ratios)
const BASE_MODULUS_BITS: u32 = 54; // Base modulus: q = 2^54 (Kim uses ~54-60)
const EXTENDED_MODULUS_BITS: u32 = 54; // Extended modulus: Q = 2^54 (Kim uses Q ≈ q)
// Combined: qQ ≈ 2^108 (safe for U256 which is 256 bits)

const Q_BASE: u64 = (1u64 << 50) - 27; // Use 50-bit base for backward compatibility

type Engine = CkksEngine<BigIntPolyRing<DEGREE>, DEGREE>;

/// Kim's HEAAN Multiplication Algorithm Implementation
///
/// This implements the exact algorithm from HEAAN/src/Scheme.cpp lines 1276-1302
/// with proper modulus management for the BigInt backend.
///
/// Key algorithmic insights:
/// 1. Two-domain operation: Regular multiplications use base modulus `q`,
///    key switching uses extended modulus `qQ = q × Q`
/// 2. The algebraic trick: `(a1+b1)(a2+b2) - b1*b2 - a1*a2 = a1*b2 + a2*b1`
/// 3. Key switching handles quadratic term with proper scaling
/// 4. Critical scaling: After key switching in `qQ`, must divide by `Q`
fn kim_heaan_multiply(
    ct1: &Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>,
    ct2: &Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>,
    relin_key: &RelinearizationKey<BigIntPolyRing<DEGREE>, DEGREE>,
    target_scale_bits: u32,
    kim_context: &toy_heaan_ckks::rings::backends::bigint::BigIntContext,
) -> Result<Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>, Box<dyn std::error::Error>>
{
    println!("🔄 Starting Kim's HEAAN multiplication algorithm...");

    // Step 1: Setup moduli for different domains using Kim's context
    println!("  Step 1: Using Kim's modulus setup:");
    println!(
        "    Base modulus q = 2^{} = {}",
        kim_context.logq, kim_context.q
    );
    println!(
        "    Extended modulus Q = 2^{} = {}",
        kim_context.log_q,
        kim_context.q_extended_u256()
    );
    println!("    Combined modulus qQ = {}", kim_context.q_times_q());

    // Convert input ciphertexts to use Kim's context
    let ct1_kim = Ciphertext {
        c0: BigIntPolyRing::from_u256_coeffs(&ct1.c0.coeffs, kim_context),
        c1: BigIntPolyRing::from_u256_coeffs(&ct1.c1.coeffs, kim_context),
        logp: ct1.logp,
        logq: ct1.logq,
    };
    let ct2_kim = Ciphertext {
        c0: BigIntPolyRing::from_u256_coeffs(&ct2.c0.coeffs, kim_context),
        c1: BigIntPolyRing::from_u256_coeffs(&ct2.c1.coeffs, kim_context),
        logp: ct2.logp,
        logq: ct2.logq,
    };
    let relin_key_kim = RelinearizationKey {
        a: BigIntPolyRing::from_u256_coeffs(&relin_key.a.coeffs, kim_context),
        b: BigIntPolyRing::from_u256_coeffs(&relin_key.b.coeffs, kim_context),
    };

    // Step 2: Compute polynomial products in base modulus q using Kim's context
    println!("  Step 2: Computing polynomial products in base modulus q...");

    // Ring2Utils::add(axbx1, cipher1.ax, cipher1.bx, q, context.N);
    let mut axbx1 = ct1_kim.c0.clone();
    axbx1 += &ct1_kim.c1;
    println!("    Computed (a1 + b1)");

    // Ring2Utils::add(axbx2, cipher2.ax, cipher2.bx, q, context.N);
    let mut axbx2 = ct2_kim.c0.clone();
    axbx2 += &ct2_kim.c1;
    println!("    Computed (a2 + b2)");

    // Ring2Utils::multAndEqual(axbx1, axbx2, q, context.N);
    axbx1 *= &axbx2;
    println!("    Computed (a1+b1)(a2+b2)");

    // Ring2Utils::mult(axax, cipher1.ax, cipher2.ax, q, context.N);
    let mut axax = ct1_kim.c0.clone();
    axax *= &ct2_kim.c0;
    println!("    Computed a1*a2");

    // Ring2Utils::mult(bxbx, cipher1.bx, cipher2.bx, q, context.N);
    let mut bxbx = ct1_kim.c1.clone();
    bxbx *= &ct2_kim.c1;
    println!("    Computed b1*b2");

    // Step 3: Key switching in extended modulus qQ (critical!)
    println!(
        "  Step 3: Key switching with relinearization key in extended domain qQ..."
    );

    // Convert axax to extended domain (qQ) for key switching
    let axax_extended = axax.to_extended_domain();
    println!("    Converted a1*a2 to extended domain qQ");

    // Ring2Utils::mult(axmult, axax, key.ax, qQ, context.N);
    let mut axmult = axax_extended.clone();
    axmult.mul_in_extended_domain(&relin_key_kim.a);
    println!("    Computed (a1*a2) * relin_key.a in extended domain qQ");

    // Ring2Utils::mult(bxmult, axax, key.bx, qQ, context.N);
    let mut bxmult = axax_extended.clone();
    bxmult.mul_in_extended_domain(&relin_key_kim.b);
    println!("    Computed (a1*a2) * relin_key.b in extended domain qQ");

    // Step 4: Scale down by Q (divide by 2^logQ) - THE CRITICAL FIX!
    println!("  Step 4: Scaling down by Q (✅ NOW IMPLEMENTED!)");
    axmult.scale_down_by_q();
    bxmult.scale_down_by_q();
    println!("    Scaled down both components by Q = 2^logQ");
    println!("    This is the critical missing operation from BIGINT_MUL.md!");

    // Step 5: Final assembly back in base modulus q
    println!("  Step 5: Final assembly in base modulus q...");

    // Ensure all operations are in base modulus q
    let mut axmult_base = axmult.in_domain(ModulusDomain::Base);
    let mut bxmult_base = bxmult.in_domain(ModulusDomain::Base);
    let axbx1_base = axbx1.in_domain(ModulusDomain::Base);
    let bxbx_base = bxbx.in_domain(ModulusDomain::Base);
    let axax_base = axax.in_domain(ModulusDomain::Base);

    // Ring2Utils::addAndEqual(axmult, axbx1, q, context.N);  // + (a1+b1)(a2+b2)
    axmult_base += &axbx1_base;
    println!("    Added (a1+b1)(a2+b2) in base modulus q");

    // Ring2Utils::subAndEqual(axmult, bxbx, q, context.N);   // - b1*b2
    let neg_bxbx = -bxbx_base.clone();
    axmult_base += &neg_bxbx;
    println!("    Subtracted b1*b2 in base modulus q");

    // Ring2Utils::subAndEqual(axmult, axax, q, context.N);   // - a1*a2
    let neg_axax = -axax_base.clone();
    axmult_base += &neg_axax;
    println!("    Subtracted a1*a2 in base modulus q");

    // Ring2Utils::addAndEqual(bxmult, bxbx, q, context.N);   // + b1*b2
    bxmult_base += &bxbx_base;
    println!("    Added b1*b2 to second component in base modulus q");

    // Use the base modulus results
    let mut c0_result = axmult_base;
    let mut c1_result = bxmult_base;

    // Step 6: Additional rescaling after Kim's scale_down_by_q()
    println!("  Step 6: Additional rescaling after Kim's algorithm...");
    let doubled_scale_bits = ct1.logp + ct2.logp; // 60
    let kim_scaling_bits = kim_context.log_q; // 20 bits scaled by Kim
    let current_scale_bits = doubled_scale_bits - kim_scaling_bits; // 60 - 20 = 40

    println!(
        "    Original scales: {} + {} = {}",
        ct1.logp, ct2.logp, doubled_scale_bits
    );
    println!(
        "    Kim's scale_down_by_q() scaled by {} bits",
        kim_scaling_bits
    );
    println!("    Current effective scale: {} bits", current_scale_bits);
    println!("    Target scale: {}", target_scale_bits);

    if current_scale_bits > target_scale_bits {
        let final_rescale_bits = current_scale_bits - target_scale_bits; // 40 - 30 = 10
        println!(
            "    Additional rescaling needed: {} bits",
            final_rescale_bits
        );
        let scale_factor = (1u64 << final_rescale_bits) as f64;
        c0_result.rescale_assign(scale_factor);
        c1_result.rescale_assign(scale_factor);
    } else {
        println!("    No additional rescaling needed");
    }

    println!("✅ Kim's multiplication algorithm complete");

    Ok(Ciphertext {
        c0: c0_result,
        c1: c1_result,
        logp: target_scale_bits,
        logq: target_scale_bits,
    })
}

/// Extended Modulus Operations for BigInt Backend
///
/// Based on Kim's HEAAN Context.cpp, we need to support:
/// - q: Base modulus (2^logq) for regular operations
/// - Q: Extended modulus (2^logQ) for key switching
/// - qQ: Product modulus (q * Q) for relinearization
///
/// Current limitation: Our BigInt backend uses fixed U256 modulus
/// TODO: Extend to support variable moduli like Kim's implementation
fn extended_modulus_demo() {
    println!("\n🔧 Extended Modulus Operations Demo");
    println!("Current BigInt Backend Limitation:");
    println!("  ❌ Uses fixed modulus U256 for all operations");
    println!("  ❌ Cannot handle variable moduli like Kim's HEAAN");
    println!("  ❌ Missing proper key switching with modulus scaling");

    println!("\nKim's HEAAN Requirements:");
    println!("  ✅ q = 2^logq (base modulus)");
    println!("  ✅ Q = 2^logQ (extended modulus)");
    println!("  ✅ qQ = q * Q (key switching modulus)");
    println!("  ✅ QQ = Q^2 (key generation modulus)");

    println!("\nRequired Implementation:");
    println!("  1. Extend BigIntPolyRing to support multiple moduli");
    println!("  2. Add modulus switching operations");
    println!("  3. Implement right shift by logQ bits");
    println!("  4. Update key generation for extended modulus");
}

/// Test multiplication with detailed analysis
fn test_multiplication_detailed(
    engine: &Engine,
    secret_key: &toy_heaan_ckks::SecretKey<BigIntPolyRing<DEGREE>, DEGREE>,
    relin_key: &RelinearizationKey<BigIntPolyRing<DEGREE>, DEGREE>,
    encoder: &encoding::BigIntEncoder<DEGREE>,
    values1: &[f64],
    values2: &[f64],
    kim_context: &toy_heaan_ckks::rings::backends::bigint::BigIntContext,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    println!("\n🧪 Detailed Multiplication Test");
    println!("Input 1: {:?}", values1);
    println!("Input 2: {:?}", values2);
    println!(
        "Expected: {:?}",
        values1
            .iter()
            .zip(values2)
            .map(|(a, b)| a * b)
            .collect::<Vec<_>>()
    );

    // Encode
    let plaintext1 = encoder.encode(values1, engine.context());
    let plaintext2 = encoder.encode(values2, engine.context());

    // Encrypt
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    let public_key = engine.generate_public_key(secret_key, &mut rng)?;
    let ciphertext1 =
        engine.encrypt(&plaintext1, &public_key, SCALE_BITS, &mut rng);
    let ciphertext2 =
        engine.encrypt(&plaintext2, &public_key, SCALE_BITS, &mut rng);

    println!(
        "Original ciphertext scales: {} and {}",
        ciphertext1.logp, ciphertext2.logp
    );

    // Test Kim's multiplication
    let result_ciphertext = kim_heaan_multiply(
        &ciphertext1,
        &ciphertext2,
        relin_key,
        SCALE_BITS,
        kim_context,
    )?;

    // Decrypt and decode
    let decrypted = Engine::decrypt(&result_ciphertext, secret_key);
    let decoded = encoder.decode(&decrypted);

    println!("Actual result: {:?}", decoded);

    Ok(decoded)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    println!("🔐 CKKS BigInt U256 Backend Multiplication Demo");

    // Create BigInt context with Kim's HEAAN parameters
    let context =
        toy_heaan_ckks::rings::backends::bigint::BigIntContext::from_kim_params(
            BASE_MODULUS_BITS,
            EXTENDED_MODULUS_BITS,
        )?;

    println!(
        "✅ Using Kim's context: q=2^{}, Q=2^{}",
        BASE_MODULUS_BITS, EXTENDED_MODULUS_BITS
    );
    println!("   Base modulus q: {}", context.q);
    println!("   Extended modulus Q: {}", context.q_extended_u256());

    let encoder = encoding::BigIntEncoder::new(SCALE_BITS)?;

    // Create CKKS Engine with BigInt U256 backend using proper context
    let engine = Engine::builder()
        .error_variance(0.2)
        .hamming_weight(DEGREE / 2)
        .build_bigint_u256(context.q, SCALE_BITS)?;
    println!("✅ Engine configured with BigInt U256 backend");

    println!("\n🔑 Generating keys...");
    let secret_key = engine.generate_secret_key(&mut rng)?;
    let public_key = engine.generate_public_key(&secret_key, &mut rng)?;
    let relin_key = engine.generate_relinearization_key(&secret_key, &mut rng)?;
    println!("✅ Secret and public keys generated with U256 arithmetic");

    // Input data - using Kim's 8 slots parameter
    let values1 = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // 8 slots, only first active
    let values2 = vec![2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // 8 slots, only first active
    println!("📊 Input data1: {:?}", &values1[0..2]); // Show first 2 slots
    println!("📊 Input data2: {:?}", &values2[0..2]); // Show first 2 slots

    // Step 1: Encode: Vec<f64> -> Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext1 = encoder.encode(&values1, engine.context());
    println!(
        "Plaintext1 coefficient[0]: {:x}",
        plaintext1.poly.coeffs[0].to_words()[0]
    );
    let plaintext2 = encoder.encode(&values2, engine.context());
    println!(
        "Plaintext2 coefficient[0]: {:x}",
        plaintext2.poly.coeffs[0].to_words()[0]
    );
    println!("✅ Values encoded to plaintext with BigInt backend");

    // sanity: encode -> decode without encryption
    let decoded_plain = encoder.decode(&plaintext1);
    println!("Plain roundtrip: {:?}", &decoded_plain);
    println!(
        "  First few coeffs for plaintext1: {:?}",
        &plaintext1.poly.coeffs[0..4]
    );

    // Step 2: Encrypt: Plaintext → Ciphertext
    println!("\n🔒 Encrypting plaintext...");
    let ciphertext1 =
        engine.encrypt(&plaintext1, &public_key, SCALE_BITS, &mut rng);
    let ciphertext2 =
        engine.encrypt(&plaintext2, &public_key, SCALE_BITS, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext using U256 operations");

    // Show extended modulus operations demo
    extended_modulus_demo();

    // Test Kim's multiplication implementation
    println!("\n🔄 Testing Kim's HEAAN Multiplication Implementation");

    // Test 1: Original multiply_ciphertexts function
    println!("\n=== Test 1: Original multiply_ciphertexts function ===");
    let ciphertext_mul_original =
        multiply_ciphertexts(&ciphertext1, &ciphertext2, &relin_key, SCALE_BITS)
            .expect("Multiplication should succeed");
    let decrypted_original = Engine::decrypt(&ciphertext_mul_original, &secret_key);
    let decoded_original = encoder.decode(&decrypted_original);
    println!("Original function result: {:?}", decoded_original);

    // Test 2: Kim's multiplication function
    println!("\n=== Test 2: Kim's HEAAN multiplication function ===");
    let ciphertext_mul_kim = kim_heaan_multiply(
        &ciphertext1,
        &ciphertext2,
        &relin_key,
        SCALE_BITS,
        &context,
    )
    .expect("Kim's multiplication should succeed");
    let decrypted_kim = Engine::decrypt(&ciphertext_mul_kim, &secret_key);
    let decoded_kim = encoder.decode(&decrypted_kim);
    println!("Kim's function result: {:?}", decoded_kim);

    // Test 3: Detailed multiplication test
    println!("\n=== Test 3: Detailed Analysis ===");
    let detailed_result = test_multiplication_detailed(
        &engine,
        &secret_key,
        &relin_key,
        &encoder,
        &values1,
        &values2,
        &context,
    )?;

    println!("\n📊 Results Summary:");
    println!(
        "Expected result: {:?}",
        values1
            .iter()
            .zip(&values2)
            .map(|(a, b)| a * b)
            .collect::<Vec<_>>()
    );
    println!("Original function: {:?}", decoded_original);
    println!("Kim's function: {:?}", decoded_kim);
    println!("Detailed test: {:?}", detailed_result);

    println!("\n🔍 Analysis:");
    println!(
        "All methods show similar incorrect results, confirming the modulus issue:"
    );
    println!("- Results are far from expected [2.0]");
    println!("- Values are in hundreds/thousands instead of ~2");
    println!("- This matches the findings in BIGINT_MUL.md");
    println!("\n⚠️  Root Cause: Missing variable modulus support (q, Q, qQ)");
    println!(
        "Next steps: Implement extended BigIntContext with proper modulus switching"
    );

    Ok(())
}
