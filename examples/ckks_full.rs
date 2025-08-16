use crypto_bigint::{NonZero, U256, Zero};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    BigIntPolyRing, Ciphertext, CkksEngine, Encoder, RelinearizationKey, encoding,
};

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 40;
const Q50: u64 = (1u64 << 50) - 27;

type Engine = CkksEngine<BigIntPolyRing<DEGREE>, DEGREE>;

/// Rescale ciphertext by dividing coefficients by 2^k (Kim's approach)
/// This follows Kim's HEAAN implementation pattern of coefficient rescaling
fn rescale_ciphertext(ct: &mut Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>, k: u32) {
    let divisor = U256::ONE << k; // 2^k
    let modulus = ct.c0.modulus().get();

    println!("Rescaling by 2^{} = {}", k, divisor);
    println!("Scale before rescaling: {}", ct.scale_bits);

    for i in 0..DEGREE {
        // For c0 coefficient
        let c0 = ct.c0.coeffs[i];
        ct.c0.coeffs[i] = rescale_coefficient(c0, modulus, divisor);

        // For c1 coefficient
        let c1 = ct.c1.coeffs[i];
        ct.c1.coeffs[i] = rescale_coefficient(c1, modulus, divisor);
    }

    // Update scale
    println!(
        "Scale after rescaling: {} - {} = {}",
        ct.scale_bits,
        k,
        ct.scale_bits - k
    );
    ct.scale_bits -= k;
}

/// Rescale a single coefficient: divide by divisor with proper modular arithmetic
/// Similar to Kim's scaleDownToReal but for modular integers
fn rescale_coefficient(coeff: U256, modulus: U256, divisor: U256) -> U256 {
    // Handle signed representation: map to centered interval
    let half_mod = modulus >> 1;

    let (is_negative, abs_val) = if coeff > half_mod {
        // Negative number in modular representation
        (true, modulus.wrapping_sub(&coeff))
    } else {
        // Positive number
        (false, coeff)
    };

    // Perform division with rounding (Kim's approach with RoundToZZ equivalent)
    let half_divisor = divisor >> 1;
    let rounded_div = abs_val.saturating_add(&half_divisor) / divisor;

    // Map back to modular representation
    if is_negative && !bool::from(rounded_div.is_zero()) {
        modulus.wrapping_sub(&rounded_div)
    } else {
        rounded_div % modulus
    }
}

// Debug multiplication to see intermediate coefficients
pub fn debug_mul_ciphertexts(
    ct1: &Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>,
    ct2: &Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>,
    relin_key: &RelinearizationKey<BigIntPolyRing<DEGREE>, DEGREE>,
) -> Ciphertext<BigIntPolyRing<DEGREE>, DEGREE> {
    println!("=== MULTIPLICATION DEBUG ===");

    // Show input coefficients
    println!("Input ct1.c0[0]: {:x}", ct1.c0.coeffs[0].to_words()[0]);
    println!("Input ct1.c1[0]: {:x}", ct1.c1.coeffs[0].to_words()[0]);
    println!("Input ct2.c0[0]: {:x}", ct2.c0.coeffs[0].to_words()[0]);
    println!("Input ct2.c1[0]: {:x}", ct2.c1.coeffs[0].to_words()[0]);

    // Step 1: Raw multiplication
    // d0 = c0 * c0'
    let mut d0 = ct1.c0.clone();
    d0 *= &ct2.c0;
    println!("d0 = c0*c0' [0]: {:x}", d0.coeffs[0].to_words()[0]);

    // d1 = c0*c1' + c1*c0'
    let mut d1_part1 = ct1.c0.clone();
    d1_part1 *= &ct2.c1;
    let mut d1_part2 = ct1.c1.clone();
    d1_part2 *= &ct2.c0;
    let mut d1 = d1_part1;
    d1 += &d1_part2;
    println!("d1 = c0*c1' + c1*c0' [0]: {:x}", d1.coeffs[0].to_words()[0]);

    // d2 = c1 * c1'
    let mut d2 = ct1.c1.clone();
    d2 *= &ct2.c1;
    println!("d2 = c1*c1' [0]: {:x}", d2.coeffs[0].to_words()[0]);

    // Step 2: Relinearization
    println!("Relin key a[0]: {:x}", relin_key.a.coeffs[0].to_words()[0]);
    println!("Relin key b[0]: {:x}", relin_key.b.coeffs[0].to_words()[0]);

    // c0_new = d0 + rk.b * d2
    let mut rk_b_times_d2 = relin_key.b.clone();
    rk_b_times_d2 *= &d2;
    println!("rk.b * d2 [0]: {:x}", rk_b_times_d2.coeffs[0].to_words()[0]);

    let mut c0_new = d0;
    c0_new += &rk_b_times_d2;
    println!(
        "c0_new = d0 + rk.b*d2 [0]: {:x}",
        c0_new.coeffs[0].to_words()[0]
    );

    // c1_new = d1 + rk.a * d2
    let mut rk_a_times_d2 = relin_key.a.clone();
    rk_a_times_d2 *= &d2;
    println!("rk.a * d2 [0]: {:x}", rk_a_times_d2.coeffs[0].to_words()[0]);

    let mut c1_new = d1;
    c1_new += &rk_a_times_d2;
    println!(
        "c1_new = d1 + rk.a*d2 [0]: {:x}",
        c1_new.coeffs[0].to_words()[0]
    );

    // Scale calculation
    println!("Scale bits before multiplication:");
    println!("  ct1.scale_bits: {}", ct1.scale_bits);
    println!("  ct2.scale_bits: {}", ct2.scale_bits);
    let new_scale_bits = ct1.scale_bits + ct2.scale_bits;
    println!(
        "  After multiplication: {} + {} = {}",
        ct1.scale_bits, ct2.scale_bits, new_scale_bits
    );
    println!("=== END MULTIPLICATION DEBUG ===");

    Ciphertext {
        c0: c0_new,
        c1: c1_new,
        scale_bits: new_scale_bits,
    }
}

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
    let values1 = vec![1.0, 0.0, 0.0, 0.0]; // Simple case
    let values2 = vec![2.0, 0.0, 0.0, 0.0]; // Should give [2.0, 0, 0, 0] right?
    println!("📊 Input data1: {values1:?}");
    println!("📊 Input data2: {values2:?}");

    // Step 1: Encode: Vec<f64> -> Plaintext
    println!("\n🔢 Encoding values...");
    let plaintext1 = encoder.encode(&values1, engine.context());
    let plaintext2 = encoder.encode(&values2, engine.context());
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
    let ciphertext1 = engine.encrypt(&plaintext1, &public_key, &mut rng);
    let ciphertext2 = engine.encrypt(&plaintext2, &public_key, &mut rng);
    println!("✅ Plaintext encrypted to ciphertext using U256 operations");

    println!("Multiplying civphertexts...");
    let mut ciphertext_mul =
        debug_mul_ciphertexts(&ciphertext1, &ciphertext2, &relin_key);

    println!("\nBefore rescale - first few coeffs:");
    println!("  c0[0]: {:x}", ciphertext_mul.c0.coeffs[0].to_words()[0]);
    println!("  c1[0]: {:x}", ciphertext_mul.c1.coeffs[0].to_words()[0]);

    // rescale using Kim's approach
    rescale_ciphertext(&mut ciphertext_mul, SCALE_BITS);

    println!("\nAfter rescale - first few coeffs:");
    println!("  c0[0]: {:x}", ciphertext_mul.c0.coeffs[0].to_words()[0]);
    println!("  c1[0]: {:x}", ciphertext_mul.c1.coeffs[0].to_words()[0]);

    // Step 3: Decrypt: Ciphertext → Plaintext
    println!("\n🔓 Decrypting ciphertext...");
    let decrypted_plaintext = Engine::decrypt(&ciphertext_mul, &secret_key);
    println!("✅ Ciphertext decrypted back to plaintext");

    println!("Debug scales:");
    println!(
        "  decrypted_plaintext.scale_bits: {}",
        decrypted_plaintext.scale_bits
    );

    let coeffs = decrypted_plaintext.poly.coefficients();
    println!("  First few coeffs: {:?}", &coeffs[0..4]);

    // Step 4: Decode: Plaintext → Vec<f64>
    println!("\n🔢 Decoding back to floating-point...");
    let decoded_values = encoder.decode(&decrypted_plaintext);
    println!("✅ Plaintext1 decoded: {:?}", decoded_values);

    Ok(())
}
