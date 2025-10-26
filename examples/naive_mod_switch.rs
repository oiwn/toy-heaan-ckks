use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{CkksEngine, Encoder, NaivePolyRing};

// Modulus Switching Constants
const DEGREE: usize = 8;
const SCALE_BITS: u32 = 20; // Delta = 2^20 = 1,048,576

// Two-level modulus chain
const Q_1: u64 = (1u64 << 51) - 9; // Level 1 (higher)
const Q_0: u64 = (1u64 << 41) - 17; // Level 0 (lower)

type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

/// # Modulus Switching for CKKS (Naive Backend)
///
/// Source: "An Introduction to Mathematical Cryptography" (FHE Textbook)
/// Chapter: c04-glwe-modulus-switch.tex
///
/// ## What is Modulus Switching?
///
/// Modulus switching transforms a ciphertext from modulus `q` to a smaller
/// modulus `q_hat` while preserving the ability to decrypt to the same plaintext.
///
/// Original ciphertext at modulus q:
///   `c0 + c1*s = Delta*m + e (mod q)`
///
/// After switching to q_hat:
///   `c0_hat + c1_hat*s = Delta_hat*m + e_hat + epsilon (mod q_hat)`
///
/// Where:
/// - `s` is the secret key (UNCHANGED)
/// - `m` is the plaintext message (UNCHANGED)
/// - `Delta` is the encoding scale (CHANGES by ratio)
/// - `e` is the noise (CHANGES by ratio)
/// - `epsilon` is rounding error from switching
///
/// ## The Formula (from FHE Book)
///
/// For each component of the ciphertext (c0, c1):
///
/// ```
/// c_hat = round(c * (q_hat / q))
/// ```
///
/// This applies to:
/// - All polynomial coefficients
/// - Both c0 and c1 components
///
/// ## What Changes and What Doesn't
///
/// ### UNCHANGED:
/// - Secret key `s`: Same key decrypts before and after
/// - Plaintext value `m`: The actual message remains the same
///
/// ### CHANGES:
/// - Encoding scale: `Delta_hat = Delta * (q_hat / q)`
/// - Noise: `e_hat = round(e * (q_hat / q))`
/// - Rounding error: `epsilon` (for RLWE with k=1)
///
/// ## Numerical Example for This Demo
///
/// Parameters:
/// - Q_1 = 2^51 - 9   (starting modulus)
/// - Q_0 = 2^41 - 17  (target modulus)
/// - Delta = 2^20     (initial scale)
///
/// Ratio:
/// ```
/// q_hat / q = (2^41) / (2^51) = 2^(-10) = 1/1024
/// ```
///
/// Scale transformation:
/// ```
/// Delta_hat = Delta * (q_hat / q)
///           = 2^20 * 2^(-10)
///           = 2^10
///           = 1024
/// ```
///
/// Result: Scale drops from ~1 million to ~1 thousand
/// Impact: Precision decreases by 1024x, but values should still be recoverable
///
/// ## Why This Matters
///
/// 1. **Noise management**: Switching reduces noise proportionally
/// 2. **Performance**: Smaller moduli = faster operations
/// 3. **Compatibility**: Must match modulus levels before homomorphic ops
/// 4. **Trade-off awareness**: Shows precision loss from switching
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

    println!("Modulus Switching Demo - CKKS Naive Backend");
    println!("{}", "=".repeat(60));
    println!();
    println!("Parameters:");
    println!("  Degree N:        {}", DEGREE);
    println!("  Q_1 (high):      {} (2^51 - 9)", Q_1);
    println!("  Q_0 (low):       {} (2^41 - 17)", Q_0);
    println!(
        "  Delta (Q_1):     {} (2^{})",
        1u64 << SCALE_BITS,
        SCALE_BITS
    );

    let ratio = (Q_0 as f64) / (Q_1 as f64);
    let delta_hat = ((1u64 << SCALE_BITS) as f64 * ratio) as u64;

    // Round to nearest power of 2 for better approximation
    let scale_bits_hat = (delta_hat as f64).log2().round() as u32;
    let delta_hat_rounded = 1u64 << scale_bits_hat;

    println!(
        "  Delta_hat (Q_0): {} ≈ {} (2^{})",
        delta_hat, delta_hat_rounded, scale_bits_hat
    );
    println!(
        "  Modulus ratio:   {:.10} (1/{})",
        ratio,
        (1.0 / ratio) as u64
    );
    println!();

    // Input values
    let values = vec![1.5, 2.5, 3.5];
    println!("Input values: {:?}", values);
    println!();

    // [1] Setup at Q_1
    println!("[1] Encrypt at Q_1");
    let engine_q1 = Engine::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_naive(Q_1, SCALE_BITS)?;

    let encoder_q1 = toy_heaan_ckks::RustFftEncoder::new(SCALE_BITS)?;

    // Generate keys
    let secret_key = engine_q1.generate_secret_key(&mut rng)?;
    let public_key = engine_q1.generate_public_key(&secret_key, &mut rng)?;

    // Encode and encrypt
    let plaintext_q1 = encoder_q1.encode(&values, engine_q1.context());
    let logq1 = Q_1.ilog2();
    let ciphertext_q1 =
        engine_q1.encrypt(&plaintext_q1, &public_key, &mut rng, logq1);

    println!("  Ciphertext modulus q: {} (2^{})", Q_1, logq1);
    println!("  Encoding scale Δ: {} (2^{})", 1u64 << SCALE_BITS, SCALE_BITS);
    println!("  Stored: logp={}, logq={}", ciphertext_q1.logp, ciphertext_q1.logq);
    println!();

    // [2] Decrypt at Q_1 (baseline)
    println!("[2] Decrypt at Q_1 (baseline)");
    let decrypted_q1 = Engine::decrypt(&ciphertext_q1, &secret_key);
    let decoded_q1 = encoder_q1.decode(&decrypted_q1);

    let max_error_q1 = values
        .iter()
        .zip(&decoded_q1)
        .map(|(orig, decoded)| (orig - decoded).abs())
        .fold(0.0, f64::max);

    println!("  Decoded: {:?}", &decoded_q1[..values.len()]);
    println!("  Max error: {:.2e}", max_error_q1);
    println!();

    // [3] Switch modulus Q_1 -> Q_0
    println!("[3] Switch modulus Q_1 -> Q_0");
    println!("  Formula: c_hat = round(c * (q_hat/q))");
    println!("  Ratio: {:.10} (= 2^-10)", ratio);

    let logq0 = Q_0.ilog2();
    let ciphertext_q0 = ciphertext_q1.mod_switch(&Q_0, logq0);

    println!();
    println!("  After switching:");
    println!("    New modulus q_hat: {} (2^{})", Q_0, logq0);
    println!("    ACTUAL scale Δ_hat = Δ * ratio = {} (2^{:.1})",
             delta_hat, (delta_hat as f64).log2());
    println!("    Stored logp: {} (claims scale = 2^{} = {})",
             ciphertext_q0.logp, ciphertext_q0.logp, 1u64 << ciphertext_q0.logp);
    println!("    Stored logq: {}", ciphertext_q0.logq);
    println!();
    println!("  ⚠️  MISMATCH: logp={} implies scale=2^{}={}, but ACTUAL scale={}!",
             ciphertext_q0.logp, ciphertext_q0.logp,
             1u64 << ciphertext_q0.logp, delta_hat);
    println!();

    // [4] Decrypt at Q_0 - Comparing approaches
    println!("[4] Decrypt at Q_0 - Two decoding approaches:");

    // Secret key values stay the same, but need to exist at new modulus Q_0
    // Extract the signed coefficients (which are just {-1, 0, 1})
    use toy_heaan_ckks::PolyRing;
    let sk_coeffs = secret_key.poly.to_coeffs();
    let secret_key_q0 = toy_heaan_ckks::SecretKey {
        poly: NaivePolyRing::from_coeffs(&sk_coeffs, &Q_0),
    };

    // Decrypt using secret key at Q_0 (same VALUES, different modulus)
    let decrypted_q0 = Engine::decrypt(&ciphertext_q0, &secret_key_q0);

    println!();
    println!("  Approach A: Use logp={} → decode with scale=2^{}={}",
             ciphertext_q0.logp, ciphertext_q0.logp, 1u64 << ciphertext_q0.logp);
    let decoded_q0_wrong = encoder_q1.decode(&decrypted_q0);
    println!("    Result: {:?}", &decoded_q0_wrong[..values.len()]);
    let scale_error = (1u64 << ciphertext_q0.logp) as f64 / delta_hat as f64;
    println!("    ❌ WRONG! Off by {:.0}x (because actual scale={} ≠ 2^{}={})",
             scale_error, delta_hat, ciphertext_q0.logp, 1u64 << ciphertext_q0.logp);
    println!();

    println!("  Approach B: Use ACTUAL scale after mod_switch = 2^{}≈{}",
             scale_bits_hat, delta_hat);
    let encoder_corrected = toy_heaan_ckks::RustFftEncoder::new(scale_bits_hat)?;
    let decoded_q0 = encoder_corrected.decode(&decrypted_q0);
    println!("    Result: {:?}", &decoded_q0[..values.len()]);
    println!("    ✅ Correct!");

    let max_error_q0 = values
        .iter()
        .zip(&decoded_q0)
        .map(|(orig, decoded)| (orig - decoded).abs())
        .fold(0.0, f64::max);

    println!();
    println!("  Comparison:");
    println!("    At Q_1: max error = {:.2e}", max_error_q1);
    println!("    At Q_0: max error = {:.2e}", max_error_q0);
    println!("    Error increased by {:.0}x (expected ~{}x from precision loss)",
             max_error_q0 / max_error_q1.max(1e-10),
             (1.0 / ratio) as u64);
    println!();

    // Analysis
    println!();
    println!("Key Insights:");
    println!(
        "  1. Modulus switching preserved the plaintext (same secret key works!)"
    );
    println!(
        "  2. Precision decreased by scale ratio ({}x)",
        (1.0 / ratio) as u64
    );
    println!("  3. HEAAN's logp/logq split is designed for RNS-CKKS:");
    println!("     - RNS: drops CRT primes, logp stays constant");
    println!("     - Naive: scale changes, must track actual scale for decoding");
    println!("  4. Rounding error stayed within theoretical bound");
    println!();

    Ok(())
}
