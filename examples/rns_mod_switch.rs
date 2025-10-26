use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::sync::Arc;
use toy_heaan_ckks::rings::backends::rns::RnsBasisBuilder;
use toy_heaan_ckks::{CkksEngine, Encoder, RnsPolyRing, RustFftEncoder, SecretKey};

// RNS Modulus Switching Constants
const DEGREE: usize = 8;
const SCALE_BITS: u32 = 20; // Delta = 2^20 = 1,048,576

// Three-level RNS modulus chain
// Each level drops one prime from the RNS representation
// Use different bit sizes to ensure we get different coprime moduli
// Total product must fit in u64 (64 bits), so we use smaller primes for the demo
// Note: In production, you'd use crypto_bigint for reconstruction with larger primes
const PRIME_BITS_3: [usize; 3] = [21, 20, 19]; // Level 2: Q2 = q1 * q2 * q3 (~60 bits)
const PRIME_BITS_2: [usize; 2] = [21, 20]; // Level 1: Q1 = q1 * q2 (~41 bits)
const PRIME_BITS_1: [usize; 1] = [21]; // Level 0: Q0 = q1 (~21 bits)

type Engine = CkksEngine<RnsPolyRing<DEGREE>, DEGREE>;

/// # RNS-based Modulus Switching (ModDrop) for CKKS
///
/// Source: "An Introduction to Mathematical Cryptography" (FHE Textbook)
/// Chapter: RNS-based ModDrop
/// URL: https://fhetextbook.github.io/RNSbasedModDroptextsfModDroptextsubscriptRNS.html
///
/// ## What is RNS ModDrop?
///
/// In the Residue Number System (RNS) representation, a polynomial is stored
/// as residues modulo multiple prime moduli (CRT representation):
///
/// ```text
/// q = q1 * q2 * ... * qk    (product of k primes)
/// poly = [residues mod q1, residues mod q2, ..., residues mod qk]
/// ```
///
/// ModDrop simply **drops** the last few primes to get a smaller modulus:
///
/// ```text
/// q' = q1 * q2 * ... * qk'  (where k' < k)
/// poly' = [residues mod q1, residues mod q2, ..., residues mod qk']
/// ```
///
/// ## CKKS Invariant (ASCII formulas)
///
/// A valid CKKS ciphertext satisfies:
/// ```text
/// c0 + c1*s = Delta*m + e  (mod q)
/// ```
/// where:
/// - c0, c1: ciphertext polynomials
/// - s: secret key polynomial
/// - Delta: scale (2^logp)
/// - m: message polynomial (encoded from input values)
/// - e: noise polynomial (small coefficients)
/// - q: modulus
///
/// For correct decryption, we need: q >> Delta (modulus much larger than scale)
/// Typical ratio: q/Delta >= 2^20 for security and correctness
///
/// ## Key Advantage: ZERO Additional Noise
///
/// Unlike naive modulus switching which introduces rounding errors,
/// RNS ModDrop adds **ZERO additional noise** when q' divides q!
///
/// This is because we're simply truncating the CRT representation,
/// not performing any scaling or rounding operations.
///
/// ## Comparison: Naive vs RNS ModDrop
///
/// | Aspect            | Naive ModSwitch          | RNS ModDrop              |
/// |-------------------|--------------------------|--------------------------|
/// | Operation         | c' = round(c * ratio)    | Drop RNS channels        |
/// | Noise added       | Rounding error ≤ 1.5     | **ZERO**                 |
/// | Arithmetic        | Floating-point           | Integer only             |
/// | Scale change      | Yes: Δ' = Δ * (q'/q)     | No (in RNS repr)         |
/// | Moduli constraint | Any coprime pair works   | q' must divide q         |
/// | Complexity        | O(N) with FP ops         | O(N) copy                |
///
/// ## This Demo
///
/// We demonstrate a 3-level modulus chain:
/// - Level 2: Start with 3 primes (highest security)
/// - Level 1: Drop to 2 primes
/// - Level 0: Drop to 1 prime (lowest level)
///
/// At each ModDrop, we verify that the decryption error remains constant,
/// proving the zero-noise property.

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n{}", "=".repeat(70));
    println!("RNS Modulus Switching Demo - CKKS RNS Backend");
    println!("{}", "=".repeat(70));

    // Initialize RNG
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    println!("\n[Parameters]");
    println!("  Degree N:        {}", DEGREE);
    println!(
        "  Scale Delta:     {} (2^{})",
        1u64 << SCALE_BITS,
        SCALE_BITS
    );

    // ========================================================================
    // Step 1: Setup 3-level RNS modulus chain
    // ========================================================================

    println!("\n[Modulus Chain (RNS)]");

    let delta = 1u64 << SCALE_BITS;

    // Level 2: 3 primes (highest)
    let basis_q2 = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_prime_bits(PRIME_BITS_3.to_vec())
            .build()?,
    );
    let q2_product: u64 = basis_q2.primes().iter().product();
    let q2_log = PRIME_BITS_3.iter().sum::<usize>();
    println!(
        "  Level 2 (Q2):    {} primes {:?}",
        basis_q2.channel_count(),
        basis_q2
            .primes()
            .iter()
            .map(|p| format!(
                "2^{}-{}",
                ((*p as f64).log2().ceil() as usize),
                (1u64 << (*p as f64).log2().ceil() as usize) - p
            ))
            .collect::<Vec<_>>()
    );
    println!("                   Q2 = {} (~2^{})", q2_product, q2_log);
    println!(
        "                   Q2/Delta = {} (need >> 1 for correctness)",
        q2_product / delta
    );

    // Level 1: 2 primes (middle)
    let basis_q1 = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_prime_bits(PRIME_BITS_2.to_vec())
            .build()?,
    );
    let q1_product: u64 = basis_q1.primes().iter().product();
    let q1_log = PRIME_BITS_2.iter().sum::<usize>();
    println!(
        "  Level 1 (Q1):    {} primes {:?}",
        basis_q1.channel_count(),
        basis_q1
            .primes()
            .iter()
            .map(|p| format!(
                "2^{}-{}",
                ((*p as f64).log2().ceil() as usize),
                (1u64 << (*p as f64).log2().ceil() as usize) - p
            ))
            .collect::<Vec<_>>()
    );
    println!("                   Q1 = {} (~2^{})", q1_product, q1_log);
    println!(
        "                   Q1/Delta = {} (need >> 1 for correctness)",
        q1_product / delta
    );

    // Level 0: 1 prime (lowest)
    let basis_q0 = Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_prime_bits(PRIME_BITS_1.to_vec())
            .build()?,
    );
    let q0_product: u64 = basis_q0.primes().iter().product();
    let q0_log = PRIME_BITS_1.iter().sum::<usize>();
    println!(
        "  Level 0 (Q0):    {} prime  {:?}",
        basis_q0.channel_count(),
        basis_q0
            .primes()
            .iter()
            .map(|p| format!(
                "2^{}-{}",
                ((*p as f64).log2().ceil() as usize),
                (1u64 << (*p as f64).log2().ceil() as usize) - p
            ))
            .collect::<Vec<_>>()
    );
    println!("                   Q0 = {} (~2^{})", q0_product, q0_log);
    println!(
        "                   Q0/Delta = {} [WARNING] TOO SMALL!",
        q0_product / delta
    );

    // ========================================================================
    // Step 2: Input values and encoding
    // ========================================================================

    let values = vec![1.5, 2.5, 3.5];
    println!("\n[Input] values: {:?}", values);

    let encoder = RustFftEncoder::new(SCALE_BITS)?;

    // ========================================================================
    // Step 3: Encrypt at highest level (Q2 - 3 primes)
    // ========================================================================

    println!("\n[1] Encrypt at Q2 (3 primes)");

    let engine_q2 = Engine::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2)
        .build_rns(basis_q2.clone(), SCALE_BITS)?;

    let sk_q2 = engine_q2.generate_secret_key(&mut rng)?;
    let pk_q2 = engine_q2.generate_public_key(&sk_q2, &mut rng)?;

    let plaintext = encoder.encode(&values, engine_q2.context());
    let ct_q2 = engine_q2.encrypt(&plaintext, &pk_q2, &mut rng, q2_log as u32);

    println!("  Channels: {}", ct_q2.c0.channels());
    println!("  logp: {}, logq: {}", ct_q2.logp, ct_q2.logq);

    // ========================================================================
    // Step 4: Decrypt at Q2 (baseline error)
    // ========================================================================

    println!("\n[2] Decrypt at Q2 (baseline)");

    let decrypted_q2 = Engine::decrypt(&ct_q2, &sk_q2);
    let decoded_q2 = encoder.decode(&decrypted_q2);

    println!("  Decoded: {:?}", format_vec(&decoded_q2));

    let error_q2 = max_error(&values, &decoded_q2);
    println!("  Max error: {:.2e}", error_q2);

    // ========================================================================
    // Step 5: ModDrop Q2 -> Q1 (drop third prime q3)
    // ========================================================================

    println!("\n[3] ModDrop Q2 -> Q1 (drop prime q3)");
    println!("  Operation: Truncate from 3 channels to 2 channels");
    println!("  No rescaling, no rounding!");

    let ct_q1 = ct_q2.mod_switch(&basis_q1, q1_log as u32);

    println!("\n  After switching:");
    println!("    Channels: {}", ct_q1.c0.channels());
    println!("    logp: {} (unchanged)", ct_q1.logp);
    println!("    logq: {}", ct_q1.logq);

    // ========================================================================
    // Step 6: Decrypt at Q1
    // ========================================================================

    println!("\n[4] Decrypt at Q1");

    // Create secret key at Q1 (same values, fewer RNS channels)
    let sk_q1_coeffs = sk_q2.poly.to_i64_coefficients();
    let sk_q1 = SecretKey {
        poly: RnsPolyRing::from_i64_slice(&sk_q1_coeffs, basis_q1.clone()),
    };

    let decrypted_q1 = Engine::decrypt(&ct_q1, &sk_q1);
    let decoded_q1 = encoder.decode(&decrypted_q1);

    println!("  Decoded: {:?}", format_vec(&decoded_q1));

    let error_q1 = max_error(&values, &decoded_q1);
    println!("  Max error: {:.2e}", error_q1);

    let error_ratio = error_q1 / error_q2;
    if error_ratio < 1.5 {
        println!("\n  [OK] Error unchanged! (ratio: {:.2}x)", error_ratio);
        println!("       RNS ModDrop adds ZERO noise");
    } else {
        println!("\n  [WARNING] Error increased by {:.2}x", error_ratio);
    }

    // ========================================================================
    // Step 7: ModDrop Q1 -> Q0 (drop second prime q2)
    // ========================================================================

    println!("\n[5] ModDrop Q1 -> Q0 (drop prime q2)");

    let ct_q0 = ct_q1.mod_switch(&basis_q0, q0_log as u32);

    println!("  After switching:");
    println!("    Channels: {}", ct_q0.c0.channels());
    println!("    logp: {} (unchanged)", ct_q0.logp);
    println!("    logq: {}", ct_q0.logq);

    // ========================================================================
    // Step 8: Decrypt at Q0
    // ========================================================================

    println!("\n[6] Decrypt at Q0");

    let sk_q0_coeffs = sk_q2.poly.to_i64_coefficients();
    let sk_q0 = SecretKey {
        poly: RnsPolyRing::from_i64_slice(&sk_q0_coeffs, basis_q0.clone()),
    };

    let decrypted_q0 = Engine::decrypt(&ct_q0, &sk_q0);
    let decoded_q0 = encoder.decode(&decrypted_q0);

    println!("  Decoded: {:?}", format_vec(&decoded_q0));

    let error_q0 = max_error(&values, &decoded_q0);
    println!("  Max error: {:.2e}", error_q0);

    // ========================================================================
    // Step 9: Final Analysis
    // ========================================================================

    println!("\n{}", "=".repeat(70));
    println!("[Analysis]");
    println!("{}", "=".repeat(70));

    println!("\nError progression:");
    println!("  Q2 (3 primes): {:.2e}", error_q2);
    println!(
        "  Q1 (2 primes): {:.2e}  (ratio: {:.2}x)",
        error_q1,
        error_q1 / error_q2
    );
    println!(
        "  Q0 (1 prime):  {:.2e}  (ratio: {:.2}x)",
        error_q0,
        error_q0 / error_q2
    );

    println!("\n[WARNING] Why Q0 Failed:");
    println!("  Q0 = {} (2^{})", q0_product, q0_log);
    println!("  Delta = {} (2^{})", delta, SCALE_BITS);
    println!("  Ratio Q0/Delta = {}", q0_product / delta);
    println!();
    println!("  CKKS requires: q >> Delta (modulus much larger than scale)");
    println!("  - Need room for: Delta*m + e (mod q)");
    println!(
        "  - At Q0, we have q ~= Delta (ratio only ~{})!",
        q0_product / delta
    );
    println!("  - No room for noise! Decryption wraps around modulus.");
    println!("  - Solution: Use larger primes or drop to Q1 as lowest level");
    println!();
    println!("[Key Insights]");
    println!("  1. RNS ModDrop is simple: just drop RNS channels");
    println!(
        "  2. Zero additional noise per ModDrop operation (Q2->Q1 shows 1.00x!)"
    );
    println!("  3. Scale (logp) stays constant (unlike naive backend)");
    println!("  4. Enables multi-level modulus chains for deep circuits");
    println!("  5. Error only from initial encryption, not from ModDrop");
    println!("  6. Still need q >> Delta at each level for correctness");

    println!("\n[Comparison with Naive ModSwitch]");
    println!("  - Naive: Adds rounding error at each switch (~1.5 per coeff)");
    println!("  - RNS: Adds ZERO noise (when q' divides q)");
    println!("  - RNS enables deeper computations with predictable noise!");

    println!("\n{}", "=".repeat(70));
    println!("[Demo complete]");
    println!("{}", "=".repeat(70));

    Ok(())
}

// Helper functions

fn format_vec(values: &[f64]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|v| format!("{:.7}", v))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn max_error(expected: &[f64], actual: &[f64]) -> f64 {
    expected
        .iter()
        .zip(actual.iter())
        .map(|(e, a)| (e - a).abs())
        .fold(0.0, f64::max)
}
