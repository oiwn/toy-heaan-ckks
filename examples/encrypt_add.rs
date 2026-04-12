//! CKKS homomorphic addition demonstration.
//!
//! Full pipeline:
//!
//!   encode → encrypt → homomorphic add → decrypt → decode
//!
//! Two vectors of real values are encoded and encrypted independently.
//! Their ciphertexts are added homomorphically (no decryption needed).
//! The result is decrypted and decoded, then compared to the expected
//! element-wise sums.

use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

/// Ring degree — polynomial lives in Z[X]/(X^N + 1), N/2 real slots available.
const N: usize = 16;

/// Encoding scale: Δ = 2^SCALE_BITS. Controls precision vs coefficient size.
const SCALE_BITS: u32 = 30;

fn main() {
    println!("╔═══════════════════════════════════════╗");
    println!("║   CKKS Homomorphic Addition  —  Demo  ║");
    println!("╚═══════════════════════════════════════╝\n");

    // ── 1. Parameters ──────────────────────────────────────────────────────────
    let delta = 2f64.powi(SCALE_BITS as i32);
    println!("Ring degree N    : {N}");
    println!("Slots  (N/2)     : {}", N / 2);
    println!("Scale bits       : {SCALE_BITS}  →  Δ ≈ {delta:.3e}\n");

    // ── 2. RNS basis ───────────────────────────────────────────────────────────
    // Three 31-bit NTT-friendly primes → Q ≈ 2^93, large enough to hold
    // scaled coefficients without wrap-around.
    let primes = generate_primes(31, 3, N as u64);
    println!("RNS primes (each ≡ 1 mod 2N = {}):", 2 * N);
    for (i, &p) in primes.iter().enumerate() {
        println!("  q{i} = {p}");
    }

    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();
    println!("  logq (Σ floor(log2 qᵢ)) = {logq}\n");

    // ── 3. Engine & keys ───────────────────────────────────────────────────────
    let params = CkksParams {
        error_variance: 3.2,
        hamming_weight: N / 2,
        scale_bits: SCALE_BITS,
    };
    let engine = CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params);

    let mut rng = ChaCha20Rng::seed_from_u64(42);
    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    println!("Keys generated (secret + public).\n");

    // ── 4. Input vectors ───────────────────────────────────────────────────────
    let a: [f64; 4] = [1.0, 2.0, 3.0, 4.0];
    let b: [f64; 4] = [0.5, -1.5, 0.25, -0.75];
    let expected: Vec<f64> = a.iter().zip(b.iter()).map(|(x, y)| x + y).collect();

    println!("Input A : {:?}", a);
    println!("Input B : {:?}", b);
    println!("Expected: {:?}\n", expected);

    // ── 5. Encode ──────────────────────────────────────────────────────────────
    let encoder = CkksEncoder::<N>::new(SCALE_BITS);
    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());

    // ── 6. Encrypt ─────────────────────────────────────────────────────────────
    let ct_a = engine.encrypt(&pt_a, &pk, logq, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, logq, &mut rng);
    println!(
        "Both plaintexts encrypted (logp={}, logq={}).\n",
        ct_a.logp, ct_a.logq
    );

    // ── 7. Homomorphic add ─────────────────────────────────────────────────────
    let ct_sum = CkksEngine::<RnsPoly<N>, N>::add_ciphertexts(&ct_a, &ct_b);
    println!("Homomorphic addition done (no secret key used).");
    println!("  Result logp={}, logq={}\n", ct_sum.logp, ct_sum.logq);

    // ── 8. Decrypt ─────────────────────────────────────────────────────────────
    let pt_sum = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_sum, &sk);

    // ── 9. Decode ──────────────────────────────────────────────────────────────
    let decoded = encoder.decode(&pt_sum);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    // ── 10. Verify ─────────────────────────────────────────────────────────────
    println!(
        "─── Results ──────────────────────────────────────────────────────\n"
    );
    println!(
        "{:<6} {:>10} {:>10} {:>12} {:>10}",
        "slot", "a[i]", "b[i]", "expected", "decoded"
    );
    let mut max_err = 0f64;
    for (i, (&exp, &got)) in expected.iter().zip(result.iter()).enumerate() {
        let err = (exp - got).abs();
        max_err = max_err.max(err);
        println!(
            "{:<6} {:>10.4} {:>10.4} {:>12.6} {:>10.6}  |err|={err:.2e}",
            i, a[i], b[i], exp, got
        );
    }

    println!("\nMax absolute error : {max_err:.2e}");
    // Noise budget: Gaussian errors (σ=3.2) in e0, e1 are folded into slots via
    // the canonical embedding. Each slot error ≈ σ·√(hamming_weight·N)/Δ (from
    // the e1·s term, which dominates). Two ciphertexts add two independent
    // errors, so the expected max grows by ≈ √2. We use a 10× safety factor
    // for the small slot count.
    let sigma = 3.2f64;
    let hw = (N / 2) as f64;
    let enc_noise = 10.0 * sigma * (hw * N as f64).sqrt() / delta;
    let bound = enc_noise + 4.0 / delta; // + encoding rounding
    if max_err <= bound {
        println!("✓  within expected bound (≲ {bound:.2e})");
    } else {
        println!("✗  ERROR: {max_err:.2e} exceeds bound {bound:.2e}");
        std::process::exit(1);
    }

    println!("\nDone.");
}
