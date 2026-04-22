//! Stress-tests how many sequential rotations can be applied before noise
//! degrades decoding accuracy.
//!
//! Each rotation performs one gadget key-switch, adding noise ≈ l×√N×σ×(q/Δ).
//! With SCALE=58 and 30-bit primes (Δ >> q_prime) this is ≈ 2.5e-7 per step.
//! Noise accumulates linearly: `total ≈ k × 2.5e-7`.
//!
//! The example rotates by +1 repeatedly, checking the error at each checkpoint.

use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
use toy_heaan_ckks::table;

const N: usize = 32;
const SLOTS: usize = N / 2; // 16
const SCALE: u32 = 58;

fn main() {
    println!("╔═══════════════════════════════════════════════════════╗");
    println!("║  CKKS rotation stress test  (N={N}, {SLOTS} slots, Δ=2^{SCALE})  ║");
    println!("╚═══════════════════════════════════════════════════════╝\n");

    // ── Setup ─────────────────────────────────────────────────────────────────
    let primes = generate_primes(30, 3, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid primes"));
    let logq = basis.total_bits();

    let engine = CkksEngine::new(
        basis.clone(),
        CkksParams { error_variance: 3.2, hamming_weight: N / 2, scale_bits: SCALE },
    );
    let encoder = CkksEncoder::<N>::new(SCALE);
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    let sk = engine.generate_secret_key(&mut rng).unwrap();
    let pk = engine.generate_public_key(&sk, &mut rng).unwrap();

    // One rotation key for offset +1 — reused at every step.
    let rotk = engine.generate_gadget_rotation_key(&sk, 1, &mut rng);

    // ── Encrypt ───────────────────────────────────────────────────────────────
    let input: Vec<f64> = (1..=SLOTS).map(|i| i as f64).collect();
    let pt = encoder.encode(&input, basis.clone());
    let mut ct = engine.encrypt(&pt, &pk, logq, &mut rng);

    // Plaintext reference: cyclic left-shift by 1 each step.
    let mut expected = input.clone();

    // ── Stress loop ───────────────────────────────────────────────────────────
    let checkpoints = [1usize, 5, 10, 50, 100, 200, 400, 800];
    let mut prev = 0usize;

    let mut t = table::new(["k (rots)", "max_err", "slot[0] got", "slot[0] want", "ok"]);

    for &checkpoint in &checkpoints {
        let steps = checkpoint - prev;
        for _ in 0..steps {
            ct = CkksEngine::<RnsPoly<N>, N>::rotate_ciphertext(&ct, &rotk);
            expected = expected[1..].iter().chain(expected[..1].iter()).copied().collect();
        }
        prev = checkpoint;

        let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct, &sk);
        let decoded = encoder.decode(&pt_out);

        let max_err = expected
            .iter()
            .zip(&decoded)
            .map(|(e, g)| (e - g).abs())
            .fold(0f64, f64::max);

        t.add_row([
            checkpoint.to_string(),
            format!("{max_err:.3e}"),
            format!("{:.6}", decoded[0]),
            format!("{:.1}", expected[0]),
            if max_err < 1e-3 { "✓" } else { "✗" }.to_string(),
        ]);
    }

    println!("{t}");

    println!("\nNoise model: ≈ k × l×√N×σ×(q_prime/Δ)");
    println!("           = k × 3 × √32 × 3.2 × 2^30 / 2^58");
    println!("           ≈ k × 2.0e-7");
    println!("\nNote: multiplications are the real depth bottleneck — each");
    println!("mul+rescale consumes one prime level. Rotations are level-free.");
}
