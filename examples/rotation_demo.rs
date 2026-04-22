//! Demonstrates homomorphic slot rotation in CKKS.
//!
//! ## Pipeline
//!
//! ```text
//! encode([1..16]) → encrypt
//!   → rotate(+1)          : [2,3,...,16,1]
//!   → add(original)       : [3,5,7,...,31,17]
//!   → rotate(+2)          : [7,9,...,31,17,3,5]
//! → decrypt → decode → verify
//! ```
//!
//! ## What this demonstrates
//!
//! Rotation requires a **rotation key** to switch the ciphertext from being
//! encrypted under `s(X^{5^k})` (the automorphism-permuted secret key) back
//! to `s(X)`.  The key-switching step uses the same RNS gadget decomposition
//! as relinearization after multiplication — but encodes `s_k` instead of `s²`.
//!
//! Unlike multiplication, rotation does **not** consume a level: `logp` and
//! `logq` are unchanged.  The example therefore uses only 2 primes (one spare
//! for headroom) and never calls `rescale_ciphertext`.
//!
//! ## Noise budget
//!
//! Gadget key switching adds noise ≈ `l × √N × q_i × σ` in polynomial space.
//! After decoding (÷ Δ), this is `l × √N × σ × (q_i / Δ)`.
//! To keep decoded noise small we need **Δ >> q_i** — achieved by using a large
//! encoding scale (SCALE=58) with small 30-bit NTT primes.  The decoded
//! key-switching noise is then ≈ `3 × √32 × 3.2 × 2^30 / 2^58 ≈ 2e-7`.
//!
//! ## Parameters
//!
//! | symbol | value | meaning                                          |
//! |--------|-------|--------------------------------------------------|
//! | N      |    32 | ring degree                                      |
//! | slots  |    16 | N/2 real plaintext slots                         |
//! | SCALE  |    58 | encoding scale bits (Δ = 2^58 >> q_i = 2^30)    |
//! | primes |     3 | three 30-bit NTT-friendly primes (Q ≈ 2^90 > Δ) |

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
/// Encoding scale bits.  Must satisfy Δ = 2^SCALE >> q_prime ≈ 2^30 so that
/// key-switching noise (≈ l × √N × q_prime × σ) is negligible after ÷ Δ.
const SCALE: u32 = 58;
/// Number of 30-bit NTT primes.  Q = q0×q1×q2 ≈ 2^90 > Δ = 2^58.
const NUM_PRIMES: usize = 3;
const ROT1: i32 = 1;
const ROT2: i32 = 2;

fn rotate_vec(v: &[f64], k: i32) -> Vec<f64> {
    let n = v.len();
    let k = k.rem_euclid(n as i32) as usize;
    v[k..].iter().chain(v[..k].iter()).copied().collect()
}

fn main() {
    println!("╔══════════════════════════════════════════════════════╗");
    println!("║        CKKS slot rotation demo  (N={N}, {SLOTS} slots)       ║");
    println!("╚══════════════════════════════════════════════════════╝\n");

    // ── 1. Setup ──────────────────────────────────────────────────────────────
    let primes = generate_primes(30, NUM_PRIMES, N as u64);
    println!("NTT-friendly primes (≡ 1 mod 2N = {}):", 2 * N);
    for (i, &p) in primes.iter().enumerate() {
        println!("  q{i} = {p}  ({}-bit)", 64 - p.leading_zeros());
    }

    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid primes"));
    let logq = basis.total_bits();
    println!("  total logq = {logq}\n");

    let engine = CkksEngine::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        },
    );
    let encoder = CkksEncoder::<N>::new(SCALE);
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    // ── 2. Key generation ─────────────────────────────────────────────────────
    println!("Generating keys …");
    let sk = engine.generate_secret_key(&mut rng).expect("sk");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");
    let rotk1 = engine.generate_gadget_rotation_key(&sk, ROT1, &mut rng);
    let rotk2 = engine.generate_gadget_rotation_key(&sk, ROT2, &mut rng);
    println!("  sk, pk, rotk({ROT1:+}), rotk({ROT2:+}) — done.\n");

    // ── 3. Encode + encrypt ───────────────────────────────────────────────────
    let input: Vec<f64> = (1..=SLOTS).map(|i| i as f64).collect();
    println!("Input slots: {:?}\n", &input);

    let pt_in = encoder.encode(&input, basis.clone());
    let ct = engine.encrypt(&pt_in, &pk, logq, &mut rng);

    // ── 4. Homomorphic pipeline ───────────────────────────────────────────────

    // Step A: rotate by ROT1
    println!("Step A: rotate ciphertext by {ROT1:+} slot(s) …");
    let ct_rot1 = CkksEngine::<RnsPoly<N>, N>::rotate_ciphertext(&ct, &rotk1);
    println!(
        "  logp={}, logq={} (unchanged — rotation is level-free)\n",
        ct_rot1.logp, ct_rot1.logq
    );

    // Step B: add original and rotated ciphertexts
    println!("Step B: add original + rotated …");
    let ct_added = CkksEngine::<RnsPoly<N>, N>::add_ciphertexts(&ct, &ct_rot1);
    println!("  logp={}, logq={}\n", ct_added.logp, ct_added.logq);

    // Step C: rotate sum by ROT2
    println!("Step C: rotate sum by {ROT2:+} slot(s) …");
    let ct_result =
        CkksEngine::<RnsPoly<N>, N>::rotate_ciphertext(&ct_added, &rotk2);
    println!("  logp={}, logq={}\n", ct_result.logp, ct_result.logq);

    // ── 5. Decrypt + decode ───────────────────────────────────────────────────
    println!("Decrypting …");
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_result, &sk);
    let decoded = encoder.decode(&pt_out);

    // ── 6. Plaintext reference ────────────────────────────────────────────────
    let ref_rot1 = rotate_vec(&input, ROT1);
    let ref_added: Vec<f64> =
        input.iter().zip(&ref_rot1).map(|(a, b)| a + b).collect();
    let expected = rotate_vec(&ref_added, ROT2);

    // ── 7. Results table ──────────────────────────────────────────────────────
    println!();
    let mut t = table::new([
        "slot",
        "input",
        &format!("rot({ROT1:+})"),
        "add",
        &format!("rot({ROT2:+})"),
        "decoded",
        "expected",
    ]);
    for i in 0..SLOTS {
        t.add_row([
            i.to_string(),
            format!("{:.1}", input[i]),
            format!("{:.1}", ref_rot1[i]),
            format!("{:.1}", ref_added[i]),
            format!("{:.3}", expected[i]),
            format!("{:.6}", decoded[i]),
            format!("{:.3}", expected[i]),
        ]);
    }
    println!("{t}");

    // ── 8. Error check ────────────────────────────────────────────────────────
    let max_err = expected
        .iter()
        .zip(decoded.iter().take(SLOTS))
        .map(|(e, g)| (e - g).abs())
        .fold(0f64, f64::max);

    println!();
    println!("Max absolute error over {SLOTS} slots: {max_err:.3e}");

    let bound = 1e-4_f64;
    if max_err <= bound {
        println!("✓  All slots within error bound ({bound:.0e}).");
    } else {
        println!("✗  ERROR: {max_err:.3e} exceeds bound {bound:.0e}");
        std::process::exit(1);
    }
}
