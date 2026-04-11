//! CKKS homomorphic multiplication demonstration.
//!
//! Full pipeline:
//!
//!   encode → encrypt → homomorphic mul (gadget relin) → rescale → decrypt → decode
//!
//! Two vectors of real values are encoded and encrypted independently.
//! Their ciphertexts are multiplied homomorphically (no decryption needed).
//! After rescaling the result is decrypted and decoded, then compared to the
//! element-wise products.

use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::engine::{CkksEngine, CkksParams, RnsGadgetRelinKey};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::keys::SecretKey;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

/// Ring degree — polynomial lives in Z[X]/(X^N + 1), N/2 real slots available.
const N: usize = 16;

/// Encoding scale: Δ = 2^SCALE_BITS.
const SCALE_BITS: u32 = 30;

fn main() {
    println!("╔═══════════════════════════════════════════╗");
    println!("║  CKKS Homomorphic Multiplication  —  Demo ║");
    println!("╚═══════════════════════════════════════════╝\n");

    // ── 1. Parameters ──────────────────────────────────────────────────────────
    let delta = 2f64.powi(SCALE_BITS as i32);
    println!("Ring degree N    : {N}");
    println!("Slots  (N/2)     : {}", N / 2);
    println!("Scale bits       : {SCALE_BITS}  →  Δ ≈ {delta:.3e}\n");

    // ── 2. RNS basis ───────────────────────────────────────────────────────────
    // Four 31-bit NTT-friendly primes.  After multiplying two ciphertexts at
    // logp=30, rescaling drops the last prime (≈30 bits), leaving three primes
    // with logp back at 30.
    let primes = generate_primes(31, 4, N as u64);
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
    let rlk: RnsGadgetRelinKey<N> = engine.generate_gadget_relin_key(&sk, &mut rng);
    println!("Keys generated (secret, public, gadget-relin).\n");

    // ── 4. Input vectors ───────────────────────────────────────────────────────
    let a: [f64; 4] = [1.0, 2.0, 3.0, 4.0];
    let b: [f64; 4] = [0.5, 1.0, 1.5, 2.0];
    let expected: Vec<f64> = a.iter().zip(b.iter()).map(|(x, y)| x * y).collect();

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
    println!("Both plaintexts encrypted (logp={}, logq={}).\n", ct_a.logp, ct_a.logq);

    // ── 7. Homomorphic multiply ────────────────────────────────────────────────
    // Uses the RNS gadget relinearization key to keep the degree at 1.
    let ct_product =
        CkksEngine::<RnsPoly<N>, N>::mul_ciphertexts_gadget(&ct_a, &ct_b, &rlk);
    println!(
        "Gadget multiplication done  (logp={}, logq={}).",
        ct_product.logp, ct_product.logq
    );

    // ── 8. Rescale ─────────────────────────────────────────────────────────────
    // Drops the last RNS prime and divides coefficients by it, halving logp.
    let ct_rescaled =
        CkksEngine::<RnsPoly<N>, N>::rescale_ciphertext(&ct_product)
            .expect("rescale");
    println!(
        "Rescale done                (logp={}, logq={}).\n",
        ct_rescaled.logp, ct_rescaled.logq
    );

    // ── 9. Decrypt ─────────────────────────────────────────────────────────────
    // The rescaled ciphertext lives in a 3-prime basis; reduce the secret key to
    // match so the basis Arc::ptr_eq check in polynomial arithmetic passes.
    let reduced_basis = ct_rescaled.c0.basis().clone();
    let level = reduced_basis.channel_count();
    let sk_channels: Vec<[u64; N]> = sk.poly.channels()[..level].to_vec();
    let sk_poly_reduced =
        RnsPoly::from_channels(sk_channels, reduced_basis, false)
            .expect("reduced sk is correctly reduced");
    let sk_reduced = SecretKey { poly: sk_poly_reduced };

    let pt_product = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_rescaled, &sk_reduced);

    // ── 10. Decode ─────────────────────────────────────────────────────────────
    let decoded = encoder.decode(&pt_product);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    // ── 11. Verify ─────────────────────────────────────────────────────────────
    println!("─── Results ──────────────────────────────────────────────────────\n");
    println!(
        "{:<6} {:>8} {:>8} {:>10} {:>12} {:>10}",
        "slot", "a[i]", "b[i]", "expected", "decoded", "|err|"
    );

    let mut max_err = 0f64;
    for (i, (&exp, &got)) in expected.iter().zip(result.iter()).enumerate() {
        let err = (exp - got).abs();
        max_err = max_err.max(err);
        println!(
            "{:<6} {:>8.4} {:>8.4} {:>10.6} {:>12.6} {:>10.2e}",
            i, a[i], b[i], exp, got, err
        );
    }

    println!("\nMax absolute error : {max_err:.2e}");

    // Loose bound: encryption noise (~1e-7 each) + relin noise (~5e-8)
    // + tensor-product noise accumulation.  Use 1e-4 as a generous safety margin.
    let bound = 1e-4_f64;
    if max_err <= bound {
        println!("✓  within expected bound (≤ {bound:.0e})");
    } else {
        println!("✗  ERROR: {max_err:.2e} exceeds bound {bound:.0e}");
        std::process::exit(1);
    }

    println!("\nDone.");
}
