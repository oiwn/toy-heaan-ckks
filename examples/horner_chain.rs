//! Homomorphic evaluation of a Horner-chain over 4096 slots.
//!
//! Computes 10 operations in a row (5 multiplications + 5 additions) using
//! N = 8192 and seven 61-bit NTT-friendly primes.
//!
//! ## What is computed
//!
//! Starting from an input vector `x₀`, each iteration applies:
//!
//!   x ← x · α + β
//!
//! five times.  Written as a polynomial in the initial value:
//!
//!   x₅ = α⁵·x₀ + β·(α⁴ + α³ + α² + α + 1)
//!
//! With α = 0.8 and β = 0.1 this is a contractive affine map whose fixed
//! point is β / (1 - α) = 0.5.
//!
//! ## Why 7 primes with SCALE = 61
//!
//! In CKKS, each multiplication doubles the scale (`logp`); rescaling then
//! drops one prime and brings `logp` back down.  After k rescalings the
//! ciphertext lives in a (7 - k)-prime RNS basis.  Five multiplications
//! leave **two** primes (Q ≈ 2^122).
//!
//! The rule `SCALE = prime_bit_size` ensures that after every rescale the
//! recovered `logp` equals SCALE exactly.
//!
//! A single remaining prime (Q ≈ 2^61 = Δ) is not enough: the plaintext
//! polynomial's constant (DC) coefficient equals `mean(slots) * Δ`.  With
//! α = 0.8 and β = 0.1 the fixed point is 0.5, so after five iterations
//! the mean converges to 0.5 and the DC term ≈ Δ/2.  Since q₀ < 2^61 = Δ,
//! we have q₀/2 < Δ/2, so the DC coefficient exceeds q₀/2 and the centered
//! CRT maps it to the wrong sign — shifting every decoded slot by exactly −1.
//! Two remaining primes give Q/2 ≈ 2^121 ≫ Δ, preventing the wrap.
//!
//! ## Level management
//!
//! Every ciphertext operand (α or β) must be freshly encrypted at the
//! *current* level with the same `logq` as the running accumulator.  A
//! reduced secret key is derived from the original by truncating channels,
//! and new public / relinearization keys are generated from it before each
//! multiply.  The `logq` field is taken from the accumulator rather than
//! computed from the basis to account for the cumulative ≈1-bit drift
//! introduced by each rescaling.

use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::engine::{CkksEngine, CkksParams, RnsGadgetRelinKey};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::keys::SecretKey;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
use toy_heaan_ckks::table;

// ── compile-time constants ────────────────────────────────────────────────────

/// Ring degree.  N/2 = 4096 real slots available.
const N: usize = 8192;

/// Encoding scale.  With 61-bit primes: bits_dropped = 61 per rescale, so
/// logp returns exactly to SCALE after every mul+rescale pair.
const SCALE: u32 = 61;

/// Number of (mul + add) iterations = 10 total homomorphic operations.
const ITERS: usize = 5;

/// Primes needed: one per multiplication level, plus two to keep Q ≫ Δ at
/// the final level (a single remaining prime gives Q ≈ Δ, which is too tight
/// when mean slot values approach 0.5 — see module doc for details).
const NUM_PRIMES: usize = ITERS + 2; // 7

// ── iteration parameters ──────────────────────────────────────────────────────

/// Multiplicative scale factor applied at each iteration.
const ALPHA: f64 = 0.8;

/// Additive offset applied after each rescale.
const BETA: f64 = 0.1;

// ── helpers ───────────────────────────────────────────────────────────────────

type Ct = toy_heaan_ckks::crypto::types::Ciphertext<RnsPoly<N>, N>;
type Sk = SecretKey<RnsPoly<N>, N>;
type Basis = Arc<RnsBasis<N>>;

/// Derive a secret key valid at the level described by `basis`.
fn reduce_sk(sk: &Sk, basis: Basis) -> Sk {
    let level = basis.channel_count();
    let channels: Vec<[u64; N]> = sk.poly.channels()[..level].to_vec();
    let poly = RnsPoly::from_channels(channels, basis, false)
        .expect("reduce_sk: channels are already correctly reduced");
    SecretKey { poly }
}

/// Build a `CkksEngine` at `basis` with standard noise parameters.
fn make_engine(basis: Basis) -> CkksEngine<RnsPoly<N>, N> {
    CkksEngine::new(
        basis,
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        },
    )
}

/// Encrypt a constant (all slots = `value`) at the given `basis` and `logq`.
fn encrypt_scalar(
    value: f64,
    basis: Basis,
    logq: u32,
    encoder: &CkksEncoder<N>,
    engine: &CkksEngine<RnsPoly<N>, N>,
    pk: &toy_heaan_ckks::keys::PublicKey<RnsPoly<N>, N>,
    rng: &mut ChaCha20Rng,
) -> Ct {
    let slots = N / 2;
    let values = vec![value; slots];
    let pt = encoder.encode(&values, basis);
    engine.encrypt(&pt, pk, logq, rng)
}

// ── main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔══════════════════════════════════════════════════════════╗");
    println!(
        "║  CKKS Horner chain — {ITERS} × (mul + add) over {slots} slots  ║",
        slots = N / 2
    );
    println!("╚══════════════════════════════════════════════════════════╝\n");

    // ── parameters ───────────────────────────────────────────────────────────
    println!("Ring degree N     : {N}");
    println!("Slots (N/2)       : {}", N / 2);
    println!("Scale bits        : {SCALE}  →  Δ = 2^{SCALE}");
    println!(
        "RNS primes        : {NUM_PRIMES} × {SCALE}-bit (Q ≈ 2^{})",
        NUM_PRIMES * SCALE as usize
    );
    println!(
        "Operations        : {ITERS} × (mul·α + add·β)  =  {} total",
        2 * ITERS
    );
    println!("α (scale factor)  : {ALPHA}");
    println!("β (additive bias) : {BETA}");
    println!();

    // ── 1. RNS basis ─────────────────────────────────────────────────────────
    let primes = generate_primes(SCALE as usize, NUM_PRIMES, N as u64);
    println!("NTT-friendly primes (≡ 1 mod 2N = {}):", 2 * N);
    for (i, &p) in primes.iter().enumerate() {
        println!("  q{i} = {p}  ({}-bit)", 64 - p.leading_zeros());
    }
    let basis_top = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq_top = basis_top.total_bits();
    println!("  total logq = {logq_top}\n");

    // ── 2. Top-level engine and keys ─────────────────────────────────────────
    println!("Generating top-level keys (sk, pk, rlk) …");
    let engine_top = make_engine(basis_top.clone());
    let encoder = CkksEncoder::<N>::new(SCALE);
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    let sk_top = engine_top.generate_secret_key(&mut rng).expect("sk");
    let pk_top = engine_top
        .generate_public_key(&sk_top, &mut rng)
        .expect("pk");
    let rlk_top: RnsGadgetRelinKey<N> =
        engine_top.generate_gadget_relin_key(&sk_top, &mut rng);
    println!("  Done.\n");

    // ── 3. Encrypt the initial input ─────────────────────────────────────────
    // Each slot holds a distinct starting value: i / (N/2) ∈ (0, 1].
    let slots = N / 2;
    let x0: Vec<f64> = (0..slots).map(|i| (i + 1) as f64 / slots as f64).collect();
    println!("Encrypting input x₀ ({slots} slots, values in (0, 1]) …");
    let pt_x0 = encoder.encode(&x0, basis_top.clone());
    let mut ct_x = engine_top.encrypt(&pt_x0, &pk_top, logq_top, &mut rng);
    println!(
        "  ct_x: logp={}, logq={}, primes={}\n",
        ct_x.logp,
        ct_x.logq,
        ct_x.c0.basis().channel_count()
    );

    // Plain-text reference: follow the same chain without encryption.
    let mut x_ref: Vec<f64> = x0.clone();

    // Keys carried into the loop (start at top level, then derived per step).
    let mut sk_cur = sk_top.clone();
    let mut pk_cur = pk_top.clone();
    let mut rlk_cur = rlk_top;
    let mut engine_cur = engine_top;

    // ── 4. Horner iterations ─────────────────────────────────────────────────
    for iter in 1..=ITERS {
        println!(
            "── Iteration {iter}/{ITERS} ──────────────────────────────────────────"
        );
        let level_before = ct_x.c0.basis().channel_count();
        println!(
            "  Level before mul : {level_before} primes, logq={}",
            ct_x.logq
        );

        // ── 4a. Multiply: ct_x ← ct_x · α ───────────────────────────────────
        let ct_alpha = encrypt_scalar(
            ALPHA,
            ct_x.c0.basis().clone(),
            ct_x.logq,
            &encoder,
            &engine_cur,
            &pk_cur,
            &mut rng,
        );
        let ct_prod = CkksEngine::<RnsPoly<N>, N>::mul_ciphertexts_gadget(
            &ct_x, &ct_alpha, &rlk_cur,
        );

        // ── 4b. Rescale ───────────────────────────────────────────────────────
        ct_x = CkksEngine::<RnsPoly<N>, N>::rescale_ciphertext(&ct_prod)
            .expect("rescale");
        let level_after = ct_x.c0.basis().channel_count();
        println!(
            "  After mul+rescale: {level_after} primes, logp={}, logq={}",
            ct_x.logp, ct_x.logq
        );

        // ── 4c. Derive keys at the new (lower) level ──────────────────────────
        let basis_cur = ct_x.c0.basis().clone();
        sk_cur = reduce_sk(&sk_top, basis_cur.clone());
        engine_cur = make_engine(basis_cur.clone());
        pk_cur = engine_cur
            .generate_public_key(&sk_cur, &mut rng)
            .expect("pk");

        // ── 4d. Add: ct_x ← ct_x + β ─────────────────────────────────────────
        let ct_beta = encrypt_scalar(
            BETA,
            basis_cur.clone(),
            ct_x.logq,
            &encoder,
            &engine_cur,
            &pk_cur,
            &mut rng,
        );
        ct_x = CkksEngine::<RnsPoly<N>, N>::add_ciphertexts(&ct_x, &ct_beta);
        println!(
            "  After add β      : {level_after} primes, logp={}, logq={}",
            ct_x.logp, ct_x.logq
        );

        // Plaintext reference step.
        for v in x_ref.iter_mut() {
            *v = *v * ALPHA + BETA;
        }

        // Generate next iteration's relinearization key (skip on last iter).
        if iter < ITERS {
            rlk_cur = engine_cur.generate_gadget_relin_key(&sk_cur, &mut rng);
        }
        println!();
    }

    // ── 5. Decrypt ───────────────────────────────────────────────────────────
    let final_primes = ct_x.c0.basis().channel_count();
    println!(
        "Decrypting ({final_primes} primes remaining, Q ≈ 2^{}) …",
        final_primes as u32 * SCALE
    );
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_x, &sk_cur);

    // ── 6. Decode ────────────────────────────────────────────────────────────
    // to_coeffs() runs CRT reconstruction over 2 primes → Q ≈ 2^122, fits in u128.
    let decoded = encoder.decode(&pt_out);

    // ── 7. Verify ────────────────────────────────────────────────────────────
    println!();
    let mut t = table::new(["slot", "x₀[i]", "expected", "decoded"]);
    for i in 0..10 {
        t.add_row([
            i.to_string(),
            format!("{:.6}", x0[i]),
            format!("{:.6}", x_ref[i]),
            format!("{:.6}", decoded[i]),
        ]);
    }
    t.add_row(["…", &format!("({slots} slots total)"), "", ""]);
    println!("{t}");
    println!();

    let max_err = x_ref
        .iter()
        .zip(decoded.iter().take(slots))
        .map(|(e, g)| (e - g).abs())
        .fold(0f64, f64::max);

    let delta = 2f64.powi(SCALE as i32);
    // Each mul+rescale contributes ≈ ITERS * σ * N^0.5 / Δ of relin noise.
    // Encryption noise per slot ≈ σ * sqrt(hw * N) / Δ.
    // After 5 iterations the dominant term is still ≪ 1 / Δ.
    let bound = 1e-5_f64;

    println!("Scale Δ = 2^{SCALE} ≈ {delta:.3e}");
    println!("Max absolute error over all {slots} slots: {max_err:.3e}");
    println!("Error bound used              : {bound:.0e}");

    if max_err <= bound {
        println!("\n✓  All {slots} slots within error bound.");
    } else {
        println!("\n✗  ERROR: {max_err:.3e} exceeds bound {bound:.0e}");
        std::process::exit(1);
    }

    // Fixed-point sanity: with α=0.8, β=0.1, the map f(x)=0.8x+0.1 converges
    // to 0.5.  After 5 iterations, slots starting near 1.0 should be ≈ 0.434.
    let last_slot = decoded[slots - 1];
    println!(
        "\nSlot {}: x₀={:.4} → x₅={last_slot:.6} (expected {:.6})",
        slots - 1,
        x0[slots - 1],
        x_ref[slots - 1]
    );
    println!(
        "Fixed point of f(x)=αx+β: β/(1-α) = {:.6}",
        BETA / (1.0 - ALPHA)
    );

    println!("\nDone.");
}
