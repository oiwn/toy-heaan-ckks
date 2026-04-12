//! CKKS encode / decode demonstration.
//!
//! Walks through the textbook CKKS encoding pipeline step by step:
//!
//!   1. Pick ring degree N and scaling factor Δ = 2^scale_bits.
//!   2. Build an RNS basis of NTT-friendly primes (p_i ≡ 1 mod 2N).
//!   3. Encode real slot values → polynomial plaintext via the special IDFT.
//!   4. Decode the plaintext back → approximate slot values via the special DFT.
//!   5. Verify the round-trip error is bounded by ≈ 1/Δ per slot.
//!
//! Reference: <https://fhetextbook.github.io/EncodingandDecoding.html>

use std::sync::Arc;

use rustfft::num_complex::Complex64;
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::RnsBasis;

/// Ring degree. Polynomial lives in Z[X]/(X^N + 1).
/// N/2 complex values (or N/2 real values) can be packed as slots.
const N: usize = 16;

/// Scaling factor exponent. Δ = 2^SCALE_BITS controls encoding precision.
/// Larger values reduce rounding error at the cost of larger coefficients.
const SCALE_BITS: u32 = 30;

fn main() {
    println!("╔═══════════════════════════════════════╗");
    println!("║   CKKS Encode / Decode  —  Demo       ║");
    println!("╚═══════════════════════════════════════╝\n");

    // ── 1. Parameters ─────────────────────────────────────────────────────────
    let delta = 2f64.powi(SCALE_BITS as i32);
    println!("Ring degree N    : {N}");
    println!("Max slots (N/2)  : {}", N / 2);
    println!(
        "Scale bits       : {SCALE_BITS}  →  Δ = 2^{SCALE_BITS} ≈ {delta:.3e}"
    );
    println!(
        "Expected per-slot rounding error ≲ 1/Δ ≈ {:.2e}\n",
        1.0 / delta
    );

    // ── 2. RNS basis ──────────────────────────────────────────────────────────
    // Each prime must satisfy p ≡ 1 (mod 2N) so Z_p contains a primitive 2N-th
    // root of unity, enabling the negacyclic NTT. Three 31-bit primes give
    // Q = p0·p1·p2 ≈ 2^93, comfortably above the scaled coefficient range.
    let primes = generate_primes(31, 3, N as u64);

    println!("RNS primes  (each ≡ 1  mod  2N = {}):", 2 * N);
    for (i, &p) in primes.iter().enumerate() {
        println!("  q{i} = {p}");
    }

    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let encoder = CkksEncoder::<N>::new(SCALE_BITS);

    // ── 3 & 4. Real-value round-trip ──────────────────────────────────────────
    println!(
        "\n─── Real-value round-trip ───────────────────────────────────────\n"
    );

    // Fill all N/2 = 8 slots.
    let originals: [f64; 8] = [
        std::f64::consts::PI,     // 3.14159…
        std::f64::consts::E,      // 2.71828…
        std::f64::consts::SQRT_2, // 1.41421…
        1.0 / 3.0,                // 0.33333…
        -1.5,
        0.5,
        -0.25,
        0.125,
    ];

    println!("Input values ({} slots):", originals.len());
    for (i, v) in originals.iter().enumerate() {
        println!("  slot[{i}] = {v:.8}");
    }

    let pt = encoder.encode(&originals, basis.clone());
    let decoded = encoder.decode(&pt);

    println!("\nDecoded values:");
    let mut max_err = 0f64;
    for (i, (orig, dec)) in originals.iter().zip(decoded.iter()).enumerate() {
        let err = (orig - dec).abs();
        max_err = max_err.max(err);
        println!("  slot[{i}]  orig={orig:.8}  decoded={dec:.8}  |err|={err:.2e}");
    }
    println!("\nMax absolute error : {max_err:.2e}");
    check_bound(max_err, delta);

    // ── Complex-value round-trip ──────────────────────────────────────────────
    println!(
        "\n─── Complex-value round-trip ────────────────────────────────────\n"
    );

    let cx_in = [
        Complex64::new(1.0, 0.5),
        Complex64::new(-0.5, 0.25),
        Complex64::new(0.0, -1.0),
        Complex64::new(0.75, -0.75),
    ];

    println!("Input values ({} slots):", cx_in.len());
    for (i, c) in cx_in.iter().enumerate() {
        println!("  slot[{i}] = ({:.6}, {:.6}i)", c.re, c.im);
    }

    let pt_cx = encoder.encode_complex(&cx_in, basis.clone());
    let cx_out = encoder.decode_complex(&pt_cx);

    println!("\nDecoded values:");
    let mut max_err_cx = 0f64;
    for (i, (orig, dec)) in cx_in.iter().zip(cx_out.iter()).enumerate() {
        let err = (orig - dec).norm();
        max_err_cx = max_err_cx.max(err);
        println!(
            "  slot[{i}]  orig=({:.6},{:.6}i)  decoded=({:.6},{:.6}i)  |err|={err:.2e}",
            orig.re, orig.im, dec.re, dec.im
        );
    }
    println!("\nMax complex error  : {max_err_cx:.2e}");
    check_bound(max_err_cx, delta);

    // ── Precision vs. scale_bits ──────────────────────────────────────────────
    println!(
        "\n─── Precision vs. scale_bits (encoding π) ──────────────────────\n"
    );
    println!(
        "{:<12} {:<18} {:<18} {}",
        "scale_bits", "Δ", "decoded π", "|error|"
    );

    let pi_slice = [std::f64::consts::PI];
    for bits in [10u32, 20, 30, 40] {
        let enc = CkksEncoder::<N>::new(bits);
        let pt = enc.encode(&pi_slice, basis.clone());
        let dec = enc.decode(&pt);
        let err = (std::f64::consts::PI - dec[0]).abs();
        println!(
            "{:<12} {:<18.3e} {:<18.10} {:.3e}",
            bits,
            2f64.powi(bits as i32),
            dec[0],
            err
        );
    }

    println!("\nDone.");
}

fn check_bound(err: f64, delta: f64) {
    let bound = 1.0 / delta;
    if err <= bound * 2.0 {
        println!("✓  within expected rounding bound (≲ 2/Δ)");
    } else {
        println!("✗  ERROR: {err:.2e} exceeds 2/Δ = {:.2e}", 2.0 * bound);
    }
}
