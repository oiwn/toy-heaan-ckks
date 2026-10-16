//! CKKS key generation demonstration.
//!
//! Shows the three key types used in CKKS and verifies the algebraic
//! relations that make them correct:
//!
//!   Secret key s  : ternary polynomial (coefficients in {-1, 0, 1})
//!   Public key (a, b) : RLWE sample  →  b + a·s  ≈  0  (only noise)
//!   Relin key  (a[k], b[k]) : gadget encryption of s²
//!                       →  b[k] + a[k]·s  ≈  s²·g_k (only noise)
//!
//! Reference: <https://fhetextbook.github.io/KeyGeneration.html>,
//! key switching: <https://fhetextbook.github.io/HomomorphicKeySwitching2.html>

use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::ckks::keys::DIGIT_BITS;
use toy_heaan_ckks::ckks::keys::public_key::{PublicKey, PublicKeyParams};
use toy_heaan_ckks::ckks::keys::relin_key::RelinearizationKey;
use toy_heaan_ckks::ckks::keys::secret_key::{SecretKey, SecretKeyParams};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
use toy_heaan_ckks::rings::traits::PolyRing;

/// Ring degree — polynomial lives in Z[X]/(X^N + 1).
const N: usize = 16;
/// Fraction of non-zero coefficients in the secret key.
const HAMMING_WEIGHT: usize = N / 2;
/// Standard deviation for the Gaussian error distribution.
const ERROR_STD: f64 = 3.2;

fn main() {
    println!("╔═══════════════════════════════════════╗");
    println!("║   CKKS Key Generation  —  Demo        ║");
    println!("╚═══════════════════════════════════════╝\n");

    // ── Parameters ────────────────────────────────────────────────────────────
    println!("Ring degree N       : {N}");
    println!(
        "Hamming weight h    : {HAMMING_WEIGHT}  ({} non-zero out of {N} coeffs)",
        HAMMING_WEIGHT
    );
    println!("Error std σ         : {ERROR_STD}");

    // ── RNS basis ─────────────────────────────────────────────────────────────
    // Two 20-bit NTT-friendly primes. Q = q0·q1 ≈ 2^40 is large enough to hold
    // any intermediate product that arises during key verification (≈ N·q ≈ 2^24).
    let primes = generate_primes(20, 2, N as u64);
    println!("\nRNS primes (each ≡ 1  mod  2N={}):", 2 * N);
    for (i, &p) in primes.iter().enumerate() {
        println!("  q{i} = {p}");
    }

    let basis: Arc<RnsBasis<N>> =
        Arc::new(RnsBasis::new(primes).expect("valid NTT primes"));
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    // ── Secret key ────────────────────────────────────────────────────────────
    println!(
        "\n─── Secret key ──────────────────────────────────────────────────\n"
    );

    let sk_params = SecretKeyParams::<N>::new(HAMMING_WEIGHT).unwrap();
    let sk = SecretKey::<RnsPoly<N>, N>::generate(&sk_params, &basis, &mut rng)
        .expect("secret key");

    let sk_coeffs = sk.poly.to_coeffs();
    let nonzero: Vec<(usize, i64)> = sk_coeffs
        .iter()
        .enumerate()
        .filter(|&(_, c)| *c != 0)
        .map(|(i, &c)| (i, c))
        .collect();

    println!(
        "Hamming weight      : {} (expected {})",
        nonzero.len(),
        HAMMING_WEIGHT
    );
    println!(
        "Non-zero positions  : {:?}",
        nonzero.iter().map(|(i, _)| i).collect::<Vec<_>>()
    );
    println!(
        "Non-zero values     : {:?}",
        nonzero.iter().map(|(_, v)| v).collect::<Vec<_>>()
    );
    let all_ternary = sk_coeffs.iter().all(|&c| c == -1 || c == 0 || c == 1);
    println!("All coeffs ∈ {{-1,0,1}} : {}", check(all_ternary));

    // ── Public key ────────────────────────────────────────────────────────────
    println!(
        "\n─── Public key ──────────────────────────────────────────────────\n"
    );
    println!("Algorithm: sample uniform a, small error e, set b = -(a·s) + e");
    println!("Relation : b + a·s = e  (small noise)\n");

    let pk_params = PublicKeyParams::<N>::new(ERROR_STD).unwrap();
    let pk =
        PublicKey::<RnsPoly<N>, N>::generate(&sk, &pk_params, &basis, &mut rng)
            .expect("public key");

    // Verify: b + a·s should equal the (unknown) error e — coefficients ≈ 0.
    let mut a_times_s = pk.a.clone();
    a_times_s *= &sk.poly; // a·s  (schoolbook, coeff domain)
    let mut check_pk = pk.b.clone();
    check_pk += &a_times_s; // b + a·s  ≈  e

    let residual = check_pk.to_coeffs();
    let max_pk = residual.iter().map(|c| c.unsigned_abs()).max().unwrap_or(0);
    println!("b + a·s  (should be small error):");
    println!("  coefficients : {:?}", &residual[..8]);
    println!("  ... (showing first 8 of {N})");
    println!(
        "  max |coeff|  : {max_pk}  (expect ≲ 3σ ≈ {})",
        (3.0 * ERROR_STD) as u64
    );
    println!("  relation ok  : {}", check(max_pk < 50));

    // ── Relinearization key ───────────────────────────────────────────────────
    println!(
        "\n─── Relinearization key ─────────────────────────────────────────\n"
    );
    println!(
        "Algorithm: compute s²; per gadget term k sample uniform a[k], small error e,"
    );
    println!("           set b[k] = -(a[k]·s) + e + s²·g_k");
    println!(
        "Relation : b[k] + a[k]·s = s²·g_k + e  (g_k: digit k of base 2^{DIGIT_BITS})\n"
    );

    let rk =
        RelinearizationKey::<RnsPoly<N>, N>::generate(&sk, ERROR_STD, &mut rng);
    println!(
        "Gadget terms        : {}  (one per channel x base-2^{DIGIT_BITS} digit)",
        rk.ksk.a.len()
    );

    // Compute s² for comparison.
    let mut s_sq = sk.poly.clone();
    s_sq *= &sk.poly; // s²

    // Verify every term: b[k] + a[k]·s - s²·g_k should equal e (small noise).
    let mut max_rk = 0;
    for (k, (a, b)) in rk.ksk.a.iter().zip(&rk.ksk.b).enumerate() {
        let mut check_rk = a.clone();
        check_rk *= &sk.poly; // a[k]·s
        check_rk += b; // b[k] + a[k]·s
        check_rk -= &s_sq.gadget_term(k, DIGIT_BITS); // ... - s²·g_k  ≈  e
        let residual_rk = check_rk.to_coeffs();
        if k == 0 {
            println!("b[0] + a[0]·s - s²·g_0  (should be small error):");
            println!("  coefficients : {:?}", &residual_rk[..8]);
            println!("  ... (showing first 8 of {N})");
        }
        let term_max = residual_rk.iter().map(|c| c.unsigned_abs()).max();
        max_rk = max_rk.max(term_max.unwrap_or(0));
    }
    println!(
        "  max |coeff| over all terms : {max_rk}  (expect ≲ 3σ ≈ {})",
        (3.0 * ERROR_STD) as u64
    );
    println!("  relation ok  : {}", check(max_rk < 50));

    println!("\nDone.");
}

fn check(ok: bool) -> &'static str {
    if ok { "✓" } else { "✗  FAIL" }
}
