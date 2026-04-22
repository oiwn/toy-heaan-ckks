//! Integration tests for homomorphic multiplication.
//!
//! Tests use N=1024 (moderately large ring degree) and large NTT-friendly primes
//! to exercise the full RNS gadget multiplication pipeline.
//!
//! ## Prime-size / SCALE relationship
//!
//! In CKKS, rescaling drops one prime `q` and divides coefficients by ~q,
//! reducing `logp` by `bits_dropped = floor(log2(q)) + 1`.  For the
//! post-rescale `logp` to equal the pre-multiplication `logp = SCALE` we need:
//!   `2*SCALE - bits_dropped = SCALE`  →  `bits_dropped = SCALE`
//!
//! This means each prime should be a SCALE-bit number (generated with
//! `generate_primes(SCALE, …)`).
//!
//! ## u128 constraint
//!
//! CRT reconstruction in this library uses u128, so the total modulus
//! `Q = q₀·…·q_{L-1}` must fit in 128 bits.  The prime count is therefore
//! limited to:  `L < 128 / SCALE`.
//!
//! | SCALE | max primes | max levels |
//! |-------|-----------|------------|
//! |  62   |     2     |     1      |  ← "large-prime" tests
//! |  40   |     3     |     2      |  ← "chained-op" tests
//!
//! Scenarios covered:
//!   1. Single multiplication — large (≈62-bit) primes, N=1024
//!   2. Two chained multiplications — 40-bit primes, N=1024
//!   3. Addition followed by multiplication — large primes, N=1024
//!   4. Multiplication followed by addition — 40-bit primes, N=1024
//!   5. Full-slot (N/2 = 512 slots) single multiplication — large primes, N=1024

use std::sync::Arc;

use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::keys::SecretKey;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

/// Ring degree: Z[X]/(X^N + 1), giving N/2 = 512 real slots.
const N: usize = 1024;

/// Scale for large-prime tests (≈62-bit primes, 2 primes fit in u128 → 1 level).
const SCALE_LARGE: u32 = 62;

/// Scale for chained-op tests (40-bit primes, 3 primes fit in u128 → 2 levels).
const SCALE_CHAIN: u32 = 40;

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Reduce `sk` to the channel count of `basis` (used after each rescaling).
fn reduce_sk(
    sk: &SecretKey<RnsPoly<N>, N>,
    basis: Arc<RnsBasis<N>>,
) -> SecretKey<RnsPoly<N>, N> {
    let level = basis.channel_count();
    let channels: Vec<[u64; N]> = sk.poly.channels()[..level].to_vec();
    let poly = RnsPoly::from_channels(channels, basis, false)
        .expect("reduced sk: channels are already correctly reduced");
    SecretKey { poly }
}

/// Build a `CkksEngine` at `basis` with standard noise parameters.
fn make_engine(basis: Arc<RnsBasis<N>>, scale_bits: u32) -> CkksEngine<RnsPoly<N>, N> {
    CkksEngine::new(
        basis,
        CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits,
        },
    )
}

/// One homomorphic multiplication + rescale.
fn mul_and_rescale(
    ct_a: &toy_heaan_ckks::crypto::types::Ciphertext<RnsPoly<N>, N>,
    ct_b: &toy_heaan_ckks::crypto::types::Ciphertext<RnsPoly<N>, N>,
    rlk: &toy_heaan_ckks::crypto::engine::RnsGadgetRelinKey<N>,
) -> toy_heaan_ckks::crypto::types::Ciphertext<RnsPoly<N>, N> {
    let ct_prod = CkksEngine::<RnsPoly<N>, N>::mul_ciphertexts_gadget(ct_a, ct_b, rlk);
    CkksEngine::<RnsPoly<N>, N>::rescale_ciphertext(&ct_prod).expect("rescale")
}

/// Maximum absolute element-wise error between two equal-length slices.
fn max_abs_err(expected: &[f64], actual: &[f64]) -> f64 {
    expected
        .iter()
        .zip(actual)
        .map(|(e, g)| (e - g).abs())
        .fold(0f64, f64::max)
}

// ── Test 1: single multiplication, large primes ───────────────────────────────

/// Verify one homomorphic multiplication with ≈62-bit NTT primes decrypts
/// to element-wise products within the expected noise bound.
///
/// Two 62-bit primes → Q ≈ 2^124 (fits in u128); one prime consumed by
/// rescaling, one remains for decryption.
#[test]
fn single_multiplication_large_primes() {
    // Two 62-bit primes (SCALE=62 → bits_dropped=62, logp preserved after rescale).
    let primes = generate_primes(SCALE_LARGE as usize, 2, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();

    let engine = make_engine(basis.clone(), SCALE_LARGE);
    let encoder = CkksEncoder::<N>::new(SCALE_LARGE);
    let mut rng = ChaCha20Rng::seed_from_u64(1);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    let a: Vec<f64> = vec![0.5, -0.25, 0.75, -0.125, 0.9, -0.6, 0.3, -0.8];
    let b: Vec<f64> = vec![0.4, 0.8, -0.2, 0.6, -0.5, 0.35, -0.7, 0.15];
    let expected: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();

    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, logq, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, logq, &mut rng);

    let ct_rs = mul_and_rescale(&ct_a, &ct_b, &rlk);

    let sk_l1 = reduce_sk(&sk, ct_rs.c0.basis().clone());
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_rs, &sk_l1);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err = max_abs_err(&expected, &result);
    assert!(
        err < 1e-8,
        "single mul (large primes): max error {err:.2e} exceeds 1e-8"
    );
}

// ── Test 2: two chained multiplications ───────────────────────────────────────

/// Verify that two sequential multiplications (each followed by rescaling)
/// yield element-wise triple products within the expected accumulated noise.
///
/// Three 40-bit primes → Q ≈ 2^120 (fits in u128); two primes are consumed
/// by the two rescalings, one remains for decryption.
///
/// After the first rescale the level drops by one; a fresh keypair at that
/// level is used to encrypt the third operand so the polynomial bases match.
#[test]
fn two_chained_multiplications() {
    // Three 40-bit primes (SCALE=40 → bits_dropped=40; 2 levels of multiplication).
    let primes = generate_primes(SCALE_CHAIN as usize, 3, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();

    let engine = make_engine(basis.clone(), SCALE_CHAIN);
    let encoder = CkksEncoder::<N>::new(SCALE_CHAIN);
    let mut rng = ChaCha20Rng::seed_from_u64(2);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    // Keep inputs in (0, 1] so products stay well within the noise budget.
    let a: Vec<f64> = vec![0.9, 0.5, 0.8, 0.3, 0.7, 0.4, 0.6, 0.2];
    let b: Vec<f64> = vec![0.8, 0.6, 0.4, 0.9, 0.5, 0.7, 0.3, 0.85];
    let c: Vec<f64> = vec![0.7, 0.9, 0.3, 0.5, 0.6, 0.8, 0.4, 0.1];
    let expected: Vec<f64> = a
        .iter()
        .zip(&b)
        .zip(&c)
        .map(|((x, y), z)| x * y * z)
        .collect();

    // ── First multiplication: a × b ───────────────────────────────────────────
    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, logq, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, logq, &mut rng);

    let ct_ab = mul_and_rescale(&ct_a, &ct_b, &rlk);

    // ── Derive level-2 keys from the reduced basis ────────────────────────────
    // ct_ab.logq carries the rescaling adjustment; use it (not basis.total_bits())
    // so that the fresh ciphertext's logq field matches exactly.
    let basis_l2 = ct_ab.c0.basis().clone();
    let sk_l2 = reduce_sk(&sk, basis_l2.clone());
    let engine_l2 = make_engine(basis_l2.clone(), SCALE_CHAIN);
    let pk_l2 = engine_l2.generate_public_key(&sk_l2, &mut rng).expect("keygen l2");
    let rlk_l2 = engine_l2.generate_gadget_relin_key(&sk_l2, &mut rng);

    // ── Second multiplication: (a×b) × c ─────────────────────────────────────
    let pt_c = encoder.encode(&c, basis_l2.clone());
    let ct_c = engine_l2.encrypt(&pt_c, &pk_l2, ct_ab.logq, &mut rng);

    let ct_abc = mul_and_rescale(&ct_ab, &ct_c, &rlk_l2);

    // ── Decrypt and verify ────────────────────────────────────────────────────
    let sk_l1 = reduce_sk(&sk, ct_abc.c0.basis().clone());
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_abc, &sk_l1);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err = max_abs_err(&expected, &result);
    assert!(
        err < 1e-4,
        "chained mul (a×b×c): max error {err:.2e} exceeds 1e-4"
    );
}

// ── Test 3: addition followed by multiplication ───────────────────────────────

/// Verify that (a + b) × c decrypts correctly with large (≈62-bit) primes.
///
/// Homomorphic addition leaves the ciphertext level unchanged, so the sum and
/// c share the same basis when the multiplication is performed.
#[test]
fn add_then_multiply() {
    let primes = generate_primes(SCALE_LARGE as usize, 2, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();

    let engine = make_engine(basis.clone(), SCALE_LARGE);
    let encoder = CkksEncoder::<N>::new(SCALE_LARGE);
    let mut rng = ChaCha20Rng::seed_from_u64(3);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    let a: Vec<f64> = vec![0.3, -0.4, 0.6, -0.2, 0.8, -0.1, 0.5, -0.7];
    let b: Vec<f64> = vec![-0.1, 0.5, -0.3, 0.7, -0.4, 0.6, -0.2, 0.4];
    let c: Vec<f64> = vec![0.9, 0.7, 0.5, 0.3, 0.8, 0.6, 0.4, 0.2];
    let expected: Vec<f64> = a
        .iter()
        .zip(&b)
        .zip(&c)
        .map(|((x, y), z)| (x + y) * z)
        .collect();

    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let pt_c = encoder.encode(&c, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, logq, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, logq, &mut rng);
    let ct_c = engine.encrypt(&pt_c, &pk, logq, &mut rng);

    // Addition does not change the level.
    let ct_sum = CkksEngine::<RnsPoly<N>, N>::add_ciphertexts(&ct_a, &ct_b);

    let ct_rs = mul_and_rescale(&ct_sum, &ct_c, &rlk);

    let sk_l1 = reduce_sk(&sk, ct_rs.c0.basis().clone());
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_rs, &sk_l1);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err = max_abs_err(&expected, &result);
    assert!(
        err < 1e-8,
        "(a+b)×c: max error {err:.2e} exceeds 1e-8"
    );
}

// ── Test 4: multiplication followed by addition ───────────────────────────────

/// Verify that (a × b) + c decrypts correctly with 40-bit primes.
///
/// After multiplying a and b and rescaling (level drops from 3 to 2), c is
/// encrypted fresh at that level so the logp and logq fields match for the
/// homomorphic addition.
#[test]
fn multiply_then_add() {
    // Three 40-bit primes: one consumed by rescale, two remain for the addition.
    let primes = generate_primes(SCALE_CHAIN as usize, 3, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();

    let engine = make_engine(basis.clone(), SCALE_CHAIN);
    let encoder = CkksEncoder::<N>::new(SCALE_CHAIN);
    let mut rng = ChaCha20Rng::seed_from_u64(4);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    let a: Vec<f64> = vec![0.6, -0.3, 0.8, -0.5, 0.4, -0.7, 0.2, -0.9];
    let b: Vec<f64> = vec![0.5, 0.7, 0.3, 0.9, 0.6, 0.4, 0.8, 0.1];
    let c: Vec<f64> = vec![0.1, -0.2, 0.4, -0.3, 0.7, -0.5, 0.3, -0.6];
    let expected: Vec<f64> = a
        .iter()
        .zip(&b)
        .zip(&c)
        .map(|((x, y), z)| x * y + z)
        .collect();

    // ── Multiply a × b, then rescale ─────────────────────────────────────────
    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, logq, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, logq, &mut rng);

    let ct_ab = mul_and_rescale(&ct_a, &ct_b, &rlk);

    // ── Encrypt c at the post-rescale level ───────────────────────────────────
    // Use ct_ab.logq (not basis_l2.total_bits()) to ensure logq fields match.
    let basis_l2 = ct_ab.c0.basis().clone();
    let sk_l2 = reduce_sk(&sk, basis_l2.clone());
    let engine_l2 = make_engine(basis_l2.clone(), SCALE_CHAIN);
    let pk_l2 = engine_l2.generate_public_key(&sk_l2, &mut rng).expect("keygen l2");

    let pt_c = encoder.encode(&c, basis_l2.clone());
    let ct_c = engine_l2.encrypt(&pt_c, &pk_l2, ct_ab.logq, &mut rng);

    // ── Add: (a×b) + c ───────────────────────────────────────────────────────
    // Both ciphertexts are at level 2 with logp=SCALE_CHAIN and logq=ct_ab.logq.
    let ct_sum = CkksEngine::<RnsPoly<N>, N>::add_ciphertexts(&ct_ab, &ct_c);

    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_sum, &sk_l2);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err = max_abs_err(&expected, &result);
    assert!(
        err < 1e-4,
        "(a×b)+c: max error {err:.2e} exceeds 1e-4"
    );
}

// ── Test 5: full-slot single multiplication ───────────────────────────────────

/// Verify a single multiplication over all N/2 = 512 slots with large primes.
///
/// Values are drawn from a fixed RNG seed for full reproducibility.
#[test]
fn full_slots_single_multiplication() {
    let primes = generate_primes(SCALE_LARGE as usize, 2, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();

    let engine = make_engine(basis.clone(), SCALE_LARGE);
    let encoder = CkksEncoder::<N>::new(SCALE_LARGE);
    let mut rng = ChaCha20Rng::seed_from_u64(5);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    // All 512 slots filled with values in (-0.9, 0.9), fixed seed for reproducibility.
    let slots = N / 2;
    let mut val_rng = ChaCha20Rng::seed_from_u64(99);
    let a: Vec<f64> = (0..slots).map(|_| val_rng.random::<f64>() * 1.8 - 0.9).collect();
    let b: Vec<f64> = (0..slots).map(|_| val_rng.random::<f64>() * 1.8 - 0.9).collect();
    let expected: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();

    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, logq, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, logq, &mut rng);

    let ct_rs = mul_and_rescale(&ct_a, &ct_b, &rlk);

    let sk_l1 = reduce_sk(&sk, ct_rs.c0.basis().clone());
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_rs, &sk_l1);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(slots).collect();

    let err = max_abs_err(&expected, &result);
    assert!(
        err < 1e-6,
        "full-slot mul ({slots} slots): max error {err:.2e} exceeds 1e-6"
    );
}
