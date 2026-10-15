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
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

/// Ring degree: Z[X]/(X^N + 1), giving N/2 = 512 real slots.
const N: usize = 1024;

/// Scale for large-prime tests (≈62-bit primes, 2 primes fit in u128 → 1 level).
const SCALE_LARGE: u32 = 62;

/// Scale for chained-op tests (40-bit primes, 3 primes fit in u128 → 2 levels).
const SCALE_CHAIN: u32 = 40;

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Build a `CkksEngine` at `basis` with standard noise parameters.
fn make_engine(
    basis: Arc<RnsBasis<N>>,
    scale_bits: u32,
) -> CkksEngine<RnsPoly<N>, N> {
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
    ct_a: &toy_heaan_ckks::ckks::types::Ciphertext<RnsPoly<N>, N>,
    ct_b: &toy_heaan_ckks::ckks::types::Ciphertext<RnsPoly<N>, N>,
    rlk: &toy_heaan_ckks::ckks::engine::RnsGadgetRelinKey<N>,
) -> toy_heaan_ckks::ckks::types::Ciphertext<RnsPoly<N>, N> {
    ct_a.mul(ct_b, rlk).rescale().expect("rescale")
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
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, &mut rng);

    let ct_rs = mul_and_rescale(&ct_a, &ct_b, &rlk);

    let sk_l1 = sk.mod_switch_to(ct_rs.level()).expect("sk level");
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
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, &mut rng);

    let ct_ab = mul_and_rescale(&ct_a, &ct_b, &rlk);

    // ── Derive level-2 keys from the reduced basis ────────────────────────────
    let basis_l2 = ct_ab.c0.basis().clone();
    let sk_l2 = sk.mod_switch_to(ct_ab.level()).expect("sk level");
    let engine_l2 = make_engine(basis_l2.clone(), SCALE_CHAIN);
    let pk_l2 = engine_l2
        .generate_public_key(&sk_l2, &mut rng)
        .expect("keygen l2");
    let rlk_l2 = engine_l2.generate_gadget_relin_key(&sk_l2, &mut rng);

    // ── Second multiplication: (a×b) × c ─────────────────────────────────────
    let pt_c = encoder.encode(&c, basis_l2.clone());
    let ct_c = engine_l2.encrypt(&pt_c, &pk_l2, &mut rng);

    let ct_abc = mul_and_rescale(&ct_ab, &ct_c, &rlk_l2);

    // ── Decrypt and verify ────────────────────────────────────────────────────
    let sk_l1 = sk.mod_switch_to(ct_abc.level()).expect("sk level");
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
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, &mut rng);
    let ct_c = engine.encrypt(&pt_c, &pk, &mut rng);

    // Addition does not change the level.
    let ct_sum = &ct_a + &ct_b;

    let ct_rs = mul_and_rescale(&ct_sum, &ct_c, &rlk);

    let sk_l1 = sk.mod_switch_to(ct_rs.level()).expect("sk level");
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_rs, &sk_l1);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err = max_abs_err(&expected, &result);
    assert!(err < 1e-8, "(a+b)×c: max error {err:.2e} exceeds 1e-8");
}

// ── Test 4: subtraction ──────────────────────────────────────────────────────

/// Verify ciphertext subtraction decrypts to element-wise differences without
/// consuming a level.
#[test]
fn subtract_ciphertexts() {
    let primes = generate_primes(SCALE_LARGE as usize, 2, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));

    let engine = make_engine(basis.clone(), SCALE_LARGE);
    let encoder = CkksEncoder::<N>::new(SCALE_LARGE);
    let mut rng = ChaCha20Rng::seed_from_u64(6);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");

    let a: Vec<f64> = vec![0.8, -0.3, 0.5, -0.7, 0.2, 0.9, -0.4, 0.1];
    let b: Vec<f64> = vec![0.1, 0.4, -0.2, -0.5, 0.7, -0.6, -0.3, 0.8];
    let expected: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x - y).collect();

    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, &mut rng);

    let ct_diff = &ct_a - &ct_b;
    assert_eq!(ct_diff.logp, ct_a.logp);
    assert_eq!(ct_diff.level(), ct_a.level());

    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_diff, &sk);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err = max_abs_err(&expected, &result);
    assert!(err < 1e-8, "a-b: max error {err:.2e} exceeds 1e-8");
}

// ── Test 5: multiplication followed by addition ───────────────────────────────

/// Verify that (a × b) + c decrypts correctly with 40-bit primes.
///
/// After multiplying a and b and rescaling (level drops from 3 to 2), c is
/// encrypted fresh at that level so logp and the basis match for the
/// homomorphic addition.
#[test]
fn multiply_then_add() {
    // Three 40-bit primes: one consumed by rescale, two remain for the addition.
    let primes = generate_primes(SCALE_CHAIN as usize, 3, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));

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
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, &mut rng);

    let ct_ab = mul_and_rescale(&ct_a, &ct_b, &rlk);

    // ── Encrypt c at the post-rescale level ───────────────────────────────────
    let basis_l2 = ct_ab.c0.basis().clone();
    let sk_l2 = sk.mod_switch_to(ct_ab.level()).expect("sk level");
    let engine_l2 = make_engine(basis_l2.clone(), SCALE_CHAIN);
    let pk_l2 = engine_l2
        .generate_public_key(&sk_l2, &mut rng)
        .expect("keygen l2");

    let pt_c = encoder.encode(&c, basis_l2.clone());
    let ct_c = engine_l2.encrypt(&pt_c, &pk_l2, &mut rng);

    // ── Add: (a×b) + c ───────────────────────────────────────────────────────
    // Both ciphertexts are at level 2 with logp=SCALE_CHAIN.
    let ct_sum = &ct_ab + &ct_c;

    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_sum, &sk_l2);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err = max_abs_err(&expected, &result);
    assert!(err < 1e-4, "(a×b)+c: max error {err:.2e} exceeds 1e-4");
}

// ── Test 7: ciphertext x plaintext multiplication ─────────────────────────────

/// Verify `mul_plain` (no relinearization / key switch) decrypts to
/// element-wise products within the expected noise bound, and that its error
/// is no worse than the ct x ct gadget path on the same inputs.
#[test]
fn mul_plain_matches_ct_mul() {
    let primes = generate_primes(SCALE_LARGE as usize, 2, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));

    let engine = make_engine(basis.clone(), SCALE_LARGE);
    let encoder = CkksEncoder::<N>::new(SCALE_LARGE);
    let mut rng = ChaCha20Rng::seed_from_u64(7);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    let a: Vec<f64> = vec![0.5, -0.25, 0.75, -0.125, 0.9, -0.6, 0.3, -0.8];
    let b: Vec<f64> = vec![0.4, 0.8, -0.2, 0.6, -0.5, 0.35, -0.7, 0.15];
    let expected: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();

    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);

    // ct x pt: multiply, then rescale to bring logp back to SCALE.
    let ct_pt_prod = ct_a.mul_plain(&pt_b);
    assert_eq!(ct_pt_prod.logp, ct_a.logp + pt_b.scale_bits);
    let ct_rs = ct_pt_prod.rescale().expect("rescale");

    let sk_l1 = sk.mod_switch_to(ct_rs.level()).expect("sk level");
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_rs, &sk_l1);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(a.len()).collect();

    let err_pt = max_abs_err(&expected, &result);

    // Reference: the ct x ct gadget path on identical inputs. Both rescales
    // land on the same chain basis, so the same reduced key decrypts it.
    let ct_b = engine.encrypt(&pt_b, &pk, &mut rng);
    let ct_ct = mul_and_rescale(&ct_a, &ct_b, &rlk);
    let pt_out_ct = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_ct, &sk_l1);
    let decoded_ct = encoder.decode(&pt_out_ct);
    let result_ct: Vec<f64> = decoded_ct.into_iter().take(a.len()).collect();
    let err_ct = max_abs_err(&expected, &result_ct);

    assert!(
        err_pt < 1e-8,
        "ct x pt mul: max error {err_pt:.2e} exceeds 1e-8"
    );
    assert!(
        err_pt <= err_ct * 10.0,
        "ct x pt error {err_pt:.2e} unexpectedly larger than ct x ct {err_ct:.2e}"
    );
}

// ── Test 6: full-slot single multiplication ───────────────────────────────────

/// Verify a single multiplication over all N/2 = 512 slots with large primes.
///
/// Values are drawn from a fixed RNG seed for full reproducibility.
#[test]
fn full_slots_single_multiplication() {
    let primes = generate_primes(SCALE_LARGE as usize, 2, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));

    let engine = make_engine(basis.clone(), SCALE_LARGE);
    let encoder = CkksEncoder::<N>::new(SCALE_LARGE);
    let mut rng = ChaCha20Rng::seed_from_u64(5);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    // All 512 slots filled with values in (-0.9, 0.9), fixed seed for reproducibility.
    let slots = N / 2;
    let mut val_rng = ChaCha20Rng::seed_from_u64(99);
    let a: Vec<f64> = (0..slots)
        .map(|_| val_rng.random::<f64>() * 1.8 - 0.9)
        .collect();
    let b: Vec<f64> = (0..slots)
        .map(|_| val_rng.random::<f64>() * 1.8 - 0.9)
        .collect();
    let expected: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();

    let pt_a = encoder.encode(&a, basis.clone());
    let pt_b = encoder.encode(&b, basis.clone());
    let ct_a = engine.encrypt(&pt_a, &pk, &mut rng);
    let ct_b = engine.encrypt(&pt_b, &pk, &mut rng);

    let ct_rs = mul_and_rescale(&ct_a, &ct_b, &rlk);

    let sk_l1 = sk.mod_switch_to(ct_rs.level()).expect("sk level");
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_rs, &sk_l1);
    let decoded = encoder.decode(&pt_out);
    let result: Vec<f64> = decoded.into_iter().take(slots).collect();

    let err = max_abs_err(&expected, &result);
    assert!(
        err < 1e-6,
        "full-slot mul ({slots} slots): max error {err:.2e} exceeds 1e-6"
    );
}

// ── Test 8: mod switch ────────────────────────────────────────────────────────

/// Verify `mod_switch_to` / `match_level` drop primes without changing the
/// slots or `logp`, land on the rescaled ciphertext's basis, and reject
/// levels outside the chain.
#[test]
fn mod_switch_keeps_slots_and_meets_rescaled_level() {
    let primes = generate_primes(SCALE_CHAIN as usize, 3, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));

    let engine = make_engine(basis.clone(), SCALE_CHAIN);
    let encoder = CkksEncoder::<N>::new(SCALE_CHAIN);
    let mut rng = ChaCha20Rng::seed_from_u64(8);

    let sk = engine.generate_secret_key(&mut rng).expect("keygen");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("keygen");
    let rlk = engine.generate_gadget_relin_key(&sk, &mut rng);

    let a: Vec<f64> = vec![0.6, -0.3, 0.8, -0.5, 0.4, -0.7, 0.2, -0.9];
    let ct_a = engine.encrypt(&encoder.encode(&a, basis.clone()), &pk, &mut rng);

    // Mod switch to level 2: one prime fewer, same slots, same logp.
    let ct_switched = ct_a.mod_switch_to(2).expect("mod switch");
    assert_eq!(ct_switched.level(), 2);
    assert_eq!(ct_switched.logp, ct_a.logp);
    assert!(ct_switched.logq() < ct_a.logq());

    let sk_l2 = sk.mod_switch_to(2).expect("sk level");
    let pt_out = CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_switched, &sk_l2);
    let decoded: Vec<f64> =
        encoder.decode(&pt_out).into_iter().take(a.len()).collect();
    let err = max_abs_err(&a, &decoded);
    assert!(err < 1e-6, "mod switch: max error {err:.2e} exceeds 1e-6");

    // match_level meets a rescaled ciphertext on the very same basis Arc,
    // so the two can be added.
    let ct_sq = mul_and_rescale(&ct_a, &ct_a, &rlk);
    let ct_matched = ct_a.match_level(&ct_sq).expect("match level");
    assert!(Arc::ptr_eq(ct_matched.c0.basis(), ct_sq.c0.basis()));
    let _sum = &ct_sq + &ct_matched;

    // Levels outside the chain are errors, not panics.
    assert!(ct_a.mod_switch_to(0).is_err());
    assert!(ct_a.mod_switch_to(4).is_err());
    assert!(ct_sq.match_level(&ct_a).is_err());
}
