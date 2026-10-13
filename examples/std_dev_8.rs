//! Plan for homomorphic standard deviation over 8 CKKS slots.
//!
//! This example is intentionally a scaffold.  It documents the encrypted
//! pipeline we want to build, and keeps a plaintext reference calculation in
//! `main` so the target values are concrete before the missing homomorphic
//! helper operations are implemented.
//!
//! ## Goal
//!
//! Start with one ciphertext containing 8 real-valued slots:
//!
//! ```text
//! ct = Enc([x0, x1, x2, x3, x4, x5, x6, x7])
//! ```
//!
//! Produce one ciphertext where every slot contains the same standard
//! deviation:
//!
//! ```text
//! Enc([std_dev, std_dev, std_dev, std_dev,
//!      std_dev, std_dev, std_dev, std_dev])
//! ```
//!
//! ## Plain formula
//!
//! ```text
//! mean = (x0 + x1 + ... + x7) / 8
//!
//! variance = ((x0 - mean)^2
//!           + (x1 - mean)^2
//!           + ...
//!           + (x7 - mean)^2) / 8
//!
//! std_dev = sqrt(variance)
//! ```
//!
//! ## Encrypted pipeline
//!
//! 1. Encrypt 8 input values into one ciphertext.
//!
//! ```text
//! ct = Enc([x0, x1, x2, x3, x4, x5, x6, x7])
//! ```
//!
//! 2. Sum all slots with a binary rotation tree.
//!
//! For 8 slots, the required rotations are 1, 2, and 4:
//!
//! ```text
//! sum_ct = ct
//! sum_ct = sum_ct + rotate(sum_ct, 1)
//! sum_ct = sum_ct + rotate(sum_ct, 2)
//! sum_ct = sum_ct + rotate(sum_ct, 4)
//! ```
//!
//! After this step:
//!
//! ```text
//! sum_ct = Enc([sum, sum, sum, sum, sum, sum, sum, sum])
//! sum = x0 + x1 + ... + x7
//! ```
//!
//! Required operations:
//!
//! ```text
//! rotate_ciphertext(ct, offset)
//! add_ciphertexts(ct1, ct2)
//! rotation keys for offsets 1, 2, 4
//! ```
//!
//! 3. Compute the encrypted mean.
//!
//! ```text
//! mean_ct = mul_plain_scalar(sum_ct, 1.0 / 8.0)
//! ```
//!
//! After this step:
//!
//! ```text
//! mean_ct = Enc([mean, mean, mean, mean, mean, mean, mean, mean])
//! ```
//!
//! Required operations:
//!
//! ```text
//! mul_plain_scalar(ct, scalar)
//! rescale
//! ```
//!
//! 4. Subtract the mean from the original input.
//!
//! ```text
//! centered_ct = sub_ciphertexts(ct, mean_ct)
//! ```
//!
//! After this step:
//!
//! ```text
//! centered_ct = Enc([x0 - mean, x1 - mean, ..., x7 - mean])
//! ```
//!
//! Required operations:
//!
//! ```text
//! neg_ciphertext(ct)
//! add_ciphertexts(ct1, ct2)
//! sub_ciphertexts(ct1, ct2)
//! level/scale alignment before subtraction
//! ```
//!
//! 5. Square each centered value.
//!
//! ```text
//! sq_diff_ct = mul_ciphertexts_gadget(centered_ct, centered_ct, rlk)
//! sq_diff_ct = rescale(sq_diff_ct)
//! ```
//!
//! After this step:
//!
//! ```text
//! sq_diff_ct = Enc([
//!   (x0 - mean)^2,
//!   (x1 - mean)^2,
//!   ...
//!   (x7 - mean)^2
//! ])
//! ```
//!
//! Required operations:
//!
//! ```text
//! mul_ciphertexts_gadget(ct1, ct2, rlk)
//! rescale
//! ```
//!
//! 6. Sum squared deviations with rotations 1, 2, and 4.
//!
//! ```text
//! var_sum_ct = sq_diff_ct
//! var_sum_ct = var_sum_ct + rotate(var_sum_ct, 1)
//! var_sum_ct = var_sum_ct + rotate(var_sum_ct, 2)
//! var_sum_ct = var_sum_ct + rotate(var_sum_ct, 4)
//! ```
//!
//! After this step:
//!
//! ```text
//! var_sum_ct = Enc([sum_sq_diff, sum_sq_diff, ..., sum_sq_diff])
//! sum_sq_diff = (x0 - mean)^2 + ... + (x7 - mean)^2
//! ```
//!
//! 7. Divide by 8 to get encrypted variance.
//!
//! ```text
//! variance_ct = mul_plain_scalar(var_sum_ct, 1.0 / 8.0)
//! ```
//!
//! This is the first implementation milestone:
//!
//! ```text
//! variance_ct = Enc([variance, variance, ..., variance])
//! ```
//!
//! 8. Later, approximate square root.
//!
//! ```text
//! std_dev_ct = eval_poly(variance_ct, sqrt_coeffs)
//! ```
//!
//! The square-root approximation should be decided separately.  Bounding
//! inputs to a small interval helps, but `sqrt(x)` remains difficult near
//! zero, so the polynomial must be fitted for the expected variance range.
//!
//! ## Minimum operation set for encrypted variance
//!
//! ```text
//! add_ciphertexts
//! neg_ciphertext
//! sub_ciphertexts
//! rotate_ciphertext
//! sum_slots_power_of_two
//! mul_plain_scalar
//! mul_ciphertexts_gadget
//! rescale
//! level/scale alignment helper
//! ```
//!
//! ## Depth budget
//!
//! ```text
//! mean = sum / 8                1 level
//! squared deviations            1 level
//! variance = sum_sq_diff / 8    1 level
//! ```
//!
//! Encrypted variance needs about 3 levels.  A degree-3 polynomial square
//! root would add about 2 more levels, for about 5 levels total.

use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::engine::{CkksEngine, CkksParams, RnsGadgetRelinKey};
use toy_heaan_ckks::crypto::pipeline::{
    RotationKeyMap, align_level, decode_wide, generate_rotation_key_set,
    mul_plain_scalar, reduce_secret_key, sum_slots,
};
use toy_heaan_ckks::crypto::types::Ciphertext;
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::keys::SecretKey;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

const N: usize = 16;
const SLOTS: usize = 8;
const SCALE: u32 = 61;
const SCALAR_SCALE: u32 = SCALE;
const NUM_PRIMES: usize = 5;
const INPUT: [f64; SLOTS] = [0.12, 0.18, 0.21, 0.27, 0.31, 0.38, 0.42, 0.49];

type Ct = Ciphertext<RnsPoly<N>, N>;
type Sk = SecretKey<RnsPoly<N>, N>;
type Basis = Arc<RnsBasis<N>>;
type RotKeyMap = RotationKeyMap<N>;

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

fn make_rotkeys(
    engine: &CkksEngine<RnsPoly<N>, N>,
    sk: &Sk,
    rng: &mut ChaCha20Rng,
) -> RotKeyMap {
    generate_rotation_key_set(engine, sk, SLOTS, rng)
}

fn mul_and_rescale(ct_a: &Ct, ct_b: &Ct, rlk: &RnsGadgetRelinKey<N>) -> Ct {
    let product =
        CkksEngine::<RnsPoly<N>, N>::mul_ciphertexts_gadget(ct_a, ct_b, rlk);
    CkksEngine::<RnsPoly<N>, N>::rescale_ciphertext(&product)
        .expect("mul_and_rescale: rescale")
}

fn main() {
    let mean = INPUT.iter().sum::<f64>() / INPUT.len() as f64;
    let sum = INPUT.iter().sum::<f64>();
    let variance = INPUT
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>()
        / INPUT.len() as f64;
    let std_dev = variance.sqrt();

    println!("CKKS std-dev over 8 slots: encrypted variance demo");
    println!("input     : {INPUT:?}");
    println!("sum       : {sum:.8}");
    println!("mean      : {mean:.8}");
    println!("variance  : {variance:.8}");
    println!("std-dev   : {std_dev:.8}");
    println!();

    let primes = generate_primes(SCALE as usize, NUM_PRIMES, N as u64);
    let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
    let logq = basis.total_bits();
    let engine = make_engine(basis.clone());
    let encoder = CkksEncoder::<N>::new(SCALE);
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    let sk = engine.generate_secret_key(&mut rng).expect("sk");
    let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");
    let rotkeys = make_rotkeys(&engine, &sk, &mut rng);

    let pt_input = encoder.encode(&INPUT, basis);
    let ct_input = engine.encrypt(&pt_input, &pk, logq, &mut rng);
    let ct_sum = sum_slots(&ct_input, SLOTS, &rotkeys);
    println!(
        "After sum_slots: logp={}, logq={}, primes={}",
        ct_sum.logp,
        ct_sum.logq,
        ct_sum.c0.basis().channel_count()
    );
    let pt_sum = CkksEngine::<RnsPoly<N>, N>::decrypt(
        &ct_sum,
        &reduce_secret_key(&sk, ct_sum.c0.basis().clone()),
    );
    let decoded_sum = decode_wide(&pt_sum);
    println!(
        "  decoded sum (wide): {:.8}, expected {sum:.8}, err={:.2e}",
        decoded_sum[0],
        (decoded_sum[0] - sum).abs()
    );
    println!();

    let ct_mean = mul_plain_scalar(&ct_sum, 1.0 / SLOTS as f64, SCALAR_SCALE);
    println!(
        "After mean = sum * 1/{SLOTS}: logp={}, logq={}, primes={}",
        ct_mean.logp,
        ct_mean.logq,
        ct_mean.c0.basis().channel_count()
    );
    let pt_mean = CkksEngine::<RnsPoly<N>, N>::decrypt(
        &ct_mean,
        &reduce_secret_key(&sk, ct_mean.c0.basis().clone()),
    );
    let decoded_mean = decode_wide(&pt_mean);
    println!(
        "  decoded mean (wide): {:.8}, expected {mean:.8}, err={:.2e}",
        decoded_mean[0],
        (decoded_mean[0] - mean).abs()
    );
    println!();

    let ct_input_mean_level =
        align_level(&ct_input, ct_mean.c0.basis().clone(), ct_mean.logq);
    let ct_centered = CkksEngine::<RnsPoly<N>, N>::sub_ciphertexts(
        &ct_input_mean_level,
        &ct_mean,
    );
    println!(
        "After centered = input - mean: logp={}, logq={}, primes={}",
        ct_centered.logp,
        ct_centered.logq,
        ct_centered.c0.basis().channel_count()
    );

    let sk_centered = reduce_secret_key(&sk, ct_centered.c0.basis().clone());
    let engine_centered = make_engine(ct_centered.c0.basis().clone());
    let rlk_centered =
        engine_centered.generate_gadget_relin_key(&sk_centered, &mut rng);
    let ct_sq_diff = mul_and_rescale(&ct_centered, &ct_centered, &rlk_centered);
    println!(
        "After squared differences: logp={}, logq={}, primes={}",
        ct_sq_diff.logp,
        ct_sq_diff.logq,
        ct_sq_diff.c0.basis().channel_count()
    );

    let sk_sq = reduce_secret_key(&sk, ct_sq_diff.c0.basis().clone());
    let engine_sq = make_engine(ct_sq_diff.c0.basis().clone());
    let rotkeys_sq = make_rotkeys(&engine_sq, &sk_sq, &mut rng);
    let ct_var_sum = sum_slots(&ct_sq_diff, SLOTS, &rotkeys_sq);
    let ct_variance =
        mul_plain_scalar(&ct_var_sum, 1.0 / SLOTS as f64, SCALAR_SCALE);
    println!(
        "After variance = sum_sq_diff * 1/{SLOTS}: logp={}, logq={}, primes={}",
        ct_variance.logp,
        ct_variance.logq,
        ct_variance.c0.basis().channel_count()
    );

    let sk_variance = reduce_secret_key(&sk, ct_variance.c0.basis().clone());
    let pt_variance =
        CkksEngine::<RnsPoly<N>, N>::decrypt(&ct_variance, &sk_variance);
    let decoded_variance = decode_wide(&pt_variance);

    let max_variance_err = decoded_variance
        .iter()
        .take(SLOTS)
        .map(|got| (got - variance).abs())
        .fold(0.0f64, f64::max);

    println!("Encrypted variance result:");
    println!("  plaintext variance : {variance:.8}");
    println!("  plaintext std-dev  : {std_dev:.8}");
    for (slot, got) in decoded_variance.iter().take(SLOTS).enumerate() {
        println!(
            "  slot {slot}: decoded={got:.8}, expected={variance:.8}, err={:.2e}",
            (got - variance).abs()
        );
    }
    println!("max variance error: {max_variance_err:.2e}");
}
