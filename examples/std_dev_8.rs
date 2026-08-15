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

use std::collections::HashMap;
use std::sync::Arc;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use rustfft::num_complex::Complex64;
use toy_heaan_ckks::crypto::engine::{
    CkksEngine, CkksParams, RnsGadgetRelinKey, RnsGadgetRotationKey,
};
use toy_heaan_ckks::crypto::types::{Ciphertext, Plaintext};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::encoding::special_fft::{VandermondeTables, special_dft};
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
type RotKeyMap = HashMap<i32, RnsGadgetRotationKey<N>>;

fn reduce_sk(sk: &Sk, basis: Basis) -> Sk {
    let level = basis.channel_count();
    let channels: Vec<[u64; N]> = sk.poly.channels()[..level].to_vec();
    let poly = RnsPoly::from_channels(channels, basis, false)
        .expect("reduce_sk: channels are already reduced");
    SecretKey { poly }
}

fn reduce_ct_to_basis(ct: &Ct, basis: Basis, logq: u32) -> Ct {
    let level = basis.channel_count();
    let c0_channels: Vec<[u64; N]> = ct.c0.channels()[..level].to_vec();
    let c1_channels: Vec<[u64; N]> = ct.c1.channels()[..level].to_vec();
    Ciphertext {
        c0: RnsPoly::from_channels(
            c0_channels,
            basis.clone(),
            ct.c0.is_ntt_domain(),
        )
        .expect("reduce c0: channels are already reduced"),
        c1: RnsPoly::from_channels(c1_channels, basis, ct.c1.is_ntt_domain())
            .expect("reduce c1: channels are already reduced"),
        logp: ct.logp,
        logq,
    }
}

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
    [1, 2, 4]
        .into_iter()
        .map(|rotation| {
            (
                rotation,
                engine.generate_gadget_rotation_key(sk, rotation, rng),
            )
        })
        .collect()
}

fn sum_slots_power_of_two(ct: &Ct, slots: usize, rotkeys: &RotKeyMap) -> Ct {
    assert!(slots.is_power_of_two(), "slots must be a power of two");
    assert!(slots <= N / 2, "slots cannot exceed available CKKS slots");

    let mut acc = Ciphertext {
        c0: ct.c0.clone(),
        c1: ct.c1.clone(),
        logp: ct.logp,
        logq: ct.logq,
    };

    let mut step = 1usize;
    while step < slots {
        let rotation = i32::try_from(step).expect("rotation offset fits in i32");
        let rotkey = rotkeys
            .get(&rotation)
            .expect("missing rotation key for sum_slots_power_of_two");
        let rotated = CkksEngine::<RnsPoly<N>, N>::rotate_ciphertext(&acc, rotkey);
        acc = CkksEngine::<RnsPoly<N>, N>::add_ciphertexts(&acc, &rotated);
        step *= 2;
    }

    acc
}

fn mul_plain_scalar(ct: &Ct, scalar: f64, scalar_scale_bits: u32) -> Ct {
    let scalar_encoder = CkksEncoder::<N>::new(scalar_scale_bits);
    let scalar_values = vec![scalar; SLOTS];
    let pt_scalar = scalar_encoder.encode(&scalar_values, ct.c0.basis().clone());

    let mut c0 = ct.c0.clone();
    c0 *= &pt_scalar.poly;

    let mut c1 = ct.c1.clone();
    c1 *= &pt_scalar.poly;

    let product = Ciphertext {
        c0,
        c1,
        logp: ct.logp + scalar_scale_bits,
        logq: ct.logq,
    };

    CkksEngine::<RnsPoly<N>, N>::rescale_ciphertext(&product)
        .expect("mul_plain_scalar: rescale")
}

fn mul_and_rescale(ct_a: &Ct, ct_b: &Ct, rlk: &RnsGadgetRelinKey<N>) -> Ct {
    let product =
        CkksEngine::<RnsPoly<N>, N>::mul_ciphertexts_gadget(ct_a, ct_b, rlk);
    CkksEngine::<RnsPoly<N>, N>::rescale_ciphertext(&product)
        .expect("mul_and_rescale: rescale")
}

fn decode_wide(pt: &Plaintext<RnsPoly<N>, N>) -> Vec<f64> {
    let mut poly = pt.poly.clone();
    poly.to_coeff_domain();

    let coeffs: Vec<Complex64> = (0..N)
        .map(|idx| {
            let residues: Vec<u64> =
                poly.channels().iter().map(|channel| channel[idx]).collect();
            Complex64::new(
                reconstruct_centered_coeff_f64(&residues, poly.basis().moduli()),
                0.0,
            )
        })
        .collect();

    let tables = VandermondeTables::<N>::new();
    let scale = 2f64.powi(pt.scale_bits as i32);
    special_dft::<N>(&coeffs, &tables)
        .into_iter()
        .take(pt.slots)
        .map(|slot| (slot / scale).re)
        .collect()
}

fn reconstruct_centered_coeff_f64(residues: &[u64], primes: &[u64]) -> f64 {
    let mut value = BigUint::from(residues[0]);
    let mut prefix_product = BigUint::from(primes[0]);

    for (idx, &prime) in primes.iter().enumerate().skip(1) {
        let current_mod = (&value % prime)
            .to_u64()
            .expect("value mod u64 prime always fits in u64");
        let delta = sub_mod(residues[idx], current_mod, prime);
        if delta != 0 {
            let prefix_mod = (&prefix_product % prime)
                .to_u64()
                .expect("prefix_product mod u64 prime always fits in u64");
            let prefix_inv = mod_inverse(prefix_mod, prime);
            let adjustment =
                ((delta as u128 * prefix_inv as u128) % prime as u128) as u64;
            if adjustment != 0 {
                value += &prefix_product * BigUint::from(adjustment);
            }
        }
        prefix_product *= prime;
    }

    let half_modulus = &prefix_product >> 1usize;
    if value > half_modulus {
        let magnitude = prefix_product - value;
        -magnitude
            .to_f64()
            .expect("centered coefficient magnitude must fit in f64")
    } else {
        value
            .to_f64()
            .expect("centered coefficient must fit in f64")
    }
}

fn sub_mod(value: u64, subtrahend: u64, modulus: u64) -> u64 {
    if value >= subtrahend {
        (value - subtrahend) % modulus
    } else {
        (value + modulus - subtrahend) % modulus
    }
}

fn mod_inverse(value: u64, modulus: u64) -> u64 {
    let (gcd, x, _) = extended_gcd(value as i128, modulus as i128);
    assert_eq!(gcd, 1, "mod_inverse: values must be coprime");
    ((x % modulus as i128 + modulus as i128) % modulus as i128) as u64
}

fn extended_gcd(value: i128, modulus: i128) -> (i128, i128, i128) {
    if value == 0 {
        (modulus, 0, 1)
    } else {
        let (gcd, x, y) = extended_gcd(modulus % value, value);
        (gcd, y - (modulus / value) * x, x)
    }
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
    let ct_sum = sum_slots_power_of_two(&ct_input, SLOTS, &rotkeys);
    println!(
        "After sum_slots: logp={}, logq={}, primes={}",
        ct_sum.logp,
        ct_sum.logq,
        ct_sum.c0.basis().channel_count()
    );
    println!(
        "  skipped high-level sum decode: current plaintext coeff API is i64-bound"
    );
    println!();

    let ct_mean = mul_plain_scalar(&ct_sum, 1.0 / SLOTS as f64, SCALAR_SCALE);
    println!(
        "After mean = sum * 1/{SLOTS}: logp={}, logq={}, primes={}",
        ct_mean.logp,
        ct_mean.logq,
        ct_mean.c0.basis().channel_count()
    );

    println!(
        "  skipped high-level mean decode: current plaintext coeff API is i64-bound"
    );
    println!();

    let ct_input_mean_level =
        reduce_ct_to_basis(&ct_input, ct_mean.c0.basis().clone(), ct_mean.logq);
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

    let sk_centered = reduce_sk(&sk, ct_centered.c0.basis().clone());
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

    let sk_sq = reduce_sk(&sk, ct_sq_diff.c0.basis().clone());
    let engine_sq = make_engine(ct_sq_diff.c0.basis().clone());
    let rotkeys_sq = make_rotkeys(&engine_sq, &sk_sq, &mut rng);
    let ct_var_sum = sum_slots_power_of_two(&ct_sq_diff, SLOTS, &rotkeys_sq);
    let ct_variance =
        mul_plain_scalar(&ct_var_sum, 1.0 / SLOTS as f64, SCALAR_SCALE);
    println!(
        "After variance = sum_sq_diff * 1/{SLOTS}: logp={}, logq={}, primes={}",
        ct_variance.logp,
        ct_variance.logq,
        ct_variance.c0.basis().channel_count()
    );

    let sk_variance = reduce_sk(&sk, ct_variance.c0.basis().clone());
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
