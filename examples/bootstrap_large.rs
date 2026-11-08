//! Bootstrap at large ring degrees (up to N=8192) to study its behaviour.
//!
//! Large degrees need a sparse secret (the mod-raise overflow grows with the
//! secret weight) and multi-level transforms (a one-level transform needs
//! N/2-1 rotation keys). The chain is sized from the schedule. This is not
//! a secure parameter set: N=8192 supports only ~218 modulus bits at
//! 128-bit security, and this chain has over a thousand.
//!
//! Weight, message ratio and double angles default to
//! `BootstrapConfig::reference(N)`, which picks them from N: weight
//! min(N/2, 64), five double angles above weight 32, and a ratio that drops
//! from 2^17 at small N as N grows. Key-switching noise (20-bit gadget
//! digits, no special modulus) grows about N^1.5 in the output, and
//! bootstrap amplifies it by q0/message_scale; a smaller ratio trades that
//! for sine linearization error, which stays small because large-N message
//! coefficients are small. N=8192 needs ~7.5 GB, mostly rotation keys.
//!
//! ```text
//! cargo run --release --example bootstrap_large -- [--degree 16..8192]
//!     [--weight W] [--cts L] [--stc L] [--ratio R] [--poly 31] [--angles A]
//!     [--cycles 2] [--plan-only]
//! ```
use std::sync::Arc;
use std::time::Instant;

use num_complex::Complex64;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::PolyRing;
use toy_heaan_ckks::ckks::algorithms::bootstrap::{BootstrapConfig, BootstrapPlan};
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::ckks::keys::{DIGIT_BITS, RelinearizationKey};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

type Error = Box<dyn std::error::Error>;

struct Options {
    degree: usize,
    weight: Option<usize>,
    cts: Option<usize>,
    stc: Option<usize>,
    ratio: Option<u32>,
    poly: usize,
    angles: Option<u32>,
    cycles: usize,
    plan_only: bool,
}

fn parse() -> Result<Options, Error> {
    let mut options = Options {
        degree: 8192,
        weight: None,
        cts: None,
        stc: None,
        ratio: None,
        poly: 31,
        angles: None,
        cycles: 2,
        plan_only: false,
    };
    let mut args = std::env::args().skip(1);
    while let Some(flag) = args.next() {
        if flag == "--plan-only" {
            options.plan_only = true;
            continue;
        }
        let value = args.next().ok_or(format!("{flag} needs a value"))?;
        match flag.as_str() {
            "--degree" => options.degree = value.parse()?,
            "--weight" => options.weight = Some(value.parse()?),
            "--cts" => options.cts = Some(value.parse()?),
            "--stc" => options.stc = Some(value.parse()?),
            "--ratio" => options.ratio = Some(value.parse()?),
            "--poly" => options.poly = value.parse()?,
            "--angles" => options.angles = Some(value.parse()?),
            "--cycles" => options.cycles = value.parse()?,
            _ => return Err(format!("unknown option {flag}").into()),
        }
    }
    Ok(options)
}

fn mib(bytes: usize) -> f64 {
    bytes as f64 / f64::from(1 << 20)
}

fn run<const N: usize>(options: &Options) -> Result<(), Error> {
    let stages = (N / 2).ilog2() as usize;
    // Defaults: CtS groups of two butterfly stages (6 rotations per level,
    // keys at the large top levels stay few), StC groups of three (14
    // rotations, but at small levels).
    let reference = BootstrapConfig::reference(N)?;
    let config = BootstrapConfig {
        max_secret_weight: options.weight.unwrap_or(reference.max_secret_weight),
        message_ratio_bits: options.ratio.unwrap_or(reference.message_ratio_bits),
        polynomial_degree: options.poly,
        double_angles: options.angles.unwrap_or(reference.double_angles),
        coeffs_to_slots_levels: options.cts.unwrap_or(stages.div_ceil(2)),
        slots_to_coeffs_levels: options.stc.unwrap_or(stages.div_ceil(3)),
        ..reference
    };
    let polynomial_depth = config.polynomial_degree.ilog2() as usize + 1;
    let depth = config.coeffs_to_slots_levels
        + polynomial_depth
        + config.double_angles as usize
        + config.slots_to_coeffs_levels;
    // One spare level after the output, for the squaring between cycles.
    let primes = depth + 2;

    let started = Instant::now();
    let basis =
        Arc::new(RnsBasis::<N>::new(generate_primes(61, primes, N as u64))?);
    let plan = BootstrapPlan::<RnsPoly<N>, N>::new(basis.clone(), config.clone())?;
    let schedule = plan.schedule();
    println!(
        "N={N} slots={} weight={} radius={} ratio=2^{} message_scale=2^{}",
        N / 2,
        config.max_secret_weight,
        config.radius(),
        config.message_ratio_bits,
        schedule.message_scale().log2()
    );
    println!(
        "budgets cts={} stc={} depth={} primes={primes} modulus_bits~{} \
         output_level={} plan={:?}",
        config.coeffs_to_slots_levels,
        config.slots_to_coeffs_levels,
        schedule.depth(),
        61 * primes,
        schedule.output_level(),
        started.elapsed()
    );
    for stage in schedule.stages() {
        println!(
            "  {:<18} level {:>2} -> {:>2}  log2(scale)={:.2}",
            format!("{:?}", stage.stage),
            stage.input_level,
            stage.output_level,
            stage.output_scale.log2()
        );
    }

    // Each key holds two polynomials per gadget term at its level.
    let required = schedule.required_keys();
    let key_bytes = |level: usize| -> Result<usize, Error> {
        let context = basis.at_level(level).ok_or("level outside the chain")?;
        Ok(2 * RnsPoly::<N>::gadget_len(&context, DIGIT_BITS) * N * level * 8)
    };
    let mut total = 0;
    for &level in &required.relinearization {
        total += key_bytes(level)?;
    }
    for &(level, _) in &required.rotations {
        total += key_bytes(level)?;
    }
    for &level in &required.conjugation {
        total += key_bytes(level)?;
    }
    println!(
        "keys: relin={} rotation={} conjugation={} estimated_memory={:.0} MiB",
        required.relinearization.len(),
        required.rotations.len(),
        required.conjugation.len(),
        mib(total)
    );
    if options.plan_only {
        return Ok(());
    }

    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: config.max_secret_weight,
            scale: schedule.message_scale(),
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let secret = engine.generate_secret_key(&mut rng)?;
    let public = engine.generate_public_key(&secret, &mut rng)?;
    let started = Instant::now();
    let keys = plan.generate_keys(&secret, 3.2f64.sqrt(), &mut rng)?;
    println!("key generation {:?}", started.elapsed());

    let encoder = CkksEncoder::new(schedule.message_scale());
    let mut expected: Vec<_> = (0..N / 2)
        .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
        .collect();
    let mut current = engine.encrypt(
        &encoder.encode_complex(&expected, basis),
        &public,
        &mut rng,
    );
    for cycle in 0..options.cycles {
        let started = Instant::now();
        let refreshed = plan.bootstrap(&current, &keys)?;
        let elapsed = started.elapsed();
        let sk = secret.mod_switch_to(refreshed.level())?;
        let actual = CkksEncoder::new(refreshed.scale)
            .decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
        let errors: Vec<_> = actual
            .iter()
            .zip(&expected)
            .map(|(a, b)| (*a - *b).norm())
            .collect();
        let max = errors.iter().copied().fold(0.0, f64::max);
        let mean = errors.iter().sum::<f64>() / errors.len() as f64;
        println!(
            "cycle={cycle} time={elapsed:?} level={} log2(scale)={:.2} \
             max_error={max:.3e} mean_error={mean:.3e} precision_bits={:.1} {}",
            refreshed.level(),
            refreshed.scale.log2(),
            -max.log2(),
            if max < 1e-3 {
                "pass"
            } else {
                "FAIL (gate 1e-3)"
            }
        );
        if cycle + 1 < options.cycles {
            let relin = RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
            current = refreshed.mul(&refreshed, &relin).rescale()?;
            for x in &mut expected {
                *x *= *x;
            }
        }
    }
    Ok(())
}

fn main() -> Result<(), Error> {
    let options = parse()?;
    match options.degree {
        16 => run::<16>(&options),
        32 => run::<32>(&options),
        64 => run::<64>(&options),
        128 => run::<128>(&options),
        256 => run::<256>(&options),
        512 => run::<512>(&options),
        1024 => run::<1024>(&options),
        2048 => run::<2048>(&options),
        4096 => run::<4096>(&options),
        8192 => run::<8192>(&options),
        other => Err(format!("unsupported degree {other}").into()),
    }
}
