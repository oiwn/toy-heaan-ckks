//! Bootstrap at large ring degrees (up to N=8192) to study its behaviour.
//!
//! ```bash
//! cargo run --release --example bootstrap_large -- [--degree 16..8192]
//!     [--weight W] [--cts L] [--stc L] [--ratio R] [--poly 31] [--angles A]
//!     [--cycles 2] [--plan-only]
//! ```
//!
//! A narrated experiment: derive the bootstrap schedule for a large ring,
//! estimate the key memory, then run refresh cycles and report the error
//! and precision. Large degrees need a sparse secret (the ModRaise overflow
//! grows with the weight) and multi-level transforms (a one-level transform
//! needs N/2 - 1 rotation keys). Not a secure parameter set: N=8192 allows
//! about 218 modulus bits at 128-bit security; this chain has over 1000.
//! N=8192 needs about 7.5 GB, mostly rotation keys; start with
//! `--degree 256` or `--plan-only`.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
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

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, dim, elapsed, formula, heading, ok, param, say, step, table};

type Error = Box<dyn std::error::Error>;

const BOUND: f64 = 1e-3;
const SEED: u64 = 0;

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

fn intro() {
    println!("{}\n", heading("Bootstrapping at large ring degrees"));
    say(
        "bootstrap_demo refreshes N=16/64. Real parameters use N in the\n\
         thousands, and three things change with N:",
    );
    say(&format!(
        "  {}  the slot transforms are split into several levels,\n\
         \x20             each with a few rotations, instead of N/2 - 1 keys\n  \
         {}  the secret is sparse (at most 64 nonzeros), so the\n\
         \x20             ModRaise overflow I stays small\n  \
         {}  key-switch noise grows about N^1.5 and the refresh\n\
         \x20             amplifies it by q0 / message scale, so that ratio\n\
         \x20             shrinks as N grows",
        formula("transforms"),
        formula("secret    "),
        formula("noise     ")
    ));
    say("\nDefaults come from BootstrapConfig::reference(N); flags override.\n");
}

fn run<const N: usize>(options: &Options) -> Result<bool, Error> {
    let stages = (N / 2).ilog2() as usize;
    // CtS groups of two butterfly stages (few keys at the large top
    // levels), StC groups of three (more rotations, at small levels).
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

    println!("{}", heading("Parameters"));
    param("N               ", N, &format!("{} complex slots", N / 2));
    param(
        "secret weight   ",
        config.max_secret_weight,
        &format!("overflow bound |I| <= {}", config.interval_bound()),
    );
    param(
        "message ratio   ",
        format!("2^{}", config.message_ratio_bits),
        &format!("message scale 2^{:.0}", schedule.message_scale().log2()),
    );
    param(
        "transform levels",
        format!(
            "{} + {}",
            config.coeffs_to_slots_levels, config.slots_to_coeffs_levels
        ),
        "CoeffsToSlots + SlotsToCoeffs",
    );
    param(
        "primes          ",
        format!("{primes} x 61 bits"),
        &format!("depth {} + 2 spare", schedule.depth()),
    );
    param(
        "seed            ",
        SEED,
        "fixed, so every run prints the same numbers",
    );
    println!(
        "  {}\n",
        dim(&format!("plan derived in {}", elapsed(started.elapsed())))
    );

    println!("{}", heading("Schedule"));
    let rows: Vec<Vec<String>> = schedule
        .stages()
        .iter()
        .map(|stage| {
            vec![
                format!("{:?}", stage.stage),
                format!("{} -> {}", stage.input_level, stage.output_level),
                format!("2^{:.2}", stage.output_scale.log2()),
            ]
        })
        .collect();
    table(&["stage", "levels", "output scale"], &rows);

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
    println!("\n{}", heading("Keys"));
    step("relinearization keys", required.relinearization.len());
    step("rotation keys       ", required.rotations.len());
    step("conjugation keys    ", required.conjugation.len());
    step("estimated memory    ", format!("{:.0} MiB", mib(total)));
    if options.plan_only {
        println!("\n{}", dim("--plan-only: stopping before key generation."));
        return Ok(true);
    }

    let engine = CkksEngine::<RnsPoly<N>, N>::new(
        basis.clone(),
        CkksParams {
            error_variance: 3.2,
            hamming_weight: config.max_secret_weight,
            scale: schedule.message_scale(),
        },
    );
    let mut rng = ChaCha20Rng::seed_from_u64(SEED);
    let secret = engine.generate_secret_key(&mut rng)?;
    let public = engine.generate_public_key(&secret, &mut rng)?;
    let started = Instant::now();
    let keys = plan.generate_keys(&secret, 3.2f64.sqrt(), &mut rng)?;
    println!(
        "  {}",
        dim(&format!("key generation {}", elapsed(started.elapsed())))
    );

    println!("\n{}", heading("Refresh cycles"));
    say("Bootstrap, then square (one level) and bootstrap again.\n");
    let encoder = CkksEncoder::new(schedule.message_scale());
    let mut expected: Vec<_> = (0..N / 2)
        .map(|j| Complex64::from_polar(0.8, j as f64 * 0.7))
        .collect();
    let mut current = engine.encrypt(
        &encoder.encode_complex(&expected, basis),
        &public,
        &mut rng,
    );
    let mut passed = true;
    for cycle in 0..options.cycles {
        let started = Instant::now();
        let refreshed = plan.bootstrap(&current, &keys)?;
        let time = elapsed(started.elapsed());
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
        step(
            &format!("cycle {cycle}"),
            format!(
                "{time}, level {}, scale 2^{:.2}, mean error {mean:.3e}, \
                 {:.1} bits",
                refreshed.level(),
                refreshed.scale.log2(),
                -max.log2()
            ),
        );
        passed &= check(&format!("  cycle {cycle} slots"), max, BOUND);
        if cycle + 1 < options.cycles {
            let relin = RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
            current = refreshed.mul(&refreshed, &relin).rescale()?;
            for x in &mut expected {
                *x *= *x;
            }
        }
    }
    Ok(passed)
}

fn main() -> Result<(), Error> {
    let options = parse()?;
    intro();
    let passed = match options.degree {
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
    }?;

    println!("\n{}", heading("Takeaways"));
    say(
        "Precision falls as N grows: key-switch noise grows with N and the\n\
         refresh amplifies it. Hybrid key switching and a special modulus\n\
         (OpenFHE, Lattigo) are the usual fix; see\n\
         specs/research/large-degree-bootstrap.md.",
    );
    if !passed {
        return Err(format!("a refresh missed the {BOUND:.0e} gate").into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
