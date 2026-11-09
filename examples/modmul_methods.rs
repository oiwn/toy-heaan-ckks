//! Three ways to compute `a * b mod q`: u128 division, Barrett, Shoup.
//!
//! ```bash
//! cargo run --release --example modmul_methods
//! ```
//!
//! A narrated walk-through: the program prints why each method works, shows
//! one 61-bit product step by step, checks that all three agree on random
//! inputs and times them. The library versions live in
//! `toy_heaan_ckks::math::modular`; precise timings are in
//! `cargo bench --bench arith` and `specs/benchmarks.md`.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.

use std::hint::black_box;
use std::time::Instant;

use owo_colors::{OwoColorize, Stream};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::math::modular::{Barrett, mul_shoup, shoup_companion};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{dim, formula, heading, ok, say, step, subheading, value};

/// Products per timing loop.
const LOOP_LEN: usize = 1 << 20;
/// Timing repetitions; the fastest is reported.
const REPEATS: usize = 5;
/// Random pairs per modulus in the agreement check.
const CHECKS: usize = 100_000;

/// The reference: widen to 128 bits and divide.
fn mul_u128(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(modulus)) as u64
}

/// High 128 bits of the 256-bit product `x * y`, from 64-bit words.
fn mul_high_u128(x: u128, y: u128) -> u128 {
    let (x_hi, x_lo) = (x >> 64, x & u128::from(u64::MAX));
    let (y_hi, y_lo) = (y >> 64, y & u128::from(u64::MAX));
    let low_low = (x_lo * y_lo) >> 64;
    let low_high = x_lo * y_hi;
    let high_low = x_hi * y_lo;
    let middle = low_low
        + (low_high & u128::from(u64::MAX))
        + (high_low & u128::from(u64::MAX));
    x_hi * y_hi + (low_high >> 64) + (high_low >> 64) + (middle >> 64)
}

fn main() {
    let modulus = generate_primes(61, 1, 1024)[0];
    intro();
    worked_example(modulus);
    summary();
    agreement();
    timing();
}

/// The problem all three methods solve.
fn intro() {
    println!(
        "{}\n",
        heading("Modular multiplication: u128, Barrett, Shoup")
    );
    say(
        "Every CKKS operation is built from products of residues modulo NTT\n\
         primes q of up to 61 bits. The product a*b then needs up to 122 bits,\n\
         and reducing it with % calls a software 128-by-64 division routine:\n\
         Rust exposes no CPU instruction for it, and the routine is slow.",
    );
    say(
        "\nThe trick of both fast methods: precompute an approximation of 1/q,\n\
         multiply by it to estimate the quotient floor(a*b / q), subtract\n\
         quotient*q, and fix the small error with at most one subtraction.\n\
         No division happens per product.\n",
    );
}

/// One product, each method narrated and spelled out.
fn worked_example(modulus: u64) {
    let a = modulus - 12_345_678_901;
    let b = modulus / 3 + 7;
    let product = u128::from(a) * u128::from(b);

    println!("{}", heading("One product modulo a 61-bit NTT prime"));
    step("q    ", modulus);
    step("a    ", a);
    step("b    ", b);
    step("a * b", product);
    println!(
        "  {}",
        dim(&format!(
            "(a * b needs {} bits)",
            128 - product.leading_zeros()
        ))
    );

    // u128: one 128-by-64 division in software.
    let reference = mul_u128(a, b, modulus);
    println!("\n{}", subheading("u128: the reference"));
    say(
        "Widen to 128 bits and take %. Simple and obviously correct, so it is\n\
         the oracle the other two are checked against, and the slowest.",
    );
    step("a * b mod q", reference);

    barrett_part(modulus, a, b, product, reference);
    shoup_part(modulus, a, b, reference);

    println!(
        "\n{} {}",
        ok("All three give"),
        value(reference).if_supports_color(Stream::Stdout, |t| t.green())
    );
    say(
        "\nLimit: both fast methods need q < 2^62, so the intermediate values\n\
         fit their words. That is why RnsBasis rejects larger primes.\n",
    );
}

/// Barrett: narration, a decimal toy, then the 61-bit product.
fn barrett_part(modulus: u64, a: u64, b: u64, product: u128, reference: u64) {
    println!("\n{}", subheading("Barrett: one constant per modulus"));
    say(&format!(
        "Precompute once per prime {}, i.e. 1/q scaled\n\
         by 2^128. Then {}.\n\
         The estimate is never too big and at most 1 too small, so\n\
         a*b - quotient*q lies in [0, 2q) and one subtraction fixes it.\n\
         It works for any a and b, so the library uses it for pointwise\n\
         ring products and pow.",
        formula("ratio = floor(2^128 / q)"),
        formula("quotient = floor(a*b * ratio / 2^128)")
    ));

    // Decimal toy: q = 97, 10^4 stands in for 2^128.
    let toy_ratio = 10_000 / 97;
    let toy_quotient = 3500 * toy_ratio / 10_000;
    println!(
        "{}",
        dim("\nToy in decimal: q = 97, 10^4 in place of 2^128, a*b = 3500.")
    );
    step("ratio = floor(10^4 / 97)        ", toy_ratio);
    step("quotient = floor(3500*103/10^4) ", toy_quotient);
    step("3500 - 36*97                    ", 3500 - toy_quotient * 97);
    assert_eq!(3500 - toy_quotient * 97, 3500 % 97);

    let ratio = u128::MAX / u128::from(modulus); // = floor(2^128 / q) for odd q
    let quotient = mul_high_u128(product, ratio) as u64;
    let remainder = (product - u128::from(quotient) * u128::from(modulus)) as u64;
    let barrett = if remainder >= modulus {
        remainder - modulus
    } else {
        remainder
    };
    println!("{}", dim("\nThe 61-bit product:"));
    step("ratio = floor(2^128 / q)             ", ratio);
    step("quotient = floor(a*b * ratio / 2^128)", quotient);
    step(
        "exact floor(a*b / q)                 ",
        product / u128::from(modulus),
    );
    step("a*b - quotient*q   (in [0, 2q))      ", remainder);
    step("after one conditional subtract       ", barrett);
    assert_eq!(barrett, reference);
    assert_eq!(Barrett::new(modulus).mul(a, b), reference);
}

/// Shoup: narration, a decimal toy, then the 61-bit product.
fn shoup_part(modulus: u64, a: u64, b: u64, reference: u64) {
    println!("\n{}", subheading("Shoup: one companion per fixed operand"));
    say(&format!(
        "When one factor w is known in advance, fold w/q into one word:\n\
         {}.\n\
         Then {} and\n\
         r = a*w - estimate*q. The true r is below 2q < 2^64, so plain\n\
         wrapping 64-bit arithmetic computes it exactly; one subtraction\n\
         finishes. The library uses it for NTT twiddles, twist factors and\n\
         rescale constants: each is fixed and gets its companion once.",
        formula("companion = floor(w * 2^64 / q)"),
        formula("estimate = high 64 bits of a*companion")
    ));

    // Decimal toy: q = 97, w = 30, 1000 stands in for 2^64.
    let toy_companion = 30 * 1000 / 97;
    let toy_estimate = 50 * toy_companion / 1000;
    println!(
        "{}",
        dim("\nToy in decimal: q = 97, w = 30, 1000 in place of 2^64, a = 50.")
    );
    step("companion = floor(30*1000 / 97) ", toy_companion);
    step("estimate = floor(50*309 / 1000) ", toy_estimate);
    step(
        "50*30 - 15*97                   ",
        50 * 30 - toy_estimate * 97,
    );
    assert_eq!(50 * 30 - toy_estimate * 97, 50 * 30 % 97);

    let companion = shoup_companion(b, modulus);
    let estimate = ((u128::from(a) * u128::from(companion)) >> 64) as u64;
    let shoup_raw = a
        .wrapping_mul(b)
        .wrapping_sub(estimate.wrapping_mul(modulus));
    let shoup = if shoup_raw >= modulus {
        shoup_raw - modulus
    } else {
        shoup_raw
    };
    println!("{}", dim("\nThe 61-bit product, with w = b:"));
    step("companion = floor(b * 2^64 / q)      ", companion);
    step("estimate = high 64 bits of a*comp.   ", estimate);
    step("a*b - estimate*q mod 2^64 (in [0,2q))", shoup_raw);
    step("after one conditional subtract       ", shoup);
    assert_eq!(shoup, reference);
    assert_eq!(mul_shoup(a, b, companion, modulus), reference);

    say(
        "\nWhy Shoup cannot replace Barrett: when both factors change every\n\
         time, building the companion costs a 128-bit division, exactly\n\
         what we are avoiding. So the library keeps both.",
    );
}

/// When each method applies.
fn summary() {
    let rows = [
        ["method", "precomputed", "applies to", "per product"],
        ["u128", "nothing", "anything", "1 mul + software div"],
        [
            "Barrett",
            "floor(2^128/q) per prime",
            "any a, b",
            "6 multiplies",
        ],
        [
            "Shoup",
            "floor(w*2^64/q) per w",
            "a * fixed w",
            "3 multiplies",
        ],
    ];
    println!("{}", heading("Which one when"));
    for (index, [method, precomputed, applies, cost]) in rows.iter().enumerate() {
        println!("  {method:<8} | {precomputed:<26} | {applies:<22} | {cost}");
        if index == 0 {
            println!("  {}", "-".repeat(78));
        }
    }
    println!(
        "  {}\n",
        dim("Barrett: a*b, 4 for the 256-bit estimate, quotient*q.")
    );
}

/// Random residues: all three methods must give the same result.
fn agreement() {
    println!(
        "{}",
        heading(&format!("Agreement on {CHECKS} random pairs per modulus"))
    );
    let mut rng = ChaCha20Rng::seed_from_u64(7);
    for bits in [30, 50, 61] {
        let modulus = generate_primes(bits, 1, 1024)[0];
        let barrett = Barrett::new(modulus);
        for _ in 0..CHECKS {
            let a = rng.random::<u64>() % modulus;
            let b = rng.random::<u64>() % modulus;
            let expected = mul_u128(a, b, modulus);
            assert_eq!(barrett.mul(a, b), expected);
            let companion = shoup_companion(b, modulus);
            assert_eq!(mul_shoup(a, b, companion, modulus), expected);
        }
        println!(
            "  {bits}-bit q = {}: {}",
            value(modulus),
            ok("u128, Barrett and Shoup agree")
        );
    }
    println!();
}

/// Fastest of `REPEATS` runs of `body`, in nanoseconds per product.
fn time_per_product(mut body: impl FnMut()) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..REPEATS {
        let started = Instant::now();
        body();
        best = best.min(started.elapsed().as_nanos() as f64);
    }
    best / LOOP_LEN as f64
}

/// Loops of pointwise products, as in a ring multiplication or an NTT.
fn timing() {
    println!(
        "{}",
        heading(&format!(
            "Rough timing: {LOOP_LEN} products per loop, best of {REPEATS}"
        ))
    );
    println!(
        "  {:>6} | {:>10} | {:>10} | {:>10} | {:>8} | {:>8}",
        "q bits", "u128", "Barrett", "Shoup", "Barrett", "Shoup"
    );
    println!(
        "  {:>6} | {:>10} | {:>10} | {:>10} | {:>8} | {:>8}",
        "", "", "", "", "speedup", "speedup"
    );
    let mut rng = ChaCha20Rng::seed_from_u64(11);
    for bits in [30, 50, 61] {
        let modulus = generate_primes(bits, 1, 1024)[0];
        let barrett = Barrett::new(modulus);
        let left: Vec<u64> = (0..LOOP_LEN)
            .map(|_| rng.random::<u64>() % modulus)
            .collect();
        let right: Vec<u64> = (0..LOOP_LEN)
            .map(|_| rng.random::<u64>() % modulus)
            .collect();
        // Shoup needs a fixed operand, like an NTT twiddle.
        let fixed = right[0];
        let companion = shoup_companion(fixed, modulus);
        let mut out = vec![0u64; LOOP_LEN];

        let u128_ns = time_per_product(|| {
            for ((o, &a), &b) in out.iter_mut().zip(&left).zip(&right) {
                *o = mul_u128(a, b, modulus);
            }
            black_box(&out);
        });
        let barrett_ns = time_per_product(|| {
            for ((o, &a), &b) in out.iter_mut().zip(&left).zip(&right) {
                *o = barrett.mul(a, b);
            }
            black_box(&out);
        });
        let shoup_ns = time_per_product(|| {
            for (o, &a) in out.iter_mut().zip(&left) {
                *o = mul_shoup(a, fixed, companion, modulus);
            }
            black_box(&out);
        });

        // Pad before coloring, so escape codes do not break the columns.
        let barrett_speedup = format!("{:>7.1}x", u128_ns / barrett_ns);
        let shoup_speedup = format!("{:>7.1}x", u128_ns / shoup_ns);
        println!(
            "  {bits:>6} | {u128_ns:>7.2} ns | {barrett_ns:>7.2} ns | \
             {shoup_ns:>7.2} ns | {} | {}",
            ok(&barrett_speedup),
            ok(&shoup_speedup)
        );
    }
    println!(
        "\n{}",
        dim(
            "On 30-bit primes the 128-bit product has a small high word, so the\n\
             division routine takes a fast path and u128 looks better than it\n\
             is for the 50-61 bit primes CKKS actually uses. Run in --release."
        )
    );
}
