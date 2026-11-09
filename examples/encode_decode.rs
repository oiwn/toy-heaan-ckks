//! Encoding and decoding: numbers in slots, a polynomial in the ring.
//!
//! ```bash
//! cargo run --example encode_decode
//! ```
//!
//! A narrated walk-through: why CKKS encodes vectors as polynomials, how the
//! scale turns real numbers into integers, a round trip of real and complex
//! slots, and how the scale sets the precision. No encryption yet.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::f64::consts::{E, PI, SQRT_2};
use std::sync::Arc;

use rustfft::num_complex::Complex64;
use toy_heaan_ckks::PolyRing;
use toy_heaan_ckks::ckks::Scale;
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::types::Plaintext;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, dim, formula, heading, ok, param, say, step, value};

const N: usize = 16;
const SCALE_BITS: i32 = 30;

type Pt = Plaintext<RnsPoly<N>, N>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Encoding: from slots to a polynomial"));
    say(
        "CKKS computes on vectors of real or complex numbers, but encrypts\n\
         polynomials with integer coefficients in Z_Q[X]/(X^N + 1). The\n\
         encoder bridges the two. A vector of N/2 complex slots is the list\n\
         of values of a polynomial m(X) at N/2 complex roots of X^N + 1, so\n\
         an inverse FFT (the special IDFT) finds the coefficients.",
    );
    say(&format!(
        "\nThose coefficients are real, not integers. Multiplying by the\n\
         scale and rounding makes them integers:\n  {}\n  {}\n\
         Rounding loses at most 1/2 per coefficient, so the slot error is\n\
         about 1/scale. Adding and multiplying polynomials then adds and\n\
         multiplies the slots: one operation acts on every slot at once.\n",
        formula("encode: m = round(scale * IDFT(z))"),
        formula("decode: z ~ DFT(m) / scale")
    ));

    println!("{}", heading("Parameters"));
    let scale = Scale::power_of_two(SCALE_BITS)?;
    let primes = generate_primes(31, 3, N as u64);
    param("N     ", N, "ring degree; N/2 = 8 slots");
    param(
        "scale ",
        format!("2^{SCALE_BITS}"),
        "precision about 2^-30 = 9e-10",
    );
    param(
        "primes",
        "3 x 31 bits",
        "Q ~ 2^93 holds the scaled coefficients",
    );
    println!(
        "  {}\n",
        dim("each prime is 1 mod 2N, so the ring has the roots NTT needs")
    );
    let basis = Arc::new(RnsBasis::<N>::new(primes)?);
    let encoder = CkksEncoder::<N>::new(scale);
    let bound = 2.0 / scale.value();
    let mut passed = true;

    println!("{}", heading("Step 1: real slots"));
    let originals = [PI, E, SQRT_2, 1.0 / 3.0, -1.5, 0.5, -0.25, 0.125];
    let pt: Pt = encoder.encode(&originals, basis.clone());
    let coefficients = pt.poly.to_centered_f64();
    say("The encoded polynomial: integers of about scale * |z| = 2^31.");
    step(
        "m[0..4]",
        format!(
            "{:.0}, {:.0}, {:.0}, {:.0}, ...",
            coefficients[0], coefficients[1], coefficients[2], coefficients[3]
        ),
    );
    let decoded = encoder.decode(&pt);
    println!(
        "\n  {:>4} | {:>12} | {:>12} | {:>9}",
        "slot", "input", "decoded", "error"
    );
    println!("  {}", "-".repeat(46));
    let mut worst: f64 = 0.0;
    for (j, (a, b)) in originals.iter().zip(&decoded).enumerate() {
        let error = (a - b).abs();
        worst = worst.max(error);
        println!(
            "  {j:>4} | {a:>12.8} | {} | {error:>9.2e}",
            value(format!("{b:>12.8}"))
        );
    }
    passed &= check("real round trip", worst, bound);

    println!("\n{}", heading("Step 2: complex slots"));
    say(
        "The same encoder packs complex numbers. Unused slots are zero:\n\
         here 4 of the 8 slots are filled.",
    );
    let complex_in = [
        Complex64::new(1.0, 0.5),
        Complex64::new(-0.5, 0.25),
        Complex64::new(0.0, -1.0),
        Complex64::new(0.75, -0.75),
    ];
    let pt: Pt = encoder.encode_complex(&complex_in, basis.clone());
    let complex_out = encoder.decode_complex(&pt);
    let mut worst: f64 = 0.0;
    for (j, (a, b)) in complex_in.iter().zip(&complex_out).enumerate() {
        worst = worst.max((a - b).norm());
        println!(
            "  {} {}",
            formula(&format!("slot {j}: {:>12} ->", format!("{a:.2}"))),
            value(format!("{b:.8}"))
        );
    }
    passed &= check("complex round trip", worst, bound);

    println!("\n{}", heading("Step 3: the scale sets the precision"));
    say("Encode pi with growing scales: each 10 bits buy 3 decimal digits.\n");
    println!("  {:>10} | {:>14} | {:>9}", "scale", "decoded pi", "error");
    println!("  {}", "-".repeat(40));
    for bits in [10, 20, 30, 40] {
        let encoder = CkksEncoder::<N>::new(Scale::power_of_two(bits)?);
        let pt: Pt = encoder.encode(&[PI], basis.clone());
        let decoded = encoder.decode(&pt)[0];
        println!(
            "  {:>10} | {} | {:>9.2e}",
            format!("2^{bits}"),
            value(format!("{decoded:>14.10}")),
            (PI - decoded).abs()
        );
    }
    println!(
        "  {}",
        dim("At 2^40 the f64 FFT (about 1e-12 here), not the scale, limits it.")
    );

    println!("\n{}", heading("Takeaways"));
    say(
        "Slots are the user's view, coefficients the ring's view, and the\n\
         scale is the fixed-point factor between them. A larger scale gives\n\
         more precision but needs a larger modulus. Next: keys, then\n\
         encrypt_add.",
    );
    println!(
        "{}",
        dim("Reference: fhetextbook.github.io/EncodingandDecoding.html")
    );
    if !passed {
        return Err("encoding error above 2/scale".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
