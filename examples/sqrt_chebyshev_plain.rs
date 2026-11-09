//! Approximating sqrt(x) with a Chebyshev polynomial, in plain f64.
//!
//! ```bash
//! cargo run --example sqrt_chebyshev_plain
//! ```
//!
//! A narrated math sandbox, no encryption: CKKS can only add and multiply,
//! so a function like sqrt must become a polynomial. This example fits
//! sqrt(x) on an interval at Chebyshev nodes, evaluates the polynomial with
//! Clenshaw's recurrence, and measures how the error falls with the degree.
//! The fitter here is written out in full; the library's
//! `ckks::algorithms::chebyshev` does the same for encrypted inputs.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use std::f64::consts::PI;

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{check, dim, formula, heading, ok, param, say, step, table};

/// sqrt is fitted on [LOWER, UPPER]: positive, like a variance.
const LOWER: f64 = 0.01;
const UPPER: f64 = 0.2;
const DEGREE: usize = 7;
const GRID_POINTS: usize = 20_001;
/// A variance from the std-dev demos, inside the interval.
const CURRENT_VARIANCE: f64 = 0.014_093_75;
/// Max error of the degree-7 fit on the interval.
const BOUND: f64 = 1e-3;

/// Maps x in [lower, upper] to t in [-1, 1].
fn to_chebyshev_domain(x: f64, lower: f64, upper: f64) -> f64 {
    (2.0 * x - (upper + lower)) / (upper - lower)
}

/// Maps t in [-1, 1] back to x in [lower, upper].
fn from_chebyshev_domain(t: f64, lower: f64, upper: f64) -> f64 {
    0.5 * ((upper - lower) * t + upper + lower)
}

/// Coefficients c_k of sum c_k * T_k(t) interpolating sqrt at the
/// degree + 1 Chebyshev nodes t_j = cos(pi * (j + 1/2) / (degree + 1)).
fn fit_sqrt_chebyshev(lower: f64, upper: f64, degree: usize) -> Vec<f64> {
    let node_count = degree + 1;
    let theta = |j: usize| PI * (j as f64 + 0.5) / node_count as f64;
    let values: Vec<f64> = (0..node_count)
        .map(|j| from_chebyshev_domain(theta(j).cos(), lower, upper).sqrt())
        .collect();
    (0..node_count)
        .map(|k| {
            let sum: f64 = values
                .iter()
                .enumerate()
                .map(|(j, value)| value * (k as f64 * theta(j)).cos())
                .sum();
            if k == 0 {
                sum / node_count as f64
            } else {
                2.0 * sum / node_count as f64
            }
        })
        .collect()
}

/// Clenshaw's recurrence: evaluates sum c_k * T_k(t) without forming T_k.
fn eval_chebyshev_clenshaw(t: f64, coeffs: &[f64]) -> f64 {
    let mut b_next = 0.0;
    let mut b_next_next = 0.0;
    for &coeff in coeffs.iter().skip(1).rev() {
        let b = 2.0 * t * b_next - b_next_next + coeff;
        b_next_next = b_next;
        b_next = b;
    }
    t * b_next - b_next_next + coeffs[0]
}

fn approx_sqrt(x: f64, coeffs: &[f64]) -> f64 {
    eval_chebyshev_clenshaw(to_chebyshev_domain(x, LOWER, UPPER), coeffs)
}

/// Largest |approx - sqrt| on a fine grid, and where it occurs.
fn max_grid_error(coeffs: &[f64]) -> (f64, f64) {
    (0..GRID_POINTS)
        .map(|i| {
            let x = LOWER + (UPPER - LOWER) * i as f64 / (GRID_POINTS - 1) as f64;
            ((approx_sqrt(x, coeffs) - x.sqrt()).abs(), x)
        })
        .fold(
            (0.0, LOWER),
            |worst, item| {
                if item.0 > worst.0 { item } else { worst }
            },
        )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("sqrt(x) as a Chebyshev polynomial"));
    say(&format!(
        "A ciphertext supports +, - and *, so every function must be a\n\
         polynomial. A Chebyshev fit on an interval spreads the error evenly\n\
         instead of letting it explode at the ends (as a Taylor series does):\n  \
         {}\n  {}\n\
         Outside the interval the polynomial is meaningless, so the caller\n\
         must know where its inputs lie.\n",
        formula("t = (2x - (hi + lo)) / (hi - lo)       maps x to [-1, 1]"),
        formula(
            "p(x) = sum_k c_k * T_k(t)              T_k: Chebyshev polynomials"
        )
    ));

    println!("{}", heading("Parameters"));
    param(
        "interval",
        format!("[{LOWER}, {UPPER}]"),
        "variances of the std-dev demos",
    );
    param("degree  ", DEGREE, "depth 3 when evaluated on ciphertexts");
    param(
        "grid    ",
        GRID_POINTS,
        "points where the error is measured",
    );
    println!();

    println!("{}", heading("Step 1: error against degree"));
    say("Each extra degree costs multiplications (and levels in CKKS).\n");
    let rows: Vec<Vec<String>> = [3usize, 5, 7, 9, 11]
        .iter()
        .map(|&degree| {
            let (error, x) =
                max_grid_error(&fit_sqrt_chebyshev(LOWER, UPPER, degree));
            vec![
                degree.to_string(),
                format!("{error:.3e}"),
                format!("{x:.5}"),
            ]
        })
        .collect();
    table(&["degree", "max error", "at x"], &rows);
    say(
        "\nThe worst error sits at the lower end: sqrt is steepest near 0,\n\
         so a wider interval toward 0 needs a higher degree.",
    );

    println!(
        "\n{}",
        heading(&format!("Step 2: the degree-{DEGREE} coefficients"))
    );
    let coeffs = fit_sqrt_chebyshev(LOWER, UPPER, DEGREE);
    for (k, c) in coeffs.iter().enumerate() {
        step(&format!("c_{k}"), format!("{c:>12.8}"));
    }
    println!(
        "  {}",
        dim("They shrink fast: sqrt is smooth on this interval.")
    );

    println!("\n{}", heading("Step 3: samples"));
    let samples = [LOWER, CURRENT_VARIANCE, 0.02, 0.05, 0.10, 0.16, UPPER];
    let rows: Vec<Vec<String>> = samples
        .iter()
        .map(|&x| {
            let got = approx_sqrt(x, &coeffs);
            vec![
                format!("{x:.6}"),
                format!("{:.8}", x.sqrt()),
                format!("{got:.8}"),
                format!("{:.2e}", (got - x.sqrt()).abs()),
            ]
        })
        .collect();
    table(&["x", "sqrt(x)", "p(x)", "error"], &rows);
    let (error, at) = max_grid_error(&coeffs);
    println!();
    step("worst x", format!("{at:.5}"));
    let passed = check(&format!("degree {DEGREE} on the interval"), error, BOUND);

    println!("\n{}", heading("Takeaways"));
    say(
        "Pick the interval from what the data can be, then the lowest degree\n\
         that meets the error budget. sensor_std_plain applies this to the\n\
         std-dev pipeline; sensor_std_dev runs it encrypted.",
    );
    if !passed {
        return Err("sqrt approximation above bound".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
