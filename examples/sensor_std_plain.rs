//! Plaintext reference for the sensor std-dev design (no encryption).
//!
//! The encrypted demo (`examples/sensor_std_dev/`) normalizes readings with
//! a range contract and approximates sqrt with a Chebyshev polynomial fitted
//! on the contract's std envelope. This example checks that design in plain
//! f64 using the same library fitter and evaluator.
//!
//! Three sections:
//!   1. Envelope x degree sweep: max std error in degC over the envelope
//!   2. The demo pipeline in plaintext (normalize -> variance -> sqrt ->
//!      restore)
//!   3. Behavior outside the envelope (why the client checks it)

use toy_heaan_ckks::ckks::algorithms::chebyshev::{
    eval_chebyshev_clenshaw, fit_chebyshev, to_chebyshev_domain,
};

// ── Contract (matches `SensorContract::demo()`) ──────────────────────────

const RANGE_LO: f64 = -40.0;
const RANGE_HI: f64 = 60.0;
const STD_MIN: f64 = 1.0;
const STD_MAX: f64 = 8.0;
const DEGREE: usize = 7;

const WIDTH: f64 = RANGE_HI - RANGE_LO;

// ── Sweep parameters ─────────────────────────────────────────────────────

const GRID_POINTS: usize = 2_001;
const ENVELOPES: [(f64, f64); 5] =
    [(1.0, 10.0), (1.0, 8.0), (0.5, 8.0), (2.0, 10.0), (1.0, 6.0)];
const DEGREES: [usize; 4] = [5, 7, 9, 11];

/// 32 readings, one per 45 min: a fixed diurnal day (std ~2.9 degC).
const READINGS: [f64; 32] = [
    11.4, 11.0, 10.9, 11.2, 11.6, 12.4, 13.2, 14.3, 15.2, 16.3, 17.1, 18.0, 18.6,
    19.0, 19.2, 19.1, 18.8, 18.1, 17.5, 16.6, 15.7, 14.8, 13.9, 13.2, 12.6, 12.1,
    11.8, 11.5, 11.3, 11.2, 11.3, 11.4,
];

// ── Helpers ──────────────────────────────────────────────────────────────

/// Population std (divide by n), matching the encrypted pipeline.
fn population_std(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    (values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt()
}

fn normalize(reading: f64) -> f64 {
    (reading.clamp(RANGE_LO, RANGE_HI) - RANGE_LO) / WIDTH
}

/// Normalized variance interval the sqrt is fitted on.
fn variance_interval(std_min: f64, std_max: f64) -> (f64, f64) {
    ((std_min / WIDTH).powi(2), (std_max / WIDTH).powi(2))
}

/// Approximate std in degC from a std in degC, through the normalized
/// variance and the Chebyshev sqrt (what the encrypted pipeline computes).
fn approx_std(std_dev: f64, interval: (f64, f64), coeffs: &[f64]) -> f64 {
    let (lo, hi) = interval;
    let variance_normalized = (std_dev / WIDTH).powi(2);
    let t = to_chebyshev_domain(variance_normalized, lo, hi);
    eval_chebyshev_clenshaw(t, coeffs) * WIDTH
}

/// Max |approx - true| std error over the envelope, and where it occurs.
fn max_envelope_error(std_min: f64, std_max: f64, degree: usize) -> (f64, f64) {
    let interval = variance_interval(std_min, std_max);
    let coeffs = fit_chebyshev(f64::sqrt, interval.0, interval.1, degree);
    (0..GRID_POINTS)
        .map(|i| {
            let std_dev =
                std_min + (std_max - std_min) * i as f64 / (GRID_POINTS - 1) as f64;
            (
                (approx_std(std_dev, interval, &coeffs) - std_dev).abs(),
                std_dev,
            )
        })
        .fold(
            (0.0, std_min),
            |worst, item| {
                if item.0 > worst.0 { item } else { worst }
            },
        )
}

// ── Section 1: envelope x degree sweep ───────────────────────────────────

fn section1() {
    println!("=== Section 1: envelope x degree sweep (W = {WIDTH}) ===\n");
    println!("Max std error (degC) over the envelope. Accuracy depends only on");
    println!("std_max / std_min and the degree, not on the range width.\n");

    print!("{:>14} {:>7}", "envelope", "ratio");
    for degree in DEGREES {
        print!("  {:>10}", format!("deg {degree}"));
    }
    println!();

    for (std_min, std_max) in ENVELOPES {
        let ratio = (std_max / std_min).powi(2);
        print!("{:>14} {ratio:>7.0}", format!("[{std_min}, {std_max}]"));
        for degree in DEGREES {
            let (error, _) = max_envelope_error(std_min, std_max, degree);
            print!("  {error:>10.3e}");
        }
        println!();
    }

    let (error, at) = max_envelope_error(STD_MIN, STD_MAX, DEGREE);
    println!(
        "\nDemo choice: [{STD_MIN}, {STD_MAX}] degC, degree {DEGREE} -> max error \
         {error:.3} degC at std = {at:.2} degC\n"
    );
}

// ── Section 2: demo pipeline in plaintext ────────────────────────────────

fn section2() {
    println!("=== Section 2: demo pipeline in plaintext ===\n");

    let interval = variance_interval(STD_MIN, STD_MAX);
    let coeffs = fit_chebyshev(f64::sqrt, interval.0, interval.1, DEGREE);
    let max_coeff = coeffs.iter().fold(0.0f64, |m, c| m.max(c.abs()));

    let normalized: Vec<f64> = READINGS.iter().map(|&x| normalize(x)).collect();
    let std_normalized = population_std(&normalized);
    let std_true = population_std(&READINGS);

    let t = to_chebyshev_domain(std_normalized.powi(2), interval.0, interval.1);
    let std_approx = eval_chebyshev_clenshaw(t, &coeffs) * WIDTH;

    println!(
        "Contract:     range [{RANGE_LO}, {RANGE_HI}] degC, envelope [{STD_MIN}, {STD_MAX}] degC"
    );
    println!(
        "Fit interval: [{:.2e}, {:.2e}] (normalized variance)",
        interval.0, interval.1
    );
    println!("Max |coeff|:  {max_coeff:.4} (small: no encoding limits in play)");
    println!(
        "Affine slope: {:.1} (4 / interval width)\n",
        4.0 / (interval.1 - interval.0)
    );
    println!("Readings:     32 x 45 min, std {std_true:.6} degC");
    println!("Normalized:   std {std_normalized:.6e} (= std / {WIDTH})");
    println!(
        "Chebyshev:    std {std_approx:.6} degC, error {:.2e} degC\n",
        (std_approx - std_true).abs()
    );
}

// ── Section 3: outside the envelope ──────────────────────────────────────

fn section3() {
    println!("=== Section 3: behavior outside the envelope ===\n");

    let interval = variance_interval(STD_MIN, STD_MAX);
    let coeffs = fit_chebyshev(f64::sqrt, interval.0, interval.1, DEGREE);

    println!(
        "{:>12}  {:>12}  {:>10}  note",
        "true std", "approx std", "error"
    );
    for (std_dev, note) in [
        (0.0, "constant readings"),
        (0.3, "indoor sensor"),
        (STD_MIN, "envelope lower bound"),
        (3.0, "typical outdoor day"),
        (STD_MAX, "envelope upper bound"),
        (10.0, "above envelope"),
        (15.0, "far above envelope"),
    ] {
        let approx = approx_std(std_dev, interval, &coeffs);
        println!(
            "{std_dev:>12.2}  {approx:>12.4}  {:>10.2e}  {note}",
            (approx - std_dev).abs()
        );
    }
    println!(
        "\nBelow the envelope the result is biased high; above it the \
         polynomial diverges.\nThe client checks its plaintext std against \
         the envelope before encrypting.\n"
    );
}

fn main() {
    section1();
    section2();
    section3();
}
