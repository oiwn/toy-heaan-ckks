//! The sensor std-dev design checked in plain f64 (no encryption).
//!
//! ```bash
//! cargo run --example sensor_std_plain
//! ```
//!
//! A narrated design check for `sensor_std_dev`: the encrypted demo
//! normalizes readings with a range contract and approximates sqrt with a
//! Chebyshev polynomial fitted on the contract's std_dev range. This
//! example sweeps ranges and degrees, runs the pipeline in plain f64 with
//! the same library fitter, and shows what happens outside the range.
//!
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
use toy_heaan_ckks::ckks::algorithms::chebyshev::{
    eval_chebyshev_clenshaw, fit_chebyshev, to_chebyshev_domain,
};

#[path = "support/narrate.rs"]
mod narrate;
use narrate::{bad, check, formula, heading, ok, param, say, step, table};

// Contract, matching `SensorContract::demo()`.
const RANGE_LO: f64 = -40.0;
const RANGE_HI: f64 = 60.0;
const STD_MIN: f64 = 1.0;
const STD_MAX: f64 = 8.0;
const DEGREE: usize = 7;
const WIDTH: f64 = RANGE_HI - RANGE_LO;
/// The demo's std_dev budget in degC.
const BOUND: f64 = 0.1;

const GRID_POINTS: usize = 2_001;
const ENVELOPES: [(f64, f64); 5] =
    [(1.0, 10.0), (1.0, 8.0), (0.5, 8.0), (2.0, 10.0), (1.0, 6.0)];
const DEGREES: [usize; 4] = [5, 7, 9, 11];

/// 32 readings, one per 45 min: a fixed diurnal day (std_dev ~2.9 degC).
const READINGS: [f64; 32] = [
    11.4, 11.0, 10.9, 11.2, 11.6, 12.4, 13.2, 14.3, 15.2, 16.3, 17.1, 18.0, 18.6,
    19.0, 19.2, 19.1, 18.8, 18.1, 17.5, 16.6, 15.7, 14.8, 13.9, 13.2, 12.6, 12.1,
    11.8, 11.5, 11.3, 11.2, 11.3, 11.4,
];

/// Population std_dev (divide by n), matching the encrypted pipeline.
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

/// The std_dev the encrypted pipeline would return for a true std_dev, in
/// degC: normalized variance -> Chebyshev sqrt -> back to degC.
fn approx_std(std_dev: f64, interval: (f64, f64), coeffs: &[f64]) -> f64 {
    let (lo, hi) = interval;
    let t = to_chebyshev_domain((std_dev / WIDTH).powi(2), lo, hi);
    eval_chebyshev_clenshaw(t, coeffs) * WIDTH
}

/// Max |approx - true| std_dev error over the range, and where it occurs.
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

fn sweep() -> bool {
    println!("{}", heading("Step 1: which range and degree?"));
    say(
        "Max std_dev error (degC) over each range. sqrt is homogeneous, so\n\
         the error depends only on std_max / std_min and the degree, not on\n\
         the reading range width.\n",
    );
    let rows: Vec<Vec<String>> = ENVELOPES
        .iter()
        .map(|&(std_min, std_max)| {
            let mut row = vec![
                format!("[{std_min}, {std_max}]"),
                format!("{:.0}", (std_max / std_min).powi(2)),
            ];
            for degree in DEGREES {
                row.push(format!(
                    "{:.3e}",
                    max_envelope_error(std_min, std_max, degree).0
                ));
            }
            row
        })
        .collect();
    table(
        &[
            "std_dev range",
            "var ratio",
            "deg 5",
            "deg 7",
            "deg 9",
            "deg 11",
        ],
        &rows,
    );
    let (error, at) = max_envelope_error(STD_MIN, STD_MAX, DEGREE);
    println!();
    step(
        "demo choice",
        format!(
            "[{STD_MIN}, {STD_MAX}] degC, degree {DEGREE}, worst at {at:.2} degC"
        ),
    );
    check("demo range and degree", error, BOUND)
}

fn pipeline() -> bool {
    println!("\n{}", heading("Step 2: the demo pipeline in plain f64"));
    say(&formula(
        "normalize -> std_dev -> variance -> Chebyshev sqrt -> x width",
    ));
    let interval = variance_interval(STD_MIN, STD_MAX);
    let coeffs = fit_chebyshev(f64::sqrt, interval.0, interval.1, DEGREE);
    let max_coeff = coeffs.iter().fold(0.0f64, |m, c| m.max(c.abs()));
    let normalized: Vec<f64> = READINGS.iter().map(|&x| normalize(x)).collect();
    let std_normalized = population_std(&normalized);
    let std_true = population_std(&READINGS);
    let t = to_chebyshev_domain(std_normalized.powi(2), interval.0, interval.1);
    let std_approx = eval_chebyshev_clenshaw(t, &coeffs) * WIDTH;
    step(
        "fit interval (normalized variance)",
        format!("[{:.2e}, {:.2e}]", interval.0, interval.1),
    );
    step(
        "max |coefficient|                 ",
        format!("{max_coeff:.4}"),
    );
    step(
        "affine slope 2 / interval width   ",
        format!("{:.1}", 2.0 / (interval.1 - interval.0)),
    );
    step(
        "readings std_dev                  ",
        format!("{std_true:.6} degC"),
    );
    step(
        "normalized std_dev                ",
        format!("{std_normalized:.6e}"),
    );
    step(
        "Chebyshev std_dev                 ",
        format!("{std_approx:.6} degC"),
    );
    println!();
    check("pipeline std_dev", (std_approx - std_true).abs(), BOUND)
}

fn outside() {
    println!("\n{}", heading("Step 3: outside the range"));
    say("The polynomial is only fitted on [1, 8] degC.\n");
    let interval = variance_interval(STD_MIN, STD_MAX);
    let coeffs = fit_chebyshev(f64::sqrt, interval.0, interval.1, DEGREE);
    let rows: Vec<Vec<String>> = [
        (0.0, "constant readings"),
        (0.3, "indoor sensor"),
        (STD_MIN, "lower bound"),
        (3.0, "typical outdoor day"),
        (STD_MAX, "upper bound"),
        (10.0, "above the range"),
        (15.0, "far above the range"),
    ]
    .iter()
    .map(|&(std_dev, note)| {
        let approx = approx_std(std_dev, interval, &coeffs);
        vec![
            format!("{std_dev:.2}"),
            format!("{approx:.4}"),
            format!("{:.2e}", (approx - std_dev).abs()),
            note.to_string(),
        ]
    })
    .collect();
    table(&["true std_dev", "approx", "error", "note"], &rows);
    say(&format!(
        "\n{} below the range the result is biased high; above it\n\
         the polynomial diverges. The client checks its plaintext std_dev\n\
         against the range before encrypting.",
        bad("Shown on purpose:")
    ));
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("{}\n", heading("Sensor std_dev design, in plain f64"));
    say(
        "sensor_std_dev computes a standard deviation on encrypted data.\n\
         Its only approximate step is the square root, a Chebyshev\n\
         polynomial fitted on a std_dev range the client promises. This\n\
         example checks that design without encryption: same fitter, same\n\
         normalization, exact arithmetic.\n",
    );
    println!("{}", heading("Parameters"));
    param(
        "reading range",
        format!("[{RANGE_LO}, {RANGE_HI}] degC"),
        "normalized to [0, 1]",
    );
    param(
        "std_dev range",
        format!("[{STD_MIN}, {STD_MAX}] degC"),
        "the sqrt fit interval",
    );
    param("degree       ", DEGREE, "depth 3 on ciphertexts");
    param(
        "budget       ",
        format!("{BOUND} degC"),
        "allowed std_dev error",
    );
    println!();

    let mut passed = sweep();
    passed &= pipeline();
    outside();

    println!("\n{}", heading("Takeaways"));
    say(
        "Agreeing on the input range first lets a low-degree polynomial do\n\
         an exact function's job. Next: sensor_std_dev runs this encrypted.",
    );
    if !passed {
        return Err("design check missed the budget".into());
    }
    println!("\n{}", ok("All checks passed."));
    Ok(())
}
