//! Plaintext sensor std-dev pipeline with Chebyshev sqrt approximation.
//!
//! Answers: can Chebyshev degree-7 sqrt on [LOWER, 1.0] serve a sensor
//! std-dev pipeline, and what normalization strategy works?
//!
//! Three sections:
//!   1. Chebyshev fit quality sweep (LOWER x DEGREE, probe at key x-values)
//!   2. Full sensor std-dev pipeline (multiple normalization bounds)
//!   3. Edge cases (zero variance, extreme variance, tiny jitter)

// ── Sensor parameters ────────────────────────────────────────────────────

const SENSOR_LO: f64 = -50.0;
const SENSOR_HI: f64 = 100.0;
const SENSOR_COUNT: usize = 16;

const SENSOR_DATA: [f64; SENSOR_COUNT] = [
    14.0, 14.0, 15.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 20.0, 20.0, 17.0,
    16.0, 15.0, 14.0,
];

/// Theoretical maximum variance: (hi - lo)^2 / 4  (population / N denominator).
const MAX_VAR_POP: f64 = {
    let range = SENSOR_HI - SENSOR_LO;
    range * range / 4.0
};

/// Theoretical maximum sample variance: N/(N-1) × population max.
const MAX_VAR_SAMPLE: f64 = {
    let range = SENSOR_HI - SENSOR_LO;
    let n = SENSOR_COUNT as f64;
    (n / (n - 1.0)) * (range * range / 4.0)
};

// ── Chebyshev parameters ─────────────────────────────────────────────────

const GRID_POINTS: usize = 20_001;
const SWEEP_LOWERS: [f64; 7] = [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1];
const SWEEP_DEGREES: [usize; 5] = [3, 5, 7, 9, 11];
const PROBE_X: [f64; 4] = [0.001, 0.01, 0.1, 0.5];

/// Normalization bounds to test in Section 2.
const NORM_BOUNDS: [f64; 4] = [10.0, 50.0, 100.0, MAX_VAR_SAMPLE];

/// Default Chebyshev interval for the pipeline.
const CHEBY_LO_DEFAULT: f64 = 0.01;
const CHEBY_HI_DEFAULT: f64 = 1.0;
const CHEBY_DEGREE_DEFAULT: usize = 7;

const ERROR_THRESHOLD: f64 = 1e-3;

// ── Statistics ────────────────────────────────────────────────────────────

fn sample_mean(data: &[f64]) -> f64 {
    data.iter().sum::<f64>() / data.len() as f64
}

fn sample_variance(data: &[f64]) -> f64 {
    let n = data.len() as f64;
    let mean = sample_mean(data);
    let sum_sq = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>();
    sum_sq / (n - 1.0)
}

fn sample_std_dev(data: &[f64]) -> f64 {
    sample_variance(data).sqrt()
}

// ── Chebyshev math ───────────────────────────────────────────────────────

fn to_chebyshev_domain(x: f64, lower: f64, upper: f64) -> f64 {
    (2.0 * x - (upper + lower)) / (upper - lower)
}

fn from_chebyshev_domain(t: f64, lower: f64, upper: f64) -> f64 {
    0.5 * ((upper - lower) * t + upper + lower)
}

fn fit_sqrt_chebyshev(lower: f64, upper: f64, degree: usize) -> Vec<f64> {
    let node_count = degree + 1;
    let mut values = Vec::with_capacity(node_count);

    for j in 0..node_count {
        let theta = std::f64::consts::PI * (j as f64 + 0.5) / node_count as f64;
        let t = theta.cos();
        let x = from_chebyshev_domain(t, lower, upper);
        values.push(x.sqrt());
    }

    let mut coeffs = vec![0.0; node_count];
    for (k, coeff) in coeffs.iter_mut().enumerate() {
        let sum = values
            .iter()
            .enumerate()
            .map(|(j, value)| {
                let theta =
                    std::f64::consts::PI * (j as f64 + 0.5) / node_count as f64;
                value * (k as f64 * theta).cos()
            })
            .sum::<f64>();

        *coeff = if k == 0 {
            sum / node_count as f64
        } else {
            2.0 * sum / node_count as f64
        };
    }

    coeffs
}

fn eval_chebyshev_clenshaw(t: f64, coeffs: &[f64]) -> f64 {
    assert!(!coeffs.is_empty(), "coefficients must not be empty");

    let mut b_next = 0.0;
    let mut b_next_next = 0.0;

    for &coeff in coeffs.iter().skip(1).rev() {
        let b = 2.0 * t * b_next - b_next_next + coeff;
        b_next_next = b_next;
        b_next = b;
    }

    t * b_next - b_next_next + coeffs[0]
}

fn approx_sqrt(x: f64, lower: f64, upper: f64, coeffs: &[f64]) -> f64 {
    let t = to_chebyshev_domain(x, lower, upper);
    eval_chebyshev_clenshaw(t, coeffs)
}

fn max_grid_error(lower: f64, upper: f64, coeffs: &[f64]) -> (f64, f64) {
    let mut max_error = 0.0;
    let mut max_x = lower;

    for i in 0..GRID_POINTS {
        let ratio = i as f64 / (GRID_POINTS - 1) as f64;
        let x = lower + (upper - lower) * ratio;
        let got = approx_sqrt(x, lower, upper, coeffs);
        let error = (got - x.sqrt()).abs();

        if error > max_error {
            max_error = error;
            max_x = x;
        }
    }

    (max_error, max_x)
}

fn point_error(x: f64, lower: f64, upper: f64, coeffs: &[f64]) -> f64 {
    (approx_sqrt(x, lower, upper, coeffs) - x.sqrt()).abs()
}

// ── Pipeline helpers ─────────────────────────────────────────────────────

struct SqrtApproximator {
    lower: f64,
    upper: f64,
    coeffs: Vec<f64>,
}

impl SqrtApproximator {
    fn new(lower: f64, upper: f64, degree: usize) -> Self {
        let coeffs = fit_sqrt_chebyshev(lower, upper, degree);
        Self {
            lower,
            upper,
            coeffs,
        }
    }

    fn eval(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let clipped = if x < self.lower { self.lower } else { x };
        approx_sqrt(clipped, self.lower, self.upper, &self.coeffs)
    }
}

fn pipeline_std(
    data: &[f64],
    norm_bound: f64,
    approx: &SqrtApproximator,
) -> (f64, f64) {
    let var = sample_variance(data);
    let var_norm = var / norm_bound;
    let sqrt_norm = approx.eval(var_norm);
    let std_approx = sqrt_norm * norm_bound.sqrt();
    let std_true = sample_std_dev(data);
    (std_approx, std_true)
}

// ── Section 1: Chebyshev fit quality sweep ───────────────────────────────

fn section1() {
    println!("=== Section 1: Chebyshev fit quality sweep ===\n");
    println!("Sweep: LOWER x DEGREE, max error + probe points\n");

    let degrees = &SWEEP_DEGREES;
    let lowers = &SWEEP_LOWERS;
    let probes = &PROBE_X;

    for &lower in lowers.iter() {
        println!("LOWER = {lower:.0e}");

        // header
        print!("{:>7} ", "degree");
        print!("{:>12} ", "max_err");
        for &px in probes.iter() {
            print!("err@{px:.3}  ");
        }
        println!();

        for &deg in degrees.iter() {
            if lower <= 0.0 && lower != 0.0 {
                continue; // negative lower makes no sense for sqrt
            }
            let actual_lower = if lower == 0.0 { 1e-12 } else { lower };
            let coeffs = fit_sqrt_chebyshev(actual_lower, CHEBY_HI_DEFAULT, deg);
            let (max_err, _) =
                max_grid_error(actual_lower, CHEBY_HI_DEFAULT, &coeffs);

            print!("{deg:>7} ",);
            print!("{max_err:>12.3e} ",);
            for &px in probes.iter() {
                if px < actual_lower {
                    print!("{:>10} ", "N/A");
                } else {
                    let err =
                        point_error(px, actual_lower, CHEBY_HI_DEFAULT, &coeffs);
                    print!("{err:>10.3e} ",);
                }
            }
            println!();
        }
        println!();
    }

    println!("(N/A = probe point below fit interval lower bound)\n");
}

// ── Section 2: Sensor std-dev pipeline ───────────────────────────────────

fn section2() {
    println!("=== Section 2: Sensor std-dev pipeline ===\n");

    let data: Vec<f64> = SENSOR_DATA.to_vec();
    let mean = sample_mean(&data);
    let var = sample_variance(&data);
    let std_true = sample_std_dev(&data);

    println!("Input ({SENSOR_COUNT} sensors, range [{SENSOR_LO}, {SENSOR_HI}]):");
    println!("  {SENSOR_DATA:?}\n");
    println!("  true mean     = {mean:.8}");
    println!("  true variance = {var:.8}  (sample, N-1 denominator)");
    println!("  true std_dev  = {std_true:.8}\n");

    println!("Theoretical bounds:");
    println!("  max variance = {MAX_VAR_POP:.1} (population)");
    println!("  max variance = {MAX_VAR_SAMPLE:.3} (sample, N={SENSOR_COUNT})");
    println!("  max std_dev  = {:.3}", MAX_VAR_POP.sqrt());
    println!("  max std_dev  = {:.3} (sample)\n", MAX_VAR_SAMPLE.sqrt());

    println!(
        "Chebyshev fit: [{CHEBY_LO_DEFAULT}, {CHEBY_HI_DEFAULT}], \
         degree {CHEBY_DEGREE_DEFAULT}\n",
    );

    let approx = SqrtApproximator::new(
        CHEBY_LO_DEFAULT,
        CHEBY_HI_DEFAULT,
        CHEBY_DEGREE_DEFAULT,
    );

    println!(
        "{:>8}  {:>12}  {:>12}  {:>12}  {:>16}  {:>10}",
        "B", "var_norm", "sqrt_norm_c", "sqrt_norm_t", "std_approx", "err"
    );
    println!(
        "{:>8}  {:>12}  {:>12}  {:>12}  {:>16}  {:>10}",
        "", "", "(cheb)", "(true)", "", ""
    );

    for &b in NORM_BOUNDS.iter() {
        let var_norm = var / b;
        let sqrt_norm_ch = approx.eval(var_norm);
        let sqrt_norm_tr = var_norm.sqrt();
        let std_approx = sqrt_norm_ch * b.sqrt();
        let err = (std_approx - std_true).abs();

        // Check if var_norm is inside fit interval
        let clipped = if var_norm < CHEBY_LO_DEFAULT {
            format!("(clipped @ {CHEBY_LO_DEFAULT}) {var_norm:.6e}")
        } else if var_norm > CHEBY_HI_DEFAULT {
            format!("(out of range!) {var_norm:.2e}")
        } else {
            format!("{var_norm:.6}")
        };

        let status = if err < ERROR_THRESHOLD { "OK" } else { "FAIL" };

        println!(
            "{b:>8.1}  {clipped:>12}  {sqrt_norm_ch:>12.8}  \
             {sqrt_norm_tr:>12.8}  {std_approx:>16.10}  {err:>8.1e}  {status}"
        );
    }
    println!();
}

// ── Section 3: Edge cases ────────────────────────────────────────────────

fn section3() {
    println!("=== Section 3: Edge cases ===\n");

    let approx = SqrtApproximator::new(
        CHEBY_LO_DEFAULT,
        CHEBY_HI_DEFAULT,
        CHEBY_DEGREE_DEFAULT,
    );

    // Case A: all same (variance = 0)
    {
        let data: Vec<f64> = vec![20.0; SENSOR_COUNT];
        let var = sample_variance(&data);
        let std_true = sample_std_dev(&data);
        println!("Case A: all sensors = 20.0");
        println!("  true variance = {var:.8}");
        println!("  true std_dev  = {std_true:.8}");
        for &b in NORM_BOUNDS.iter() {
            let (std_approx, _) = pipeline_std(&data, b, &approx);
            let err = (std_approx - std_true).abs();
            let status = if err < ERROR_THRESHOLD { "OK" } else { "FAIL" };
            println!(
                "  B={b:>8.1}: std_approx={std_approx:.6}, \
                 err={err:.2e}  {status}"
            );
        }
        println!();
    }

    // Case B: extreme (half at min, half at max)
    {
        let half = SENSOR_COUNT / 2;
        let other_half = SENSOR_COUNT - half;
        let mut data = vec![SENSOR_LO; half];
        data.extend(vec![SENSOR_HI; other_half]);
        let var = sample_variance(&data);
        let std_true = sample_std_dev(&data);
        println!(
            "Case B: extreme [{SENSOR_LO} x {half}, \
             {SENSOR_HI} x {other_half}]"
        );
        println!("  true variance = {var:.4}");
        println!("  true std_dev  = {std_true:.4}");
        for &b in NORM_BOUNDS.iter() {
            let (std_approx, _) = pipeline_std(&data, b, &approx);
            let err = (std_approx - std_true).abs();
            let status = if err < ERROR_THRESHOLD { "OK" } else { "FAIL" };
            let var_norm = var / b;
            println!(
                "  B={b:>8.1}: var_norm={var_norm:.4e}, std_approx={std_approx:.4}, \
                 err={err:.2e}  {status}"
            );
        }
        println!();
    }

    // Case C: tiny jitter around 20
    {
        let mut data = Vec::with_capacity(SENSOR_COUNT);
        for i in 0..SENSOR_COUNT {
            if i % 2 == 0 {
                data.push(20.0);
            } else {
                data.push(20.1);
            }
        }
        let var = sample_variance(&data);
        let std_true = sample_std_dev(&data);
        println!("Case C: tiny jitter [20.0, 20.1, 20.0, 20.1, ...]");
        println!("  true variance = {var:.8}");
        println!("  true std_dev  = {std_true:.8}");
        for &b in NORM_BOUNDS.iter() {
            let (std_approx, _) = pipeline_std(&data, b, &approx);
            let err = (std_approx - std_true).abs();
            let var_norm = var / b;
            let clipped = if var_norm < CHEBY_LO_DEFAULT {
                "[clipped]"
            } else {
                ""
            };
            let status = if err < ERROR_THRESHOLD { "OK" } else { "FAIL" };
            println!(
                "  B={b:>8.1}: var_norm={var_norm:.6e} {clipped}, \
                 std_approx={std_approx:.6}, err={err:.2e}  {status}"
            );
        }
        println!();
    }
}

// ── main ──────────────────────────────────────────────────────────────────

fn main() {
    section1();
    section2();
    section3();
}
