//! Plaintext Chebyshev approximation demo for `sqrt(x)`.
//!
//! This does not use encryption. It is a math sandbox for the std-dev demo:
//! fit `sqrt(x)` on a variance interval, evaluate the polynomial in the
//! Chebyshev basis with Clenshaw's recurrence, and inspect the approximation
//! error before porting anything into CKKS.

const LOWER: f64 = -1.00;
const UPPER: f64 = 1.00;
const DEGREE: usize = 7;
const GRID_POINTS: usize = 20_001;
const CURRENT_VARIANCE: f64 = 0.014_093_75;

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

fn print_degree_sweep() {
    println!("Degree sweep on [{LOWER}, {UPPER}]");
    println!("{:>6}  {:>14}  {:>12}", "degree", "max_error", "at_x");

    for degree in [3usize, 5, 7, 9, 11] {
        let coeffs = fit_sqrt_chebyshev(LOWER, UPPER, degree);
        let (error, x) = max_grid_error(LOWER, UPPER, &coeffs);
        println!("{degree:>6}  {error:>14.6e}  {x:>12.8}");
    }
}

fn print_coefficients(coeffs: &[f64]) {
    println!();
    println!("Chebyshev coefficients for degree {DEGREE}");
    println!("p(x) = sum c_k * T_k(t), t = (2*x - (hi + lo)) / (hi - lo)");
    println!("const SQRT_CHEB_COEFFS: [f64; {}] = [", coeffs.len());
    for coeff in coeffs {
        println!("    {coeff:.17e},");
    }
    println!("];");
}

fn print_samples(coeffs: &[f64]) {
    let samples = [LOWER, CURRENT_VARIANCE, 0.02, 0.05, 0.10, 0.16, UPPER];

    println!();
    println!("Sample approximation errors");
    println!(
        "{:>12}  {:>14}  {:>14}  {:>14}",
        "x", "sqrt(x)", "approx", "abs_error"
    );

    for x in samples {
        let expected = x.sqrt();
        let got = approx_sqrt(x, LOWER, UPPER, coeffs);
        let error = (got - expected).abs();
        println!("{x:>12.8}  {expected:>14.10}  {got:>14.10}  {error:>14.6e}");
    }
}

fn main() {
    println!("Plaintext Chebyshev sqrt approximation");
    println!("Interval: [{LOWER}, {UPPER}]");
    println!("Selected degree: {DEGREE}");
    println!();

    print_degree_sweep();

    let coeffs = fit_sqrt_chebyshev(LOWER, UPPER, DEGREE);
    let (max_error, max_x) = max_grid_error(LOWER, UPPER, &coeffs);

    print_coefficients(&coeffs);
    print_samples(&coeffs);

    println!();
    println!(
        "Selected degree {DEGREE}: max grid error {max_error:.6e} at x={max_x:.8}"
    );
}
