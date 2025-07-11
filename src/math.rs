use rand::{Rng, seq::SliceRandom};
use rand_distr::{Distribution, Normal};

// Sample uniformly random integer coefficients in range [0, max_value)
pub fn uniform_coefficients<const DEGREE: usize, R: Rng + ?Sized>(
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    let mut coeffs = [0u64; DEGREE];
    for i in 0..DEGREE {
        coeffs[i] = rng.random_range(0..max_value);
    }
    coeffs
}

/// Sample Gaussian noise as signed integers, then convert to unsigned mod max_value
pub fn gaussian_coefficients<const DEGREE: usize, R: Rng + ?Sized>(
    std_dev: f64,
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    let normal = Normal::new(0.0, std_dev).expect("Invalid Gaussian std_dev");
    let mut coeffs = [0u64; DEGREE];

    for i in 0..DEGREE {
        let sample = normal.sample(rng);
        let noise_int = sample.round() as i64;

        // Convert to unsigned: handle negative values properly
        coeffs[i] = if noise_int < 0 {
            let abs_val = noise_int.unsigned_abs() % max_value;
            if abs_val == 0 { 0 } else { max_value - abs_val }
        } else {
            (noise_int as u64) % max_value
        };
    }

    coeffs
}

/// Sample a ternary polynomial (coeffs in \{-1,0,1\}) with given Hamming weight.
pub fn ternary_coefficients<const DEGREE: usize, R: Rng + ?Sized>(
    hamming_weight: usize,
    rng: &mut R,
) -> [i64; DEGREE] {
    let mut out = [0i64; DEGREE];
    // shuffle indices and pick the first `hamming_weight`
    let mut indices: Vec<usize> = (0..DEGREE).collect();
    indices.shuffle(rng);
    for &idx in indices.iter().take(hamming_weight) {
        out[idx] = if rng.random_bool(0.5) { 1 } else { -1 };
    }
    out
}
