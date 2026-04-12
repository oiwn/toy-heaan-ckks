use rand::{Rng, RngExt, seq::SliceRandom};
use rand_distr::{Distribution, Normal, uniform::Uniform};

/// Samples uniform integer coefficients in `[0, max_value)`.
///
/// # Panics
///
/// Panics if `max_value == 0`.
pub fn uniform_coefficients<const DEGREE: usize, R: Rng + ?Sized>(
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    let distribution = Uniform::new(0, max_value).unwrap_or_else(|_| {
        panic!(
            "uniform_coefficients: invalid range [0, {max_value}), \
            max_value must be positive"
        )
    });
    let mut coeffs = [0u64; DEGREE];
    for coeff in &mut coeffs {
        *coeff = distribution.sample(rng);
    }
    coeffs
}

/// Samples rounded Gaussian integers and maps them into `[0, max_value)`.
///
/// # Panics
///
/// Panics if `std_dev` is not finite and positive, or if `max_value == 0`.
pub fn gaussian_coefficients<const DEGREE: usize, R: Rng + ?Sized>(
    std_dev: f64,
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    assert!(
        std_dev.is_finite() && std_dev > 0.0,
        "gaussian_coefficients: std_dev must be finite and positive"
    );
    assert!(
        max_value > 0,
        "gaussian_coefficients: max_value must be positive"
    );
    let normal = Normal::new(0.0, std_dev)
        .expect("gaussian_coefficients: failed to create Normal distribution");
    let mut coeffs = [0u64; DEGREE];

    for coeff in &mut coeffs {
        let sample = normal.sample(rng);
        let noise_int = sample.round() as i64;

        // Convert to unsigned: handle negative values properly
        *coeff = if noise_int < 0 {
            let abs_val = noise_int.unsigned_abs() % max_value;
            if abs_val == 0 { 0 } else { max_value - abs_val }
        } else {
            (noise_int as u64) % max_value
        };
    }

    coeffs
}

/// Samples a ternary vector with coefficients in `{-1, 0, 1}`.
///
/// Exactly `hamming_weight` entries are non-zero.
///
/// # Panics
///
/// Panics if `hamming_weight > DEGREE`.
pub fn ternary_coefficients<const DEGREE: usize, R: Rng + ?Sized>(
    hamming_weight: usize,
    rng: &mut R,
) -> [i64; DEGREE] {
    assert!(
        hamming_weight <= DEGREE,
        "ternary_coefficients: hamming_weight must be <= DEGREE"
    );
    let mut out = [0i64; DEGREE];
    // Shuffle indices and assign signs on the selected support.
    let mut indices: Vec<usize> = (0..DEGREE).collect();
    indices.shuffle(rng);
    for &idx in indices.iter().take(hamming_weight) {
        out[idx] = if rng.random_bool(0.5) { 1 } else { -1 };
    }
    out
}

#[cfg(test)]
mod tests {
    use super::{
        gaussian_coefficients, ternary_coefficients, uniform_coefficients,
    };
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn uniform_coefficients_stay_in_range() {
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let coeffs = uniform_coefficients::<128, _>(17, &mut rng);
        for &coeff in &coeffs {
            assert!(coeff < 17);
        }
    }

    #[test]
    #[should_panic(
        expected = "uniform_coefficients: invalid range [0, 0), max_value must be positive"
    )]
    fn uniform_coefficients_panics_on_zero_max_value() {
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let _ = uniform_coefficients::<8, _>(0, &mut rng);
    }

    #[test]
    fn uniform_coefficients_are_roughly_balanced() {
        const DEGREE: usize = 8192;
        const MODULUS: usize = 8;
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let coeffs = uniform_coefficients::<DEGREE, _>(MODULUS as u64, &mut rng);

        let mut buckets = [0usize; MODULUS];
        for &coeff in &coeffs {
            buckets[coeff as usize] += 1;
        }

        let expected = DEGREE as f64 / MODULUS as f64;
        for &count in &buckets {
            let deviation = (count as f64 - expected).abs();
            assert!(
                deviation <= expected * 0.30,
                "bucket count {count} too far from expected {expected}"
            );
        }
    }

    #[test]
    fn gaussian_coefficients_stay_in_range() {
        let mut rng = ChaCha20Rng::seed_from_u64(11);
        let coeffs = gaussian_coefficients::<256, _>(3.2, 97, &mut rng);
        for &coeff in &coeffs {
            assert!(coeff < 97);
        }
    }

    #[test]
    #[should_panic(
        expected = "gaussian_coefficients: std_dev must be finite and positive"
    )]
    fn gaussian_coefficients_panics_on_non_positive_std_dev() {
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let _ = gaussian_coefficients::<8, _>(0.0, 17, &mut rng);
    }

    #[test]
    #[should_panic(
        expected = "gaussian_coefficients: std_dev must be finite and positive"
    )]
    fn gaussian_coefficients_panics_on_non_finite_std_dev() {
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let _ = gaussian_coefficients::<8, _>(f64::NAN, 17, &mut rng);
    }

    #[test]
    #[should_panic(expected = "gaussian_coefficients: max_value must be positive")]
    fn gaussian_coefficients_panics_on_zero_max_value() {
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let _ = gaussian_coefficients::<8, _>(2.0, 0, &mut rng);
    }

    #[test]
    fn gaussian_coefficients_have_reasonable_mean_and_variance() {
        const DEGREE: usize = 16_384;
        let std_dev = 3.2;
        let max_value = 1_000_003u64;
        let mut rng = ChaCha20Rng::seed_from_u64(99);
        let coeffs =
            gaussian_coefficients::<DEGREE, _>(std_dev, max_value, &mut rng);

        let centered: Vec<f64> = coeffs
            .iter()
            .map(|&x| {
                if x <= max_value / 2 {
                    x as f64
                } else {
                    x as f64 - max_value as f64
                }
            })
            .collect();

        let mean = centered.iter().sum::<f64>() / DEGREE as f64;
        let variance = centered
            .iter()
            .map(|&x| {
                let diff = x - mean;
                diff * diff
            })
            .sum::<f64>()
            / DEGREE as f64;

        let expected_variance = std_dev * std_dev;
        assert!(mean.abs() <= 0.25, "mean too far from 0: {mean}");
        assert!(
            (variance - expected_variance).abs() <= expected_variance * 0.35,
            "variance {variance} too far from expected {expected_variance}"
        );
    }

    #[test]
    fn ternary_coefficients_have_exact_hamming_weight() {
        let mut rng = ChaCha20Rng::seed_from_u64(123);
        let coeffs = ternary_coefficients::<256, _>(31, &mut rng);
        let non_zero = coeffs.iter().filter(|&&x| x != 0).count();
        assert_eq!(non_zero, 31);
    }

    #[test]
    fn ternary_coefficients_are_in_expected_set() {
        let mut rng = ChaCha20Rng::seed_from_u64(321);
        let coeffs = ternary_coefficients::<256, _>(64, &mut rng);
        for &x in &coeffs {
            assert!(x == -1 || x == 0 || x == 1);
        }
    }

    #[test]
    fn ternary_coefficients_handle_weight_extremes() {
        let mut rng = ChaCha20Rng::seed_from_u64(999);
        let all_zero = ternary_coefficients::<64, _>(0, &mut rng);
        assert!(all_zero.iter().all(|&x| x == 0));

        let full = ternary_coefficients::<64, _>(64, &mut rng);
        assert!(full.iter().all(|&x| x == -1 || x == 1));
    }

    #[test]
    #[should_panic(
        expected = "ternary_coefficients: hamming_weight must be <= DEGREE"
    )]
    fn ternary_coefficients_panics_on_oversized_hamming_weight() {
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let _ = ternary_coefficients::<8, _>(9, &mut rng);
    }
}
