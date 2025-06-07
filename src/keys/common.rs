use crate::PolyRing;
use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Generate a polynomial with uniformly random coefficients
pub fn generate_random_poly<R: Rng>(
    n: usize,
    modulus: u64,
    ring_dim: usize,
    rng: &mut R,
) -> PolyRing {
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Use gen_range for unbiased sampling
        let random_val = rng.random_range(0..modulus);
        coeffs.push(random_val);
    }

    PolyRing::from_coeffs(&coeffs, modulus, ring_dim)
}

/// Generate a polynomial with coefficients from discrete Gaussian distribution
pub fn generate_error_poly<R: Rng>(
    n: usize,
    modulus: u64,
    variance: f64,
    ring_dim: usize,
    rng: &mut R,
) -> PolyRing {
    let std_dev = variance.sqrt();
    let normal = Normal::new(0.0, std_dev).expect("Invalid standard deviation");
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Sample from Gaussian distribution and round to nearest integer
        let noise = normal.sample(rng);
        let noise_int = noise.round() as i64;

        // Convert to unsigned, correctly handling negative values with modular arithmetic
        let coeff = if noise_int < 0 {
            modulus - (noise_int.unsigned_abs() % modulus)
        } else {
            noise_int as u64 % modulus
        };

        coeffs.push(coeff);
    }

    PolyRing::from_coeffs(&coeffs, modulus, ring_dim)
}

/// sparsity controls what fraction of coefficients will be non-zero (0.0-1.0)
pub fn generate_ternary_poly<T: Rng>(
    n: usize,
    modulus: u64,
    sparsity: f64,
    rng: &mut T,
) -> PolyRing {
    // Calculate hamming weight from sparsity
    let hamming_weight = (n as f64 * sparsity).round() as usize;

    // Initialize all coefficients to zero
    let mut coeffs = vec![0u64; n];
    let mut count = 0;

    // Add exactly hamming_weight non-zero coefficients
    while count < hamming_weight {
        let idx = rng.random_range(0..n);

        // Only set if position is still zero
        if coeffs[idx] == 0 {
            // Generate either 1 or -1 with equal probability
            coeffs[idx] = if rng.random_bool(0.5) {
                1
            } else {
                modulus - 1 // -1 mod q
            };
            count += 1;
        }
    }

    PolyRing::from_coeffs(&coeffs, modulus, n)
}

/// Compute the negative of a polynomial in the ring
pub fn negate_poly(poly: &PolyRing, modulus: u64) -> PolyRing {
    // Create a polynomial with coefficients (q - a_i) for each coefficient a_i
    let mut neg_coeffs =
        Vec::with_capacity((poly.polynomial_degree() + 1) as usize);

    for coeff in poly {
        if *coeff == 0 {
            neg_coeffs.push(0);
        } else {
            neg_coeffs.push(modulus - *coeff);
        }
    }

    PolyRing::from_coeffs(&neg_coeffs, modulus, poly.ring_dim())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_uniform_distribution() {
        let mut rng = ChaCha20Rng::seed_from_u64(0);
        let modulus = 17;
        let poly_len = 8;
        let samples_count = 100000;

        let mut counts = vec![0usize; modulus as usize];

        for _ in 0..samples_count {
            let poly = generate_random_poly(poly_len, modulus, poly_len, &mut rng);
            for &coeff in poly.into_iter() {
                counts[coeff as usize] += 1;
            }
        }

        let expected = samples_count as f64 * poly_len as f64 / modulus as f64;
        let tolerance_ratio = 1.2;

        for &count in &counts {
            assert!(
                (count as f64 - expected).abs() / expected < tolerance_ratio,
                "Distribution not uniform enough: count {}, expected {}",
                count,
                expected
            );
        }
    }

    #[test]
    fn test_gaussian_distribution() {
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let modulus = 17;
        let poly_len = 8;
        let variance = 3.0;
        let samples_count = 10000;

        let mut values = Vec::with_capacity(samples_count * poly_len);

        for _ in 0..samples_count {
            let poly = generate_error_poly(
                poly_len, modulus, variance, poly_len, &mut rng,
            );
            for &coeff in poly.into_iter() {
                let val = coeff as i64;
                let val = if val > modulus as i64 / 2 {
                    val - modulus as i64
                } else {
                    val
                };
                values.push(val as f64);
            }
        }

        let mean: f64 = values.iter().copied().sum::<f64>() / values.len() as f64;
        let std_dev = (values.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / values.len() as f64)
            .sqrt();

        let expected_std_dev = variance.sqrt();
        let tolerance = 0.2;

        assert!(
            (std_dev - expected_std_dev).abs() < tolerance,
            "Standard deviation deviates too much: actual {}, expected {}",
            std_dev,
            expected_std_dev
        );
    }

    #[test]
    fn test_ternary_poly_hamming_weight() {
        let n = 1000;
        let modulus = 65537;
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Test different sparsity levels
        let sparsity_levels = [0.1, 0.3, 0.5, 0.7, 0.9];

        for &sparsity in &sparsity_levels {
            let expected_weight = (n as f64 * sparsity).round() as usize;
            let poly = generate_ternary_poly(n, modulus, sparsity, &mut rng);

            // Count non-zero coefficients
            let actual_weight = poly.into_iter().filter(|&&x| x != 0).count();

            assert_eq!(
                actual_weight, expected_weight,
                "Hamming weight doesn't match for sparsity {}",
                sparsity
            );
        }
    }

    #[test]
    fn test_ternary_poly_distribution() {
        let n = 10000;
        let modulus = 65537;
        let sparsity = 0.6;
        let mut rng = ChaCha20Rng::seed_from_u64(123);

        let poly = generate_ternary_poly(n, modulus, sparsity, &mut rng);

        // Count occurrences of each value
        let ones_count = poly.into_iter().filter(|&&x| x == 1).count();
        let neg_ones_count =
            poly.into_iter().filter(|&&x| x == modulus - 1).count();
        let zeros_count = poly.into_iter().filter(|&&x| x == 0).count();

        // Calculate hamming weight
        let expected_weight = (n as f64 * sparsity).round() as usize;
        let actual_weight = ones_count + neg_ones_count;

        // Check total non-zero count
        assert_eq!(actual_weight, expected_weight);

        // Check zeros count
        assert_eq!(zeros_count, n - expected_weight);

        // Check distribution of 1 and -1 (should be approximately equal)
        let ratio = ones_count as f64 / actual_weight as f64;
        assert!(
            (ratio - 0.5).abs() < 0.1,
            "1 and -1 distribution ratio is {}, expected close to 0.5",
            ratio
        );
    }

    #[test]
    fn test_deterministic_generation() {
        let n = 100;
        let modulus = 65537;
        let sparsity = 0.4;

        // Create two RNGs with the same seed
        let mut rng1 = ChaCha20Rng::seed_from_u64(42);
        let mut rng2 = ChaCha20Rng::seed_from_u64(42);

        let poly1 = generate_ternary_poly(n, modulus, sparsity, &mut rng1);
        let poly2 = generate_ternary_poly(n, modulus, sparsity, &mut rng2);

        // Polynomials should be identical with same seed
        for i in 0..n {
            assert_eq!(
                poly1.into_iter().nth(i),
                poly2.into_iter().nth(i),
                "Coefficients at position {} differ",
                i
            );
        }
    }

    #[test]
    fn test_edge_cases() {
        let n = 100;
        let modulus = 65537;
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Test zero sparsity (all zeros)
        let zero_poly = generate_ternary_poly(n, modulus, 0.0, &mut rng);
        let non_zero_count = zero_poly.into_iter().filter(|&&x| x != 0).count();
        assert_eq!(non_zero_count, 0, "Expected all zeros with sparsity 0.0");

        // Test full sparsity (all non-zero)
        let full_poly = generate_ternary_poly(n, modulus, 1.0, &mut rng);
        let zero_count = full_poly.into_iter().filter(|&&x| x == 0).count();
        assert_eq!(zero_count, 0, "Expected no zeros with sparsity 1.0");
    }
}
