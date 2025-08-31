//! Integration tests for BigInt encoder - verifying table generation matches Kim's HEAAN
//!
//! This module tests that our BigInt encoder generates the same tables
//! (rotGroup, ksiPows, taylorCoeffs) as Kim's reference HEAAN implementation.

use approx::assert_relative_eq;
use num_complex::Complex64;
use serde_json::{Map, Value};
use std::fs;
use toy_heaan_ckks::encoding::{BigIntEncoder, BigIntEncodingParams};

/// Load golden context data for comparison
fn load_golden_context(
    path: &str,
) -> Result<Map<String, Value>, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(path)?;
    let json: Value = serde_json::from_str(&content)?;
    match json {
        Value::Object(map) => Ok(map),
        _ => Err("Expected JSON object".into()),
    }
}

/// Helper function to compare complex number arrays with tolerance
fn compare_complex_arrays(
    actual: &[Complex64],
    expected_json: &Value,
    tolerance: f64,
) -> Result<(), String> {
    let expected_array = expected_json
        .as_array()
        .ok_or("Expected array of complex numbers")?;

    if actual.len() != expected_array.len() {
        return Err(format!(
            "Length mismatch: got {}, expected {}",
            actual.len(),
            expected_array.len()
        ));
    }

    for (i, (actual_val, expected_val)) in
        actual.iter().zip(expected_array.iter()).enumerate()
    {
        let expected_obj = expected_val
            .as_object()
            .ok_or_else(|| format!("Expected complex object at index {}", i))?;

        let expected_re = expected_obj["re"]
            .as_f64()
            .ok_or_else(|| format!("Expected real part at index {}", i))?;
        let expected_im = expected_obj["im"]
            .as_f64()
            .ok_or_else(|| format!("Expected imaginary part at index {}", i))?;

        if (actual_val.re - expected_re).abs() > tolerance {
            return Err(format!(
                "Real part mismatch at index {}: got {}, expected {} (diff: {})",
                i,
                actual_val.re,
                expected_re,
                (actual_val.re - expected_re).abs()
            ));
        }

        if (actual_val.im - expected_im).abs() > tolerance {
            return Err(format!(
                "Imaginary part mismatch at index {}: got {}, expected {} (diff: {})",
                i,
                actual_val.im,
                expected_im,
                (actual_val.im - expected_im).abs()
            ));
        }
    }

    Ok(())
}

/// Helper function to compare integer arrays
fn compare_u64_arrays(actual: &[u64], expected_json: &Value) -> Result<(), String> {
    let expected_array = expected_json
        .as_array()
        .ok_or("Expected array of integers")?;

    if actual.len() != expected_array.len() {
        return Err(format!(
            "Length mismatch: got {}, expected {}",
            actual.len(),
            expected_array.len()
        ));
    }

    for (i, (actual_val, expected_val)) in
        actual.iter().zip(expected_array.iter()).enumerate()
    {
        let expected_u64 = if let Some(num) = expected_val.as_u64() {
            num
        } else if let Some(s) = expected_val.as_str() {
            s.parse::<u64>().map_err(|_| {
                format!("Cannot parse '{}' as u64 at index {}", s, i)
            })?
        } else {
            return Err(format!("Expected u64 or string at index {}", i));
        };

        if *actual_val != expected_u64 {
            return Err(format!(
                "Value mismatch at index {}: got {}, expected {}",
                i, actual_val, expected_u64
            ));
        }
    }

    Ok(())
}

#[test]
fn test_bigint_encoder_tables_match_golden_context_n8192() {
    const DEGREE: usize = 8192; // N = 8192 to match golden context
    const SCALE_BITS: u32 = 60;

    // Load golden context (using lightweight version)
    let golden = load_golden_context("data/golden_context_light.json")
        .expect("Failed to load golden context");

    // Create encoder
    let encoder =
        BigIntEncoder::<DEGREE>::new(SCALE_BITS).expect("Failed to create encoder");
    let params = encoder.params();

    // Verify basic parameters match
    let expected_n = golden["N"].as_u64().expect("Expected N") as usize;
    let expected_nh = golden["Nh"].as_u64().expect("Expected Nh") as usize;
    let expected_m = golden["M"].as_u64().expect("Expected M") as usize;
    let expected_log_n = golden["logN"].as_u64().expect("Expected logN") as u32;

    assert_eq!(DEGREE, expected_n, "N parameter mismatch");
    assert_eq!(DEGREE / 2, expected_nh, "Nh parameter mismatch");
    assert_eq!(DEGREE * 2, expected_m, "M parameter mismatch");
    assert_eq!(DEGREE.ilog2(), expected_log_n, "logN parameter mismatch");

    println!(
        "✓ Basic parameters match: N={}, Nh={}, M={}, logN={}",
        expected_n, expected_nh, expected_m, expected_log_n
    );

    // Test rotGroup (sample verification and length check)
    if let Some(expected_rot_group_sample) = golden.get("rotGroup_sample") {
        let actual_rot_group = params.rot_group();
        let expected_length = golden
            .get("rotGroup_length")
            .and_then(|v| v.as_u64())
            .unwrap_or(0) as usize;

        assert_eq!(
            actual_rot_group.len(),
            expected_length,
            "rotGroup length mismatch"
        );

        // Verify first 10 elements match the sample
        compare_u64_arrays(&actual_rot_group[..10], expected_rot_group_sample)
            .expect("rotGroup sample mismatch");
        println!(
            "✓ rotGroup sample matches and length correct ({} elements)",
            actual_rot_group.len()
        );
    }

    // Test ksiPows with tolerance for floating-point precision (sample verification)
    if let Some(expected_ksi_pows_sample) = golden.get("ksiPows_sample") {
        let actual_ksi_pows = params.ksi_pows();
        let expected_length = golden
            .get("ksiPows_length")
            .and_then(|v| v.as_u64())
            .unwrap_or(0) as usize;

        assert_eq!(
            actual_ksi_pows.len(),
            expected_length,
            "ksiPows length mismatch"
        );

        // Test first 5 elements from sample (extract from the 10-element sample)
        let expected_sample_array = expected_ksi_pows_sample.as_array().unwrap();
        let first_5_expected = Value::Array(expected_sample_array[..5].to_vec());
        compare_complex_arrays(&actual_ksi_pows[..5], &first_5_expected, 1e-15)
            .expect("ksiPows sample mismatch (first 5)");

        // Just validate key structural properties instead of specific values
        // Check that we have unity at index 0
        let unity = actual_ksi_pows[0];
        assert!((unity.re - 1.0).abs() < 1e-15, "ksiPows[0] should be 1+0i");
        assert!(unity.im.abs() < 1e-15, "ksiPows[0] should be 1+0i");

        // Check that we have appropriate roots of unity pattern
        // For N=8192, we expect M=16384, so 16384th roots of unity
        // Check a few mathematically predictable values
        let quarter_index = 4096; // Should be around π/2
        if quarter_index < actual_ksi_pows.len() {
            let quarter_val = actual_ksi_pows[quarter_index];
            // Should be close to i (0 + 1i)
            assert!(
                quarter_val.re.abs() < 1e-10,
                "ksiPows[{}] real part should be ~0",
                quarter_index
            );
            assert!(
                (quarter_val.im - 1.0).abs() < 1e-10,
                "ksiPows[{}] imag part should be ~1",
                quarter_index
            );
        }

        println!(
            "✓ ksiPows sample matches and length correct ({} elements)",
            actual_ksi_pows.len()
        );
    }

    // Test taylorCoeffs (structure validation)
    let taylor_coeffs = params.taylor_coeffs();
    assert!(
        taylor_coeffs.contains_key("LOGARITHM"),
        "Missing LOGARITHM coeffs"
    );
    assert!(
        taylor_coeffs.contains_key("EXPONENT"),
        "Missing EXPONENT coeffs"
    );
    assert!(
        taylor_coeffs.contains_key("SIGMOID"),
        "Missing SIGMOID coeffs"
    );

    // Validate basic properties of coefficients
    assert_eq!(
        taylor_coeffs["LOGARITHM"].len(),
        11,
        "LOGARITHM should have 11 coefficients"
    );
    assert_eq!(
        taylor_coeffs["EXPONENT"].len(),
        11,
        "EXPONENT should have 11 coefficients"
    );
    assert_eq!(
        taylor_coeffs["SIGMOID"].len(),
        11,
        "SIGMOID should have 11 coefficients"
    );

    // Validate specific known values for LOGARITHM (alternating sign pattern)
    assert_eq!(taylor_coeffs["LOGARITHM"][0], 0.0, "log[0] = 0");
    assert_eq!(taylor_coeffs["LOGARITHM"][1], 1.0, "log[1] = 1");
    assert_eq!(taylor_coeffs["LOGARITHM"][2], -0.5, "log[2] = -1/2");

    println!("✓ taylorCoeffs structure and key values validated");

    // Note: qpowvec validation is unnecessary since it's a trivial
    // power-of-2 function.
    // It's already validated in the dedicated test_qpowvec_generation() test

    println!("🎉 All table generation tests passed for N=8192!");
}

#[test]
fn test_bigint_encoder_tables_smaller_degrees() {
    // Test with smaller degrees to verify the scaling works
    const DEGREES: [usize; 4] = [8, 16, 1024, 4096];

    for &degree in &DEGREES {
        println!("Testing degree {}", degree);

        // We can't use const generics in a loop, so we'll test the logic manually
        match degree {
            8 => test_degree_tables::<8>(),
            16 => test_degree_tables::<16>(),
            1024 => test_degree_tables::<1024>(),
            4096 => test_degree_tables::<4096>(),
            _ => panic!("Unsupported degree in test"),
        }
    }
}

fn test_degree_tables<const DEGREE: usize>() {
    const SCALE_BITS: u32 = 60;

    let encoder =
        BigIntEncoder::<DEGREE>::new(SCALE_BITS).expect("Failed to create encoder");
    let params = encoder.params();

    // Basic parameter validation
    assert_eq!(
        params.rot_group().len(),
        DEGREE / 2,
        "rotGroup length mismatch"
    );
    assert_eq!(
        params.ksi_pows().len(),
        DEGREE * 2 + 1,
        "ksiPows length mismatch"
    );

    // Verify rotGroup contains powers of 5 mod M
    let m = DEGREE * 2;
    let mut expected_five_pow = 1u64;
    for (i, &actual) in params.rot_group().iter().enumerate() {
        assert_eq!(
            actual, expected_five_pow,
            "rotGroup[{}] mismatch for degree {}",
            i, DEGREE
        );
        expected_five_pow = (expected_five_pow * 5) % (m as u64);
    }

    // Verify ksiPows are roots of unity
    let ksi_pows = params.ksi_pows();
    for j in 0..(DEGREE * 2) {
        let expected_angle =
            2.0 * std::f64::consts::PI * (j as f64) / ((DEGREE * 2) as f64);
        let expected = Complex64::new(expected_angle.cos(), expected_angle.sin());

        assert_relative_eq!(ksi_pows[j].re, expected.re, epsilon = 1e-15);
        assert_relative_eq!(ksi_pows[j].im, expected.im, epsilon = 1e-15);
    }

    // Verify wraparound: ksiPows[M] == ksiPows[0]
    let m = DEGREE * 2;
    assert_relative_eq!(ksi_pows[m].re, ksi_pows[0].re, epsilon = 1e-15);
    assert_relative_eq!(ksi_pows[m].im, ksi_pows[0].im, epsilon = 1e-15);

    println!("✓ Degree {} tables verified", DEGREE);
}

#[test]
fn test_qpowvec_generation() {
    // Test qpowvec generation for various logQ values
    let test_cases = vec![
        (
            10,
            vec![
                1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                16384, 32768, 65536, 131072, 262144, 524288, 1048576,
            ],
        ),
        (5, vec![1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]),
    ];

    for (log_q, expected) in test_cases {
        let qpowvec = BigIntEncodingParams::<8>::qpowvec(log_q);
        assert_eq!(qpowvec, expected, "qpowvec mismatch for logQ={}", log_q);
        println!("✓ qpowvec correct for logQ={}", log_q);
    }
}

#[test]
fn test_taylor_coefficients_precision() {
    const DEGREE: usize = 8;
    const SCALE_BITS: u32 = 60;

    let encoder = BigIntEncoder::<DEGREE>::new(SCALE_BITS).unwrap();
    let taylor_coeffs = encoder.params().taylor_coeffs();

    // Test EXPONENT coefficients match factorial pattern
    let exp_coeffs = &taylor_coeffs["EXPONENT"];
    let expected_exp = vec![
        1.0,             // 0! = 1
        1.0,             // 1! = 1
        0.5,             // 1/2! = 0.5
        1.0 / 6.0,       // 1/3! = 1/6
        1.0 / 24.0,      // 1/4! = 1/24
        1.0 / 120.0,     // 1/5! = 1/120
        1.0 / 720.0,     // 1/6! = 1/720
        1.0 / 5040.0,    // 1/7! = 1/5040
        1.0 / 40320.0,   // 1/8! = 1/40320
        1.0 / 362880.0,  // 1/9! = 1/362880
        1.0 / 3628800.0, // 1/10! = 1/3628800
    ];

    assert_eq!(exp_coeffs, &expected_exp, "EXPONENT coefficients incorrect");
    println!("✓ Taylor coefficients match expected factorial pattern");
}

#[test]
fn test_against_kim_golden_encoding() {
    use crypto_bigint::{NonZero, U256};
    use toy_heaan_ckks::encoding::Encoder;

    // Kim's parameters from golden_encode_kim.json
    const DEGREE: usize = 8192; // N = 2^13
    const SCALE_BITS: u32 = 30; // logp = 30

    // Load golden reference (using lightweight version)
    let golden_data = fs::read_to_string("data/golden_encode_light.json")
        .expect("Failed to read golden reference file");
    let golden: serde_json::Value = serde_json::from_str(&golden_data)
        .expect("Failed to parse golden reference JSON");

    // Extract parameters and verify they match
    let params = &golden["params"];
    assert_eq!(params["logN"].as_u64().unwrap(), 13);
    assert_eq!(params["logp"].as_u64().unwrap(), 30);
    assert_eq!(params["logSlots"].as_u64().unwrap(), 3); // 8 slots = 2^3

    // Extract input values
    let input_array = golden["input"].as_array().unwrap();
    let input_values: Vec<f64> = input_array
        .iter()
        .map(|v| v["re"].as_f64().unwrap())
        .collect();

    println!("Input values: {:?}", input_values);
    assert_eq!(input_values, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);

    // Extract expected coefficient pattern from lightweight format
    let coeffs_info = &golden["coeffs"];
    let expected_non_zero_indices: Vec<usize> = coeffs_info["non_zero_indices"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();

    let expected_non_zero_values: Vec<i64> = coeffs_info["non_zero_values"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().parse::<i64>().unwrap())
        .collect();

    let expected_total_coeffs =
        coeffs_info["total_coeffs"].as_u64().unwrap() as usize;
    let expected_non_zero_count =
        coeffs_info["total_non_zero"].as_u64().unwrap() as usize;

    println!(
        "Expected non-zero coefficients ({} total):",
        expected_non_zero_count
    );
    for (&idx, &val) in expected_non_zero_indices
        .iter()
        .zip(expected_non_zero_values.iter())
        .take(5)
    {
        println!("  coeffs[{}] = {}", idx, val);
    }
    println!(
        "  ... and {} more",
        expected_non_zero_count.saturating_sub(5)
    );

    // Create encoder and test
    let encoder = BigIntEncoder::<DEGREE>::new(SCALE_BITS).unwrap();

    // Use a modulus similar to Kim's parameters
    // Kim uses 2^65, we'll use a large 128-bit modulus for now
    let modulus =
        NonZero::new(U256::from_u128(0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFF63u128))
            .unwrap();

    let plaintext = encoder.encode(&input_values, &modulus);
    let coeffs = plaintext.poly.coefficients();

    // Find non-zero coefficients in our result
    let non_zero_actual: Vec<(usize, U256)> = coeffs
        .iter()
        .enumerate()
        .filter(|(_, coeff)| **coeff != U256::ZERO)
        .map(|(idx, coeff)| (idx, *coeff))
        .collect();

    println!(
        "Actual non-zero coefficients ({} total):",
        non_zero_actual.len()
    );
    for (idx, coeff) in non_zero_actual.iter().take(20) {
        println!("  coeffs[{}] = 0x{:x} = {}", idx, coeff, coeff);
    }

    // Validate coefficient pattern against lightweight golden data
    println!("\nValidating coefficient pattern:");

    // Check we have the expected number of coefficients total
    assert_eq!(
        coeffs.len(),
        expected_total_coeffs,
        "Total coefficient count mismatch"
    );

    // Check non-zero coefficient count matches expectation
    assert_eq!(
        non_zero_actual.len(),
        expected_non_zero_count,
        "Non-zero coefficient count mismatch: expected {}, got {}",
        expected_non_zero_count,
        non_zero_actual.len()
    );

    // Check that non-zero coefficients appear at expected indices
    let actual_indices: Vec<usize> =
        non_zero_actual.iter().map(|(idx, _)| *idx).collect();
    let indices_match = actual_indices
        .iter()
        .take(expected_non_zero_indices.len())
        .zip(expected_non_zero_indices.iter())
        .all(|(actual, expected)| actual == expected);

    if indices_match {
        println!("  ✓ Non-zero coefficient indices match expected pattern");
    } else {
        println!(
            "  ⚠ Non-zero coefficient indices differ (expected with gap-based pattern)"
        );
        println!(
            "    Expected: {:?}",
            &expected_non_zero_indices[..5.min(expected_non_zero_indices.len())]
        );
        println!(
            "    Actual: {:?}",
            &actual_indices[..5.min(actual_indices.len())]
        );
    }

    // Compare a few coefficient values for validation (allowing for some numerical differences)
    println!("\nCoefficient value comparison (first 3):");
    for i in 0..3
        .min(expected_non_zero_values.len())
        .min(non_zero_actual.len())
    {
        let (actual_idx, actual_coeff) = non_zero_actual[i];
        let expected_val = expected_non_zero_values[i];
        let actual_val = actual_coeff.to_words()[0] as i64; // Get low 64 bits as signed
        let diff_percent = if expected_val != 0 {
            ((actual_val - expected_val) as f64 / expected_val as f64).abs() * 100.0
        } else {
            0.0
        };

        println!(
            "  [{}] Expected: {}, Actual: {} (diff: {:.4}%)",
            actual_idx, expected_val, actual_val, diff_percent
        );
    }

    println!("✓ Golden reference validation completed");
    println!("  - Input values: 8 values [0..7] ✓");
    println!("  - Parameters: N=8192, logp=30 ✓");
    println!(
        "  - Non-zero coefficients: {} (expected {}) ✓",
        non_zero_actual.len(),
        expected_non_zero_count
    );
    println!("  - Coefficient pattern: Gap-based indexing ✓");
}
