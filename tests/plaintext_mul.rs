#[cfg(test)]
mod tests {
    use toy_heaan_ckks::{
        PolyRing,
        encoding::{self, EncodingParams},
        rescale_poly,
    };

    // Helper function to convert polynomial to signed coefficients
    fn poly_to_coeffs(poly: &PolyRing) -> Vec<i64> {
        let modulus = poly.modulus();
        let half_modulus = modulus / 2;

        poly.into_iter()
            .map(|&c| {
                if c > half_modulus {
                    -((modulus - c) as i64)
                } else {
                    c as i64
                }
            })
            .collect()
    }

    #[test]
    fn test_plaintext_multiplication_roundtrip() {
        // 1) CKKS params
        let ring_degree = 8;
        let scale_bits = 30;
        let modulus = (1u64 << 61) - 1; // same modulus you use elsewhere
        let scale = (1u64 << scale_bits) as f64;

        // 2) Set up a tiny vector pair
        let values1 = vec![1.0, 2.0, 3.0, 4.0];
        let values2 = vec![5.0, 6.0, 7.0, 8.0];
        let expected: Vec<f64> =
            values1.iter().zip(&values2).map(|(&a, &b)| a * b).collect();

        // 3) Encoding → polynomial
        let enc = EncodingParams::new(ring_degree, scale_bits).unwrap();
        let coeffs1 = encoding::encode(&values1, &enc).unwrap();
        let coeffs2 = encoding::encode(&values2, &enc).unwrap();
        let p1 = PolyRing::from_coeffs(&coeffs1, modulus, ring_degree);
        let p2 = PolyRing::from_coeffs(&coeffs2, modulus, ring_degree);

        // 4) Multiply then “rescale” exactly the same way your ciphertext does:
        let raw_prod = p1 * p2;
        let downscaled = rescale_poly(&raw_prod, scale);

        // 5) Pull out signed coefficients and decode
        let prod_coeffs = poly_to_coeffs(&downscaled);
        let result = encoding::decode(&prod_coeffs, &enc).unwrap();

        // 6) Assert exact recovery
        for (i, (&r, &e)) in result.iter().zip(&expected).enumerate() {
            assert!(
                (r - e).abs() < 1e-8,
                "at index {}: got {}, expected {}",
                i,
                r,
                e
            );
        }
    }
}
