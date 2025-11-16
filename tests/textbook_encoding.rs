use std::sync::Arc;

use approx::assert_relative_eq;
use rustfft::num_complex::Complex64;

use toy_heaan_ckks::encoding::{Encoder, TextbookEncoder, TextbookEncodingParams};
use toy_heaan_ckks::rings::backends::rns::params::toy_basis;

const SCALE_BITS: u32 = 30;

#[test]
fn textbook_encode_decode_roundtrip_degree8() {
    const DEGREE: usize = 8;
    let basis = toy_basis::<DEGREE>().expect("toy basis");
    let params = TextbookEncodingParams::<DEGREE>::new(
        SCALE_BITS,
        Arc::new(basis.primes().clone()),
    )
    .expect("params");
    let encoder = TextbookEncoder::new(params);

    let values = vec![
        Complex64::new(0.5, 0.0),
        Complex64::new(-1.25, 0.75),
        Complex64::new(0.0, -0.125),
    ];

    let plaintext = encoder.encode_complex(&values, &basis).expect("encode");
    let decoded = encoder.decode_complex(&plaintext).expect("decode");

    for (expected, actual) in values.iter().zip(decoded.iter()) {
        assert_relative_eq!(expected.re, actual.re, epsilon = 1e-6);
        assert_relative_eq!(expected.im, actual.im, epsilon = 1e-6);
    }
}

#[test]
fn textbook_real_roundtrip_via_trait() {
    const DEGREE: usize = 16;
    let basis = toy_basis::<DEGREE>().expect("toy basis");
    let params = TextbookEncodingParams::<DEGREE>::new(
        SCALE_BITS,
        Arc::new(basis.primes().clone()),
    )
    .expect("params");
    let encoder = TextbookEncoder::new(params);

    let reals = vec![0.125, -0.5, 1.75, 0.0, -0.0625];
    let plaintext = encoder.encode(&reals, &basis);
    let decoded = encoder.decode(&plaintext);

    for (expected, actual) in reals.iter().zip(decoded.iter()) {
        assert_relative_eq!(expected, actual, epsilon = 1e-6);
    }
}
