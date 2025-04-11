use criterion::{
    BenchmarkId, Criterion, black_box, criterion_group, criterion_main,
};
use toy_heaan_ckks::PolyRing;

// Helper function to create random polynomials for benchmarking
fn create_random_poly(degree: usize, modulus: u64) -> PolyRing {
    let coeffs = (0..=degree)
        .map(|i| (i as u64 * 17 + 11) % 1231231237)
        .collect::<Vec<_>>();
    PolyRing::from_coeffs(coeffs.as_slice(), modulus, 8)
}

fn bench_addition(c: &mut Criterion) {
    let mut group = c.benchmark_group("polynomial_addition");
    let modulus = 1231231237; // Large prime for testing

    // Test different polynomial degrees
    for degree in [100, 1000].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            degree,
            |b, &degree| {
                let p1 = create_random_poly(degree, modulus);
                let p2 = create_random_poly(degree, modulus);

                b.iter(|| {
                    let _result = black_box(p1.clone()) + black_box(p2.clone());
                });
            },
        );
    }
    group.finish();
}

fn bench_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("polynomial_multiplication");
    let modulus = 1231231237;

    // Test different polynomial degrees
    for degree in [100, 1000].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            degree,
            |b, &degree| {
                let p1 = create_random_poly(degree, modulus);
                let p2 = create_random_poly(degree, modulus);

                b.iter(|| {
                    let _result = black_box(p1.clone()) * black_box(p2.clone());
                });
            },
        );
    }
    group.finish();
}

fn bench_modular_reduction(c: &mut Criterion) {
    let mut group = c.benchmark_group("modular_reduction");

    // Test different moduli sizes
    let moduli = [17, 65537, 1231231237];

    for (i, &modulus) in moduli.iter().enumerate() {
        let degree = 100; // Fixed degree for comparing different moduli
        group.bench_with_input(BenchmarkId::new("modulus_bits", i), &i, |b, _| {
            let mut p = create_random_poly(degree, modulus);

            b.iter(|| {
                p.reduce_coeffs();
                black_box(&mut p);
            });
        });
    }
    group.finish();
}

fn bench_complex_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("complex_operations");
    let modulus = 1231231237;
    let degree = 100;

    group.bench_function("add_mul_sequence", |b| {
        let p1 = create_random_poly(degree, modulus);
        let p2 = create_random_poly(degree, modulus);
        let p3 = create_random_poly(degree, modulus);

        b.iter(|| {
            // (p1 + p2) * p3
            let _result = black_box(black_box(p1.clone()) + black_box(p2.clone()))
                * black_box(p3.clone());
        });
    });

    group.bench_function("mul_add_sequence", |b| {
        let p1 = create_random_poly(degree, modulus);
        let p2 = create_random_poly(degree, modulus);
        let p3 = create_random_poly(degree, modulus);

        b.iter(|| {
            // (p1 * p2) + p3
            let _result = black_box(black_box(p1.clone()) * black_box(p2.clone()))
                + black_box(p3.clone());
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_addition,
    bench_multiplication,
    bench_modular_reduction,
    bench_complex_operations
);
criterion_main!(benches);
