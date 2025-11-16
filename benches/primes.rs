use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use toy_heaan_ckks::math::{generate_primes, is_prime, is_prime_reference};

fn bench_prime_checks(c: &mut Criterion) {
    let inputs = [2_147_483_647u64, 9_223_372_036_854_775_783u64];
    let mut group = c.benchmark_group("is_prime");

    for &n in &inputs {
        group.bench_with_input(BenchmarkId::new("miller_rabin", n), &n, |b, &n| {
            b.iter(|| black_box(is_prime(black_box(n))));
        });
        group.bench_with_input(BenchmarkId::new("reference", n), &n, |b, &n| {
            b.iter(|| black_box(is_prime_reference(black_box(n))));
        });
    }

    group.finish();
}

fn bench_prime_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("generate_primes");
    let configs = [(30usize, 5usize, 8192u64), (61, 5, 8192)];

    for &(bits, count, degree) in &configs {
        let label = format!("{bits}b_{count}x_deg{degree}");
        group.bench_with_input(
            BenchmarkId::from_parameter(label),
            &(bits, count, degree),
            |b, &(bits, count, degree)| {
                b.iter(|| {
                    black_box(generate_primes(bits, count, degree));
                });
            },
        );
    }

    group.finish();
}

criterion_group!(primes, bench_prime_checks, bench_prime_generation);
criterion_main!(primes);
