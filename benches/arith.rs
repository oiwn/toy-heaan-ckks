//! Modular arithmetic primitives over NTT-friendly primes.
//!
//! Variants share bench ids (`<op>/<variant>/<bits>b`): `u128` is the
//! reference, `barrett` reduces general products and `shoup` multiplies by a
//! fixed operand with a precomputed companion (task 0052). Each variant is
//! checked against the reference once before timing.

mod common;

use std::hint::black_box;

use criterion::{
    BenchmarkId, Criterion, Throughput, criterion_group, criterion_main,
};
use rand::RngExt;
use toy_heaan_ckks::math::modular::{
    Barrett, add_mod, mul_mod, mul_shoup, pow_mod_reference, shoup_companion,
    sub_mod,
};

const BITS: [usize; 3] = [30, 50, 61];
const BATCH: usize = 4096;

/// One prime of `bits` bits and two reduced operand vectors of `len`.
fn operands(bits: usize, len: usize) -> (u64, Vec<u64>, Vec<u64>) {
    let modulus = common::primes(bits, 1, BATCH)[0];
    let mut rng = common::rng(bits as u64);
    let mut draw = || (0..len).map(|_| rng.random::<u64>() % modulus).collect();
    let a: Vec<u64> = draw();
    let b: Vec<u64> = draw();
    (modulus, a, b)
}

fn scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("modmul");
    for bits in BITS {
        let (q, a, b) = operands(bits, 1);
        let (x, y) = (a[0], b[0]);
        let barrett = Barrett::new(q);
        let y_companion = shoup_companion(y, q);
        assert_eq!(barrett.mul(x, y), mul_mod(x, y, q));
        assert_eq!(mul_shoup(x, y, y_companion, q), mul_mod(x, y, q));
        assert_eq!(barrett.pow(x, q - 2), pow_mod_reference(x, q - 2, q));

        let id = |op: &str, variant: &str| {
            BenchmarkId::new(op, format!("{variant}/{bits}b"))
        };
        group.bench_function(id("mul", "u128"), |bench| {
            bench.iter(|| mul_mod(black_box(x), black_box(y), black_box(q)));
        });
        group.bench_function(id("mul", "barrett"), |bench| {
            bench.iter(|| black_box(&barrett).mul(black_box(x), black_box(y)));
        });
        group.bench_function(id("mul", "shoup"), |bench| {
            bench.iter(|| {
                mul_shoup(
                    black_box(x),
                    black_box(y),
                    black_box(y_companion),
                    black_box(q),
                )
            });
        });
        group.bench_function(id("add", "u128"), |bench| {
            bench.iter(|| add_mod(black_box(x), black_box(y), black_box(q)));
        });
        group.bench_function(id("sub", "u128"), |bench| {
            bench.iter(|| sub_mod(black_box(x), black_box(y), black_box(q)));
        });
        group.bench_function(id("pow", "u128"), |bench| {
            bench.iter(|| {
                pow_mod_reference(black_box(x), black_box(q - 2), black_box(q))
            });
        });
        group.bench_function(id("pow", "barrett"), |bench| {
            bench.iter(|| black_box(&barrett).pow(black_box(x), black_box(q - 2)));
        });
    }
    group.finish();
}

fn batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("modmul_batch");
    group.throughput(Throughput::Elements(BATCH as u64));
    for bits in BITS {
        let (q, a, b) = operands(bits, BATCH);
        let barrett = Barrett::new(q);
        // A fixed multiplier, as for NTT twiddles and rescale constants.
        let (w, w_companion) = (b[0], shoup_companion(b[0], q));
        let mut out = vec![0u64; BATCH];
        for (x, y) in a.iter().zip(&b) {
            assert_eq!(barrett.mul(*x, *y), mul_mod(*x, *y, q));
            assert_eq!(mul_shoup(*x, w, w_companion, q), mul_mod(*x, w, q));
        }

        let id = |op: &str, variant: &str| {
            BenchmarkId::new(op, format!("{variant}/{bits}b"))
        };
        group.bench_function(id("pointwise_mul", "u128"), |bench| {
            bench.iter(|| {
                for ((o, x), y) in out.iter_mut().zip(&a).zip(&b) {
                    *o = mul_mod(*x, *y, q);
                }
                black_box(&out);
            });
        });
        group.bench_function(id("pointwise_mul", "barrett"), |bench| {
            bench.iter(|| {
                for ((o, x), y) in out.iter_mut().zip(&a).zip(&b) {
                    *o = barrett.mul(*x, *y);
                }
                black_box(&out);
            });
        });
        group.bench_function(id("pointwise_mul_const", "u128"), |bench| {
            bench.iter(|| {
                for (o, x) in out.iter_mut().zip(&a) {
                    *o = mul_mod(*x, w, q);
                }
                black_box(&out);
            });
        });
        group.bench_function(id("pointwise_mul_const", "shoup"), |bench| {
            bench.iter(|| {
                for (o, x) in out.iter_mut().zip(&a) {
                    *o = mul_shoup(*x, w, w_companion, q);
                }
                black_box(&out);
            });
        });
    }
    group.finish();
}

criterion_group! {
    name = arith;
    config = common::config();
    targets = scalar, batch
}
criterion_main!(arith);
