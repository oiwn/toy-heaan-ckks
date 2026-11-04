//! Modular arithmetic primitives over NTT-friendly primes.
//!
//! Variants share bench ids (`<op>/<variant>/<bits>b`), so task 0052 adds
//! `shoup` and `barrett` next to the `u128` reference and Criterion compares
//! them directly.

mod common;

use std::hint::black_box;

use criterion::{
    BenchmarkId, Criterion, Throughput, criterion_group, criterion_main,
};
use rand::RngExt;
use toy_heaan_ckks::math::modular::{add_mod, mul_mod, pow_mod, sub_mod};

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
        let id = |op: &str| BenchmarkId::new(op, format!("u128/{bits}b"));
        group.bench_function(id("mul"), |bench| {
            bench.iter(|| mul_mod(black_box(x), black_box(y), black_box(q)));
        });
        group.bench_function(id("add"), |bench| {
            bench.iter(|| add_mod(black_box(x), black_box(y), black_box(q)));
        });
        group.bench_function(id("sub"), |bench| {
            bench.iter(|| sub_mod(black_box(x), black_box(y), black_box(q)));
        });
        group.bench_function(id("pow"), |bench| {
            bench.iter(|| pow_mod(black_box(x), black_box(q - 2), black_box(q)));
        });
    }
    group.finish();
}

fn batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("modmul_batch");
    group.throughput(Throughput::Elements(BATCH as u64));
    for bits in BITS {
        let (q, a, b) = operands(bits, BATCH);
        let mut out = vec![0u64; BATCH];
        let id = BenchmarkId::new("pointwise_mul", format!("u128/{bits}b"));
        group.bench_function(id, |bench| {
            bench.iter(|| {
                for ((o, x), y) in out.iter_mut().zip(&a).zip(&b) {
                    *o = mul_mod(*x, *y, q);
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
