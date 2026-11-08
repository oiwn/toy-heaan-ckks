//! Ring operations (`PolyRing`) per backend, degree, prime width and level,
//! plus the RNS NTT kernel on one channel.

#[macro_use]
mod common;

use std::hint::black_box;

use common::Backend;
use criterion::{
    BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main,
};
use rand::RngExt;
use toy_heaan_ckks::ckks::keys::DIGIT_BITS;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

/// Every ring operation at a low and a high level of a 61-bit chain.
fn ring_ops<P: Backend<N>, const N: usize>(c: &mut Criterion) {
    let levels: [usize; 2] = if common::full() { [3, 20] } else { [3, 12] };
    let top = P::chain(&common::primes(61, levels[1], N));
    let bottom = common::context_at::<P, N>(&top, 1);
    let mut rng = common::rng(N as u64);
    let mut group = c.benchmark_group("ring");

    for level in levels {
        let context = common::context_at::<P, N>(&top, level);
        let a = P::sample_uniform(&context, &mut rng);
        let b = P::sample_uniform(&context, &mut rng);
        let low = P::sample_uniform(&bottom, &mut rng);
        // `to_coeffs` returns i64, so the conversion benches use a message-
        // sized polynomial (|c| < 2^40), as encode/decode do.
        let coeffs: Vec<i64> = (0..N)
            .map(|_| (rng.random::<u64>() >> 23) as i64 - (1 << 40))
            .collect();
        let small = P::from_coeffs(&coeffs, &context);
        assert_eq!(
            small.to_coeffs().as_slice(),
            coeffs.as_slice(),
            "coeff round trip"
        );
        let id = |op: &str| common::id::<P, N>(op, 61, level);

        group.bench_function(id("add"), |bench| {
            bench.iter_batched(
                || a.clone(),
                |mut x| {
                    x += &b;
                    x
                },
                BatchSize::LargeInput,
            );
        });
        group.bench_function(id("mul"), |bench| {
            bench.iter_batched(
                || a.clone(),
                |mut x| {
                    x *= &b;
                    x
                },
                BatchSize::LargeInput,
            );
        });
        group.bench_function(id("neg"), |bench| {
            bench.iter_batched(|| a.clone(), |x| -x, BatchSize::LargeInput);
        });
        group.bench_function(id("automorphism"), |bench| {
            bench.iter(|| black_box(a.automorphism(5)));
        });
        group.bench_function(id("decompose"), |bench| {
            bench.iter(|| black_box(a.decompose(DIGIT_BITS)));
        });
        group.bench_function(id("gadget_term"), |bench| {
            bench.iter(|| black_box(a.gadget_term(0, DIGIT_BITS)));
        });
        group.bench_function(id("rescale"), |bench| {
            bench.iter(|| black_box(a.rescale().expect("level >= 2")));
        });
        group.bench_function(id("mod_switch"), |bench| {
            bench.iter(|| black_box(a.mod_switch_to(1).expect("level 1")));
        });
        group.bench_function(id("mod_raise"), |bench| {
            bench.iter(|| black_box(low.mod_raise(&context).expect("same chain")));
        });
        group.bench_function(id("from_coeffs"), |bench| {
            bench.iter(|| black_box(P::from_coeffs(&coeffs, &context)));
        });
        group.bench_function(id("to_coeffs"), |bench| {
            bench.iter(|| black_box(small.to_coeffs()));
        });
    }
    group.finish();
}

/// `mul` and `rescale` with 50- and 61-bit primes at the same level.
fn prime_width<P: Backend<N>, const N: usize>(c: &mut Criterion) {
    let level = 3;
    let mut rng = common::rng(N as u64);
    let mut group = c.benchmark_group("prime_width");
    for bits in [50, 61] {
        let context = P::chain(&common::primes(bits, level, N));
        let a = P::sample_uniform(&context, &mut rng);
        let b = P::sample_uniform(&context, &mut rng);
        group.bench_function(common::id::<P, N>("mul", bits, level), |bench| {
            bench.iter_batched(
                || a.clone(),
                |mut x| {
                    x *= &b;
                    x
                },
                BatchSize::LargeInput,
            );
        });
        group.bench_function(common::id::<P, N>("rescale", bits, level), |bench| {
            bench.iter(|| black_box(a.rescale().expect("level >= 2")));
        });
    }
    group.finish();
}

/// Forward and inverse negacyclic NTT on one 61-bit RNS channel.
fn ntt<const N: usize>(c: &mut Criterion) {
    let basis = RnsBasis::<N>::new(common::primes(61, 1, N)).expect("primes");
    let table = basis.ntt_table(0);
    let mut rng = common::rng(N as u64);
    let values: Vec<u64> = (0..N)
        .map(|_| rng.random::<u64>() % table.modulus)
        .collect();
    let input: Box<[u64; N]> =
        values.into_boxed_slice().try_into().expect("N values");

    // Correctness once, outside timing: the transforms round trip.
    let mut check = input.clone();
    table.forward(&mut check);
    table.inverse(&mut check);
    assert_eq!(check, input, "NTT round trip");

    let mut group = c.benchmark_group("ntt");
    let id = |op: &str| BenchmarkId::new(op, format!("rns_ntt/N{N}/61b"));
    group.bench_function(id("forward"), |bench| {
        bench.iter_batched_ref(
            || input.clone(),
            |v| table.forward(v),
            BatchSize::LargeInput,
        );
    });
    group.bench_function(id("inverse"), |bench| {
        bench.iter_batched_ref(
            || input.clone(),
            |v| table.inverse(v),
            BatchSize::LargeInput,
        );
    });
    group.finish();
}

fn benches(c: &mut Criterion) {
    ntt::<1024>(c);
    ntt::<2048>(c);
    ntt::<4096>(c);
    for_degrees!(c, ring_ops, RnsPoly, [1024, 2048, 4096]);
    for_degrees!(c, prime_width, RnsPoly, [1024, 4096]);
    if common::full() {
        ntt::<8192>(c);
        ntt::<16384>(c);
        for_degrees!(c, ring_ops, RnsPoly, [8192, 16384]);
        for_degrees!(c, prime_width, RnsPoly, [16384]);
    }
}

criterion_group! {
    name = ring;
    config = common::config();
    targets = benches
}
criterion_main!(ring);
