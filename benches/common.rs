//! Shared bench harness: backends, seeded setup and bench ids.
//!
//! Every bench is generic over a [`Backend`], so a new ring backend is one
//! `impl` here plus one line per bench file. `BENCH_FULL=1` turns on the full
//! parameter matrix; the default run stays small.
// Each bench file uses a different subset of the helpers.
#![allow(dead_code, unused_macros)]

use std::sync::Arc;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::PolyRing;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

/// A ring backend the benches can build from a prime chain.
pub trait Backend<const N: usize>: PolyRing<N> {
    /// Short name used in bench ids.
    const NAME: &'static str;

    /// The top-level context of a new chain over `primes`.
    fn chain(primes: &[u64]) -> Self::Context;
}

impl<const N: usize> Backend<N> for RnsPoly<N> {
    const NAME: &'static str = "rns_ntt";

    fn chain(primes: &[u64]) -> Self::Context {
        Arc::new(RnsBasis::new(primes.to_vec()).expect("NTT-friendly primes"))
    }
}

/// `true` when `BENCH_FULL` is set to anything but `0`.
pub fn full() -> bool {
    std::env::var("BENCH_FULL").is_ok_and(|value| value != "0")
}

/// Criterion settings. The default run is a quick survey of ~130 benches:
/// short warm-up and measurement, ten samples, no plots. `BENCH_FULL`
/// restores Criterion's defaults (3 s + 5 s, 100 samples, HTML plots) for
/// precise comparisons; combine it with a filter to stay fast.
/// Command-line flags still override both.
pub fn config() -> Criterion {
    if full() {
        Criterion::default()
    } else {
        Criterion::default()
            .warm_up_time(Duration::from_millis(300))
            .measurement_time(Duration::from_secs(1))
            // Criterion's minimum: slow ops (decompose ~90 ms) need
            // ~n^2/2 iterations under linear sampling.
            .sample_size(10)
            .without_plots()
    }
}

/// `count` NTT-friendly primes of `bits` bits for degree `n`.
pub fn primes(bits: usize, count: usize, n: usize) -> Vec<u64> {
    generate_primes(bits, count, n as u64)
}

/// The context at `level` of the chain whose top context is `top`.
///
/// Derive every level from one `top`: a second `Backend::chain` call
/// builds a separate chain, and operations across chains are rejected.
pub fn context_at<P: Backend<N>, const N: usize>(
    top: &P::Context,
    level: usize,
) -> P::Context {
    P::zero(top)
        .mod_switch_to(level)
        .expect("level inside the chain")
        .context()
        .clone()
}

/// Seeded RNG, so every run benches the same inputs.
pub fn rng(seed: u64) -> ChaCha20Rng {
    ChaCha20Rng::seed_from_u64(seed)
}

/// Bench id `backend/N<n>/<bits>b/L<level>`, e.g. `rns_ntt/N1024/61b/L3`.
pub fn id<P: Backend<N>, const N: usize>(
    op: &str,
    bits: usize,
    level: usize,
) -> BenchmarkId {
    BenchmarkId::new(op, format!("{}/N{N}/{bits}b/L{level}", P::NAME))
}

/// Calls `$bench::<$backend<N>, N>($criterion)` for every listed degree.
/// Bench files load it with `#[macro_use] mod common;`.
macro_rules! for_degrees {
    ($criterion:expr, $bench:ident, $backend:ident, [$($degree:literal),* $(,)?]) => {
        $( $bench::<$backend<$degree>, $degree>($criterion); )*
    };
}
