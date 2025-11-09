use super::poly::{RnsBasis, RnsBasisBuilder, RnsError, RnsResult};
use std::sync::Arc;

/// Small pool of hand-picked NTT-friendly primes for toy parameter sets.
///
/// Each entry satisfies `p ≡ 1 (mod 2 * degree)` for at least one CKKS degree
/// used in tests/examples (degrees are powers of two).
const TOY_PRIME_POOL: [u64; 8] = [97, 193, 257, 769, 1153, 3329, 7681, 12289];

/// Build a minimal RNS basis with two CRT primes suitable for quick tests.
pub fn toy_basis<const DEGREE: usize>() -> RnsResult<Arc<RnsBasis>> {
    toy_basis_with_channels::<DEGREE>(2)
}

/// Build a toy RNS basis with the requested number of CRT channels.
pub fn toy_basis_with_channels<const DEGREE: usize>(
    channels: usize,
) -> RnsResult<Arc<RnsBasis>> {
    let requirement = 2 * DEGREE as u64;
    let mut selected = Vec::with_capacity(channels);

    for &candidate in TOY_PRIME_POOL.iter() {
        if (candidate - 1).is_multiple_of(requirement) {
            selected.push(candidate);
        }
        if selected.len() == channels {
            break;
        }
    }

    if selected.len() < channels {
        return Err(RnsError::MissingToyPrimes {
            degree: DEGREE,
            needed: channels,
            found: selected.len(),
        });
    }

    Ok(Arc::new(
        RnsBasisBuilder::new(DEGREE)
            .with_custom_primes(selected)
            .build()?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_basis_provides_two_channels() {
        const DEGREE: usize = 16;
        let basis = toy_basis::<DEGREE>().expect("toy basis");
        assert_eq!(basis.channel_count(), 2);
        for &prime in basis.primes() {
            assert_eq!((prime - 1) % (2 * DEGREE as u64), 0);
        }
    }

    #[test]
    fn toy_basis_with_custom_channels() {
        const DEGREE: usize = 32;
        let basis = toy_basis_with_channels::<DEGREE>(1)
            .expect("toy basis with one channel");
        assert_eq!(basis.channel_count(), 1);
    }
}
