//! Slot aggregation: sum packed slots with a binary rotate-and-add tree.

use crate::PolyRing;
use crate::ckks::evaluator::RotationKeyMap;
use crate::ckks::types::Ciphertext;

/// Rotation offsets needed by [`sum_slots`] for `slots` packed values:
/// `1, 2, 4, ..., slots/2`.
///
/// # Panics
/// Panics if `slots` is not a power-of-two >= 2.
#[must_use]
pub fn required_rotation_offsets(slots: usize) -> Vec<i32> {
    assert!(slots.is_power_of_two(), "slots must be a power of two");
    assert!(slots >= 2, "summing fewer than two slots is a no-op");

    (0..slots.ilog2()).map(|shift| 1i32 << shift).collect()
}

/// Sums windows of `slots` consecutive slots via a binary rotation tree.
///
/// Slot `i` of the result holds `x_i + ... + x_{i+slots-1}` (indices mod
/// `DEGREE / 2`, since rotations wrap over all CKKS slots). So
/// `Enc([x0, ..., x_{slots-1}])` becomes `Enc([sum, ..., sum])` in every
/// slot only if `slots == DEGREE / 2` or the input is tiled (slot `i` holds
/// `x_{i mod slots}`). With unfilled slots left at zero, only slot 0 holds
/// the full sum. Consumes no levels; the level and basis are unchanged.
///
/// # Panics
/// Panics if `slots` is not a power of two, exceeds `DEGREE / 2`, or a
/// required rotation key is missing from `rotkeys`.
pub fn sum_slots<P: PolyRing<DEGREE>, const DEGREE: usize>(
    ct: &Ciphertext<P, DEGREE>,
    slots: usize,
    rotkeys: &RotationKeyMap<P, DEGREE>,
) -> Ciphertext<P, DEGREE> {
    assert!(
        slots <= DEGREE / 2,
        "slots cannot exceed available CKKS slots"
    );

    let mut acc = ct.clone();

    for offset in required_rotation_offsets(slots) {
        let rotkey = rotkeys
            .get(&offset)
            .unwrap_or_else(|| panic!("missing rotation key for offset {offset}"));
        acc = &acc + &acc.rotate(rotkey);
    }

    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn offsets_are_powers_of_two_below_the_window() {
        assert_eq!(required_rotation_offsets(2), vec![1]);
        assert_eq!(required_rotation_offsets(8), vec![1, 2, 4]);
    }
}
