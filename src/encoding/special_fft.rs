//! Helpers for the textbook CKKS embedding transforms. The final implementation
//! will mirror the FHE Textbook Vandermonde/special FFT flow while keeping the
//! logic reusable for future optimizations.

use std::f64::consts::PI;

use rustfft::num_complex::Complex64;

fn pow_mod(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut acc = 1u64;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = (acc * base) % modulus;
        }
        base = (base * base) % modulus;
        exp >>= 1;
    }
    acc
}

fn j_function(h: usize, degree: usize, is_plus: bool) -> u64 {
    let modulus = (2 * degree) as u64;
    let base = pow_mod(5, h as u64, modulus);
    if is_plus {
        base
    } else {
        (modulus - base) % modulus
    }
}

/// Accepted input views for the slot builder. Keeps encode callers ergonomic for
/// both purely real and complex inputs without forcing extra allocations.
#[derive(Clone, Copy, Debug)]
pub enum SlotInput<'a> {
    Real(&'a [f64]),
    Complex(&'a [Complex64]),
}

impl<'a> SlotInput<'a> {
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            SlotInput::Real(values) => values.len(),
            SlotInput::Complex(values) => values.len(),
        }
    }

    #[inline]
    pub fn get(&self, idx: usize) -> Complex64 {
        match self {
            SlotInput::Real(values) => Complex64::new(values[idx], 0.0),
            SlotInput::Complex(values) => values[idx],
        }
    }
}

impl<'a> From<&'a [f64]> for SlotInput<'a> {
    fn from(values: &'a [f64]) -> Self {
        SlotInput::Real(values)
    }
}

impl<'a> From<&'a [Complex64]> for SlotInput<'a> {
    fn from(values: &'a [Complex64]) -> Self {
        SlotInput::Complex(values)
    }
}

/// Cached Vandermonde data derived from the canonical CKKS roots of unity.
#[derive(Debug, Clone)]
pub struct VandermondeTables<const N: usize> {
    pub degree: usize,
    pub psi: Complex64,
    pub psi_inv: Complex64,
    pub omega: Complex64,
    pub omega_inv: Complex64,
    pub slot_roots: Vec<Complex64>,
    pub slot_roots_inv: Vec<Complex64>,
    pub psi_pows: Vec<Complex64>,
    pub psi_inv_pows: Vec<Complex64>,
}

impl<const N: usize> VandermondeTables<N> {
    pub fn new() -> Self {
        assert!(N.is_power_of_two(), "Degree must be a power of two");
        let angle_unit = PI / N as f64; // angle for primitive 2N-th root
        let psi = Complex64::from_polar(1.0, angle_unit);
        let psi_inv = psi.conj();
        let omega = psi * psi;
        let omega_inv = omega.conj();

        let mut slot_roots = Vec::with_capacity(N);
        let mut slot_roots_inv = Vec::with_capacity(N);
        let mut exponents = Vec::with_capacity(N);

        // Enumerate exponents via textbook J-function ordering.
        for h in 0..(N / 2) {
            exponents.push(j_function(h, N, true));
        }
        for h in (0..(N / 2)).rev() {
            exponents.push(j_function(h, N, false));
        }

        for exponent in exponents {
            let root = psi.powu(exponent as u32);
            slot_roots.push(root);
            slot_roots_inv.push(root.conj());
        }

        let mut psi_pows = Vec::with_capacity(N);
        let mut psi_inv_pows = Vec::with_capacity(N);
        let mut current = Complex64::new(1.0, 0.0);
        let mut current_inv = Complex64::new(1.0, 0.0);
        for _ in 0..N {
            psi_pows.push(current);
            psi_inv_pows.push(current_inv);
            current *= psi;
            current_inv *= psi_inv;
        }

        Self {
            degree: N,
            psi,
            psi_inv,
            omega,
            omega_inv,
            slot_roots,
            slot_roots_inv,
            psi_pows,
            psi_inv_pows,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.degree
    }
}

impl<const N: usize> Default for VandermondeTables<N> {
    fn default() -> Self {
        Self::new()
    }
}

/// Builds a conjugate-symmetric slot vector of length `N`, zero-padding the
/// remaining entries. Supports both real (`&[f64]`) and complex inputs.
pub fn build_conjugate_slots<const N: usize>(
    values: SlotInput<'_>,
) -> Vec<Complex64> {
    let mut slots = vec![Complex64::new(0.0, 0.0); N];
    let max_slots = N / 2;
    assert!(values.len() <= max_slots, "input exceeds slot capacity");

    for idx in 0..values.len() {
        let value = values.get(idx);
        slots[idx] = value;
        if idx > 0 {
            slots[N - idx] = value.conj();
        }
    }

    slots
}

#[inline]
pub fn build_slots_from_complex<const N: usize>(
    values: &[Complex64],
) -> Vec<Complex64> {
    build_conjugate_slots::<N>(SlotInput::from(values))
}

#[inline]
pub fn build_slots_from_real<const N: usize>(values: &[f64]) -> Vec<Complex64> {
    build_conjugate_slots::<N>(SlotInput::from(values))
}

/// Textbook inverse embedding (slot vector -> coefficient vector) using direct
/// Vandermonde evaluation with the conjugate roots of unity.
pub fn special_idft<const N: usize>(
    values: &[Complex64],
    tables: &VandermondeTables<N>,
) -> Vec<Complex64> {
    assert_eq!(values.len(), tables.len(), "expected N slot values");
    let mut permuted = vec![Complex64::new(0.0, 0.0); N];
    for (dst, src) in permuted.iter_mut().zip(values.iter().rev()) {
        *dst = *src;
    }
    let mut coeffs = vec![Complex64::new(0.0, 0.0); N];
    let inv_n = 1.0 / N as f64;

    for (slot_idx, value) in permuted.iter().enumerate() {
        let mut power = Complex64::new(1.0, 0.0);
        let root = tables.slot_roots[slot_idx];
        for coeff in coeffs.iter_mut() {
            *coeff += *value * power;
            power *= root;
        }
    }

    for coeff in coeffs.iter_mut() {
        *coeff *= inv_n;
    }

    coeffs
}

/// Textbook forward embedding (coefficients -> slot vector) using direct
/// Vandermonde evaluation.
pub fn special_dft<const N: usize>(
    coeffs: &[Complex64],
    tables: &VandermondeTables<N>,
) -> Vec<Complex64> {
    assert_eq!(coeffs.len(), tables.len(), "expected N coefficients");
    let mut slots = vec![Complex64::new(0.0, 0.0); N];

    for (slot_idx, slot) in slots.iter_mut().enumerate() {
        let mut power = Complex64::new(1.0, 0.0);
        let root = tables.slot_roots_inv[slot_idx];
        for coeff in coeffs {
            *slot += coeff * power;
            power *= root;
        }
    }

    slots.reverse();
    slots
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn conjugate_slots_build_symmetry() {
        const N: usize = 8;
        let input = vec![
            Complex64::new(1.0, 0.5),
            Complex64::new(-0.25, 0.75),
            Complex64::new(0.0, -1.0),
        ];
        let slots = build_slots_from_complex::<N>(&input);
        assert_eq!(slots.len(), N);
        assert_eq!(slots[0], input[0]);
        assert_eq!(slots[1], input[1]);
        assert_eq!(slots[2], input[2]);
        assert_eq!(slots[N - 1], input[1].conj());
        assert_eq!(slots[N - 2], input[2].conj());
    }

    #[test]
    fn conjugate_slots_real_inputs() {
        const N: usize = 8;
        let input = [1.5_f64, -2.25, 0.75];
        let slots = build_slots_from_real::<N>(&input);
        assert_eq!(slots[0], Complex64::new(1.5, 0.0));
        assert_eq!(slots[1], Complex64::new(-2.25, 0.0));
        assert_eq!(slots[2], Complex64::new(0.75, 0.0));
        assert_eq!(slots[N - 1], Complex64::new(-2.25, 0.0));
        assert_eq!(slots[N - 2], Complex64::new(0.75, 0.0));
    }

    #[test]
    #[should_panic(expected = "input exceeds slot capacity")]
    fn conjugate_slots_rejects_overflow() {
        const N: usize = 8;
        let input = vec![Complex64::new(0.0, 0.0); N];
        // Needs N/2 + 1 entries to trigger the panic.
        let _ = build_slots_from_complex::<N>(&input[..(N / 2 + 1)]);
    }

    #[test]
    fn vandermonde_table_consistency() {
        const N: usize = 8;
        let tables = VandermondeTables::<N>::new();
        let product = tables.psi * tables.psi;
        assert_relative_eq!(product.re, tables.omega.re, epsilon = 1e-12);
        assert_relative_eq!(product.im, tables.omega.im, epsilon = 1e-12);
        for idx in 0..N {
            let conj = tables.slot_roots[idx].conj();
            assert_relative_eq!(
                conj.re,
                tables.slot_roots_inv[idx].re,
                epsilon = 1e-12
            );
            assert_relative_eq!(
                conj.im,
                tables.slot_roots_inv[idx].im,
                epsilon = 1e-12
            );
            assert_relative_eq!(
                (tables.psi_pows[idx] * tables.psi_inv_pows[idx]).re,
                1.0,
                epsilon = 1e-12
            );
            assert_relative_eq!(
                (tables.psi_pows[idx] * tables.psi_inv_pows[idx]).im,
                0.0,
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn vandermonde_roundtrip() {
        const N: usize = 8;
        let tables = VandermondeTables::<N>::new();
        let coeffs: Vec<Complex64> = (0..N)
            .map(|idx| Complex64::new(idx as f64 / 7.0, -(idx as f64) / 11.0))
            .collect();
        let slots = special_dft(&coeffs, &tables);
        let recovered = special_idft(&slots, &tables);

        for (expected, actual) in coeffs.iter().zip(recovered.iter()) {
            assert_relative_eq!(expected.re, actual.re, epsilon = 1e-9);
            assert_relative_eq!(expected.im, actual.im, epsilon = 1e-9);
        }
    }
}
