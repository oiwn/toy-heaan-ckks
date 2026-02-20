#!/usr/bin/env python3
"""
Reference CKKS encoding/decoding implementation based on FHE textbook.

This serves as ground truth for validating the Rust implementation.
Uses explicit matrix operations to match the textbook algorithm exactly.

References:
- https://fhetextbook.github.io/EncodingandDecoding.html
- https://github.com/fhetextbook/fhetextbook.github.io/blob/main/source%20code/ckks.py
"""

import numpy as np
from typing import List, Tuple

# Parameters for toy example
N = 8  # Polynomial degree (must be power of 2)
M = 2 * N  # Order for primitive root
SCALE_BITS = 30
DELTA = 2**SCALE_BITS  # Scaling factor

print(f"=== CKKS Encoding Reference (N={N}) ===\n")


def j_function(h: int, degree: int, is_plus: bool) -> int:
    """
    J-function for root ordering: J(h) = 5^h mod 2N

    Args:
        h: index (0 to N/2-1)
        degree: polynomial degree N
        is_plus: True for J(h), False for -J(h) mod 2N

    Returns:
        Exponent for primitive 2N-th root
    """
    modulus = 2 * degree
    base = pow(5, h, modulus)
    if is_plus:
        return base
    else:
        return (modulus - base) % modulus


def build_vandermonde_matrix(degree: int) -> np.ndarray:
    """
    Build the Vandermonde matrix W for CKKS canonical embedding.

    Row i corresponds to root ψ^{exponent[i]}
    Column j contains powers: [1, root^1, root^2, ..., root^{N-1}]

    Returns:
        N×N complex matrix W
    """
    psi = np.exp(1j * np.pi / degree)  # Primitive 2N-th root

    # Compute exponents using J-function ordering
    exponents = []
    for h in range(degree // 2):
        exponents.append(j_function(h, degree, True))  # First half: J(h)
    for h in range(degree // 2 - 1, -1, -1):
        exponents.append(j_function(h, degree, False))  # Second half: -J(h) reversed

    print(f"Root exponents (J-function ordering): {exponents}")

    # Build Vandermonde matrix
    W = np.zeros((degree, degree), dtype=complex)
    for row_idx, exp in enumerate(exponents):
        root = psi ** exp
        for col_idx in range(degree):
            W[row_idx, col_idx] = root ** col_idx

    return W


def build_conjugate_slots(values: List[complex]) -> np.ndarray:
    """
    Build conjugate-symmetric slot vector from N/2 input values.

    Implements the textbook pattern: [z0, z1, ..., z_{N/2-1}, conj(z0), ..., conj(z_{N/2-1})]

    Args:
        values: List of N/2 complex numbers

    Returns:
        N-length array with conjugate symmetry
    """
    assert len(values) <= N // 2, f"Too many values: {len(values)} > {N//2}"

    # Pad to N/2 if needed
    padded = list(values) + [0.0+0.0j] * (N // 2 - len(values))

    # Append conjugates sequentially (textbook pattern)
    result = np.array(padded + [np.conj(z) for z in padded], dtype=complex)

    print(f"\nConjugate-symmetric slots ({len(values)} input values):")
    for i, val in enumerate(result):
        print(f"  slot[{i}] = {val.real:+.6f} {val.imag:+.6f}i")

    return result


def encode(values: List[complex], W: np.ndarray, scale: float = DELTA) -> Tuple[np.ndarray, np.ndarray]:
    """
    Encode slot vector into polynomial coefficients.

    Algorithm:
        1. Build conjugate-symmetric vector (N values)
        2. Apply anti-diagonal flip: reverse the vector
        3. Multiply by Vandermonde transpose: (1/N) * W^T @ flipped
        4. Scale by Δ and round to integers

    Args:
        values: List of N/2 complex numbers
        W: N×N Vandermonde matrix
        scale: Scaling factor Δ

    Returns:
        (unscaled_coeffs, scaled_int_coeffs)
    """
    # Step 1: Build conjugate-symmetric slots
    slots = build_conjugate_slots(values)

    # Step 2: Apply anti-diagonal flip (reverse)
    anti_diag = np.eye(N)[::-1]  # Anti-diagonal identity matrix
    flipped = anti_diag @ slots
    print(f"\nAfter anti-diagonal flip:")
    for i, val in enumerate(flipped[:4]):  # Print first 4
        print(f"  flipped[{i}] = {val.real:+.6f} {val.imag:+.6f}i")
    print("  ...")

    # Step 3: Vandermonde transform with 1/N scaling
    # CRITICAL: Use W^T (transpose), not W!
    W_T = W.T
    coeffs_unscaled = (1.0 / N) * (W_T @ flipped)

    print(f"\nCoefficients before scaling:")
    max_imag = max(abs(c.imag) for c in coeffs_unscaled)
    print(f"  Max imaginary part: {max_imag:.2e}")
    for i, c in enumerate(coeffs_unscaled[:4]):
        print(f"  coeff[{i}] = {c.real:+.6f} {c.imag:+.6f}i")
    print("  ...")

    # Check if coefficients are real (within tolerance)
    IMAG_TOL = 1e-6
    if max_imag > IMAG_TOL:
        print(f"  ⚠ WARNING: Imaginary parts exceed tolerance {IMAG_TOL}")
    else:
        print(f"  ✓ Coefficients are nearly real (max imag: {max_imag:.2e})")

    # Step 4: Scale and round
    coeffs_scaled = coeffs_unscaled * scale
    coeffs_int = np.round(coeffs_scaled.real).astype(int)

    print(f"\nScaled integer coefficients (Δ = 2^{SCALE_BITS}):")
    for i, c in enumerate(coeffs_int[:4]):
        print(f"  coeff_int[{i}] = {c}")
    print("  ...")

    return coeffs_unscaled, coeffs_int


def decode(coeffs_int: np.ndarray, W: np.ndarray, scale: float = DELTA) -> np.ndarray:
    """
    Decode polynomial coefficients back to slot vector.

    Algorithm:
        1. Unscale: divide by Δ
        2. Multiply by W (NOT transpose): evaluate polynomial at roots
        3. Reverse the result (undo anti-diagonal)
        4. Extract first N/2 values

    Args:
        coeffs_int: N integer coefficients
        W: N×N Vandermonde matrix
        scale: Scaling factor Δ

    Returns:
        Decoded slot vector (first N/2 values)
    """
    # Step 1: Unscale and convert to complex
    coeffs_complex = coeffs_int.astype(complex) / scale

    # Step 2: Multiply by W (NOT W^T!)
    # Since encoding uses W^T, decoding uses W
    slots_full = W @ coeffs_complex

    # Step 3: Reverse (undo anti-diagonal)
    slots_unreversed = slots_full[::-1]

    print(f"\nDecoded slots (full):")
    for i, val in enumerate(slots_unreversed[:4]):
        print(f"  slot[{i}] = {val.real:+.6f} {val.imag:+.6f}i")
    print("  ...")

    # Step 4: Extract first N/2 values (original slots)
    slots_decoded = slots_unreversed[:N//2]

    return slots_decoded


def test_roundtrip():
    """Test encode-decode roundtrip with simple inputs."""
    print("\n" + "="*60)
    print("TEST: Encode-Decode Roundtrip")
    print("="*60)

    # Test input: a few complex values
    test_values = [
        1.0 + 0.0j,
        -0.5 + 0.25j,
        0.75 + 0.0j,
        0.0 + 0.5j
    ]

    print(f"\nInput slots ({len(test_values)} values):")
    for i, val in enumerate(test_values):
        print(f"  input[{i}] = {val.real:+.6f} {val.imag:+.6f}i")

    # Build Vandermonde matrix
    W = build_vandermonde_matrix(N)

    print(f"\nVandermonde matrix W (first 4 rows, 4 cols):")
    for i in range(min(4, N)):
        row_str = "  "
        for j in range(min(4, N)):
            row_str += f"{W[i,j].real:+.3f}{W[i,j].imag:+.3f}i  "
        print(row_str)
    print("  ...")

    # Encode
    print("\n" + "-"*60)
    print("ENCODING")
    print("-"*60)
    coeffs_unscaled, coeffs_int = encode(test_values, W)

    # Decode
    print("\n" + "-"*60)
    print("DECODING")
    print("-"*60)
    decoded = decode(coeffs_int, W)

    # Compare
    print("\n" + "-"*60)
    print("ROUNDTRIP ERROR")
    print("-"*60)
    print(f"\nComparison (input vs decoded):")
    errors = []
    for i, (orig, dec) in enumerate(zip(test_values, decoded)):
        re_err = abs(orig.real - dec.real)
        im_err = abs(orig.imag - dec.imag)
        max_err = max(re_err, im_err)
        errors.append(max_err)
        print(f"  slot[{i}]:")
        print(f"    Input:   {orig.real:+.6f} {orig.imag:+.6f}i")
        print(f"    Decoded: {dec.real:+.6f} {dec.imag:+.6f}i")
        print(f"    Error:   {max_err:.2e}")

    max_error = max(errors)
    print(f"\nMax absolute error: {max_error:.2e}")

    if max_error < 1e-6:
        print("✓ PASS: Roundtrip error within tolerance")
    else:
        print(f"✗ FAIL: Roundtrip error {max_error:.2e} exceeds 1e-6")

    return max_error < 1e-6


def verify_vandermonde_orthogonality():
    """Verify that W^T @ W = N * I (orthogonality property)."""
    print("\n" + "="*60)
    print("VERIFICATION: Vandermonde Orthogonality")
    print("="*60)

    W = build_vandermonde_matrix(N)
    W_T = W.T.conj()
    product = W_T @ W

    # Should equal N * Identity
    expected = N * np.eye(N, dtype=complex)
    diff = product - expected
    max_diff = np.max(np.abs(diff))

    print(f"\nW^T @ W - N*I max error: {max_diff:.2e}")

    if max_diff < 1e-10:
        print("✓ PASS: Vandermonde matrix is orthogonal")
        return True
    else:
        print(f"✗ FAIL: Orthogonality error {max_diff:.2e} too large")
        return False


def main():
    """Run all tests and print summary."""
    print("\nRunning CKKS encoding reference implementation...\n")

    # Verify Vandermonde properties
    ortho_ok = verify_vandermonde_orthogonality()

    # Test roundtrip
    roundtrip_ok = test_roundtrip()

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Vandermonde orthogonality: {'PASS' if ortho_ok else 'FAIL'}")
    print(f"Encode-decode roundtrip:   {'PASS' if roundtrip_ok else 'FAIL'}")

    if ortho_ok and roundtrip_ok:
        print("\n✓ All checks passed!")
        print("\nThis reference can be used to validate the Rust implementation.")
    else:
        print("\n✗ Some checks failed - review algorithm!")


if __name__ == "__main__":
    main()
