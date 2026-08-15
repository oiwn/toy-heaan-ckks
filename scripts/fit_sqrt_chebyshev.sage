# Fit polynomial approximations for sqrt(x) on variance intervals.
#
# Usage:
#   sage scripts/fit_sqrt_chebyshev.sage
#
# If Sage tries to write under ~/.sage and your environment blocks it, run:
#   SAGE_DOT_SAGE=/private/tmp/sage-dot-sage sage scripts/fit_sqrt_chebyshev.sage

R = RealField(120)


def fit_cheb_interpolant(lo, hi, degree):
    """Return power-basis coefficients for a Chebyshev-node interpolant."""
    point_count = degree + 1
    xs = []
    ys = []
    pi_value = R.pi()

    for idx in range(point_count):
        node = R(cos((2 * idx + 1) * pi_value / (2 * point_count)))
        x_value = (R(lo) + R(hi)) / 2 + (R(hi) - R(lo)) / 2 * node
        xs.append(x_value)
        ys.append(sqrt(x_value))

    system = matrix(
        R,
        [[x_value**power for power in range(point_count)] for x_value in xs],
    )
    return list(system.solve_right(vector(R, ys)))


def max_grid_error(coeffs, lo, hi, sample_count=20001):
    """Estimate max absolute approximation error on a dense grid."""
    max_error = R(0)
    max_x = R(lo)

    for idx in range(sample_count):
        x_value = R(lo) + (R(hi) - R(lo)) * R(idx) / R(sample_count - 1)
        approx = sum(
            coeff * x_value**power for power, coeff in enumerate(coeffs)
        )
        error = abs(approx - sqrt(x_value))
        if error > max_error:
            max_error = error
            max_x = x_value

    return max_error, max_x


def print_fit(lo, hi, degree):
    coeffs = fit_cheb_interpolant(lo, hi, degree)
    max_error, max_x = max_grid_error(coeffs, lo, hi)

    print("interval:", n(R(lo), digits=12), n(R(hi), digits=12))
    print("degree:", degree)
    print("power basis coefficients:")
    for idx, coeff in enumerate(coeffs):
        print(f"  a{idx} = {n(coeff, digits=20)}")
    print("max grid error:", n(max_error, digits=14))
    print("max error at x:", n(max_x, digits=14))
    print()


def print_rust_coeff_array(lo, hi, degree):
    coeffs = fit_cheb_interpolant(lo, hi, degree)
    max_error, max_x = max_grid_error(coeffs, lo, hi)

    print(
        "Rust coefficients for sqrt(x), interval ["
        + str(n(R(lo), digits=12))
        + ", "
        + str(n(R(hi), digits=12))
        + "], degree "
        + str(degree)
    )
    print("const SQRT_COEFFS: [f64; " + str(degree + 1) + "] = [")
    for coeff in coeffs:
        print("    " + str(n(coeff, digits=20)) + ",")
    print("];")
    print("max grid error:", n(max_error, digits=14))
    print("max error at x:", n(max_x, digits=14))
    print()


def print_degree_sweep(lo, hi, degrees):
    print("interval:", n(R(lo), digits=12), n(R(hi), digits=12))
    print("degree sweep:")
    print("  degree | horner_muls | balanced_depth_est | max_grid_error | at_x")
    for degree in degrees:
        coeffs = fit_cheb_interpolant(lo, hi, degree)
        max_error, max_x = max_grid_error(coeffs, lo, hi)
        if degree <= 1:
            balanced_depth = degree
        else:
            balanced_depth = ceil(log(R(degree), R(2)))
        print(
            "  "
            + str(degree).rjust(6)
            + " | "
            + str(degree).rjust(11)
            + " | "
            + str(balanced_depth).rjust(18)
            + " | "
            + str(n(max_error, digits=12)).rjust(14)
            + " | "
            + str(n(max_x, digits=12))
        )
    print()


HI = R("0.25")
DEFAULT_DEGREE = 3
DEGREES = [3, 5, 7, 9, 11]
LOW_BOUNDS = [
    R("0"),
    R("0.0001"),
    R("0.001"),
    R("0.005"),
    R("0.01"),
]

print("Detailed coefficients for the default degree")
print("============================================")
print()

for low in LOW_BOUNDS:
    print_fit(low, HI, DEFAULT_DEGREE)

print("Degree/error sweep")
print("==================")
print()

for low in LOW_BOUNDS:
    print_degree_sweep(low, HI, DEGREES)

print("Recommended demo coefficients")
print("=============================")
print()

print_rust_coeff_array(R("0.01"), HI, 7)
print_rust_coeff_array(R("0.01"), HI, 9)
