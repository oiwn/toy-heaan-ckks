"""Plaintext EvalMod experiment; no encrypted-accuracy claim. Run with Python 3."""

import math


def fit(function, degree):
    count = degree + 1
    angles = [math.pi * (j + 0.5) / count for j in range(count)]
    values = [function(math.cos(angle)) for angle in angles]
    return [
        (1 if k == 0 else 2)
        * sum(value * math.cos(k * angle) for value, angle in zip(values, angles))
        / count
        for k in range(count)
    ]


def evaluate(coefficients, value):
    next_term = previous = 0.0
    for coefficient in reversed(coefficients[1:]):
        term = 2 * value * next_term - previous + coefficient
        previous, next_term = next_term, term
    return value * next_term - previous + coefficients[0]


for degree_ring in (16, 64):
    weight = degree_ring // 2
    bound = (weight + 2) // 2
    radius = bound + 1
    ratio = 2.0 ** -21
    for degree in (15, 31, 63):
        for doubling in (3, 4, 5):
            coefficients = fit(
                lambda t: math.cos(2 * math.pi * radius * t / 2**doubling),
                degree,
            )
            worst = 0.0
            for integer in range(-bound, bound + 1):
                for sample in range(101):
                    message = -1 + 2 * sample / 100
                    value = integer + ratio * message
                    result = evaluate(coefficients, (value - 0.25) / radius)
                    for _ in range(doubling):
                        result = 2 * result * result - 1
                    restored = result / (2 * math.pi * ratio)
                    worst = max(worst, abs(restored - message))
            print(
                f"N={degree_ring} degree={degree} double={doubling} "
                f"wrap_bound={bound} scalar_error={worst:.6g} "
                f"slot_bound={degree_ring * worst:.6g}"
            )
