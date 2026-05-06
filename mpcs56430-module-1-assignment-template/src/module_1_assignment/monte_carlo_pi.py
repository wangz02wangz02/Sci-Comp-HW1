"""Monte Carlo estimation of pi via random points in the unit square."""

from __future__ import annotations

import random
from typing import List, Tuple


Point = Tuple[float, float]


def simulate_pi(num_points: int, seed: int | None = None) -> tuple[float, List[Point], int]:
    """Throw ``num_points`` darts at the unit square and estimate pi.

    Returns ``(pi_estimate, points, hits)`` where ``points`` is the list of
    sampled (x, y) coordinates and ``hits`` is the count that fell inside the
    quarter circle of radius 1.
    """
    rng = random.Random(seed)
    points: list[Point] = []
    hits = 0

    for _ in range(num_points):
        x = rng.uniform(0, 1)
        y = rng.uniform(0, 1)
        points.append((x, y))
        if x * x + y * y <= 1:
            hits += 1

    pi_estimate = 4 * hits / num_points if num_points else 0.0
    return pi_estimate, points, hits


def plot_pi_simulation(points: List[Point]) -> None:
    """Scatter the sampled points, coloring those inside the unit circle red."""
    import matplotlib.pyplot as plt

    inside_x: list[float] = []
    inside_y: list[float] = []
    outside_x: list[float] = []
    outside_y: list[float] = []

    for x, y in points:
        if x * x + y * y <= 1:
            inside_x.append(x)
            inside_y.append(y)
        else:
            outside_x.append(x)
            outside_y.append(y)

    plt.figure(figsize=(6, 6))
    plt.scatter(inside_x, inside_y, color="red", s=20)
    plt.scatter(outside_x, outside_y, color="blue", s=20)
    plt.xlim(0, 1.05)
    plt.ylim(0, 1.05)
    plt.gca().set_aspect("equal")
    plt.show()
