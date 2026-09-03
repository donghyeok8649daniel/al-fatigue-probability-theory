"""Verify that identical one-point P0/F0 do not uniquely determine later P.

Two finite generalized-LJ chains are initialized from permutations of the same
spacing/rate pairs. Therefore their one-point empirical P0(lambda) and
F0(lambda,c) are exactly identical. Their ordered neighbor structure differs,
so the deterministic accelerations and later spacing distributions differ.

This is a reduced-information audit, not a stochastic model and not a new
physical probability assumption.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


M_EXP = 12.19
N_EXP = 6.0
DT = 1.0e-4
END_TIME = 1.0

SPACING_A = np.array([0.98, 0.99, 1.00, 1.01, 1.02, 1.005, 0.995, 1.015])
PERMUTATION = np.array([0, 4, 2, 6, 1, 7, 3, 5])
SPACING_B = SPACING_A[PERMUTATION]
RATE_A = np.zeros_like(SPACING_A)
RATE_B = np.zeros_like(SPACING_B)


def normalized_lj_force(stretch: np.ndarray) -> np.ndarray:
    return (
        stretch ** (-N_EXP - 1.0) - stretch ** (-M_EXP - 1.0)
    ) / (M_EXP - N_EXP)


def spacing_acceleration(spacing: np.ndarray, q: float = 0.0) -> np.ndarray:
    """Spacing acceleration for the active fixed-left / force-right chain."""
    dphi = normalized_lj_force(spacing)
    out = np.empty_like(spacing)
    out[0] = dphi[1] - dphi[0]
    out[1:-1] = dphi[2:] - 2.0 * dphi[1:-1] + dphi[:-2]
    out[-1] = q + dphi[-2] - 2.0 * dphi[-1]
    return out


def node_force(x: np.ndarray, q: float = 0.0) -> np.ndarray:
    spacing = np.diff(x)
    dphi = normalized_lj_force(spacing)
    force = np.zeros_like(x)
    force[1:-1] = dphi[1:] - dphi[:-1]
    force[-1] = -dphi[-1] + q
    return force


def state_from_spacing_and_rate(spacing: np.ndarray, rate: np.ndarray):
    x = np.concatenate(([0.0], np.cumsum(spacing)))
    velocity = np.concatenate(([0.0], np.cumsum(rate)))
    return x, velocity


def simulate(spacing0: np.ndarray, rate0: np.ndarray):
    x, velocity = state_from_spacing_and_rate(spacing0, rate0)
    force = node_force(x)
    n_steps = int(round(END_TIME / DT))

    for _ in range(n_steps):
        velocity[1:] += 0.5 * DT * force[1:]
        x[1:] += DT * velocity[1:]
        new_force = node_force(x)
        velocity[1:] += 0.5 * DT * new_force[1:]
        force = new_force

    return np.diff(x), np.diff(velocity)


def empirical_ks(sample_a: np.ndarray, sample_b: np.ndarray) -> float:
    grid = np.unique(np.concatenate((sample_a, sample_b)))
    cdf_a = np.array([np.mean(sample_a <= value) for value in grid])
    cdf_b = np.array([np.mean(sample_b <= value) for value in grid])
    return float(np.max(np.abs(cdf_a - cdf_b)))


def second_moment_curvature(spacing: np.ndarray, rate: np.ndarray) -> float:
    # For g(lambda)=lambda^2:
    # d2 <g>/dt2 = mean(2*c^2 + 2*lambda*lambda_ddot).
    acceleration = spacing_acceleration(spacing)
    return float(np.mean(2.0 * rate**2 + 2.0 * spacing * acceleration))


def main() -> None:
    assert np.array_equal(np.sort(SPACING_A), np.sort(SPACING_B))
    assert np.array_equal(np.sort(RATE_A), np.sort(RATE_B))

    final_a, final_rate_a = simulate(SPACING_A, RATE_A)
    final_b, final_rate_b = simulate(SPACING_B, RATE_B)

    result = {
        "classification": "deterministic finite-chain reduced-initial-data counterexample",
        "model": {
            "m": M_EXP,
            "n": N_EXP,
            "dt": DT,
            "end_time": END_TIME,
            "external_force": 0.0,
            "number_of_spacings": int(SPACING_A.size),
        },
        "initial_data": {
            "spacing_multiset_identical": True,
            "rate_multiset_identical": True,
            "one_point_P0_identical": True,
            "one_point_F0_identical": True,
            "ordered_neighbor_structure_identical": False,
            "spacing_order_A": SPACING_A.tolist(),
            "spacing_order_B": SPACING_B.tolist(),
        },
        "analytic_local_check": {
            "second_derivative_of_mean_lambda_squared_A": second_moment_curvature(
                SPACING_A, RATE_A
            ),
            "second_derivative_of_mean_lambda_squared_B": second_moment_curvature(
                SPACING_B, RATE_B
            ),
        },
        "final_state": {
            "empirical_KS_distance_between_spacing_distributions": empirical_ks(
                final_a, final_b
            ),
            "mean_spacing_A": float(np.mean(final_a)),
            "mean_spacing_B": float(np.mean(final_b)),
            "variance_spacing_A": float(np.var(final_a)),
            "variance_spacing_B": float(np.var(final_b)),
            "max_sorted_spacing_difference": float(
                np.max(np.abs(np.sort(final_a) - np.sort(final_b)))
            ),
            "spacing_A": final_a.tolist(),
            "spacing_B": final_b.tolist(),
            "rate_A": final_rate_a.tolist(),
            "rate_B": final_rate_b.tolist(),
        },
        "interpretation": [
            "P0 alone is not sufficient because it discards rates and spatial ordering.",
            "Even one-point F0(lambda,c) is not sufficient: the two chains have exactly the same multiset of spacing/rate pairs but different neighbor ordering.",
            "The future one-point spacing distributions differ under the same deterministic generalized-LJ dynamics and the same loading.",
            "The ordered spacing/rate state, equivalently Gamma0 for the fixed-left chain, remains an exact sufficient deterministic initialization.",
        ],
    }

    output = Path("results/data/initial_data_sufficiency/summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
