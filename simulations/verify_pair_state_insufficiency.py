"""Verify that identical P0, F0, and nearest-neighbor pair data need not determine later P.

The two chains below have the same one-point spacing multiset, zero spacing-rate
multiset, the same directed nearest-neighbor pair multiset, and the same first
and last spacing. Their triplet ordering differs. Under the same deterministic
1D generalized-LJ dynamics their later one-point spacing distributions diverge.

This is an exact reduced-information counterexample. It does not use a random
force, a named PDF family, or a statistical closure.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import numpy as np

M_EXP = 12.19
N_EXP = 6.0
DT = 1.0e-3
END_TIME = 2.0
SNAPSHOT_TIMES = (0.1, 0.2, 0.5, 1.0, 2.0)

A = 0.99
B = 1.01
C = 1.02

SPACING_A = np.array([A, A, A, B, A, C], dtype=float)
SPACING_B = np.array([A, A, B, A, A, C], dtype=float)
RATE_A = np.zeros_like(SPACING_A)
RATE_B = np.zeros_like(SPACING_B)


def normalized_lj_force(stretch: np.ndarray) -> np.ndarray:
    return (
        stretch ** (-N_EXP - 1.0) - stretch ** (-M_EXP - 1.0)
    ) / (M_EXP - N_EXP)


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


def pair_counter(values: np.ndarray) -> Counter:
    return Counter(zip(values[:-1].tolist(), values[1:].tolist()))


def triplet_counter(values: np.ndarray) -> Counter:
    return Counter(
        zip(values[:-2].tolist(), values[1:-1].tolist(), values[2:].tolist())
    )


def empirical_ks(sample_a: np.ndarray, sample_b: np.ndarray) -> float:
    grid = np.unique(np.concatenate((sample_a, sample_b)))
    cdf_a = np.array([np.mean(sample_a <= value) for value in grid])
    cdf_b = np.array([np.mean(sample_b <= value) for value in grid])
    return float(np.max(np.abs(cdf_a - cdf_b)))


def simulate(spacing0: np.ndarray, rate0: np.ndarray):
    x, velocity = state_from_spacing_and_rate(spacing0, rate0)
    force = node_force(x)
    snapshots: dict[float, np.ndarray] = {}
    targets = list(SNAPSHOT_TIMES)
    target_index = 0
    n_steps = int(round(END_TIME / DT))

    for step in range(n_steps):
        velocity[1:] += 0.5 * DT * force[1:]
        x[1:] += DT * velocity[1:]
        new_force = node_force(x)
        velocity[1:] += 0.5 * DT * new_force[1:]
        force = new_force

        t = (step + 1) * DT
        while target_index < len(targets) and t + 1.0e-12 >= targets[target_index]:
            snapshots[targets[target_index]] = np.diff(x).copy()
            target_index += 1

    return snapshots


def counter_to_json(counter: Counter) -> dict[str, int]:
    return {str(tuple(key)): int(value) for key, value in sorted(counter.items())}


def main() -> None:
    assert np.array_equal(np.sort(SPACING_A), np.sort(SPACING_B))
    assert np.array_equal(RATE_A, RATE_B)
    assert SPACING_A[0] == SPACING_B[0]
    assert SPACING_A[-1] == SPACING_B[-1]
    assert pair_counter(SPACING_A) == pair_counter(SPACING_B)
    assert triplet_counter(SPACING_A) != triplet_counter(SPACING_B)

    snapshots_a = simulate(SPACING_A, RATE_A)
    snapshots_b = simulate(SPACING_B, RATE_B)

    time_results = []
    for t in SNAPSHOT_TIMES:
        sample_a = snapshots_a[t]
        sample_b = snapshots_b[t]
        time_results.append(
            {
                "time": t,
                "empirical_KS_distance": empirical_ks(sample_a, sample_b),
                "max_sorted_spacing_difference": float(
                    np.max(np.abs(np.sort(sample_a) - np.sort(sample_b)))
                ),
                "mean_spacing_A": float(np.mean(sample_a)),
                "mean_spacing_B": float(np.mean(sample_b)),
                "variance_spacing_A": float(np.var(sample_a)),
                "variance_spacing_B": float(np.var(sample_b)),
                "spacing_A": sample_a.tolist(),
                "spacing_B": sample_b.tolist(),
            }
        )

    result = {
        "classification": "deterministic nearest-neighbor pair-state insufficiency counterexample",
        "model": {
            "m": M_EXP,
            "n": N_EXP,
            "dt": DT,
            "end_time": END_TIME,
            "external_force": 0.0,
            "number_of_spacings": int(SPACING_A.size),
            "boundary_condition": "fixed-left / zero-force-right",
        },
        "initial_state_A": SPACING_A.tolist(),
        "initial_state_B": SPACING_B.tolist(),
        "initial_checks": {
            "P0_identical": True,
            "F0_identical": True,
            "directed_nearest_neighbor_pair_measure_identical": True,
            "first_spacing_identical": True,
            "last_spacing_identical": True,
            "triplet_measure_identical": False,
            "directed_pair_counts": counter_to_json(pair_counter(SPACING_A)),
            "triplet_counts_A": counter_to_json(triplet_counter(SPACING_A)),
            "triplet_counts_B": counter_to_json(triplet_counter(SPACING_B)),
        },
        "time_results": time_results,
        "main_result": {
            "KS_at_time_1": next(
                item["empirical_KS_distance"] for item in time_results if item["time"] == 1.0
            ),
            "max_sorted_spacing_difference_at_time_1": next(
                item["max_sorted_spacing_difference"]
                for item in time_results
                if item["time"] == 1.0
            ),
        },
        "interpretation": [
            "Nearest-neighbor pair information is enough to form the instantaneous one-point mean-neighbor force terms, but it is not an autonomous state.",
            "The pair transport equation needs outside neighbors and therefore triplet information.",
            "Identical P0, F0, directed pair data, and boundary endpoint spacings do not guarantee identical later one-point P.",
            "Nearest-neighbor mechanics produces a local correlation hierarchy F1 <- F2 <- F3 <- ... rather than an exact pair closure.",
        ],
    }

    output = Path("results/data/pair_state_insufficiency/summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
