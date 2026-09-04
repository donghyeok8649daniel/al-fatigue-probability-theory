# === 한국어 파일 안내 시작 ===
# - 파일 역할: cycle-periodic snapshot P가 cumulative first passage를 유일하게 결정하는지 반례로 검사한다.
# - 핵심: 매 cycle 같은 label이 threshold를 넘는 경우와, 다른 label이 돌아가며 넘는 경우는 P가 완전히 같지만 cumulative first-passage가 다르다.
# - 이 파일은 물리 calibration이 아니라 identifiability/math audit이다.
# === 한국어 파일 안내 끝 ===
"""Audit whether periodic snapshot P determines cumulative first passage.

This is a synthetic mathematical counterexample, not an aluminum calibration.
"""
from __future__ import annotations

from pathlib import Path
import json
import math

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data" / "periodic_p_first_passage"

M = 100
LAMBDA_LOW = 1.000
LAMBDA_HIGH = 1.020
LAMBDA_C = 1.010
PHASES = np.linspace(0.0, 1.0, 101)
CYCLES_TO_REPORT = [1, 2, 5, 10, 25, 50, 100, 150]


def excursion(phase: float) -> float:
    """Smooth one-cycle excursion with one threshold crossing interval."""
    return LAMBDA_LOW + (LAMBDA_HIGH - LAMBDA_LOW) * math.sin(math.pi * phase) ** 2


def snapshot(active_label: int, phase: float) -> np.ndarray:
    values = np.full(M, LAMBDA_LOW, dtype=float)
    values[active_label] = excursion(phase)
    return values


def check_snapshot_identity() -> float:
    """Compare sorted marginals for different active labels at every phase."""
    max_diff = 0.0
    for phase in PHASES:
        a = np.sort(snapshot(0, float(phase)))
        b = np.sort(snapshot(37, float(phase)))
        max_diff = max(max_diff, float(np.max(np.abs(a - b))))
    return max_diff


def deterministic_first_passage(rotating: bool, cycles: int) -> float:
    crossed: set[int] = set()
    for n in range(cycles):
        active = n % M if rotating else 0
        if max(excursion(float(p)) for p in PHASES) >= LAMBDA_C:
            crossed.add(active)
    return len(crossed) / M


def stochastic_rare_trial(p_cycle: float, cycles: int) -> float:
    """Independent conditional rare-trial diagnostic; extra stochastic physics."""
    return 1.0 - (1.0 - p_cycle) ** cycles


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)

    max_snapshot_difference = check_snapshot_identity()

    deterministic_rows = []
    for n in CYCLES_TO_REPORT:
        deterministic_rows.append(
            {
                "cycles": n,
                "same_label_first_passage_fraction": deterministic_first_passage(False, n),
                "rotating_label_first_passage_fraction": deterministic_first_passage(True, n),
            }
        )

    p_cycle = 1.0e-5
    stochastic_cycles = [1, 10, 100, 1000, 10000, 100000, 1000000]
    stochastic_rows = [
        {
            "cycles": n,
            "cumulative_probability": stochastic_rare_trial(p_cycle, n),
        }
        for n in stochastic_cycles
    ]

    payload = {
        "classification": "periodic-P first-passage identifiability audit",
        "status": "exact synthetic counterexample; not material calibration",
        "M": M,
        "lambda_low": LAMBDA_LOW,
        "lambda_high_peak": LAMBDA_HIGH,
        "lambda_threshold": LAMBDA_C,
        "max_sorted_snapshot_difference_between_histories": max_snapshot_difference,
        "deterministic_histories": deterministic_rows,
        "stochastic_rare_trial_diagnostic": {
            "per_cycle_probability": p_cycle,
            "rows": stochastic_rows,
            "warning": "Independent per-cycle trials are extra stochastic physics and are not adopted by this audit."
        },
        "verdict": [
            "The same cycle-periodic instantaneous spacing marginal can coexist with different cumulative first-passage histories.",
            "If the same label repeats the threshold-crossing trajectory every cycle, cumulative first passage saturates after the first cycle.",
            "If the crossing label rotates while the one-point marginal stays identical, cumulative first passage grows until all labels have crossed.",
            "Therefore snapshot P alone is insufficient to determine fatigue first passage; survivor/path information is required.",
            "Permanent drift of the normalized P shape is not mathematically necessary if survivor escape continues by mixing, stochasticity, or another hidden-state mechanism."
        ]
    }

    (DATA / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
