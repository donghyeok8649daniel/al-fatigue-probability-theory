# === 한국어 파일 안내 시작 ===
# - 파일 역할: omega*M을 고정한 동적 유사성 조건에서 represented system size에 따른 closure mismatch를 검사한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: dynamically_matched_omega, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Dynamically matched represented-system-size sweep for the 1D layer-LJ closure.

A naive fixed-frequency M sweep is not a clean finite-size test because the
acoustic transit count changes with chain length. This script instead keeps

    omega * M = constant

where M is the number of represented layer spacings. With a linear acoustic
speed of order unity in the normalized chain, this keeps the loading period
relative to the chain transit time approximately fixed.

The comparison remains a controlled numerical diagnostic, not a proof of an
M->infinity limit.
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from theory.normal_lj_chain import NormalLJParameters, simulate_normal_lj_chain
from theory.normal_lj_closure_validation import compare_snapshot_to_closure

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIG = ROOT / "results" / "figures"

ATOM_COUNTS = (32, 64, 128, 256)
REFERENCE_SPACINGS = 31
REFERENCE_OMEGA = 0.02
OMEGA_TIMES_M = REFERENCE_OMEGA * REFERENCE_SPACINGS
FORCE_AMPLITUDE = 0.03
SAMPLE_CYCLE = 2


def dynamically_matched_omega(represented_spacings: int) -> float:
    """Angular frequency for the controlled scaling omega*M=constant."""
    if represented_spacings < 2:
        raise ValueError("represented_spacings must be at least 2")
    return OMEGA_TIMES_M / represented_spacings


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    rows = []
    for atoms in ATOM_COUNTS:
        represented_spacings = atoms - 1
        omega = dynamically_matched_omega(represented_spacings)
        result = simulate_normal_lj_chain(
            NormalLJParameters(
                force_amplitude=FORCE_AMPLITUDE,
                omega=omega,
            ),
            atoms=atoms,
            cycles=SAMPLE_CYCLE,
            record_stride=10_000_000,
        )
        values = result.cycle_snapshots[SAMPLE_CYCLE]
        comparison = compare_snapshot_to_closure(
            values,
            closure_quadrature_order=640,
            cdf_quadrature_order=128,
        )
        target_time = SAMPLE_CYCLE * result.period
        row = {
            "atoms": atoms,
            "represented_spacings": represented_spacings,
            "omega": omega,
            "omega_times_M": omega * represented_spacings,
            "target_time_star": target_time,
            "target_time_over_M": target_time / represented_spacings,
            **asdict(comparison),
        }
        rows.append(row)

    with (DATA / "normal_lj_closure_system_size.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    inv_m = 1.0 / np.asarray([row["represented_spacings"] for row in rows], dtype=float)
    ks = np.asarray([row["kolmogorov_distance"] for row in rows], dtype=float)
    skew_sim = np.asarray([row["empirical_skewness"] for row in rows], dtype=float)
    skew_closure = np.asarray([row["closure_skewness"] for row in rows], dtype=float)

    ks_fit = np.polyfit(inv_m, ks, 1)
    skew_sim_fit = np.polyfit(inv_m, skew_sim, 1)
    skew_closure_fit = np.polyfit(inv_m, skew_closure, 1)

    summary = {
        "classification": "controlled numerical represented-system-size diagnostic",
        "dynamic_similarity_rule": "omega * M = constant",
        "omega_times_M": OMEGA_TIMES_M,
        "force_amplitude": FORCE_AMPLITUDE,
        "sample_cycle": SAMPLE_CYCLE,
        "rows": rows,
        "exploratory_linear_in_inverse_M_extrapolation": {
            "classification": "EXPLORATORY NUMERICAL EXTRAPOLATION; not a theorem",
            "KS_intercept_M_to_infinity": float(ks_fit[1]),
            "empirical_skewness_intercept_M_to_infinity": float(skew_sim_fit[1]),
            "closure_skewness_intercept_M_to_infinity": float(skew_closure_fit[1]),
        },
    }
    (DATA / "normal_lj_closure_system_size.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    m_values = np.asarray([row["represented_spacings"] for row in rows], dtype=float)

    plt.figure(figsize=(7.5, 5.0))
    plt.plot(m_values, ks, marker="o")
    plt.xlabel("Represented spacings M")
    plt.ylabel("Kolmogorov distance")
    plt.title("1D layer-LJ closure: dynamically matched system-size sweep")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_closure_system_size_ks.svg")
    plt.close()

    plt.figure(figsize=(7.5, 5.0))
    plt.plot(m_values, skew_sim, marker="o", label="deterministic 1D layer-LJ")
    plt.plot(m_values, skew_closure, marker="o", label="same mean + energy closure")
    plt.xlabel("Represented spacings M")
    plt.ylabel("Skewness")
    plt.title("Skewness mismatch versus represented system size")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_closure_system_size_skewness.svg")
    plt.close()

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
