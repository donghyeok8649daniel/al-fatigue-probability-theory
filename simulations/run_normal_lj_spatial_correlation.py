# === 한국어 파일 안내 시작 ===
# - 파일 역할: 동적으로 matched된 여러 chain size에서 spacing spatial correlation을 계산하고 CSV/JSON/figure를 생성한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Spatial-correlation sweep for dynamically matched 1D layer-LJ chains."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from simulations.run_normal_lj_closure_system_size import (
    ATOM_COUNTS,
    FORCE_AMPLITUDE,
    SAMPLE_CYCLE,
    dynamically_matched_omega,
)
from theory.normal_lj_chain import NormalLJParameters, simulate_normal_lj_chain
from theory.normal_lj_spatial_correlation import (
    correlation_profile,
    random_permutation_expected_rho,
    summarize_spatial_correlation,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIG = ROOT / "results" / "figures"
SCALED_LAGS = (0.02, 0.05, 0.10, 0.20, 0.30, 0.40)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    profile_rows = []
    scaled_rows = []

    for atoms in ATOM_COUNTS:
        m_count = atoms - 1
        omega = dynamically_matched_omega(m_count)
        result = simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=FORCE_AMPLITUDE, omega=omega),
            atoms=atoms,
            cycles=SAMPLE_CYCLE,
            record_stride=10_000_000,
        )
        values = result.cycle_snapshots[SAMPLE_CYCLE]
        summary = summarize_spatial_correlation(values)
        lags, covariance, rho = correlation_profile(values)
        summary_rows.append({
            "atoms": atoms,
            "M": m_count,
            "omega": omega,
            "omega_times_M": omega * m_count,
            "mean_stretch": float(np.mean(values)),
            **asdict(summary),
            "random_permutation_expected_rho_k": random_permutation_expected_rho(m_count),
        })
        for k, c_k, rho_k in zip(lags, covariance, rho):
            profile_rows.append({
                "M": m_count,
                "omega": omega,
                "lag_k": int(k),
                "scaled_lag_eta": float(k / m_count),
                "C_k": float(c_k),
                "rho_k": float(rho_k),
            })
        for eta in SCALED_LAGS:
            k = max(1, int(round(eta * m_count)))
            scaled_rows.append({
                "M": m_count,
                "omega": omega,
                "eta_target": eta,
                "lag_k": k,
                "eta_actual": k / m_count,
                "rho_k": float(rho[k]),
            })

    def write_csv(path: Path, rows):
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(DATA / "normal_lj_spatial_correlation_summary.csv", summary_rows)
    write_csv(DATA / "normal_lj_spatial_correlation_profile.csv", profile_rows)
    write_csv(DATA / "normal_lj_spatial_correlation_scaled.csv", scaled_rows)

    payload = {
        "classification": "controlled numerical spatial-correlation diagnostic in the active 1D normal layer-LJ model",
        "protocol": {
            "force_amplitude": FORCE_AMPLITUDE,
            "sample_cycle": SAMPLE_CYCLE,
            "dynamic_similarity_rule": "omega * M = 0.62",
        },
        "exact_permutation_null": "For any nonzero lag k under a uniformly random permutation of a centered finite sample, E[rho_k] = -1/(M-1).",
        "summary_rows": summary_rows,
        "scaled_lag_rows": scaled_rows,
    }
    (DATA / "normal_lj_spatial_correlation.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    plt.figure(figsize=(7.5, 5.0))
    for row in summary_rows:
        m_count = row["M"]
        selected = [r for r in profile_rows if r["M"] == m_count]
        plt.plot([r["scaled_lag_eta"] for r in selected], [r["rho_k"] for r in selected], label=f"M={m_count}")
    plt.axhline(0.0, linestyle="--")
    plt.xlim(0.0, 0.5)
    plt.xlabel("Scaled lag eta = k/M")
    plt.ylabel("Normalized spacing correlation rho_k")
    plt.title("1D layer-LJ spatial correlation collapse")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_spatial_correlation_scaled.svg")
    plt.close()

    plt.figure(figsize=(7.5, 5.0))
    m_values = [row["M"] for row in summary_rows]
    plt.plot(m_values, [row["rho1"] for row in summary_rows], marker="o", label="deterministic rho_1")
    plt.plot(m_values, [row["random_permutation_expected_rho_k"] for row in summary_rows], marker="o", label="random-permutation expectation")
    plt.xlabel("Represented spacings M")
    plt.ylabel("Nearest-neighbor normalized correlation")
    plt.title("Ordering information absent from one-point P(lambda,t)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_spatial_correlation_rho1.svg")
    plt.close()

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
