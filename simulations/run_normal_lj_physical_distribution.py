# === 한국어 파일 안내 시작 ===
# - 파일 역할: current layer-LJ calibration에서 force-biased stable/barrier point와 metastable thermal P의 dimensionless 진단값을 생성한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: main
# - 출력: results/data/normal_lj_physical_distribution.csv, .json 및 results/reports/NORMAL_LJ_PHYSICAL_DISTRIBUTION.md
# - 주의: chi는 dimensionless diagnostic 값이며 A0가 물리적으로 결정되기 전에는 실제 Al 온도/수명 예측으로 해석하지 않는다.
# === 한국어 파일 안내 끝 ===
"""Generate dimensionless physical-statistical diagnostics for Milestone 12."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from theory.normal_lj_chain import critical_dimensionless_force, critical_stretch
from theory.normal_lj_physical_distribution import (
    metastable_barrier_height,
    metastable_gibbs_density,
    metastable_stationary_points,
    metastable_tail_probability,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
REPORTS = ROOT / "results" / "reports"


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    fc = critical_dimensionless_force()
    lam_c = critical_stretch()
    force_fractions = (0.10, 0.30, 0.50, 0.70, 0.90, 0.97)
    chi_values = (40.0, 100.0, 250.0, 500.0)

    rows: list[dict[str, float]] = []
    barriers: list[dict[str, float]] = []

    for fraction in force_fractions:
        f = fraction * fc
        stable, barrier = metastable_stationary_points(f)
        delta_w = metastable_barrier_height(f)
        barriers.append(
            {
                "force_fraction_of_fc": fraction,
                "dimensionless_force": f,
                "stable_spacing": stable,
                "barrier_spacing": barrier,
                "barrier_height_delta_w": delta_w,
            }
        )

        lower = min(0.72, 0.85 * stable)
        grid = np.linspace(lower, barrier * (1.0 - 1.0e-8), 20001)
        for chi in chi_values:
            density = metastable_gibbs_density(grid, f, chi)
            mean = float(np.trapezoid(grid * density, grid))
            variance = float(np.trapezoid((grid - mean) ** 2 * density, grid))
            tail = metastable_tail_probability(grid, density)
            rows.append(
                {
                    "force_fraction_of_fc": fraction,
                    "dimensionless_force": f,
                    "chi": chi,
                    "stable_spacing": stable,
                    "lambda_c": lam_c,
                    "barrier_spacing": barrier,
                    "barrier_height_delta_w": delta_w,
                    "mean_spacing": mean,
                    "variance_spacing": variance,
                    "metastable_tail_Qc": tail,
                }
            )

    csv_path = DATA / "normal_lj_physical_distribution.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "classification": "dimensionless model diagnostic under metastable local-equilibrium assumption",
        "m": 12.19,
        "n": 6.0,
        "critical_dimensionless_force": fc,
        "lambda_c": lam_c,
        "warning": "chi is not mapped to an aluminum temperature until representative layer area A0 is physically fixed",
        "barrier_curve": barriers,
        "metastable_density_diagnostics": rows,
    }
    (DATA / "normal_lj_physical_distribution.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    report = f"""# Normal-LJ Physical Statistical Distribution Diagnostic

## Classification

**DIMENSIONLESS MODEL DIAGNOSTIC / CONTROLLED METASTABLE APPROXIMATION.**

This report evaluates the full nonlinear reduced layer-LJ potential. It does not assign a physical aluminum temperature because the representative coarse-grained area $A_0$ has not yet been fixed.

The critical reduced force and tangent-instability spacing are

$$
f_c={fc:.12g},
\qquad
\lambda_c={lam_c:.12g}.
$$

For each $f/f_c$, the stable point $\lambda_s$, unstable barrier point $\lambda_b$, and

$$
\Delta w=w_f(\lambda_b)-w_f(\lambda_s)
$$

are computed without a Taylor expansion. Metastable densities are conditioned on $0<\lambda<\lambda_b$ and evaluated for illustrative dimensionless $\chi$ values.

The CSV and JSON files contain the complete numerical table. The important qualitative checks are:

- $\Delta w>0$ for every tested $0<f<f_c$;
- the barrier decreases as $f\to f_c^-$;
- increasing $\chi$ concentrates the intact-basin density near $\lambda_s$;
- the reported $Q_c$ is an instantaneous basin population, not an escape rate or fatigue life.

---

# 한국어 번역 — Normal-LJ 물리 통계분포 진단

## 분류

**DIMENSIONLESS MODEL DIAGNOSTIC / CONTROLLED METASTABLE APPROXIMATION.**

이 보고서는 full nonlinear reduced layer-LJ potential을 계산한다. representative coarse-grained area $A_0$가 아직 물리적으로 정해지지 않았으므로 실제 aluminum temperature를 부여하지 않는다.

critical reduced force와 tangent-instability spacing은

$$
f_c={fc:.12g},
\qquad
\lambda_c={lam_c:.12g}
$$

이다.

각 $f/f_c$에 대해 stable point $\lambda_s$, unstable barrier point $\lambda_b$ 및

$$
\Delta w=w_f(\lambda_b)-w_f(\lambda_s)
$$

를 Taylor expansion 없이 계산한다. metastable density는 $0<\lambda<\lambda_b$에 조건부로 두고 illustrative dimensionless $\chi$ 값에서 계산한다.

전체 수치표는 CSV와 JSON에 저장한다. 핵심 qualitative check는 다음과 같다.

- 모든 tested $0<f<f_c$에서 $\Delta w>0$;
- $f\to f_c^-$이면 barrier가 감소;
- $\chi$가 증가하면 intact-basin density가 $\lambda_s$ 근처로 집중;
- $Q_c$는 instantaneous basin population이며 escape rate나 fatigue life가 아님.
"""
    (REPORTS / "NORMAL_LJ_PHYSICAL_DISTRIBUTION.md").write_text(
        report, encoding="utf-8"
    )


if __name__ == "__main__":
    main()
