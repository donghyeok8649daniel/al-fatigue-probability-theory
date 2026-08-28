# === 한국어 파일 안내 시작 ===
# - 파일 역할: push-forward 관련 과거 harmonic/Taylor 진단 결과를 재현한다. active 전역 분포 모델은 아니다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: empirical_kolmogorov_distance, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Evaluate the governing-equation push-forward clue against existing 1D snapshots."""
from __future__ import annotations

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from theory.normal_lj_chain import NormalLJParameters, simulate_normal_lj_chain
from theory.normal_lj_closure_validation import compare_snapshot_to_closure
from theory.normal_lj_pushforward import (
    arcsine_cdf,
    arcsine_density,
    two_harmonic_max_abs_skewness,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIG = ROOT / "results" / "figures"

CASES = (
    {
        "name": "slow_t10T",
        "parameters": NormalLJParameters(force_amplitude=0.03, omega=0.01),
        "cycles": 12,
        "sample_cycle": 10,
    },
    {
        "name": "dynamic_t2T",
        "parameters": NormalLJParameters(force_amplitude=0.03, omega=0.02),
        "cycles": 3,
        "sample_cycle": 2,
    },
)


def empirical_kolmogorov_distance(values, model_cdf) -> float:
    x = np.sort(np.asarray(values, dtype=float))
    F = np.asarray(model_cdf(x), dtype=float)
    M = len(x)
    upper = np.arange(1, M + 1, dtype=float) / M
    lower = np.arange(0, M, dtype=float) / M
    return float(max(np.max(np.abs(upper - F)), np.max(np.abs(lower - F))))


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rows = []
    skew_bound = two_harmonic_max_abs_skewness()

    for case in CASES:
        result = simulate_normal_lj_chain(case["parameters"], cycles=case["cycles"])
        values = result.cycle_snapshots[case["sample_cycle"]]
        mean = float(np.mean(values))
        variance = float(np.var(values))
        centered = values - mean
        skewness = float(np.mean(centered ** 3) / variance ** 1.5)
        amplitude = float(np.sqrt(2.0 * variance))
        arcsine_ks = empirical_kolmogorov_distance(
            values,
            lambda x, mu=mean, A=amplitude: arcsine_cdf(x, mu, A),
        )
        closure = compare_snapshot_to_closure(
            values,
            closure_quadrature_order=640,
            cdf_quadrature_order=128,
        )
        row = {
            "case": case["name"],
            "represented_spacings": len(values),
            "mean_stretch": mean,
            "variance": variance,
            "skewness": skewness,
            "single_mode_amplitude_from_variance": amplitude,
            "single_mode_arcsine_KS": arcsine_ks,
            "mean_energy_exponential_closure_KS": closure.kolmogorov_distance,
            "two_harmonic_max_abs_skewness": skew_bound,
            "exceeds_two_harmonic_skewness_bound": abs(skewness) > skew_bound,
        }
        rows.append(row)

        lo = mean - 1.15 * amplitude
        hi = mean + 1.15 * amplitude
        grid = np.linspace(lo, hi, 2000)
        plt.figure(figsize=(8, 5))
        plt.hist(
            values,
            bins=12,
            density=True,
            histtype="step",
            linewidth=2,
            label="deterministic 1D layer-LJ",
        )
        density = arcsine_density(grid, mean, amplitude)
        finite = np.isfinite(density)
        plt.plot(grid[finite], density[finite], label="single linear-mode push-forward")
        plt.xlabel("Normalized layer spacing lambda")
        plt.ylabel("Density")
        plt.title(case["name"])
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(FIG / f"normal_lj_pushforward_{case['name']}.svg")
        plt.close()

    with (DATA / "normal_lj_pushforward_clue.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "classification": "governing-equation push-forward structural clue plus numerical falsification of single-mode sufficiency",
        "single_mode_origin": "linearized 1D layer-LJ normal mode sampled uniformly in spatial phase",
        "two_harmonic_skewness_bound": skew_bound,
        "rows": rows,
        "conclusion": "The one-point density is exactly a spatial push-forward. A single linear mode gives an arcsine density, but the tested driven states require richer spatial mode content; the slow snapshot is too skewed even for a first-plus-second-harmonic waveform.",
    }
    (DATA / "normal_lj_pushforward_clue.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
