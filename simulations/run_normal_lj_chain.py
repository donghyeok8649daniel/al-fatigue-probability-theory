# === 한국어 파일 안내 시작 ===
# - 파일 역할: 보존적 1D layer-LJ chain을 실행하고 cycle history, energy balance, instability 관련 결과를 저장한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Run the main normal-opening generalized-LJ proof-of-principle cases."""

from __future__ import annotations

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from theory.normal_lj_chain import (
    NormalLJParameters,
    atomic_time_scale,
    critical_dimensionless_force,
    critical_stretch,
    dimensionless_omega_from_frequency,
    normalized_lj_force,
    physical_frequency_from_dimensionless_omega,
    simulate_normal_lj_chain,
    stress_to_dimensionless_force,
)

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures"
DATA = ROOT / "results" / "data"


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    youngs_modulus = 69.0e9
    stress_100mpa = 100.0e6
    f_100 = stress_to_dimensionless_force(stress_100mpa, youngs_modulus)

    cases = {
        "100MPa_equivalent": simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=f_100, omega=0.02),
            cycles=12,
        ),
        "subcritical_slow": simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=0.03, omega=0.01),
            cycles=12,
        ),
        "subcritical_dynamic": simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=0.03, omega=0.02),
            cycles=5,
        ),
    }

    sweep = {}
    for omega in (0.01, 0.02, 0.05, 0.10):
        result = simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=0.03, omega=omega),
            cycles=5,
        )
        sweep[str(omega)] = (
            result.first_instability.cycle
            if result.first_instability is not None
            else None
        )

    lam_c = critical_stretch()
    force_c = critical_dimensionless_force()

    atomic_mass_al = 26.9815385 * 1.66053906660e-27
    a0 = 2.8627442948e-10
    area0 = 6.0338e-20
    t0 = atomic_time_scale(atomic_mass_al, a0, youngs_modulus, area0)

    summary = {
        "lambda_c": lam_c,
        "dimensionless_static_critical_force": force_c,
        "static_critical_stress_pa_if_sigma_over_E_mapping": force_c * youngs_modulus,
        "dimensionless_100MPa_amplitude": f_100,
        "frequency_sweep_force_amplitude": 0.03,
        "frequency_sweep_first_instability_cycle": sweep,
        "100MPa_first_instability": (
            None
            if cases["100MPa_equivalent"].first_instability is None
            else cases["100MPa_equivalent"].first_instability.__dict__
        ),
        "100MPa_energy_balance_relative_error":
            cases["100MPa_equivalent"].energy_balance_relative_error,
        "omega_star_0p02_physical_frequency_hz":
            physical_frequency_from_dimensionless_omega(0.02, t0),
        "20Hz_dimensionless_omega_star":
            dimensionless_omega_from_frequency(20.0, t0),
    }

    (DATA / "normal_lj_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    with (DATA / "normal_lj_cycle_history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.writer(f)
        writer.writerow(
            ["case", "cycle", "mean_spacing", "variance", "max_spacing", "min_spacing"]
        )
        for name, result in cases.items():
            for idx in range(len(result.cycle_mean_spacing)):
                writer.writerow(
                    [
                        name,
                        idx + 1,
                        result.cycle_mean_spacing[idx],
                        result.cycle_variance_spacing[idx],
                        result.cycle_max_spacing[idx],
                        result.cycle_min_spacing[idx],
                    ]
                )

    lam = np.linspace(0.94, 1.18, 500)
    plt.figure(figsize=(7, 5))
    plt.plot(lam, normalized_lj_force(lam))
    plt.axvline(lam_c, linestyle="--", label=f"stability loss λc={lam_c:.4f}")
    plt.axhline(f_100, linestyle=":", label="100 MPa / 69 GPa")
    plt.xlabel("Normal stretch λ = a/a0")
    plt.ylabel("Dimensionless normal stress σ/E")
    plt.title("Generalized-LJ normal traction–stretch relation")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_traction_stretch.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    for name, result in cases.items():
        n = np.arange(1, len(result.cycle_max_spacing) + 1)
        if len(n):
            plt.plot(n, result.cycle_max_spacing, marker="o", label=name)
    plt.axhline(lam_c, linestyle="--", label="LJ local stability limit")
    plt.xlabel("Cycle number")
    plt.ylabel("Maximum bond stretch max(a_i/a0)")
    plt.title("Pure normal generalized-LJ chain: cycle-end state")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_cycle_max_spacing.png", dpi=180)
    plt.close()

    result = cases["subcritical_dynamic"]
    plt.figure(figsize=(8, 5))
    for cycle in sorted(result.cycle_snapshots):
        if cycle <= 2:
            values = result.cycle_snapshots[cycle]
            plt.hist(
                values,
                bins=12,
                density=True,
                histtype="step",
                linewidth=1.8,
                label=f"cycle {cycle}",
            )
    plt.axvline(lam_c, linestyle="--", label="λc")
    plt.xlabel("Local normal spacing λ_i = a_i/a0")
    plt.ylabel("Finite empirical density P_N")
    plt.title("Normal-spacing distribution before dynamic instability")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_spacing_distribution.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    x = np.asarray([float(k) for k in sweep])
    y = np.asarray([np.nan if sweep[k] is None else sweep[k] for k in sweep])
    plt.scatter(x, y)
    for omega, cycle in zip(x, y):
        if np.isnan(cycle):
            plt.annotate(
                "no crossing in 5 cycles",
                (omega, 4.6),
                ha="center",
                rotation=30,
            )
    plt.xlabel("Dimensionless angular frequency ω*")
    plt.ylabel("First λc crossing (cycle)")
    plt.title("Sub-static-critical normal forcing: frequency dependence")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_frequency_sweep.png", dpi=180)
    plt.close()

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
