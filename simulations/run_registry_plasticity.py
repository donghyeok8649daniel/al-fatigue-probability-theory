# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: _pulse, _write_history, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Generate the dimensionless active ideal-registry plasticity demonstration."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from theory.registry_lattice import (
    RegistryLattice,
    preferred_registry,
    registry_energy,
    registry_energy_derivative,
)
from theory.registry_plasticity import RegistryTransportConfig, solve_registry


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data" / "registry_plasticity"
FIGURES = ROOT / "results" / "figures" / "registry_plasticity"


def _pulse(time: np.ndarray) -> np.ndarray:
    force = np.zeros_like(time)
    ramp_up = (time >= 2.0) & (time < 4.0)
    force[ramp_up] = 0.55 * (time[ramp_up] - 2.0) / 2.0
    force[(time >= 4.0) & (time < 8.0)] = 0.55
    ramp_down = (time >= 8.0) & (time < 12.0)
    force[ramp_down] = 0.55 * (12.0 - time[ramp_down]) / 4.0
    return force


def _write_history(path: Path, history) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "time_reduced",
                "resolved_force_reduced",
                "mean_registry_s_over_b",
                "mean_well_index",
                "mean_intrawell_registry",
                "variance",
                "mean_lattice_energy_over_epsilon_c",
                "work_over_epsilon_c",
                "entropy_production_reduced",
                "edge_probability",
            ]
        )
        for row in zip(
            history.time,
            history.generalized_force,
            history.mean_registry,
            history.mean_well_index,
            history.mean_intrawell_registry,
            history.variance,
            history.mean_lattice_energy,
            history.work,
            history.entropy_production,
            history.boundary_probability,
        ):
            writer.writerow([float(value) for value in row])


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    lattice = RegistryLattice(normal_ratio=1.0, bessel_modes=20)
    config = RegistryTransportConfig(
        lattice=lattice,
        inverse_temperature=20.0,
        u_min=-6.0,
        u_max=7.0,
        cells=390,
    )
    time = np.linspace(0.0, 30.0, 301)
    force = _pulse(time)
    pulse = solve_registry(time, force, config, max_dt=0.025)
    # Six complete, exactly antisymmetric cycles on the same 0--30 interval.
    cyclic_force = 0.55 * np.sin(2.0 * np.pi * time / 5.0)
    symmetric = solve_registry(time, cyclic_force, config, max_dt=0.025)

    phase = np.linspace(0.0, 1.0, 2001)
    energy = np.asarray(registry_energy(phase, lattice))
    traction = np.asarray(registry_energy_derivative(phase, lattice))
    critical_force = float(np.max(np.abs(traction)))
    _write_history(DATA / "resolved_shear_pulse.csv", pulse)
    _write_history(DATA / "symmetric_cycle.csv", symmetric)
    summary = {
        "status": "dimensionless ideal single-registry mechanism demonstration",
        "m": lattice.m,
        "n": lattice.n,
        "a_over_b": lattice.normal_ratio,
        "sigma_LJ_over_b": lattice.sigma_ratio,
        "inverse_temperature_epsilon_c_over_kBT": config.inverse_temperature,
        "bessel_modes": lattice.bessel_modes,
        "preferred_registry_s_over_b": preferred_registry(lattice),
        "ideal_reduced_registry_force": critical_force,
        "pulse_peak_reduced_force": float(np.max(force)),
        "pulse_final_mean_well_index": float(pulse.mean_well_index[-1]),
        "pulse_final_mean_intrawell_registry": float(
            pulse.mean_intrawell_registry[-1]
        ),
        "pulse_final_work_over_epsilon_c": float(pulse.work[-1]),
        "pulse_max_edge_probability": float(np.max(pulse.boundary_probability)),
        "symmetric_final_mean_well_index": float(
            symmetric.mean_well_index[-1]
        ),
        "physical_plastic_shear_mapping": "gamma_p=(b/h_slip)*mean_well_index",
        "physical_axial_mapping": "epsilon_p=M_schmid*gamma_p",
        "limitations": [
            "not calibrated to aluminum",
            "one ideal registry; no dislocation storage or hardening",
            "normal separation is prescribed, not coupled dynamically",
            "mobility and representative interface area require atomistic calibration",
        ],
    }
    (DATA / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.6), constrained_layout=True)
    axes[0, 0].plot(phase, energy - np.min(energy), color="black")
    axes[0, 0].plot(phase, -0.55 * phase + energy - np.min(energy), color="tab:red")
    axes[0, 0].set(
        xlabel=r"registry $s/b$",
        ylabel=r"energy $/\varepsilon_c$",
        title="Exact Fourier--Bessel landscape and subcritical tilt",
    )
    axes[0, 0].legend(["zero load", "tilted by peak resolved force"])

    mesh = axes[0, 1].pcolormesh(
        pulse.time,
        pulse.registry,
        pulse.density.T,
        shading="auto",
        cmap="magma",
    )
    axes[0, 1].set(
        xlabel="reduced time",
        ylabel=r"unwrapped registry $s/b$",
        title="Probability crosses lattice wells",
    )
    figure.colorbar(mesh, ax=axes[0, 1], label="density")

    axes[1, 0].plot(pulse.time, pulse.generalized_force, label="resolved force")
    axes[1, 0].plot(pulse.time, pulse.mean_well_index, label=r"$\langle z\rangle$")
    axes[1, 0].plot(
        pulse.time,
        pulse.mean_intrawell_registry,
        label="mean intrawell registry",
    )
    axes[1, 0].set(
        xlabel="reduced time",
        title="Unloading recovers intrawell motion but not well population",
    )
    axes[1, 0].legend()

    axes[1, 1].plot(
        pulse.mean_registry, pulse.generalized_force, color="tab:blue", label="biased pulse"
    )
    axes[1, 1].plot(
        symmetric.mean_registry,
        symmetric.generalized_force,
        color="tab:orange",
        alpha=0.8,
        label="symmetric cycles",
    )
    axes[1, 1].set(
        xlabel=r"mean unwrapped registry $\langle s/b\rangle$",
        ylabel="reduced resolved force",
        title="Registry hysteresis and residual translation",
    )
    axes[1, 1].legend()
    figure.savefig(FIGURES / "active_registry_plasticity.png", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
