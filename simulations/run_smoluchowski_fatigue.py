# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: history, cycle_areas, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Generate the uncalibrated dimensionless Smoluchowski demonstration."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from theory.smoluchowski_escape import TransportConfig, solve


ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures" / "smoluchowski"
DATA = ROOT / "results" / "data" / "smoluchowski"


def history(period: float, boundary: str, cycles: int = 6):
    time = np.linspace(0.0, cycles * period, cycles * 100 + 1)
    force = 0.008 + 0.007 * np.sin(2 * np.pi * time / period)
    initiation = "tangent_instability" if boundary == "absorbing" else "fixed_coordinate"
    return solve(time, force, TransportConfig(inverse_temperature=2000, cells=180,
                 lambda_max=1.34, boundary=boundary, initiation_definition=initiation),
                 max_dt=min(0.03, period/100))


def cycle_areas(h, period):
    boundaries = np.arange(0, h.time[-1] + 0.5 * period, period)
    return np.diff(np.interp(boundaries, h.time, h.work))


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True); DATA.mkdir(parents=True, exist_ok=True)
    reflecting = history(12.0, "reflecting")
    absorbing = history(12.0, "absorbing")

    fig, axes = plt.subplots(3, 2, figsize=(11, 10), constrained_layout=True)
    axes[0, 0].plot(reflecting.time, reflecting.force)
    axes[0, 0].set(xlabel=r"reduced time $t/t_r$", ylabel=r"force $F a_0/E_0$")
    image = axes[0, 1].pcolormesh(reflecting.time, reflecting.stretch,
                                  reflecting.density.T, shading="auto")
    fig.colorbar(image, ax=axes[0, 1], label=r"$p(\lambda,t)$")
    axes[0, 1].set(xlabel="reduced time", ylabel=r"$\lambda$")
    axes[1, 0].plot(reflecting.time, reflecting.mean, label="mean")
    axes[1, 0].plot(reflecting.time, reflecting.variance, label="variance")
    axes[1, 0].plot(reflecting.time, reflecting.skewness, label="skewness")
    axes[1, 0].legend(); axes[1, 0].set_xlabel("reduced time")
    axes[1, 1].plot(reflecting.time, reflecting.mean_shifted_energy, label="mean U")
    axes[1, 1].plot(reflecting.time, reflecting.tail_conditional, label="tail")
    axes[1, 1].legend(); axes[1, 1].set_xlabel("reduced time")
    axes[2, 0].plot(reflecting.mean[-200:], reflecting.force[-200:])
    axes[2, 0].set(xlabel=r"$\bar\lambda$", ylabel="reduced force", title="last two loops")
    areas = cycle_areas(reflecting, 12.0)
    axes[2, 1].plot(np.arange(1, len(areas)+1), areas, "o-")
    axes[2, 1].set(xlabel="cycle", ylabel=r"$H_k/E_0$", title="loop work")
    fig.suptitle("Dimensionless reflecting demonstration (not Al lifetime prediction)")
    fig.savefig(FIG / "reflecting_overview.png", dpi=180); plt.close(fig)

    # Equal force (mean force) on the rising and falling branches of cycle 6.
    load_i, unload_i = 525, 575
    fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    ax.plot(reflecting.stretch, reflecting.density[load_i], label="loading")
    ax.plot(reflecting.stretch, reflecting.density[unload_i], "--", label="unloading")
    ax.set(xlabel=r"$\lambda$", ylabel="density",
           title="Different distributions at the same force")
    ax.legend(); fig.savefig(FIG / "same_force_phase_lag.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    axes[0, 0].plot(absorbing.time, absorbing.survival)
    axes[0, 0].set(title="survival", xlabel="reduced time", ylabel="$S$")
    axes[0, 1].plot(absorbing.time, absorbing.initiation)
    axes[0, 1].set(title="cumulative initiation", xlabel="reduced time", ylabel="$1-S$")
    axes[1, 0].plot(absorbing.time, absorbing.hazard)
    axes[1, 0].set(title="hazard", xlabel="reduced time", ylabel="$h t_r$")
    axes[1, 1].plot(absorbing.time, absorbing.outflux)
    axes[1, 1].set(title="first-passage outflux", xlabel="reduced time",
                   ylabel="$J_{out}t_r$")
    fig.suptitle(r"Crack initiation: first passage through $\lambda_c$")
    fig.savefig(FIG / "escape_observables.png", dpi=180); plt.close(fig)

    periods = np.asarray([0.5, 2, 6, 12, 30, 80, 200.0])
    loop, escape = [], []
    for period in periods:
        ref = history(float(period), "reflecting", cycles=4)
        esc = history(float(period), "absorbing", cycles=4)
        loop.append(abs(cycle_areas(ref, float(period))[-1]))
        escape.append(esc.initiation[-1])
    fig, ax1 = plt.subplots(figsize=(7, 4.7), constrained_layout=True)
    ax1.loglog(1/periods, loop, "o-", label="loop area")
    ax1.set(xlabel=r"reduced frequency $f t_r$", ylabel="loop area")
    ax2 = ax1.twinx(); ax2.semilogx(1/periods, escape, "s--", color="tab:red")
    ax2.set_ylabel("four-cycle escape probability", color="tab:red")
    ax1.set_title("Frequency response: phase lag versus irreversible escape")
    fig.savefig(FIG / "frequency_sweep.png", dpi=180); plt.close(fig)

    with (DATA / "absorbing_history.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp); writer.writerow(["time_reduced", "force_reduced", "mean_lambda",
            "variance", "skewness", "mean_shifted_energy", "tail_conditional", "survival",
            "hazard_reduced", "initiation"])
        writer.writerows(zip(absorbing.time, absorbing.force, absorbing.mean,
            absorbing.variance, absorbing.skewness, absorbing.mean_shifted_energy,
            absorbing.tail_conditional, absorbing.survival, absorbing.hazard,
            absorbing.initiation))
    summary = {"status": "dimensionless demonstration; not calibrated aluminum life",
               "m": 12.19, "n": 6.0, "inverse_temperature": 2000.0,
               "initiation_definition": "first passage through tangent-instability stretch",
               "periods": periods.tolist(), "frequencies": (1/periods).tolist(),
               "loop_areas": loop, "four_cycle_escape": escape,
               "same_force_density_l1_difference": float(np.sum(np.abs(
                   reflecting.density[load_i]-reflecting.density[unload_i]))
                   * (reflecting.stretch[1]-reflecting.stretch[0])),
               "final_mean_offset": float(absorbing.mean[-1]-absorbing.mean[0]),
               "final_variance_change": float(absorbing.variance[-1]-absorbing.variance[0]),
               "final_survival": float(absorbing.survival[-1])}
    (DATA / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
