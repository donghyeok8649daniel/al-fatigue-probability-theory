# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: write_csv, render_report, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Milestone 15: quasistatic-limit diagnosis of the current correlation protocol."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt

from theory.normal_lj_chain import NormalLJParameters, simulate_normal_lj_chain
from theory.normal_lj_quasistatic_protocol import (
    cycle_boundary_force,
    residual_snapshot_metrics,
    stable_stretch_for_tensile_force,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIG = ROOT / "results" / "figures"
REPORT = ROOT / "results" / "reports"

REPRESENTED_SPACINGS = (31, 63)
OMEGA_M_VALUES = (0.62, 0.31, 0.155, 0.0775)
FORCE_AMPLITUDE = 0.03
SAMPLE_CYCLE = 2
INTEGRATION_CYCLES = 2.01
DT = 0.05


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def render_report(rows: list[dict]) -> str:
    lines = [
        "# Normal-LJ Quasistatic Protocol Diagnostic",
        "",
        "## Classification",
        "",
        "**EXACT STATIC RESULT + CONTROLLED NUMERICAL DYNAMICAL DIAGNOSTIC**",
        "",
        "For the homogeneous open chain under a constant tensile end force $f$,",
        "",
        "$$",
        "\\Pi=\\sum_{i=1}^{M}[\\phi(\\lambda_i)-f\\lambda_i].",
        "$$",
        "",
        "On the stable branch $\\phi''>0$, stationarity gives $\\phi'(\\lambda_i)=f$ for every spacing, hence all spacings equal the unique stable root $\\lambda_s(f)$. Therefore the exact zero-temperature quasistatic empirical distribution is",
        "",
        "$$",
        "P_M^{\\rm qs}(\\lambda\\mid f)=\\delta[\\lambda-\\lambda_s(f)].",
        "$$",
        "",
        "The Milestone 13 snapshots are taken at integer cycle 2 for a zero-mean sine load. At that exact phase the applied force is zero, so the quasistatic reference is $\\lambda_i=1$ and $C_0=0$. Any nonzero variance in those snapshots is therefore dynamical residual structure, not static material randomness.",
        "",
        "## Numerical sweep",
        "",
        f"The sweep uses force amplitude `{FORCE_AMPLITUDE}`, sample cycle `{SAMPLE_CYCLE}`, and decreases the protocol parameter $\\alpha=\\omega M$ while keeping the model strictly one-dimensional.",
        "",
        "| M | alpha=omega M | C0 | sqrt(C0) | |mean-1| | rho1 | M_eff^(+) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['M']} | {row['omega_times_M']:.6g} | {row['variance_c0']:.6e} | "
            f"{row['rms_nonuniformity']:.6e} | {row['mean_offset_from_quasistatic']:.6e} | "
            f"{row['rho1']:.6f} | {row['m_eff_positive_window']:.6f} |"
        )
    lines += [
        "",
        "The key diagnostic is that the fluctuation amplitude collapses toward the exact quasistatic state as $\\omega M$ is reduced, while the normalized correlation shape and its positive-window effective count remain of order three. A normalized correlation length can therefore remain apparently system-scale even while the field it normalizes is disappearing.",
        "",
        "This **supersedes the interpretation**, not the arithmetic, of the earlier $M_{\\rm eff}^{(+)}\\approx3$ result. That number remains a valid shape diagnostic for the selected deterministic residual snapshot, but it is not evidence for a finite material statistical-cell count or a material correlation length.",
        "",
        "## Consequence for P(lambda,t)",
        "",
        "A nontrivial fatigue probability distribution cannot be obtained from the adiabatic limit of one perfectly homogeneous deterministic zero-temperature chain alone. A physically broad $P$ requires a justified ensemble source, for example finite-temperature microstates, physically specified initial-condition uncertainty, or independently justified material heterogeneity. This does not authorize an arbitrary fitted distribution.",
        "",
        "The next 1D target is therefore to separate",
        "",
        "$$",
        "P_{\\rm spatial}^{\\rm traj}(\\lambda,t)=\\frac1M\\sum_i\\delta(\\lambda-\\lambda_i(t))",
        "$$",
        "",
        "from an ensemble-averaged physical probability state and test the latter under cyclic driving.",
        "",
        "---",
        "",
        "# Normal-LJ 준정적 프로토콜 진단",
        "",
        "## 분류",
        "",
        "**정확한 정적 결과 + 통제된 수치 동역학 진단**",
        "",
        "균질한 open chain에 일정한 인장 end force $f$를 가하면 spacing 좌표에서",
        "",
        "$$",
        "\\Pi=\\sum_{i=1}^{M}[\\phi(\\lambda_i)-f\\lambda_i]",
        "$$",
        "",
        "이다. 안정 branch에서는 $\\phi''>0$이므로 정지조건 $\\phi'(\\lambda_i)=f$의 안정해가 유일하고 모든 spacing은 같은 $\\lambda_s(f)$가 된다. 따라서 zero-temperature 준정적 empirical distribution은",
        "",
        "$$",
        "P_M^{\\rm qs}(\\lambda\\mid f)=\\delta[\\lambda-\\lambda_s(f)]",
        "$$",
        "",
        "이다.",
        "",
        "Milestone 13 snapshot은 zero-mean sine load의 정수주기 cycle 2에서 저장된다. 이 정확한 위상에서 applied force는 0이므로 준정적 기준은 $\\lambda_i=1$, $C_0=0$이다. 따라서 기존 snapshot의 nonzero variance는 정적 물질 확률분포가 아니라 동적 잔류구조다.",
        "",
        "## 수치 sweep",
        "",
        "위 표와 동일한 수치에서 $\\alpha=\\omega M$을 낮추면 variance와 mean offset은 준정적 값 0으로 급격히 감소하지만, normalized correlation shape와 positive-window $M_{\\rm eff}^{(+)}$는 약 3 수준을 유지한다.",
        "",
        "따라서 기존 $M_{\\rm eff}^{(+)}\\approx3$의 **계산 자체가 틀린 것은 아니지만 해석은 수정해야 한다**. 그것은 선택한 deterministic residual snapshot의 normalized shape 진단값이지, 물질 고유 통계셀 개수나 물질 고유 correlation length의 증거가 아니다.",
        "",
        "## P(lambda,t)에 대한 결과",
        "",
        "완전히 균질한 deterministic zero-temperature chain의 adiabatic limit만으로는 폭을 가진 fatigue probability distribution이 생기지 않는다. 물리적인 broad $P$를 만들려면 finite-temperature microstate, 물리적으로 정의된 initial-condition ensemble, 또는 독립적으로 정당화된 material heterogeneity 같은 ensemble source가 필요하다. 그렇다고 임의의 분포 fitting을 허용하는 것은 아니다.",
        "",
        "다음 1D 목표는 한 trajectory의 spatial empirical measure와 실제 ensemble-averaged probability state를 명확히 분리하고, 후자를 cyclic loading에서 검증하는 것이다.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for m_count in REPRESENTED_SPACINGS:
        for alpha in OMEGA_M_VALUES:
            omega = alpha / m_count
            parameters = NormalLJParameters(
                mean_force=0.0,
                force_amplitude=FORCE_AMPLITUDE,
                omega=omega,
                ramp_cycles=2,
            )
            boundary_force = cycle_boundary_force(parameters, SAMPLE_CYCLE)
            quasistatic_stretch = stable_stretch_for_tensile_force(boundary_force)
            result = simulate_normal_lj_chain(
                parameters,
                atoms=m_count + 1,
                dt=DT,
                cycles=INTEGRATION_CYCLES,
                record_stride=10_000_000,
            )
            values = result.cycle_snapshots[SAMPLE_CYCLE]
            metrics = residual_snapshot_metrics(
                values,
                quasistatic_stretch=quasistatic_stretch,
            )
            rows.append({
                "M": m_count,
                "omega_times_M": alpha,
                "omega": omega,
                "sample_cycle": SAMPLE_CYCLE,
                "applied_force_at_exact_cycle_boundary": boundary_force,
                "quasistatic_stretch": quasistatic_stretch,
                **asdict(metrics),
            })

    for m_count in REPRESENTED_SPACINGS:
        group = [r for r in rows if r["M"] == m_count]
        baseline = group[0]
        for row in group:
            row["variance_reduction_vs_alpha_0p62"] = (
                baseline["variance_c0"] / row["variance_c0"]
                if row["variance_c0"] > 0.0 else float("inf")
            )
            row["rms_reduction_vs_alpha_0p62"] = (
                baseline["rms_nonuniformity"] / row["rms_nonuniformity"]
                if row["rms_nonuniformity"] > 0.0 else float("inf")
            )

    write_csv(DATA / "normal_lj_quasistatic_protocol.csv", rows)
    payload = {
        "classification": "EXACT STATIC RESULT + CONTROLLED NUMERICAL DYNAMICAL DIAGNOSTIC",
        "exact_static_result": "For homogeneous force-controlled 1D layer-LJ spacings on the stable branch, phi'(lambda_i)=f has one stable root, so all lambda_i=lambda_s(f) and the quasistatic empirical P is a delta distribution.",
        "snapshot_phase_correction": "The existing Milestone 13 sample_cycle=2 snapshot for zero-mean sine loading is at exact zero applied force. Its nonzero variance is therefore a dynamical residual, not a static equilibrium distribution.",
        "protocol": {
            "M_values": list(REPRESENTED_SPACINGS),
            "omega_times_M_values": list(OMEGA_M_VALUES),
            "force_amplitude": FORCE_AMPLITUDE,
            "sample_cycle": SAMPLE_CYCLE,
            "integration_cycles": INTEGRATION_CYCLES,
            "dt": DT,
        },
        "interpretation": "As omega*M is reduced, residual variance and mean offset collapse toward the exact homogeneous quasistatic state while normalized correlation shape diagnostics remain O(1) and M_eff^(+) remains near three. Therefore normalized residual correlation alone cannot define a material statistical-cell length.",
        "rows": rows,
        "next_target": "Construct and validate a physically justified 1D ensemble probability P_ens(lambda,t), distinct from a single-trajectory spatial empirical measure, under cyclic loading.",
    }
    (DATA / "normal_lj_quasistatic_protocol.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )

    plt.figure(figsize=(7.5, 5.0))
    for m_count in REPRESENTED_SPACINGS:
        group = [r for r in rows if r["M"] == m_count]
        plt.loglog(
            [r["omega_times_M"] for r in group],
            [r["variance_c0"] for r in group],
            marker="o",
            label=f"M={m_count}",
        )
    plt.xlabel("Protocol parameter alpha = omega M")
    plt.ylabel("Residual spacing variance C0")
    plt.title("Residual variance collapses toward quasistatic homogeneity")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_quasistatic_residual_variance.svg")
    plt.close()

    plt.figure(figsize=(7.5, 5.0))
    for m_count in REPRESENTED_SPACINGS:
        group = [r for r in rows if r["M"] == m_count]
        plt.semilogx(
            [r["omega_times_M"] for r in group],
            [r["m_eff_positive_window"] for r in group],
            marker="o",
            label=f"M={m_count}",
        )
    plt.xlabel("Protocol parameter alpha = omega M")
    plt.ylabel("Positive-window M_eff diagnostic")
    plt.title("Normalized correlation shape survives while amplitude collapses")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_quasistatic_meff.svg")
    plt.close()

    (REPORT / "NORMAL_LJ_QUASISTATIC_PROTOCOL.md").write_text(
        render_report(rows),
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
