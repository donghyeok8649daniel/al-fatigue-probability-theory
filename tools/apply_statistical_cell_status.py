#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Milestone 13/13A의 1D 통계셀 상태를 README/OPEN_PROBLEMS/manifest에 동기화한다."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def insert_after_once(text: str, anchor: str, marker: str, block: str) -> str:
    if marker in text:
        return text
    if anchor not in text:
        raise RuntimeError(f"anchor not found: {anchor!r}")
    return text.replace(anchor, anchor + "\n\n" + block.rstrip(), 1)


def update_readme() -> None:
    path = ROOT / "README.md"
    text = path.read_text(encoding="utf-8")
    en = r"""## Active 1D statistical-cell dependence scale

<!-- STATISTICAL_CELL_STATUS_EN -->

Probability aggregation is now explicitly separated into complete identical dependence, partial dependence, and true independence. For a second-order stationary 1D spacing process, the exact variance identity is

$$
\operatorname{Var}(\bar\lambda_M)
=\frac{\sigma_\lambda^2}{M}\tau_M,
\qquad
\tau_M=1+2\sum_{k=1}^{M-1}\left(1-\frac{k}{M}\right)\rho_k.
$$

This defines the variance-equivalent independent count and axial statistical length

$$
M_{\rm eff}=\frac{M}{\tau_M},
\qquad
\ell_{\rm stat}^{(2)}=a_0\tau_M.
$$

A single deterministic finite snapshot cannot use the all-lag population formula directly because sample-mean centering gives an exact weighted zero-sum identity. Such snapshots therefore use a separately labeled first-positive-lobe estimator. In the dynamically matched sweep $M=31,63,127,255$, the corrected estimate gives $M_{\rm eff}^{(+)}\approx2.93,2.99,3.03,3.05$ while $\ell_{\rm stat}^{(2,+)}/a_0\approx10.58,21.10,41.93,83.49$. Thus the tested protocol retains system-scale coherence rather than converging to a local material correlation length.

The mechanical calibration area $A_0$ is not identified with a transverse statistical independence area. The active scope remains strictly one-dimensional.
"""
    ko = r"""## 활성 1D 통계셀 종속성 척도

<!-- STATISTICAL_CELL_STATUS_KO -->

확률을 합칠 때 완전 동일 종속, 부분 종속, 진짜 독립을 명시적으로 구분한다. second-order stationary 1D spacing process에서는

$$
\operatorname{Var}(\bar\lambda_M)
=\frac{\sigma_\lambda^2}{M}\tau_M,
\qquad
\tau_M=1+2\sum_{k=1}^{M-1}\left(1-\frac{k}{M}\right)\rho_k
$$

가 정확하고, 이에 따라

$$
M_{\rm eff}=\frac{M}{\tau_M},
\qquad
\ell_{\rm stat}^{(2)}=a_0\tau_M
$$

를 정의한다.

하지만 하나의 deterministic finite snapshot을 자기 sample mean으로 center하면 모든 lag를 넣은 weighted correlation sum이 정확히 0이 되는 finite-sample identity가 있으므로 population 식을 그대로 plug-in하면 안 된다. 따라서 snapshot에는 별도로 표시한 first-positive-lobe estimator를 사용한다. Dynamically matched $M=31,63,127,255$ sweep에서는 corrected estimate가 $M_{\rm eff}^{(+)}\approx2.93,2.99,3.03,3.05$이고 $\ell_{\rm stat}^{(2,+)}/a_0\approx10.58,21.10,41.93,83.49$이다. 즉 tested protocol은 local material correlation length로 수렴하지 않고 system-scale coherence를 유지한다.

mechanical calibration area $A_0$를 transverse statistical independence area와 동일시하지 않는다. 활성 범위는 계속 엄격한 1D다.
"""
    text = insert_after_once(text, "# Al Fatigue Probability Theory", "<!-- STATISTICAL_CELL_STATUS_EN -->", en)
    text = insert_after_once(text, "# 한국어 번역", "<!-- STATISTICAL_CELL_STATUS_KO -->", ko)
    path.write_text(text, encoding="utf-8")


def update_open_problems() -> None:
    path = ROOT / "docs" / "OPEN_PROBLEMS.md"
    text = path.read_text(encoding="utf-8")
    en = r"""## 1D statistical-cell and dependence priority

<!-- STATISTICAL_CELL_OPEN_EN -->

The current dynamically matched chain behaves, in a variance-equivalent sense, like roughly three independent axial probability blocks even as $M$ increases from 31 to 255. The estimated axial statistical length grows with system size instead of converging.

Open questions:

1. determine whether this system-scale coherence is caused by the present boundary/loading protocol or survives in a physically appropriate 1D ensemble;
2. do not assign a fixed local mini-cell length until a size-independent axial scale converges;
3. keep the variance-equivalent length separate from a crack-tail/first-passage clustering length;
4. keep $A_0$ (energy calibration area) separate from any future transverse statistical area.
"""
    ko = r"""## 1D 통계셀 및 종속성 우선문제

<!-- STATISTICAL_CELL_OPEN_KO -->

현재 dynamically matched chain은 variance-equivalent 관점에서 $M=31$에서 255로 커져도 대략 3개의 독립 axial probability block처럼 행동한다. 추정 axial statistical length는 수렴하지 않고 system size와 함께 증가한다.

미해결 문제는 다음과 같다.

1. 이 system-scale coherence가 현재 boundary/loading protocol 때문인지, 물리적으로 적절한 1D ensemble에서도 유지되는지 확인한다.
2. system-size-independent axial scale이 수렴하기 전에는 fixed local mini-cell length를 지정하지 않는다.
3. variance-equivalent length와 crack-tail/first-passage clustering length를 분리한다.
4. $A_0$(energy calibration area)와 미래의 transverse statistical area를 분리한다.
"""
    text = insert_after_once(text, "# Open Problems — Active 1D Layer-LJ Mainline", "<!-- STATISTICAL_CELL_OPEN_EN -->", en)
    text = insert_after_once(text, "# 한국어 번역 — 활성 1D Layer-LJ Mainline 미해결 문제", "<!-- STATISTICAL_CELL_OPEN_KO -->", ko)
    path.write_text(text, encoding="utf-8")


def update_manifest() -> None:
    path = ROOT / "results" / "data" / "result_manifest.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["active_mainline"] = (
        "Strictly 1D normal layer-LJ in continuous physical time t, with full nonlinear mechanics, "
        "exact transport, explicitly labeled statistical-mechanical limits, and correlation-derived axial statistical-cell scales"
    )
    data["scope_freeze"] = (
        "Remain strictly in one-dimensional normal tension until P(lambda,t) is physically identified and validated. "
        "2D/3D FEM, multiaxial criteria, specimen hazard integration, and 3D mini-mesh coupling are deferred."
    )
    data["state"] = (
        "P(lambda,t) has exact transport and physical equilibrium constraints. Probability aggregation distinguishes "
        "complete identical dependence, partial dependence, and true independence. The exact variance-equivalent "
        "statistical length belongs to the true stationary correlation sequence; deterministic finite snapshots use "
        "a separately labeled first-positive-lobe estimator because the all-lag sample-mean-centered plug-in sum has an exact zero-sum artifact."
    )
    active = data.setdefault("active_files", {})
    active["statistical_cell_theory"] = "theory/normal_lj_statistical_cells.py"
    active["statistical_cell_test"] = "tests/test_normal_lj_statistical_cells.py"
    active["statistical_cell_derivation"] = "docs/MILESTONE13_1D_STATISTICAL_CELL.md"
    active["statistical_cell_finite_snapshot_correction"] = "docs/MILESTONE13A_FINITE_SNAPSHOT_CORRECTION.md"
    active["statistical_cell_variables"] = "docs/VARIABLE_DEFINITIONS_STATISTICAL_CELL.md"
    active["statistical_cell_report"] = "results/reports/NORMAL_LJ_STATISTICAL_CELL.md"
    data["current_result"] = (
        "For true stationary 1D correlations, Var(mean_M)=sigma_lambda^2*tau_M/M exactly and "
        "M_eff=M/tau_M, ell_stat_2=a0*tau_M. For one deterministic sample-mean-centered snapshot the all-lag "
        "weighted covariance sum is exactly zero, so a first-positive-lobe estimator is used instead. In the "
        "dynamically matched M=31,63,127,255 sweep, M_eff_hat_plus is approximately 2.929, 2.986, 3.029, 3.054, "
        "while ell_stat_2_hat_plus/a0 is approximately 10.584, 21.101, 41.934, 83.489. Thus the tested chain retains "
        "system-scale coherent dependence and does not yet provide a converged local material correlation length. "
        "Complete identity remains distinct from partial dependence, and event independence still requires joint factorization."
    )
    data["next_target"] = (
        "Stay in 1D. Determine whether the observed system-scale coherence is caused by the present boundary/loading "
        "protocol or persists in a physically appropriate 1D ensemble. Continue identifying P(lambda,t), and later "
        "derive an event/first-passage clustering length separately from the second-moment statistical length. Do not "
        "assign a fixed local mini-cell or identify A0 with a statistical independence area without an independent physical derivation."
    )
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    update_readme()
    update_open_problems()
    update_manifest()
    print("Milestone 13 statistical-cell status synchronized.")


if __name__ == "__main__":
    main()
