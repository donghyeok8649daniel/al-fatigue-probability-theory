#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Milestone 12의 물리 통계역학 P 상태를 README/문서/manifest에 동기화한다."""
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
    en = """## Active physical-statistical constraint on $P$\n\n<!-- PHYSICAL_P_STATUS_EN -->\n\nThe current mechanics now gives a physically grounded hierarchy of possible spacing distributions without fitting a named family. With\n\n$$\nU(a)=E_0\phi(a/a_0)+C,\n$$\n\nthe elastic calibration $E=(a_0/A_0)U''(a_0)$ and $\phi''(1)=1$ imply\n\n$$\n\boxed{E_0=EA_0a_0},\n\qquad\n\boxed{\chi=\frac{EA_0a_0}{k_BT}}.\n$$\n\nAt zero-temperature homogeneous quasistatic equilibrium, the distribution is not broad:\n\n$$\n\boxed{P(\lambda\mid f)=\delta[\lambda-\lambda_s(f)]}.\n$$\n\nAt fixed total normalized length, canonical equilibrium gives the exact finite-$M$ marginal\n\n$$\n\boxed{\nP_M(\lambda\mid L,\chi)\n=\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}.\n}\n$$\n\nFor constant tensile force $f>0$, the full-domain Gibbs integral diverges because $\phi(\lambda)-f\lambda\to-\infty$. For $0<f<f_c$, a Gibbs density conditioned on the intact basin $0<\lambda<\lambda_b(f)$ is therefore only a **controlled metastable/local-equilibrium approximation**, not a global equilibrium law or a fatigue-life law. See `docs/MILESTONE12_PHYSICAL_STATISTICAL_P.md`.\n\nThe representative layer area $A_0$ is now an explicit physical bottleneck: no numerical aluminum thermal distribution is claimed until $A_0$ is defined consistently with the coarse-grained layer interaction.\n"""
    ko = """## $P$에 대한 활성 물리 통계역학 제약\n\n<!-- PHYSICAL_P_STATUS_KO -->\n\n현재 mechanics에서는 named distribution을 fitting하지 않고도 가능한 spacing distribution의 물리적 hierarchy를 얻는다.\n\n$$\nU(a)=E_0\phi(a/a_0)+C\n$$\n\n로 두면 elastic calibration $E=(a_0/A_0)U''(a_0)$와 $\phi''(1)=1$에서\n\n$$\n\boxed{E_0=EA_0a_0},\n\qquad\n\boxed{\chi=\frac{EA_0a_0}{k_BT}}\n$$\n\n가 나온다.\n\nzero-temperature homogeneous quasistatic equilibrium에서는 broad distribution이 아니라\n\n$$\n\boxed{P(\lambda\mid f)=\delta[\lambda-\lambda_s(f)]}\n$$\n\n이다.\n\nfixed total normalized length의 canonical equilibrium에서는 exact finite-$M$ marginal\n\n$$\n\boxed{\nP_M(\lambda\mid L,\chi)\n=\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}\n}\n$$\n\n을 얻는다.\n\nconstant tensile force $f>0$에서는 $\phi(\lambda)-f\lambda\to-\infty$이므로 full-domain Gibbs integral이 발산한다. 따라서 $0<f<f_c$에서 intact basin $0<\lambda<\lambda_b(f)$에 조건부로 둔 Gibbs density는 **controlled metastable/local-equilibrium approximation**일 뿐 global equilibrium law나 fatigue-life law가 아니다. `docs/MILESTONE12_PHYSICAL_STATISTICAL_P.md`를 현재 물리 통계역학 기준으로 사용한다.\n\nrepresentative layer area $A_0$가 이제 명시적인 physical bottleneck이다. coarse-grained layer interaction과 일관된 $A_0$가 정해지기 전에는 numerical aluminum thermal distribution을 주장하지 않는다.\n"""
    text = insert_after_once(text, "# Al Fatigue Probability Theory", "<!-- PHYSICAL_P_STATUS_EN -->", en)
    text = insert_after_once(text, "# 한국어 번역", "<!-- PHYSICAL_P_STATUS_KO -->", ko)
    path.write_text(text, encoding="utf-8")


def update_open_problems() -> None:
    path = ROOT / "docs" / "OPEN_PROBLEMS.md"
    text = path.read_text(encoding="utf-8")
    en = """## Current physical priority — identify the ensemble before closing $P$\n\n<!-- PHYSICAL_P_OPEN_EN -->\n\nMilestone 12 shows that established physical theory already separates the candidate forms of $P$. Zero-temperature stable equilibrium gives a delta measure; fixed-length thermal equilibrium gives an exact canonical marginal; constant tensile force has no normalizable full-domain Gibbs equilibrium and therefore requires either a length constraint or a metastable intact-basin interpretation.\n\nThe highest-priority unresolved physical questions are now:\n\n1. define the representative layer-patch area $A_0$ consistently with the calibrated effective layer potential, because $E_0=EA_0a_0$ and $\chi=E_0/(k_BT)$;\n2. determine whether the driven conservative chain samples an approximately microcanonical state, a local-equilibrium state, or a strongly nonequilibrium coherent state;\n3. compare the physically derived fixed-length/microcanonical/metastable forms directly with deterministic spacing statistics before adding any further closure;\n4. if equilibrium forms fail, derive the nonequilibrium correction from the exact $F_1/P_2$ transport hierarchy rather than fitting a distribution.\n"""
    ko = """## 현재 물리 우선순위 — $P$를 닫기 전에 ensemble부터 결정\n\n<!-- PHYSICAL_P_OPEN_KO -->\n\nMilestone 12에서 기존 물리이론만으로도 $P$ 후보 함수형이 regime별로 갈린다는 것이 드러났다. zero-temperature stable equilibrium은 delta measure를 주고, fixed-length thermal equilibrium은 exact canonical marginal을 주며, constant tensile force에는 normalizable full-domain Gibbs equilibrium이 없다. 따라서 length constraint 또는 metastable intact-basin 해석이 필요하다.\n\n현재 최우선 미해결 물리문제는 다음이다.\n\n1. calibrated effective layer potential과 일관되게 representative layer-patch area $A_0$를 정의한다. $E_0=EA_0a_0$, $\chi=E_0/(k_BT)$이기 때문이다.\n2. driven conservative chain이 approximately microcanonical state인지, local-equilibrium state인지, strongly nonequilibrium coherent state인지 판별한다.\n3. 추가 closure를 넣기 전에 physically derived fixed-length/microcanonical/metastable form을 deterministic spacing statistics와 직접 비교한다.\n4. equilibrium form이 실패하면 distribution fitting 대신 exact $F_1/P_2$ transport hierarchy에서 nonequilibrium correction을 유도한다.\n"""
    text = insert_after_once(text, "# Open Problems — Active 1D Layer-LJ Mainline", "<!-- PHYSICAL_P_OPEN_EN -->", en)
    text = insert_after_once(text, "# 한국어 번역 — 활성 1D Layer-LJ Mainline 미해결 문제", "<!-- PHYSICAL_P_OPEN_KO -->", ko)
    path.write_text(text, encoding="utf-8")


def update_assumptions() -> None:
    path = ROOT / "docs" / "ASSUMPTIONS.md"
    text = path.read_text(encoding="utf-8")
    en = """## Statistical-mechanical ensemble rule\n\n<!-- PHYSICAL_P_ASSUMPTIONS_EN -->\n\nA smooth thermal $P$ is not assumed merely because the coordinate is treated probabilistically. The following physical distinctions are mandatory:\n\n- the zero-temperature homogeneous quasistatic state is a delta distribution on the stable mechanical branch;\n- a microcanonical distribution requires an isolated conservative equilibrium interpretation and includes the kinetic density-of-states factor after momenta are integrated out;\n- a canonical fixed-length distribution requires a thermal reservoir and fixed total length;\n- a tensile intact-basin Gibbs density is only a metastable/local-equilibrium approximation and requires intrabasin equilibration to be fast relative to loading/escape;\n- no Kramers or Arrhenius escape prefactor is introduced without an independently derived bath/friction/phonon time scale;\n- $A_0$ is a physical coarse-graining input and may not be tuned to obtain a desired tail probability or fatigue life.\n"""
    ko = """## 통계역학 ensemble 규칙\n\n<!-- PHYSICAL_P_ASSUMPTIONS_KO -->\n\ncoordinate를 probabilistic하게 다룬다는 이유만으로 smooth thermal $P$를 가정하지 않는다. 다음 물리적 구분을 반드시 유지한다.\n\n- zero-temperature homogeneous quasistatic state는 stable mechanical branch 위의 delta distribution이다.\n- microcanonical distribution은 isolated conservative equilibrium 해석을 요구하며 momentum 적분 뒤 kinetic density-of-states factor를 포함한다.\n- canonical fixed-length distribution은 thermal reservoir와 fixed total length를 요구한다.\n- tensile intact-basin Gibbs density는 metastable/local-equilibrium approximation일 뿐이며 intrabasin equilibration이 loading/escape보다 빨라야 한다.\n- independently derived bath/friction/phonon time scale 없이 Kramers 또는 Arrhenius escape prefactor를 넣지 않는다.\n- $A_0$는 physical coarse-graining input이며 원하는 tail probability나 fatigue life를 얻도록 tuning하면 안 된다.\n"""
    text = insert_after_once(text, "# Assumptions and Approximations", "<!-- PHYSICAL_P_ASSUMPTIONS_EN -->", en)
    text = insert_after_once(text, "# 한국어 번역 — 가정과 근사", "<!-- PHYSICAL_P_ASSUMPTIONS_KO -->", ko)
    path.write_text(text, encoding="utf-8")


def update_manifest() -> None:
    path = ROOT / "results" / "data" / "result_manifest.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["active_mainline"] = (
        "1D normal layer-LJ in continuous physical time t, with full nonlinear mechanics, "
        "exact transport, and explicitly labeled statistical-mechanical ensemble limits"
    )
    data["state"] = (
        "P(lambda,t) has exact driven-transport constraints and physically derived equilibrium limits: "
        "T=0 quasistatic delta state, fixed-length canonical/microcanonical marginals, and a controlled "
        "metastable intact-basin Gibbs form under subcritical tension"
    )
    active = data.setdefault("active_files", {})
    active["physical_distribution_theory"] = "theory/normal_lj_physical_distribution.py"
    active["physical_distribution_test"] = "tests/test_normal_lj_physical_distribution.py"
    active["physical_distribution_simulation"] = "simulations/run_normal_lj_physical_distribution.py"
    active["physical_distribution_derivation"] = "docs/MILESTONE12_PHYSICAL_STATISTICAL_P.md"
    active["physical_distribution_variables"] = "docs/VARIABLE_DEFINITIONS_PHYSICAL_P.md"
    active["physical_distribution_data"] = "results/data/normal_lj_physical_distribution.json"
    active["physical_distribution_table"] = "results/data/normal_lj_physical_distribution.csv"
    active["physical_distribution_report"] = "results/reports/NORMAL_LJ_PHYSICAL_DISTRIBUTION.md"
    data["current_result"] = (
        "Physical theory now fixes several non-arbitrary forms of P. The elastic calibration gives "
        "E0=E*A0*a0 and chi=E0/(k_B*T). Zero-temperature homogeneous quasistatic equilibrium is a delta "
        "measure at the stable force-balance spacing. An isolated equilibrium chain has a microcanonical "
        "configurational measure weighted by remaining kinetic phase volume. At fixed total length, the exact "
        "canonical one-spacing marginal is exp[-chi*phi(lambda)]*Z_{M-1}(L-lambda)/Z_M(L). A global tensile "
        "force-controlled Gibbs density does not normalize; for 0<f<fc an intact-basin Gibbs density is only a "
        "controlled metastable/local-equilibrium approximation."
    )
    data["next_target"] = (
        "Physically define the representative layer area A0, then test whether deterministic cyclic states are "
        "consistent with microcanonical/local-equilibrium statistics. Compare physically derived P forms directly "
        "with mechanics; if they fail, derive the nonequilibrium correction from the exact F1/P2 hierarchy."
    )
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    update_readme()
    update_open_problems()
    update_assumptions()
    update_manifest()
    print("Milestone 12 physical-statistical status synchronized.")


if __name__ == "__main__":
    main()
