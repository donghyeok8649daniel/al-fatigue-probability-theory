#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Milestone 15 준정적/protocol 교정 상태를 README, OPEN_PROBLEMS, manifest에 동기화한다."""
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
    en_anchor = "The mechanical calibration area $A_0$ is not identified with a transverse statistical independence area. The active scope remains strictly one-dimensional."
    en = r"""## Active correction — quasistatic limit of the deterministic correlation snapshot

<!-- QUASISTATIC_PROTOCOL_STATUS_EN -->

The earlier $M_{\rm eff}^{(+)}\approx3$ arithmetic remains a valid normalized-shape diagnostic for the selected deterministic snapshot, but its physical interpretation is now corrected. The committed Milestone 13 snapshot is taken at integer cycle 2 under zero-mean sinusoidal end loading, so the exact applied force at that phase is zero.

For the homogeneous force-controlled chain,

$$
\Pi=\sum_i[\phi(\lambda_i)-f\lambda_i],
\qquad
\phi'(\lambda_i)=f.
$$

On the stable branch $\phi''>0$, the root is unique and therefore every spacing is identical:

$$
\lambda_i=\lambda_s(f)\quad\forall i.
$$

At the sampled zero-force phase this gives $\lambda_i=1$ and $C_0=0$ exactly in the quasistatic state. A new $\alpha=\omega M$ sweep shows that the deterministic residual variance collapses rapidly toward zero as the drive is slowed, while the normalized correlation shape and positive-window $M_{\rm eff}^{(+)}$ remain near three. Therefore normalized residual correlation alone cannot define a material statistical-cell length.

The project now explicitly distinguishes a single-trajectory spatial empirical measure $P_{M,\mathrm{spatial}}^{\mathrm{traj}}$ from a future physically specified ensemble probability $P_{\rm ens}$. The next active target is a justified 1D initial ensemble under the same nonlinear cyclic mechanics, not a fitted named distribution.
"""
    ko_anchor = "mechanical calibration area $A_0$를 transverse statistical independence area와 동일시하지 않는다. 활성 범위는 계속 엄격한 1D다."
    ko = r"""## 활성 교정 — deterministic correlation snapshot의 준정적 극한

<!-- QUASISTATIC_PROTOCOL_STATUS_KO -->

기존 $M_{\rm eff}^{(+)}\approx3$ 계산 자체는 선택한 deterministic snapshot의 normalized-shape 진단값으로 유효하지만 물리적 해석은 교정되었다. Milestone 13 snapshot은 zero-mean sinusoidal end loading에서 정수 cycle 2에 저장되므로 그 정확한 위상에서 applied force는 0이다.

균질 force-controlled chain에서는

$$
\Pi=\sum_i[\phi(\lambda_i)-f\lambda_i],
\qquad
\phi'(\lambda_i)=f
$$

이고 안정 branch에서 $\phi''>0$이므로 안정 root가 유일하다. 따라서

$$
\lambda_i=\lambda_s(f)\quad\forall i
$$

이다.

현재 snapshot의 zero-force 위상에서는 준정적 상태가 정확히 $\lambda_i=1$, $C_0=0$이다. 새 $\alpha=\omega M$ sweep에서는 drive를 느리게 할수록 deterministic residual variance가 0으로 급격히 감소하지만 normalized correlation shape와 positive-window $M_{\rm eff}^{(+)}$는 약 3을 유지한다. 따라서 normalized residual correlation만으로 물질 고유 statistical-cell length를 정의할 수 없다.

이제 한 trajectory의 spatial empirical measure $P_{M,\mathrm{spatial}}^{\mathrm{traj}}$와 물리적으로 정의해야 할 ensemble probability $P_{\rm ens}$를 명시적으로 구분한다. 다음 활성 목표는 임의 named distribution fitting이 아니라 같은 nonlinear cyclic mechanics 위에서 물리적으로 정당한 1D initial ensemble을 만드는 것이다.
"""
    text = insert_after_once(text, en_anchor, "<!-- QUASISTATIC_PROTOCOL_STATUS_EN -->", en)
    text = insert_after_once(text, ko_anchor, "<!-- QUASISTATIC_PROTOCOL_STATUS_KO -->", ko)
    path.write_text(text, encoding="utf-8")


def update_open_problems() -> None:
    path = ROOT / "docs" / "OPEN_PROBLEMS.md"
    text = path.read_text(encoding="utf-8")
    en_anchor = "4. keep $A_0$ (energy calibration area) separate from any future transverse statistical area."
    en = r"""## Quasistatic-protocol correction and new probability target

<!-- QUASISTATIC_PROTOCOL_OPEN_EN -->

Milestone 15 identifies the committed deterministic cycle-2 correlation snapshot as a zero-applied-force residual-dynamics state. In the exact homogeneous quasistatic force-controlled solution all spacings are identical, so the spatial variance is zero. Slowing the drive collapses the residual variance while leaving the normalized correlation shape approximately intact.

Open questions are therefore updated:

1. define a physically justified 1D initial phase-space ensemble $\Gamma_0$;
2. evolve the ensemble with the same full nonlinear conservative layer-LJ dynamics under cyclic loading;
3. distinguish ensemble broadening from single-trajectory residual waves;
4. test one-point $P_{\rm ens}$, pair dependence, and first-passage/tail convergence with system size and ensemble size;
5. only after those tests define an event-relevant statistical-cell or clustering scale.
"""
    ko_anchor = "4. $A_0$(energy calibration area)와 미래의 transverse statistical area를 분리한다."
    ko = r"""## 준정적 프로토콜 교정과 새로운 확률 목표

<!-- QUASISTATIC_PROTOCOL_OPEN_KO -->

Milestone 15에서는 기존 deterministic cycle-2 correlation snapshot이 applied force가 0인 동적 잔류상태임을 확인했다. 정확한 균질 force-controlled 준정적 해에서는 모든 spacing이 동일하므로 spatial variance는 0이다. drive를 느리게 하면 residual variance는 0으로 무너지지만 normalized correlation shape는 거의 유지된다.

따라서 미해결 문제를 다음처럼 갱신한다.

1. 물리적으로 정당한 1D initial phase-space ensemble $\Gamma_0$를 정의한다.
2. 같은 full nonlinear conservative layer-LJ dynamics와 cyclic loading으로 ensemble을 진화시킨다.
3. ensemble broadening과 single-trajectory residual wave를 분리한다.
4. one-point $P_{\rm ens}$, pair dependence, first-passage/tail의 system-size 및 ensemble-size 수렴을 검사한다.
5. 그 검증 뒤에만 event-relevant statistical-cell 또는 clustering scale을 정의한다.
"""
    text = insert_after_once(text, en_anchor, "<!-- QUASISTATIC_PROTOCOL_OPEN_EN -->", en)
    text = insert_after_once(text, ko_anchor, "<!-- QUASISTATIC_PROTOCOL_OPEN_KO -->", ko)
    path.write_text(text, encoding="utf-8")


def update_manifest() -> None:
    path = ROOT / "results" / "data" / "result_manifest.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["state"] = (
        "Exact nonlinear transport remains active, but Milestone 15 corrects the interpretation of the deterministic spatial-correlation snapshot. "
        "For a homogeneous force-controlled 1D layer-LJ chain the exact stable quasistatic state has lambda_i=lambda_s(f) for all i. "
        "The existing cycle-2 zero-mean-sine snapshot is taken at zero applied force, so its nonzero variance is residual dynamics. "
        "A drive-rate sweep shows residual variance collapsing toward zero while normalized M_eff^(+) stays near three. "
        "Single-trajectory spatial empirical P is therefore separated from a future physically specified ensemble probability P_ens."
    )
    active = data.setdefault("active_files", {})
    active["quasistatic_protocol_theory"] = "theory/normal_lj_quasistatic_protocol.py"
    active["quasistatic_protocol_test"] = "tests/test_normal_lj_quasistatic_protocol.py"
    active["quasistatic_protocol_derivation"] = "docs/MILESTONE15_QUASISTATIC_PROTOCOL_COHERENCE.md"
    active["quasistatic_protocol_variables"] = "docs/VARIABLE_DEFINITIONS_QUASISTATIC_PROTOCOL.md"
    active["quasistatic_protocol_simulation"] = "simulations/run_normal_lj_quasistatic_protocol.py"
    active["quasistatic_protocol_data"] = "results/data/normal_lj_quasistatic_protocol.json"
    active["quasistatic_protocol_table"] = "results/data/normal_lj_quasistatic_protocol.csv"
    active["quasistatic_protocol_report"] = "results/reports/NORMAL_LJ_QUASISTATIC_PROTOCOL.md"
    data["probability_semantics"] = {
        "single_trajectory_spatial_empirical": "P_M_spatial_traj(lambda,t)=(1/M) sum_i delta(lambda-lambda_i(t)); exact trajectory measure but not by itself a physical randomness model.",
        "ensemble_target": "P_ens(lambda,t)=E_{Gamma0}[(1/M) sum_i delta(lambda-lambda_i(t;Gamma0))]; requires a physically specified initial phase-space ensemble before physical fatigue interpretation."
    }
    data["current_result"] = (
        "Milestone 15 proves that the homogeneous stable force-controlled 1D chain has an exactly uniform quasistatic spacing state. "
        "The prior Milestone 13 cycle-2 zero-mean-sine correlation snapshot occurs at exact zero applied force and is therefore a residual-dynamics snapshot. "
        "In the new M=31,63 sweep, reducing alpha=omega*M from 0.62 to 0.0775 drives C0 and |mean(lambda)-1| strongly toward zero, while the positive-window normalized M_eff diagnostic remains near three. "
        "Thus M_eff~3 is a residual normalized-shape diagnostic, not evidence for a material local correlation length."
    )
    data["next_target"] = (
        "Stay strictly in 1D. Define a physically justified initial phase-space ensemble Gamma0, evolve it under the same full nonlinear cyclic layer-LJ mechanics, "
        "and compare the resulting ensemble-averaged P_ens(lambda,t) with single-trajectory spatial empirical P. Test ensemble-size and system-size convergence of one-point shape, pair dependence, and tail/first-passage statistics before defining any event-relevant statistical-cell length."
    )
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    update_readme()
    update_open_problems()
    update_manifest()
    print("Milestone 15 quasistatic/protocol status synchronized.")


if __name__ == "__main__":
    main()
