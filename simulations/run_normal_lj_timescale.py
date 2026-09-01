# === 한국어 파일 안내 시작 ===
# - 파일 역할: 물리 주파수와 reduced-model time scale의 대응을 계산해 reference 결과를 생성한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Reproduce the normal-LJ conservative time-scale falsification calculation."""
from __future__ import annotations

import json
from pathlib import Path

from theory.normal_lj_chain import atomic_time_scale, critical_stretch, normalized_lj_stiffness
from theory.normal_lj_timescale import (
    homogeneous_stretch_for_dimensionless_stress,
    local_small_oscillation_frequency_hz,
    moving_atoms_for_target_frequency,
    near_critical_distance_for_target_local_frequency,
    normalized_lj_third_derivative,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "data"

MASS_AL = 26.9815385 * 1.66053906660e-27
A0_M = 2.8627442948e-10
E_PA = 69.0e9
REFERENCE_AREA_M2 = 6.0338e-20


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    t0 = atomic_time_scale(MASS_AL, A0_M, E_PA, REFERENCE_AREA_M2)
    target_hz = 20.0
    atoms = moving_atoms_for_target_frequency(target_hz, t0)
    length_m = atoms * A0_M

    lam100 = homogeneous_stretch_for_dimensionless_stress(100.0e6 / E_PA)
    k100 = float(normalized_lj_stiffness(lam100))
    f100 = local_small_oscillation_frequency_hz(lam100, t0)

    lam_c = critical_stretch()
    delta = near_critical_distance_for_target_local_frequency(target_hz, t0)

    result = {
        "classification": "falsification / time-scale diagnostic",
        "target_frequency_hz": target_hz,
        "atomic_time_scale_s": t0,
        "implied_lattice_speed_m_per_s": A0_M / t0,
        "moving_atoms_for_lowest_mode_at_20Hz": atoms,
        "chain_length_for_lowest_mode_at_20Hz_m": length_m,
        "lambda_c": lam_c,
        "phi3_at_lambda_c": normalized_lj_third_derivative(lam_c),
        "required_distance_from_lambda_c_for_20Hz_local_softening": delta,
        "lambda_at_100MPa": lam100,
        "tangent_stiffness_at_100MPa": k100,
        "local_frequency_at_100MPa_hz": f100,
    }

    (OUT / "normal_lj_timescale_summary.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
