# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: test_stable_branch_endpoints, test_quasistatic_cycle_is_closed
#   test_frequency_is_absent_from_phase_path, test_rejects_unstable_or_compressive_cycles
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Regression tests for the exact quasistatic homogeneous LJ branch."""

import math

import numpy as np
import pytest

from theory.normal_lj_chain import critical_dimensionless_force, critical_stretch
from theory.normal_lj_quasistatic_limit import (
    branch_retrace_error,
    closed_cycle_residual,
    quasistatic_cycle,
    stable_stretch_from_force,
)


def test_stable_branch_endpoints():
    fc = critical_dimensionless_force()
    assert stable_stretch_from_force(0.0) == pytest.approx(1.0)
    assert stable_stretch_from_force(fc) == pytest.approx(critical_stretch())


def test_quasistatic_cycle_is_closed():
    residual = closed_cycle_residual(mean_force=0.02, force_amplitude=0.01)
    assert residual < 1.0e-12
    assert branch_retrace_error(mean_force=0.02, force_amplitude=0.01) < 1.0e-12


def test_frequency_is_absent_from_phase_path():
    phase, force, stretch = quasistatic_cycle(
        mean_force=0.02,
        force_amplitude=0.01,
        samples=361,
    )
    # A change of laboratory frequency changes t=theta/(2*pi*f), not the
    # equilibrium path as a function of phase theta.
    for frequency_hz in (1.0, 20.0, 100.0):
        time = phase / (2.0 * math.pi * frequency_hz)
        recovered_phase = 2.0 * math.pi * frequency_hz * time
        assert np.allclose(recovered_phase, phase)
        _, force_again, stretch_again = quasistatic_cycle(
            mean_force=0.02,
            force_amplitude=0.01,
            samples=361,
        )
        assert np.allclose(force_again, force)
        assert np.allclose(stretch_again, stretch)


def test_rejects_unstable_or_compressive_cycles():
    fc = critical_dimensionless_force()
    with pytest.raises(ValueError):
        quasistatic_cycle(mean_force=0.005, force_amplitude=0.01)
    with pytest.raises(ValueError):
        quasistatic_cycle(mean_force=fc, force_amplitude=0.001)
