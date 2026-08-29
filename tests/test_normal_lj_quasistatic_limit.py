# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D normal layer-LJ의 준정적 안정 branch가 단일값이며 닫힌 인장 cycle에서 경로를 되짚는지 검증한다.
# - 주요 클래스: 없음
# - 주요 함수: test_stable_branch_endpoints, test_quasistatic_cycle_is_closed, test_frequency_is_absent_from_phase_path, test_rejects_unstable_or_compressive_cycles
# - 주의: 이 테스트는 준정적 homogeneous equilibrium branch의 성질만 검증하며 실제 피로 누적 법칙을 가정하지 않는다.
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
