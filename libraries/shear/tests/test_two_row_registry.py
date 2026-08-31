# 이 파일은 비활성 두 원자열 registry 수학의 convention, 합, 주기성, 미분을 검증한다.
# 주요 검사는 q=2 닫힌형, well-depth 변환, registry 최소점의 거리 의존성이다.
# 이 테스트는 archive 검산이며 활성 normal tensile solver의 물리 검증을 대신하지 않는다.

from __future__ import annotations

import math

import pytest

from libraries.shear.theory.two_row_registry import (
    inverse_power_direct,
    inverse_power_q2_closed,
    pair_equilibrium_ratio,
    pair_potential_coefficient,
    registry_force_per_repeat_direct,
    two_row_cross_energy_per_repeat_direct,
    well_depth_scale,
)


def test_well_depth_and_coefficient_conventions_are_exactly_convertible() -> None:
    m, n = 12.19, 6.0
    epsilon_well = 2.3
    epsilon_coefficient = well_depth_scale(m, n) * epsilon_well
    r_e = pair_equilibrium_ratio(m, n)
    value = pair_potential_coefficient(
        r_e,
        epsilon_coefficient,
        1.0,
        m,
        n,
    )
    assert value == pytest.approx(-epsilon_well, rel=2.0e-14)


@pytest.mark.parametrize("delta,eta", [(0.0, 0.7), (0.23, 0.7), (0.5, 1.3)])
def test_shifted_q2_sum_converges_to_independent_closed_form(
    delta: float,
    eta: float,
) -> None:
    exact = inverse_power_q2_closed(delta, eta)
    coarse = abs(inverse_power_direct(2.0, delta, eta, 100) - exact)
    fine = abs(inverse_power_direct(2.0, delta, eta, 1000) - exact)
    assert fine < 0.101 * coarse


def test_exact_q2_form_is_periodic_and_even() -> None:
    delta, eta = 0.217, 0.83
    reference = inverse_power_q2_closed(delta, eta)
    assert inverse_power_q2_closed(delta + 1.0, eta) == pytest.approx(reference)
    assert inverse_power_q2_closed(-delta, eta) == pytest.approx(reference)


def test_registry_force_is_energy_derivative() -> None:
    args = dict(
        a=0.93,
        b=1.0,
        epsilon_coefficient=1.0,
        sigma=1.0,
        m=12.19,
        n=6.0,
        half_width=2000,
    )
    s = 0.19
    step = 1.0e-6
    plus = two_row_cross_energy_per_repeat_direct(s=s + step, **args)
    minus = two_row_cross_energy_per_repeat_direct(s=s - step, **args)
    finite_difference = (plus - minus) / (2.0 * step)
    exact_derivative = registry_force_per_repeat_direct(s=s, **args)
    assert finite_difference == pytest.approx(exact_derivative, rel=2.0e-8)


def test_preferred_registry_depends_on_normal_separation() -> None:
    common = dict(
        b=1.0,
        epsilon_coefficient=1.0,
        sigma=1.0,
        m=12.19,
        n=6.0,
        half_width=2000,
    )
    close_aligned = two_row_cross_energy_per_repeat_direct(a=1.0, s=0.0, **common)
    close_offset = two_row_cross_energy_per_repeat_direct(a=1.0, s=0.5, **common)
    far_aligned = two_row_cross_energy_per_repeat_direct(a=1.5, s=0.0, **common)
    far_offset = two_row_cross_energy_per_repeat_direct(a=1.5, s=0.5, **common)
    assert close_offset < close_aligned
    assert far_aligned < far_offset


def test_invalid_inputs_are_rejected() -> None:
    with pytest.raises(ValueError):
        well_depth_scale(6.0, 6.0)
    with pytest.raises(ValueError):
        inverse_power_direct(1.0, 0.0, 1.0, 10)
    with pytest.raises(ValueError):
        inverse_power_q2_closed(0.0, 0.0)

