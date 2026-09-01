# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: test_active_bessel_identity_matches_independent_direct_sum
#   test_active_bessel_energy_is_periodic_and_its_force_is_exact_derivative
#   test_preferred_registry_is_computed_not_assumed, test_signed_schmid_projection_for_fcc_example
#   _transport_config, test_zero_load_preserves_symmetry_and_probability
#   test_discrete_gibbs_state_has_zero_registry_current
#   test_resolved_shear_pulse_leaves_an_unwrapped_residual_slip_population
#   test_symmetric_zero_mean_cycles_do_not_create_directed_plastic_drift
#   test_registry_grid_and_timestep_refinement
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
from __future__ import annotations

import math

import numpy as np
import pytest

from libraries.shear.theory.two_row_registry import (
    inverse_power_direct,
    inverse_power_q2_closed,
    registry_force_per_repeat_direct,
)
from theory.registry_lattice import (
    RegistryLattice,
    preferred_registry,
    registry_energy,
    registry_energy_derivative,
    schmid_factor,
    shifted_inverse_power_bessel,
    h_q_delta_derivative_bessel,
)
from theory.registry_plasticity import (
    RegistryTransportConfig,
    registry_grid,
    registry_step,
    solve_registry,
)


def test_active_bessel_identity_matches_independent_direct_sum() -> None:
    for q in (2.0, 6.0, 12.19):
        for delta, eta in ((0.0, 0.7), (0.23, 0.9), (0.5, 1.3)):
            bessel = shifted_inverse_power_bessel(q, delta, eta, modes=28)
            # q=2 has a slowly decaying O(1/N) direct-sum tail, so compare it
            # with the independent exact hyperbolic closed form.  Higher q
            # converges rapidly enough for a direct real-space check.
            direct = (
                inverse_power_q2_closed(delta, eta)
                if q == 2.0
                else inverse_power_direct(q, delta, eta, half_width=20_000)
            )
            tolerance = 3.0e-13 if q == 2.0 else 2.0e-12
            assert bessel == pytest.approx(direct, abs=tolerance)


def test_active_bessel_energy_is_periodic_and_its_force_is_exact_derivative() -> None:
    lattice = RegistryLattice(normal_ratio=0.93, bessel_modes=24)
    delta = 0.19
    assert registry_energy(delta + 1.0, lattice) == pytest.approx(
        registry_energy(delta, lattice), abs=2.0e-14
    )
    direct_derivative = (
        lattice.sigma_ratio**lattice.m
        * h_q_delta_derivative_bessel(
            lattice.m, delta, lattice.normal_ratio, 28, 64
        )
        - lattice.sigma_ratio**lattice.n
        * h_q_delta_derivative_bessel(
            lattice.n, delta, lattice.normal_ratio, 28, 64
        )
    )
    assert registry_energy_derivative(delta, lattice) == pytest.approx(
        direct_derivative, rel=2.0e-13
    )


def test_preferred_registry_is_computed_not_assumed() -> None:
    assert preferred_registry(RegistryLattice(normal_ratio=1.0)) == pytest.approx(0.5)
    assert preferred_registry(RegistryLattice(normal_ratio=1.5)) == pytest.approx(0.0)


def test_signed_schmid_projection_for_fcc_example() -> None:
    factor = schmid_factor([0, 0, 1], [1, 1, 1], [1, 0, -1])
    assert factor == pytest.approx(-1.0 / math.sqrt(6.0))
    with pytest.raises(ValueError):
        schmid_factor([0, 0, 1], [1, 1, 1], [1, 1, 0])


def _transport_config() -> RegistryTransportConfig:
    return RegistryTransportConfig(
        lattice=RegistryLattice(normal_ratio=1.0, bessel_modes=16),
        inverse_temperature=20.0,
        u_min=-6.0,
        u_max=7.0,
        cells=390,
    )


def test_zero_load_preserves_symmetry_and_probability() -> None:
    time = np.linspace(0.0, 4.0, 41)
    history = solve_registry(time, np.zeros_like(time), _transport_config(), max_dt=0.025)
    dx = history.registry[1] - history.registry[0]
    mass = np.sum(history.density, axis=1) * dx
    assert np.max(np.abs(mass - 1.0)) < 2.0e-12
    assert abs(history.mean_well_index[-1]) < 2.0e-10
    assert np.all(history.density >= 0.0)
    assert np.all(history.entropy_production >= 0.0)


def test_discrete_gibbs_state_has_zero_registry_current() -> None:
    config = RegistryTransportConfig(
        lattice=RegistryLattice(normal_ratio=1.0, bessel_modes=16),
        inverse_temperature=20.0,
        u_min=-2.0,
        u_max=3.0,
        cells=200,
    )
    registry, dx = registry_grid(config)
    energy = np.asarray(registry_energy(registry, config.lattice))
    density = np.exp(-config.inverse_temperature * (energy - np.min(energy)))
    density /= np.sum(density) * dx
    updated, currents = registry_step(density, energy, dx, 0.1, 0.0, config)
    assert np.max(np.abs(updated - density)) < 8.0e-15
    assert np.max(np.abs(currents)) < 8.0e-15


def test_resolved_shear_pulse_leaves_an_unwrapped_residual_slip_population() -> None:
    time = np.linspace(0.0, 30.0, 301)
    force = np.zeros_like(time)
    force[(time >= 2.0) & (time < 8.0)] = 0.55
    ramp = (time >= 8.0) & (time < 12.0)
    force[ramp] = 0.55 * (12.0 - time[ramp]) / 4.0
    history = solve_registry(time, force, _transport_config(), max_dt=0.025)

    # The external force has been zero for 18 reduced time units.  Intrawell
    # registry relaxes, while the signed population transferred into adjacent
    # wells leaves a nonzero residual mean z.  This, not a well crossing by
    # itself, is the operational reduced-plasticity statement.
    assert force[-1] == 0.0
    assert history.mean_well_index[-1] > 0.60
    assert abs(history.mean_intrawell_registry[-1]) < 2.0e-6
    assert history.boundary_probability.max() < 5.0e-8
    assert history.work[-1] > 0.0


def test_symmetric_zero_mean_cycles_do_not_create_directed_plastic_drift() -> None:
    time = np.linspace(0.0, 20.0, 201)
    force = 0.55 * np.sin(2.0 * math.pi * time / 4.0)
    history = solve_registry(time, force, _transport_config(), max_dt=0.025)
    assert abs(history.mean_well_index[-1]) < 2.0e-3
    assert history.boundary_probability.max() < 5.0e-8


def test_registry_grid_and_timestep_refinement() -> None:
    time = np.linspace(0.0, 20.0, 201)
    force = np.zeros_like(time)
    ramp_up = (time >= 2.0) & (time < 4.0)
    force[ramp_up] = 0.55 * (time[ramp_up] - 2.0) / 2.0
    force[(time >= 4.0) & (time < 8.0)] = 0.55
    ramp_down = (time >= 8.0) & (time < 12.0)
    force[ramp_down] = 0.55 * (12.0 - time[ramp_down]) / 4.0
    lattice = RegistryLattice(normal_ratio=1.0, bessel_modes=16)
    coarse = solve_registry(
        time,
        force,
        RegistryTransportConfig(
            lattice=lattice,
            inverse_temperature=20.0,
            u_min=-6.0,
            u_max=7.0,
            cells=260,
        ),
        max_dt=0.05,
    )
    fine = solve_registry(
        time,
        force,
        RegistryTransportConfig(
            lattice=lattice,
            inverse_temperature=20.0,
            u_min=-6.0,
            u_max=7.0,
            cells=520,
        ),
        max_dt=0.025,
    )
    assert abs(coarse.mean_well_index[-1] - fine.mean_well_index[-1]) < 0.007
    assert abs(coarse.work[-1] - fine.work[-1]) < 0.004
    assert max(coarse.boundary_probability.max(), fine.boundary_probability.max()) < 1.0e-8
