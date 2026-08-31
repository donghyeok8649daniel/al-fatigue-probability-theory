# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: test_reflecting_probability_and_equilibrium_preservation
#   test_absorbing_mass_equals_integrated_outflux
#   test_tangent_instability_is_exact_operational_initiation_boundary
#   test_finite_rate_distribution_lag_and_quasistatic_loop_collapse
#   test_grid_timestep_refinement_and_free_energy_balance
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import numpy as np
import pytest

from theory.smoluchowski_escape import TransportConfig, domain_max, solve
from theory.normal_lj_chain import critical_stretch, normalized_lj_energy


def test_reflecting_probability_and_equilibrium_preservation():
    t = np.linspace(0, 2, 81)
    h = solve(t, np.zeros_like(t), TransportConfig(cells=120), max_dt=0.04)
    assert np.max(np.abs(h.survival-1)) < 2e-13
    assert np.max(np.abs(h.mean-h.mean[0])) < 2e-12
    assert np.all(h.entropy_production >= 0)


def test_absorbing_mass_equals_integrated_outflux():
    t = np.linspace(0, 3, 151)
    h = solve(t, np.full_like(t, 0.018),
              TransportConfig(cells=140, boundary="absorbing", lambda_max=1.20),
              max_dt=0.02)
    cumulative = np.cumsum(h.outflux[1:]*np.diff(t))
    assert np.max(np.abs((1-h.survival[1:])-cumulative)) < 2e-11
    assert np.all(np.diff(h.survival) <= 2e-13)
    assert np.all(h.hazard >= 0)


def test_tangent_instability_is_exact_operational_initiation_boundary():
    config = TransportConfig(cells=140, boundary="absorbing",
                             initiation_definition="tangent_instability",
                             lambda_max=9.0)
    assert domain_max(config) == pytest.approx(critical_stretch())
    t = np.linspace(0, 1, 51)
    h = solve(t, np.full_like(t, 0.01), config, max_dt=0.02)
    assert h.stretch[-1] < critical_stretch()
    assert np.all(h.tail_conditional == 0.0)


def test_finite_rate_distribution_lag_and_quasistatic_loop_collapse():
    def run(period):
        t = np.linspace(0, 4*period, 401)
        f = 0.006 + 0.004*np.sin(2*np.pi*t/period)
        return solve(t, f, TransportConfig(cells=100, inverse_temperature=100),
                     max_dt=min(0.03, period/80))
    # Dissipation first rises away from both limiting rates; compare a
    # finite-rate loop to the genuinely slow side of that response curve.
    fast = run(12.0)
    slow = run(200.0)
    fast_area = abs(fast.work[-1]-fast.work[300])
    slow_area = abs(slow.work[-1]-slow.work[300])
    assert fast_area > slow_area
    # Same force at quarter/three-quarter phases, different transport history.
    assert np.linalg.norm(fast.density[325]-fast.density[375]) > 1e-3


def test_grid_timestep_refinement_and_free_energy_balance():
    t = np.linspace(0, 8, 401)
    f = 0.01*(1-np.exp(-t/0.3))
    coarse = solve(t, f, TransportConfig(cells=80), max_dt=0.04)
    fine = solve(t, f, TransportConfig(cells=160), max_dt=0.01)
    assert abs(coarse.mean[-1]-fine.mean[-1]) < 2e-4
    assert abs(coarse.variance[-1]-fine.variance[-1]) < 3e-5

    c = TransportConfig(cells=160)
    dx = fine.stretch[1]-fine.stretch[0]
    phi = normalized_lj_energy(fine.stretch, c.m, c.n)
    free = np.asarray([np.sum(row*(phi-f[k]*(fine.stretch-1)
        + np.log(row)/c.inverse_temperature))*dx
        for k, row in enumerate(fine.density)])
    load_power = np.sum(-0.5*((fine.mean[1:]-1)+(fine.mean[:-1]-1))*np.diff(f))
    dissipation = np.trapezoid(fine.entropy_production, t)
    residual = (free[-1]-free[0])-load_power+dissipation
    assert abs(residual) < 5e-6
