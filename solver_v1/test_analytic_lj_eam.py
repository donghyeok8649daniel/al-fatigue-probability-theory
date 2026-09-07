import numpy as np
import pytest
from scipy.optimize import root

from solver_v1.analytic_lj_eam import AnalyticLJEAM, SquareRootEmbedding
from solver_v1.configurational_landscape import bound_configurational_barrier
from solver_v1.energy_surface import LocalEnergySurface
from solver_v1.exponential_density import ExponentialDensityParams
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_2d import (
    CyclicLoad2D,
    Grid2DParams,
    PDETimeParams,
    build_grid,
    energy_grid,
    run_probability_pde_2d,
)
from solver_v1.reduced_fast_a import bound_basin_conditional, conditional_fast_a


def _params() -> ModelParams:
    return ModelParams(
        n_cells=1,
        kT=0.02,
        mobility_a=1.0,
        mobility_s=0.05,
        chi_axial_projection=0.20,
    )


def _hybrid(*, amplitude: float = 0.08, rho0: float = 1.0) -> AnalyticLJEAM:
    return AnalyticLJEAM(
        _params(),
        density_params=ExponentialDensityParams(
            rho0=rho0, beta_rho=3.0, r_e=1.0, tol=1.0e-14
        ),
        embedding=SquareRootEmbedding(amplitude=amplitude, rho_ref=1.0),
    )


@pytest.fixture(scope="module")
def prepared_hybrid() -> AnalyticLJEAM:
    model = _hybrid()
    model._build_opening_table()
    return model


def test_square_root_embedding_derivatives_are_analytic() -> None:
    embedding = SquareRootEmbedding(amplitude=0.7, rho_ref=1.3)
    rho = np.array([0.3, 1.0, 4.0])
    h = 2.0e-5
    first_fd = (embedding.value(rho + h) - embedding.value(rho - h)) / (2.0 * h)
    second_fd = (
        embedding.value(rho + h)
        - 2.0 * embedding.value(rho)
        + embedding.value(rho - h)
    ) / h**2
    np.testing.assert_allclose(
        embedding.first_derivative(rho), first_fd, rtol=7.0e-10, atol=2.0e-10
    )
    np.testing.assert_allclose(
        embedding.second_derivative(rho), second_fd, rtol=7.0e-6, atol=2.0e-6
    )


def test_hybrid_local_gradient_and_hessian_match_finite_differences() -> None:
    model = _hybrid()
    h_gradient = 2.0e-6
    h_hessian = 1.0e-4
    for a, s in ((0.72, -0.18), (0.84, 0.11), (1.1, 0.31)):
        grad = np.array([model.local_deda(a, s), model.local_deds(a, s)])
        energy_gradient = np.array(
            [
                (model.local_energy(a + h_gradient, s) - model.local_energy(a - h_gradient, s))
                / (2.0 * h_gradient),
                (model.local_energy(a, s + h_gradient) - model.local_energy(a, s - h_gradient))
                / (2.0 * h_gradient),
            ]
        )
        np.testing.assert_allclose(grad, energy_gradient, rtol=2.0e-8, atol=2.0e-8)

        hessian = model.local_hessian(a, s)
        fd_a = (
            np.array(
                [model.local_deda(a + h_hessian, s), model.local_deds(a + h_hessian, s)]
            )
            - np.array(
                [model.local_deda(a - h_hessian, s), model.local_deds(a - h_hessian, s)]
            )
        ) / (2.0 * h_hessian)
        fd_s = (
            np.array(
                [model.local_deda(a, s + h_hessian), model.local_deds(a, s + h_hessian)]
            )
            - np.array(
                [model.local_deda(a, s - h_hessian), model.local_deds(a, s - h_hessian)]
            )
        ) / (2.0 * h_hessian)
        fd_hessian = np.column_stack((fd_a, fd_s))
        np.testing.assert_allclose(hessian, fd_hessian, rtol=8.0e-6, atol=8.0e-6)
        np.testing.assert_allclose(hessian, hessian.T, rtol=0.0, atol=2.0e-14)


def test_zero_embedding_amplitude_and_zero_density_recover_lj() -> None:
    lj = TwoRowLJ(_params())
    for hybrid in (_hybrid(amplitude=0.0), _hybrid(rho0=0.0)):
        assert hybrid.a0 == lj.a0
        assert hybrid.sigma_over_E_force_scale() == lj.sigma_over_E_force_scale()
        for a, s in ((0.7, -0.2), (0.9, 0.17), (1.3, 0.4)):
            assert hybrid.local_energy(a, s) == lj.local_energy(a, s)
            assert hybrid.local_deda(a, s) == lj.local_deda(a, s)
            assert hybrid.local_deds(a, s) == lj.local_deds(a, s)
            np.testing.assert_array_equal(
                hybrid.local_hessian(a, s), lj.local_hessian(a, s)
            )


def test_hybrid_equilibrium_hessian_and_periodicity() -> None:
    model = _hybrid()
    assert abs(model.local_deda(model.a0, 0.0)) < 1.0e-10
    hessian = model.local_hessian(model.a0, 0.0)
    np.testing.assert_allclose(hessian, hessian.T, rtol=0.0, atol=1.0e-14)
    assert np.all(np.linalg.eigvalsh(hessian) > 0.0)
    for a, s in ((0.75, -0.21), (1.1, 0.33)):
        np.testing.assert_allclose(
            model.local_energy(a, s + model.p.b),
            model.local_energy(a, s),
            rtol=0.0,
            atol=2.0e-13,
        )
        np.testing.assert_allclose(
            model.local_deds(a, s + model.p.b),
            model.local_deds(a, s),
            rtol=0.0,
            atol=2.0e-12,
        )
    reference = float(model.local_energy(model.a0, 0.0))
    for a, s in ((0.74, 0.11), (0.95, -0.22)):
        assert float(model.local_energy_referenced(a, s)) == pytest.approx(
            model.local_energy(a, s) - reference, abs=2.0e-15
        )


def test_hybrid_relaxed_kappa_gives_unit_small_stress_tangent() -> None:
    model = _hybrid()
    hessian = model.local_hessian(model.a0, 0.0)
    c = np.array([1.0, model.p.chi_axial_projection])
    for reduced_stress in (-1.0e-5, 1.0e-5):
        force = float(model.force_from_sigma_over_E(reduced_stress))
        initial = np.array([model.a0, 0.0]) + np.linalg.solve(hessian, force * c)
        stationary = root(
            lambda q: np.array(
                [
                    model.local_deda(float(q[0]), float(q[1])) - force,
                    model.local_deds(float(q[0]), float(q[1]))
                    - force * model.p.chi_axial_projection,
                ]
            ),
            initial,
            options={"xtol": 1.0e-11},
        )
        strain = (
            stationary.x[0] - model.a0
            + model.p.chi_axial_projection * stationary.x[1]
        ) / model.a0
        assert stationary.success
        assert np.linalg.norm(stationary.fun, ord=np.inf) < 1.0e-9
        assert abs(strain / reduced_stress - 1.0) < 2.0e-3


def test_hybrid_implements_energy_surface_and_probability_energy_grid() -> None:
    model = _hybrid()
    assert isinstance(model, LocalEnergySurface)
    np.testing.assert_array_equal(
        model.grad(0.8, 0.1),
        [model.local_deda(0.8, 0.1), model.local_deds(0.8, 0.1)],
    )
    np.testing.assert_array_equal(model.hessian(0.8, 0.1), model.local_hessian(0.8, 0.1))
    assert model.energy(0.8, 0.1) == model.local_energy(0.8, 0.1)
    grid = build_grid(
        model, Grid2DParams(n_a=7, n_s=9, s_wells=3, a_upper=1.4)
    )
    force = 0.4
    values = energy_grid(model, grid, force)
    for ia in (0, 3, 6):
        for is_ in (0, 4, 8):
            expected = model.local_energy(grid.a[ia], grid.s[is_]) - force * (
                (grid.a[ia] - model.a0)
                + model.p.chi_axial_projection * grid.s[is_]
            )
            np.testing.assert_allclose(values[ia, is_], expected, rtol=0.0, atol=2.0e-13)


def test_fixed_boundary_hybrid_free_energy_gradient_identity() -> None:
    model = _hybrid()
    s = 0.12
    force = 0.35
    kwargs = dict(a_lower=0.69, a_upper=1.35, quadrature_order=160)
    conditional = conditional_fast_a(model, s, force, **kwargs)
    h = 2.0e-5
    plus = conditional_fast_a(model, s + h, force, **kwargs)
    minus = conditional_fast_a(model, s - h, force, **kwargs)
    finite_difference = (plus.F_eff - minus.F_eff) / (2.0 * h)
    np.testing.assert_allclose(
        conditional.dF_eff_ds,
        conditional.mean_dG_ds,
        rtol=0.0,
        atol=2.0e-14,
    )
    np.testing.assert_allclose(
        conditional.dF_eff_ds, finite_difference, rtol=2.0e-7, atol=2.0e-8
    )


def test_hybrid_opening_table_and_full_pde_interface(
    prepared_hybrid: AnalyticLJEAM,
) -> None:
    lj = TwoRowLJ(_params())
    lj._build_opening_table()
    assert prepared_hybrid._fc_interp(0.0) > lj._fc_interp(0.0)
    assert np.min(prepared_hybrid._fc_grid) > np.min(lj._fc_grid)
    result = run_probability_pde_2d(
        prepared_model=prepared_hybrid,
        grid_params=Grid2DParams(
            n_a=7, n_s=9, s_wells=3, a_upper=1.4
        ),
        time_params=PDETimeParams(
            max_dt=1.0e-3,
            record_interval=1.0e-3,
            integrator="implicit",
        ),
        load=CyclicLoad2D(
            force_min=0.0, force_max=0.0, period=1.0, cycles=0.0
        ),
        preload_force=0.0,
    )
    assert result["model"] is prepared_hybrid
    assert abs(float(result["intact_probability_mass"][-1]) - 1.0) < 1.0e-12
    assert float(result["cumulative_absorbed_mass"][-1]) == 0.0


def test_hypothetical_embedding_sensitivity_raises_force_four_barriers(
    prepared_hybrid: AnalyticLJEAM,
) -> None:
    barrier = bound_configurational_barrier(
        prepared_hybrid, 4.0, n_s=101, quadrature_order=40
    )
    assert barrier is not None
    assert 5.0e-3 < barrier.barrier < 1.0e-2
    assert barrier.opening_barrier_at_saddle > 9.0e-2


def test_hybrid_moving_bound_basin_keeps_leibniz_correction(
    prepared_hybrid: AnalyticLJEAM,
) -> None:
    s = 0.17
    force = 3.5
    step = 1.0e-5
    center = bound_basin_conditional(
        prepared_hybrid, s, force, quadrature_order=192
    )
    plus = bound_basin_conditional(
        prepared_hybrid, s + step, force, quadrature_order=192
    )
    minus = bound_basin_conditional(
        prepared_hybrid, s - step, force, quadrature_order=192
    )
    finite_difference = (plus.F_eff - minus.F_eff) / (2.0 * step)
    assert abs(center.boundary_correction_s) > 1.0e-8
    assert finite_difference == pytest.approx(
        center.dF_eff_ds, rel=3.0e-5, abs=3.0e-8
    )
