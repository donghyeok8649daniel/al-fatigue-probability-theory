import numpy as np
from scipy.optimize import root, root_scalar

from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_2d import (
    CyclicLoad2D,
    Grid2DParams,
    PDETimeParams,
    build_grid,
    configurational_interface_fluxes,
    configurational_well_observables,
    cyclic_load_from_sigma_over_E,
    energy_grid,
    initial_gibbs_density,
    observables,
    run_probability_pde_2d,
)


def _params() -> ModelParams:
    return ModelParams(
        n_cells=1,
        kT=0.02,
        mobility_a=1.0,
        mobility_s=0.05,
        chi_axial_projection=0.20,
    )


def test_gibbs_initial_density_normalizes_without_sampling():
    model = TwoRowLJ(_params())
    model._build_opening_table()
    grid = build_grid(
        model,
        Grid2DParams(n_a=25, n_s=41, s_wells=3, a_upper=1.6),
    )
    density = initial_gibbs_density(model, grid)
    mass = float(np.sum(density) * grid.cell_volume)
    assert abs(mass - 1.0) < 1.0e-12
    assert np.all(density >= 0.0)


def test_reference_spacing_satisfies_pristine_stationarity_to_high_precision():
    model = TwoRowLJ(_params())
    assert abs(model.local_deda(model.a0, 0.0)) < 1.0e-10


def test_pristine_hessian_is_symmetric_and_stable():
    model = TwoRowLJ(_params())
    hessian = model.local_hessian(model.a0, 0.0)
    np.testing.assert_allclose(hessian, hessian.T, rtol=0.0, atol=1.0e-14)
    assert abs(float(hessian[0, 1])) < 1.0e-12
    assert np.all(np.linalg.eigvalsh(hessian) > 0.0)


def test_frozen_normal_force_scale_matches_stationarity_equation():
    model = TwoRowLJ(_params())
    for reduced_stress in (-1.0e-5, 1.0e-5):
        force = model.frozen_normal_sigma_over_E_force_scale() * reduced_stress
        result = root_scalar(
            lambda a: model.local_deda(float(a), 0.0) - force,
            bracket=(0.97 * model.a0, 1.03 * model.a0),
            method="brentq",
            xtol=5.0e-15,
            rtol=4.0 * np.finfo(float).eps,
        )
        strain = (float(result.root) - model.a0) / model.a0
        assert result.converged
        assert abs(strain / reduced_stress - 1.0) < 2.0e-3


def test_sigma_over_E_force_scale_matches_relaxed_total_axial_strain():
    model = TwoRowLJ(_params())
    hessian = model.local_hessian(model.a0, 0.0)
    c = np.array([1.0, model.p.chi_axial_projection], dtype=float)
    for reduced_stress in (-1.0e-5, 1.0e-5):
        force = float(model.force_from_sigma_over_E(reduced_stress))
        initial = np.array([model.a0, 0.0]) + np.linalg.solve(hessian, force * c)
        result = root(
            lambda q: np.array(
                [
                    model.local_deda(float(q[0]), float(q[1])) - force,
                    model.local_deds(float(q[0]), float(q[1]))
                    - force * model.p.chi_axial_projection,
                ]
            ),
            initial,
            method="hybr",
            options={"xtol": 1.0e-11},
        )
        total_strain = (
            (float(result.x[0]) - model.a0)
            + model.p.chi_axial_projection * float(result.x[1])
        ) / model.a0
        assert result.success
        assert np.linalg.norm(result.fun, ord=np.inf) < 1.0e-10
        assert abs(total_strain / reduced_stress - 1.0) < 2.0e-3
    assert model.sigma_over_E_force_scale() > 0.0
    assert (
        model.sigma_over_E_force_scale()
        < model.frozen_normal_sigma_over_E_force_scale()
    )


def test_reduced_stress_history_is_not_used_as_raw_lj_force():
    model = TwoRowLJ(_params())
    load = cyclic_load_from_sigma_over_E(
        model,
        sigma_over_E_min=-0.01,
        sigma_over_E_max=0.02,
        period=1.0,
        cycles=1.0,
    )
    scale = model.sigma_over_E_force_scale()
    assert np.isclose(load.force_min, -0.01 * scale)
    assert np.isclose(load.force_max, 0.02 * scale)


def test_strain_decomposition_contains_signed_well_index_plasticity():
    model = TwoRowLJ(_params())
    model._build_opening_table()
    grid = build_grid(
        model,
        Grid2DParams(n_a=25, n_s=61, s_wells=3, a_upper=1.6),
    )
    density = np.zeros((grid.a.size, grid.s.size), dtype=float)
    ia = int(np.argmin(np.abs(grid.a - model.a0)))
    candidates = np.where(grid.s > 0.55 * model.p.b)[0]
    assert candidates.size
    is_ = int(candidates[0])
    density[ia, is_] = 1.0 / grid.cell_volume

    obs = observables(density, model, grid, force=0.0)
    assert obs["plastic_strain"] > 0.0
    assert obs["mean_well_index"] > 0.5
    assert obs["plastic_well_activity"] > 0.5
    assert np.isclose(
        obs["strain"],
        obs["normal_strain"] + obs["intrawell_strain"] + obs["plastic_strain"],
        rtol=0.0,
        atol=1.0e-14,
    )


def test_principal_well_intrawell_strain_does_not_create_plastic_strain():
    model = TwoRowLJ(_params())
    model._build_opening_table()
    grid = build_grid(
        model,
        Grid2DParams(n_a=25, n_s=61, s_wells=3, a_upper=1.6),
    )
    density = np.zeros((grid.a.size, grid.s.size), dtype=float)
    ia = int(np.argmin(np.abs(grid.a - model.a0)))
    candidates = np.where(
        (grid.s > 0.05 * model.p.b) & (grid.s < 0.45 * model.p.b)
    )[0]
    assert candidates.size
    is_ = int(candidates[0])
    density[ia, is_] = 1.0 / grid.cell_volume

    obs = observables(density, model, grid, force=0.0)
    assert obs["intrawell_strain"] > 0.0
    assert abs(obs["plastic_strain"]) < 1.0e-14
    assert np.isclose(
        obs["strain"],
        obs["normal_strain"] + obs["intrawell_strain"] + obs["plastic_strain"],
        rtol=0.0,
        atol=1.0e-14,
    )


def test_zero_load_probability_mass_is_conserved():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=21, n_s=31, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(max_dt=1.0e-3, cfl=0.40, record_interval=0.01),
        load=CyclicLoad2D(force_min=0.0, force_max=0.0, period=1.0, cycles=0.03),
    )
    assert abs(float(result["survival"][-1]) - 1.0) < 1.0e-9
    assert float(result["initiation_probability"][-1]) < 1.0e-9
    assert np.allclose(
        result["strain"],
        result["normal_strain"] + result["intrawell_strain"] + result["plastic_strain"],
        rtol=0.0,
        atol=2.0e-12,
        equal_nan=True,
    )


def test_compression_does_not_create_tensile_opening_first_passage():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=21, n_s=31, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(max_dt=1.0e-3, cfl=0.40, record_interval=0.01),
        load=CyclicLoad2D(force_min=-1.0, force_max=-1.0, period=1.0, cycles=0.02),
    )
    assert float(result["initiation_probability"][-1]) < 1.0e-9


def test_above_opening_spinodal_absorbs_probability_mass():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=21, n_s=31, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(max_dt=1.0e-3, cfl=0.40, record_interval=0.01),
        load=CyclicLoad2D(force_min=6.0, force_max=6.0, period=1.0, cycles=0.0),
    )
    assert float(result["survival"][-1]) < 1.0e-12
    assert float(result["initiation_probability"][-1]) > 1.0 - 1.0e-12


def test_survival_is_monotone_nonincreasing():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=21, n_s=31, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(max_dt=1.0e-3, cfl=0.40, record_interval=0.01),
        load=CyclicLoad2D(force_min=0.0, force_max=4.2, period=0.2, cycles=0.5),
    )
    survival = np.asarray(result["survival"], dtype=float)
    assert np.all(np.diff(survival) <= 1.0e-10)


def test_physical_first_passage_bookkeeping_matches_discrete_absorption_flux():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=11, n_s=18, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(
            max_dt=2.0e-3,
            cfl=0.40,
            record_interval=4.0e-3,
            integrator="implicit",
        ),
        load=CyclicLoad2D(force_min=4.2, force_max=4.2, period=1.0, cycles=0.008),
    )
    absorbed = np.asarray(result["cumulative_absorbed_mass"], dtype=float)
    survival = np.asarray(result["survival"], dtype=float)
    time = np.asarray(result["time"], dtype=float)
    flux = np.asarray(result["first_passage_flux"], dtype=float)
    integrated_from_records = float(np.sum(flux[1:] * np.diff(time)))

    assert absorbed[-1] > 1.0e-8
    assert np.all(np.diff(absorbed) >= -1.0e-15)
    np.testing.assert_allclose(survival, 1.0 - absorbed, rtol=0.0, atol=0.0)
    assert np.isclose(
        integrated_from_records,
        float(result["integrated_first_passage_flux"][-1]),
        rtol=2.0e-14,
        atol=1.0e-18,
    )
    assert np.isclose(
        absorbed[-1],
        float(result["initial_absorbed_mass"][-1]) + integrated_from_records,
        rtol=2.0e-14,
        atol=1.0e-18,
    )
    assert np.max(np.abs(result["flux_consistency_residual"])) < 1.0e-14
    assert np.max(np.abs(result["mass_balance_residual"])) < 1.0e-12


def test_aligned_well_populations_obey_interwell_flux_balance():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=11, n_s=18, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(
            max_dt=2.0e-3,
            cfl=0.40,
            record_interval=4.0e-3,
            integrator="implicit",
        ),
        load=CyclicLoad2D(force_min=0.5, force_max=0.5, period=1.0, cycles=0.008),
    )
    populations = np.asarray(result["well_populations"], dtype=float)
    net = np.asarray(result["interwell_net_flux"], dtype=float)
    gross = np.asarray(result["interwell_gross_flux"], dtype=float)

    assert populations.shape[1] == 3
    assert net.shape[1] == 2
    assert np.all(gross + 1.0e-30 >= np.abs(net))
    assert np.max(np.abs(result["interwell_boundary_alignment_error"])) < 1.0e-14
    assert np.max(np.abs(result["well_population_balance_residual"])) < 1.0e-12
    np.testing.assert_allclose(
        np.sum(populations, axis=1),
        result["intact_probability_mass"],
        rtol=0.0,
        atol=2.0e-14,
    )
    assert np.max(np.abs(result["registry_moment_balance_residual"])) < 1.0e-12
    factor = (
        result["model"].p.chi_axial_projection
        * result["model"].p.b
        / result["model"].a0
    )
    np.testing.assert_allclose(
        result["plastic_strain"],
        factor
        * result["unnormalized_registry_moment"]
        / result["intact_probability_mass"],
        rtol=0.0,
        atol=2.0e-14,
    )


def test_initial_crack_absorption_does_not_create_interwell_plastic_flow():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=11, n_s=18, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(
            max_dt=2.0e-3,
            cfl=0.40,
            record_interval=0.01,
            integrator="implicit",
        ),
        load=CyclicLoad2D(force_min=6.0, force_max=6.0, period=1.0, cycles=0.0),
        preload_force=0.0,
    )
    assert result["cumulative_absorbed_mass"][-1] > 1.0 - 1.0e-12
    assert result["accumulated_net_registry_transfer"][-1] == 0.0
    assert result["net_interwell_registry_rate"][-1] == 0.0


def test_configurational_interface_flux_is_zero_for_discrete_gibbs_state():
    model = TwoRowLJ(_params())
    model._build_opening_table()
    grid = build_grid(
        model,
        Grid2DParams(n_a=17, n_s=30, s_wells=3, a_upper=1.55),
    )
    density = initial_gibbs_density(
        model, grid, preload_force=0.0, principal_well_only=False
    )
    energy = energy_grid(model, grid, force=0.0)
    signed, gross = configurational_interface_fluxes(density, energy, model, grid)
    well = configurational_well_observables(density, energy, model, grid)

    assert np.max(np.abs(signed)) < 1.0e-12
    assert np.all(gross >= 0.0)
    assert np.max(np.abs(well["interwell_net_flux"])) < 1.0e-12
