import numpy as np
from scipy.optimize import minimize_scalar

from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_2d import (
    CyclicLoad2D,
    Grid2DParams,
    PDETimeParams,
    build_grid,
    cyclic_load_from_sigma_over_E,
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


def test_sigma_over_E_force_scale_matches_fast_normal_young_limit():
    model = TwoRowLJ(_params())
    reduced_stress = 1.0e-5
    force = float(model.force_from_sigma_over_E(reduced_stress))
    result = minimize_scalar(
        lambda a: model.local_energy(float(a), 0.0) - force * (float(a) - model.a0),
        bounds=(0.97 * model.a0, 1.03 * model.a0),
        method="bounded",
        options={"xatol": 1.0e-13},
    )
    strain = (float(result.x) - model.a0) / model.a0
    assert result.success
    assert abs(strain / reduced_stress - 1.0) < 2.0e-3
    assert model.sigma_over_E_force_scale() > 1.0


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
    ia = int(np.argmin(np.abs(grid.a - model.a0)) )
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
