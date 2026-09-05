import numpy as np

from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_2d import (
    CyclicLoad2D,
    Grid2DParams,
    PDETimeParams,
    build_grid,
    initial_gibbs_density,
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


def test_zero_load_probability_mass_is_conserved():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=21, n_s=31, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(max_dt=1.0e-3, cfl=0.40, record_interval=0.01),
        load=CyclicLoad2D(force_min=0.0, force_max=0.0, period=1.0, cycles=0.03),
    )
    assert abs(float(result["survival"][-1]) - 1.0) < 1.0e-9
    assert float(result["initiation_probability"][-1]) < 1.0e-9


def test_compression_does_not_create_tensile_opening_first_passage():
    result = run_probability_pde_2d(
        model_params=_params(),
        grid_params=Grid2DParams(n_a=21, n_s=31, s_wells=3, a_upper=1.55),
        time_params=PDETimeParams(max_dt=1.0e-3, cfl=0.40, record_interval=0.01),
        load=CyclicLoad2D(force_min=-1.0, force_max=-1.0, period=1.0, cycles=0.02),
    )
    assert float(result["initiation_probability"][-1]) < 1.0e-9


def test_above_opening_spinodal_absorbs_probability_mass():
    # Static force above the local LJ opening strength should remove the
    # initially intact probability mass through the opening criterion without
    # requiring any random trajectory sampling.
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
