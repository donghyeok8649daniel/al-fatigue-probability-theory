import numpy as np

from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_4d import (
    CyclicLoad4D,
    Grid4DParams,
    PDE4DTimeParams,
    build_grid_4d,
    initial_gibbs_density_4d,
    observables_4d,
    run_probability_pde_4d,
)


def _params() -> ModelParams:
    return ModelParams(
        n_cells=2,
        kT=0.025,
        mobility_a=1.0,
        mobility_s=0.04,
        chi_axial_projection=0.20,
    )


def _grid() -> Grid4DParams:
    return Grid4DParams(n_a=7, n_s=9, s_wells=3, a_upper=1.50)


def _time() -> PDE4DTimeParams:
    return PDE4DTimeParams(max_dt=5.0e-4, cfl=0.30, record_interval=0.004)


def test_n2_gibbs_density_normalizes_and_is_nonnegative():
    model = TwoRowLJ(_params())
    model._build_opening_table()
    grid = build_grid_4d(model, _grid())
    density = initial_gibbs_density_4d(model, grid)
    mass = float(np.sum(density) * grid.cell_volume)
    assert abs(mass - 1.0) < 1.0e-12
    assert np.all(density >= 0.0)


def test_n2_interaction_produces_non_product_joint_density():
    model = TwoRowLJ(_params())
    model._build_opening_table()
    grid = build_grid_4d(model, _grid())
    density = initial_gibbs_density_4d(model, grid)
    obs = observables_4d(density, model, grid, 0.0)
    # The full interacting two-cell energy should not collapse to a product of
    # independent one-cell joint densities.  Coarse grids can weaken the
    # measured discrepancy, so only require a clear nonzero signal.
    assert obs["product_closure_l1_error"] > 1.0e-8


def test_n2_zero_load_mass_is_conserved():
    result = run_probability_pde_4d(
        model_params=_params(),
        grid_params=_grid(),
        time_params=_time(),
        load=CyclicLoad4D(force_min=0.0, force_max=0.0, period=1.0, cycles=0.01),
    )
    assert abs(float(result["survival"][-1]) - 1.0) < 5.0e-8
    assert float(result["initiation_probability"][-1]) < 5.0e-8


def test_n2_compression_does_not_create_tensile_opening_first_passage():
    result = run_probability_pde_4d(
        model_params=_params(),
        grid_params=_grid(),
        time_params=_time(),
        load=CyclicLoad4D(force_min=-0.8, force_max=-0.8, period=1.0, cycles=0.006),
    )
    assert float(result["initiation_probability"][-1]) < 5.0e-8


def test_n2_survival_is_monotone_nonincreasing():
    result = run_probability_pde_4d(
        model_params=_params(),
        grid_params=_grid(),
        time_params=_time(),
        load=CyclicLoad4D(force_min=0.0, force_max=4.2, period=0.08, cycles=0.1),
    )
    survival = np.asarray(result["survival"], dtype=float)
    assert np.all(np.diff(survival) <= 2.0e-9)
