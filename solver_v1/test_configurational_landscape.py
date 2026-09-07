import pytest

from solver_v1.configurational_landscape import bound_configurational_barrier
from solver_v1.model import ModelParams, TwoRowLJ


@pytest.fixture(scope="module")
def model() -> TwoRowLJ:
    value = TwoRowLJ(
        ModelParams(
            n_cells=1,
            kT=0.02,
            mobility_a=1.0,
            mobility_s=0.05,
            chi_axial_projection=0.20,
        )
    )
    value._build_opening_table()
    return value


def test_bound_basin_configurational_barrier_decreases_toward_spinodal(
    model: TwoRowLJ,
) -> None:
    barriers = [
        bound_configurational_barrier(
            model, force, n_s=101, quadrature_order=40
        )
        for force in (2.0, 3.5, 4.0)
    ]
    assert all(value is not None for value in barriers)
    values = [value.barrier for value in barriers if value is not None]
    assert values[0] > values[1] > values[2] > 0.0
    for value in barriers:
        assert value is not None
        assert value.minimum_curvature > 0.0
        assert value.saddle_curvature < 0.0
        assert value.opening_barrier_at_minimum > 0.0
        assert value.opening_barrier_at_saddle > 0.0


def test_bound_configurational_minimum_and_saddle_merge_near_force_4p03(
    model: TwoRowLJ,
) -> None:
    before = bound_configurational_barrier(
        model, 4.02, n_s=151, quadrature_order=48
    )
    after = bound_configurational_barrier(
        model, 4.03, n_s=151, quadrature_order=48
    )
    assert before is not None
    assert before.barrier < 5.0e-5
    assert after is None
