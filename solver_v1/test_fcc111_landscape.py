import numpy as np
import pytest

from solver_v1.energy_model_registry import (
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    TWO_ROW_LJ_REFERENCE,
    build_energy_model,
)
from solver_v1.fcc111_full_energy import FullFCC111StackEnergy
from solver_v1.fcc111_landscape import (
    configurational_spinodal,
    loaded_energy,
    opening_barrier,
    opening_spinodal,
    relaxed_registry_barrier,
)


@pytest.fixture(scope="module")
def full_lj():
    return FullFCC111StackEnergy.from_reduced_model(
        build_energy_model(TWO_ROW_LJ_REFERENCE)
    )


def test_loaded_energy_uses_existing_generalized_force_convention(full_lj):
    a, s, force = full_lj.a0 + 0.02, 0.11, 0.3
    expected = full_lj.energy(a, s) - force * (
        a - full_lj.a0 + full_lj.p.chi_axial_projection * s
    )
    assert loaded_energy(full_lj, a, s, force) == pytest.approx(expected)


def test_full_lj_static_barriers_and_spinodals_are_distinct(full_lj):
    opening_force, opening_a = opening_spinodal(full_lj, 0.0)
    registry = configurational_spinodal(full_lj)
    assert registry is not None
    assert 0.0 < registry.force < opening_force
    assert registry.residual_norm < 1e-8
    assert opening_a > full_lj.a0
    barrier = relaxed_registry_barrier(full_lj, 0.0, samples=25)
    opening = opening_barrier(full_lj, 0.0, 0.0)
    assert barrier is not None and barrier.barrier > 0.0
    assert opening is not None and opening.barrier > barrier.barrier


def test_hypothetical_hybrid_has_finite_reference_spinodals():
    model = FullFCC111StackEnergy.from_reduced_model(
        build_energy_model(ANALYTIC_LJ_EAM_HYPOTHETICAL)
    )
    registry = configurational_spinodal(model)
    opening_force, _ = opening_spinodal(model, 0.0)
    assert registry is not None
    assert np.isfinite(registry.force)
    assert 0.0 < registry.force < opening_force
