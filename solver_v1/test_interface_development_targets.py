"""The new development loss must not silently reuse held-out claims."""
import numpy as np
import pytest

from .interface_development_targets import development_observations
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import UNITS


@pytest.fixture(scope='module')
def targets():
    source, old, states = source_and_targets()
    rows, provenance = development_observations(source, old, states)
    return source, old, rows, provenance


def test_inspected_failures_reclassified_and_new_tests_excluded(targets):
    _, old, rows, provenance = targets
    by_name = {o.name: o for o in rows}
    for o in old:
        if o.role == 'heldout':
            assert by_name[o.name].role == 'fit'
            assert 'development' in provenance[o.name]
    test_rows = [o for o in rows if o.role == 'heldout']
    assert len(test_rows) == 30
    assert all(o.name.startswith('new_') for o in test_rows)
    assert len(by_name) == len(rows)


def test_source_jets_and_declared_units_are_consistent(targets):
    source, _, rows, _ = targets
    for obs in rows:
        if obs.state is None:
            continue
        v = source.evaluate(obs.state)
        jet = np.r_[v.energy, v.gradient, v.hessian[0], v.hessian[1, 1:], v.hessian[2, 2]]
        assert float(np.asarray(obs.jet_weights)@jet) == pytest.approx(obs.target, abs=1e-14)
        assert obs.scale > 0
    force = next(o for o in rows if o.name == 'opening_1.3h_normal_force')
    assert UNITS.force_to_traction_mpa(force.scale) == pytest.approx(250.)


def test_old_observations_not_mutated_and_no_specimen_or_kinetic_loss(targets):
    _, old, rows, _ = targets
    assert sum(o.role == 'heldout' for o in old) == 9
    assert len([o for o in rows if o.bulk_index is not None]) == 5
    assert all(not any(key in o.name for key in ('yield', 'fatigue', 'mobility', 'A_c')) for o in rows)
