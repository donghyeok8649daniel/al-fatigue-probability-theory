import numpy as np
import pytest

from .static_bloch_targets import source_bloch_targets,append_bloch_problem
from .reference_eam_targets import MishinRigidFCCReference
from .quartic_angular_material import QuarticSymmetryTailBulkBasis


def test_predeclared_source_spatial_targets_have_independent_roles_and_units():
    source=MishinRigidFCCReference();targets=source_bloch_targets(source)
    assert len(targets)==14
    fit={t.fractional_cubic for t in targets if t.observation.role=='fit'}
    hold={t.fractional_cubic for t in targets if t.observation.role=='heldout'}
    assert not fit&hold
    assert all(t.observation.units=='eV/L0^2' for t in targets)
    assert all(t.observation.scale==.1*t.observation.target for t in targets)
    with pytest.raises(ValueError):source_bloch_targets(MishinRigidFCCReference(lattice_angstrom=4.065))


def test_spatial_projection_preserves_full_coefficients_and_tail_bound():
    targets=source_bloch_targets(MishinRigidFCCReference())
    first=QuarticSymmetryTailBulkBasis([2.8,5.6,11.9],radius=12.)
    second=QuarticSymmetryTailBulkBasis([2.8,5.6,11.9],radius=16.)
    a,obs,tail=append_bloch_problem(np.zeros((0,10)),[],first,targets)
    b,_,_=append_bloch_problem(np.zeros((0,10)),[],second,targets)
    assert np.all(abs(a-b)<=tail+1e-12)
    assert np.array_equal(a[:,-1],np.zeros(14)) # quartic moment has no harmonic bulk term
    assert len(obs)==14 and a.shape==(14,10)
