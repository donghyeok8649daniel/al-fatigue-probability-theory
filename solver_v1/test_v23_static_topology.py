from types import SimpleNamespace

import numpy as np

from .fcc111_geometry import fcc111_geometry
from .run_v23_static_topology import registry_basis_reduced
from .vector_interface_reference import MishinVectorInterfaceReference
from .vector_registry_audit import stationary_state


def test_source_and_candidate_registry_equivalence_use_same_reduced_units():
    physical = fcc111_geometry(4.05)
    length = physical.b
    source = MishinVectorInterfaceReference(SimpleNamespace(
        interpolation='cubic', geometry=physical, h=physical.h111), length)
    candidate = SimpleNamespace(geometry=fcc111_geometry(np.sqrt(2)))
    actual = registry_basis_reduced(source)
    np.testing.assert_allclose(actual, registry_basis_reduced(candidate), atol=2e-16)
    np.testing.assert_allclose(np.linalg.solve(actual, actual@np.array([1., -2.])), [1., -2.])
    partial = np.linalg.solve(actual, np.array([.5, np.sqrt(3)/6]))
    assert np.linalg.norm(partial-np.rint(partial)) > .4


def test_positive_Hxx_does_not_exclude_coupled_index_one_saddle():
    H = np.array([[8., 0., 0.], [0., 1., 2.], [0., 2., 1.]])
    q0 = np.array([1., 0., 0.])
    def evaluate(q):
        delta = np.asarray(q)-q0
        return SimpleNamespace(energy=float(.5*delta@H@delta), gradient=H@delta, hessian=H)
    result = stationary_state(SimpleNamespace(h=1., evaluate=evaluate), q0, expected_index=1)
    assert result['valid'] and result['morse_index'] == 1
    assert result['evaluation'].hessian[1, 1] > 0
    np.testing.assert_allclose(result['eigenvalues'], [-1., 3., 8.], atol=2e-15)
