"""Geometry/normalization tests; optional engine controls run in the audit CLI."""
import numpy as np
import pytest

from .silicon_atomistic_reference import restricted_cell, force_metrics, periodic_basis_audit


def test_restricted_cell_preserves_energy_gradient_and_volume():
    rng = np.random.default_rng(38501)
    cell = np.array([[4., 1., .5], [.2, 5., 1.], [.1, -.5, 6.]])
    b, q = restricted_cell(cell)
    np.testing.assert_allclose(cell @ q, b, atol=2e-15)
    np.testing.assert_allclose(q.T @ q, np.eye(3), atol=5e-16)
    assert np.linalg.det(q) == pytest.approx(1.)
    assert np.linalg.det(b) == pytest.approx(np.linalg.det(cell))
    assert np.all(np.diag(b) > 0)
    r, f, dr = rng.normal(size=(3, 13, 3))
    np.testing.assert_allclose(np.sum(f*dr), np.sum((f @ q)*(dr @ q)), atol=1e-14)
    np.testing.assert_allclose((f @ q) @ q.T, f, atol=1e-15)
    np.testing.assert_allclose(np.linalg.norm(r @ q, axis=1), np.linalg.norm(r, axis=1))


@pytest.mark.parametrize('cell', [np.zeros((3, 3)), np.diag([-1., 1., 1.]),
                                      np.full((3, 3), np.nan), np.eye(2)])
def test_restricted_cell_rejects_invalid_geometry(cell):
    with pytest.raises(ValueError):
        restricted_cell(cell)


def test_force_metric_component_weighting():
    first = force_metrics(np.array([[3., 0., 0.]]), np.zeros((1, 3)))
    second = force_metrics(np.zeros((2, 3)), np.zeros((2, 3)))
    assert first['component_rmse_eV_A'] == pytest.approx(np.sqrt(3.))
    assert first['max_vector_error_eV_A'] == 3.
    combined = np.sqrt((first['squared_error']+second['squared_error']) /
                       (first['components']+second['components']))
    assert combined == 1.


def test_force_metric_rejects_nonfinite_data():
    with pytest.raises(ValueError):
        force_metrics(np.full((2, 3), np.nan), np.zeros((2, 3)))


def test_left_handed_periodic_reindex_preserves_lattice_vectors():
    c = np.array([[4., .1, 0.], [.2, 5., 0.], [.3, .4, -8.]])
    changed, signs = periodic_basis_audit(c, (True, True, True))
    shifts = np.array([[1, 3, -2], [-2, 0, 1], [0, 0, 0]])
    np.testing.assert_allclose(shifts @ c, (shifts*signs) @ changed)
    assert signs.tolist() == [1, 1, -1]
    with pytest.raises(ValueError):
        periodic_basis_audit(c, (True, True, False))
