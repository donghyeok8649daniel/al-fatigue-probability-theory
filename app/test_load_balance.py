import numpy as np
import pytest
from .specimen_mesh import cylinder
from .load_balance import face_geometry, stress_traction, resultant, correction_operator, balanced_traction


def test_uniform_stress_closed_body_balances():
    mesh = cylinder()
    stress = np.array([[10, 2, 3], [2, 20, 4], [3, 4, 30]])
    traction = stress_traction(mesh, np.arange(len(mesh.faces)), stress)
    np.testing.assert_allclose(resultant(mesh, traction), 0, atol=2e-10)


def test_top_and_bottom_tension_have_opposite_force():
    mesh = cylinder()
    _, normals = face_geometry(mesh)
    top = np.flatnonzero(normals[:, 2] > .99)
    bottom = np.flatnonzero(normals[:, 2] < -.99)
    stress = np.diag([0, 0, 10])
    ft = resultant(mesh, stress_traction(mesh, top, stress))
    fb = resultant(mesh, stress_traction(mesh, bottom, stress))
    assert ft[2] > 0 and fb[2] < 0
    np.testing.assert_allclose(ft+fb, 0, atol=1e-10)


def test_correction_balances_arbitrary_time_history():
    mesh = cylinder()
    _, normals = face_geometry(mesh)
    top = np.flatnonzero(normals[:, 2] > .99)
    bottom = np.flatnonzero(normals[:, 2] < -.99)
    correction = correction_operator(mesh, bottom)
    for t in [0, .013, .4, 1.7]:
        stress = np.array([[0, 0, 3*np.cos(t)], [0, 0, 0], [3*np.cos(t), 0, 7+np.sin(3*t)]])
        traction = stress_traction(mesh, top, stress)
        before = traction.copy()
        np.testing.assert_allclose(resultant(mesh, balanced_traction(mesh, traction, correction)), 0, atol=1e-9)
        np.testing.assert_array_equal(traction, before)


def test_rank_deficient_correction_and_asymmetric_stress_refused():
    mesh = cylinder()
    with pytest.raises(ValueError): correction_operator(mesh, [0])
    with pytest.raises(ValueError): stress_traction(mesh, [0], [[1,2,0],[0,1,0],[0,0,1]])
