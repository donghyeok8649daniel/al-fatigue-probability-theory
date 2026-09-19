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


def test_unassigned_faces_and_empty_loads_have_exactly_zero_external_force():
    from .load_balance import applied_traction, balance_report
    from .load_workflow import FaceLoad
    mesh = cylinder()
    ids = np.flatnonzero(face_geometry(mesh)[1][:,2] > .99)
    load = FaceLoad('top',1200.,1500.,0.,0.,tuple(ids),mesh.areas[ids].sum(),
                    '0,0,0;0,0,0;0,0,normal_mean+normal_amp*sin(2*pi*f*t)',25.)
    for t in (0.,.01,.03):
        empty = applied_traction(mesh, [], t)
        np.testing.assert_array_equal(empty, 0.)
        assert balance_report(mesh, empty)['unloaded']
        value = applied_traction(mesh, [load], t)
        unassigned = np.ones(len(mesh.faces), bool); unassigned[ids] = False
        np.testing.assert_array_equal(value[unassigned], 0.)
        assert not balance_report(mesh, value)['balanced']


def test_preflight_finds_sinusoidal_imbalance_even_when_initial_force_is_zero():
    from .load_balance import preflight_balance, LoadBalanceError
    from .load_workflow import FaceLoad
    mesh = cylinder()
    ids = np.flatnonzero(face_geometry(mesh)[1][:,2] > .99)
    load = FaceLoad('top',0.,10.,0.,0.,tuple(ids),mesh.areas[ids].sum(),
                    '0,0,0;0,0,0;0,0,normal_amp*sin(2*pi*f*t)',25.)
    with pytest.raises(LoadBalanceError) as caught:
        preflight_balance(mesh, [load], None, .04)
    data = caught.value.ui_error_data
    assert data['key'] == 'solid.unbalanced_detail'
    assert float(data['values']['time']) > 0
    assert '-' in data['values']['opposite_force']
    bottom = np.flatnonzero(face_geometry(mesh)[1][:,2] < -.99)
    preflight_balance(mesh,[load],correction_operator(mesh,bottom),.04)
