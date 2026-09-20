import numpy as np
import pytest

from .load_balance import stress_traction, face_geometry, correction_operator, balanced_traction
from .specimen_mesh import cylinder
from .solid_mechanics import SolidMechanics, FaceStressHistory, check_balance, tensor_from_voigt
from .load_workflow import FaceLoad


@pytest.fixture(scope='module')
def solid():
    return SolidMechanics(cylinder(radius_mm=2, length_mm=8, target_mm=3), 69000., .33, 2.)


def test_all_six_stresses_patch_test_energy_and_rigid_gauge(solid):
    tensor = np.array([[120., 11., -8.], [11., 30., 7.], [-8., 7., 80.]])
    mesh = solid.mesh
    traction = stress_traction(mesh, np.arange(len(mesh.faces)), tensor)
    result = solid.solve(traction)
    expected = np.broadcast_to(tensor, (len(solid.cells), 3, 3))
    np.testing.assert_allclose(tensor_from_voigt(result['stress_voigt_mpa']), expected, rtol=1e-9, atol=1e-8)
    assert result['relative_residual'] < 1e-10
    assert np.linalg.norm(result['gauge_reaction_n']) < 1e-8
    displacement = result['displacement_mm'].ravel()
    np.testing.assert_allclose(solid.rigid.T @ displacement, 0, atol=1e-12)
    assert 2*result['strain_energy_n_mm'] == pytest.approx(
        np.sum(result['displacement_mm']*result['nodal_force_n']), rel=1e-10)


def test_stress_recovery_respects_rigid_motion_before_thin_element_cancellation():
    from .solid_mechanics import elastic_matrix, element_stress
    vertices = np.array([[[0., 0., 0.], [1., .3, 0.], [.4, 1., 0.], [.6, .7, 1e-9]]])
    grad = np.linalg.inv((vertices[:, 1:]-vertices[:, :1]).transpose(0, 2, 1))
    gradients = np.concatenate((-grad.sum(axis=1)[:, None, :], grad), axis=1)
    c = elastic_matrix(69000., .33)
    translated = np.broadcast_to([1000., -2000., 3000.], vertices.shape)
    np.testing.assert_array_equal(element_stress(gradients, translated, c), 0.)
    # A linear displacement has the prescribed six strain components even in
    # a thin tetrahedron; this checks the physical gradient/Voigt convention.
    strain = np.array([.001, -.0002, .0003, .0001, -.0004, .0002])
    xx, yy, zz, yz, xz, xy = strain
    affine = np.array([[xx, xy/2, xz/2], [xy/2, yy, yz/2], [xz/2, yz/2, zz]])
    actual = element_stress(gradients, vertices @ affine.T, c)
    np.testing.assert_allclose(actual[0], c @ strain, atol=2e-5, rtol=1e-7)


def test_sparse_gauge_matches_full_neumann_constraints_for_auxiliary_loads(solid):
    from scipy import sparse
    from scipy.sparse.linalg import splu
    constraint = sparse.csc_matrix(solid.rigid)
    augmented = sparse.bmat([[solid.k/solid.young_mpa, constraint],
                              [constraint.T, sparse.csc_matrix((6, 6))]], format='csc')
    factor = splu(augmented)
    # An arbitrary unbalanced traction also exercises the auxiliary load bases
    # used for superposition. Projecting it must reproduce the old saddle solve.
    traction = np.random.default_rng(512).normal(size=(len(solid.mesh.faces), 3))
    result = solid.solve(traction, require_balance=False)
    force = result['nodal_force_n'].ravel()
    reference = factor.solve(np.r_[force/solid.young_mpa, np.zeros(6)])
    np.testing.assert_allclose(result['displacement_mm'].ravel(), reference[:-6], rtol=1e-8, atol=1e-11)
    np.testing.assert_allclose(result['gauge_reaction_n'].ravel(),
        solid.rigid @ (reference[-6:]*solid.young_mpa), rtol=1e-9, atol=1e-10)
    np.testing.assert_allclose(result['residual_n']+result['gauge_reaction_n'], 0, atol=1e-9)


def test_chunked_assembly_preserves_stiffness(solid, monkeypatch):
    from . import solid_mechanics as sm
    monkeypatch.setattr(sm, 'ASSEMBLY_CHUNK_CELLS', 17)
    chunked = sm.assemble_stiffness(solid.gradients, solid.c, solid.volumes, solid.cells, solid.k.shape[0])
    difference = chunked-solid.k
    assert np.max(np.abs(difference.data), initial=0) < 1e-8


def test_projected_multigrid_matches_direct_neumann_stress_and_auxiliary_gauge(solid, monkeypatch):
    from . import solid_mechanics as sm
    monkeypatch.setattr(sm, 'DIRECT_SOLVER_MAX_NODES', 0)
    iterative = sm.SolidMechanics(solid.mesh, solid.young_mpa, solid.poisson, 2.)
    np.testing.assert_array_equal(iterative.cells, solid.cells)
    tensor = np.array([[120., 11., -8.], [11., 30., 7.], [-8., 7., 80.]])
    patch = stress_traction(solid.mesh, np.arange(len(solid.mesh.faces)), tensor)
    auxiliary = np.random.default_rng(817).normal(size=patch.shape)
    for traction, balanced in ((patch, True), (auxiliary, False)):
        actual = iterative.solve(traction, require_balance=balanced)
        expected = solid.solve(traction, require_balance=balanced)
        np.testing.assert_allclose(actual['displacement_mm'], expected['displacement_mm'], rtol=1e-7, atol=1e-11)
        np.testing.assert_allclose(actual['stress_voigt_mpa'], expected['stress_voigt_mpa'], rtol=1e-7, atol=1e-7)
        np.testing.assert_allclose(actual['gauge_reaction_n'], expected['gauge_reaction_n'], atol=1e-10)
        assert 0 < actual['linear_iterations'] < 2000
        np.testing.assert_allclose(iterative.rigid.T @ actual['displacement_mm'].ravel(), 0, atol=1e-12)
    # A failed iterative solve cannot be presented as a completed mechanics run.
    monkeypatch.setattr(sm, 'cg', lambda *a, **kw: (np.zeros(iterative.k.shape[0]), 2000))
    with pytest.raises(ValueError, match='solid.iterative_failed'): iterative.solve(patch)


def test_multiple_balanced_end_loads_and_unbalanced_rejection(solid):
    mesh = solid.mesh
    _, normals = face_geometry(mesh)
    top = np.flatnonzero(normals[:, 2] > .99)
    bottom = np.flatnonzero(normals[:, 2] < -.99)
    expression = '0,0,0;0,0,0;0,0,normal_mean+normal_amp*sin(2*pi*f*t)'
    loads = [FaceLoad('custom', 100, 30, 0, 0, tuple(ids), mesh.areas[ids].sum(), expression, 25)
             for ids in (top, bottom)]
    history = FaceStressHistory(solid, loads)
    from .surface_setup import encode, decode
    restored_loads, _ = decode(encode(mesh, loads, None), mesh)
    assert restored_loads == loads  # numpy face indices from geometry are portable
    for time in (0., .0037, .01, .037):
        field = history.face_stress(time)
        np.testing.assert_allclose(field[:, 2], 100+30*np.sin(2*np.pi*25*time), atol=1e-8)
        np.testing.assert_allclose(field[:, [0, 1, 3, 4, 5]], 0, atol=1e-8)
    with pytest.raises(ValueError, match='solid.unbalanced'):
        FaceStressHistory(solid, loads[:1]).face_stress(0.)


def test_balanced_force_but_nonzero_torque_is_rejected(solid):
    _, normals = face_geometry(solid.mesh)
    top, bottom = normals[:, 2] > .99, normals[:, 2] < -.99
    traction = np.zeros((len(normals), 3))
    traction[top, 0] = 10; traction[bottom, 0] = -10
    with pytest.raises(ValueError, match='solid.unbalanced'):
        solid.solve(traction)


def test_empty_face_history_is_unloaded_and_does_not_solve_fictitious_bases(solid, monkeypatch):
    monkeypatch.setattr(solid, 'solve', lambda *_a, **_k: pytest.fail('no force basis exists'))
    history = FaceStressHistory(solid, [])
    assert history.cell_basis.shape == (0,len(solid.cells),6)
    for time in (0.,.013,.1):
        c,net,balance,residual = history.coefficients(time)
        assert len(c) == 0 and balance == residual == 0.
        np.testing.assert_array_equal(net,0.)
        np.testing.assert_array_equal(history.face_stress(time),0.)


def test_unloaded_3d_pde_ignores_unassigned_draft_stress_and_roundtrips(tmp_path):
    from .solid_mechanics import run_solid_probability, validate_solid_result, solid_field
    from .solver_adapter import UIAnalysisConfig
    from .project_file import save_bundle, load_bundle
    mesh = cylinder(radius_mm=2.,length_mm=8.,target_mm=3.)
    config = UIAnalysisConfig(stress_mean_mpa=1200.,stress_amplitude_mpa=1500.,
                              model_frequency=1000.,cycles=.01,steps_per_cycle=10)
    reference,result = run_solid_probability(config,mesh,[],None,poisson=.33,target_mm=2.,
                                             direction=[0,0,1],sample_count=2)
    assert reference['initial_stress_mpa'] == 0.
    np.testing.assert_array_equal(reference['applied_stress_mpa'],0.)
    np.testing.assert_array_equal(solid_field(result,'stress',-1),0.)
    np.testing.assert_array_equal(result['force_torque_balance'],0.)
    validate_solid_result(result,mesh)
    save_bundle(tmp_path/'unloaded.ftgsim',result)
    restored=load_bundle(tmp_path/'unloaded.ftgsim')
    validate_solid_result(restored,mesh)
    np.testing.assert_array_equal(restored['local_initiation_probability'],result['local_initiation_probability'])


def test_3d_neck_uses_actual_cross_section_for_stress():
    from .test_axial_specimen import necked_specimen
    mesh = necked_specimen(12)
    system = SolidMechanics(mesh, 69000., .33, 4.)
    _, normals = face_geometry(mesh)
    ends = np.flatnonzero(np.abs(normals[:, 2]) > .99)
    result = system.solve(stress_traction(mesh, ends, np.diag([0., 0., 100.])))
    z = (system.nodes[system.cells].mean(axis=1)+system.origin)[:, 2]
    for mask, expected in (((z > 25) & (z < 35), 400.), ((z > 1) & (z < 5), 100.)):
        mean = np.average(result['stress_voigt_mpa'][mask, 2], weights=system.volumes[mask])
        assert mean == pytest.approx(expected, rel=.005)
    assert result['relative_residual'] < 1e-9


def test_steiner_allowance_counts_added_nodes_and_reports_actual_limit(monkeypatch):
    import tetgen
    from . import solid_mechanics as sm
    mesh = cylinder(radius_mm=2, length_mm=8, target_mm=3)
    monkeypatch.setattr(sm, 'MAX_SOLID_NODES', 200)
    class ExhaustingMesher:
        def __init__(self, vertices, _faces, _markers):
            self.input_nodes = len(vertices)
        def tetrahedralize(self, **options):
            total = self.input_nodes + options['steinerleft']
            return np.zeros((total, 3)), np.zeros((10, 4), dtype=int), None, None
    monkeypatch.setattr(tetgen, 'TetGen', ExhaustingMesher)
    with pytest.raises(sm.SolidMeshLimitError) as error:
        sm.SolidMechanics(mesh, 69000., .33, 3.)
    values = error.value.ui_error_data['values']
    assert values['nodes'] == values['max_nodes'] == 200
    assert values['cells'] == 10 and values['target'] == '3'


def test_approved_time_dependent_correction_is_actually_solved(solid):
    mesh = solid.mesh
    _, normals = face_geometry(mesh)
    top = np.flatnonzero(normals[:, 2] > .99)
    bottom = np.flatnonzero(normals[:, 2] < -.99)
    expression = '0,0,3*cos(2*pi*f*t);0,0,0;3*cos(2*pi*f*t),0,100+30*sin(2*pi*f*t)'
    load = FaceLoad('custom', 100, 30, 0, 0, tuple(top), mesh.areas[top].sum(), expression, 25)
    correction = correction_operator(mesh, bottom)
    history = FaceStressHistory(solid, [load], correction)
    for time in (0., .0037, .019):
        coef, net, balance, linear = history.coefficients(time)
        traction = np.einsum('b,bfi->fi', coef, history.traction_basis)
        direct = solid.solve(traction)
        np.testing.assert_allclose(history.face_stress(time), direct['face_stress_voigt_mpa'], atol=1e-8)
        assert balance < 1e-10 and linear < 1e-10
        np.testing.assert_allclose(net, 0, atol=1e-9)


def test_real_3d_field_drives_actual_pde_and_persists(tmp_path):
    from .solid_mechanics import run_solid_probability, validate_solid_result, solid_field
    from .solver_adapter import UIAnalysisConfig, run_ui_analysis
    from .project_file import save_bundle, load_bundle
    mesh = cylinder(radius_mm=2, length_mm=8, target_mm=3)
    expression = 'normal_mean+normal_amp*sin(2*pi*f*t),15,0;15,80,0;0,0,30'
    load = FaceLoad('all', 500, 100, 0, 0, tuple(range(len(mesh.faces))),
                    mesh.areas.sum(), expression, 1000.)
    config = UIAnalysisConfig(model_frequency=1000., cycles=.2, steps_per_cycle=10,
                               stress_mean_mpa=500., stress_amplitude_mpa=100.)
    records = []
    reference, result = run_solid_probability(config, mesh, [load], None,
        poisson=.33, target_mm=2., direction=[1, 0, 0], sample_count=4,
        record_callback=records.append)
    assert len(records) > 1
    for key in ('model_time', 'local_initiation_probability', 'applied_stress_mpa'):
        np.testing.assert_array_equal([row[key] for row in records], reference[key])
    validate_solid_result(result, mesh)
    assert len(result['sample_cells']) == 1  # homogeneous stress, not element volumes
    scalar = run_ui_analysis(config)
    np.testing.assert_allclose(result['local_initiation_probability'][:, 0],
        np.interp(result['model_time'], scalar['model_time'], scalar['local_initiation_probability']),
        rtol=1e-8, atol=1e-28)
    assert reference['initial_stress_mpa'] == pytest.approx(500)
    for field, value in [('yy', 80), ('xy', 15), ('zz', 30)]:
        np.testing.assert_allclose(solid_field(result, field, 0), value, atol=1e-8)
    assert np.max(result['relative_linear_residual']) < 1e-10
    assert not result['probability_resolution_certified']
    save_bundle(tmp_path/'solid.ftgsim', result)
    restored = load_bundle(tmp_path/'solid.ftgsim')
    validate_solid_result(restored, mesh)
    np.testing.assert_array_equal(solid_field(restored, 'initiation', -1),
                                  solid_field(result, 'initiation', -1))
    # Invalid saved volumes must not silently corrupt the new volume risk maps.
    for key, value in (
        ('cell_volumes_mm3', -result['cell_volumes_mm3']),
        ('cell_volumes_mm3', 2*result['cell_volumes_mm3']),
        ('tetrahedra', result['tetrahedra']+len(result['nodes_mm'])),
        ('local_survival_probability', np.zeros_like(result['local_survival_probability'])),
    ):
        with pytest.raises(ValueError, match='axial.incomplete'):
            validate_solid_result(dict(result, **{key: value}))
