"""Independent section geometry, force balance, and actual local-PDE checks."""
from dataclasses import replace

import numpy as np
import pytest

from .axial_specimen import (prepare_axial_sections, run_axial_probability,
                            validate_spatial_result)
from .solver_adapter import UIAnalysisConfig, run_ui_analysis
from .specimen_mesh import SurfaceMesh, cylinder, refine_selected, refine_surface


def necked_specimen(n=24):
    """Polygonal solid of revolution, r=5 at grips and r=2.5 in the neck."""
    z = np.array([0., 10., 20., 40., 50., 60.])
    radius = np.array([5., 5., 2.5, 2.5, 5., 5.])
    angle = np.arange(n)*2*np.pi/n
    vertices = [[r*np.cos(a), r*np.sin(a), x]
                for x, r in zip(z, radius) for a in angle]
    vertices += [[0, 0, z[0]], [0, 0, z[-1]]]
    faces = []
    for j in range(len(z)-1):
        for i in range(n):
            a, b = j*n+i, j*n+(i+1) % n
            faces.extend([[a, b, b+n], [a, b+n, a+n]])
    for i in range(n):
        faces.extend([[len(vertices)-2, (i+1) % n, i],
                      [len(vertices)-1, (len(z)-1)*n+i, (len(z)-1)*n+(i+1) % n]])
    mesh = SurfaceMesh(vertices, faces, 'synthetic_necked_specimen')
    return mesh


def bottom(mesh):
    return np.flatnonzero(np.all(mesh.vertices[mesh.faces, 2] == mesh.vertices[:, 2].min(), axis=1))


def test_constant_bar_and_neck_force_conservation():
    bar = cylinder(radius_mm=5, length_mm=60, target_mm=3)
    straight = prepare_axial_sections(bar, bottom(bar), 12)
    np.testing.assert_allclose(straight.stress_factors, 1, rtol=2e-14)
    mesh = necked_specimen()
    section = prepare_axial_sections(mesh, bottom(mesh), 12)
    radius = np.interp(section.positions_mm, [0, 10, 20, 40, 50, 60], [5, 5, 2.5, 2.5, 5, 5])
    polygon_factor = 24*np.sin(2*np.pi/24)/2
    np.testing.assert_allclose(section.areas_mm2, polygon_factor*radius**2, rtol=1e-14)
    np.testing.assert_allclose(section.stress_factors, (5/radius)**2, rtol=1e-14)
    stress = 100*section.stress_factors
    np.testing.assert_allclose(stress*section.areas_mm2, 100*section.loaded_area_mm2, rtol=1e-14)
    assert stress[5] == pytest.approx(400)
    assert stress[0] == pytest.approx(100)


def test_rotation_translation_scaling_and_triangle_refinement():
    mesh = necked_specimen()
    original = prepare_axial_sections(mesh, bottom(mesh), 12)
    axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
    cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rotation = np.eye(3)+np.sin(.7)*cross+(1-np.cos(.7))*(cross@cross)
    moved = SurfaceMesh(mesh.vertices @ rotation.T * .5 + 1e5, mesh.faces, 'moved')
    actual = prepare_axial_sections(moved, bottom(mesh), 12)
    np.testing.assert_allclose(actual.areas_mm2, original.areas_mm2/4, rtol=1e-10)
    np.testing.assert_allclose(actual.stress_factors, original.stress_factors, rtol=1e-10)
    refined, _ = refine_selected(mesh, [0, 1, 30])
    fine = prepare_axial_sections(refined, bottom(refined), 12)
    np.testing.assert_allclose(fine.areas_mm2, original.areas_mm2, rtol=1e-14)
    # The same physical stations are also present on the 3x finer axial grid.
    finer = prepare_axial_sections(mesh, bottom(mesh), 36)
    np.testing.assert_allclose(finer.areas_mm2[1::3], original.areas_mm2, rtol=1e-14)


def test_reject_open_side_partial_end_and_offcenter_sections():
    mesh = necked_specimen()
    with pytest.raises(ValueError, match='axial.mesh'):
        prepare_axial_sections(SurfaceMesh(mesh.vertices, mesh.faces[:-1], 'open'), bottom(mesh))
    for ids in (bottom(mesh)[:-1], np.array([0, 1]), np.arange(len(mesh.faces))):
        with pytest.raises(ValueError, match='axial.end'):
            prepare_axial_sections(mesh, ids)
    shifted = mesh.vertices.copy(); shifted[48:96, 0] += 1
    with pytest.raises(ValueError, match='axial.centered'):
        prepare_axial_sections(SurfaceMesh(shifted, mesh.faces, 'bent'), bottom(mesh))


def test_actual_pdes_receive_stress_not_patch_area_and_keep_floors(tmp_path):
    from .project_file import save_bundle, load_bundle
    mesh = necked_specimen()
    sections = prepare_axial_sections(mesh, bottom(mesh), 6)
    config = UIAnalysisConfig(stress_mean_mpa=500, stress_amplitude_mpa=100,
                              model_frequency=1000, cycles=.2, steps_per_cycle=10)
    result = run_axial_probability(config, sections)
    validate_spatial_result(result, mesh)
    # Independently run the existing PDE at the 4x neck load, not a scaled P.
    neck = run_ui_analysis(replace(config, stress_mean_mpa=2000, stress_amplitude_mpa=400))
    np.testing.assert_allclose(result['local_initiation_probability'][:, 2],
                               np.interp(result['model_time'], neck['model_time'],
                                         neck['local_initiation_probability']), rtol=1e-10, atol=1e-30)
    assert result['local_initiation_probability'][-1, 2] > result['local_initiation_probability'][-1, 0]
    assert result['force_balance_residual_n'] < 1e-10
    assert not result['probability_resolution_certified']
    assert not result['spatial_resolution_certified']
    assert np.max(np.abs(result['mass_balance_residual'])) < 1e-10
    assert np.all(result['local_rare_event_floor'] >= 0)
    path = tmp_path/'spatial.ftgsim'
    save_bundle(path, result)
    restored = load_bundle(path)
    validate_spatial_result(restored, mesh)
    np.testing.assert_array_equal(restored['local_initiation_probability'], result['local_initiation_probability'])
    altered = SurfaceMesh(mesh.vertices*2, mesh.faces, 'changed')
    with pytest.raises(ValueError, match='axial.stale'):
        validate_spatial_result(restored, altered)


def test_stop_never_returns_a_partial_spatial_map():
    mesh = necked_specimen()
    sections = prepare_axial_sections(mesh, bottom(mesh), 6)
    with pytest.raises(InterruptedError, match='axial.cancelled'):
        run_axial_probability(UIAnalysisConfig(), sections, stop_requested=lambda: True)
