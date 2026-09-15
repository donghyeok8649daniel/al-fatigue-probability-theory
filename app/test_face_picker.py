import numpy as np
from app.face_picker import ray_triangle_pick, planar_patch
from app.specimen_mesh import cylinder


def test_ray_picks_visible_surface_and_misses_empty_space():
    mesh = cylinder()
    i = ray_triangle_pick(mesh.vertices, mesh.faces, np.array([1., 0., 50.]), np.array([0., 0., -1.]))
    assert i is not None
    assert np.all(mesh.vertices[mesh.faces[i], 2] == 30)
    assert ray_triangle_pick(mesh.vertices, mesh.faces, np.array([30., 0., 50.]), np.array([0., 0., -1.])) is None
    patch = planar_patch(mesh, i)
    assert len(patch) > 1
    assert np.all(mesh.vertices[mesh.faces[patch], 2] == 30)


def test_bottom_patch_is_not_top_and_patch_does_not_wrap_cylinder():
    mesh = cylinder()
    bottom = ray_triangle_pick(mesh.vertices, mesh.faces, np.array([1., 0., -10.]), np.array([0., 0., 1.]))
    assert np.all(mesh.vertices[mesh.faces[planar_patch(mesh, bottom)], 2] == 0)
    assert len(planar_patch(mesh, 0)) < len(mesh.faces)/2
