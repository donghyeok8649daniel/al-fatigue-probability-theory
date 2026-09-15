import numpy as np
import pytest
from .specimen_mesh import cylinder, refine_selected, SurfaceMesh
from .mesh_probability_view import extrapolate_patches


def test_local_refinement_conserves_geometry_and_conformity():
    mesh = cylinder()
    refined, ids = refine_selected(mesh, [0,1])
    assert refined.closed
    assert len(refined.faces) > len(mesh.faces)
    assert len(refined.faces) < 4*len(mesh.faces)
    assert len(ids) == 12
    assert refined.areas.sum() == pytest.approx(mesh.areas.sum(), rel=1e-14)
    assert refined.enclosed_volume_mm3 == pytest.approx(mesh.enclosed_volume_mm3, rel=1e-14)
    assert mesh.enclosed_volume_mm3 == pytest.approx(np.pi*5**2*30, rel=.05)
    shifted = SurfaceMesh(mesh.vertices+1e6,mesh.faces,'translated')
    assert shifted.enclosed_volume_mm3 == pytest.approx(mesh.enclosed_volume_mm3,rel=1e-9)


def test_open_and_inconsistently_oriented_mesh_refuse_volume():
    mesh = cylinder()
    with pytest.raises(ValueError):
        _ = SurfaceMesh(mesh.vertices,mesh.faces[:-1],'open').enclosed_volume_mm3
    faces = mesh.faces.copy(); faces[0] = faces[0,::-1]
    with pytest.raises(ValueError):
        _ = SurfaceMesh(mesh.vertices,faces,'bad_orientation').enclosed_volume_mm3


@pytest.mark.parametrize('p',[0.,1e-16,.1,1.])
def test_patch_extrapolation_preserves_total_on_refinement(p):
    mesh = cylinder(); finer,_ = refine_selected(mesh,[0,1])
    s,p_init,log_s = extrapolate_patches(p,mesh.areas,.1)
    sf,pf,log_sf = extrapolate_patches(p,finer.areas,.1)
    assert np.exp(log_s.sum()) == pytest.approx(np.exp(log_sf.sum()),abs=1e-15)
    np.testing.assert_allclose(s+p_init,1)
    assert np.all(pf>=0)
    # No floor/certification is manufactured by this mathematical function.
    if p == 1e-16: assert p_init.max() > 0


def test_patch_hazard_extreme_ratio_does_not_create_nan():
    s,p,_ = extrapolate_patches(0.,[1.],1e-320)
    np.testing.assert_array_equal(s,[1.]); np.testing.assert_array_equal(p,[0.])
    s,p,_ = extrapolate_patches(1e-300,[1.],1e-300)
    assert s[0] == pytest.approx(np.exp(-1))
