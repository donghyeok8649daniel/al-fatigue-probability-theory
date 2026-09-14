import numpy as np
import pytest

from app.specimen_mesh import cylinder, refine_surface, load_surface


def test_default_cylinder_closed_dimensions_and_outward_volume():
    m = cylinder()
    assert m.closed
    np.testing.assert_allclose(np.ptp(m.vertices, axis=0)[2], 30)
    t = m.vertices[m.faces]
    volume = np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2])).sum()/6
    assert volume > 0
    assert volume == pytest.approx(np.pi*25*30, rel=.05)


def test_refinement_preserves_geometry_area_and_closure():
    m = cylinder()
    r = refine_surface(m, 2)
    assert r.closed and len(r.faces) > len(m.faces)
    assert sum(r.areas) == pytest.approx(sum(m.areas), rel=1e-12)
    t = r.vertices[r.faces]
    assert max(np.linalg.norm(t[:, i]-t[:, (i+1)%3], axis=1).max() for i in range(3)) <= 2
    assert not r.vertices.flags.writeable


def test_obj_units_negative_indices_and_triangle_validation(tmp_path):
    p = tmp_path/'test.obj'
    p.write_text('v 0 0 0\nv 1 0 0\nv 0 1 0\nf -3 -2 -1\n')
    m = load_surface(p, 1000)
    assert sum(m.areas) == 500000
    assert not m.closed
    p.write_text('v 0 0 0\nf 1 1 1\n')
    with pytest.raises(ValueError):
        load_surface(p)


def test_ascii_binary_stl_equivalence(tmp_path):
    import struct
    p = tmp_path/'test.stl'
    p.write_text('solid t\n facet normal 0 0 1\n outer loop\n vertex 0 0 0\n vertex 1 0 0\n vertex 0 1 0\n endloop\n endfacet\nendsolid t')
    a = load_surface(p)
    p.write_bytes(b' '*80+struct.pack('<I', 1)+struct.pack('<12fH', 0,0,1,0,0,0,1,0,0,0,1,0,0))
    b = load_surface(p)
    np.testing.assert_array_equal(a.vertices, b.vertices)
    np.testing.assert_array_equal(a.faces, b.faces)


@pytest.mark.parametrize('size', [0, -1, float('nan'), float('inf')])
def test_invalid_mesh_size_rejected(size):
    with pytest.raises(ValueError):
        refine_surface(cylinder(), size)


def test_large_mesh_request_rejected_without_partial_mutation():
    m = cylinder()
    with pytest.raises(ValueError, match='mesh_limit'):
        refine_surface(m, 1e-10)
    assert len(m.faces) < 20000
