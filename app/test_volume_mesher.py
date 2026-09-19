import numpy as np
import pytest

from .specimen_mesh import cylinder
from .volume_mesher import tetrahedralize, SolidMeshError


@pytest.mark.parametrize('failure', ['engine', 'degenerate', 'boundary'])
def test_invalid_first_mesh_retries_without_changing_input_or_validation(monkeypatch, failure):
    import tetgen
    original = tetgen.TetGen
    options = []
    mesh = cylinder(radius_mm=2., length_mm=8., target_mm=3.)
    vertices, faces = mesh.vertices.copy(), mesh.faces.copy()
    class Trial:
        def __init__(self, *args): self.real = original(*args)
        def tetrahedralize(self, **kwargs):
            options.append(kwargs)
            self.first = not kwargs['nobisect']
            if self.first and failure == 'engine': raise RuntimeError('fixture engine failure')
            output = list(self.real.tetrahedralize(**kwargs))
            if self.first and failure == 'degenerate':
                output[1] = output[1].copy(); output[1][0] = output[1][0,0]
            return output
        @property
        def triface_markers(self):
            markers = self.real.triface_markers.copy()
            if self.first and failure == 'boundary':
                i = np.flatnonzero(markers > 0)[0]
                markers[i] = 1 if markers[i] != 1 else 2
            return markers
        def __getattr__(self, name): return getattr(self.real, name)
    monkeypatch.setattr(tetgen, 'TetGen', Trial)
    result = tetrahedralize(mesh, 2., mesh.vertices.mean(axis=0))
    assert [o['nobisect'] for o in options] == [False, True]
    assert all(o['quality'] and o['maxvolume'] == 8/6 for o in options)
    assert all(o['steinerleft'] > 0 for o in options)
    assert result['meshing_info']['mode'] == 'preserved_input_facets'
    assert result['volumes'].min() > 0
    assert result['volumes'].sum() == pytest.approx(mesh.enclosed_volume_mm3, rel=1e-10)
    mapped = np.bincount(result['parents'], weights=result['boundary_areas'], minlength=len(faces))
    np.testing.assert_allclose(mapped, mesh.areas, rtol=1e-10)
    np.testing.assert_array_equal(mesh.vertices, vertices)
    np.testing.assert_array_equal(mesh.faces, faces)


def test_engine_failure_retains_detail_and_no_invalid_mesh_is_returned(monkeypatch):
    import tetgen
    class Failure:
        def __init__(self, *_): pass
        def tetrahedralize(self, **_): raise RuntimeError('fixture-engine')
    monkeypatch.setattr(tetgen, 'TetGen', Failure)
    mesh = cylinder()
    with pytest.raises(SolidMeshError) as caught:
        tetrahedralize(mesh, 3., mesh.vertices.mean(axis=0))
    assert caught.value.ui_error_data == dict(key='solid.mesh_engine', values=dict(detail='fixture-engine'))
