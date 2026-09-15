import pytest
from .surface_setup import encode, decode
from .load_workflow import FaceLoad
from .specimen_mesh import cylinder
from .tensor_load import default_tensor_expressions


def test_setup_roundtrip_and_mesh_binding():
    mesh = cylinder()
    load = FaceLoad('custom', 100, 20, 0, 0, (0,1), float(mesh.areas[:2].sum()), default_tensor_expressions(), 25)
    text = encode(mesh, [load], None)
    loads, correction = decode(text, mesh)
    assert loads == [load] and correction is None
    with pytest.raises(ValueError): decode(text, cylinder(radius_mm=6))
    with pytest.raises(ValueError): decode(text.replace('MPa', 'Pa'), mesh)
    with pytest.raises(ValueError): decode(text.replace('normal_mean+', '__import__+'), mesh)
