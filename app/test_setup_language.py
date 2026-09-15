from dataclasses import replace
import pytest
from .setup_language import emit, compile_setup, analyze_setup
from .load_workflow import FaceLoad
from .specimen_mesh import cylinder
from .tensor_load import default_tensor_expressions


def test_roundtrip_and_metadata():
    mesh = cylinder()
    load = FaceLoad('custom', 100., 20., 0., 0., (0,1), float(mesh.areas[:2].sum()), default_tensor_expressions(), 25.)
    text = emit(mesh, [load, replace(load,frequency=30)])
    loads, correction = compile_setup(text, mesh)
    assert loads == [load,replace(load,frequency=30)] and correction is None
    assert not analyze_setup(mesh, loads, correction)['spatial_solver_connected']
    for bad in [text.replace('MPa','Pa'), text.replace('end','exec code',1),
                text.replace('load L2','load L1'), text.replace('frequency 25.0','frequency nan')]:
        with pytest.raises(ValueError): compile_setup(bad,mesh)


def test_expression_budget_rejects_unbounded_draft_arithmetic():
    from .tensor_load import compile_tensor_matrix, evaluate_tensor
    with pytest.raises(ValueError): compile_tensor_matrix('1e999,0,0;0,0,0;0,0,0')
    with pytest.raises(ValueError): compile_tensor_matrix('True,0,0;0,0,0;0,0,0')
    with pytest.raises(ValueError): compile_tensor_matrix('9'*400+',0,0;0,0,0;0,0,0')
    with pytest.raises(ValueError): compile_tensor_matrix('+'.join(['1']*200)+',0,0;0,0,0;0,0,0')
    compiled = compile_tensor_matrix('9**999999999,0,0;0,0,0;0,0,0')
    with pytest.raises(ValueError):
        evaluate_tensor(compiled,t=0,frequency=1,normal_mean=0,normal_amplitude=0,shear_mean=0,shear_amplitude=0)
