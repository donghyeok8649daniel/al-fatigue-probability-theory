import numpy as np
import pytest

from .core_farfield_load_audit import affine_row_traction,nominal_shear_audit
from .current_material_rows import CurrentMaterialScrewCore
from .test_current_material_rows import model,reference_far
from .run_source_core_reference import load_source_material,build_source_core
from .isolated_screw_core import ScrewFarField


def test_affine_energy_traction_derivative_and_units(model):
    core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=.85,ring=7,burgers_sign=0)
    gamma=np.array([.0003,.001]);value=affine_row_traction(core,gamma);step=2e-7
    for j in range(2):
        delta=np.eye(2)[j]*step
        derivative=(affine_row_traction(core,gamma+delta)['energy_change_per_atom_eV']
                    -affine_row_traction(core,gamma-delta)['energy_change_per_atom_eV'])/(2*step)
        assert derivative==pytest.approx(value['Piola_yz_columns_eV_L0cubed'][0,j]
            *value['atomic_volume_over_L0cubed'],abs=2e-8,rel=2e-6)
    tiny=nominal_shear_audit(core,.01,length_scale_m=4.05/np.sqrt(2)*1e-10)
    assert tiny['sigma_xz_error_MPa']==pytest.approx(0.,abs=2e-8)
    assert abs(tiny['actual_traction_columns_MPa'][0,0])<2e-8
    assert tiny['nonlinear_stress_control'] is False


def test_source_nominal_load_uses_own_tensor_without_changing_pressure():
    model,tensor,meta=load_source_material()
    far=ScrewFarField(tensor,model.geometry.b,(np.sqrt(3)/12,model.h/2))
    core=build_source_core(model,far,free_radius=.85,ring=5,burgers_sign=0)
    result=nominal_shear_audit(core,.01,length_scale_m=meta['length_scale_m'])
    assert abs(result['sigma_xz_error_MPa'])<2e-8
    assert abs(result['zero_load_columns_MPa'][2,1])>1e-8
    np.testing.assert_allclose(result['actual_traction_columns_MPa']-result['zero_load_columns_MPa'],
                               result['traction_increment_columns_MPa'],atol=0,rtol=0)


def test_bad_affine_strain_and_scale_refused(model):
    core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=.85,ring=1,burgers_sign=0)
    with pytest.raises(ValueError): affine_row_traction(core,[0.,0.,0.])
    with pytest.raises(ValueError): nominal_shear_audit(core,50.,length_scale_m=0.)
