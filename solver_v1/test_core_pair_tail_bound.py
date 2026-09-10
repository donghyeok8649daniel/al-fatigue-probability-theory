import numpy as np
import pytest

from .core_pair_tail_bound import pair_force_tail_bound,explicit_pair_gradient
from .current_core_coefficient_basis import current_core_coefficients
from .current_material_rows import CurrentMaterialScrewCore
from .test_current_material_rows import model,reference_far


def test_pair_force_independent_site_counting(model):
    core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=.85,ring=2)
    q=core.initial+.005*np.sin(np.arange(core.initial.size)).reshape(core.initial.shape)
    basis=current_core_coefficients(core,q)
    actual=explicit_pair_gradient(core,q,u=model.coefficients[0],v=model.coefficients[1])
    np.testing.assert_allclose(actual,basis['gradient'][:,:,:2]@model.coefficients[:2],atol=3e-12,rtol=2e-12)


def test_infinite_bound_encloses_actual_pair_refinement(model):
    values=[]
    for ring in (2,4,8):
        core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=.85,ring=ring)
        q=core.initial+.01*np.cos(np.arange(core.initial.size)).reshape(core.initial.shape)
        values.append(explicit_pair_gradient(core,q,u=model.coefficients[0],v=model.coefficients[1]))
        bound=pair_force_tail_bound(ring=ring,b=core.rows.b,h=core.rows.h,
            u=model.coefficients[0],v=model.coefficients[1],
            maximum_transverse_displacement=float(np.max(np.linalg.norm(q[:,1:],axis=1))))
        if ring==2: coarse=bound['pair_gradient_norm_bound_eV_L0']
        if ring==4: medium=bound['pair_gradient_norm_bound_eV_L0']
        assert bound['environmental_force_tail_certified'] is False
    assert np.max(np.linalg.norm(values[0]-values[-1],axis=1))<coarse
    assert np.max(np.linalg.norm(values[1]-values[-1],axis=1))<medium<coarse


def test_shell_majorant_matches_large_positive_sum():
    from scipy.special import gammaln
    b=1.;h=np.sqrt(2/3);u=.01;v=.04;U=.03;ring=4
    result=pair_force_tail_bound(ring=ring,b=b,h=h,u=u,v=v,maximum_transverse_displacement=U)
    alpha=result['row_map_minimum_singular_value'];k=np.arange(ring+1,10001.)
    radius=alpha*k-2*U;finite=0.
    for c,m in ((u,12),(v,6)):
        nu=(m+2)/2;integral=np.sqrt(np.pi)*np.exp(gammaln(nu-.5)-gammaln(nu))/b
        finite+=np.sum(8*k*np.hypot(b/2,2*U)*c*m*(m+1)*(2*radius**(-m-2)+integral*radius**(-m-1)))
    assert result['pair_gradient_norm_bound_eV_L0']==pytest.approx(finite,rel=3e-14)
    assert result['force_modified'] is False


def test_unsafe_tail_geometry_refused():
    with pytest.raises(ValueError):
        pair_force_tail_bound(ring=2,b=1.,h=1.,u=.1,v=.1,maximum_transverse_displacement=10.)
    with pytest.raises(ValueError):
        pair_force_tail_bound(ring=2,b=1.,h=1.,u=-.1,v=.1,maximum_transverse_displacement=0.)
