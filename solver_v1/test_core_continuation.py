"""Initialization/provenance tests, not measured Al core or kinetics data."""
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest

from .core_continuation import continue_core_field
from .isolated_screw_core import ScrewFarField
from .nonlocal_interface_elasticity import cubic_elastic_tensor,rotate_elastic_tensor
from .fcc111_geometry import fcc111_geometry_from_b


def fixture(radius=2.,shear=0.):
    geometry=fcc111_geometry_from_b(1.)
    tensor=rotate_elastic_tensor(cubic_elastic_tensor(3.,1.,1.),geometry.plane_basis_in_stacked_cubic_axes())
    far=ScrewFarField(tensor,1.,(.1,.1))
    indices=np.array([[0,0],[1,0],[2,0]])
    xyz=np.array([[0.,0.,0.],[0.,.9,0.],[0.,1.8,0.]])
    length=2.8e-10
    initial=far.displacement(xyz,shear_traction=shear*1e6*length**3/1.602176634e-19)
    core=SimpleNamespace(far_field=far,free_radius=radius,indices=indices,
        free_ids=np.arange(3),xyz=xyz,initial=initial.copy())
    saved={(0,0):np.array([.17,.02,-.01]),(1,0):np.array([.81,-.03,.02])}
    prior=dict(parameter_sha256='test-potential',length_scale_m=length,
        center_over_L0=(.1,.1),radius_over_L0=1.,free_sites=2,shear_traction_MPa=0.)
    return core,prior,saved,dict(parameter_sha256='test-potential',length_scale_m=length,shear_mpa=shear)


def test_expand_preserves_rows_and_far_field_without_interpolation():
    core,prior,saved,kw=fixture(); before=deepcopy(saved)
    before_initial=core.initial.copy()
    result,meta=continue_core_field(core,prior,saved,allow_expansion=True,**kw)
    np.testing.assert_array_equal(core.initial,before_initial)
    np.testing.assert_array_equal(result[:2],list(saved.values()))
    np.testing.assert_array_equal(result[2],core.far_field.displacement(core.xyz)[2])
    for key in saved:
        np.testing.assert_array_equal(saved[key],before[key])
    assert meta['inherited_free_rows']==2 and meta['new_free_rows']==1
    assert meta['initialization_only'] and not meta['domain_convergence_certified']


def test_load_change_applies_only_physical_affine_increment_to_saved_rows():
    core,prior,saved,kw=fixture(shear=25.)
    result,_=continue_core_field(core,prior,saved,allow_expansion=True,**kw)
    affine=core.far_field.displacement(core.xyz,shear_traction=25e6*kw['length_scale_m']**3/1.602176634e-19,burgers_sign=0)
    np.testing.assert_allclose(result[:2],np.array(list(saved.values()))+affine[:2],atol=1e-15)


@pytest.mark.parametrize('change', ['parameter','length','center','missing','implicit','unstable'])
def test_refuse_mismatched_or_unverified_state(change):
    core,prior,saved,kw=fixture(); controls=dict(allow_expansion=True)
    if change=='parameter': kw['parameter_sha256']='different'
    if change=='length': kw['length_scale_m']*=1.01
    if change=='center': prior['center_over_L0']=(.2,.1)
    if change=='missing': saved.pop((1,0))
    if change=='implicit': controls['allow_expansion']=False
    if change=='unstable': controls.update(require_stable=True,prior_summary={})
    with pytest.raises(ValueError):
        continue_core_field(core,prior,saved,**controls,**kw)


def test_same_disk_legacy_copy_and_explicit_stable_gate():
    core,prior,saved,kw=fixture()
    prior['radius_over_L0']=core.free_radius
    saved[(2,0)]=np.array([.96,.01,0.]); prior['free_sites']=3
    summary=dict(final_stability_probe=dict(stable_on_tested_fixed_boundary=True))
    result,meta=continue_core_field(core,prior,saved,require_stable=True,prior_summary=summary,**kw)
    np.testing.assert_array_equal(result,np.array(list(saved.values())))
    assert not meta['domain_expansion'] and meta['new_free_rows']==0


def test_analytic_newton_and_lbfgs_recover_same_stable_perfect_state():
    from .isolated_screw_core import IsolatedScrewCore
    from .run_isolated_screw_core import material
    from .core_stability import lowest_core_mode
    surface,tensor,_=material('historical')
    far=ScrewFarField(tensor,1.,(np.sqrt(3)/12,surface.h/2))
    core=IsolatedScrewCore(surface,far,free_radius=.75,ring=1,burgers_sign=0)
    seed=core.initial+.004*np.sin(np.arange(core.initial.size)).reshape(core.initial.shape)
    results=[core.relax(seed,method=method,force_tolerance=3e-7,max_iterations=65)
             for method in ('lbfgs','newton_cg')]
    for result in results:
        assert result['force_converged']
        assert lowest_core_mode(core,result['field'])['stable_on_tested_fixed_boundary']
    np.testing.assert_allclose(results[0]['field'],results[1]['field'],atol=2e-7)
    assert abs(results[0]['energy']-results[1]['energy'])<2e-13
