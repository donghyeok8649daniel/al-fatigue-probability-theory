import numpy as np
import pytest

from .current_core_coefficient_basis import current_core_coefficients
from .core_force_identifiability import force_compatibility_lower_bound,constrained_force_profile
from .current_material_rows import CurrentMaterialScrewCore
from .frozen_core_force_target import FrozenCoreForceTarget
from .test_current_material_rows import model,reference_far


def test_full_current_energy_and_force_coefficient_identity(model):
    core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=.85,ring=2)
    q=core.initial+.004*np.cos(np.arange(core.initial.size)*.33).reshape(core.initial.shape)
    columns=current_core_coefficients(core,q)
    value=core.evaluate(q)
    assert columns['energy']@model.coefficients == pytest.approx(value['energy'],abs=3e-12)
    np.testing.assert_allclose(columns['gradient']@model.coefficients,value['gradient'],atol=2e-11,rtol=2e-11)
    v=np.sin(np.arange(q.size)*.7).reshape(q.shape); v/=np.linalg.norm(v)
    step=1e-6
    plus=current_core_coefficients(core,q+step*v)
    minus=current_core_coefficients(core,q-step*v)
    np.testing.assert_allclose((plus['energy']-minus['energy'])/(2*step),
        np.einsum('nic,ni->c',columns['gradient'],v),atol=3e-8,rtol=3e-7)


def test_force_lower_bound_preserves_exact_constraints_and_reports_failure():
    design=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[1.,1.,1.]])
    target=np.array([0.,0.,0.,10.])
    exact=np.array([[1.,1.,0.]])
    baseline=np.array([.3,.7,2.])
    result=force_compatibility_lower_bound(design,target,exact,[1.],baseline,scales=[.1,2.,100.])
    assert result['exact_rank']==1 and result['null_dimension']==2
    assert exact@result['coefficients']==pytest.approx([1.])
    assert result['unconstrained_lower_bound_force_rms']>0
    assert result['unconstrained_lower_bound_force_rms']<result['initial_force_rms']
    assert result['least_squares_normal_residual']<1e-10
    assert not result['material_accepted'] and not result['physical_inequalities_enforced']


def test_identifiable_exact_problem_and_bad_shapes():
    result=force_compatibility_lower_bound(np.eye(2),[1.,2.],np.eye(2),[1.,2.],[.9,1.9],scales=[1.,1.])
    assert result['null_dimension']==0 and result['unconstrained_lower_bound_force_rms']<1e-14
    with pytest.raises(ValueError):
        force_compatibility_lower_bound(np.eye(2),[1.,2.],np.eye(2),[1.,2.],[1.,2.],scales=[0.,1.])
    with pytest.raises(ValueError,match='inconsistent'):
        force_compatibility_lower_bound(np.eye(2),[1.,2.],[[1.,0.],[1.,0.]],[1.,2.],
                                       [1.,2.],scales=[1.,1.])


def test_sign_constrained_core_fit_cannot_beat_relaxed_lower_bound():
    design=np.eye(3); target=np.array([.8,.2,-2.])
    exact=np.array([[1.,1.,0.]]); rhs=np.array([1.]); baseline=np.array([.5,.5,1.])
    lower=force_compatibility_lower_bound(design,target,exact,rhs,baseline,scales=np.ones(3))
    fit=constrained_force_profile(design,target,exact,rhs,nonnegative=(0,1,2))
    np.testing.assert_allclose(fit['coefficients'],[.8,.2,0.],atol=1e-10)
    assert fit['force_rms']>=lower['unconstrained_lower_bound_force_rms']
    assert fit['exact_residual']<1e-10 and fit['kkt_residual']<1e-7
    assert not fit['material_accepted'] and not fit['relaxed_core_recomputed']


def test_frozen_interior_force_matches_full_source_atomic_field(model):
    core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=2.,ring=2)
    field=core.initial+.003*np.sin(np.arange(core.initial.size)*.17).reshape(core.initial.shape)
    full=current_core_coefficients(core,field)
    probe=FrozenCoreForceTarget(core,field,radius=.85)
    small=probe.coefficient_matrix(model,ring=2)
    lookup={tuple(i):k for k,i in enumerate(core.indices[core.free_ids])}
    selected=[lookup[tuple(i)] for i in small['logical_rows']]
    np.testing.assert_allclose(small['design'],full['gradient'][selected].reshape(-1,10),
                               atol=4e-12,rtol=1e-11)
    np.testing.assert_allclose(small['design']@model.coefficients,small['target'],atol=2e-12)
