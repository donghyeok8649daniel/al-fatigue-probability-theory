"""Independent conservative-force and identifiable-amplitude controls."""
import numpy as np
import pytest
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_environment_research import environment_jet
from .silicon_sw_material_fit import (force_features,amplitude_parameters,fit_amplitudes,
                                      fit_with_bulk_energy_guard)


@pytest.fixture
def parameters():
    return source_parameters()[0]


def test_full_force_derivatives_and_force_torque_balance(parameters):
    x=np.array([[10.,10.,10.],[12.4,10.2,10.1],[9.4,12.2,10.3],[9.5,9.6,12.3]])
    cell=np.eye(3)*30;pbc=np.zeros(3,bool)
    e,f=force_features(x,cell,pbc,parameters)
    h=2e-5
    for atom,axis in [(0,0),(1,2),(2,1),(3,0)]:
        d=np.zeros_like(x);d[atom,axis]=h
        plus,minus=[force_features(q,cell,pbc,parameters)[0] for q in [x+d,x-d]]
        np.testing.assert_allclose(-(plus-minus)/(2*h),f[atom,axis],atol=2e-8,rtol=3e-7)
    np.testing.assert_allclose(f.sum(axis=0),0.,atol=2e-14)
    for channel in range(3):
        np.testing.assert_allclose(np.cross(x-x.mean(axis=0),f[:,:,channel]).sum(axis=0),0.,atol=2e-13)


def test_centered_triple_and_nonself_periodic_images(parameters):
    # Cell repeats must preserve per-cell energy and atom force, even when the
    # cell is shorter than twice cutoff and an atom sees its own images.
    cell=np.diag([3.1,20.,20.]);x=np.array([[1.,10.,10.],[2.5,10.6,10.2]])
    pbc=np.array([True,False,False])
    e,f=force_features(x,cell,pbc,parameters)
    repeated=np.r_[x,x+cell[0]];large=cell.copy();large[0]*=2
    ee,ff=force_features(repeated,large,pbc,parameters)
    np.testing.assert_allclose(ee,2*e,atol=3e-13)
    np.testing.assert_allclose(ff,np.r_[f,f],atol=3e-13)


def test_component_mapping_preserves_energy_and_force(parameters):
    x=np.array([[10.,10.,10.],[12.4,10.2,10.1],[9.4,12.2,10.3]])
    cell=np.eye(3)*30;pbc=np.zeros(3,bool);a=np.array([.9,1.2,1.4])
    e,f=force_features(x,cell,pbc,parameters)
    mapped=amplitude_parameters(parameters,a)
    ee,ff=force_features(x,cell,pbc,mapped)
    np.testing.assert_allclose(ee.sum(),e@a,atol=2e-13)
    np.testing.assert_allclose(ff.sum(axis=2),np.einsum('nck,k->nc',f,a),atol=2e-13)
    # A separately differentiated explicit site implementation checks total E.
    direct=sum(environment_jet(np.delete(x,center,axis=0)-x[center],mapped,'direct').value
               for center in range(len(x)))
    np.testing.assert_allclose(ee.sum(),direct,atol=2e-13)


def test_positive_fit_and_rank_gate():
    x=np.array([[1.,2.,0.],[0.,1.,3.],[2.,0.,1.],[3.,-1.,2.]])
    truth=np.array([.8,1.1,1.3]);a,record=fit_amplitudes(x,x@truth)
    np.testing.assert_allclose(a,truth,atol=2e-14)
    assert record['rank']==3 and record['kkt_max_abs'] < 1e-12
    with pytest.raises(ValueError,match='rank-deficient'):
        fit_amplitudes(np.ones((5,3)),np.ones(5))


def test_boundary_fit_cannot_become_an_admissible_candidate(parameters):
    a,record=fit_amplitudes(np.eye(3),np.array([-1.,2.,3.]))
    assert a[0] == 0. and record['active_mask'][0] == -1
    with pytest.raises(ValueError,match='positive'):
        amplitude_parameters(parameters,a)


def test_bulk_guard_constrained_optimum_and_infeasibility():
    a,record=fit_with_bulk_energy_guard(np.eye(3),np.ones(3)*.5,
                                       np.array([[1.,1.,0.]]),np.array([2.]),.01)
    np.testing.assert_allclose(a,[.95,.95,.5],atol=3e-8)
    assert record['bulk_guard_active'] and record['bulk_guard_multiplier']>0
    assert record['kkt_max_abs'] < 2e-7
    with pytest.raises(ValueError,match='cannot satisfy'):
        fit_with_bulk_energy_guard(np.eye(3),np.ones(3),np.zeros((1,3)),np.array([2.]),.01)


def test_bulk_guard_boundary_dual_and_inactive_cases():
    # Exact active bound with the first amplitude at its nonnegative boundary.
    a,r=fit_with_bulk_energy_guard(np.eye(3),np.array([-2.,.5,.5]),
                                   np.array([[0.,1.,1.]]),np.array([2.]),.01)
    np.testing.assert_allclose(a,[0.,.95,.95],atol=2e-12)
    assert r['active_mask']==[-1,0,0] and r['kkt_max_abs']<1e-11
    assert r['complementarity_residual']<1e-11
    a,r=fit_with_bulk_energy_guard(np.eye(3),np.array([1.,2.,3.]),
                                   np.eye(3),np.array([1.,2.,3.]),.01)
    np.testing.assert_allclose(a,[1.,2.,3.],atol=1e-14)
    assert not r['bulk_guard_active'] and r['bulk_guard_multiplier']==0.


def test_periodic_strain_work_and_rotational_symmetry(parameters):
    cell=np.array([[3.1,0.,0.],[.2,3.4,0.],[.1,.3,4.2]])
    x=np.array([[.3,.4,.5],[1.7,1.5,2.1]])
    pbc=np.ones(3,bool);h=3e-6
    energy,force,strain=force_features(x,cell,pbc,parameters,strain_derivative=True)
    np.testing.assert_allclose(strain,strain.transpose(1,0,2),atol=3e-13)
    for i,j in [(0,0),(1,1),(2,2),(0,1),(1,2),(2,0)]:
        mode=np.zeros((3,3));mode[i,j]=1.
        eplus,eminus=[force_features(x@(np.eye(3)+s*h*mode).T,
                    cell@(np.eye(3)+s*h*mode).T,pbc,parameters)[0] for s in [1,-1]]
        np.testing.assert_allclose((eplus-eminus)/(2*h),strain[i,j],rtol=2e-7,atol=2e-7)
    # Wrapping atom coordinates changes neither periodic virial nor forces.
    _,ff,ss=force_features(x+np.array([[2,-1,1],[-1,3,0]])@cell,cell,pbc,parameters,
                           strain_derivative=True)
    np.testing.assert_allclose(ff,force,atol=5e-12)
    np.testing.assert_allclose(ss,strain,atol=5e-12)


def test_nonzero_self_images_carry_stress_when_atomic_force_vanishes(parameters):
    cell=np.diag([3.1,20.,20.]);x=np.array([[1.,10.,10.]])
    pbc=np.array([True,False,False]);h=2e-6
    e,f,s=force_features(x,cell,pbc,parameters,strain_derivative=True)
    np.testing.assert_allclose(f,0.,atol=1e-14)
    assert np.linalg.norm(s)>1e-4
    ep,em=[force_features(x*np.array([1+sign*h,1,1]),
                       cell*np.array([1+sign*h,1,1]),pbc,parameters)[0] for sign in [1,-1]]
    np.testing.assert_allclose((ep-em)/(2*h),s[0,0],atol=1e-8,rtol=1e-7)
