import numpy as np
import pytest

from .density_angular_cross_material import (
    CrossDensityAngularInterface,CrossDensityAngularCache,CrossDensityAngularBulk,
    cross_site_jet,cross_sector_profile,nonlinear_sector_matrix,
)
from .quartic_angular_material import QuarticSymmetryInterface
from .vector_material_calibration import MaterialObservation


def test_cross_is_per_site_not_product_of_global_sums():
    rho=np.zeros((2,10));rho[:,0]=[.1,.3]
    tensor=np.zeros((2,10,2));tensor[:,0]=[[1,0],[0,2]]
    jet=cross_site_jet(rho,tensor)
    assert jet[0]==pytest.approx(.1+1.2)
    assert jet[0]!=pytest.approx(.4*5)
    assert np.array_equal(jet[1:],np.zeros(9))


def test_cross_recovers_old_surface_and_has_no_pristine_harmonic_part():
    decays=(2.8,5.6,11.9);c=np.array([.33,1.5,10.,16.,.3,4.,8.,.8,.26,1e5])
    old=QuarticSymmetryInterface(decays,c)
    new=CrossDensityAngularInterface(decays,np.r_[c,0.])
    state=(new.h*1.06,.21,.10)
    assert new.evaluate(state).energy==old.evaluate(state).energy
    assert np.array_equal(new.evaluate(state).hessian,old.evaluate(state).hessian)
    jet,_=new.cross_jet((new.h,0.,0.))
    assert np.max(abs(jet))<1e-20
    cols,tails=CrossDensityAngularBulk(decays).evaluate((.2,.1,.05))
    assert np.array_equal(cols[-1],np.zeros((3,3))) and tails[-1]==0


def test_cross_analytic_jet_periodicity_and_nonlinear_sector_bound():
    decays=(2.8,5.6,11.9);c=np.array([.33,1.5,10.,16.,.3,4.,8.,.8,.26,1e5,200.])
    model=CrossDensityAngularInterface(decays,c)
    q=np.array([model.h*1.07,.23,.12]);value=model.evaluate(q)
    gradients=[];hessians=[];step=2e-5
    for axis in range(3):
        delta=np.eye(3)[axis]*step
        a=model.evaluate(q+delta);b=model.evaluate(q-delta)
        gradients.append((a.energy-b.energy)/(2*step));hessians.append((a.gradient-b.gradient)/(2*step))
    assert np.allclose(gradients,value.gradient,rtol=2e-6,atol=5e-7)
    assert np.allclose(np.array(hessians).T,value.hessian,rtol=3e-6,atol=2e-5)
    periodic=model.evaluate(q+[0.,1.,0.])
    assert periodic.energy==pytest.approx(value.energy,abs=1e-12)
    with pytest.raises(ValueError):CrossDensityAngularInterface(decays,np.r_[c[:10],1000.])
    for y in (-.8,-.1,.2,1.):
        for I in (0.,1e-4,.01,1.):
            assert np.array([y,I])@nonlinear_sector_matrix(c)@np.array([y,I])>=0


def test_convex_sector_cut_solves_known_projection_without_relaxing_constraint():
    targets=np.ones(11);targets[10]=6.
    obs=[MaterialObservation(str(i),float(y),1.,'test','exact' if i<2 else 'fit') for i,y in enumerate(targets)]
    columns=np.zeros((1,11,3,3));columns[0,0]=np.eye(3)
    result=cross_sector_profile(np.eye(11),obs,columns,np.zeros((1,11)),nonnegative=(0,1,2,4,5,8,9))
    assert result['coefficients'][4]==pytest.approx(7/3,abs=1e-8)
    assert result['coefficients'][9]==pytest.approx(7/3,abs=1e-8)
    assert result['coefficients'][10]==pytest.approx(14/3,abs=1e-8)
    assert result['environment_sector_verified'] and result['kkt_residual']<1e-6


def test_cross_observation_is_actual_energy_column_not_surrogate():
    decays=(2.8,5.6,11.9);c=np.array([.33,1.5,10.,16.,.3,4.,8.,.8,.26,1e5,200.])
    model=CrossDensityAngularInterface(decays,c)
    state=(model.h*1.07,.23,.12)
    obs=[MaterialObservation('energy',0,1,'eV','fit',state,tuple(np.eye(10)[0]))]
    matrix=CrossDensityAngularCache(obs).matrix(decays)
    assert float((matrix@c)[0])==pytest.approx(model.evaluate(state).energy,abs=1e-11)
