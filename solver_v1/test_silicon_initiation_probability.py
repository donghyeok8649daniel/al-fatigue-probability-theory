import numpy as np
import pytest
from scipy.linalg import expm

from solver_v1.silicon_charge_dynamics import sg_chain
from solver_v1.silicon_initiation_probability import (
    check_generator, initiation_collectors, evolve_first_passage,
    basin_committor, first_passage_moments)


def generator(n,edges):
    g=np.zeros((n,n))
    for source,target,rate in edges:g[target,source]+=rate;g[source,source]-=rate
    return g


def test_competing_causes_against_closed_form():
    g=generator(3,[(0,1,.3),(0,2,.7),(1,0,100),(2,0,200)])
    model=initiation_collectors(g,{'surface':[1],'interface':[2]})
    time=np.array([0.,.01,.5,1.,5.])
    result=evolve_first_passage(model,[1.],time)
    assert result['survival']==pytest.approx(np.exp(-time),abs=3e-15)
    assert result['cause_cumulative']==pytest.approx((1-np.exp(-time))[:,None]*[.3,.7],abs=3e-15)
    assert np.max(abs(result['mass_residual']))<3e-15
    assert result['minimum_probability']>=0
    moments=first_passage_moments(model)
    assert moments['mean_first_passage_time']==pytest.approx([1.])
    assert moments['eventual_cause_probability'][:,0]==pytest.approx([.3,.7])


def test_returns_from_intermediate_are_not_counted_as_initiation():
    g=generator(5,[(0,1,2),(1,0,3),(0,2,1),(2,0,4),(1,3,.4),(2,4,.5),(3,1,7),(4,2,9)])
    h=basin_committor(g,intact_core=[0],initiated_basins=[3,4])
    assert h['committor']==pytest.approx([0,.4/3.4,.5/4.5,1,1])
    assert np.max(abs(h['interior_residual']))<2e-16
    true=initiation_collectors(g,{'surface':[3],'interface':[4]})
    early=initiation_collectors(g,{'transient_opening':[1,2,3,4]})
    x=evolve_first_passage(true,[1,0,0],[.1,1])
    y=evolve_first_passage(early,[1],[.1,1])
    assert np.all(x['survival']>y['survival'])
    tau=(1+2/3.4+1/4.5)/(2*.4/3.4+1*.5/4.5)
    assert first_passage_moments(true)['mean_first_passage_time'][0]==pytest.approx(tau)


def test_sg_boundary_absorption_against_independent_dense_exponential():
    x=np.linspace(-1,1,31)
    g=sg_chain(.7*(x*x-.36)**2,thermal_energy=.2,spacing=x[1]-x[0],face_mobility=.4)
    model=initiation_collectors(g,{'left_surface':[0],'right_interface':[30]})
    p=np.exp(-40*x[1:-1]**2);p/=p.sum()
    result=evolve_first_passage(model,p,[0,.1,.7,3.])
    full=np.r_[p,0.,0.]
    for i,t in enumerate([0,.1,.7,3.]):
        independent=expm(model['generator'].toarray()*t)@full
        assert result['transient_mass'][i]==pytest.approx(independent[:-2],abs=1e-14)
        assert result['cause_cumulative'][i]==pytest.approx(independent[-2:],abs=1e-14)
    assert np.max(abs(result['cause_cumulative'][:,0]-result['cause_cumulative'][:,1]))<1e-14
    assert np.max(abs(result['mass_residual']))<3e-14
    assert np.all(np.diff(result['survival'])<=0)
    assert np.all(np.diff(result['cause_cumulative'],axis=0)>=0)


def test_dividing_a_cause_into_subcauses_preserves_total_probability():
    g=generator(4,[(0,1,1),(1,0,2),(0,2,.2),(1,3,.4)])
    separate=initiation_collectors(g,{'surface':[2],'interface':[3]})
    together=initiation_collectors(g,{'any':[2,3]})
    a=evolve_first_passage(separate,[1,0],[0,.3,2,20])
    b=evolve_first_passage(together,[1,0],[0,.3,2,20])
    assert a['survival']==pytest.approx(b['survival'],abs=1e-14)
    assert a['cause_cumulative'].sum(axis=1)==pytest.approx(b['cause_cumulative'][:,0],abs=1e-14)


def test_frozen_chemical_mixture_preserves_conditional_first_passage():
    # Two frozen chemical configurations; the specimen never switches between them.
    # These rates are dimensionless synthetic controls, not Si or dopant estimates.
    rates=np.array([.02,2.]);weights=np.array([.9,.1]);time=np.array([0.,.1,1.,10.])
    g=generator(4,[(0,2,rates[0]),(1,3,rates[1])])
    model=initiation_collectors(g,{'first_crack':[2,3]})
    result=evolve_first_passage(model,weights,time)
    conditional=weights[None,:]*np.exp(-time[:,None]*rates)
    exact=conditional.sum(axis=1)
    assert result['survival']==pytest.approx(exact,abs=3e-15)
    assert result['cause_cumulative'][:,0]==pytest.approx(1-exact,abs=3e-15)
    hazard=result['cause_flux'][:,0]/result['survival']
    assert hazard==pytest.approx((conditional@rates)/exact,abs=3e-15)
    assert np.all(np.diff(hazard)<0)
    assert exact[-1]>3*np.exp(-time[-1]*(weights@rates))
    assert weights@first_passage_moments(model)['mean_first_passage_time']==pytest.approx(weights@(1/rates))


def test_uniform_diffusion_committor_mfpt_and_continuum_refinement():
    diffusion,time=.2,2.
    k=np.arange(60);odd=2*k+1
    survival_exact=4/np.pi*np.sum((-1.)**k/odd*np.exp(-diffusion*(odd*np.pi/2)**2*time))
    errors=[]
    for count in (31,61,121):
        x=np.linspace(-1,1,count)
        g=sg_chain(np.zeros(count),thermal_energy=1.,spacing=x[1]-x[0],face_mobility=diffusion)
        model=initiation_collectors(g,{'left':[0],'right':[count-1]})
        h=basin_committor(g,intact_core=[0],initiated_basins=[count-1])['committor']
        assert h==pytest.approx((x+1)/2,abs=3e-14)
        tau=first_passage_moments(model)['mean_first_passage_time']
        assert tau==pytest.approx((1-x[1:-1]**2)/(2*diffusion),abs=3e-12)
        p=np.zeros(count-2);p[(count-2)//2]=1
        got=evolve_first_passage(model,p,[time])['survival'][0]
        errors.append(abs(got-survival_exact))
    assert 3.8<errors[0]/errors[1]<4.2
    assert 3.8<errors[1]/errors[2]<4.2


def test_closed_class_is_not_regularized_into_a_rate():
    g=generator(4,[(0,1,1),(1,0,1),(2,3,1)])
    model=initiation_collectors(g,{'crack':[3]})
    assert not model['all_states_can_reach_B']
    with pytest.raises(ValueError,match='closed class'):first_passage_moments(model)
    with pytest.raises(ValueError,match='unable'):basin_committor(g,intact_core=[2],initiated_basins=[3])


def test_unresolved_rare_absorption_is_not_given_a_finite_clock():
    g=generator(3,[(0,1,1),(1,0,1),(1,2,1e-16)])
    model=initiation_collectors(g,{'crack':[2]})
    assert model['all_states_can_reach_B']
    with pytest.raises(ValueError,match='numerically unresolved'):
        first_passage_moments(model)


def test_overlapping_cause_labels_are_rejected():
    g=generator(3,[(0,1,1),(0,2,1)])
    with pytest.raises(ValueError,match='disjoint'):
        initiation_collectors(g,{'surface':[1],'interface':[1,2]})


@pytest.mark.parametrize('g',[np.array([[-1,0],[.5,0]]),np.array([[-1,-.1],[1,.1]])])
def test_invalid_generator_rejected(g):
    with pytest.raises(ValueError):check_generator(g)


def test_conservation_validation_is_invariant_to_time_units():
    good=generator(3,[(0,1,.3),(1,0,.7),(1,2,.2)])
    bad=good.copy();bad[0,0]*=.5
    for scale in (1e-15,1.,1e15):
        checked=check_generator(scale*good)
        assert checked.toarray()==pytest.approx(scale*good)
        with pytest.raises(ValueError,match='conserve'):
            check_generator(scale*bad)
    # Zero-rate states have exactly zero residual and remain valid.
    assert check_generator(np.zeros((2,2))).nnz==0


def test_fast_columns_do_not_hide_missing_probability_in_slow_columns():
    good=generator(4,[(0,1,1.),(1,0,2.),(2,3,1e-15),(3,2,2e-15)])
    check_generator(good)
    bad=good.copy();bad[2,2]*=.5
    with pytest.raises(ValueError,match='conserve'):
        check_generator(bad)
