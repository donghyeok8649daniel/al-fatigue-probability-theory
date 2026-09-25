"""Detailed balance, fast-limit projection and quenched-mixture counterexample."""
import numpy as np
import pytest
from scipy.linalg import expm
from numpy.testing import assert_allclose

from .silicon_charge_dynamics import charge_memory_blocks, coupled_charge_generator, reversible_propagator, sg_chain


def setup(cells=48, rate=30.):
    dx = 2./cells
    q = -1+dx*(np.arange(cells)+.5)
    f = np.array([.5*(q+.35)**2, .8*(q-.25)**2+.03])
    r = coupled_charge_generator(f, excess_counts=[0, 1], chemical_potential=.02,
                                 thermal_energy=.12, spacing=dx, mobilities=[1., .3], charge_rate=rate)
    initial = np.exp(-((q+.15)/.2)**2)
    initial /= sum(initial)
    return q, initial, r


def test_joint_probability_is_conserved_positive_reversible_and_gibbs_stationary():
    _, p, r = setup()
    l, pi = r['generator'], r['equilibrium']
    assert_allclose(np.asarray(l.sum(axis=0)), 0, atol=1e-12)
    assert_allclose(l@pi, 0, atol=2e-14)
    conductance = l.toarray()*pi
    assert_allclose(conductance, conductance.T, atol=2e-14)
    t = reversible_propagator(l, pi, .3)
    state = t@(r['lift']@p)
    assert state.min() > 0
    assert_allclose(sum(state), 1, atol=1e-12)
    assert_allclose(t, expm(.3*l.toarray()), rtol=2e-9, atol=1e-12)


def test_projection_preserves_the_exact_finite_grid_reversible_generator():
    _, _, r = setup()
    assert_allclose((r['projection']@r['lift']).toarray(), np.eye(48), atol=3e-16)
    assert_allclose((r['projection']@r['charge']).toarray(), 0, atol=1e-13)
    assert_allclose(r['coarse']@r['coarse_equilibrium'], 0, atol=2e-14)
    assert_allclose(r['continuum']@r['coarse_equilibrium'], 0, atol=2e-14)
    # Both share Gibbs equilibrium; their finite-grid rates need not coincide.
    assert np.max(abs((r['coarse']-r['continuum']).data)) > 1e-3


def test_fast_charge_limit_approaches_projected_dynamics_without_changing_mobility():
    errors = []
    for rate in [30., 300., 3000.]:
        _, p, r = setup(cells=24, rate=rate)
        actual = r['projection']@reversible_propagator(r['generator'], r['equilibrium'], .3)@(r['lift']@p)
        limit = reversible_propagator(r['coarse'], r['coarse_equilibrium'], .3)@p
        errors.append(sum(abs(actual-limit)))
    assert errors[1] < errors[0]/5
    assert errors[2] < errors[1]/5


def test_continuum_grand_potential_limit_requires_spatial_refinement_too():
    errors = []
    for cells in [24, 48, 96]:
        _, p, r = setup(cells=cells)
        exact_grid = reversible_propagator(r['coarse'], r['coarse_equilibrium'], .3)@p
        continuum_grid = reversible_propagator(r['continuum'], r['coarse_equilibrium'], .3)@p
        errors.append(sum(abs(exact_grid-continuum_grid)))
    assert errors[1] < errors[0]/3
    assert errors[2] < errors[1]/3


def test_quenched_mixture_is_not_the_semigroup_of_the_average_generator():
    cells, kt = 24, .12
    q = -1+2/cells*(np.arange(cells)+.5)
    potentials = [.5*(q+.4)**2, .5*(q-.4)**2]
    generators = [sg_chain(f, thermal_energy=kt, spacing=2/cells, face_mobility=1.).toarray()
                  for f in potentials]
    def mix(t):
        return .5*expm(t*generators[0])+.5*expm(t*generators[1])
    p = np.exp(-(q/.2)**2); p /= sum(p)
    one = mix(.6)@p
    reset_hidden_state = mix(.3)@mix(.3)@p
    averaged = expm(.3*(generators[0]+generators[1]))@p
    assert sum(abs(one-reset_hidden_state)) > .01
    assert sum(abs(one-averaged)) > .02
    assert_allclose([sum(one), sum(reset_hidden_state), sum(averaged)], 1, atol=1e-12)


def test_zero_charge_rate_keeps_sectors_quenched_and_one_sector_is_original_sg():
    _, _, r = setup(rate=0.)
    assert r['charge'].nnz == 0 or np.max(abs(r['charge'].data)) == 0
    f = np.linspace(-.1, .2, 10)
    single = coupled_charge_generator(f[None], excess_counts=[1], chemical_potential=.05,
               thermal_energy=.12, spacing=.2, mobilities=[.3], charge_rate=10.)
    expected = sg_chain(f, thermal_energy=.12, spacing=.2, face_mobility=.3)
    assert_allclose(single['generator'].toarray(), expected.toarray(), atol=1e-14)


@pytest.mark.parametrize('counts,rate', [([0, .2], 1), ([1, 1], 1), ([0, 1], -1)])
def test_fractional_duplicate_sectors_and_negative_rates_are_rejected(counts, rate):
    with pytest.raises(ValueError):
        coupled_charge_generator(np.zeros((2, 8)), excess_counts=counts, chemical_potential=0,
            thermal_energy=.1, spacing=.2, mobilities=[1, 1], charge_rate=rate)


def test_sector_dependent_absorption_conserves_probability_and_has_weighted_fast_limit():
    errors=[]
    for rate in [30.,300.,3000.]:
        q,p,r=setup(cells=24,rate=rate)
        sink=np.vstack([.5*(q>.25),2.*(q>.25)])
        average=np.sum(sink*r['weights'],axis=0)
        joint=np.zeros((49,49));joint[:48,:48]=r['generator'].toarray()-np.diag(sink.ravel())
        joint[-1,:48]=sink.ravel()
        coarse=np.zeros((25,25));coarse[:24,:24]=r['coarse'].toarray()-np.diag(average)
        coarse[-1,:24]=average
        a=expm(.3*joint)@np.r_[r['lift']@p,0.]
        b=expm(.3*coarse)@np.r_[p,0.]
        assert a.min()>=0
        assert_allclose(a.sum(),1,atol=1e-12,rtol=0)
        errors.append(abs(a[-1]-b[-1]))
    assert errors[1]<errors[0]/5 and errors[2]<errors[1]/5


@pytest.mark.parametrize('sectors',[2,3])
def test_exact_hidden_charge_coordinates_and_memory_resolvent(sectors):
    from scipy import sparse
    q=np.linspace(-1,1,12)
    f=np.array([.5*(q+.2*j)**2+.01*j for j in range(sectors)])
    r=coupled_charge_generator(f,excess_counts=np.arange(sectors),chemical_potential=.02,
        thermal_energy=.12,spacing=.1,mobilities=np.linspace(.3,1,sectors),charge_rate=3.)
    m=charge_memory_blocks(r)
    transform=sparse.hstack([r['lift'],m['hidden_lift']]).toarray()
    inverse=sparse.vstack([r['projection'],m['hidden_projection']]).toarray()
    assert_allclose(inverse@transform,np.eye(12*sectors),atol=5e-16)
    blocks=sparse.bmat([[m['A'],m['B']],[m['C'],m['D']]]).toarray()
    assert_allclose(transform@blocks@inverse,r['generator'].toarray(),atol=1e-13)
    z=.4+1.7j
    a,b,c,d=[m[k].toarray() for k in ['A','B','C','D']]
    memory=b@np.linalg.solve(z*np.eye(len(d))-d,c)
    reduced=np.linalg.inv(z*np.eye(len(a))-a-memory)
    full=r['projection']@np.linalg.solve(z*np.eye(len(blocks))-r['generator'].toarray(),r['lift'].toarray())
    assert_allclose(reduced,full,atol=2e-14,rtol=2e-13)
    assert np.max(abs(memory))>1e-3


@pytest.mark.parametrize('scale',[1e-15,1.,1e15])
def test_reversible_validation_is_invariant_under_time_units(scale):
    generator=np.array([[-2.,1.],[2.,-1.]])*scale
    pi=np.array([1/3,2/3])
    expected=expm(generator*(.7/scale))
    assert_allclose(reversible_propagator(generator,pi,.7/scale),expected,atol=1e-14)
    leaking=generator.copy();leaking[1,1]-=.1*scale
    with pytest.raises(ValueError,match='conserve'):
        reversible_propagator(leaking,pi,.7/scale)
    with pytest.raises(ValueError,match='stationary'):
        reversible_propagator(generator,np.array([.5,.5]),.7/scale)


@pytest.mark.parametrize('scale',[1e-15,1.,1e15])
def test_stationary_irreversible_cycle_is_rejected_at_every_time_scale(scale):
    cycle=scale*np.array([[-1.,0.,1.],[1.,-1.,0.],[0.,1.,-1.]])
    with pytest.raises(ValueError,match='reversible'):
        reversible_propagator(cycle,np.ones(3)/3,.7/scale)


def test_fast_sector_cannot_hide_invalid_slow_sector():
    from scipy.linalg import block_diag
    slow=1e-15*np.array([[-1.,1.],[1.,-1.1]])
    fast=np.array([[-1.,1.],[1.,-1.]])
    with pytest.raises(ValueError,match='conserve'):
        reversible_propagator(block_diag(slow,fast),np.ones(4)/4,.5)


def test_zero_generator_and_unnormalized_equilibrium():
    assert_allclose(reversible_propagator(np.zeros((2,2)),np.array([.2,.8]),1e15),np.eye(2))
    with pytest.raises(ValueError,match='normalized'):
        reversible_propagator(np.zeros((2,2)),np.array([.2,.8000001]),1.)
