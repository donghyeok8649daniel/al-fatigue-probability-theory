import numpy as np

from solver_v1.silicon_crack_research import AtomicEvaluation, RelaxedCoordinates
from solver_v1.silicon_conditional_research import AMU_EV_PS2_A2
from solver_v1.silicon_thermal_dynamics import (
    reflecting_drift, constrained_verlet, stationary_correlation,
)


def test_box_drift_multiple_reflections_and_time_reversal():
    x,v=np.array([.1,.7,-.6]),np.array([1.,-3.,8.])
    moved,speed,count=reflecting_drift(x,v,1.7,1.)
    assert count>3
    np.testing.assert_allclose(abs(speed),abs(v),atol=0)
    restored,reverse,_=reflecting_drift(moved,-speed,1.7,1.)
    np.testing.assert_allclose(restored,x,atol=2e-15)
    np.testing.assert_allclose(reverse,-v,atol=0)


class HarmonicAtoms:
    def evaluate(self,r):
        return AtomicEvaluation(float(np.sum(r*r)),2*r,np.sum(r*r,axis=1))


def test_constrained_verlet_preserves_gap_grips_energy_order_and_reversal():
    reference=np.array([[0.,0.,0.],[0.,-.2,0.],[0.,.2,0.]])
    coords=RelaxedCoordinates(reference,[True,False,False],bond=[1,2])
    initial=np.array([.1,.0,-.03,.07,.02])
    velocity=np.array([.03,.02,-.04,.01,.05])
    errors=[]
    for dt in [.04,.02,.01]:
        run=constrained_verlet(HarmonicAtoms(),coords,initial,velocity,gap=.4,
            dt_ps=dt,steps=round(2/dt),halfwidth=2.,mass_amu=1/AMU_EV_PS2_A2)
        errors.append(run['max_energy_residual_eV'])
        positions=run['final_positions_A']
        np.testing.assert_allclose(positions[0],reference[0],atol=0)
        np.testing.assert_allclose(positions[2,1]-positions[1,1],.4,atol=2e-15)
        assert run['reflection_count']==0
    assert 3.8<errors[0]/errors[1]<4.2
    assert 3.8<errors[1]/errors[2]<4.2
    reverse=constrained_verlet(HarmonicAtoms(),coords,run['final_bath_coordinates_A'],
        -run['final_bath_velocities_A_ps'],gap=.4,dt_ps=.01,steps=200,halfwidth=2.,
        mass_amu=1/AMU_EV_PS2_A2)
    np.testing.assert_allclose(reverse['final_bath_coordinates_A'],initial,atol=2e-14)
    np.testing.assert_allclose(reverse['final_bath_velocities_A_ps'],-velocity,atol=2e-14)


def test_fft_correlation_equals_direct_time_origins_with_external_mean():
    x=np.random.default_rng(119).standard_normal(257)+.7
    for center in [None,.6]:
        mean=x.mean() if center is None else center
        actual=stationary_correlation(x,100,center=center)
        expected=np.array([np.mean((x[:len(x)-lag]-mean)*(x[lag:]-mean)) for lag in range(101)])
        np.testing.assert_allclose(actual,expected,atol=7e-16)
