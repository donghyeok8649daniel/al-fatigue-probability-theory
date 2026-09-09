import numpy as np
import pytest

from .activation_stress_sensitivity import barrier_and_derivative, apparent_rate_slope_from_barrier
from .interface_static_scenarios import InterfaceUnits
from .vector_interface_reference import VectorInterfaceEvaluation


class AnalyticTestLandscape:
    h=1.
    def evaluate(self,q):
        a,x,y=q; angle=2*np.pi*x
        return VectorInterfaceEvaluation(float((a-1)**2+1-np.cos(angle)+y*y),
            np.array([2*(a-1),2*np.pi*np.sin(angle),2*y]),
            np.diag([2.,(2*np.pi)**2*np.cos(angle),2.]),{}, {})


def test_stationary_envelope_derivative_uses_atomic_cell_not_specimen_count():
    model=AnalyticTestLandscape();units=InterfaceUnits(3e-10,8e-20)
    value=barrier_and_derivative(model,[1.,0.,0.],[1.,.5,0.],
        traction_MPa=[2.,4.,1.],units=units)
    for axis in range(3):
        energies=[]
        for sign in (-1,1):
            load=np.array([2.,4.,1.]);load[axis]+=sign*.01
            energies.append(barrier_and_derivative(model,value['minimum']['q'],
                value['saddle']['q'],traction_MPa=load,units=units)['barrier_eV_cell'])
        derivative=-(energies[1]-energies[0])/.02
        assert derivative==pytest.approx(value['minus_barrier_traction_derivative_eV_per_MPa'][axis],abs=3e-13)
    assert value['physical_local_activation_event_validated'] is False
    assert value['cell_mobility'] is None


def test_barrier_only_slope_retains_sign_and_does_not_provide_rate():
    got=apparent_rate_slope_from_barrier(np.array([-1,0,1])*1e-29,293.)
    assert got[0]==-got[2] and got[1]==0
    assert got[2]==pytest.approx(1e-23/(1.380649e-23*293))
    with pytest.raises(ValueError):apparent_rate_slope_from_barrier(1.,0.)
