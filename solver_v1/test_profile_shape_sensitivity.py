import numpy as np
import pytest
from .core_interface_compatibility import minimax_compatibility
from .profile_shape_sensitivity import envelope_gradient,difference_stencil
from .vector_material_calibration import MaterialObservation


@pytest.mark.parametrize('target,theta',[(0.,1.),(2.,1.),(2.,3.)])
def test_dual_derivative_matches_independently_resolved_profile(target,theta):
    obs=[MaterialObservation('exact',1.,.7,'test','exact'),
         MaterialObservation('fit',target,3.,'test','fit')]
    def evaluate(x):
        M=np.array([[np.exp(.2*x)],[x*np.exp(.2*x)]])
        return M,minimax_compatibility(M,obs,[1],nonnegative=(0,))
    M,p=evaluate(theta)
    dM=np.exp(.2*theta)*np.array([[.2],[1+.2*theta]])
    value=envelope_gradient(dM[None],obs,p)[0]
    h=1e-5
    numerical=(evaluate(theta+h)[1]['minimax_normalized_error']-evaluate(theta-h)[1]['minimax_normalized_error'])/(2*h)
    assert value==pytest.approx(numerical,abs=2e-10)
    assert value==pytest.approx(np.sign(theta-target)/3,abs=1e-12)


@pytest.mark.parametrize('x',[-1.,0.,1.])
def test_bounded_stencil_and_column_gauge(x):
    stencil=difference_stencil([x],0,[-1.],[1.],1e-4)
    assert all(-1<=v[0]<=1 for v,w in stencil)
    assert sum(w*v[0]**2 for v,w in stencil)==pytest.approx(2*x,abs=2e-10)
    obs=[MaterialObservation('exact',1.,1.,'test','exact'),MaterialObservation('fit',0.,1.,'test','fit')]
    M=np.exp(x)*np.array([[1.],[2.]])
    profile=minimax_compatibility(M,obs,[1],nonnegative=(0,))
    assert abs(envelope_gradient(M[None],obs,profile)[0])<1e-12


def test_reject_uncertified_or_malformed_sensitivity():
    with pytest.raises(ValueError):envelope_gradient(np.zeros((1,1,1)),[],{'completed':False})
    with pytest.raises(ValueError):difference_stencil([2.],0,[-1.],[1.],1e-4)
