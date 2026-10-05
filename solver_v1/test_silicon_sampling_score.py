import numpy as np
import pytest

from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_sampling_score import nonlinear_correction_gradient


def target(temperature):
    r = np.array([[0.,0.,0.],[5.,0.,0.]])
    basis = np.zeros((6,3))
    basis[3:] = np.array([[.6,-.8,0.],[.8,.6,0.],[0.,0.,1.]])
    h = np.array([[2.,.3,.1],[.3,3.,-.2],[.1,-.2,1.5]])
    return FixedGripTarget(r,basis,h,[.05,-.07,.02],temperature,r,1.,1.)


@pytest.mark.parametrize('temperature',[80.,300.,600.])
def test_actual_nonlinear_energy_difference_in_rotated_coupled_coordinates(temperature):
    t = target(temperature)
    u = np.array([.3,-.8,.2])
    def energy(v):
        x,y,z = t.coordinates(v)[1]-t.positions[1]
        return .7*x*x+1.3*y*y+.9*z*z+.6*x*y+.2*x**4+.1*y*y*z+.04*x
    r = t.coordinates(u)[1]-t.positions[1]
    x,y,z = r
    forces = np.array([[1e8,-2e8,3e8],
                      [-1.4*x-.6*y-.8*x**3-.04,
                       -2.6*y-.6*x-.2*y*z, -1.8*z-.1*y*y]])
    gradient = nonlinear_correction_gradient(t,u,forces)
    step = 2e-5
    finite = np.array([(t.correction(u+step*e,energy(u+step*e),0.)
                       -t.correction(u-step*e,energy(u-step*e),0.))/(2*step)
                      for e in np.eye(3)])
    np.testing.assert_allclose(gradient,finite,rtol=0,atol=2e-10)
    forces[0] = 0.
    np.testing.assert_array_equal(gradient,nonlinear_correction_gradient(t,u,forces))


def test_exact_reference_energy_leaves_zero_correction_gradient():
    t = target(300.)
    h = t.reference.lower@t.reference.lower.T
    for u in (np.zeros(3),np.array([1.,-.4,.7])):
        d = t.reference.transform(u)
        forces = -(t.basis@(h@(d-t.reference.mean))).reshape(2,3)
        np.testing.assert_allclose(nonlinear_correction_gradient(t,u,forces),0.,atol=5e-16)


def test_score_rejects_invalid_forces_and_coordinates():
    t = target(300.)
    with pytest.raises(ValueError): nonlinear_correction_gradient(t,[1.,2.],np.zeros((2,3)))
    with pytest.raises(ValueError): nonlinear_correction_gradient(t,[1.,2.,3.],np.full((2,3),np.nan))
