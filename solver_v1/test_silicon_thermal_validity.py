import numpy as np
import pytest
from scipy.constants import h, electron_volt
from scipy.integrate import quad

from .silicon_thermal_validity import KB, oscillator_statistics, validate_fixed_cartesian


def test_classical_variance_against_direct_gaussian_integral():
    k, temperature = .31, 300.
    density = lambda x: np.exp(-.5*k*x*x/(KB*temperature))
    variance = quad(lambda x:x*x*density(x),-np.inf,np.inf)[0]/quad(density,-np.inf,np.inf)[0]
    assert oscillator_statistics([k],temperature)['classical_variance_A2'][0] == pytest.approx(variance, rel=1e-11)


def test_quantum_against_finite_oscillator_partition_sum():
    stats = oscillator_statistics([1.7],300.)
    energy = stats['quantum_energy_eV'][0]
    levels = np.arange(500)+.5
    weights = np.exp(-levels*energy/(KB*300)); weights /= weights.sum()
    # <k x^2> = total oscillator energy by the virial theorem.
    variance = np.dot(weights, levels*energy)/1.7
    assert stats['quantum_variance_A2'][0] == pytest.approx(variance,rel=1e-12)
    assert energy == pytest.approx(h/electron_volt*1e12*stats['omega_rad_ps'][0]/(2*np.pi))


def test_classical_high_temperature_and_quantum_low_temperature_limits():
    high=oscillator_statistics([.1,1,10],1e8)
    np.testing.assert_allclose(high['variance_ratio'],1,rtol=1e-9)
    low=oscillator_statistics([1.7],1e-3)
    assert low['quantum_variance_A2'][0] == pytest.approx(low['quantum_energy_eV'][0]/(2*1.7))


def test_rotation_covariance_and_fixed_grips():
    free=np.array([True,False]);basis=np.vstack([np.eye(3),np.zeros((3,3))])
    hessian=np.array([[2.,.5,.2],[.5,1.,-.1],[.2,-.1,3.]])
    q,_=np.linalg.qr(np.array([[1.,2.,4.],[3.,1.,5.],[2.,8.,1.]]))
    h1,_=validate_fixed_cartesian(hessian,basis,free)
    h2,_=validate_fixed_cartesian(q.T@hessian@q,basis@q,free)
    original=basis@np.linalg.inv(h1)@basis.T
    rotated=basis@q@np.linalg.inv(h2)@q.T@basis.T
    np.testing.assert_allclose(original,rotated,atol=1e-14)
    assert np.max(abs(original[3:])) == 0


@pytest.mark.parametrize('values,temperature,mass', [([0],300,28),([-1],300,28),([np.nan],300,28),([1],0,28),([1],300,0)])
def test_reject_nonphysical_oscillators(values,temperature,mass):
    with pytest.raises(ValueError): oscillator_statistics(values,temperature,mass)


def test_reject_moving_grip_basis():
    basis=np.vstack([np.eye(3),np.zeros((3,3))]);basis[3,0]=.1
    with pytest.raises(ValueError): validate_fixed_cartesian(np.eye(3),basis,np.array([True,False]))
