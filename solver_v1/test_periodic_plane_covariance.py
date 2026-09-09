import numpy as np
import pytest

from .periodic_plane_covariance import (
    plane_mode_equipartition,empirical_plane_modes,KB_J_K,SourceEAMBlochHessian,
)


def test_periodic_nearest_plane_chain_has_exact_equipartition_scaling():
    n=12;count=576;T=300.;spring=np.diag([17.,31.,43.])
    weights=4*np.sin(np.pi*np.arange(1,n)/n)**2
    modes,local=plane_mode_equipartition(weights[:,None,None]*spring,
        planes=n,atoms_per_plane=count,temperature_K=T)
    expected=KB_J_K*T/count*np.linalg.inv(spring)
    assert np.allclose(modes,expected,rtol=2e-15,atol=0)
    assert np.allclose(local,(1-1/n)*expected,rtol=2e-15,atol=0)
    doubled=plane_mode_equipartition(weights[:,None,None]*spring,
        planes=n,atoms_per_plane=2*count,temperature_K=T)[1]
    assert np.array_equal(doubled,local/2)
    bad=weights[:,None,None]*spring;bad[3]*=-1
    with pytest.raises(ValueError):plane_mode_equipartition(bad,planes=n,atoms_per_plane=count,temperature_K=T)


def test_unitary_dft_covariance_matches_actual_real_space_without_independence_claim():
    t=np.arange(90)[:,None,None];p=np.arange(12)[None,:,None];axes=np.arange(3)[None,None,:]
    q=1e-12*np.sin(.17*t+2*np.pi*p/12+.71*axes)
    C,local=empirical_plane_modes(q)
    assert np.allclose(C.sum(axis=0).real/12,local,rtol=2e-14,atol=1e-39)
    assert np.max(abs(C.sum(axis=0).imag))<1e-39
    assert np.max(abs(C[0]))<1e-50


def test_published_source_bloch_hessian_is_even_and_translation_invariant():
    from .reference_eam_targets import MishinRigidFCCReference
    source=MishinRigidFCCReference(lattice_angstrom=4.065)
    operator=SourceEAMBlochHessian(source)
    assert np.array_equal(operator.evaluate([0,0,0]),np.zeros((3,3)))
    q=np.array([.2,.2,.2]);H=operator.evaluate(q)
    assert np.allclose(H,operator.evaluate(-q),rtol=0,atol=1e-14)
    assert np.linalg.eigvalsh(H)[0]>0


@pytest.mark.parametrize('mode',[1,3,6])
def test_source_bloch_derivative_matches_independent_site_energy(mode):
    from .reference_eam_targets import MishinRigidFCCReference
    source=MishinRigidFCCReference(lattice_angstrom=4.065)
    operator=SourceEAMBlochHessian(source)
    n=12;q=np.ones(3)*2*np.pi*mode/(n*4.065)
    H=operator.evaluate(q)
    for direction in np.eye(3):
        baseline=operator.direct_sinusoidal_energy(planes=n,mode=mode,
            polarization_plane=direction,amplitude_angstrom=0.)
        expected=direction@H@direction
        errors=[]
        for amplitude in (4e-4,2e-4):
            energy=operator.direct_sinusoidal_energy(planes=n,mode=mode,
                polarization_plane=direction,amplitude_angstrom=amplitude)
            mean_cosine2=np.mean(np.cos(2*np.pi*mode*np.arange(n)/n)**2)
            measured=2*(energy-baseline)/(amplitude**2*mean_cosine2)
            errors.append(abs(measured-expected))
        assert errors[-1]<3e-5*expected
        assert errors[-1]<=errors[0]+1e-7
