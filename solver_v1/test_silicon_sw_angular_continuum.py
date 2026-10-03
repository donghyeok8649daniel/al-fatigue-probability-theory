"""Polynomial root/normalization controls used by continuous angular fitting."""
import numpy as np
from .run_silicon_sw_angular_continuum import real_roots,squared_polynomial


def test_quartic_squared_residual_and_derivative_roots():
    # R(beta)=[beta^2-.04, beta-.2]; its exact minimum is beta=.2.
    vectors=np.array([[-.04,-.2],[0.,1.],[1.,0.]])
    loss=squared_polynomial(vectors)
    roots=real_roots(np.polynomial.polynomial.polyder(loss),-.375,.375)
    assert any(abs(r-.2)<1e-12 for r in roots)
    for beta in [-.375,-.07,0.,.13,.375]:
        residual=vectors[0]+beta*vectors[1]+beta**2*vectors[2]
        np.testing.assert_allclose(np.polynomial.polynomial.polyval(beta,loss),residual@residual,atol=2e-16)


def test_guard_roots_and_eos_mean_preserve_physical_scaling():
    roots=real_roots([-.04,0.,1.,0.,0.],-.375,.375)
    np.testing.assert_allclose(sorted(roots),[-.2,.2],atol=1e-14)
    assert real_roots([1.],-.375,.375)==[]
    assert real_roots([.04,0.,1.],-.375,.375)==[]
    eos=np.arange(12,dtype=float)/1000
    normalized=eos/(.02*np.sqrt(len(eos)))
    np.testing.assert_allclose(normalized@normalized,np.mean((eos/.02)**2),atol=3e-17)
