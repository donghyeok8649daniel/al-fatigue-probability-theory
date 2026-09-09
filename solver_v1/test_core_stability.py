"""Synthetic quadratic unit tests, not material/MD calibration data."""
import numpy as np
import pytest

from .core_stability import lowest_core_mode


class QuadraticTest:
    def __init__(self,diagonal):
        self.diagonal=np.array(diagonal).reshape((1,3))

    def evaluate(self,q):
        return dict(energy=float(np.sum(self.diagonal*q*q)/2),gradient=self.diagonal*q)

    def linearize(self,q):
        return self.evaluate(q),lambda v:self.diagonal*v


def test_zero_force_saddle_is_rejected_by_independent_energy_curvature():
    check=lowest_core_mode(QuadraticTest([-2.,3.,4.]),np.zeros((1,3)))
    assert check['maximum_free_force']==0
    assert check['minimum_eigenvalue']==pytest.approx(-2.)
    assert not check['stable_on_tested_fixed_boundary']
    assert check['eigenpair_residual']<1e-12
    assert check['eigenvector'][0,0]>0
    for fd in check['energy_curvature_checks']:
        assert fd['curvature']==pytest.approx(-2.)
        assert fd['plus_energy_change']<0 and fd['minus_energy_change']<0


def test_positive_curvature_still_requires_stationary_force():
    core=QuadraticTest([2.,3.,4.])
    assert lowest_core_mode(core,np.zeros((1,3)))['stable_on_tested_fixed_boundary']
    out=lowest_core_mode(core,np.ones((1,3))*.01)
    assert out['positive_hessian_on_tested_fixed_boundary']
    assert not out['stable_on_tested_fixed_boundary']
