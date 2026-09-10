"""Numerical stationary roots distinguish minima/saddles from dynamics."""
import numpy as np
import pytest

from .core_stationary_point import stationary_core


class QuarticTestLandscape:
    def evaluate(self,q):
        x,y = np.asarray(q).ravel()
        return dict(energy=(x*x-1)**2/4+y*y,
                    gradient=np.array([[x*(x*x-1),2*y]]))

    def linearize(self,q):
        x,y = np.asarray(q).ravel()
        matrix = np.diag([3*x*x-1,2.])
        def action(v):
            return (matrix@np.asarray(v).ravel()).reshape(1,2)
        action.explicit_matrix=lambda:matrix
        return self.evaluate(q),action


def test_stationary_saddle_root_not_a_minimization_or_time_law():
    core = QuarticTestLandscape()
    result = stationary_core(core,np.array([[.08,.12]]),force_tolerance=1e-11)
    assert result['force_converged'] and result['negative_eigenvalues']==1
    assert result['unresolved_eigenvalues']==0
    np.testing.assert_allclose(result['field'],0.,atol=1e-11)
    assert result['energy']==pytest.approx(.25)
    assert result['physical_dynamics'] is False


def test_stationary_minimum_and_unconverged_budget_are_separate():
    core = QuarticTestLandscape()
    result = stationary_core(core,np.array([[.8,.2]]),force_tolerance=1e-11)
    assert result['force_converged'] and result['negative_eigenvalues']==0
    assert result['energy']<1e-20
    incomplete = stationary_core(core,np.array([[.8,.2]]),force_tolerance=1e-14,max_iterations=1,max_step=.001)
    assert not incomplete['force_converged']
