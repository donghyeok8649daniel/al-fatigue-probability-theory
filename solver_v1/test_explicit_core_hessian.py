"""Dense and matrix-free analytic core Hessians must be the SAME operator."""
import numpy as np
import pytest

from .core_stability import lowest_core_mode
from .current_material_rows import CurrentMaterialScrewCore
from .isolated_screw_core import IsolatedScrewCore
from .test_current_material_rows import model, reference_far


@pytest.mark.parametrize('current,radius,ring', [(True,.85,1),(True,1.3,2),(False,.85,1),(False,1.3,2)])
def test_explicit_matrix_equals_full_analytic_action(model,current,radius,ring):
    cls, surface = (CurrentMaterialScrewCore,model) if current else (IsolatedScrewCore,model.base.surface)
    core = cls(surface,reference_far(model),free_radius=radius,ring=ring)
    q = core.initial+.002*np.sin(np.arange(core.initial.size)*.4).reshape(core.initial.shape)
    _, H = core.linearize(q)
    matrix = H.explicit_matrix()
    assert np.max(abs(matrix-matrix.T)) < 2e-11
    for vector in (np.sin(np.arange(q.size)*.73),np.cos(np.arange(q.size)*1.13)):
        np.testing.assert_allclose(matrix@vector,H(vector.reshape(q.shape)).ravel(),atol=3e-10,rtol=3e-11)
    a = lowest_core_mode(core,q,check_energy=False,method='operator')
    b = lowest_core_mode(core,q,check_energy=False,method='assembled')
    assert a['minimum_eigenvalue'] == pytest.approx(b['minimum_eigenvalue'],abs=2e-10,rel=2e-10)
    assert b['eigenpair_residual'] < 1e-9


def test_no_silent_hessian_method_replacement(model):
    core = CurrentMaterialScrewCore(model,reference_far(model),free_radius=.85,ring=1)
    with pytest.raises(ValueError,match='Morse'):
        lowest_core_mode(core,core.initial,method='guess')
