"""Source-only core checks; the hybrid LJ/Bessel model is never replaced."""
from types import SimpleNamespace

import numpy as np
import pytest

from .isolated_screw_core import ScrewFarField
from .reference_eam_targets import MishinRigidFCCReference
from .source_core_reference import SourceRowKernel, MishinScrewCoreReference, source_elastic_tensor
from .vector_material_calibration import LENGTH_M


@pytest.fixture(scope='module')
def setup():
    source = MishinRigidFCCReference()
    tensor = source_elastic_tensor(source, LENGTH_M)
    b, h = source.geometry.b/(LENGTH_M/1e-10), source.h/(LENGTH_M/1e-10)
    return source, ScrewFarField(tensor, b, (np.sqrt(3)*b/12, h/2))


def test_source_row_units_and_derivatives(setup):
    source, far = setup
    rows = SimpleNamespace(source=source, length_angstrom=LENGTH_M/1e-10, b=far.b)
    kernel = SourceRowKernel(rows)
    points = np.array([[.12, .88, .13], [-.24, .33, 1.15]])
    exact = kernel.evaluate(points, order=2)
    for axis in range(3):
        step = np.eye(3)[axis]*1e-6
        plus, minus = kernel.evaluate(points+step), kernel.evaluate(points-step)
        np.testing.assert_allclose((plus['value']-minus['value'])/2e-6,
            exact['gradient'][..., axis], atol=3e-8, rtol=3e-7)
        np.testing.assert_allclose((plus['gradient']-minus['gradient'])/2e-6,
            exact['hessian'][..., axis], atol=5e-6, rtol=3e-6)
    assert exact['source_cutoff_exact'] and not exact['reciprocal_Bessel_model']


def test_source_core_counts_forces_hessian_and_bulk(setup):
    source, far = setup
    core = MishinScrewCoreReference(source, far, length_scale_m=LENGTH_M,
        free_radius=.85, ring=4, burgers_sign=0)
    zero = core.evaluate(core.initial)
    assert abs(zero['energy']) < 1e-14 and np.max(abs(zero['gradient'])) < 1e-13
    np.testing.assert_allclose(core.affine_antiplane_hessian(), far.matrix, atol=1e-11, rtol=1e-11)
    q = core.initial+.002*np.cos(np.arange(core.initial.size)*.71).reshape(core.initial.shape)
    out, H = core.linearize(q)
    v = np.sin(np.arange(q.size)+.4).reshape(q.shape); v /= np.linalg.norm(v)
    np.testing.assert_allclose(H.explicit_matrix()@v.ravel(),H(v).ravel(),atol=1e-12,rtol=1e-12)
    step = 1e-6
    plus, minus = core.evaluate(q+step*v), core.evaluate(q-step*v)
    assert abs((plus['energy']-minus['energy'])/(2*step)-np.sum(out['gradient']*v)) < 1e-8
    np.testing.assert_allclose((plus['gradient']-minus['gradient'])/(2*step), H(v), atol=5e-6, rtol=3e-6)
    assert np.max(abs(out['all_site_gradient'].sum(axis=0))) < 1e-12
    assert core.reference_only and not core.production_eligible


def test_source_cutoff_guard_does_not_silently_omit_neighbors(setup):
    source, far = setup
    core = MishinScrewCoreReference(source, far, length_scale_m=LENGTH_M,
        free_radius=.85, ring=1, burgers_sign=0)
    with pytest.raises(ValueError, match='enclose'):
        core.evaluate(core.initial)


def test_source_periodic_row_shift_is_not_a_strength_factor(setup):
    source, far = setup
    core = MishinScrewCoreReference(source, far, length_scale_m=LENGTH_M,
        free_radius=.85, ring=4)
    q = core.initial
    first = core.evaluate(q)
    q[0, 0] += core.rows.b
    second = core.evaluate(q)
    assert first['energy'] == pytest.approx(second['energy'], abs=1e-12)
    np.testing.assert_allclose(first['gradient'], second['gradient'], atol=3e-12, rtol=3e-12)
