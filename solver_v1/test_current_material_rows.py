"""Representation bridge checks, not material/kinetic/yield adoption."""
from itertools import product
from pathlib import Path
import json

import numpy as np
import pytest

from .coordination_screening import CoordinationScreenedInterface
from .current_material_rows import CurrentMaterialRowKernel, CurrentMaterialSiteLaw, CurrentMaterialScrewCore
from .fcc111_geometry import fcc111_geometry_from_b
from .isolated_screw_core import IsolatedScrewCore, ScrewFarField, row_environment
from .nonlinear_fcc_screw import stf_basis
from .nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor
from .quadrupole_saturation import CubicInterfaceMomentSites
from .vector_interface_reference import HESSIAN_INDICES
from .vector_fcc_validation import direct_row_channels


ROOT = Path(__file__).resolve().parents[1]
V20 = ROOT/'results/fcc111_active_interface/tangent_calibration_v20/shape_refinement/old_family/calibration.json'


@pytest.fixture(scope='module')
def model():
    saved = json.loads(V20.read_bytes())
    return CoordinationScreenedInterface(saved['best']['shape'], saved['best']['coefficients'], law=saved['screening_law'])


def reference_far(model):
    plane = fcc111_geometry_from_b(1.).plane_basis_in_stacked_cubic_axes()
    # Proper eV/L0^3 units, fixed independent cubic elastic combinations.
    length = 4.05/np.sqrt(2)*1e-10
    C = rotate_elastic_tensor(cubic_elastic_tensor(114., 62., 32.), plane)*1e9*length**3/1.602176634e-19
    return ScrewFarField(C, 1., (np.sqrt(3)/12, model.h/2))


def site_interface_channels(model, state):
    density = model.base._density(np.asarray(state))[0]
    moments = {m.invariant.rank: m for _, m in model.base.moments}
    jets = [CubicInterfaceMomentSites(model.base, moments[r]).site_jets(tuple(state))[0]
            for r in (1, 2, 3)]
    jets.append(CubicInterfaceMomentSites(model.base, model.moment.parts[0]).site_jets(tuple(state))[0])
    count = max(len(density), *(len(j) for j in jets))
    out = np.zeros((count, 10, 22))
    out[:len(density), :, 1] = density
    for start, jet in zip((2, 5, 10, 17), jets):
        out[:len(jet), :, start:start+jet.shape[-1]] = jet
    return out


def test_same_interface_energy_gradient_hessian_all_site_terms(model):
    law = CurrentMaterialSiteLaw(model)
    for state in ([model.h, 0, 0], [1.08*model.h, .19, -.037], [1.47*model.h, .41, .09]):
        channels = site_interface_channels(model, state)
        E, g, H = law.evaluate(channels[:, 0])
        first = channels[:, 1:4]
        second = channels[:, HESSIAN_INDICES]
        grad = np.einsum('nc,nic->i', g, first)
        hess = np.einsum('ncd,nic,njd->ij', H, first, first)+np.einsum('nc,nijc->ij', g, second)
        jet = 2*np.r_[E.sum(), grad, hess[0], hess[1, 1:], hess[2, 2]]
        old = model.evaluate(state)
        expected = sum(v for k, v in old.components.items() if not k.startswith('pair_'))
        np.testing.assert_allclose(jet, expected, atol=3e-10, rtol=2e-10)


def test_site_law_fd_and_nonzero_density_moment_cross_hessian(model):
    # Exercise the existing power screening law away from its zero extension.
    shape = model.screened_shape.copy(); shape[5] = -.3
    screened = CoordinationScreenedInterface(shape, model.coefficients, law='power')
    law = CurrentMaterialSiteLaw(screened)
    z = .001*np.sin(np.arange(44)).reshape(2, 22)
    z[:, 1] = -.13*law.rho_ref
    E, g, H = law.evaluate(z)
    assert np.max(abs(H[:, 1, 2:5])) > 0
    assert np.max(abs(H[:, 5:10, 17:22])) > 0
    np.testing.assert_allclose(H, H.swapaxes(-1, -2), atol=1e-13, rtol=1e-13)
    v = np.cos(np.arange(44)*.71).reshape(z.shape); v[:, 1] *= law.rho_ref
    step = 2e-7
    plus = law.evaluate(z+step*v); minus = law.evaluate(z-step*v)
    np.testing.assert_allclose((plus[0]-minus[0])/(2*step), np.sum(g*v, axis=1), atol=2e-8, rtol=2e-7)
    errors = []
    exact = np.einsum('nij,nj->ni', H, v)
    # The fixed Eg amplitude gauge makes third derivatives large. Verify
    # actual central-difference refinement, not a looser material tolerance.
    for step in (2e-7, 1e-7, 5e-8):
        plus = law.evaluate(z+step*v); minus = law.evaluate(z-step*v)
        fd = (plus[1]-minus[1])/(2*step)
        errors.append(float(np.max(abs(fd-exact))))
    assert errors[1] < .35*errors[0] and errors[2] < .35*errors[1]
    np.testing.assert_allclose(fd, exact, atol=3e-6, rtol=3e-7)


def test_extra_quadrupole_infinite_row_vs_independent_real_atoms(model):
    rows = row_environment(model.base.surface)
    kernel = CurrentMaterialRowKernel(rows)
    vectors = np.array([[.17, .69, .31], [-.22, -.57, .83], [.37, 2.17, -.67]])
    values = kernel.evaluate(vectors, order=2)
    inv = rows.surface.angular
    expected = []
    for x, y, z in vectors:
        atoms = np.column_stack([np.arange(-192, 193)*rows.b+x, np.full(385, y), np.full(385, z)])
        weight = inv.amplitude*np.exp(-inv.kappa*np.linalg.norm(atoms, axis=1))
        raw = np.array([np.sum(weight*np.prod(atoms[:, axes], axis=1)) for axes in product(range(3), repeat=2)])
        expected.append(raw@stf_basis(2))
    np.testing.assert_allclose(values['value'][:, 17:], expected, atol=2e-13, rtol=2e-11)
    for axis in range(3):
        d = np.eye(3)[axis]*2e-6
        plus, minus = kernel.evaluate(vectors+d), kernel.evaluate(vectors-d)
        np.testing.assert_allclose((plus['gradient']-minus['gradient'])/4e-6,
                                   values['hessian'][..., axis], atol=2e-6, rtol=3e-6)
    reverse = kernel.evaluate(-vectors, order=2)
    for field, order in [('value', 0), ('gradient', 1), ('hessian', 2)]:
        parity = kernel.parity.reshape((1, 22)+(1,)*order)*(-1)**order
        np.testing.assert_allclose(reverse[field], values[field]*parity, atol=2e-12, rtol=2e-12)


def test_all_current_channels_independent_real_atom_derivatives_and_tail(model):
    rows = row_environment(model.base.surface)
    kernel = CurrentMaterialRowKernel(rows)
    vectors = np.array([[.13, .74, -.43], [-.31, .26, 1.12], [.48, 1.81, 2.12]])
    analytic = kernel.evaluate(vectors, order=2)
    for i, vector in enumerate(vectors):
        coarse = direct_row_channels(rows, vector, images=32, include_odd_quadrupole=True)
        fine = direct_row_channels(rows, vector, images=128, include_odd_quadrupole=True)
        for key in ('value', 'gradient', 'hessian'):
            np.testing.assert_allclose(fine[key], analytic[key][i], atol=2e-11, rtol=3e-10)
        # Uncorrected finite LJ sums have a nonzero algebraic tail. The Bessel
        # reference is infinite; this real-space evaluator is validation only.
        error32 = abs(coarse['value'][0]-analytic['value'][i, 0])
        error128 = abs(fine['value'][0]-analytic['value'][i, 0])
        assert error128 < error32/100


def test_same_material_bulk_antiplane_identity(model):
    core = CurrentMaterialScrewCore(model, reference_far(model), free_radius=.85, ring=4)
    np.testing.assert_allclose(core.affine_antiplane_hessian(), core.far_field.matrix,
                               atol=7e-9, rtol=3e-9)


def test_full_core_chain_rule_symmetry_and_exact_zero_state(model):
    core = CurrentMaterialScrewCore(model, reference_far(model), free_radius=.85, ring=1, burgers_sign=0)
    perfect = core.evaluate(core.initial)
    assert abs(perfect['energy']) < 1e-14
    assert np.max(abs(perfect['gradient'])) < 1e-12
    q = core.initial+.002*np.sin(np.arange(core.initial.size)).reshape(core.initial.shape)
    out, H = core.linearize(q)
    v = np.cos(np.arange(q.size)).reshape(q.shape); v /= np.linalg.norm(v)
    w = np.sin(np.arange(q.size)*.6+.2).reshape(q.shape)
    assert abs(np.sum(w*H(v))-np.sum(v*H(w))) < 2e-10
    step = 2e-6
    plus, minus = core.evaluate(q+step*v), core.evaluate(q-step*v)
    assert abs((plus['energy']-minus['energy'])/(2*step)-np.sum(out['gradient']*v)) < 1e-7
    np.testing.assert_allclose((plus['gradient']-minus['gradient'])/(2*step), H(v), atol=2e-6, rtol=3e-6)
    assert np.max(abs(out['all_site_gradient'].sum(axis=0))) < 1e-11
    shifted = q.copy(); shifted[0, 0] += core.rows.b
    after = core.evaluate(shifted)
    assert after['energy'] == pytest.approx(out['energy'], abs=3e-12)
    np.testing.assert_allclose(after['gradient'], out['gradient'], atol=2e-11)


def test_nested_zero_extension_recovers_historical_core(model):
    shape = model.screened_shape.copy(); shape[3] = shape[5] = 0.
    coefficients = model.coefficients.copy(); coefficients[8:] = 0.
    simple = CoordinationScreenedInterface(shape, coefficients, law='power')
    far = reference_far(simple)
    old = IsolatedScrewCore(simple.base.surface, far, free_radius=.85, ring=1)
    new = CurrentMaterialScrewCore(simple, far, free_radius=.85, ring=1)
    q = old.initial+.001*np.cos(np.arange(old.initial.size)).reshape(old.initial.shape)
    a, A = old.linearize(q); b, B = new.linearize(q)
    assert a['energy'] == pytest.approx(b['energy'], abs=2e-12)
    np.testing.assert_allclose(a['gradient'], b['gradient'], atol=2e-11)
    np.testing.assert_allclose(A(q), B(q), atol=2e-10)


def test_unsafe_density_and_wrong_model_refused(model):
    law = CurrentMaterialSiteLaw(model)
    z = np.zeros((1, 22)); z[:, 1] = -law.rho_bulk
    with pytest.raises(ValueError, match='density'):
        law.evaluate(z)
    with pytest.raises(TypeError):
        CurrentMaterialSiteLaw(model.base.surface)
