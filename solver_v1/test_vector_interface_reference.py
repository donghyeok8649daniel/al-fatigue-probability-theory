"""Independent vector registry identities; no material/kinetic certification."""
import copy

import numpy as np
import pytest

from .fcc111_geometry import DIRECT_110, SHOCKLEY_112, registry_path
from .run_low_stress_cyclic_diagnostic import build_surface
from .reference_eam_targets import DEFAULT_CACHE, MishinRigidFCCReference
from .vector_interface_reference import FullRegistryInterface, MishinVectorInterfaceReference
from .vector_registry_audit import stationary_state, relax_at_x, saddle_downhill_endpoints


@pytest.fixture(scope='module')
def analytic():
    surface, units, metadata = build_surface(tolerance=2e-11)
    return FullRegistryInterface(surface), units, metadata


@pytest.fixture(scope='module')
def source(analytic):
    if not DEFAULT_CACHE.exists():
        pytest.skip('optional checksum-verified NIST source not downloaded')
    return MishinVectorInterfaceReference(MishinRigidFCCReference(), analytic[1].length_scale_m/1e-10)


@pytest.mark.parametrize('path', [DIRECT_110, SHOCKLEY_112])
def test_recovers_verified_scalar_direction(analytic, path):
    model, _, _ = analytic
    surface, _, _ = build_surface(path, tolerance=2e-11)
    direction = registry_path(path).direction(); a = 1.041*model.h; s = .219
    out = model.evaluate([a, *(s*direction)])
    tangent = np.r_[0., direction]
    packed = np.array([out.energy, out.gradient[0], out.gradient@tangent,
                      out.hessian[0, 0], out.hessian[0]@tangent, tangent@out.hessian@tangent])
    np.testing.assert_allclose(packed, surface.packed(a, s), rtol=3e-10, atol=2e-10)


@pytest.mark.parametrize('qfrac', [(1.017, .173, .079), (1.21, -.217, .341)])
def test_full_analytic_gradient_hessian(analytic, qfrac):
    model = analytic[0]; q = np.array(qfrac); q[0] *= model.h
    out = model.evaluate(q); steps = (2e-5, 1e-5)
    errors = []
    for step in steps:
        plus = [model.evaluate(q+step*axis) for axis in np.eye(3)]
        minus = [model.evaluate(q-step*axis) for axis in np.eye(3)]
        fdg = np.array([(p.energy-m.energy)/(2*step) for p, m in zip(plus, minus)])
        fdh = np.column_stack([(p.gradient-m.gradient)/(2*step) for p, m in zip(plus, minus)])
        np.testing.assert_allclose(out.gradient, fdg, rtol=8e-7, atol=7e-8)
        np.testing.assert_allclose(out.hessian, fdh, rtol=8e-7, atol=8e-7)
        errors.append(np.max(abs(out.hessian-fdh)))
    assert errors[1] < errors[0]*.4
    np.testing.assert_array_equal(out.hessian, out.hessian.T)


def test_registry_lattice_and_threefold_covariance(analytic):
    model = analytic[0]; q = np.array([1.1*model.h, .174, .083])
    out = model.evaluate(q)
    for lattice in (model.geometry.a1, model.geometry.a2):
        translated = model.evaluate(q+np.r_[0., lattice])
        np.testing.assert_allclose(translated.energy, out.energy, atol=5e-12)
        np.testing.assert_allclose(translated.gradient, out.gradient, atol=5e-11)
        np.testing.assert_allclose(translated.hessian, out.hessian, atol=5e-10)
    angle = 2*np.pi/3
    rotation = np.eye(3); rotation[1:, 1:] = [[np.cos(angle), -np.sin(angle)],
                                            [np.sin(angle), np.cos(angle)]]
    rotated = model.evaluate(rotation@q)
    np.testing.assert_allclose(rotated.energy, out.energy, atol=5e-12)
    np.testing.assert_allclose(rotated.gradient, rotation@out.gradient, atol=5e-11)
    np.testing.assert_allclose(rotated.hessian, rotation@out.hessian@rotation.T, atol=5e-10)


def test_jet_components_are_per_site_and_sum_exactly(analytic):
    model = analytic[0]; out = model.evaluate([model.h, .213, .13])
    jet = sum(out.components.values())
    np.testing.assert_array_equal(jet[1:4], out.gradient)
    assert jet[0] == out.energy
    assert set(out.components) == {'pair_6', 'pair_3', 'embedding', 'angular_1', 'angular_2', 'angular_3'}


def test_direct_sum_arbitrary_vector_and_tail_refinement(analytic):
    model = analytic[0]; q = np.array([1.09*model.h, .237, .087])
    face = copy.copy(model.face)
    # Independent real-space validator already proven for the scalar projection;
    # only the explicitly supplied constant registry offset changes here.
    face._active_delta = lambda k, s: face._baseline_delta(k)+q[1:]+s*face.direction
    out = model.evaluate(q); errors = []
    for radial, layers in ((24, 24), (48, 48)):
        val = face.direct_reference(q[0], 0., radial_index=radial, layers=layers)
        energy = val.energy
        for amplitude, plane in model.moments:
            moment = copy.copy(plane.invariant); moment.interface = face
            energy += amplitude*moment.direct_value(q[0], 0., radius=radial, layers=layers)
        errors.append(abs(energy-out.energy))
    assert errors[1] < errors[0]/6
    assert errors[1] < 3e-6  # independent algebraic finite real-space tail, not roundoff


def test_reciprocal_depth_tolerance_refinement(analytic):
    loose = analytic[0]
    tight = FullRegistryInterface(build_surface(tolerance=2e-13)[0])
    q = [1.13*loose.h, .307, -.071]
    a, b = loose.evaluate(q), tight.evaluate(q)
    np.testing.assert_allclose(a.energy, b.energy, atol=3e-11, rtol=0)
    np.testing.assert_allclose(a.gradient, b.gradient, atol=4e-10, rtol=0)
    np.testing.assert_allclose(a.hessian, b.hessian, atol=8e-10, rtol=0)


def test_direct_path_has_omitted_force_not_true_stationary_saddle(analytic):
    model = analytic[0]; v = model.evaluate([model.h, .5, 0.])
    assert abs(v.gradient[1]) < 1e-10
    assert abs(v.gradient[2]) > 1.
    assert abs(v.gradient[0]) > 1.
    # Numerical nonzero and a maximum in ONE path are not full-state evidence.


def test_partial_saddle_full_index_connectivity(analytic):
    model = analytic[0]
    saddle = stationary_state(model, [1.03*model.h, *(.6*model.geometry.tau)], expected_index=1)
    assert saddle['valid']
    endpoints = saddle_downhill_endpoints(model, saddle)
    energies = sorted(e['evaluation'].energy for e in endpoints)
    assert abs(energies[0]) < 1e-9
    assert .04 < energies[1] < .08
    assert all(e['morse_index'] == 0 for e in endpoints)


def test_transverse_schur_derivative(analytic):
    model = analytic[0]; x = .18; step = 2e-5
    central = relax_at_x(model, x)
    z = central['q'][[0, 2]]
    plus, minus = (relax_at_x(model, x+sign*step, z) for sign in (1, -1))
    fd = (plus['evaluation'].gradient[1]-minus['evaluation'].gradient[1])/(2*step)
    assert fd == pytest.approx(central['schur_curvature'], rel=5e-6, abs=1e-7)


def test_source_matches_existing_path_and_full_derivatives(source):
    q = [1.04*source.h, .197, .081]; out = source.evaluate(q); step = 1e-5
    fd = np.column_stack([(source.evaluate(np.array(q)+step*v).gradient-
                          source.evaluate(np.array(q)-step*v).gradient)/(2*step) for v in np.eye(3)])
    np.testing.assert_allclose(out.hessian, fd, rtol=3e-6, atol=3e-6)
    s = .231; path = registry_path(SHOCKLEY_112); t = np.r_[0., path.direction()]
    out = source.evaluate([q[0], *(s*path.direction())])
    old = source.reference.interface_derivatives(q[0]*source.length, s*source.length, path_id=SHOCKLEY_112)
    actual = [out.energy, out.gradient[0], out.gradient@t, out.hessian[0, 0],
              out.hessian[0]@t, t@out.hessian@t]
    np.testing.assert_allclose(actual, old*np.r_[1., [source.length]*2, [source.length**2]*3], atol=2e-12)


def test_stationary_guard_rejects_flat_cutoff_region(source):
    result = stationary_state(source, [1.3*source.h, .5, -.28], expected_index=0)
    assert not result['valid']  # this is a high-index maximum, not an intact minimum


def test_static_mixed_load_unload_and_no_time_claim(analytic):
    model, units, metadata = analytic
    pristine = stationary_state(model, [model.h, 0., 0.], expected_index=0)
    load = units.traction_mpa_to_force([50., 25., 15.])
    loaded = stationary_state(model, pristine['q'], force=load, expected_index=0)
    unloaded = stationary_state(model, loaded['q'], expected_index=0)
    assert all(r['valid'] for r in (pristine, loaded, unloaded))
    np.testing.assert_allclose(unloaded['q'], pristine['q'], atol=2e-9)
    assert not metadata['physical_hz'] and not metadata['material_calibration_accepted']


@pytest.mark.parametrize('q', [[0., 0., 0.], [1., np.nan, 0.], [1., 0.]])
def test_invalid_vector_state_rejected(analytic, q):
    with pytest.raises(ValueError):
        analytic[0].evaluate(q)
