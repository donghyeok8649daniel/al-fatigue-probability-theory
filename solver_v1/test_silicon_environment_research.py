"""Independent static/derivative checks; no Si material or clock certification."""
import numpy as np
import pytest
from scipy.optimize import brentq, root

from results.silicon_wafer_feasibility.run_static_probe import (
    RigidSlab, analytic_shuffle_hessian, pair, source_parameters,
)
from .silicon_environment_research import (
    DiamondCell, Jet, RigidInterface, environment_jet, relaxed_hessian,
)


@pytest.fixture(scope='module')
def reference():
    p, _ = source_parameters()
    bond = brentq(lambda x: float(pair(x, p, derivative=True)), 2.2, 2.5, xtol=1e-14)
    return p, bond, 4*bond/np.sqrt(3)


def test_jet_mixed_chain_rule_against_closed_form():
    x, y = (Jet.variable(v, 2, i) for i, v in enumerate([.3, -.4]))
    f = (x*y).exp() + 1/(2+x-y)
    z, e, v = 2.7, np.exp(-.12), np.array([1., -1.])
    expected_g = e*np.array([-.4, .3])-v/z**2
    expected_h = e*np.array([[.16, .88], [.88, .09]])+2*np.outer(v, v)/z**3
    np.testing.assert_allclose(f.gradient, expected_g, atol=1e-14)
    np.testing.assert_allclose(f.hessian, expected_h, atol=1e-14)


@pytest.mark.parametrize('scale', [.88, 1., 1.13])
def test_full_cartesian_force_hessian_identity_rotation_and_torque(reference, scale):
    p, bond, _ = reference
    tetra = np.array([[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1.]])*bond/np.sqrt(3)
    deformation = np.array([[scale, .11, 0.], [0., 1., -.04], [.02, 0., 1.]])
    vectors = tetra@deformation.T
    a, b = (environment_jet(vectors, p, method) for method in ['direct', 'moments'])
    assert abs(a.value-b.value) < 1e-12
    np.testing.assert_allclose(a.gradient, b.gradient, atol=2e-12)
    np.testing.assert_allclose(a.hessian, b.hessian, atol=3e-11)
    np.testing.assert_allclose(a.hessian, a.hessian.T, atol=1e-13)
    rotation, _ = np.linalg.qr(np.array([[1,2,3.], [-2,1,4.], [3,4,-1.]]))
    rotated = environment_jet(vectors@rotation.T, p)
    transform = np.kron(np.eye(4), rotation)
    np.testing.assert_allclose(rotated.gradient, transform@a.gradient, atol=2e-12)
    np.testing.assert_allclose(rotated.hessian, transform@a.hessian@transform.T, atol=3e-11)
    torque = np.cross(vectors, -a.gradient.reshape(-1, 3)).sum(axis=0)
    np.testing.assert_allclose(torque, 0., atol=2e-13)


def test_cutoff_and_zero_neighbor_derivatives(reference):
    p, _, _ = reference
    cutoff = p['sigma']*p['a']
    for vectors in [np.empty((0,3)), np.array([[cutoff+.1, 0., 0.]])]:
        j = environment_jet(vectors, p)
        assert j.value == 0.
        assert np.all(j.gradient == 0.) and np.all(j.hessian == 0.)
    with pytest.raises(ValueError, match='coincident'):
        environment_jet(np.zeros((1,3)), p)


def test_diamond_internal_relaxation_and_direct_finite_strain(reference):
    p, bond, lattice = reference
    cell = DiamondCell(p, lattice)
    base = cell.evaluate(np.zeros(9))
    direct = cell.evaluate(np.zeros(9), 'direct')
    assert abs(base.value-4*float(pair(bond, p))) < 2e-12
    np.testing.assert_allclose(base.gradient, 0., atol=5e-12)
    np.testing.assert_allclose(base.hessian, direct.hessian, atol=3e-11)
    effective, response = relaxed_hessian(base, range(6))
    assert np.linalg.eigvalsh(base.hessian[6:,6:])[0] > 0
    assert effective[5,5] < .6*base.hessian[5,5]
    np.testing.assert_allclose(effective[:3,:3], base.hessian[:3,:3], atol=1e-12)
    errors = []
    for gamma in [2e-3, 1e-3, 5e-4]:
        state = np.zeros(9); state[5] = gamma
        start = response[:,5]*gamma
        def fun(u):
            q = state.copy(); q[6:] = u
            j = cell.evaluate(q, 'direct')
            return j.gradient[6:], j.hessian[6:,6:]
        solved = root(fun, start, jac=True, options={'xtol': 1e-10})
        assert np.max(abs(solved.fun)) < 1e-11
        state[6:] = solved.x
        energy = cell.evaluate(state, 'direct', derivatives=False)
        curvature = 2*(energy-base.value)/gamma**2
        errors.append(abs(curvature-effective[5,5]))
        np.testing.assert_allclose(solved.x/gamma, response[:,5], atol=2e-5)
    assert errors[-1] < 2e-5
    assert errors[-1] < errors[0]/8


@pytest.mark.parametrize('kind', ['shuffle', 'glide'])
def test_interface_derivatives_against_independent_slab_and_mesh_extent(reference, kind):
    p, bond, lattice = reference
    old = RigidSlab(p, lattice, kind, repeats=1, periods=3, image_shell=3)
    interface = RigidInterface(old, image_shell=5)
    q = np.array([.3, .2*old.period, 0.])
    j = interface.evaluate(q)
    direct = interface.evaluate(q, 'direct')
    assert abs(j.value-old.energy(*q[:2])) < 2e-11
    np.testing.assert_allclose(j.gradient, direct.gradient, atol=2e-12)
    np.testing.assert_allclose(j.hessian, direct.hessian, atol=2e-11)
    fd_errors = []
    for h in [1e-3, 5e-4, 2.5e-4]:
        directions = np.eye(2)*h
        grad = [(old.energy(*(q[:2]+d))-old.energy(*(q[:2]-d)))/(2*h)
                for d in directions]
        fd_errors.append(float(np.max(abs(grad-j.gradient[:2]))))
    assert fd_errors[-1] < 2e-6
    assert fd_errors[-1] < fd_errors[0]/8
    def independent_energy(state):
        positions = old.positions.copy()
        positions[old.upper] += state@interface.frame
        return (old.total_energy(positions)-old.reference)/old.cells
    hessian_errors = []
    for h in [1e-3, 5e-4]:
        basis = np.eye(3)*h
        hessian = np.empty((3,3))
        for a in range(3):
            for b in range(3):
                if a == b:
                    hessian[a,a] = (independent_energy(q+basis[a])-2*independent_energy(q)
                                    + independent_energy(q-basis[a]))/h**2
                else:
                    hessian[a,b] = sum(sa*sb*independent_energy(q+sa*basis[a]+sb*basis[b])
                                      for sa in [-1,1] for sb in [-1,1])/(4*h**2)
        hessian_errors.append(float(np.max(abs(hessian-j.hessian))))
    assert hessian_errors[-1] < 3e-5
    assert hessian_errors[-1] < hessian_errors[0]/2
    bigger = RigidInterface(RigidSlab(p, lattice, kind, repeats=2, periods=4, image_shell=3), image_shell=5)
    bj = bigger.evaluate(q)
    assert abs(bj.value-j.value) < 2e-11
    np.testing.assert_allclose(bj.gradient, j.gradient, atol=2e-11)
    np.testing.assert_allclose(bj.hessian, j.hessian, atol=2e-10)
    if kind == 'shuffle':
        equilibrium = interface.evaluate(np.zeros(3))
        np.testing.assert_allclose(equilibrium.hessian[:2,:2], analytic_shuffle_hessian(bond, p), atol=2e-11)
        assert abs(j.gradient[2]) < 1e-12
    else:
        assert abs(j.gradient[2]) > 1.  # Omitted direction is not stationary.


def test_insufficient_images_and_unstable_elimination_are_rejected(reference):
    p, _, lattice = reference
    with pytest.raises(ValueError, match='shell'):
        DiamondCell(p, lattice, image_shell=1).evaluate(np.zeros(9))
    j = Jet(0., np.zeros(2), np.diag([1., -1.]))
    with pytest.raises(ValueError, match='positive curvature'):
        relaxed_hessian(j, [0])
    with pytest.raises(ValueError, match='distinct'):
        relaxed_hessian(j, [0, 0])
