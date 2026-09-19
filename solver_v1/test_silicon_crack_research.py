"""Independent checks for spatial energy, boundaries and local coordinates."""
import numpy as np
import pytest

from results.silicon_wafer_feasibility.run_static_probe import pair, source_parameters
from .silicon_environment_research import environment_jet
from .silicon_crack_research import (
    RelaxedCoordinates, SpatialSW, diamond_crack_strip, relax_atoms,
)


@pytest.fixture
def material():
    from scipy.optimize import minimize_scalar
    p, _ = source_parameters()
    bond = minimize_scalar(lambda r: float(pair(r, p)), bounds=(2.2, 2.5),
                           method='bounded', options={'xatol':1e-14}).x
    return p, 4*bond/np.sqrt(3)


def test_spatial_force_matches_independent_site_jets(material):
    p, _ = material
    positions = np.array([[0., 0., 0.], [1.4, 1.5, 1.3], [1.6, -1.5, -1.2],
                          [-1.5, 1.2, -1.7], [-1.4, -1.4, 1.5]])
    expected_energy, expected_gradient = 0., np.zeros_like(positions)
    for i in range(len(positions)):
        ids = np.array([j for j in range(len(positions)) if j != i])
        jet = environment_jet(positions[ids]-positions[i], p, 'direct')
        expected_energy += jet.value
        expected_gradient[ids] += jet.gradient.reshape(-1, 3)
        expected_gradient[i] -= jet.gradient.reshape(-1, 3).sum(axis=0)
    model = SpatialSW(p)
    current = model.evaluate(positions)
    np.testing.assert_allclose(current.energy, expected_energy, atol=2e-13)
    np.testing.assert_allclose(current.gradient, expected_gradient, atol=2e-13)
    np.testing.assert_allclose(model.evaluate(positions, representation='direct').gradient,
                               expected_gradient, atol=2e-13)
    np.testing.assert_allclose(current.gradient.sum(axis=0), 0, atol=1e-13)
    np.testing.assert_allclose(np.cross(positions, current.gradient).sum(axis=0), 0, atol=1e-13)
    theta = .731
    rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta), np.cos(theta), 0], [0, 0, 1.]])
    rotated = model.evaluate(positions@rotation.T+[9, -3, 2])
    np.testing.assert_allclose(rotated.energy, current.energy, atol=2e-13)
    np.testing.assert_allclose(rotated.gradient, current.gradient@rotation.T, atol=3e-13)


def test_periodic_images_skin_rebuild_and_new_neighbors(material):
    p, _ = material
    positions = np.array([[0., 0., .2], [2.2, .1, 2.1], [-2.1, .5, 3.8], [0., 4.5, 1.]])
    period = 4.6
    model = SpatialSW(p, front_period=period, skin=.7)
    for delta in [0., .1, -.8, -2.5, -.02]:
        displaced = positions.copy()
        displaced[-1, 1] += delta
        current = model.evaluate(displaced)
        fresh = SpatialSW(p, front_period=period, skin=0).evaluate(displaced, representation='direct')
        np.testing.assert_allclose(current.energy, fresh.energy, atol=2e-13)
        np.testing.assert_allclose(current.gradient, fresh.gradient, atol=2e-12)
    shifted = displaced.copy()
    shifted[1, 2] += 5*period
    shifted[3, 2] -= 3*period
    wrapped = model.evaluate(shifted)
    np.testing.assert_allclose(wrapped.energy, current.energy, atol=2e-12)
    np.testing.assert_allclose(wrapped.gradient, current.gradient, atol=2e-12)
    assert model.rebuilds >= 3


def test_sparse_hessian_keeps_fixed_center_contributions(material):
    p, _ = material
    positions = np.array([[0., 0., 0.], [1.4, 1.5, 1.3], [1.6, -1.5, -1.2],
                          [-1.5, 1.2, -1.7], [-1.4, -1.4, 1.5]])
    model = SpatialSW(p)
    free = np.array([1, 2, 4])
    hessian = model.hessian(positions, free_atoms=free).toarray()
    direction = np.sin(np.arange(9)+.4)
    delta = np.zeros_like(positions)
    delta[free] = direction.reshape(-1, 3)
    errors = []
    for h in [2e-3, 1e-3, 5e-4]:
        derivative = (model.evaluate(positions+h*delta).gradient[free]
                      - model.evaluate(positions-h*delta).gradient[free]).ravel()/(2*h)
        errors.append(np.max(abs(derivative-hessian@direction)))
    assert errors[-1] < 1e-4
    assert errors[0]/errors[-1] > 12
    np.testing.assert_allclose(hessian, hessian.T, atol=2e-13)
    full = model.hessian(positions).toarray()
    dofs = (free[:, None]*3+np.arange(3)).ravel()
    np.testing.assert_allclose(hessian, full[np.ix_(dofs, dofs)], atol=2e-13)


def test_strip_and_rigid_separation_normalization(material):
    p, lattice = material
    strip = diamond_crack_strip(lattice, nx=5, ny=3, nz=2)
    assert len(strip.positions) == 5*3*2*12
    assert len(strip.crossing_bonds) > 0
    np.testing.assert_allclose(strip.frame@strip.frame.T, np.eye(3), atol=2e-15)
    assert np.linalg.det(strip.frame) > 0
    model = SpatialSW(p, front_period=strip.front_period)
    initial = model.evaluate(strip.positions)
    separated = strip.positions.copy()
    separated[separated[:, 1] > 0, 1] += 5.
    opened = model.evaluate(separated)
    # At ideal shuffle geometry all original SW angles are tetrahedral. Each
    # severed normal bond costs exactly -phi(r0), with no area multiplier.
    bond = lattice*np.sqrt(3)/4
    expected = -len(strip.crossing_bonds)*float(pair(bond, p))
    np.testing.assert_allclose(opened.energy-initial.energy, expected, atol=2e-10)
    healed = model.evaluate(strip.positions)
    np.testing.assert_allclose(healed.energy, initial.energy, atol=2e-11)
    np.testing.assert_allclose(healed.gradient, initial.gradient, atol=2e-11)


@pytest.mark.parametrize('coherent', [False, True])
def test_local_coordinates_exact_constraint_and_virtual_work(material, coherent):
    _, lattice = material
    strip = diamond_crack_strip(lattice, nx=4, ny=3, nz=3, grip_width=3.)
    bonds = strip.crossing_bonds
    eligible = bonds[~np.any(strip.fixed[bonds], axis=1)]
    bond = eligible[len(eligible)//2]
    coordinates = RelaxedCoordinates(strip.positions, strip.fixed,
        groups=strip.front_groups if coherent else None, bond=bond, components=(1, 0))
    q = np.array([3.1, .12])
    x = np.sin(np.arange(coordinates.dimension))*.01
    positions = coordinates.positions(x, q)
    np.testing.assert_allclose((positions[bond[1]]-positions[bond[0]])[[1, 0]], q, atol=2e-15)
    np.testing.assert_array_equal(positions[strip.fixed], strip.positions[strip.fixed])
    np.testing.assert_allclose(coordinates.encode(positions), x, atol=3e-15)
    gradient = np.cos(np.arange(positions.size)).reshape(positions.shape)
    gx, gq = coordinates.pullback(gradient)
    tangent = coordinates.tangent_matrix()
    np.testing.assert_allclose(tangent.T@gradient.ravel(), gx, atol=2e-15)
    dx, dq = np.cos(np.arange(len(x))), np.array([.4, -.3])
    h = 1e-4
    work = np.sum(gradient*(coordinates.positions(x+h*dx, q+h*dq)-positions))/h
    np.testing.assert_allclose(work, gx@dx+gq@dq, atol=2e-9)
    np.testing.assert_allclose((coordinates.positions(x+h*dx, q)-positions).ravel()/h,
                               tangent@dx, atol=4e-11)
    if not coherent:
        same_x = np.isclose(strip.positions[bonds[:, 0], 0], strip.positions[bond[0], 0])
        gaps = strip.bond_gaps(positions)[same_x]
        assert np.ptp(gaps) > .5  # A single tip bond can move without the entire front.


def test_force_checked_relaxation_does_not_hide_iteration_failure(material):
    p, _ = material
    positions = np.array([[0., 0., 0.], [2.6, .2, .1]])
    coordinates = RelaxedCoordinates(positions, [True, False])
    model = SpatialSW(p)
    relaxed, info = relax_atoms(model, coordinates, tolerance=1e-8)
    assert info['force_converged']
    assert not info['stability_certified']  # A dimer has two rotational zero modes.
    assert info['physical_time_seconds'] is None
    assert np.max(abs(model.evaluate(relaxed).gradient)) < 1e-8
    _, unfinished = relax_atoms(model, coordinates, tolerance=1e-12, maxiter=1)
    assert not unfinished['force_converged']
    assert not unfinished['optimizer_success']


def test_invalid_inputs_are_rejected(material):
    p, lattice = material
    with pytest.raises(ValueError):
        SpatialSW(p, front_period=-1)
    with pytest.raises(ValueError):
        SpatialSW(p).evaluate(np.zeros((2, 3)))
    with pytest.raises(ValueError):
        diamond_crack_strip(lattice, nx=1, ny=1, nz=1, grip_width=20)
    with pytest.raises(ValueError):
        RelaxedCoordinates(np.zeros((2, 3)), [True, False], groups=[0, 0])
    with pytest.raises(ValueError):
        RelaxedCoordinates(np.zeros((2, 3)), [True, False], bond=(0, 1))
