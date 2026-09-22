import numpy as np
import pytest
from numpy.testing import assert_allclose

from .silicon_concentration_research import (
    chemical_concentration_cm3, constrained_energy_derivative,
    deterministic_substitution_order, rigid_opening_geometry,
    independent_site_configuration_coverage)


def test_reference_material_volume_preserves_count_under_vacuum_opening():
    cell = np.diag([10., 12., 20.]); r = np.array([[1., 2., 1.], [1., 2., 19.]])
    density = chemical_concentration_cm3(2, reference_material_volume_A3=np.linalg.det(cell))
    for q in [0., 1., 20.]:
        new_r, new_c, dr, dc = rigid_opening_geometry(r, cell, q)
        assert_allclose(new_r, r)
        assert_allclose(new_c[2, 2], 20+q)
        assert_allclose(density, chemical_concentration_cm3(2, reference_material_volume_A3=2400.), rtol=1e-14)
        assert_allclose(dr, 0); assert dc.sum() == 1


def test_nonaffine_work_matches_independent_periodic_harmonic_bond():
    # Bond crosses the z image boundary: E = k/2 (L+z0-z1-l0)^2.
    c = np.diag([2., 3., 10.]); r = np.array([[.3, .4, .7], [.3, .4, 8.6]])
    k, natural = 3.4, 1.5
    d = c[2, 2]+r[0, 2]-r[1, 2]
    derivative = k*(d-natural)
    f = np.array([[0., 0., -derivative], [0., 0., derivative]])
    stress = np.zeros((3, 3)); stress[2, 2] = derivative*d/np.linalg.det(c)
    dc = np.zeros((3, 3)); dc[2, 2] = 1
    value = constrained_energy_derivative(r, f, c, stress, position_tangent=np.zeros_like(r), cell_tangent=dc)
    assert_allclose(value, derivative, atol=1e-14)
    assert abs(6*stress[2, 2]-derivative) > 1
    for h in [1e-3, 1e-4]:
        fd = (.5*k*(d+h-natural)**2-.5*k*(d-h-natural)**2)/(2*h)
        assert_allclose(value, fd, rtol=1e-10)


def test_general_cell_work_is_coordinate_rotation_invariant():
    rng = np.random.default_rng(47)
    c = np.array([[4., .1, .2], [.3, 5., .4], [.1, .2, 6.]])
    r, f, dr, dc = rng.normal(size=(4, 3)), rng.normal(size=(4, 3)), rng.normal(size=(4, 3)), rng.normal(size=(3, 3))
    s = rng.normal(size=(3, 3)); s = (s+s.T)/2
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    a = constrained_energy_derivative(r, f, c, s, position_tangent=dr, cell_tangent=dc)
    b = constrained_energy_derivative(r@rotation, f@rotation, c@rotation, rotation.T@s@rotation,
                                     position_tangent=dr@rotation, cell_tangent=dc@rotation)
    assert_allclose(a, b, rtol=1e-13, atol=1e-13)


def test_skew_cell_and_nonaffine_motion_against_direct_bond_length_derivative():
    c = np.array([[4., .1, .2], [.3, 5., .4], [.1, .2, 6.]])
    r = np.array([[.4, .5, .6], [1.3, 2.1, 3.2]])
    dr = np.array([[.2, -.1, .3], [-.4, .2, .1]])
    dc = np.array([[.1, .2, -.3], [.4, -.2, .5], [.2, -.3, .1]])
    image = np.array([1., -1., 0.])
    d = r[1]+image@c-r[0]
    length = np.linalg.norm(d); spring, natural = 2.3, 3.2
    g = spring*(1-natural/length)*d
    forces = np.array([g, -g])
    stress = np.outer(g, d)/np.linalg.det(c)
    actual = constrained_energy_derivative(r, forces, c, stress, position_tangent=dr, cell_tangent=dc)
    expected = g@(dr[1]+image@dc-dr[0])
    assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)
    def energy(q):
        vector = (r[1]+q*dr[1])+image@(c+q*dc)-(r[0]+q*dr[0])
        return .5*spring*(np.linalg.norm(vector)-natural)**2
    assert_allclose(actual, (energy(1e-5)-energy(-1e-5))/2e-5, rtol=1e-8)


def test_configuration_selection_keeps_distinct_nested_sites_and_explicit_interface():
    r = np.array([[x, y, z] for x in [1., 4., 7.] for y in [1., 4.] for z in [1., 5., 9.]])
    c = np.diag([9., 8., 10.])
    for arrangement in ['interface_cluster', 'interface_spread', 'bulk_spread']:
        order = deterministic_substitution_order(r, c, arrangement=arrangement)
        assert len(set(order)) == len(order)
        if arrangement == 'interface_spread': assert np.all(np.isin(r[order, 2], [1, 9]))
        if arrangement == 'bulk_spread': assert r[order[0], 2] == 5


def test_small_negative_periodic_coordinate_needs_explicit_modulo():
    # ASE wrap(eps=1e-7) can retain a negative Cartesian y near this size.
    c = np.diag([11.6, 13.4, 28.4])
    r = np.array([[1., -8.5e-7, 1.2], [2., 3., 27.]])
    with pytest.raises(ValueError): rigid_opening_geometry(r, c, .5)
    wrapped = r % np.diag(c)
    changed, _, _, _ = rigid_opening_geometry(wrapped, c, .5)
    assert np.all(changed >= 0)
    assert_allclose((changed-r)/np.diag(c), np.rint((changed-r)/np.diag(c)), atol=1e-15)


@pytest.mark.parametrize('count', [-1, .5, True])
def test_fractional_or_invalid_chemical_atoms_are_rejected(count):
    with pytest.raises(ValueError): chemical_concentration_cm3(count, reference_material_volume_A3=100.)


@pytest.mark.parametrize('x',[0.,1e-9,.2,.8,1.])
def test_complete_chemical_pattern_enumeration_conserves_probability(x):
    from itertools import combinations
    n=6
    patterns=[p for k in range(n+1) for p in combinations(range(n),k)]
    result=independent_site_configuration_coverage(n,x,patterns)
    assert_allclose(result['covered_mass'],1,atol=5e-15,rtol=0)
    assert_allclose(result['missing_mass'],0,atol=5e-15,rtol=0)
    expected_count=sum(len(p)*w for p,w in zip(result['unique_initial_patterns'],result['weights']))
    assert_allclose(expected_count,n*x,atol=1e-14,rtol=1e-14)


def test_selected_chemical_patterns_keep_missing_mass_and_remove_exact_duplicates():
    r=independent_site_configuration_coverage(3,.25,[[],[0],[0]])
    assert r['duplicates_removed']==1
    assert_allclose(r['covered_mass'],.75**3+.25*.75**2,atol=1e-15)
    assert_allclose(r['missing_mass'],1-r['covered_mass'],atol=1e-15)
    assert r['weights'].sum()<1
