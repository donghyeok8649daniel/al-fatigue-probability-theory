"""Independent identities and failure cases for the doping research bridge."""
import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from .silicon_doping_research import (
    DopantPopulation, charge_balance_cm3, directional_young_GPa,
    dopant_site_statistics, equilibrated_charge_branches,
    excess_electrons_for_material_volume, interpolate_elastic_constants,
    load_elastic_samples, silicon_site_density_cm3,
)
from .silicon_specimen_research import AnisotropicModeI, cubic_tensor


def test_site_count_density_and_volume_are_consistent():
    a = 5.431
    density = silicon_site_density_cm3(a)
    assert_allclose(density*a**3*1e-24, 8)
    s = dopant_site_statistics(density/8, lattice_A=a, site_count=8)
    assert_allclose(s['expected_count'], 1)
    assert_allclose(s['material_volume_A3'], a**3)
    assert_allclose(s['one_dopant_concentration_cm3'], density/8)


def test_exact_binomial_does_not_round_dopant_atoms():
    x, n, a = .19, 4, 5.431
    s = dopant_site_statistics(x*silicon_site_density_cm3(a), lattice_A=a, site_count=n)
    distribution = np.array([math.comb(n, j)*x**j*(1-x)**(n-j) for j in range(n+1)])
    assert_allclose(s['probability_at_least_one'], sum(distribution[1:]))
    assert_allclose(s['expected_count'], distribution@np.arange(n+1))
    assert_allclose(s['variance_count'], distribution@(np.arange(n+1)-n*x)**2)


def test_dilute_and_full_occupancy_limits():
    density = silicon_site_density_cm3(5.431)
    dilute = dopant_site_statistics(density*1e-20, lattice_A=5.431, site_count=100)
    assert_allclose(dilute['probability_at_least_one'], 1e-18, rtol=1e-14)
    for x in (0., 1.):
        s = dopant_site_statistics(density*x, lattice_A=5.431, site_count=100)
        assert s['probability_at_least_one'] == x


@pytest.mark.parametrize('density,count', [(-1, 4), (np.inf, 4), (1e24, 4), (1e19, 0), (1e19, 4.5), (1e19, True)])
def test_invalid_chemical_occupancy_is_rejected(density, count):
    with pytest.raises(ValueError):
        dopant_site_statistics(density, lattice_A=5.431, site_count=count)


def test_unknown_ionization_does_not_become_full_ionization():
    with pytest.raises(ValueError, match='unknown'):
        charge_balance_cm3([DopantPopulation('B', 1e19)], electrons_cm3=0, holes_cm3=1e19)
    with pytest.raises(ValueError):
        DopantPopulation('B', 1e19, 1.1)


def test_compensation_preserves_chemical_atoms():
    ntype = DopantPopulation('P', 2e19, .75)
    ptype = DopantPopulation('B', 1e19, .5)
    assert charge_balance_cm3([ntype, ptype], electrons_cm3=1e19, holes_cm3=0) == 0
    assert sum(d.chemical_density_cm3 for d in [ntype, ptype]) == 3e19
    assert charge_balance_cm3([ntype, ptype], electrons_cm3=0, holes_cm3=0) == 1e19


def test_charge_number_uses_material_volume_and_sign():
    for sign in [-1, 1]:
        assert_allclose(excess_electrons_for_material_volume(sign*1e21, material_volume_A3=100), sign*.1)
    with pytest.raises(ValueError):
        excess_electrons_for_material_volume(1e19, material_volume_A3=0)


def test_measured_sample_source_values_and_temperature_units():
    samples = load_elastic_samples()
    assert len(samples) == 7
    p = samples['P7.5']
    assert p.carrier_nominal_cm3 == 7.5e19
    assert p.carrier_average_cm3 == 7.47e19
    assert_allclose(p.constants_GPa(25), [161.4, 66.1, 78.5])
    assert_allclose(p.constants_GPa(35)[0], 161.4*(1-30.7e-6*10-78e-9*100))
    for t in [-40.001, 85.001, np.nan]:
        with pytest.raises(ValueError):
            p.constants_GPa(t)


def test_interpolation_is_explicit_same_species_and_bounded():
    samples = load_elastic_samples()
    lo, hi = samples['P4.1'], samples['P4.7']
    result = interpolate_elastic_constants('P', carrier_density_cm3=(lo.carrier_average_cm3+hi.carrier_average_cm3)/2,
                                          temperature_C=25, method='linear_in_carrier_density')
    assert_allclose(result['constants_GPa'], (lo.constants_GPa(25)+hi.constants_GPa(25))/2)
    assert_allclose(result['upper_weight'], .5)
    for species, n in [('P', 1e19), ('P', 1e21), ('Sb', 1e19), ('B', 7e19)]:
        with pytest.raises(ValueError):
            interpolate_elastic_constants(species, carrier_density_cm3=n, temperature_C=25,
                                          method='linear_in_carrier_density')
    with pytest.raises(TypeError):
        interpolate_elastic_constants('P', carrier_density_cm3=5e19, temperature_C=25)


def test_young_modulus_matches_independent_tensor_compliance():
    c = np.array([163., 65.4, 79.2])
    tensor = cubic_tensor(*c, frame=np.eye(3))
    # Use orthonormal symmetric (Kelvin) tensors, avoiding Voigt shear factors.
    basis = np.zeros((6, 3, 3))
    basis[np.arange(3), np.arange(3), np.arange(3)] = 1
    for k, (i, j) in enumerate([(1, 2), (0, 2), (0, 1)], start=3):
        basis[k, i, j] = basis[k, j, i] = 1/np.sqrt(2)
    kelvin = np.einsum('aij,ijkl,bkl->ab', basis, tensor, basis)
    for v in [[1, 0, 0], [1, 1, 0], [1, 1, 1], [2, -3, 7]]:
        direction = np.asarray(v)/np.linalg.norm(v)
        traction = np.einsum('aij,i,j->a', basis, direction, direction)
        expected = 1/(traction@np.linalg.solve(kelvin, traction))
        assert_allclose(directional_young_GPa(*c, v), expected, rtol=1e-13)
    e100, e110, e111 = [directional_young_GPa(*c, v) for v in [[1, 0, 0], [1, 1, 0], [1, 1, 1]]]
    assert_allclose(4/e110, 1/e100+3/e111)


def test_measured_elasticity_J_and_traction_at_temperature_endpoints():
    for sample in load_elastic_samples().values():
        for t in (sample.temperature_min_C, 25., sample.temperature_max_C):
            field = sample.crack_field(t)
            assert_allclose(field.contour_J(.6, radius_A=40, points=2048), field.energy_release_J_m2(.6), rtol=2e-10)
            stress = field.field([[10., 0]], .6)[2]
            assert_allclose(stress[0, 1], 60/np.sqrt(20*np.pi), rtol=1e-12)


def branch_result(q, mu=.04):
    matrices = np.array([[[2., .3], [.3, 1.]], [[1., -.2], [-.2, 3.]]])
    linear = np.array([[.2, -.1], [-.3, .25]])
    values = np.array([0., .03])+linear@q+.5*np.einsum('i,aij,j->a', q, matrices, q)
    gradients = linear+np.einsum('aij,j->ai', matrices, q)
    return equilibrated_charge_branches(values, gradients, matrices,
                                        excess_electron_counts=[0., 1.], mu_eV=mu, kBT_eV=.025)


def test_electronic_grand_potential_force_and_full_hessian():
    q, step = np.array([.08, -.03]), 1e-5
    r = branch_result(q)
    for i, axis in enumerate(np.eye(2)):
        plus, minus = branch_result(q+step*axis), branch_result(q-step*axis)
        assert_allclose((plus['grand_potential_eV']-minus['grand_potential_eV'])/(2*step), r['gradient'][i], rtol=2e-7)
        assert_allclose((plus['gradient']-minus['gradient'])/(2*step), r['hessian'][:, i], rtol=2e-7)
    numerical_mu = (branch_result(q, .04+step)['grand_potential_eV']-branch_result(q, .04-step)['grand_potential_eV'])/(2*step)
    assert_allclose(numerical_mu, -r['mean_excess_electrons'], rtol=1e-7)


def test_charge_fluctuations_change_curvature_and_energy_shift_is_harmless():
    args = dict(gradients=[[-1.], [1.]], hessians=[[[2.]], [[2.]]],
                excess_electron_counts=[0, 1], mu_eV=0, kBT_eV=.5)
    r = equilibrated_charge_branches([0., 0.], **args)
    assert_allclose(r['hessian'], [[0.]], atol=1e-15)
    shifted = equilibrated_charge_branches([1e8, 1e8], **args)
    assert_allclose(shifted['branch_probability'], [.5, .5], atol=1e-15)
    single = equilibrated_charge_branches([2.], [[3.]], [[[4.]]], excess_electron_counts=[-1.], mu_eV=.5, kBT_eV=.02)
    assert_allclose(single['grand_potential_eV'], 2.5)
    assert_allclose(single['hessian'], [[4.]])


def test_invalid_electronic_data_rejected():
    with pytest.raises(ValueError):
        equilibrated_charge_branches([0.], [[0.]], [[[0.]]], excess_electron_counts=[0.], mu_eV=0, kBT_eV=0)
    with pytest.raises(ValueError):
        equilibrated_charge_branches([0.], [[0., 0.]], [[[1., 2.], [3., 4.]]], excess_electron_counts=[0.], mu_eV=0, kBT_eV=.1)


@pytest.mark.parametrize('counts', [[0., .0016], [1., 1.]])
def test_fractional_DFT_sampling_and_duplicate_sectors_are_not_partition_states(counts):
    with pytest.raises(ValueError, match='distinct integer'):
        equilibrated_charge_branches([0., 1.], [[0.], [0.]], [[[1.]], [[1.]]],
                                     excess_electron_counts=counts, mu_eV=0., kBT_eV=.025)
