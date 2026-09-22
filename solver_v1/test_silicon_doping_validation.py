"""Source-bound support tests and independent electronic/elastic references."""
import csv
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.optimize import brentq

from .silicon_doping_research import (
    ELASTIC_SOURCE, equilibrated_charge_branches, interpolate_elastic_constants, load_elastic_samples,
)
from .silicon_doping_validation import (
    fourier_surface_compliance, two_orbital_canonical, two_orbital_product,
)


def test_v8_retains_every_original_elastic_table_value():
    root = Path(__file__).resolve().parents[1]
    with (root/'results/silicon_doping_v7/sources/jaakkola_2014_elastic.csv').open() as f:
        old = list(csv.DictReader(f))
    with ELASTIC_SOURCE.open() as f:
        new = list(csv.DictReader(f))
    assert len(old) == len(new) == 7
    for before, after in zip(old, new):
        assert all(after[k] == v for k, v in before.items())


def test_arsenic_source_temperature_exception_and_exact_endpoint():
    samples = load_elastic_samples()
    low, high = samples['As1.7'], samples['As2.5']
    assert low.temperature_max_C == 80  # Jaakkola Sec.III, explicit exception.
    low.constants_GPa(80)
    high.constants_GPa(85)
    for t in [80.0001, 85]:
        with pytest.raises(ValueError, match='As1.7'):
            low.constants_GPa(t)
    options = dict(species='As', temperature_C=85, method='linear_in_carrier_density')
    with pytest.raises(ValueError, match='As1.7'):
        interpolate_elastic_constants(carrier_density_cm3=(low.carrier_average_cm3+high.carrier_average_cm3)/2,
                                      **options)
    result = interpolate_elastic_constants(carrier_density_cm3=high.carrier_average_cm3, **options)
    assert result['source_samples'] == ['As2.5']
    assert_allclose(result['constants_GPa'], high.constants_GPa(85))


@pytest.mark.parametrize('sample_id', ['B3', 'B0.6', 'As1.7', 'As2.5', 'P4.1', 'P4.7', 'P7.5'])
def test_crack_coefficient_from_independent_displacement_halfspace(sample_id):
    sample = load_elastic_samples()[sample_id]
    for t in [sample.temperature_min_C, 25, sample.temperature_max_C]:
        reference = fourier_surface_compliance(*sample.constants_GPa(t))
        assert_allclose(reference['H_GPa_inv'], sample.crack_field(t).H_GPa_inv, rtol=3e-12)
        assert reference['subspace_residual'] < 1e-12


def test_halfspace_isotropic_repeated_root_is_resolved_as_subspace():
    young, nu = 170., .23
    mu = young/(2*(1+nu))
    lam = young*nu/((1+nu)*(1-2*nu))
    result = fourier_surface_compliance(lam+2*mu, lam, mu)
    assert_allclose(result['H_GPa_inv'], (1-nu**2)/young, rtol=1e-12)


@pytest.mark.parametrize('q', [[-.1, .05], [0., 0.], [.12, -.08]])
@pytest.mark.parametrize('mu', [-.1, 0., .1])
@pytest.mark.parametrize('kt', [.01, .025, .08])
def test_charge_sector_sum_matches_independent_fermi_product(q, mu, kt):
    sectors = two_orbital_canonical(q, kt)
    actual = equilibrated_charge_branches(*sectors, excess_electron_counts=[0, 1, 2], mu_eV=mu, kBT_eV=kt)
    reference = two_orbital_product(q, mu, kt)
    for key in ['grand_potential_eV', 'gradient', 'hessian', 'mean_excess_electrons']:
        assert_allclose(actual[key], reference[key], rtol=1e-12, atol=1e-14)
    weights = actual['branch_probability']
    variance = weights@(np.arange(3)-actual['mean_excess_electrons'])**2
    assert_allclose(variance/kt, reference['dmean_dmu'], rtol=2e-12)


def test_fixed_mean_charge_legendre_curvature_has_schur_correction():
    # Fixed MEAN N is a Legendre constraint, not an exact canonical N sector.
    q, mu, kt, step = np.array([.03, -.02]), .01, .025, 1e-5
    r = two_orbital_product(q, mu, kt)
    dn = r['dmean_dq']
    expected = r['hessian']+np.outer(dn, dn)/r['dmean_dmu']
    def constrained(point):
        chosen = brentq(lambda m: two_orbital_product(point, m, kt)['mean_excess_electrons']
                       -r['mean_excess_electrons'], -.8, .8, xtol=1e-14)
        return two_orbital_product(point, chosen, kt)
    for i, axis in enumerate(np.eye(2)):
        derivative = (constrained(q+step*axis)['gradient']-constrained(q-step*axis)['gradient'])/(2*step)
        assert_allclose(derivative, expected[:, i], rtol=1e-7, atol=1e-8)
    assert np.linalg.norm(expected-r['hessian']) > .1


def test_normalized_truncated_sector_sum_can_still_be_wrong():
    kt = .025
    # Two zero-energy independent orbitals have sector degeneracy 1,2,1.
    f = [0., -kt*np.log(2), 0.]
    full = equilibrated_charge_branches(f, np.zeros((3, 1)), np.zeros((3, 1, 1)),
                                       excess_electron_counts=[0, 1, 2], mu_eV=0, kBT_eV=kt)
    cut = equilibrated_charge_branches(f[:2], np.zeros((2, 1)), np.zeros((2, 1, 1)),
                                      excess_electron_counts=[0, 1], mu_eV=0, kBT_eV=kt)
    assert_allclose(sum(cut['branch_probability']), 1)
    assert_allclose(cut['mean_excess_electrons'], 2/3)
    assert_allclose(full['mean_excess_electrons'], 1)
    assert_allclose(cut['grand_potential_eV']-full['grand_potential_eV'], kt*np.log(4/3))
