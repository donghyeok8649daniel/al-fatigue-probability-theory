import math

import numpy as np
import pytest

from theory.fcc_dislocation_initiation import (
    AluminumSlipInitiationParameters,
    axial_tmw_initiation_life,
    empirical_first_passage_shape,
    fcc_axial_tmw_lives,
    tmw_cycles_to_initiation,
)


def test_tmw_formula_has_expected_inverse_square_stress_dependence():
    args = dict(
        shear_modulus_pa=26.0e9,
        poisson_ratio=0.33,
        surface_energy_j_m2=1.14,
        burgers_vector_m=0.286e-9,
    )
    low = tmw_cycles_to_initiation(8.0e6, **args)
    high = tmw_cycles_to_initiation(16.0e6, **args)
    assert np.isclose(low / high, 4.0)


def test_tmw_threshold_returns_infinite_life_for_this_mechanism_only():
    life = tmw_cycles_to_initiation(
        2.0e6,
        shear_modulus_pa=26.0e9,
        poisson_ratio=0.33,
        surface_energy_j_m2=1.14,
        burgers_vector_m=0.286e-9,
        friction_stress_pa=1.0e6,
    )
    assert math.isinf(life)


def test_all_twelve_fcc_systems_are_retained_and_fastest_is_first():
    lives = fcc_axial_tmw_lives(0, 0, 1, 20.0, 69.0e9, 0.33)
    assert len(lives) == 12
    assert lives[0].cycles_to_initiation <= lives[-1].cycles_to_initiation
    assert np.isclose(lives[0].system.schmid_factor, 1.0 / np.sqrt(6.0))


def test_high_purity_al_reference_scale_is_not_an_sn_curve_fit():
    """At 4 MPa RSS amplitude, physical defaults give the observed HCF scale."""
    material = AluminumSlipInitiationParameters()
    axial_amplitude = 4.0 / 0.4898979485566356
    life = axial_tmw_initiation_life(
        2, 5, -1, axial_amplitude, 69.0e9, 0.33, material
    )
    assert np.isclose(life.system.schmid_factor, 0.4898979485566356)
    assert 4.0e6 < life.cycles_to_initiation < 6.0e6


def test_empirical_first_passage_shape_preserves_probability_and_censoring():
    shape = empirical_first_passage_shape(
        np.array([1.0, 2.0, 2.0, 4.0, np.nan]), observation_end_cycles=5.0
    )
    np.testing.assert_allclose(shape.probability_mass, [0.2, 0.4, 0.2])
    np.testing.assert_allclose(shape.cumulative_probability, [0.2, 0.6, 0.8])
    assert shape.quantile_multiplier(0.5) == 1.0
    assert np.isnan(shape.quantile_multiplier(0.9))
    assert np.isclose(shape.censored_probability, 0.2)


@pytest.mark.parametrize("field", ["lattice_parameter_nm", "surface_energy_j_m2"])
def test_positive_material_inputs_are_required(field):
    values = {field: 0.0}
    with pytest.raises(ValueError):
        AluminumSlipInitiationParameters(**values).validate()
