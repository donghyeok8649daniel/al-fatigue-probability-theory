import numpy as np
import pytest

from solver_v1.configurational_plasticity import (
    crystallographic_slip_projection,
    registry_decomposition,
    registry_rate_terms,
)


def test_registry_decomposition_reconstructs_s_in_half_open_wells() -> None:
    b = 1.7
    s = np.array([-3.4, -2.55, -1.7, -0.85, 0.0, 0.85, 1.7, 2.55])
    n, xi = registry_decomposition(s, b)
    np.testing.assert_allclose(s, b * n + xi, rtol=0.0, atol=2.0e-16)
    assert np.all(xi >= -0.5 * b)
    assert np.all(xi < 0.5 * b)


def test_well_assignment_at_exact_boundaries_uses_lower_closed_upper_open() -> None:
    b = 2.0
    boundaries = np.array([-3.0, -1.0, 1.0, 3.0])
    n, xi = registry_decomposition(boundaries, b)
    np.testing.assert_array_equal(n, np.array([-1, 0, 1, 2]))
    np.testing.assert_allclose(xi, -0.5 * b, rtol=0.0, atol=0.0)


def test_registry_rate_separates_interwell_flow_from_selective_absorption() -> None:
    terms = registry_rate_terms(
        well_indices=[-1, 0, 1],
        well_populations=[0.1, 0.7, 0.1],
        interwell_net_flux=[0.02, -0.01],
        interwell_gross_flux=[0.08, 0.07],
        opening_absorption_rate_by_well=[0.03, 0.0, 0.0],
    )
    assert np.isclose(terms.net_registry_flow_rate, 0.01)
    assert np.isclose(terms.gross_configurational_activity_rate, 0.15)
    assert np.isclose(terms.opening_registry_loss_rate, -0.03)
    assert terms.gross_configurational_activity_rate > abs(
        terms.net_registry_flow_rate
    )
    # Selectively removing n=-1 increases the conditional surviving mean, but
    # that contribution is not interwell plastic flow.
    assert terms.conditional_selective_opening_rate > 0.0
    assert np.isclose(
        terms.conditional_net_registry_flow_rate,
        terms.net_registry_flow_rate / 0.9,
    )


def test_symmetric_bidirectional_hopping_can_have_zero_net_registry_flow() -> None:
    terms = registry_rate_terms(
        well_indices=[-1, 0, 1],
        well_populations=[0.1, 0.8, 0.1],
        interwell_net_flux=[0.0, 0.0],
        interwell_gross_flux=[0.2, 0.2],
        opening_absorption_rate_by_well=[0.0, 0.0, 0.0],
    )
    assert terms.net_registry_flow_rate == 0.0
    assert terms.gross_configurational_activity_rate == 0.4
    assert terms.conditional_mean_well_index == 0.0


def test_crystallographic_projection_matches_single_slip_tensor_identity() -> None:
    e = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    m = np.array([1.0, 0.0, 0.0])
    normal = np.array([0.0, 1.0, 0.0])
    result = crystallographic_slip_projection(
        e, m, normal, a0=0.8, slip_plane_spacing=0.4
    )
    assert np.isclose(result.axial_shear_projection, 0.5)
    assert np.isclose(result.chi_candidate, 1.0)


def test_crystallographic_projection_rejects_nonunit_or_invalid_slip_geometry() -> None:
    with pytest.raises(ValueError, match="normalized"):
        crystallographic_slip_projection(
            [2.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            a0=1.0,
            slip_plane_spacing=1.0,
        )
    with pytest.raises(ValueError, match="slip plane"):
        crystallographic_slip_projection(
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            a0=1.0,
            slip_plane_spacing=1.0,
        )
