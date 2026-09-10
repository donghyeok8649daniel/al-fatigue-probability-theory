"""Signed slip-area -> specimen strain, NOT a new plastic-flow/material law.

The caller supplies geometrically resolved swept surfaces from a deterministic
defect calculation or observation. No emission rate, source density, activation
length, specimen correlation area, or physical clock is inferred here.
Recoverable bowing contributes slip strain too; persistence is a separate test.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from .finite_source_reference import pinned_source_branch


def _geometry(burgers_vectors_m, plane_normals, specimen_volume_m3):
    b = np.asarray(burgers_vectors_m, float)
    n = np.asarray(plane_normals, float)
    if (b.ndim != 2 or b.shape[1] != 3 or len(b) == 0 or n.shape != b.shape
            or np.any(~np.isfinite([b, n]))
            or not np.isfinite(specimen_volume_m3) or specimen_volume_m3 <= 0):
        raise ValueError('finite (systems,3) Burgers vectors/normals and positive specimen volume required')
    lengths = np.linalg.norm(b, axis=1)
    if (np.any(lengths <= 0) or not np.allclose(np.linalg.norm(n, axis=1), 1., atol=1e-12, rtol=0)
            or np.any(abs(np.einsum('ki,ki->k', b/lengths[:, None], n)) > 1e-12)):
        raise ValueError('nonzero tangential Burgers vectors and unit plane normals required')
    return b, n


def slip_distortion(signed_area_m2, *, burgers_vectors_m, plane_normals, specimen_volume_m3):
    r"""beta_slip = sum_k A_k b_k tensor n_k / V, with exact SI units.

    A has shape (..., systems) and can be negative. Positive A and b encode
    the declared slip sense. At small strain epsilon_slip=sym(beta_slip).
    All surfaces belong to the SAME explicitly defined averaging volume.
    This identity does not assert that their motions are independent.
    """
    b, n = _geometry(burgers_vectors_m, plane_normals, specimen_volume_m3)
    area = np.asarray(signed_area_m2, float)
    if area.ndim < 1 or area.shape[-1] != len(b) or np.any(~np.isfinite(area)):
        raise ValueError('finite signed areas with one final-axis entry per slip system required')
    beta = np.einsum('...k,ki,kj->...ij', area, b, n)/specimen_volume_m3
    if np.any(~np.isfinite(beta)):
        raise ValueError('slip distortion exceeds finite numerical range')
    return beta


def symmetric_slip_strain(signed_area_m2, **geometry):
    beta = slip_distortion(signed_area_m2, **geometry)
    return (beta + np.swapaxes(beta, -1, -2))/2


def axial_slip_strain(strain_tensor, loading_direction):
    e = np.asarray(loading_direction, float)
    strain = np.asarray(strain_tensor, float)
    if (e.shape != (3,) or np.any(~np.isfinite(e))
            or not np.isclose(np.linalg.norm(e), 1., rtol=0, atol=1e-12)
            or strain.shape[-2:] != (3, 3) or np.any(~np.isfinite(strain))
            or not np.allclose(strain, np.swapaxes(strain, -1, -2), rtol=1e-12, atol=0)):
        raise ValueError('unit loading direction and symmetric finite strain tensors required')
    return np.einsum('i,...ij,j->...', e, strain, e)


def transport_history(signed_area_history_m2, **geometry):
    """Changes relative to initial slip state; separate net and gross motion.

    Forward/backward are positive/negative changes of RESOLVED swept area,
    not abs(net) and not nearest-neighbor SG recrossing counts. With sampled
    states they omit unobserved reversals between saves; time refinement is
    required. Gross area is never converted into signed plastic strain.
    """
    area = np.asarray(signed_area_history_m2, float)
    if area.ndim != 2 or len(area) < 2 or np.any(~np.isfinite(area)):
        raise ValueError('at least two finite area snapshots of shape (steps, systems) required')
    strain = symmetric_slip_strain(area - area[0], **geometry)
    increments = np.diff(area, axis=0)
    forward = np.vstack([np.zeros(area.shape[1]), np.cumsum(np.maximum(increments, 0), axis=0)])
    backward = np.vstack([np.zeros(area.shape[1]), np.cumsum(np.maximum(-increments, 0), axis=0)])
    return dict(relative_signed_area_m2=area-area[0], forward_area_m2=forward,
                backward_area_m2=backward, gross_area_m2=forward+backward,
                signed_slip_strain_tensor=strain,
                area_balance_residual_m2=forward-backward-(area-area[0]),
                residual_plasticity_certified=False)


def uniform_stress_work(stress_midpoint_pa, area_increments_m2, **geometry):
    """Independent tensor and traction evaluations of virtual work [J].

    For uniform symmetric sigma, V sigma:Delta epsilon =
    sum_k (b_k . sigma . n_k) Delta A_k. Midpoint sampling is a work quadrature,
    not an exact time integral when stress changes nonlinearly within a step.
    Nonuniform stress needs the local surface integral, not this shortcut.
    """
    b, n = _geometry(**geometry)
    area = np.asarray(area_increments_m2, float)
    strain = symmetric_slip_strain(area, **geometry)
    stress = np.asarray(stress_midpoint_pa, float)
    if (stress.shape[-2:] != (3, 3) or np.any(~np.isfinite(stress))
            or not np.allclose(stress, np.swapaxes(stress, -1, -2), atol=0, rtol=1e-12)):
        raise ValueError('finite symmetric uniform stress tensor required')
    tensor = geometry['specimen_volume_m3']*np.sum(stress*strain, axis=(-2, -1))
    traction = np.sum(np.einsum('ki,...ij,kj->...k', b, stress, n)*area, axis=-1)
    return dict(tensor_work_J=tensor, slip_surface_work_J=traction, residual_J=tensor-traction)


def pinned_swept_area(line, endpoint_angle, *, span_m, outer_log_ratio, relative_tolerance=1e-9):
    r"""Minor-arc area 2 int y(theta) Q'(theta)/p dtheta, in square metres.

    Adaptive quadrature of the analytic angular primitives avoids mistaking
    the coarse plotted polygon for a resolved swept area near the fold.
    The energy/core restrictions of pinned_source_branch remain unchanged.
    """
    if not np.isfinite(relative_tolerance) or not 1e-13 <= relative_tolerance <= 1e-3:
        raise ValueError('explicit relative quadrature tolerance in [1e-13,1e-3] required')
    branch = pinned_source_branch(line, endpoint_angle, span_m=span_m,
                                  outer_log_ratio=outer_log_ratio, points=17)
    gamma = lambda t: float(line.evaluate(t)*outer_log_ratio)
    derivative = lambda t: float(line.evaluate(t, 1)*outer_log_ratio)
    primitive = lambda t: -gamma(t)*np.cos(t)+derivative(t)*np.sin(t)
    p = branch['applied_shear_Pa']*np.linalg.norm(line.burgers_m)
    top = primitive(endpoint_angle)
    def integrand(t):
        y = (top-primitive(t))/p/span_m
        dx = float(line.stiffness(t)*outer_log_ratio)*np.cos(t)/p/span_m
        return 2*y*dx
    normalized, error = quad(integrand, 0., endpoint_angle, epsabs=0., epsrel=relative_tolerance)
    if normalized <= 0 or error > relative_tolerance*normalized:
        raise ArithmeticError('positive swept area not quadrature-resolved')
    return dict(area_m2=normalized*span_m**2, area_over_span_squared=normalized,
                estimated_quadrature_error_m2=error*span_m**2,
                applied_shear_Pa=branch['applied_shear_Pa'],
                critical_shear_outer_only_Pa=branch['critical_shear_outer_only_Pa'])


def pinned_population_strain_budget(*, area_per_source_m2, span_m, burgers_m,
                                    schmid_factor, axial_criterion):
    r"""Necessary swept-area budget, not a fitted source population or yield law.

    Identical same-sense sources give epsilon=m b rho_line A/L, where rho_line
    means INITIAL pinned length/V. eta=rho_line L^2=N_sources L^3/V is exposed
    as an overlap diagnostic. eta<<1 is a sufficient spacing heuristic for
    dilute uniformly distributed sources, not a universal interaction bound.
    Solving for rho is demand accounting; it must NOT become a fitted input.
    Recoverable bow area does not establish residual plasticity even if a
    loading-branch offset criterion is reached.
    """
    values = [area_per_source_m2, span_m, burgers_m, schmid_factor, axial_criterion]
    if np.any(~np.isfinite(values)) or min(values) <= 0 or schmid_factor > .5+1e-12:
        raise ValueError('positive finite geometry/criterion and Schmid factor in (0,1/2] required')
    required = axial_criterion*span_m/(schmid_factor*burgers_m*area_per_source_m2)
    return dict(required_initial_pinned_line_density_m2=required,
                required_source_number_density_m3=required/span_m,
                required_overlap_parameter=required*span_m**2,
                axial_strain_per_unit_overlap=schmid_factor*burgers_m*area_per_source_m2/span_m**3,
                population_calibrated=False, experimental_yield_MPa=None)
