"""Physical dislocation-line kinetics and a CONDITIONAL metric derivation.

No source coefficient is silently inserted into the a/s PDE or static energy.
B_line v=tau b has B_line [Pa s], not the cell friction [J s/m^2].
"""
import json
from pathlib import Path

import numpy as np


REFERENCE_FILE=Path(__file__).with_name('data')/'aluminum_line_drag_references.json'


def load_line_drag_references():
    data=json.loads(REFERENCE_FILE.read_bytes())
    rows=data['entries'];ids=[r['id'] for r in rows]
    if len(ids)!=len(set(ids)) or data['production_time_calibration_available']:
        raise ValueError('source-only unique line references required')
    return {r['id']:r for r in rows}


def source_drag_pa_s(reference_id, temperature_K):
    """Reproduce the published coefficient only within its declared conditions."""
    source=load_line_drag_references()[reference_id]
    T=float(temperature_K)
    if not np.isfinite(T) or T<=0:raise ValueError('physical positive temperature required')
    if source['quantity']=='B_line/T':
        if not source['temperature_min_K']<=T<=source['temperature_max_K']:
            raise ValueError('outside published temperature span; no automatic extrapolation')
        return float(T*source['value_SI'])
    if T!=source['temperature_K']:
        raise ValueError('this is a measured endpoint, not an invented interpolation law')
    return float(source['value_SI'])


def source_line_velocity_over_b(reference_id, *, shear_mpa, temperature_K):
    """v/b [1/s] of the published low-speed fit, NOT solver frequency in Hz.

    Avoid silently mixing a 0K lattice length with a finite-T source Burgers
    vector. Multiplying by a condition-matched b requires separate geometry.
    """
    source=load_line_drag_references()[reference_id]
    stress=np.asarray(shear_mpa,float)
    if np.any(~np.isfinite(stress)):raise ValueError('finite physical shear required')
    bound=source.get('maximum_abs_shear_over_T_MPa_K')
    if bound is None:raise ValueError('published fit range required; endpoint drag alone is not a velocity law')
    drag=source_drag_pa_s(reference_id,temperature_K)
    if np.max(abs(stress),initial=0)/temperature_K>bound:
        raise ValueError('outside the source low-velocity fit range')
    return stress*1e6/drag


def registry_profile_metric(coordinate_m, slip_m):
    """Integral |ds/dx|^2 dx [m] for the explicit piecewise-linear profile.

    This is a geometric friction metric. Refining a continuum reconstruction
    is not proof of atomistic-core convergence. A constant offset is irrelevant.
    """
    x=np.asarray(coordinate_m,float);s=np.asarray(slip_m,float)
    if x.ndim!=1 or x.shape!=s.shape or len(x)<2 or np.any(~np.isfinite(x+s)) or np.any(np.diff(x)<=0):
        raise ValueError('finite ordered physical x and matching slip required')
    return float(np.sum(np.diff(s)**2/np.diff(x)))


def conditional_registry_metric(*, line_drag_pa_s, profile_metric_m, atomic_cell_area_m2):
    r"""Derived local-friction hypothesis, explicitly NOT a calibration.

    Per-line Rayleigh R/ell=1/2 int eta_s (ds/dt)^2 dx. Rigid translation
    s(x-X(t)) gives B_line=eta_s int(s')^2 dx. A cell then has
    zeta_cell=eta_s A_atomic_cell, M_cell=1/zeta_cell.

    Thus M_cell=I/(B_line A_atomic_cell) [m^2/(J s)], BUT only if uniform
    local slip friction adequately represents the measured line dissipation
    and the coordinate/core/temperature/material all match. Far-field phonon
    drag is not automatically local core dissipation. No t0 is returned.
    """
    values=np.array([line_drag_pa_s,profile_metric_m,atomic_cell_area_m2],float)
    if np.any(~np.isfinite(values)) or np.any(values<=0):
        raise ValueError('positive finite drag, measured/derived profile metric and atomic cell area required')
    drag,metric,area=values
    eta=drag/metric
    return dict(slip_friction_per_area_J_s_m4=float(eta),
        conditional_cell_mobility_m2_J_s=float(metric/(drag*area)),
        assumptions_validated=False,production_calibration_available=False,
        normal_mobility_available=False,t0_seconds=None,
        interpretation='conditional uniform-local-friction metric, NOT accepted Al a/s mobility')
