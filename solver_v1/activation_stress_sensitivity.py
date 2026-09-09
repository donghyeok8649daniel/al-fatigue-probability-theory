"""Static barrier traction derivative; not an atomistic activation-volume fit.

For stationary q_m and q_dagger in G=W-A_atomic_cell L0 T.q/eV_J,
-d(Delta G_phys)/dT = A_atomic_cell L0 (q_dagger-q_m) [m^3].
This is an exact envelope identity on smooth stationary branches. It does not
make an infinite-plane energy per cell a validated localized activated event.
"""
import numpy as np

from .aluminum_calibration import EV_J
from .vector_registry_audit import stationary_state


def barrier_and_derivative(model, minimum_initial, saddle_initial, *, traction_MPa, units):
    traction=np.asarray(traction_MPa,float)
    if traction.shape!=(3,) or np.any(~np.isfinite(traction)):
        raise ValueError('finite local normal/two-shear tractions required')
    force=units.traction_mpa_to_force(traction)
    minimum=stationary_state(model,minimum_initial,force=force,expected_index=0)
    saddle=stationary_state(model,saddle_initial,force=force,expected_index=1)
    if not minimum['valid'] or not saddle['valid']:
        raise ArithmeticError('stationary branch not verified; no barrier derivative reported')
    jump=saddle['q']-minimum['q']
    energy=saddle['evaluation'].energy-minimum['evaluation'].energy-force@jump
    if energy<=0:
        raise ArithmeticError('positive forward barrier required; branch identity must be audited')
    volume=units.atomic_cell_area_m2*units.length_scale_m*jump
    return dict(minimum=minimum,saddle=saddle,barrier_eV_cell=float(energy),
        minus_barrier_traction_derivative_m3=volume,
        minus_barrier_traction_derivative_eV_per_MPa=volume*1e6/EV_J,
        physical_local_activation_event_validated=False,
        kinetic_prefactor_available=False,cell_mobility=None)


def apparent_rate_slope_from_barrier(volume_m3, temperature_K):
    """Barrier-only d ln(rate)/dT [1/MPa] under CONSTANT prefactor hypothesis.

    Does not create a rate or use it in the PDE. Experimental apparent slopes
    also include prefactor/source-population stress dependence.
    """
    if not np.isfinite(temperature_K) or temperature_K<=0:
        raise ValueError('positive explicit temperature required')
    volume=np.asarray(volume_m3,float)
    if np.any(~np.isfinite(volume)):raise ValueError('finite signed derivative required')
    return volume/(1.380649e-23*temperature_K)*1e6
