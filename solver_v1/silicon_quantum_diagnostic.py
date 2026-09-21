"""Local harmonic nuclear correction diagnostic, not a Si quantum PMF model."""
from __future__ import annotations

import numpy as np
from scipy.constants import hbar,electron_volt

from .silicon_conditional_research import KB_EV_K,AMU_EV_PS2_A2,SI_MASS_AMU


HBAR_EV_PS=hbar/electron_volt/1e-12


def quantum_harmonic_difference(values,reference_values,temperature_K,*,mass_amu=SI_MASS_AMU):
    """Difference of same-dimensional stable oscillator free energies in eV.

    These are conditional bath oscillators at fixed classical q. Their sum
    does not include q tunneling, anharmonic nuclear motion or electronic DFT.
    Applying it to an empirical potential is a sensitivity diagnostic only.
    """
    a,b=np.asarray(values,float),np.asarray(reference_values,float)
    if (a.ndim!=1 or a.shape!=b.shape or not len(a) or np.any(a<=0) or np.any(b<=0)
            or not np.isfinite(a).all() or not np.isfinite(b).all()
            or not np.isfinite(temperature_K) or temperature_K<0
            or not np.isfinite(mass_amu) or mass_amu<=0):
        raise ValueError('equal positive bath spectra and nonnegative temperature required')
    energies=HBAR_EV_PS*np.sqrt(a/(mass_amu*AMU_EV_PS2_A2))
    reference=HBAR_EV_PS*np.sqrt(b/(mass_amu*AMU_EV_PS2_A2))
    zero_point=float(.5*np.sum(energies-reference))
    if temperature_K==0:
        return zero_point
    thermal=KB_EV_K*temperature_K
    return zero_point+float(thermal*np.sum(np.log(-np.expm1(-energies/thermal))
                                          -np.log(-np.expm1(-reference/thermal))))
