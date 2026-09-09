"""Low-frequency correlation-integral audit; no production clock inference.

For a declared equilibrium harmonic projection, C0=kBT H^-1. A linear GLE
has H integral_0^inf C(t)dt = Gamma(0) C0, whether its short-time dynamics
is inertial or not. Thus Gamma(0)=kBT C0^-1 K C0^-1 and M0=Gamma(0)^-1.
This is a ZERO-FREQUENCY response, not proof of a Markov process at all lags.
No atomic mass is introduced into the production overdamped equation here.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .public_aluminum_kinetics import lag_covariances


def _covariance_roots(covariance):
    C=np.asarray(covariance,float)
    if (C.ndim!=2 or C.shape[0]!=C.shape[1] or np.any(~np.isfinite(C))
            or not np.allclose(C,C.T,rtol=1e-12,atol=0)):
        raise ValueError('finite symmetric covariance required')
    e,V=np.linalg.eigh(C)
    if np.min(e)<=0:raise ValueError('positive covariance required; no regularization')
    return (V*np.sqrt(e))@V.T,(V/np.sqrt(e))@V.T


def formal_zero_frequency_response(covariance_m2, integral_m2_seconds, thermal_energy_J):
    """Formal reciprocal matrices for an ALREADY justified complete integral.

    This algebraic helper does not certify data, a Markov window, coordinate
    equivalence, thermostat independence or any physical PDE time conversion.
    Finite-tail estimation and numerical positivity are audited separately.
    """
    root,white=_covariance_roots(covariance_m2)
    K=np.asarray(integral_m2_seconds,float)
    if (K.shape!=root.shape or np.any(~np.isfinite(K))
            or not np.allclose(K,K.T,rtol=1e-12,atol=0)
            or not np.isfinite(thermal_energy_J) or thermal_energy_J<=0):
        raise ValueError('matching symmetric integral and positive thermal energy required')
    time_matrix=white@K@white
    if np.min(np.linalg.eigvalsh(time_matrix))<=0:
        raise ValueError('positive converged integral required; do not clip a negative integral')
    friction=thermal_energy_J*white@time_matrix@white
    mobility=root@np.linalg.inv(time_matrix)@root/thermal_energy_J
    return dict(friction_J_seconds_per_m2=friction,mobility_m2_per_J_second=mobility,
        whitened_integral_seconds=time_matrix,zero_frequency_only=True,
        production_calibration_available=False)


def correlation_integral_audit(coordinates, *, frame_seconds, max_lag, blocks=4):
    """All lag cutoffs, time blocks and stride quadrature retained, no fit window.

    The maximum block discrepancy is a sensitivity, NOT a confidence interval.
    A late plateau with errors below the integral is necessary but not
    sufficient for material kinetics: this function cannot authorize seconds.
    """
    q=np.asarray(coordinates,float)
    if (q.ndim!=3 or np.any(~np.isfinite(q)) or int(blocks)!=blocks or blocks<2
            or int(max_lag)!=max_lag or max_lag<2
            or not np.isfinite(frame_seconds) or frame_seconds<=0):
        raise ValueError('finite trajectory and positive declared lag/time/block settings required')
    width=len(q)//blocks
    if width<2*(max_lag+1):raise ValueError('each block must contain at least two lag windows')
    lags=np.arange(max_lag+1);time=lags*frame_seconds
    C=lag_covariances(q,lags);_,white=_covariance_roots(C[0])
    Cblocks=np.array([lag_covariances(q[j*width:(j+1)*width],lags) for j in range(blocks)])
    K=cumulative_trapezoid(C,time,axis=0,initial=0)
    Kblocks=cumulative_trapezoid(Cblocks,time,axis=1,initial=0)
    normalized=white@K@white;normalized_blocks=white@Kblocks@white
    symmetric=(normalized+normalized.swapaxes(-1,-2))/2
    symmetric_blocks=(normalized_blocks+normalized_blocks.swapaxes(-1,-2))/2
    block_error=np.max(np.linalg.norm(symmetric_blocks-symmetric,ord=2,axis=(-2,-1)),axis=0)
    antisymmetry=np.linalg.norm(normalized-normalized.swapaxes(-1,-2),ord=2,axis=(-2,-1))/2
    # Same measured C, every-other-lag trapezoid: quadrature refinement only,
    # not a second independent trajectory or a justification for decimation.
    coarse=cumulative_trapezoid(C[::2],time[::2],axis=0,initial=0)
    change=white@(coarse-K[::2])@white
    quadrature_change=np.linalg.norm(change,ord=2,axis=(-2,-1))
    eigen=np.linalg.eigvalsh(symmetric)[::2]
    floor=block_error[::2]+antisymmetry[::2]+quadrature_change
    # Cutoff dependence is separately visible, never hidden by chosen windows.
    half_cutoff_change=np.array([np.linalg.norm(symmetric[i]-symmetric[i//2],2)
                                 for i in range(0,max_lag+1,2)])
    return dict(cutoff_seconds=time[::2],covariance_m2=C[0],
        integral_m2_seconds=K[::2],whitened_integral_seconds=symmetric[::2],
        integral_eigenvalues_seconds=eigen,block_sensitivity_seconds=block_error[::2],
        antisymmetric_seconds=antisymmetry[::2],quadrature_change_seconds=quadrature_change,
        empirical_floor_seconds=floor,half_cutoff_change_seconds=half_cutoff_change,
        positive_integral_above_floor=eigen[:,0]>floor,
        all_components_positive=np.all(eigen>0,axis=1),
        blocks=int(blocks),block_frames=width,unused_block_tail_frames=len(q)-width*blocks,
        production_calibration_available=False,statistical_confidence_claimed=False,
        interpretation='Finite correlation integrals; positive cutoffs alone are not a converged zero-frequency mobility.')
