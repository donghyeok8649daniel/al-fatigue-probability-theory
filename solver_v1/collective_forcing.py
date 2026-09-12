"""Conjugate forcing and linear response of a measured periodic plane gap.

Reference MD only: this is NOT a mobility conversion for the production PDE.
q=(mean u_upper-mean u_lower).e, G_external=-F(t)q.
"""
import numpy as np
from scipy.integrate import trapezoid


def conjugate_atom_weights(labels, lower=0, upper=1):
    labels=np.asarray(labels)
    if labels.ndim!=1 or lower==upper or np.any(labels!=np.rint(labels)):
        raise ValueError('distinct represented planes and integer labels required')
    counts=[np.count_nonzero(labels==p) for p in (lower,upper)]
    if min(counts)==0:raise ValueError('both planes must contain atoms')
    return (labels==upper)/counts[1]-(labels==lower)/counts[0]


def harmonic_response(time_ps, displacement_angstrom, frequency_per_ps,
                      signed_force_eV_A):
    """Fit offset + cos/sin; Re[chi F exp(i omega t)] convention.

    No detrending, amplitude repair, selected cycle or confidence claim.
    Caller compares temporal blocks and signed/amplitude controls separately.
    """
    t=np.asarray(time_ps,float);q=np.asarray(displacement_angstrom,float)
    if (t.ndim!=1 or q.shape!=t.shape or len(t)<8 or np.any(np.diff(t)<=0)
            or not np.all(np.isfinite(t+q)) or frequency_per_ps<=0
            or not np.isfinite(frequency_per_ps+signed_force_eV_A)
            or signed_force_eV_A==0):raise ValueError('finite ordered response required')
    X=np.stack([np.ones_like(t),np.cos(2*np.pi*frequency_per_ps*t),
                np.sin(2*np.pi*frequency_per_ps*t)],axis=1)
    coefficient,_,rank,_=np.linalg.lstsq(X,q,rcond=None)
    if rank!=3:raise ValueError('unresolved sinusoidal design')
    chi=(coefficient[1]-1j*coefficient[2])/signed_force_eV_A
    return dict(real_A2_eV=float(chi.real),imag_A2_eV=float(chi.imag),
                magnitude_A2_eV=float(abs(chi)),phase_rad=float(np.angle(chi)),
                offset_A=float(coefficient[0]),residual_rms_A=float(np.sqrt(np.mean((q-X@coefficient)**2))))


def susceptibility_from_covariance(time_seconds,covariance_m2,frequency_hz,kBT_J):
    """Canonical equilibrium FDT, not proof that finite NVE is canonical.

    With exp(+i omega t) loading, chi=(C0-i omega Ctilde)/kBT if C(inf)=0.
    Also return the endpoint term for the explicitly truncated -C' integral.
    Their difference and cutoff dependence must not be hidden.
    """
    t=np.asarray(time_seconds,float);c=np.asarray(covariance_m2,float)
    if (t.ndim!=1 or c.shape!=t.shape or len(t)<2 or t[0]!=0
            or np.any(np.diff(t)<=0) or not np.all(np.isfinite(t+c))
            or not np.isfinite(frequency_hz+kBT_J) or frequency_hz<0 or kBT_J<=0):
        raise ValueError('finite ordered covariance and physical units required')
    omega=2*np.pi*frequency_hz
    integral=trapezoid(c*np.exp(-1j*omega*t),t)
    chi=(c[0]-1j*omega*integral)/kBT_J
    endpoint=c[-1]*np.exp(-1j*omega*t[-1])/kBT_J
    return chi,chi-endpoint
