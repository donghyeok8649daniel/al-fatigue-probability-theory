"""Exact linear-coordinate/energy changes; no local-material clock approval."""
import numpy as np
from scipy.linalg import helmert


def plane_zero_sum_basis(planes,components=3):
    if int(planes)!=planes or planes<2 or int(components)!=components or components<1:
        raise ValueError('integer dimensions required')
    return np.kron(helmert(planes,full=False).T,np.eye(components))


def harmonic_generator(H,M,thermal_energy):
    H,M=np.asarray(H,float),np.asarray(M,float)
    if (H.ndim!=2 or H.shape[0]!=H.shape[1] or M.shape!=H.shape
        or not np.all(np.isfinite(H+M)) or np.linalg.norm(H-H.T)>1e-12*np.linalg.norm(H)
        or np.linalg.norm(M-M.T)>1e-12*np.linalg.norm(M) or np.linalg.eigvalsh(H).min()<=0
        or np.linalg.eigvalsh(M).min()<=0 or not np.isfinite(thermal_energy) or thermal_energy<=0):
        raise ValueError('stable symmetric Hessian/mobility and positive thermal energy required')
    return M@H,thermal_energy*M


def covariance_mobility(C,K,kBT):
    """Full matched-coordinate low-band proxy, not scalar inverse susceptibility."""
    C,K=np.asarray(C,float),np.asarray(K,float)
    harmonic_generator(C,K,kBT)  # same SPD requirements, no unit reinterpretation
    if np.linalg.matrix_rank(C)<len(C) or np.linalg.matrix_rank(K)<len(K):
        raise ValueError('numerically rank-deficient covariance/integral; no inverse mobility')
    # Declared 1% numerical inversion-error budget, NOT a physical noise floor.
    gamma=len(K)*np.finfo(float).eps/(1-len(K)*np.finfo(float).eps)
    if max(np.linalg.cond(C),np.linalg.cond(K))*gamma>.01:
        raise ValueError('condition*gamma_n exceeds 1% roundoff budget; no inverse mobility')
    H=kBT*np.linalg.inv(C)
    drag=kBT*np.linalg.solve(C,np.linalg.solve(C,K).T).T
    mobility=np.linalg.inv(drag)
    return H,(mobility+mobility.T)/2


def energy_rescale(H,M,kBT,sites):
    if not np.isfinite(sites) or sites<=0:raise ValueError('positive normalization required')
    harmonic_generator(H,M,kBT)
    return np.asarray(H)/sites,np.asarray(M)*sites,kBT/sites


def projected_susceptibility(H,M,projection,omega):
    """Exact harmonic response to conjugate force, exp(+i omega t) convention.

    This research diagnostic does not assert a Markovian scalar reduction.
    H, M and omega must use one consistent coordinate/energy/time convention.
    """
    H,M=np.asarray(H,float),np.asarray(M,float)
    drift,_=harmonic_generator(H,M,1.)
    v=np.asarray(projection,float)
    if v.shape!=(len(H),) or not np.all(np.isfinite(v)) or not np.isfinite(omega) or omega<0:
        raise ValueError('finite projection and nonnegative angular frequency required')
    return complex(v@np.linalg.solve(drift+1j*omega*np.eye(len(H)),M@v))
