"""Mode-resolved damped-correlation inference from reference MD, not PDE time.

The equilibrium position correlation solves c''+2*g*c'+w0^2*c=0,
c(0)=1,c'(0)=0. With wd^2=w0^2-g^2 its underdamped solution is
exp(-g*t)*(cos(wd*t)+g/wd*sin(wd*t)). Integral c dt=2*g/w0^2.
Fitted source timestamps supply units, never an invented atomic-mass clock.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


def damped_position_correlation(time, damping, damped_angular_frequency):
    t=np.asarray(time,float);g=float(damping);w=float(damped_angular_frequency)
    if (np.any(~np.isfinite(t)) or np.any(t<0) or not np.isfinite(g+w)
            or g<0 or w<=0):
        raise ValueError('nonnegative time/damping and positive frequency required')
    return np.exp(-g*t)*(np.cos(w*t)+g*t*np.sinc(w*t/np.pi))


def mode_correlations(coordinates_m, max_lag):
    """Unitary spatial DFT; real ordered modal autocovariance [m^2].

    Conjugate modes are deliberately NOT counted as independent samples.
    FFT time zero-padding computes the unbiased finite-origin lag average.
    """
    q=np.asarray(coordinates_m,float)
    if (q.ndim!=3 or q.shape[2]!=3 or q.shape[1]<3 or np.any(~np.isfinite(q))
            or not isinstance(max_lag,int) or not 1<=max_lag<len(q)):
        raise ValueError('finite time/plane/3-coordinate array and integer lag required')
    z=np.fft.fft(q-q.mean(axis=0,keepdims=True),axis=1,norm='ortho')
    length=1 << (2*len(q)-1).bit_length()
    spectrum=np.fft.fft(z,n=length,axis=0)
    covariance=np.fft.ifft(abs(spectrum)**2,axis=0)[:max_lag+1].real
    covariance/=np.arange(len(q),len(q)-max_lag-1,-1)[:,None,None]
    return covariance


def fit_damped_correlation(time_ps, measured, frequency_seed_rad_ps):
    """Deterministic multi-start two-parameter DHO, plus overdamped control.

    Equal normalized-correlation residuals, not independent Gaussian errors.
    Jacobian singular values are numerical sensitivity, not confidence limits.
    Every start and bound is returned. A successful optimizer is NOT a
    successful physical calibration. Held-out/window tests belong to caller.
    """
    t=np.asarray(time_ps,float);y=np.asarray(measured,float)
    if (t.ndim!=1 or y.shape!=t.shape or len(t)<8 or t[0]!=0
            or np.any(~np.isfinite(t+y)) or np.any(np.diff(t)<=0)
            or not np.isfinite(frequency_seed_rad_ps) or frequency_seed_rad_ps<=0):
        raise ValueError('finite ordered time/correlation and positive spectral seed required')
    nyquist=np.pi/np.min(np.diff(t))
    low=np.log([1e-6, .02]);high=np.log([20.,.98*nyquist])
    starts=[];best=None
    for damping in (.01,.1,1.):
        for factor in (.9,1.,1.1):
            seed=np.array([damping,frequency_seed_rad_ps*factor])
            if not np.all((np.log(seed)>low)&(np.log(seed)<high)):
                continue
            fit=least_squares(lambda x:damped_position_correlation(t,*np.exp(x))-y,
                np.log(seed),bounds=(low,high),ftol=1e-11,xtol=1e-11,gtol=1e-11,
                max_nfev=600)
            record=dict(seed=seed,parameters=np.exp(fit.x),cost=float(2*fit.cost),
                success=bool(fit.success),status=int(fit.status),nfev=int(fit.nfev),
                active_mask=fit.active_mask)
            starts.append(record)
            if best is None or fit.cost<best.cost:best=fit
    if best is None:raise ValueError('no admissible spectral initialization')
    g,w=np.exp(best.x);w0=np.hypot(g,w)
    # Same normalized lag data, no oscillatory sign repair in the control.
    control=least_squares(lambda x:np.exp(-np.exp(x[0])*t)-y,
                         [0.],bounds=([-14.],[9.]),max_nfev=300)
    return dict(damping_per_ps=float(g),damped_angular_frequency_rad_ps=float(w),
        natural_angular_frequency_rad_ps=float(w0),envelope_time_ps=float(1/g),
        formal_integral_time_ps=float(2*g/w0**2),
        training_rmse=float(np.sqrt(2*best.cost/len(t))),
        prediction=damped_position_correlation(t,g,w),
        overdamped_rate_per_ps=float(np.exp(control.x[0])),
        overdamped_training_rmse=float(np.sqrt(2*control.cost/len(t))),
        singular_values_log_parameters=np.linalg.svd(best.jac,compute_uv=False),
        optimizer_success=bool(best.success),active_mask=best.active_mask,
        starts=starts,bounds_log_parameters=[low,high],
        production_calibration_available=False)


def formal_modal_mobility(variance_m2, integral_time_seconds, temperature_K):
    """Conditional modal mobility, not per-atom or a/s mobility.

    Coordinate normalization determines the conjugate energy. For a complex
    pair, caller must use a real normalized cosine/sine coordinate variance;
    E|z|^2 equals each sqrt(2)*Re/Im variance ONLY under phase equivalence.
    """
    values=np.asarray([variance_m2,integral_time_seconds,temperature_K],float)
    if np.any(~np.isfinite(values)) or np.min(values)<=0:
        raise ValueError('positive variance, resolved integral time and temperature required')
    return float(variance_m2/(1.380649e-23*temperature_K*integral_time_seconds))


def periodic_difference_symbol(planes):
    """Exact eigenvalues of B B^T for q_l=u_(l+1)-u_l, excluding k=0."""
    if int(planes)!=planes or planes<3:
        raise ValueError('integer periodic plane count >=3 required')
    return 4*np.sin(np.pi*np.arange(1,planes)/planes)**2


def gradient_drag_hypothesis(mode_indices, damping_per_ps, *, planes, fit_modes):
    """Test g_k=D |B_k|^2, not fit an arbitrary production mobility.

    This follows from a local Rayleigh function in relative plane velocities,
    R=(Np eta/2) sum_l |qdot_l|^2, with inertia Np m for u_l in reference MD.
    It conserves uniform-translation momentum; local drag on u_l would not.
    D=eta/(2m) has units inverse time. No atomic mass is supplied or inferred
    here. Fitting modes and held-out modes must be explicitly disjoint.
    """
    k=np.asarray(mode_indices);g=np.asarray(damping_per_ps,float)
    symbol=periodic_difference_symbol(planes)
    if (k.ndim!=1 or g.shape!=k.shape or np.any(k!=np.rint(k))
            or np.any(k<1) or np.any(k>=planes) or len(np.unique(k))!=len(k)
            or np.any(~np.isfinite(g)) or np.any(g<=0)):
        raise ValueError('unique nonzero modes and positive measured damping required')
    selected=np.isin(k,np.asarray(fit_modes))
    if not np.any(selected) or np.all(selected):
        raise ValueError('nonempty distinct fit and held-out modes required')
    w=symbol[k.astype(int)-1]
    coefficient=float(w[selected]@g[selected]/(w[selected]@w[selected]))
    predicted=coefficient*w
    return dict(coefficient_per_ps=coefficient,predicted_damping_per_ps=predicted,
        fit_mask=selected,per_mode_coefficient_per_ps=g/w,
        heldout_rmse_per_ps=float(np.sqrt(np.mean((predicted[~selected]-g[~selected])**2))),
        production_calibration_available=False)
