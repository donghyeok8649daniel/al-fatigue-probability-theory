"""Finite-temperature constrained Hamiltonian reference, never a fitted clock."""
from __future__ import annotations

import numpy as np

from .silicon_conditional_research import AMU_EV_PS2_A2, SI_MASS_AMU


def reflecting_drift(position, velocity, duration, halfwidth):
    """Exact free motion and elastic reflections in an orthogonal fixed box.

    Returns crossed-wall count, including multiple crossings within one drift.
    Equality with a wall has measure zero; the right-cell convention sets it.
    """
    x,v=np.asarray(position,float),np.asarray(velocity,float)
    if (x.shape!=v.shape or x.ndim!=1 or not np.isfinite(x).all() or not np.isfinite(v).all()
            or not np.isfinite(duration) or duration<=0 or not np.isfinite(halfwidth)
            or halfwidth<=0 or np.any(abs(x)>halfwidth)):
        raise ValueError('finite vectors in a positive orthogonal box required')
    raw=x+duration*v
    cell=np.floor((raw+halfwidth)/(2*halfwidth)).astype(np.int64)
    folded=(raw+halfwidth)%(4*halfwidth)
    result=halfwidth-abs(folded-2*halfwidth)
    new_velocity=np.where(cell%2==0,v,-v)
    return result,new_velocity,int(np.sum(abs(cell)))


def constrained_verlet(model, coordinates, initial, velocities, *, gap, dt_ps, steps,
                       halfwidth, stride=1, mass_amu=SI_MASS_AMU, progress=None):
    """Physical atom-mass dynamics on the orthonormal fixed-gap bath.

    Orthogonality gives kinetic energy m |xdot|^2/2. Reflections, if reached,
    implement the SAME declared canonical box; q and all grip atoms are exact.
    No thermostat is used. Correlation of constraint forces is a nonlinear
    frozen-q diagnostic, not automatically the exact Mori-Zwanzig memory.
    """
    x,v=np.asarray(initial,float).copy(),np.asarray(velocities,float).copy()
    if (x.shape!=(coordinates.dimension,) or v.shape!=x.shape
            or not np.isfinite(x).all() or not np.isfinite(v).all()
            or not np.isfinite(dt_ps) or dt_ps<=0 or int(steps)!=steps or steps<1
            or int(stride)!=stride or stride<1 or steps%stride
            or not np.isfinite(mass_amu) or mass_amu<=0):
        raise ValueError('valid orthonormal initial state and integral time controls required')
    mass=mass_amu*AMU_EV_PS2_A2
    result=model.evaluate(coordinates.positions(x,[gap]))
    gradient,reaction=coordinates.pullback(result.gradient)
    sites0=result.site_energy.copy()
    kinetic0=.5*mass*(v@v)
    times,observations,potential,kinetic,residual=[],[],[],[],[]
    max_error,reflections=0.,0
    for step in range(int(steps)+1):
        if step:
            v-=.5*dt_ps*gradient/mass
            x,v,count=reflecting_drift(x,v,dt_ps,halfwidth)
            reflections+=count
            result=model.evaluate(coordinates.positions(x,[gap]))
            gradient,reaction=coordinates.pullback(result.gradient)
            v-=.5*dt_ps*gradient/mass
        delta_u=float(np.sum(result.site_energy-sites0))
        k=float(.5*mass*(v@v))
        error=delta_u+k-kinetic0
        if not np.isfinite(error):
            raise FloatingPointError('nonfinite constrained reference dynamics')
        max_error=max(max_error,abs(error))
        if step%stride==0:
            times.append(step*dt_ps)
            observations.append([float(reaction[0]),float(np.max(abs(x))),float(np.mean(x*x))])
            potential.append(delta_u);kinetic.append(k);residual.append(error)
        if progress is not None and step and step%10000==0:
            progress(step)
    positions=coordinates.positions(x,[gap])
    return dict(time_ps=np.asarray(times),observations=np.asarray(observations),
        potential_difference_eV=np.asarray(potential),kinetic_energy_eV=np.asarray(kinetic),
        energy_residual_eV=np.asarray(residual),max_energy_residual_eV=max_error,
        initial_kinetic_energy_eV=float(kinetic0),reflection_count=reflections,
        final_bath_coordinates_A=x,final_bath_velocities_A_ps=v,final_positions_A=positions,
        dt_ps=dt_ps,steps=steps,stride=stride,mass_amu=mass_amu,
        bath_box_halfwidth_A=halfwidth,production_t0_seconds=None)


def stationary_correlation(history, maximum_lag, *, center=None):
    """Time-origin average with the actual N-lag denominator, using an FFT.

    A supplied independent ensemble mean avoids subtracting the trajectory's
    own zero-frequency component. Both center choices should be audited.
    Overlapping origins are NOT treated as independent replicas.
    """
    values=np.asarray(history,float)
    if (values.ndim!=1 or not np.isfinite(values).all() or len(values)<2
            or int(maximum_lag)!=maximum_lag or not 0<=maximum_lag<len(values)):
        raise ValueError('finite scalar history and represented lags required')
    mean=float(values.mean()) if center is None else float(center)
    if not np.isfinite(mean):
        raise ValueError('finite center required')
    delta=values-mean
    size=1<<(2*len(values)-1).bit_length()
    spectrum=np.fft.rfft(delta,n=size)
    raw=np.fft.irfft(spectrum.conj()*spectrum,n=size)[:maximum_lag+1]
    return raw/(len(values)-np.arange(maximum_lag+1))
