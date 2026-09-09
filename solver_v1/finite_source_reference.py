"""Long-wave finite pinned-line diagnostic from the SAME bulk elastic tensor.

This is NOT a new empirical plasticity law, an atomistic core solution, a
single-arm source model, an activation barrier, or a specimen yield predictor.
There is no default source length, core radius, kinetic scale, or yield factor.
The leading logarithmic energy coefficient comes from the two-half-space
Schur kernel already derived from the infinite LJ/Bessel energy.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.linalg import solve_banded

from .nonlocal_interface_elasticity import halfspace_impedance


@dataclass(frozen=True)
class ElasticLineCoefficient:
    """Pi-periodic k(theta) [J/m] in E_line=k ln(R/r_core)+E_core.

    Theta is the angle of the tangent measured from in-plane x. Burgers vector
    is explicit. No core term or length is inferred from this static Hessian.
    Angular derivatives are of a converged Fourier series, not energy FD.
    """
    coefficients: np.ndarray
    frequencies: np.ndarray
    samples: int
    burgers_m: np.ndarray
    interpolation_error_J_m: float
    impedance_hermitian_residual: float

    def evaluate(self, theta, derivative=0):
        if derivative not in (0,1,2):
            raise ValueError('only first two angular derivatives supported')
        theta=np.asarray(theta,float)
        if np.any(~np.isfinite(theta)):
            raise ValueError('finite angle required')
        value=np.exp(1j*theta[...,None]*self.frequencies)@(
            self.coefficients*(1j*self.frequencies)**derivative)
        return value.real

    def stiffness(self, theta):
        return self.evaluate(theta)+self.evaluate(theta,2)


def line_energy_coefficient(tensor_pa, burgers_m, *, samples=64):
    """k(theta)=b^T Re[K0(n_perp)] b/(2 pi), from step disregistry.

    K(q)=|q| K0. Both +q and -q enter the real elastic energy. Isotropic
    limits: screw=mu b^2/(4pi), edge=mu b^2/[4pi(1-nu)].
    """
    b=np.array(burgers_m,dtype=float,copy=True)
    if b.shape!=(3,) or np.any(~np.isfinite(b)) or np.linalg.norm(b)==0 or b[2]!=0:
        raise ValueError('explicit nonzero in-plane Burgers vector [m] required')
    if samples<16 or samples%2:
        raise ValueError('even angular sample count >=16 required')
    angles=np.arange(samples)*math.pi/samples
    errors=[]
    def at(t):
        impedance=halfspace_impedance(tensor_pa,[-math.sin(t),math.cos(t)])
        errors.append(impedance.hermitian_residual)
        return float(b@impedance.jump_per_wave_number.real@b/(2*math.pi))
    values=np.array([at(t) for t in angles])
    frequencies=2*np.fft.fftfreq(samples,d=1/samples)
    coeff=np.fft.fft(values)/samples
    test=angles+math.pi/(2*samples)
    approx=(np.exp(1j*test[:,None]*frequencies)@coeff).real
    interpolation=float(np.max(abs(approx-np.array([at(t) for t in test]))))
    for a in (coeff,frequencies,b): a.setflags(write=False)
    return ElasticLineCoefficient(coeff,frequencies,samples,b,interpolation,max(errors))


def pinned_source_branch(line, endpoint_angle, *, span_m, outer_log_ratio, points=129):
    r"""Symmetric double-ended source, LEADING-LOG OUTER elasticity only.

    E[y]=int gamma(theta) dl - tau |b| int y dx,
    gamma(theta)=k(theta)*outer_log_ratio; core and nonlocal finite parts absent.
    From first variation: p=tau|b|, (gamma+gamma'') curvature=p.
    For positive line stiffness and reflection symmetry,
      p L/2=gamma(theta_e) sin(theta_e)+gamma'(theta_e) cos(theta_e).
    Endpoint reaches pi/2 at the local-line critical branch. Outputs are an
    athermal bow-out diagnostic, never a calibrated 0.2% proof stress.
    R/r_core is fixed DURING the variation; R=L is an externally declared
    asymptotic convention, not an independent fitted material parameter.
    """
    vals=np.array([endpoint_angle,span_m,outer_log_ratio],float)
    if np.any(~np.isfinite(vals)) or not(0<endpoint_angle<=math.pi/2) or min(span_m,outer_log_ratio)<=0:
        raise ValueError('positive span/log ratio and angle in (0,pi/2] required')
    if points<9:
        raise ValueError('at least nine curve samples required')
    theta=np.linspace(0,endpoint_angle,points)
    check=np.linspace(-math.pi/2,math.pi/2,4*line.samples+1)
    shape_scale=float(np.max(abs(line.evaluate(check))))
    symmetry=float(np.max(abs(line.evaluate(check)-line.evaluate(-check))))
    arithmetic_floor=64*np.finfo(float).eps*line.samples**2*shape_scale
    if symmetry>max(arithmetic_floor,4*line.interpolation_error_J_m):
        raise ValueError('symmetric pinned branch requires a mirror-symmetric line energy')
    stiffness=line.stiffness(check)
    if np.min(stiffness)<=max(arithmetic_floor,4*line.samples**2*line.interpolation_error_J_m):
        raise ValueError('line stiffness not positively resolved; faceted branch needs separate theory')
    energy=line.evaluate(theta)*outer_log_ratio
    derivative=line.evaluate(theta,1)*outer_log_ratio
    conjugate=energy*np.sin(theta)+derivative*np.cos(theta)
    # Exact angular primitives: q'=T cos(theta), r'=T sin(theta).
    vertical_primitive=-energy*np.cos(theta)+derivative*np.sin(theta)
    pressure=2*conjugate[-1]/span_m
    x=conjugate/pressure
    y=(vertical_primitive[-1]-vertical_primitive)/pressure
    b=float(np.linalg.norm(line.burgers_m))
    return dict(theta=theta,x_m=x,y_m=y,span_m=span_m,outer_log_ratio=outer_log_ratio,
        applied_shear_Pa=float(pressure/b),critical_shear_outer_only_Pa=float(
            2*line.evaluate(math.pi/2)*outer_log_ratio/(b*span_m)),
        line_stiffness_min_J_m=float(np.min(stiffness)*outer_log_ratio),
        symmetry_error_J_m=symmetry,angular_interpolation_error_J_m=line.interpolation_error_J_m,
        core_energy_included=False,finite_part_included=False,finite_source_atomistically_validated=False,
        observable='double_ended_pinned_source_outer_elastic_bowout_NOT_yield',
        experimental_yield_prediction=None,physical_time_available=False)


def solve_pinned_graph(line, *, span_m, outer_log_ratio, shear_Pa, segments=128,
                       force_tolerance=2e-10):
    """Independent conservative piecewise-linear variation of the same energy.

    Unknown y/L, dimensionless energy E/(gamma_ref L); endpoints are fixed.
    Exact segment derivatives: E'=q; E''=(gamma+gamma'')cos(theta)^3/dx.
    Newton iterations are static solves, NOT physical-time dynamics. Returns
    actual convergence, not a stress-based declaration of source activation.
    """
    if segments<8 or segments%2 or not np.isfinite(shear_Pa) or shear_Pa<0:
        raise ValueError('nonnegative shear and even segment count >=8 required')
    probe=pinned_source_branch(line,math.pi/2,span_m=span_m,
        outer_log_ratio=outer_log_ratio,points=17)
    if shear_Pa>=probe['critical_shear_outer_only_Pa']:
        raise ValueError('no subcritical smooth graph branch certified at/above fold')
    if not np.isfinite(force_tolerance) or force_tolerance<=0:
        raise ValueError('positive dimensionless force tolerance required')
    reference=float(line.evaluate(0)*outer_log_ratio)
    dx=1/segments
    pressure=shear_Pa*np.linalg.norm(line.burgers_m)*span_m/reference
    y=np.zeros(segments+1)
    def evaluate(y):
        slopes=np.diff(y)/dx; theta=np.arctan(slopes)
        g=line.evaluate(theta)*outer_log_ratio/reference
        dg=line.evaluate(theta,1)*outer_log_ratio/reference
        stiff=line.stiffness(theta)*outer_log_ratio/reference
        conjugate=g*np.sin(theta)+dg*np.cos(theta)
        gradient=conjugate[:-1]-conjugate[1:]-pressure*dx
        edge=stiff*np.cos(theta)**3/dx
        bands=np.zeros((3,segments-1)); bands[1]=edge[:-1]+edge[1:]
        bands[0,1:]=-edge[1:-1]; bands[2,:-1]=-edge[1:-1]
        energy=float(np.sum(g*np.hypot(dx,np.diff(y)))-pressure*dx*np.sum(y[1:-1]))
        return energy,gradient,bands
    converged=False
    for iteration in range(100):
        energy,gradient,bands=evaluate(y)
        if np.max(abs(gradient))<=force_tolerance:
            converged=True; break
        step=solve_banded((1,1),bands,-gradient)
        accepted=False
        for damping in 2.**(-np.arange(24)):
            trial=y.copy(); trial[1:-1]+=damping*step
            trial_energy=evaluate(trial)[0]
            allowance=16*np.finfo(float).eps*max(1,abs(energy))
            if trial_energy<=energy+1e-4*damping*(gradient@step)+allowance:
                y=trial; accepted=True; break
        if not accepted: break
    energy,gradient,bands=evaluate(y)
    return dict(converged=converged,iterations=iteration,segments=segments,
        x_m=np.linspace(-span_m/2,span_m/2,segments+1),y_m=y*span_m,
        swept_area_m2=float(dx*np.sum(y[1:-1])*span_m**2),
        dimensionless_force_residual=float(np.max(abs(gradient))),
        energy_J=energy*reference*span_m,shear_Pa=shear_Pa,
        experimental_yield_prediction=None,physical_time_available=False)


def swept_area_to_shear_strain(area_m2, burgers_m, specimen_volume_m3):
    """Kinematic volume average b*slipped_area/V for a specified specimen.

    Actual specimen volume, not a fitted activation volume; no probability
    aggregation. One bow-out/first slip is not the 0.2% strain criterion.
    """
    a,b,v=np.asarray([area_m2,burgers_m,specimen_volume_m3],float)
    if np.any(~np.isfinite([a,b,v])) or a<0 or min(b,v)<=0:
        raise ValueError('finite nonnegative area, positive b and specimen volume required')
    return float(a*b/v)
