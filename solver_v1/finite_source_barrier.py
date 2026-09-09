"""Finite pinned-line OUTER elastic barrier, not an atomistic Al activation law.

The same infinite-LJ/Bessel-derived anisotropic line coefficient supplies
gamma(theta)=k(theta)*log(R/r_core). The span and outer/core ratio are explicit
geometry hypotheses, never invented defaults or fitted yield parameters.
Core energy, nonlocal finite parts, obstacles, entropy, source populations
and kinetics remain missing. No Arrhenius rate or PDE clock is returned.
"""
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

from .finite_source_reference import pinned_source_branch


def pinned_line_barrier(line,*,span_m,outer_log_ratio,shear_Pa,
                        quadrature_points=64,relative_tolerance=1e-9):
    r"""Compare metastable and overhanging saddle branches of the SAME line.

    Put p=tau b, T(theta)=gamma+gamma'', Q=gamma sin(theta)+gamma' cos(theta).
    On the symmetric convex branch, Q(theta_0)=pL/2. Its saddle has endpoint
    angle pi-theta_0. The finite positive barrier/area differences are

      DeltaG=2/p int_(theta0)^(pi-theta0) [Q(theta)-Q(theta0)] T sin(theta) dtheta,
      DeltaArea=2/p^2 int_(theta0)^(pi-theta0) Q(theta) T sin(theta) dtheta.

    This avoids subtracting large almost equal total energies near the fold.
    The envelope theorem gives -dDeltaG/dtau=b DeltaArea, with units m^3.
    The quantity is a DERIVED stress derivative, not an inserted event volume.
    """
    if (not np.isfinite(shear_Pa) or shear_Pa<=0 or quadrature_points<16
            or int(quadrature_points)!=quadrature_points
            or not np.isfinite(relative_tolerance) or relative_tolerance<=0):
        raise ValueError('positive subcritical traction and declared quadrature accuracy required')
    check=pinned_source_branch(line,np.pi/2,span_m=span_m,
                              outer_log_ratio=outer_log_ratio,points=17)
    critical=check['critical_shear_outer_only_Pa']
    if shear_Pa>=critical:
        raise ValueError('no metastable pinned branch below an outer-only saddle at/above fold')
    b=float(np.linalg.norm(line.burgers_m));p=shear_Pa*b
    def gamma(theta,derivative=0):return line.evaluate(theta,derivative)*outer_log_ratio
    def Q(theta):return gamma(theta)*np.sin(theta)+gamma(theta,1)*np.cos(theta)
    requested=p*span_m/2
    theta0=brentq(lambda theta:float(Q(theta)-requested),0.,np.pi/2,xtol=2e-14)
    theta1=np.pi-theta0
    def integral(n):
        nodes,weights=leggauss(n)
        theta=np.pi/2+(theta1-theta0)/2*nodes
        scale=(theta1-theta0)/2
        tension=gamma(theta)+gamma(theta,2)
        conjugate=Q(theta)
        barrier=2/p*scale*np.dot(weights,(conjugate-requested)*tension*np.sin(theta))
        area=2/p**2*scale*np.dot(weights,conjugate*tension*np.sin(theta))
        # Independent total-energy difference; less stable near the fold.
        line_energy=2/p*scale*np.dot(weights,gamma(theta)*tension)
        return np.array([barrier,area,line_energy])
    coarse=integral(quadrature_points);fine=integral(2*quadrature_points)
    error=abs(fine-coarse)
    if np.any(fine<=0) or np.any(error>relative_tolerance*abs(fine)):
        raise ArithmeticError('outer barrier/area is not positively quadrature-resolved')
    # A positive numerical integral does not certify neglected core or
    # finite-part energy. Keep that physical limitation in every result.
    return dict(metastable_endpoint_angle=theta0,saddle_endpoint_angle=theta1,
        shear_Pa=float(shear_Pa),outer_only_critical_shear_Pa=float(critical),
        load_fraction=float(shear_Pa/critical),span_m=float(span_m),
        outer_log_ratio=float(outer_log_ratio),
        planar_line_minimum_morse_index=0,planar_line_saddle_morse_index=1,
        stability_scope='Dirichlet normal variations of the convex local anisotropic line, not atomistic stability',
        outer_only_barrier_J=float(fine[0]),additional_swept_area_m2=float(fine[1]),
        outer_only_stress_derivative_m3=float(b*fine[1]),
        line_energy_difference_J=float(fine[2]),
        total_energy_identity_residual_J=float(fine[2]-p*fine[1]-fine[0]),
        quadrature_points=2*quadrature_points,quadrature_change=error,
        line_interpolation_error_J_m=line.interpolation_error_J_m,
        core_energy_included=False,nonlocal_finite_part_included=False,
        finite_source_atomistically_validated=False,experimental_activation_fitted=False,
        kinetic_rate=None,physical_a_s_mobility=None,production_probability=None)


def normal_mode_second_variation(*,endpoint_angle,pressure_J_m2,mode):
    r"""Exact unit-amplitude Dirichlet mode curvature, within this line model.

    At equilibrium kappa=p/T and d ell=T/p d theta. Therefore
      delta²G=p int[-theta_e,theta_e] [(eta_theta)^2-eta^2] dtheta.
    For eta=amplitude*sin[j*pi*(theta+theta_e)/(2theta_e)], curvature is
      p theta_e [(j*pi/(2theta_e))²-1] [J/m²].
    Minor arc theta_e<pi/2: all positive. Major arc pi/2<theta_e<pi:
    exactly j=1 is negative. The fold has one zero mode. This does not prove
    atomistic/core/nonlocal stability beyond the specified local-line energy.
    """
    if (not np.isfinite(endpoint_angle+pressure_J_m2) or not 0<endpoint_angle<np.pi
            or pressure_J_m2<=0 or int(mode)!=mode or mode<1):
        raise ValueError('positive pressure, angle in (0,pi), positive integer mode required')
    return float(pressure_J_m2*endpoint_angle*((mode*np.pi/(2*endpoint_angle))**2-1))
