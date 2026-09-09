"""Stress-controlled STATIC research probes with explicit conjugate units.

No mobility, model time, probability evolution, specimen correlation area, or
empirical slip threshold enters these tests. Load continuation distinguishes
a resolved stable root from failed continuation; failure is not a spinodal.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.optimize import least_squares

from .aluminum_calibration import EV_J


@dataclass(frozen=True)
class InterfaceUnits:
    length_scale_m: float
    atomic_cell_area_m2: float

    def __post_init__(self):
        values=np.array([self.length_scale_m,self.atomic_cell_area_m2])
        if np.any(~np.isfinite(values)) or np.any(values<=0):
            raise ValueError("finite positive microscopic coordinate and cell scales required")

    def traction_to_force(self,traction_gpa):
        """Conjugate eV per reduced coordinate; positive opens/increases slip."""
        return np.asarray(traction_gpa)*1e9*self.atomic_cell_area_m2*self.length_scale_m/EV_J

    def force_to_traction(self,force):
        return np.asarray(force)*EV_J/(1e9*self.atomic_cell_area_m2*self.length_scale_m)

    def traction_mpa_to_force(self, traction_mpa):
        """Explicit MPa boundary: avoid unit-ambiguous research load arrays."""
        return self.traction_to_force(np.asarray(traction_mpa)/1000.)

    def force_to_traction_mpa(self, force):
        return self.force_to_traction(force)*1000.

    def energy_to_surface(self,energy_ev):
        return np.asarray(energy_ev)*EV_J/self.atomic_cell_area_m2


def resolved_uniaxial_tractions(stress_gpa,loading_direction,plane_normal,slip_direction):
    """Explicit crystallographic projection, not a fitted phenomenological chi."""
    vectors=np.asarray([loading_direction,plane_normal,slip_direction],dtype=float)
    if vectors.shape!=(3,3) or not np.all(np.isfinite(vectors)):
        raise ValueError("three finite three-dimensional vectors required")
    if not np.allclose(np.linalg.norm(vectors,axis=1),1.,atol=1e-12,rtol=0):
        raise ValueError("orientation vectors must be normalized")
    e,n,m=vectors
    if abs(n@m)>1e-12:
        raise ValueError("slip direction must lie in interface plane")
    return np.array([stress_gpa*(e@n)**2,stress_gpa*(e@n)*(e@m)])


def resolved_tensor_tractions(stress_tensor, plane_normal, slip_direction):
    """Project a symmetric 3D Cauchy stress into (normal, slip1, slip2).

    Units are preserved (e.g. input MPa -> output MPa). Tangent2=n cross m.
    ALL vectors/tensor must be in the same explicitly supplied physical frame.
    This is a geometry/work helper, NOT a three-state or spatial 3D PDE.
    A scalar-s model cannot silently discard the returned second shear load.
    """
    stress = np.asarray(stress_tensor, dtype=float)
    n, m = np.asarray(plane_normal, dtype=float), np.asarray(slip_direction, dtype=float)
    if (stress.shape != (3,3) or not np.all(np.isfinite(stress))
            or not np.allclose(stress, stress.T, rtol=1e-12, atol=1e-12)):
        raise ValueError("finite symmetric 3x3 stress tensor required")
    if (n.shape != (3,) or m.shape != (3,) or not np.all(np.isfinite([n,m]))
            or not np.allclose([np.linalg.norm(n),np.linalg.norm(m)],1.,atol=1e-12,rtol=0)
            or abs(n@m)>1e-12):
        raise ValueError("orthogonal unit plane normal and slip direction required")
    basis = np.array([n, m, np.cross(n,m)])
    return basis @ (stress @ n)


def packed_evaluation(interface,a,s):
    value=interface.evaluate(float(a),float(s))
    return np.array([value.energy,value.d_da,value.d_ds,value.d2_daa,value.d2_das,value.d2_dss])


def hessian_from_packed(value):
    return np.array([[value[3],value[4]],[value[4],value[5]]])


def continue_static_branch(evaluate,*,h,period,units,loads_gpa,max_traction_step=.1):
    """Follow the intact local minimum for an ordered sequence of (Tn,tau).

    Each prescribed step is subdivided in PHYSICAL traction, independent of
    the state solver. A failed or non-positive-Hessian step stops the branch:
    jumping to another well would fabricate branch continuation. Restricting
    root searches to a local numerical window does not alter the energy.
    """
    if h<=0 or period<=0 or max_traction_step<=0:
        raise ValueError("positive geometry and continuation step required")
    loads=np.asarray(loads_gpa,dtype=float)
    if loads.ndim!=2 or loads.shape[1]!=2 or np.any(~np.isfinite(loads)):
        raise ValueError("finite normal/shear load pairs required")
    q=np.array([h,0.]); previous=np.zeros(2); alive=True; rows=[]
    for index,load in enumerate(loads):
        status="stable_local_equilibrium"; residual=float("nan"); substeps=0
        if alive:
            count=max(1,int(np.ceil(np.max(np.abs(load-previous))/max_traction_step)))
            for fraction in np.linspace(1/count,1,count):
                subload=previous+fraction*(load-previous)
                conjugate=units.traction_to_force(subload)
                @lru_cache(maxsize=8)
                def value_at(a,s):
                    return np.asarray(evaluate(a,s),dtype=float)
                def fun(x):
                    return value_at(*x)[1:3]-conjugate
                def jac(x):
                    return hessian_from_packed(value_at(*x))
                solved=least_squares(fun,q,jac=jac,
                    bounds=([.55*h,q[1]-.25*period],[4*h,q[1]+.25*period]),
                    x_scale=[h,period],xtol=2e-12,ftol=2e-12,gtol=2e-12,max_nfev=80)
                val=value_at(*solved.x); eig=np.linalg.eigvalsh(hessian_from_packed(val))
                residual=float(np.max(np.abs(fun(solved.x)))); substeps+=1
                if not solved.success or residual>2e-8:
                    status="continuation_unresolved_or_branch_lost"
                    alive=False; break
                if eig[0]<=0:
                    status="stationary_root_not_stable"
                    alive=False; break
                q=solved.x
        else:
            status="not_evaluated_after_branch_loss"
        if alive:
            value=np.asarray(evaluate(*q))
            eig=np.linalg.eigvalsh(hessian_from_packed(value))
            row=dict(a=float(q[0]),s=float(q[1]),energy_ev_cell=float(value[0]),
                     loaded_energy_ev_cell=float(value[0]-units.traction_to_force(load)@np.array([q[0]-h,q[1]])),
                     min_hessian_eigenvalue=float(eig[0]),
                     normal_strain=float((q[0]-h)/h),slip_over_period=float(q[1]/period))
        else:
            row=dict(a=None,s=None,energy_ev_cell=None,loaded_energy_ev_cell=None,
                     min_hessian_eigenvalue=None,normal_strain=None,slip_over_period=None)
        rows.append(dict(step=index,normal_traction_GPa=float(load[0]),
                         shear_traction_GPa=float(load[1]),status=status,
                         force_residual_ev_coordinate=residual,substeps=substeps,**row))
        previous=load.copy()
    return rows
