"""Research joint bulk/interface calibration with matched rigid-EAM targets.

LJ stays analytic. A source EAM file supplies targets only. No tabulated source
function enters this model, its derivatives, or the production probability PDE.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import minimize_scalar, nnls, root_scalar

from .aluminum_full_fcc_calibration import FullFCCParameters, build_full_fcc_calibration_model
from .fcc111_active_interface import FCC111ActiveInterface
from .fcc111_full_energy import FullFCC111StackEnergy
from .fcc111_geometry import DIRECT_110, SHOCKLEY_112
from .full_fcc_calibration_audit import IDEAL_H, NAMES, bulk_basis, independent_bulk_targets
from .model import ModelParams
from .reference_eam_targets import MishinRigidFCCReference


INTERFACE_STATES=(
    ("intrinsic_fault_rigid",1.,1/math.sqrt(3),SHOCKLEY_112),
    ("shockley_half_partial_rigid",1.,.5/math.sqrt(3),SHOCKLEY_112),
    ("direct_half_registry_rigid",1.,.5,DIRECT_110),
    ("opening_20percent_rigid",1.2,0.,DIRECT_110),
    ("work_separation_rigid",math.inf,0.,DIRECT_110),
)
JOINT_NAMES=(*NAMES,*(row[0] for row in INTERFACE_STATES))


@dataclass(frozen=True)
class MinimalConvexEmbedding:
    """F(x)=-A sqrt(x)+B(x-1)+C(x-1)^2, x=rho/rho_ref.

    C is a one-term *research* extension; do not activate it until the baseline
    joint fit is measured. F(0)=-B+C must be included in atomization energy.
    The old embedding classes/default models remain unchanged.
    """
    sqrt_ev: float
    linear_ev: float
    quadratic_ev: float
    rho_ref: float

    def __post_init__(self):
        values=np.array([self.sqrt_ev,self.linear_ev,self.quadratic_ev,self.rho_ref])
        if np.any(~np.isfinite(values)) or self.sqrt_ev<0 or self.quadratic_ev<0 or self.rho_ref<=0:
            raise ValueError("finite coefficients, A,C>=0 and rho_ref>0 required")

    def value(self,density):
        x=np.asarray(density,dtype=float)/self.rho_ref
        if np.any(x<0):
            raise ValueError("negative environment density")
        return -self.sqrt_ev*np.sqrt(x)+self.linear_ev*(x-1)+self.quadratic_ev*(x-1)**2

    def first_derivative(self,density):
        x=np.asarray(density,dtype=float)/self.rho_ref
        if np.any(x<=0):
            raise ValueError("positive density required for derivative")
        return (-self.sqrt_ev/(2*np.sqrt(x))+self.linear_ev+2*self.quadratic_ev*(x-1))/self.rho_ref

    def second_derivative(self,density):
        x=np.asarray(density,dtype=float)/self.rho_ref
        if np.any(x<=0):
            raise ValueError("positive density required for derivative")
        return (self.sqrt_ev/(4*x**1.5)+2*self.quadratic_ev)/self.rho_ref**2


def matched_targets(reference=None):
    source=reference or MishinRigidFCCReference()
    target,scale=independent_bulk_targets()
    # Original rounded Mishin bulk constants remain fixed; matched reference
    # source calculations constrain the interface under identical rigid halves.
    values=[]
    for name,opening,slip,path in INTERFACE_STATES:
        value=(source.work_separation() if math.isinf(opening) else
               source.interface_energy(opening*source.h,slip*source.geometry.b,path_id=path))
        values.append(value)
    values=np.asarray(values)
    # Explicit 10% model-discrepancy allowance, not statistical measurement errors.
    return np.r_[target,values],np.r_[scale,.10*values]


@lru_cache(maxsize=1)
def interface_pair_basis():
    models=[]
    for sigma in (1.,2**(1/6)):
        model=FullFCC111StackEnergy(ModelParams(n_cells=1,b=1.,epsilon=.25,
            sigma_lj=sigma),find_equilibrium=False)
        model.a0=IDEAL_H
        models.append(FCC111ActiveInterface(model,tolerance=2e-11))
    rows=[]
    for name,opening,slip,path in INTERFACE_STATES:
        if math.isinf(opening):
            first,second=[model.separated_limit() for model in models]
            s6=(second-2*first)/2
            minus_s3=first-s6
        else:
            model=FCC111ActiveInterface(models[0].bulk,path_id=path,tolerance=2e-11)
            s6=model._power_cross_difference(opening*IDEAL_H,slip,6.)[0][0]
            minus_s3=-model._power_cross_difference(opening*IDEAL_H,slip,3.)[0][0]
        rows.append([s6,minus_s3])
    return np.asarray(rows)


@lru_cache(maxsize=256)
def joint_basis(decay):
    """Exact linear coefficient dependence at fixed decay and reference geometry."""
    bulk=build_full_fcc_calibration_model(FullFCCParameters(.1,1.,float(decay),1.))
    bulk.a0=IDEAL_H
    rho=bulk.full_environment_density(IDEAL_H,0)
    x_alpha=IDEAL_H*rho.d_da/rho.value
    x_gamma=IDEAL_H*rho.d_ds/rho.value
    c_column=np.array([0.,1.,2*(3*x_alpha)**2,2*x_alpha**2,2*x_gamma**2])
    bulk_columns=np.column_stack((bulk_basis(float(decay)),c_column))
    rows=[]
    for name,opening,slip,path in INTERFACE_STATES:
        interface=FCC111ActiveInterface(bulk,path_id=path,tolerance=2e-11)
        # Density cross contribution is exponentially absent here; LJ separation
        # uses its exact infinite-tail limit separately above.
        a=400*IDEAL_H if math.isinf(opening) else opening*IDEAL_H
        changes,_,_,_=interface._density_changes(a,slip)
        dx=np.array([row[0]/rho.value for row in changes])
        if np.any(1+dx<=0):
            raise ValueError("nonpositive interface density in coefficient basis")
        rows.append([-2*np.sum(np.sqrt(1+dx)-1),2*np.sum(dx),2*np.sum(dx**2)])
    return np.vstack((bulk_columns,np.column_stack((interface_pair_basis(),rows))))


@dataclass(frozen=True)
class JointFit:
    family: str
    decay: float
    coefficients: np.ndarray
    squared_loss: float
    residuals: np.ndarray
    predictions: np.ndarray


def fixed_decay_fit(decay,target,scales,*,family):
    count={"sqrt":3,"linear":4,"quadratic":5}[family]
    basis=joint_basis(float(decay))[:,:count]
    coefficients,_=nnls(basis/scales[:,None],target/scales,maxiter=500)
    predictions=basis@coefficients
    residuals=(predictions-target)/scales
    coeff=np.r_[coefficients,np.zeros(5-count)]
    return JointFit(family,float(decay),coeff,float(residuals@residuals),residuals,predictions)


def joint_profile(target,scales,*,family,grid=None):
    grid=np.geomspace(.35,8.,33) if grid is None else np.asarray(grid)
    rows=[fixed_decay_fit(d,target,scales,family=family) for d in grid]
    candidates=rows.copy()
    for i in range(1,len(rows)-1):
        if rows[i].squared_loss<=min(rows[i-1].squared_loss,rows[i+1].squared_loss):
            solved=minimize_scalar(lambda d:fixed_decay_fit(d,target,scales,family=family).squared_loss,
                bounds=(grid[i-1],grid[i+1]),method="bounded",options={"xatol":2e-6})
            candidates.append(fixed_decay_fit(solved.x,target,scales,family=family))
    eligible=[row for row in candidates if np.all(row.coefficients[:2]>0)]
    if not eligible:
        raise ValueError("joint profile found no positive-LJ candidate")
    return min(eligible,key=lambda row:row.squared_loss),candidates


def build_joint_model(fit,*,b=1.,rho_ref=None):
    u,v,A,B,C=fit.coefficients
    if u<=0 or v<=0:
        raise ValueError("LJ repulsive and attractive coefficients must be positive")
    sigma=(u/v)**(1/6); epsilon=v*v/(4*u)
    probe=build_full_fcc_calibration_model(FullFCCParameters(epsilon,sigma,fit.decay,1.),
                                         b=b,rho_ref=rho_ref)
    probe.embedding=MinimalConvexEmbedding(float(A),float(B),float(C),probe.embedding.rho_ref)
    return probe


def joint_equilibrium(fit):
    reference=build_joint_model(fit)
    rho_ref=reference.embedding.rho_ref
    def derivative(stretch):
        return 3*IDEAL_H*build_joint_model(fit,b=stretch,rho_ref=rho_ref).grad(IDEAL_H*stretch,0)[0]
    grid=np.linspace(.85,1.15,31)
    previous=derivative(grid[0]); roots=[]
    for left,right in zip(grid[:-1],grid[1:]):
        current=derivative(right)
        if previous<=0<=current:
            root=root_scalar(derivative,bracket=(left,right),xtol=1e-11).root
            roots.append(root)
        previous=current
    if not roots:
        raise ValueError("no stable FCC equilibrium near target")
    stretch=min(roots,key=lambda x:abs(x-1))
    model=build_joint_model(fit,b=stretch,rho_ref=rho_ref)
    model.a0=IDEAL_H*stretch; model.equilibrium_error=None
    return model,stretch
