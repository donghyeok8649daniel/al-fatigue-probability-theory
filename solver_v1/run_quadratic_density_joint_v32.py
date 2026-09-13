"""Deterministic material-target search of a three-parameter radial shape.

Auxiliary source density is ONLY the declared starting shape. It is not in
this objective. Angular shapes are held, amplitudes are conditionally fitted.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from .run_quadratic_density_material_v32 import interface_matrix
from .polynomial_exponential_density import QuadraticEnvelopeDensity
from .quadratic_density_interface import QuadraticScalarEnvironment, bulk_scalar_columns
from .coordination_screening import CoordinationScreenedBulk
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .core_interface_compatibility import minimax_compatibility
from .run_current_material_core import ROOT
from .run_finite_q_compatibility_v31 import points
from .run_source_core_reference import load_source_material
from .periodic_plane_covariance import SourceEAMBlochHessian
from .vector_material_calibration import LENGTH_M, MaterialObservation
from .run_vector_registry_audit import save_json
from .run_low_frequency_forcing_v29 import sha


def encode(kernel):
    return np.array([np.log(kernel.k),kernel.center,np.log(kernel.width)])


def decode(x):
    x=np.asarray(x,float)
    if x.shape!=(3,) or not np.all(np.isfinite(x)):
        raise ValueError('three finite radial shape coordinates required')
    return QuadraticEnvelopeDensity(1.,np.exp(x[0]),x[1],np.exp(x[2]))


def joint_decode(x,fixed_shape,all_ranges=False):
    x=np.asarray(x,float)
    if x.shape!=((8,) if all_ranges else (3,)):
        raise ValueError('matching three/eight optimization coordinates required')
    kernel=decode(x[:3]);shape=np.asarray(fixed_shape,float).copy()
    if all_ranges:
        shape[1:5]=np.exp(x[3:7]);shape[5]=x[7]
    if not np.all(np.isfinite(shape)) or not -1<=shape[5]<=1:
        raise ValueError('finite angular shape and inherited screening bounds required')
    if 2*shape[4]+min(shape[5],0.)*kernel.k<=0:
        raise ValueError('new scalar/screening isolated-neighbor decay condition failed')
    return kernel,shape


def run(parent,out,maxfev,all_ranges=False):
    parent,out=Path(parent),Path(out)
    if out.exists():
        raise FileExistsError('fresh study required')
    if maxfev<1:
        raise ValueError('positive evaluation budget required')
    prior=json.loads((parent/'summary.json').read_bytes())
    definition=json.loads((parent/'definition.json').read_bytes())
    if not prior['completed']:
        raise ValueError('completed fixed-shape material study required')
    shape=definition['parent_angular_shape']
    initial=QuadraticEnvelopeDensity(**prior['normalized_kernel'])
    L=LENGTH_M/1e-10
    lower=np.array([np.log(1.*L),0.,np.log(.05/L)])
    upper=np.array([np.log(6.*L),3./L,np.log(2./L)])
    start=encode(initial)
    if all_ranges:
        start=np.r_[start,np.log(shape[1:5]),shape[5]]
        lower=np.r_[lower,np.log([2.,2.,1e-5,2.]),-1.]
        upper=np.r_[upper,np.log([12.,12.,1e6,12.]),1.]
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    obs=list(impose_tangents(problem.observations,problem.source_tangents))
    original_count=len(obs)
    rows=[i for i,o in enumerate(obs) if o.role!='exact' and o.units=='eV/L0^2']
    source,_,binding=load_source_material();operator=SourceEAMBlochHessian(source.source)
    bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
    components=((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
    operators=[];fitpoints=[q for role,q in points() if role=='fit']
    for k,q in enumerate(fitpoints):
        operators.append(bulk.evaluate(q)[0])
        target=operator.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*L*L
        for i,j in components:
            w=1. if i==j else np.sqrt(2)
            obs.append(MaterialObservation(f'q{k}_H{i}{j}',float(w*target[i,j]),
                float(.05*np.linalg.norm(target)),'eV/L0^2','fit'))
    rows+=list(range(original_count,len(obs)))
    out.mkdir(parents=True);began=time.perf_counter();history=[]
    save_json(out/'definition.json',dict(parent_sha256=sha(parent/'summary.json'),source=binding,
        angular_shape_fixed=None if all_ranges else shape,initial_angular_shape=shape,
        all_existing_angular_ranges_varied=all_ranges,
        initial_kernel=vars(initial),coordinate_bounds=[lower,upper],
        bound_units='log(k*), center/L0, log(width/L0); optional log(k3),log(k2),log(alpha2),log(k1),z',
        physical_bounds='k in [1,6]/Angstrom, center in [0,3] Angstrom, width in [.05,2] Angstrom',
        bound_reason='finite exploratory domain around source auxiliary shape; not material priors',
        maxfev=maxfev,optimizer='bounded deterministic Powell',
        objective='same 115 inspected interface curvatures and four finite-q matrices; seven exact anchors',
        auxiliary_density_in_loss=False,withheld_q=[q for role,q in points() if role=='excluded'],
        fresh_states=definition['fresh_states'],production_changed=False))
    def objective(x):
        tick=time.perf_counter()
        try:
            kernel,current_shape=joint_decode(x,shape,all_ranges)
            environment=QuadraticScalarEnvironment(kernel)
            matrix=interface_matrix(problem,current_shape,environment);extra=[]
            current_bulk=CoordinationScreenedBulk(current_shape,radius=16.,law='power') if all_ranges else bulk
            current_operators=[current_bulk.evaluate(q)[0] for q in fitpoints] if all_ranges else operators
            for q,base in zip(fitpoints,current_operators):
                columns=base.copy()
                columns[2:5],_=bulk_scalar_columns(current_bulk,environment,q)
                extra.extend(columns[:,i,j]*(1 if i==j else np.sqrt(2)) for i,j in components)
            fit=minimax_compatibility(np.vstack([matrix,extra]),obs,rows,nonnegative=NONNEGATIVE)
            if not fit['completed']:
                raise ArithmeticError(f"coefficient LP failed: status={fit['lp_status']}, {fit['message']}")
            value=fit['minimax_normalized_error']
            entry=dict(completed=True,value=value,kernel=vars(environment.kernel),angular_shape=current_shape,
                coefficients=fit['coefficients'],positive_LJ=fit['strictly_positive_LJ'],duality_gap=fit['duality_gap'])
        except (ValueError,ArithmeticError) as error:
            value=np.inf;entry=dict(completed=False,error=str(error),coordinates=x)
        history.append(dict(**entry,seconds=time.perf_counter()-tick))
        save_json(out/'checkpoint.json',dict(completed=False,history=history))
        print('quadratic material',len(history),'eta',value,flush=True)
        return value
    result=minimize(objective,start,method='Powell',bounds=list(zip(lower,upper)),
        options=dict(maxfev=maxfev,xtol=2e-4,ftol=2e-5))
    eligible=[r for r in history if r['completed'] and r['positive_LJ']]
    save_json(out/'summary.json',dict(completed=True,optimizer_success=bool(result.success),
        message=str(result.message),history=history,best=min(eligible,key=lambda r:r['value']) if eligible else None,
        elapsed_seconds=time.perf_counter()-began,material_accepted=False,
        family_impossibility_proved=False,production_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--parent',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--maxfev',type=int,default=120)
    p.add_argument('--all-ranges',action='store_true')
    args=p.parse_args();run(args.parent,args.out,args.maxfev,args.all_ranges)
