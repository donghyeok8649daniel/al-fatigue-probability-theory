"""Bounded deterministic existing-shape search with finite-q AND interface jets."""
import argparse
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from .run_finite_q_compatibility_v31 import points
from .core_interface_compatibility import ExistingJetSubset,minimax_compatibility
from .coordination_screening import CoordinationScreenedBulk
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents,NONNEGATIVE
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_current_material_core import ROOT,load_current_material
from .run_source_core_reference import load_source_material
from .vector_material_calibration import LENGTH_M,MaterialObservation
from .run_vector_registry_audit import save_json
from .run_low_frequency_forcing_v29 import sha
import json


def shape_coordinates(shape, vary_screening=False):
    """Five log-positive radial variables; existing power z is signed."""
    shape=np.asarray(shape,float)
    if shape.shape!=(6,) or np.any(~np.isfinite(shape)) or np.any(shape[:5]<=0) or not -1<=shape[5]<=1:
        raise ValueError('positive existing radial shape and power exponent in [-1,1] required')
    return np.r_[np.log(shape[:5]),shape[5]] if vary_screening else np.log(shape[:5])


def decode_shape(x, fixed_screening, vary_screening=False):
    x=np.asarray(x,float)
    if x.shape!=((6,) if vary_screening else (5,)) or np.any(~np.isfinite(x)):
        raise ValueError('matching finite optimization coordinates required')
    shape=np.r_[np.exp(x[:5]),x[5] if vary_screening else fixed_screening]
    shape_coordinates(shape,vary_screening)
    if 2*shape[4]+min(shape[5],0)*shape[0]<=0:
        raise ValueError('isolated-neighbor environmental energy must decay')
    return shape


def run(out,maxfev,start_from=None,vary_screening=False):
    out=Path(out)
    if out.exists():raise FileExistsError('fresh radial search required')
    source,_,binding=load_source_material();operator=SourceEAMBlochHessian(source.source)
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    old=impose_tangents(problem.observations,problem.source_tangents)
    selected=[i for i,o in enumerate(old) if o.role=='exact' or o.units=='eV/L0^2']
    subset=ExistingJetSubset(problem,selected);obs=[old[i] for i in subset.rows]
    fitpoints=[q for role,q in points() if role=='fit'];components=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    for k,q in enumerate(fitpoints):
        H=operator.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
        for i,j in components:
            obs.append(MaterialObservation(f'v31_q{k}_H{i}{j}',float(H[i,j]*(1 if i==j else np.sqrt(2))),
                float(.05*np.linalg.norm(H)),'eV/L0^2','fit'))
    rows=[i for i,o in enumerate(obs) if o.role!='exact']
    model,_,mb=load_current_material(ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json')
    start=model.screened_shape.copy();lower=np.log([.7,2.,2.,1e-5,2.]);upper=np.log([8.,12.,12.,1e6,12.])
    start_hash=None
    if start_from is not None:
        prior=json.loads(Path(start_from).read_bytes())
        if not prior['completed'] or not prior['best']['positive_LJ']:raise ValueError('completed positive-LJ starting study required')
        proposed=np.asarray(prior['best']['shape'])
        if proposed.shape!=start.shape or (not vary_screening and proposed[-1]!=start[-1]):raise ValueError('same existing family required')
        start=proposed;start_hash=sha(start_from)
    if vary_screening:
        lower=np.r_[lower,-1.];upper=np.r_[upper,1.]
    out.mkdir(parents=True);began=time.perf_counter();history=[]
    save_json(out/'definition.json',dict(source=binding,start_binding=mb,start_shape=start,
        log_bounds=[lower,upper],maxfev=maxfev,optimizer='deterministic bounded Powell',followup_start_sha256=start_hash,
        objective='joint old115 inspected curvatures plus predeclared finite-q fit matrices',
        vary_existing_power_screening=vary_screening,
        shape_coordinate_convention='first five logarithmic, optional sixth linear exponent in inherited [-1,1]',
        withheld_points=[q for role,q in points() if role=='excluded'],
        no_spectral_admissibility_certificate=True,new_energy_terms=False,actual_yield_in_loss=False))
    replay=np.max(abs(subset.matrix(tuple(start))-problem.matrix(tuple(start))[subset.rows]))
    if replay>2e-10:raise ArithmeticError('subset/full replay failure')
    def fun(x):
        tick=time.perf_counter();shape=None
        try:
            shape=decode_shape(x,start[5],vary_screening)
            base=subset.matrix(tuple(shape));bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
            extra=[]
            for q in fitpoints:
                columns,_=bulk.evaluate(q)
                extra.extend(columns[:,i,j]*(1 if i==j else np.sqrt(2)) for i,j in components)
            fit=minimax_compatibility(np.vstack([base,extra]),obs,rows,nonnegative=NONNEGATIVE)
            if not fit['completed']:raise ArithmeticError('LP not certified')
            value=fit['minimax_normalized_error']
            entry=dict(shape=shape,value=value,coefficients=fit['coefficients'],
                positive_LJ=fit['strictly_positive_LJ'],duality_gap=fit['duality_gap'],completed=True)
        except (ValueError,ArithmeticError) as error:
            value=np.inf;entry=dict(shape=shape,completed=False,error=str(error))
        history.append(dict(**entry,seconds=time.perf_counter()-tick))
        save_json(out/'checkpoint.json',dict(completed=False,history=history))
        print('joint shape',len(history),'eta',value,flush=True)
        return value
    result=minimize(fun,shape_coordinates(start,vary_screening),method='Powell',bounds=list(zip(lower,upper)),
        options=dict(maxfev=maxfev,xtol=2e-4,ftol=2e-5))
    eligible=[r for r in history if r['completed'] and r['positive_LJ']]
    save_json(out/'summary.json',dict(completed=True,optimizer_success=bool(result.success),
        message=str(result.message),history=history,best=min(eligible,key=lambda x:x['value']) if eligible else None,
        subset_replay_error=float(replay),elapsed_seconds=time.perf_counter()-began,
        family_impossibility_proved=False,material_accepted=False,actual_yield_validated=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--maxfev',type=int,default=40)
    p.add_argument('--start-from',type=Path)
    p.add_argument('--vary-screening',action='store_true')
    a=p.parse_args();run(a.out,a.maxfev,a.start_from,a.vary_screening)
