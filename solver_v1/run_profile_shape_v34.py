"""Dual-guided calibration of unchanged existing analytic radial shapes.

Audit envelope derivatives against freshly solved coefficient LPs first.
Only calibration derivatives use finite differences; energy jets remain
analytic. All target scales, LJ coefficients and production gates are kept.
"""
import argparse
from dataclasses import asdict
from functools import lru_cache
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from .core_interface_compatibility import ExistingJetSubset,minimax_compatibility
from .profile_shape_sensitivity import envelope_gradient,difference_stencil
from .interface_tangent_calibration import TangentCalibrationProblem,NONNEGATIVE
from .coordination_screening import CoordinationScreenedBulk
from .run_joint_shape_v31 import shape_coordinates,decode_shape
from .run_finite_q_compatibility_v31 import points
from .run_current_material_core import ROOT
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json
from .vector_material_calibration import MaterialObservation


class ShapeProfile:
    def __init__(self,parent,out):
        self.out=Path(out);parent=Path(parent)
        if self.out.exists():raise FileExistsError('fresh profile study required')
        definition=json.loads((parent/'definition.json').read_bytes())
        saved=json.loads((parent/'profiles.json').read_bytes())[0]
        allobs=[MaterialObservation(**o) for o in json.loads((parent/'observations.json').read_bytes())]
        self.start=shape_coordinates(definition['shape'],True)
        self.lower=np.r_[np.log([.7,2.,2.,1e-5,2.]),-1.]
        self.upper=np.r_[np.log([8.,12.,12.,1e6,12.]),1.]
        problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
        count=len(problem.observations)
        selected=[i for i,o in enumerate(allobs[:count]) if o.role=='exact' or o.units=='eV/L0^2']
        self.subset=ExistingJetSubset(problem,selected)
        self.obs=[allobs[i] for i in self.subset.rows]
        self.fitpoints=[q for role,q in points() if role=='fit']
        qindices=[count+6*k+j for k,(role,q) in enumerate(points()) if role=='fit' for j in range(6)]
        self.obs.extend(allobs[i] for i in qindices)
        self.rows=[i for i,o in enumerate(self.obs) if o.role!='exact']
        self.out.mkdir(parents=True);self.history=[];self.gradient_history=[]
        self.began=time.perf_counter()
        save_json(self.out/'definition.json',dict(parent_sha256=sha(parent/'profiles.json'),
            shape=definition['shape'],lower=self.lower,upper=self.upper,
            coordinates='five log-positive existing shapes and existing signed exponent',
            target_scales_changed=False,new_energy_terms=False,
            derivative='LP envelope from equality/upper/lower duals; FD of analytic observable matrix',
            gradient_steps=[2e-4,1e-4],excluded_points_in_loss=False,
            finite_difference_energy_generator=False,production_changed=False))
        save_json(self.out/'observations.json',[asdict(o) for o in self.obs])
        M=self.matrix(tuple(self.start))
        replay=float(np.max(abs(M@saved['coefficients']-np.asarray(saved['predictions'])[self.subset.rows+qindices])))
        if replay>2e-9:raise ArithmeticError('parent coefficient/target/unit replay mismatch')
        self.replay=replay

    @lru_cache(maxsize=64)
    def matrix(self,x):
        shape=decode_shape(x,1.,True)
        M=self.subset.matrix(tuple(shape))
        bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
        components=((0,0),(1,1),(2,2),(0,1),(0,2),(1,2));extra=[]
        for q in self.fitpoints:
            columns,_=bulk.evaluate(q)
            extra.extend(columns[:,i,j]*(1. if i==j else np.sqrt(2.)) for i,j in components)
        return np.vstack([M,extra])

    @lru_cache(maxsize=128)
    def profile(self,x):
        tick=time.perf_counter()
        p=minimax_compatibility(self.matrix(x),self.obs,self.rows,nonnegative=NONNEGATIVE)
        if not p['completed']:raise ArithmeticError('uncertified coefficient profile')
        self.history.append(dict(x=x,shape=decode_shape(x,1.,True),value=p['minimax_normalized_error'],
            positive_LJ=p['strictly_positive_LJ'],coefficients=p['coefficients'],seconds=time.perf_counter()-tick,
            duality_gap=p['duality_gap']))
        save_json(self.out/'checkpoint.json',dict(completed=False,history=self.history,gradient_history=self.gradient_history))
        print('profile',len(self.history),'eta',p['minimax_normalized_error'],flush=True)
        return p

    def gradient(self,x,step=1e-4,reprofile=False):
        x=np.asarray(x,float);p=self.profile(tuple(x));derivatives=[];finite=[]
        for i in range(len(x)):
            stencil=difference_stencil(x,i,self.lower,self.upper,step)
            derivatives.append(sum(w*self.matrix(tuple(v)) for v,w in stencil))
            if reprofile:
                finite.append(sum(w*self.profile(tuple(v))['minimax_normalized_error'] for v,w in stencil))
        g=envelope_gradient(np.asarray(derivatives),self.obs,p)
        report=dict(x=x,step=step,envelope=g,profile_finite_difference=finite if reprofile else None)
        self.gradient_history.append(report)
        save_json(self.out/'gradient_history.json',self.gradient_history)
        return g,report


def run(parent,out,maxiter):
    if maxiter<0:raise ValueError('nonnegative iteration budget required')
    study=ShapeProfile(parent,out);audit=[]
    for step in (2e-4,1e-4):
        g,report=study.gradient(study.start,step,reprofile=True);audit.append(report)
    errors=[float(np.linalg.norm(np.asarray(a['envelope'])-a['profile_finite_difference'],np.inf)) for a in audit]
    change=float(np.linalg.norm(np.asarray(audit[0]['envelope'])-audit[1]['envelope'],np.inf))
    scale=max(1.,float(np.linalg.norm(g,np.inf)))
    verified=bool(max(errors+[change])<=1e-4*scale)
    projected=np.clip(study.start-g,study.lower,study.upper)-study.start
    save_json(Path(out)/'gradient_audit.json',dict(completed=True,audit=audit,absolute_errors=errors,
        refinement_change=change,algorithmic_relative_tolerance=1e-4,derivative_verified=verified,
        projected_descent=projected,parent_replay_error=study.replay,
        nonsmooth_warning='single dual is not a unique derivative across active-set changes'))
    result=None
    if maxiter and verified:
        result=minimize(lambda x:study.profile(tuple(x))['minimax_normalized_error'],study.start,
            jac=lambda x:study.gradient(x)[0],method='SLSQP',bounds=list(zip(study.lower,study.upper)),
            options=dict(maxiter=maxiter,ftol=1e-8))
    eligible=[r for r in study.history if r['positive_LJ']]
    save_json(Path(out)/'summary.json',dict(completed=True,derivative_verified=verified,
        elapsed_seconds=time.perf_counter()-study.began,profiles=len(study.history),
        best=min(eligible,key=lambda r:r['value']) if eligible else None,
        optimizer_executed=result is not None,optimizer_success=bool(result.success) if result is not None else None,
        message=str(result.message) if result is not None else 'audit only or derivative gate not passed',
        optimizer_iterations=int(result.nit) if result is not None else 0,
        material_accepted=False,whole_family_impossibility_proved=False,production_changed=False,
        physical_time_calibrated=False,actual_yield_validated=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(__doc__);parser.add_argument('--parent',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--maxiter',type=int,default=0)
    args=parser.parse_args();run(args.parent,args.out,args.maxiter)
