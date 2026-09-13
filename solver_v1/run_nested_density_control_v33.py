"""Nested pure-exponential limit of the SAME positive quadratic family.

Mix densities before each site's nonlinear energy, never fitted energies.
Both endpoint kernels have pristine bulk density one at the same decay.
For eta>0 this is exactly a quadratic-envelope kernel, not another family.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .polynomial_exponential_density import QuadraticEnvelopeDensity
from .quadratic_density_interface import QuadraticScalarEnvironment,bulk_scalar_columns
from .run_quadratic_density_material_v32 import interface_matrix
from .coordination_screening import CoordinationScreenedBulk
from .tail_constrained_material import normalized_amplitude
from .interface_tangent_calibration import TangentCalibrationProblem,NONNEGATIVE
from .core_interface_compatibility import minimax_compatibility
from .vector_material_calibration import MaterialObservation
from .run_current_material_core import ROOT
from .run_finite_q_compatibility_v31 import points
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json


def run(parent,out):
    parent,out=Path(parent),Path(out)
    if out.exists():raise FileExistsError('fresh nested-density control required')
    definition=json.loads((parent/'definition.json').read_bytes())
    saved=json.loads((parent/'profiles.json').read_bytes())[0]
    obs=[MaterialObservation(**o) for o in json.loads((parent/'observations.json').read_bytes())]
    shape=tuple(definition['shape']);k=shape[0]
    source_shape=json.loads((ROOT/'results/radial_channels_v32/quadratic_density_material/summary.json').read_bytes())['normalized_kernel']
    centers=(0.,source_shape['center']);width=source_shape['width']
    etas=(.0001,.001,.01,.1,.3,1.)
    out.mkdir(parents=True)
    save_json(out/'definition.json',dict(parent_sha256=sha(parent/'profiles.json'),shape=shape,
        centers=centers,width=width,etas=etas,eta_zero='exact existing exponential model',
        definition='f_eta=(1-eta) f_exp + eta f_quad, each normalized to bulk density one, SAME k',
        equivalent_width_squared='w^2+(1-eta) C_exp/(eta C_quad)',
        purpose='finite-width search omitted exact exponential closure; bounded continuation control',
        auxiliary_density_in_loss=False,excluded_q_in_loss=False,new_blind_validation=False,
        production_changed=False))
    start=time.perf_counter()
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    base=problem.matrix(shape);bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
    components=((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
    weights=np.array([1.,1.,1.,np.sqrt(2),np.sqrt(2),np.sqrt(2)])
    qcolumns=[bulk.evaluate(q)[0] for _,q in points()]
    def append_q(matrix,columns):
        return np.vstack([matrix,*[weights[:,None]*np.array([v[:,i,j] for i,j in components]) for v in columns]])
    original=append_q(base,qcolumns)
    replay=float(np.max(abs(original@saved['coefficients']-saved['predictions'])))
    if replay>2e-8:raise ArithmeticError('exact exponential baseline replay failed')
    profiles=[]
    def fit(label,matrix,metadata):
        result=minimax_compatibility(matrix,obs,saved['tested_rows'],nonnegative=NONNEGATIVE)
        if not result['completed']:raise ArithmeticError('coefficient LP failed')
        excluded=[]
        for j,(role,_) in enumerate(points()):
            if role=='excluded':
                ix=np.arange(len(base)+6*j,len(base)+6*(j+1))
                target=np.array([obs[i].target for i in ix])
                excluded.append(float(np.linalg.norm(result['predictions'][ix]-target)/np.linalg.norm(target)))
        profiles.append(dict(label=label,**metadata,**result,excluded_q_errors=excluded))
        save_json(out/'profiles.json',profiles)
        print(label,result['minimax_normalized_error'],result['strictly_positive_LJ'],flush=True)
    fit('exact_exponential',original,dict(eta=0.,baseline_replay_error=replay))
    for center in centers:
        endpoint=QuadraticScalarEnvironment(QuadraticEnvelopeDensity(1.,k,center,width)).kernel
        for eta in etas:
            mixed_width=np.sqrt(width**2+(1-eta)*normalized_amplitude(k)/(eta*endpoint.C))
            kernel=QuadraticEnvelopeDensity(eta*endpoint.C,k,center,mixed_width)
            environment=QuadraticScalarEnvironment(kernel)
            radii=np.array([.4,.8,1.,1.7,3.,6.])
            expected=(1-eta)*normalized_amplitude(k)*np.exp(-k*radii)+eta*endpoint.radial(radii)
            mix_error=float(np.max(abs(environment.kernel.radial(radii)-expected)))
            if mix_error>2e-10:raise ArithmeticError('density mixture/normalization identity failed')
            matrix=interface_matrix(problem,shape,environment);columns=[]
            for (_,q),old in zip(points(),qcolumns):
                value=old.copy();value[2:5],_=bulk_scalar_columns(bulk,environment,q);columns.append(value)
            design=append_q(matrix,columns)
            fit(f'center{center:g}_eta{eta:g}',design,dict(eta=eta,center=center,
                kernel=vars(environment.kernel),density_mixture_error=mix_error,
                maximum_baseline_matrix_change=float(np.max(abs(design-original)))))
    save_json(out/'summary.json',dict(completed=True,elapsed_seconds=time.perf_counter()-start,
        baseline_replay_error=replay,profiles=[dict(label=p['label'],eta_mixture=p['eta'],
            minimax_error=p['minimax_normalized_error'],positive_LJ=p['strictly_positive_LJ'],
            excluded_q_error=max(p['excluded_q_errors'])) for p in profiles],
        material_accepted=False,production_changed=False,physical_time_calibrated=False,
        interpretation='same quadratic family plus exact exponential endpoint; fixed parent angular shape only'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(__doc__);parser.add_argument('--parent',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();run(args.parent,args.out)
