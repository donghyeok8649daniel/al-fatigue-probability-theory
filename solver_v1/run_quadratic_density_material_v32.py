"""Joint static test of the derived positive density shape; no production use."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .polynomial_exponential_density import QuadraticEnvelopeDensity
from .quadratic_density_interface import QuadraticScalarEnvironment, bulk_scalar_columns
from .coordination_screening import CoordinationScreenedBulk, screened_site_jet
from .joint_fcc_interface_calibration import MinimalConvexEmbedding
from .vector_interface_reference import _embedding_jet
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .yield_elastic_metric import cubic_to_mode_matrix
from .run_current_material_core import ROOT
from .run_finite_q_compatibility_v31 import points
from .run_source_core_reference import load_source_material
from .periodic_plane_covariance import SourceEAMBlochHessian
from .vector_material_calibration import LENGTH_M, MaterialObservation
from .core_interface_compatibility import minimax_compatibility
from .run_radial_channels_v32 import conditional_svd
from .run_vector_registry_audit import save_json
from .run_low_frequency_forcing_v29 import sha
from .run_low_stress_cyclic_diagnostic import write_csv


def interface_matrix(problem,shape,environment):
    matrix=problem.matrix(tuple(shape)).copy()
    columns=environment.bulk_columns()
    matrix[:2,2:5]=columns[:2]
    matrix[2:5,2:5]=np.linalg.solve(cubic_to_mode_matrix(),columns[2:5])
    moments_model=problem.cache.site_model(shape[0],shape[4])
    jets={}
    for i,o in enumerate(problem.raw_observations):
        if o.bulk_index is not None:
            continue
        if o.state not in jets:
            density,_=environment.site_changes(o.state)
            _,moments,_=moments_model.site_inputs(o.state)
            scalar=np.column_stack([_embedding_jet(density,1.,MinimalConvexEmbedding(*c,1.)) for c in np.eye(3)])
            vector=2*screened_site_jet(density,moments,shape[5],law='power')
            jets[o.state]=(scalar,vector)
        scalar,vector=jets[o.state]
        matrix[i,2:5]=np.asarray(o.jet_weights)@scalar
        matrix[i,6]=np.asarray(o.jet_weights)@vector
    return matrix


def run(parent,density_fit,out):
    parent,density_fit,out=map(Path,(parent,density_fit,out))
    if out.exists():
        raise FileExistsError('fresh static study required')
    previous=json.loads((parent/'summary.json').read_bytes())
    if not previous['completed']:
        raise ValueError('completed prior calibration diagnostic required')
    shape=previous['best']['shape']
    trials=json.loads((density_fit/'fits.json').read_bytes())
    best=min((p for p in trials if p['name'].startswith('quadratic_')),key=lambda p:p['squared_loss'])
    p=best['parameters']; L=LENGTH_M/1e-10
    kernel=QuadraticEnvelopeDensity(1.,p['k']*L,p['center']/L,p['width']/L)
    out.mkdir(parents=True);began=time.perf_counter()
    save_json(out/'definition.json',dict(parent_sha256=sha(parent/'summary.json'),
        auxiliary_fit_sha256=sha(density_fit/'fits.json'),parent_angular_shape=shape,
        raw_reduced_kernel=vars(kernel),amplitude_gauge='normalize full FCC density once, never during strain',
        old_scalar_decay='held implementation placeholder; no longer an independent energy parameter',
        units='eV/interface cell, eV/L0, eV/L0^2',source_density_fit_is_NOT_material_calibration=True,
        fresh_states=[(.827,.071,-.033),(.879,.267,.113),(.941,.361,.207),(1.127,.419,.139)],
        physical_time_calibrated=False,production_changed=False))
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    environment=QuadraticScalarEnvironment(kernel)
    obs=list(impose_tangents(problem.observations,problem.source_tangents))
    base=interface_matrix(problem,shape,environment)
    oldrows=[i for i,o in enumerate(obs) if o.role!='exact' and o.units=='eV/L0^2']
    source,_,binding=load_source_material();source_op=SourceEAMBlochHessian(source.source)
    bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
    low=CoordinationScreenedBulk(shape,radius=12.,law='power')
    extra=[]; operators=[];tail_data=[]
    components=((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
    for index,(role,q) in enumerate(points()):
        H=source_op.evaluate(np.array(q)*2*np.pi/source.source.geometry.lattice_constant)*L*L
        columns,tail=bulk.evaluate(q);below,_=low.evaluate(q)
        columns[2:5],tail[2:5]=bulk_scalar_columns(bulk,environment,q)
        below[2:5],_=bulk_scalar_columns(low,environment,q)
        operators.append(columns)
        tail_data.append((tail,np.linalg.norm(columns-below,axis=(1,2))))
        for i,j in components:
            w=1. if i==j else np.sqrt(2)
            extra.append(w*columns[:,i,j])
            obs.append(MaterialObservation(f'q{index}_H{i}{j}',float(w*H[i,j]),float(.05*np.linalg.norm(H)),
                'eV/L0^2','fit' if role=='fit' else 'heldout'))
    matrix=np.vstack([base,extra])
    rows=oldrows+[i for i in range(len(base),len(obs)) if obs[i].role=='fit']
    fit=minimax_compatibility(matrix,obs,rows,nonnegative=NONNEGATIVE)
    if not fit['completed']:
        raise ArithmeticError('coefficient LP failed')
    coeff=fit['coefficients'];qrows=[]
    for k,(role,q) in enumerate(points()):
        ix=np.arange(len(base)+6*k,len(base)+6*(k+1));target=np.array([obs[i].target for i in ix])
        H=np.einsum('c,cij->ij',coeff,operators[k])
        qrows.append(dict(q=q,role=role,relative_matrix_error=np.linalg.norm(fit['predictions'][ix]-target)/np.linalg.norm(target),
            minimum_eigenvalue=np.linalg.eigvalsh(H).min(),tail_bound=tail_data[k][0]@abs(coeff),
            radius_change_bound=tail_data[k][1]@abs(coeff)))
    write_csv(out/'finite_q.csv',qrows)
    write_csv(out/'residuals.csv',[dict(name=o.name,reference=o.target,prediction=pred,scale=o.scale,
        residual=(pred-o.target)/o.scale,units=o.units,exact=o.role=='exact',selected=i in rows)
        for i,(o,pred) in enumerate(zip(obs,fit['predictions']))])
    save_json(out/'profile.json',dict(**fit,identifiability=conditional_svd(matrix,obs,rows)))
    summary=dict(completed=True,elapsed_seconds=time.perf_counter()-began,source=binding,
        eta=fit['minimax_normalized_error'],positive_LJ=fit['strictly_positive_LJ'],
        maximum_excluded_q_error=max(r['relative_matrix_error'] for r in qrows if r['role']=='excluded'),
        minimum_tested_eigenvalue=min(r['minimum_eigenvalue'] for r in qrows),
        normalized_kernel=vars(environment.kernel),reference_density=environment.reference_density,
        material_accepted=False,whole_shape_family_impossibility_proved=False,production_changed=False)
    save_json(out/'summary.json',summary);print(summary)


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__)
    for name in ('parent','density-fit','out'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(a.parent,a.density_fit,a.out)
