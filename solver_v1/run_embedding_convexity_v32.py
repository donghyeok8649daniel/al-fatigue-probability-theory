"""A coefficient-sign ablation is not permission to adopt unstable energies.

The source auxiliary embedding need not be globally convex. Test whether
inherited convexity alone explains the material mismatch, with LJ unchanged.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents,NONNEGATIVE
from .core_interface_compatibility import minimax_compatibility
from .coordination_screening import CoordinationScreenedBulk
from .run_current_material_core import ROOT
from .run_source_core_reference import load_source_material
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_finite_q_compatibility_v31 import points
from .vector_material_calibration import LENGTH_M,MaterialObservation
from .run_vector_registry_audit import save_json
from .run_low_frequency_forcing_v29 import sha
from .run_low_stress_cyclic_diagnostic import write_csv


def run(parent,out):
    parent,out=Path(parent),Path(out)
    if out.exists():
        raise FileExistsError('fresh sign audit required')
    previous=json.loads((parent/'summary.json').read_bytes())
    if not previous['completed'] or not previous['best']:
        raise ValueError('completed original-family shape required')
    shape=previous['best']['shape']
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    obs=list(impose_tangents(problem.observations,problem.source_tangents))
    base=problem.matrix(tuple(shape));bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
    source,_,binding=load_source_material();operator=SourceEAMBlochHessian(source.source)
    rows=[i for i,o in enumerate(obs) if o.role!='exact' and o.units=='eV/L0^2']
    extra=[];L=LENGTH_M/1e-10
    for k,(role,q) in enumerate(points()):
        if role!='fit':
            continue
        H=operator.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*L*L
        columns,_=bulk.evaluate(q)
        for i,j in ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2)):
            w=1. if i==j else np.sqrt(2)
            rows.append(len(obs));extra.append(w*columns[:,i,j])
            obs.append(MaterialObservation(f'q{k}_H{i}{j}',float(w*H[i,j]),
                float(.05*np.linalg.norm(H)),'eV/L0^2','fit'))
    matrix=np.vstack([base,extra])
    out.mkdir(parents=True)
    save_json(out/'definition.json',dict(parent_sha256=sha(parent/'summary.json'),source=binding,shape=shape,
        sign_controls=['inherited','C_embedding_signed','K3_signed','C_embedding_and_K3_signed'],
        physical_stability_not_implied_by_relaxed_signs=True,
        production_sign_constraints_changed=False))
    profiles=[]
    for name,removed in [('inherited',()),('C_embedding_signed',(4,)),('K3_signed',(9,)),
                         ('C_embedding_and_K3_signed',(4,9))]:
        fit=minimax_compatibility(matrix,obs,rows,nonnegative=tuple(i for i in NONNEGATIVE if i not in removed))
        profiles.append(dict(name=name,**fit))
        print(name,fit.get('minimax_normalized_error'),flush=True)
    source_rows=[]
    rho0=source.source._rho_bulk
    for x in (.05,.1,.2,.4,.6,.8,1.,1.1,1.2,1.3,1.5):
        source_rows.append(dict(normalized_density=x,
            embedding_second_eV=float(source.source._F(x*rho0,2))*rho0*rho0,
            auxiliary_source_gauge=True))
    write_csv(out/'source_embedding_curvature.csv',source_rows)
    save_json(out/'profiles.json',profiles)
    save_json(out/'summary.json',dict(completed=True,
        profiles=[dict(name=p['name'],completed=p['completed'],eta=p.get('minimax_normalized_error'),
            positive_LJ=p.get('strictly_positive_LJ'),dual_gap=p.get('duality_gap')) for p in profiles],
        material_accepted=False,full_shape_impossibility_proved=False,
        production_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--parent',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.parent,a.out)
