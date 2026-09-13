"""Validate a completed rank-four trial chosen by training eta only."""
import argparse
import json
from pathlib import Path
import numpy as np
from .coordination_screening import CoordinationScreenedInterface,CoordinationScreenedBulk
from .rank_four_environment import RankFourEnvironment,rank_four_bulk_column
from .vector_interface_reference import HESSIAN_INDICES,MishinVectorInterfaceReference
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_source_core_reference import load_source_material
from .run_finite_q_compatibility_v31 import points
from .vector_material_calibration import LENGTH_M
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def run(studies,out):
    studies,out=Path(studies),Path(out)
    if out.exists():raise FileExistsError('fresh validation output required')
    trials=[]
    for summary_path in sorted(studies.glob('*/summary.json')):
        parent=summary_path.parent;summary=json.loads(summary_path.read_bytes())
        if not (parent/'definition.json').is_file():continue
        definition=json.loads((parent/'definition.json').read_bytes())
        if not summary['completed'] or not definition.get('rank_four_energy'):continue
        for p in json.loads((parent/'profiles.json').read_bytes()):
            if p['completed'] and p['strictly_positive_LJ'] and len(p['coefficients'])==11:
                trials.append((p['minimax_normalized_error'],parent,definition,p))
    if not trials:raise ValueError('completed positive-LJ rank-four trials required')
    eta,parent,definition,profile=min(trials,key=lambda r:r[0])
    shape=definition['shape'];c=np.asarray(profile['coefficients'])
    k4=definition.get('rank_four_decay') or shape[0]
    base=CoordinationScreenedInterface(shape,c[:10],law='power')
    extra=RankFourEnvironment(k4);tight=RankFourEnvironment(k4,tolerance=2e-14)
    source,_,binding=load_source_material();L=LENGTH_M/1e-10
    ref_interface=MishinVectorInterfaceReference(source.source,L)
    ref_bulk=SourceEAMBlochHessian(source.source)
    def jet(q,model=extra):
        v=base.evaluate(q)
        return np.r_[v.energy,v.gradient,v.hessian[0],v.hessian[1,1:],v.hessian[2,2]]+c[10]*model.interface_jet(tuple(q))[0]
    rows=[]
    for state in definition['fresh_validation_states']:
        q=np.asarray(state);v=jet(q);reference=ref_interface.evaluate(q);finite=[]
        for h in (2e-5,1e-5):
            high=[jet(q+h*e) for e in np.eye(3)];low=[jet(q-h*e) for e in np.eye(3)]
            finite.append((np.array([(a[0]-b[0])/(2*h) for a,b in zip(high,low)]),
                           np.column_stack([(a[1:4]-b[1:4])/(2*h) for a,b in zip(high,low)])))
        grad,H=[(4*finite[1][i]-finite[0][i])/3 for i in (0,1)]
        rows.append(dict(state=state,energy=v[0],source_energy=reference.energy,
            gradient=v[1:4],source_gradient=reference.gradient,H=v[HESSIAN_INDICES],source_H=reference.hessian,
            relative_H_error=np.linalg.norm(v[HESSIAN_INDICES]-reference.hessian)/np.linalg.norm(reference.hessian),
            gradient_fd_error=np.max(abs(grad-v[1:4])),hessian_fd_error=np.max(abs(H-v[HESSIAN_INDICES])),
            hessian_tolerance_change=np.max(abs(jet(q,tight)[HESSIAN_INDICES]-v[HESSIAN_INDICES]))))
    coarse=CoordinationScreenedBulk(shape,radius=12.,law='power')
    fine=CoordinationScreenedBulk(shape,radius=16.,law='power');qrows=[]
    for role,q in points():
        matrices=[];bounds=[]
        for bulk in (coarse,fine):
            columns,tails=bulk.evaluate(q);H4,tail4=rank_four_bulk_column(extra,bulk,q)
            matrices.append(np.einsum('c,cij->ij',c[:10],columns)+c[10]*H4)
            bounds.append(tails@abs(c[:10])+abs(c[10])*tail4)
        reference=ref_bulk.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*L*L
        qrows.append(dict(q=q,role=role,relative_H_error=np.linalg.norm(matrices[1]-reference)/np.linalg.norm(reference),
            minimum_eigenvalue=np.linalg.eigvalsh(matrices[1]).min(),tail_bound=bounds[1],
            radius_change=np.linalg.norm(matrices[1]-matrices[0])))
    perfect=jet((base.h,0.,0.))
    out.mkdir(parents=True)
    write_csv(out/'fresh_states.csv',rows);write_csv(out/'finite_q.csv',qrows)
    save_json(out/'summary.json',dict(completed=True,selected_study=parent.name,selected_profile=profile['label'],
        selection='minimum training eta among completed positive-LJ trials; never excluded-q or fresh states',
        parent_sha256=sha(parent/'summary.json'),profile_sha256=sha(parent/'profiles.json'),source=binding,
        eta=eta,k4=k4,D4=c[10],maximum_fresh_H_error=max(r['relative_H_error'] for r in rows),
        maximum_gradient_fd_error=max(r['gradient_fd_error'] for r in rows),
        maximum_hessian_fd_error=max(r['hessian_fd_error'] for r in rows),
        maximum_hessian_tolerance_change=max(r['hessian_tolerance_change'] for r in rows),
        maximum_excluded_q_error=max(r['relative_H_error'] for r in qrows if r['role']=='excluded'),
        maximum_q_tail_bound=max(r['tail_bound'] for r in qrows),
        maximum_q_radius_change=max(r['radius_change'] for r in qrows),
        pristine_energy=perfect[0],pristine_gradient=perfect[1:4],pristine_H=perfect[HESSIAN_INDICES],
        material_accepted=False,all_q_stability_proved=False,physical_time_calibrated=False,production_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('--studies',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);args=p.parse_args();run(args.studies,args.out)
