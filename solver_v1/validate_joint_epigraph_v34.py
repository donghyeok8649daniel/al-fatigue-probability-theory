"""Excluded finite-q validation of a completed joint shape study."""
import argparse
import json
from pathlib import Path
import numpy as np
from .coordination_screening import CoordinationScreenedBulk
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_source_core_reference import load_source_material
from .run_finite_q_compatibility_v31 import points
from .run_joint_shape_v31 import decode_shape
from .vector_material_calibration import LENGTH_M
from .run_vector_registry_audit import save_json


def run(parent,out):
    parent=Path(parent);out=Path(out)
    if out.exists():raise FileExistsError('fresh validation output required')
    report=json.loads((parent/'joint_summary.json').read_bytes())
    if not report['completed']:raise ValueError('completed calibration required')
    profile=report['final_profile']
    if not profile['strictly_positive_LJ']:raise ValueError('zero LJ boundary is not eligible')
    shape=decode_shape(report['final_shape_coordinates'],1.,True)
    c=np.array(profile['coefficients'])
    source,_,binding=load_source_material()
    reference=SourceEAMBlochHessian(source.source)
    bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
    coarse=CoordinationScreenedBulk(shape,radius=12.,law='power')
    rows=[]
    for role,q in points():
        columns,tail=bulk.evaluate(q);low,_=coarse.evaluate(q)
        H=np.einsum('c,cij->ij',c,columns)
        R=reference.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
        rows.append(dict(role=role,q=q,relative_error=np.linalg.norm(H-R)/np.linalg.norm(R),
            minimum_eigenvalue=np.linalg.eigvalsh(H).min(),tail_bound=tail@abs(c),
            radius_change=np.linalg.norm(np.einsum('c,cij->ij',c,columns-low))))
    out.mkdir(parents=True)
    save_json(out/'summary.json',dict(completed=True,source=binding,shape=shape,
        training_eta=profile['minimax_normalized_error'],finite_q=rows,
        maximum_excluded_q_error=max(r['relative_error'] for r in rows if r['role']=='excluded'),
        validation_status='previously inspected excluded q, not newly blind data',
        material_accepted=False,physical_time_calibrated=False,production_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('--parent',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.parent,a.out)
