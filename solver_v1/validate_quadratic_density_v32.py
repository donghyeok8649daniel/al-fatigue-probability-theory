"""Independent fresh-state/derivative/tail replay, never calibration approval."""
import argparse
import json
from pathlib import Path
import numpy as np
from .quadratic_density_interface import QuadraticDensityInterface,bulk_scalar_columns
from .polynomial_exponential_density import QuadraticEnvelopeDensity
from .coordination_screening import CoordinationScreenedBulk
from .run_source_core_reference import load_source_material
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_finite_q_compatibility_v31 import points
from .vector_material_calibration import LENGTH_M
from .run_vector_registry_audit import save_json
from .run_low_frequency_forcing_v29 import sha
from .run_low_stress_cyclic_diagnostic import write_csv
from .vector_interface_reference import MishinVectorInterfaceReference


def run(parent,out):
    parent,out=Path(parent),Path(out)
    if out.exists():
        raise FileExistsError('fresh independent validation directory required')
    summary=json.loads((parent/'summary.json').read_bytes())
    definition=json.loads((parent/'definition.json').read_bytes())
    if not summary['completed'] or summary['best'] is None:
        raise ValueError('completed search with an eligible candidate required')
    best=summary['best'];shape=best.get('angular_shape',definition['angular_shape_fixed']);c=np.asarray(best['coefficients'])
    kernel=QuadraticEnvelopeDensity(**best['kernel'])
    model=QuadraticDensityInterface(shape,c,kernel)
    loose=QuadraticDensityInterface(shape,c,kernel,tolerance=2e-11)
    source,_,binding=load_source_material();operator=SourceEAMBlochHessian(source.source)
    source_interface=MishinVectorInterfaceReference(source.source,LENGTH_M/1e-10)
    states=list(definition['fresh_states'])
    extra_file=parent/'additional_validation_definition.json'
    extra=None
    if extra_file.exists():
        extra=json.loads(extra_file.read_bytes())
        if extra['used_in_fit'] or extra['source_values_inspected']:
            raise ValueError('additional states must be declared before source evaluation, outside the fit')
        states+=extra['states']
        if len({tuple(q) for q in states})!=len(states):
            raise ValueError('duplicate validation states')
    rows=[]
    for index,state in enumerate(states):
        q=np.asarray(state);v=model.evaluate(q);ref=source_interface.evaluate(q);low=loose.evaluate(q)
        fd=[]
        for h in (2e-5,1e-5):
            hi=[model.evaluate(q+h*e) for e in np.eye(3)]
            lo=[model.evaluate(q-h*e) for e in np.eye(3)]
            g=np.array([(a.energy-b.energy)/(2*h) for a,b in zip(hi,lo)])
            H=np.column_stack([(a.gradient-b.gradient)/(2*h) for a,b in zip(hi,lo)])
            fd.append((g,H))
        g,H=[(4*fd[1][i]-fd[0][i])/3 for i in (0,1)]
        role=('new_declared_before_final_selection' if index>=len(definition['fresh_states'])
              else 'previously_inspected_excluded' if extra else 'new_predeclared')
        rows.append(dict(state=state,validation_role=role,energy=v.energy,source_energy=ref.energy,
            gradient=v.gradient,source_gradient=ref.gradient,hessian=v.hessian,source_hessian=ref.hessian,
            relative_H_error=np.linalg.norm(v.hessian-ref.hessian)/np.linalg.norm(ref.hessian),
            gradient_fd_error=np.max(abs(g-v.gradient)),hessian_fd_error=np.max(abs(H-v.hessian)),
            energy_tolerance_change=abs(v.energy-low.energy),
            hessian_tolerance_change=np.max(abs(v.hessian-low.hessian))))
    highbulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
    lowbulk=CoordinationScreenedBulk(shape,radius=12.,law='power')
    qrows=[];L=LENGTH_M/1e-10
    for role,q in points():
        columns,tails=highbulk.evaluate(q);below,_=lowbulk.evaluate(q)
        columns[2:5],tails[2:5]=bulk_scalar_columns(highbulk,model.scalar_environment,q)
        below[2:5],_=bulk_scalar_columns(lowbulk,model.scalar_environment,q)
        H=np.einsum('c,cij->ij',c,columns)
        ref=operator.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*L*L
        qrows.append(dict(role=role,q=q,relative_matrix_error=np.linalg.norm(H-ref)/np.linalg.norm(ref),
            minimum_eigenvalue=np.linalg.eigvalsh(H).min(),tail_bound=tails@abs(c),
            radius_change=np.linalg.norm(np.einsum('c,cij->ij',c,columns-below))))
    out.mkdir(parents=True)
    write_csv(out/'fresh_states.csv',rows);write_csv(out/'finite_q.csv',qrows)
    save_json(out/'summary.json',dict(completed=True,parent_sha256=sha(parent/'summary.json'),source=binding,
        additional_validation_sha256=sha(extra_file) if extra else None,
        fresh_vs_previously_inspected_roles_in_csv=True,
        eta=best['value'],
        maximum_validation_H_relative_error=max(r['relative_H_error'] for r in rows),
        maximum_fresh_H_relative_error=max(r['relative_H_error'] for r in rows
                                         if r['validation_role'].startswith('new_')),
        fresh_state_count=sum(r['validation_role'].startswith('new_') for r in rows),
        previously_inspected_state_count=sum(r['validation_role']=='previously_inspected_excluded' for r in rows),
        maximum_gradient_fd_error=max(r['gradient_fd_error'] for r in rows),
        maximum_hessian_fd_error=max(r['hessian_fd_error'] for r in rows),
        maximum_hessian_tolerance_change=max(r['hessian_tolerance_change'] for r in rows),
        maximum_excluded_q_error=max(r['relative_matrix_error'] for r in qrows if r['role']=='excluded'),
        minimum_tested_eigenvalue=min(r['minimum_eigenvalue'] for r in qrows),
        maximum_q_tail_bound=max(r['tail_bound'] for r in qrows),
        maximum_q_radius_change=max(r['radius_change'] for r in qrows),
        material_accepted=False,full_spectral_stability_proved=False,production_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--parent',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.parent,a.out)
