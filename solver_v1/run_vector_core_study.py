"""Actual three-component STATIC core relaxation, no material refit or clock.

Each invocation owns a named size/ring case and retains checkpoints. A stopped
optimizer is not a converged state. Transverse cell vectors remain fixed;
the affine xz shear is stress controlled. Parameter/state provenance is checked.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import time

import numpy as np

from .nonlinear_fcc_screw import NonlinearScrewRows
from .report_nonlinear_screw_study import read_states
from .run_low_stress_cyclic_diagnostic import ROOT, build_surface, save_json, write_csv
from .run_nonlinear_screw_study import EV_J, OUT as OLD_OUT, topology_record
from .vector_fcc_rows import VectorPeriodicCell


OUT=ROOT/'results/fcc111_active_interface/vector_core_v8'


def save_state(path,field,gamma):
    write_csv(path,[dict(j=j,l=l,u_over_L0=field[j,l,0],v_over_L0=field[j,l,1],
                         w_over_L0=field[j,l,2],gamma=gamma) for j,l in np.ndindex(field.shape[:2])])


def read_state(path):
    with gzip.open(path,'rt',encoding='utf8') as stream:
        records=list(csv.DictReader(stream))
    shape=(max(int(r['j']) for r in records)+1,max(int(r['l']) for r in records)+1)
    field=np.empty(shape+(3,));gamma=float(records[0]['gamma'])
    for r in records:
        field[int(r['j']),int(r['l'])]=[float(r[k]) for k in ('u_over_L0','v_over_L0','w_over_L0')]
        if float(r['gamma'])!=gamma:
            raise ValueError('inconsistent saved affine shear')
    if len(records)!=int(np.prod(shape)):
        raise ValueError('incomplete saved vector state')
    return field,gamma


def run(args):
    if any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in args.suffix):
        raise ValueError('case suffix must be a simple filename label, not a path')
    began=time.perf_counter();surface,units,meta=build_surface(tolerance=2e-11)
    previous=json.loads((OLD_OUT/'metadata.json').read_text(encoding='utf8'))
    if previous['parameter_sha256']!=meta['parameter_sha256']:
        raise ValueError('saved scalar states have different material parameters')
    rows=NonlinearScrewRows(surface,tolerance=2e-12)
    cell=VectorPeriodicCell(rows,(args.size,args.size),ring=args.ring,tolerance=args.tolerance)
    output=OUT/f'{args.case}_n{args.size}_r{args.ring}{args.suffix}'
    if (output/'metadata.json').exists():
        raise FileExistsError('case exists; preserve it and use a new suffix, optionally --from-state')
    conversion=1e6*units.length_scale_m**3/EV_J
    if args.from_state:
        path=(OUT/args.from_state).resolve()
        if not path.is_relative_to(OUT.resolve()):
            raise ValueError('continuation state must be in this research results directory')
        source=json.loads((path.parent/'metadata.json').read_text(encoding='utf8'))
        if source['parameter_sha256']!=meta['parameter_sha256']:
            raise ValueError('checkpoint uses a different material candidate')
        field,gamma=read_state(path)
        if field.shape!=cell.shape+(3,):
            raise ValueError('checkpoint domain differs from the requested case')
        seed='explicit continuation from '+path.relative_to(OUT).as_posix()
    elif args.seed_separation is not None:
        scalar=cell.scalar.relax(cell.scalar.dipole_seed(args.seed_separation),shear_stress=0.)
        if not scalar['force_converged']:
            raise ArithmeticError('declared scalar seed failed its own force check')
        field=np.zeros(cell.shape+(3,));field[...,0]=scalar['displacement'];gamma=scalar['affine_shear']
        seed=f'actual scalar relaxation of declared seed separation {args.seed_separation:g} b'
    elif args.case=='dipole':
        scalar,gamma=read_states()[(args.size,'dipole',0)]
        field=np.zeros(cell.shape+(3,));field[...,0]=scalar
        seed='previous actual force-balanced scalar screw state, zero applied shear'
    else:
        field=np.zeros(cell.shape+(3,));gamma=0.;seed='perfect cubic crystal'
    metadata=dict(**meta,energy_normalization='eV per infinite-row repeat b*L0',
        state_space='three local displacement components; periodic transverse cell vectors fixed',
        initial_state=seed,line_invariant=True,finite_loop=False,production_registered=False,
        continuation_from_state=args.from_state,max_iterations=args.max_iterations,
        physical_seconds=False,physical_hz_enabled=False,atomic_rows_not_continuum_grid=True,
        loading='static E - tau*V_cell*gamma; V_cell=N*b*d*h*L0^3',
        size_j=args.size,size_l=args.size,transverse_ring=args.ring,
        row_tolerance=args.tolerance,force_tolerance=args.force_tolerance,
        registry_shear_columns_are_x_projection_only=True,
        x_winding_is_not_full_vector_defect_classification=True,
        imposed_shear_mpa=args.loads,mpa_to_ev_per_L0_cubed=conversion)
    save_json(output/'metadata.json',metadata)
    initial=cell.evaluate(field,gamma);histories=[];progress=[]
    initial_record=dict(energy_ev=initial['energy'],gamma=gamma,
        maximum_force=float(np.max(abs(initial['gradient']))),
        max_transverse_force=float(np.max(np.linalg.norm(initial['gradient'][...,1:],axis=-1))),
        **topology_record(cell.scalar,field[...,0]))
    save_json(output/'initial_state.json',initial_record)
    save_state(output/'initial_state.csv.gz',field,gamma)
    print(f"START {output.name}: E={initial['energy']:.10g}, max force={initial_record['maximum_force']:.5g}",flush=True)
    baseline=None
    for index,tau in enumerate(args.loads):
        tick=time.perf_counter()
        def callback(iteration,u,g):
            if iteration%20:
                return
            out=cell.evaluate(u,g)
            record=dict(step=index,tau_mpa=tau,iteration=iteration,energy_ev=out['energy'],
                max_force=float(np.max(abs(out['gradient']))),gamma=g,
                max_transverse_displacement=float(np.max(np.linalg.norm(u[...,1:],axis=-1))),
                **topology_record(cell.scalar,u[...,0]),
                elapsed_seconds=time.perf_counter()-tick)
            progress.append(record);write_csv(output/'optimizer_progress.csv',progress)
            save_state(output/'checkpoint.csv.gz',u,g)
            save_state(output/f'iteration_{index}_{iteration}.csv.gz',u,g)
            print(f"  step {index} iter {iteration}: E={out['energy']:.9g}, F={record['max_force']:.3g}, "
                  f"max transverse u={record['max_transverse_displacement']:.4g}",flush=True)
        result=cell.relax(field,gamma=gamma,shear_stress=tau*conversion,
            force_tolerance=args.force_tolerance,max_iterations=args.max_iterations,callback=callback)
        field=result['displacement'];gamma=result['affine_shear']
        shear=cell.scalar.registry_shear(field[...,0],gamma)
        if baseline is None:
            baseline=shear.copy()
        record=dict(step=index,applied_shear_mpa=tau,internal_shear_mpa=result['internal_shear_stress']/conversion,
            energy_ev_per_line_repeat=result['energy'],line_energy_j_per_m=result['energy']*EV_J/(rows.b*units.length_scale_m),
            **shear,affine_change_from_initial=gamma-baseline['affine_shear'],
            registry_change_from_initial=shear['registry_shear']-baseline['registry_shear'],
            maximum_force_residual=result['maximum_force_residual'],scaled_stress_residual=result['scaled_stress_residual'],
            force_converged=bool(result['force_converged']),iterations=result['iterations'],evaluations=result['evaluations'],
            optimizer_success=result['optimizer_success'],optimizer_message=result['optimizer_message'],
            newton_polish_steps=result['newton_polish_steps'],
            maximum_transverse_displacement=float(np.max(np.linalg.norm(field[...,1:],axis=-1))),
            minimum_transverse_radius=result['minimum_transverse_radius'],minimum_site_density=result['minimum_site_density'],
            reciprocal_modes=result['reciprocal_modes'],last_mode_envelope=result['last_mode_envelope'],
            **topology_record(cell.scalar,field[...,0]),elapsed_seconds=time.perf_counter()-tick)
        histories.append(record);write_csv(output/'stress_history.csv',histories)
        save_state(output/f'state_{index}.csv.gz',field,gamma);save_state(output/'checkpoint.csv.gz',field,gamma)
        write_csv(output/f'layer_registry_{index}.csv',cell.layer_registry(field,gamma))
        print(f"DONE step={index}, tau={tau:g}: E={result['energy']:.10g}, F={result['maximum_force_residual']:.3g}, "
              f"x-winding={record['positive_cores']}/{record['negative_cores']}, seconds={record['elapsed_seconds']:.2f}",flush=True)
        if not result['force_converged']:
            save_json(output/'summary.json',dict(completed=False,reason='actual force residual not converged',
                elapsed_seconds=time.perf_counter()-began))
            raise ArithmeticError('vector force balance incomplete; checkpoint and partial data retained')
        if args.stability and index in (0,len(args.loads)-1):
            stable=cell.minimum_curvature(field,gamma,stress_control=True,tolerance=3e-7)
            save_json(output/f'stability_{index}.json',stable)
            print(f"  Hessian min={stable['minimum_eigenvalue']:.8g}, residual={stable['eigen_residual']:.3g}",flush=True)
    save_json(output/'summary.json',dict(completed=True,all_force_balanced=True,
        elapsed_seconds=time.perf_counter()-began,actual_load_steps=len(histories),
        material_calibration_accepted=False,production_validated=False,
        scope='local vector forces at specified row cutoff; domain/tail/stability require separate audit'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size',type=int,default=24)
    parser.add_argument('--ring',type=int,default=6)
    parser.add_argument('--case',choices=['perfect','dipole'],default='dipole')
    parser.add_argument('--loads',type=float,nargs='+',default=[0.,4.,25.,50.,0.,-50.,0.])
    parser.add_argument('--force-tolerance',type=float,default=2e-7)
    parser.add_argument('--tolerance',type=float,default=2e-12)
    parser.add_argument('--max-iterations',type=int,default=500)
    parser.add_argument('--suffix',default='')
    parser.add_argument('--from-state',default=None,help='saved state relative to vector_core_v8; a NEW output suffix is required')
    parser.add_argument('--seed-separation',type=float,default=None,help='explicit diagnostic seed, not calibrated defect data')
    parser.add_argument('--stability',action='store_true')
    run(parser.parse_args())
