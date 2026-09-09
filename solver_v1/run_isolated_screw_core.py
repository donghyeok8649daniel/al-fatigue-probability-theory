"""Execute one isolated static core; preserve checkpoints and numerical status."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .even_moment_calibration import even_basis
from .full_fcc_calibration_audit import cubic_constants_gpa
from .core_stability import lowest_core_mode
from .isolated_screw_core import IsolatedScrewCore,ScrewFarField
from .nonlocal_interface_elasticity import cubic_elastic_tensor,rotate_elastic_tensor
from .range_resolved_material import build_range_surface
from .run_low_stress_cyclic_diagnostic import ROOT,FIT,build_surface,write_csv
from .run_vector_registry_audit import save_json


OUT=ROOT/'results/fcc111_active_interface/range_core_v12/isolated_core'
EV_J=1.602176634e-19


def material(name):
    if name=='historical':
        surface,units,meta=build_surface(tolerance=2e-11)
        fit=json.loads(FIT.read_bytes())['angular_monotone_opening']
        moduli=cubic_constants_gpa(even_basis(fit['scalar_decay'],fit['angular_decay'],convex=True)[:5]@fit['coefficients'])
        length=units.length_scale_m
    else:
        path=OUT.parent/'material/calibration.json'
        fit=json.loads(path.read_bytes())['best']['cubic_5pct']
        model=build_range_surface(fit['scalar_decay'],fit['odd_decay'],fit['quadrupole_decay'],fit['coefficients'])
        surface=model.surface; length=4.05/np.sqrt(2)*1e-10
        moduli=dict(zip(['C11_GPa','C12_GPa','C44_GPa'],fit['cubic_GPa']))
        meta=dict(model='range_core_v12_cubic_5pct_unadopted',parameter_source=path.relative_to(ROOT).as_posix(),
                  parameter_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),material_calibration_accepted=False)
    cubic=cubic_elastic_tensor(*[moduli[f'C{i}_GPa'] for i in (11,12,44)])
    plane=surface.interface.bulk.geometry.plane_basis_in_stacked_cubic_axes()
    conversion=1e9*length**3/EV_J
    tensor=rotate_elastic_tensor(cubic,plane)*conversion
    # One spelling only: Windows JSON consumers can be case-insensitive.
    meta.pop('physical_hz',None)
    meta.update(moduli,length_scale_m=length,elastic_tensor_unit='eV/L0^3',
                elastic_tensor_same_candidate=True,energy_unit='eV per straight-row repeat b*L0',
                physical_time=False,physical_Hz=False,experimental_yield_validated=False)
    return surface,tensor,meta


def execute(args):
    started=time.perf_counter()
    if not args.label or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in args.label):
        raise ValueError('simple case label required')
    out=OUT/args.label
    if out.exists():
        raise FileExistsError('preserve existing case; use a fresh label')
    surface,tensor,meta=material(args.model); b=surface.interface.bulk.geometry.b
    center=(np.sqrt(3)*b/12,surface.h/2)
    far=ScrewFarField(tensor,b,center)
    traction=args.shear_mpa*1e6*meta['length_scale_m']**3/EV_J
    core=IsolatedScrewCore(surface,far,free_radius=args.radius,ring=args.ring,
                          tolerance=args.tolerance,shear_traction=traction)
    seed=core.initial
    if args.from_case:
        if any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in args.from_case):
            raise ValueError('simple source-case label required')
        origin=OUT/args.from_case
        prior=json.loads((origin/'metadata.json').read_bytes())
        if prior['parameter_sha256']!=meta['parameter_sha256']:
            raise ValueError('continuation cannot switch the material potential')
        if prior['radius_over_L0']!=args.radius:
            raise ValueError('continuation must preserve the free disk')
        with (origin/'state.csv').open(newline='') as stream:
            saved={(int(r['j']),int(r['l'])):np.array([float(r[k]) for k in ('ux','uy','uz')]) for r in csv.DictReader(stream)}
        seed=np.array([saved[tuple(i)] for i in core.indices[core.free_ids]])
        delta=(args.shear_mpa-prior['shear_traction_MPa'])*1e6*meta['length_scale_m']**3/EV_J
        seed+=far.displacement(core.xyz[core.free_ids],shear_traction=delta,burgers_sign=0)
        meta.update(continued_from=args.from_case,
            source_state_sha256=hashlib.sha256((origin/'state.csv').read_bytes()).hexdigest())
    if args.escape_sign:
        if not args.from_case or not np.isfinite(args.escape_amplitude) or args.escape_amplitude<=0:
            raise ValueError('negative-mode exploration requires a saved state and positive declared perturbation')
        mode=lowest_core_mode(core,seed,force_tolerance=args.force_tolerance)
        if mode['minimum_eigenvalue']>=-max(mode['eigenpair_residual'],1e-9):
            raise ValueError('no resolved negative mode; do not manufacture a mode escape')
        if not all(c['curvature']<0 for c in mode['energy_curvature_checks']):
            raise ValueError('independent energy differences do not verify the negative mode')
        meta.update(negative_mode_probe={k:v for k,v in mode.items() if k!='eigenvector'},
                    escape_sign=args.escape_sign,escape_amplitude_over_L0=args.escape_amplitude)
        seed+=args.escape_sign*args.escape_amplitude*mode['eigenvector']
    meta.update(radius_over_L0=args.radius,ring=args.ring,shear_traction_MPa=args.shear_mpa,
        center_over_L0=center,neighbor_count=len(core.offsets),energy_sites=len(core.energy_ids),
        free_sites=len(core.free_ids),outer_sites=len(core.indices)-len(core.free_ids),
        reciprocal_tolerance=args.tolerance,force_tolerance_eV_L0=args.force_tolerance,
        farfield_log_energy_coefficient_eV_per_row=far.log_energy_coefficient,
        boundary='fixed anisotropic single-screw Volterra + uniform same-material traction',
        boundary_is_not_measured_pinning=True,transverse_infinite_tail_certified=False)
    save_json(out/'metadata.json',meta)
    def checkpoint(iteration,field):
        if iteration==1 or iteration%10==0:
            rows=[dict(j=int(i[0]),l=int(i[1]),ux=u[0],uy=u[1],uz=u[2])
                  for i,u in zip(core.indices[core.free_ids],field)]
            write_csv(out/'checkpoint.csv',rows)
            save_json(out/'progress.json',dict(completed=False,iteration=iteration,elapsed_seconds=time.perf_counter()-started))
            print(f'{args.label}: iteration {iteration}, elapsed {time.perf_counter()-started:.1f}s',flush=True)
    initial=core.evaluate(seed)
    result=core.relax(seed,max_iterations=args.iterations,force_tolerance=args.force_tolerance,callback=checkpoint)
    final=core.evaluate(result['field'])
    # Mandatory for new runs: a force-converged symmetric saddle is NOT a
    # stable core. Preserve it as evidence and expose the negative mode.
    mode=lowest_core_mode(core,result['field'],force_tolerance=args.force_tolerance)
    result['final_stability_probe']={k:v for k,v in mode.items() if k!='eigenvector'}
    write_csv(out/'state.csv',[dict(j=int(i[0]),l=int(i[1]),ux=u[0],uy=u[1],uz=u[2],
        gx=g[0],gy=g[1],gz=g[2]) for i,u,g in zip(core.indices[core.free_ids],result['field'],final['gradient'])])
    write_csv(out/'radial_energy.csv',core.radial_energy(result['field'],np.arange(1.,args.radius+.01,.5)))
    save_json(out/'summary.json',dict(**{k:v for k,v in result.items() if k!='field'},
        initial_energy=initial['energy'],initial_maximum_free_force=float(np.max(abs(initial['gradient']))),
        energy_decrease=initial['energy']-result['energy'],
        maximum_transverse_relaxation=float(np.max(np.linalg.norm(result['field'][:,1:],axis=1))),
        reciprocal_modes=final['reciprocal_modes'],last_mode_envelope=final['last_mode_envelope'],
        elapsed_seconds=time.perf_counter()-started,completed=True,
        domain_converged=False,physical_yield_validated=False,
        core_energy_matched_to_outer_elasticity=False))
    save_json(out/'progress.json',dict(completed=True,elapsed_seconds=time.perf_counter()-started))
    print(f"{args.label}: complete force={result['maximum_free_force']:.4g}, winding={result['boundary_winding']:.4g}",flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--label',required=True)
    p.add_argument('--model',choices=['historical','range'],default='historical')
    p.add_argument('--radius',type=float,default=3.)
    p.add_argument('--ring',type=int,default=3)
    p.add_argument('--shear-mpa',type=float,default=0.)
    p.add_argument('--from-case',help='same-material/same-disk saved state for static continuation')
    p.add_argument('--escape-sign',type=int,choices=[-1,0,1],default=0,
                   help='explicit two-sided negative-Hessian-mode exploration; no force/energy change')
    p.add_argument('--escape-amplitude',type=float,default=.02)
    p.add_argument('--iterations',type=int,default=150)
    p.add_argument('--force-tolerance',type=float,default=2e-6)
    p.add_argument('--tolerance',type=float,default=2e-12)
    execute(p.parse_args())


if __name__=='__main__':
    main()
