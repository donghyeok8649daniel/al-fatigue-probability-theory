"""Resolve a symmetry-connected STATIC core saddle under one fixed boundary.

Not automatically the full translation Peierls barrier. Energy is per straight
row repeat, never a finite-loop activation energy or a rate. Two descents can
verify the actual adjacent minima without assuming their identity from a plot.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np

from .core_stationary_point import stationary_core
from .core_stability import lowest_core_mode
from .current_core_diagnostics import twofold_core_partner, inner_field_difference, phase_winding_cells
from .report_current_material_core import restore_case
from .run_current_material_core import load_current_material, ROOT
from .run_source_core_reference import load_source_material, build_source_core
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--second-case', type=Path,
                        help='optional independently relaxed endpoint on EXACTLY the same fixed boundary')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--source-reference', action='store_true')
    parser.add_argument('--material', type=Path)
    parser.add_argument('--verify-descents', action='store_true')
    args = parser.parse_args()
    if args.source_reference and args.material:
        raise ValueError('source and candidate material inputs are separate')
    if args.out.exists():
        raise FileExistsError('fresh saddle directory required')
    began = time.perf_counter()
    extra = {}
    loader = load_current_material
    if args.source_reference:
        loader = load_source_material; extra['core_builder'] = build_source_core
    elif args.material:
        from functools import partial
        loader = partial(load_current_material,args.material)
    model,tensor,binding = loader()
    core,original,meta,previous,_ = restore_case(args.baseline,model,tensor,binding,**extra)
    if not previous['final_stability_probe']['stable_on_tested_fixed_boundary']:
        raise ValueError('verified stable zero-load endpoint required')
    if args.second_case:
        second_core,partner,second_meta,second_saved,_=restore_case(args.second_case,model,tensor,binding,**extra)
        if not second_saved['final_stability_probe']['stable_on_tested_fixed_boundary']:
            raise ValueError('the second endpoint must be independently stable')
        for key in ('radius_over_L0','ring','shear_traction_MPa','burgers_sign','reciprocal_tolerance'):
            if meta[key]!=second_meta[key]: raise ValueError('endpoint domain/numerics/load mismatch')
        if (not np.array_equal(core.indices,second_core.indices)
                or not np.allclose(core.boundary,second_core.boundary,atol=1e-14,rtol=0)):
            raise ValueError('the same prescribed atomic exterior is required')
    else:
        partner = twofold_core_partner(core,original)
    original_value = core.evaluate(original); partner_value = core.evaluate(partner)
    if not args.second_case and abs(original_value['energy']-partner_value['energy']) > 2e-10:
        raise ArithmeticError('twofold energy symmetry failed')
    metadata = dict(meta)
    metadata.update(baseline=args.baseline.resolve().relative_to(ROOT).as_posix(),
        source_state_sha256=hashlib.sha256((args.baseline/'state.csv').read_bytes()).hexdigest(),
        stationary_root_not_dynamics=True, full_translation_barrier_assumed=False)
    metadata.update(endpoint_kind='independently_supplied' if args.second_case else 'twofold_symmetry',
        second_endpoint=None if not args.second_case else args.second_case.resolve().relative_to(ROOT).as_posix(),
        second_state_sha256=None if not args.second_case else hashlib.sha256((args.second_case/'state.csv').read_bytes()).hexdigest())
    save_json(args.out/'metadata.json',metadata)
    def write_state(name,field):
        write_csv(args.out/name,[dict(j=int(i[0]),l=int(i[1]),ux=u[0],uy=u[1],uz=u[2])
            for i,u in zip(core.indices[core.free_ids],field)])
    def progress(iteration,field,row):
        write_state('checkpoint.csv',field)
        save_json(args.out/'progress.json',dict(completed=False,**row,elapsed_seconds=time.perf_counter()-began))
        print(f'{args.out.name}: stationary step{iteration}, force={row["maximum_force"]:.4g}',flush=True)
    result = stationary_core(core,.5*(original+partner),callback=progress)
    field = result['field']; write_state('state.csv',field)
    mode = lowest_core_mode(core,field,force_tolerance=5e-8,method='assembled')
    descents = []
    if args.verify_descents and result['force_converged'] and result['negative_eigenvalues']==1:
        for sign in (-1,1):
            tested = core.relax(field+sign*.02*result['minimum_mode'],method='newton_cg',
                max_iterations=120,force_tolerance=5e-7)
            curvature = lowest_core_mode(core,tested['field'],method='assembled',force_tolerance=5e-7)
            write_state(f'descent_{sign:+d}.csv',tested['field'])
            descents.append(dict(sign=sign,energy=tested['energy'],maximum_force=tested['maximum_free_force'],
                stable_fixed_boundary=curvature['stable_on_tested_fixed_boundary'],
                minimum_H=curvature['minimum_eigenvalue'],
                difference_to_original=inner_field_difference(core,original,tested['field'],radius=min(2.,core.free_radius)),
                difference_to_partner=inner_field_difference(core,partner,tested['field'],radius=min(2.,core.free_radius))))
    save_json(args.out/'projected_circulation.json',phase_winding_cells(core,field))
    write_csv(args.out/'stationary_history.csv',result['history'])
    write_csv(args.out/'hessian_eigenvalues.csv',[dict(index=i,eigenvalue=value) for i,value in enumerate(result['eigenvalues'])])
    save_json(args.out/'saddle_summary.json',dict(completed=True,force_converged=result['force_converged'],
        maximum_force=result['maximum_force'],negative_eigenvalues=result['negative_eigenvalues'],
        unresolved_eigenvalues=result['unresolved_eigenvalues'],hessian_resolution_floor=result['hessian_resolution_floor'],
        saddle_energy=result['energy'],original_energy=original_value['energy'],partner_energy=partner_value['energy'],
        barrier_eV_per_straight_row_repeat=result['energy']-original_value['energy'],
        stationary_root_is_index_one=bool(result['force_converged'] and result['negative_eigenvalues']==1),
        endpoint_kind=metadata['endpoint_kind'],partner_means_second_reference_endpoint=True,
        minimum_mode_probe={k:v for k,v in mode.items() if k!='eigenvector'},descents=descents,
        full_translation_Peierls_validated=False,finite_loop_activation_validated=False,
        material_accepted=False,physical_yield_validated=False,physical_Hz=False,
        elapsed_seconds=time.perf_counter()-began))
    print(f'{args.out.name}: index={result["negative_eigenvalues"]}, '
        f'barrier={result["energy"]-original_value["energy"]:.8g} eV/row repeat',flush=True)


if __name__=='__main__':
    main()
