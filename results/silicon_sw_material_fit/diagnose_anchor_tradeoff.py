"""Replay fixed-shape anchor conflicts using saved matrices only.

No geometry evaluation, material refit, or parameter selection. A full-rank
three-amplitude matrix fixes the amplitudes at this shape; this is not a
structural identifiability or impossibility proof for other energy families.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def main(args):
    if args.output.exists():
        raise ValueError('fresh diagnostic output required')
    root = args.package
    read = lambda name: json.loads((root / name).read_text(encoding='utf-8'))
    selection = read('final_selection.json')
    chosen_name = 'raw_runs/' + selection['selected_run'] + '/selected_fit.json'
    chosen = read(chosen_name)
    groups_name = 'raw_runs/force_span_diagnostic_recovery_results/group_force_bounds.json'
    groups = read(groups_name)
    reference = read('primary_sources/test-results/model-CASTEP_ASE-test-bulk_diamond-properties.json')
    matrix = np.array(chosen['anchor']['anchor_matrix'])
    singular = np.linalg.svd(matrix, compute_uv=False)
    if singular[-1] <= singular[0] * 1e-11:
        raise RuntimeError('fixed-shape anchor matrix is numerically unresolved')
    volume = reference['diamond_a0']**3 / 4
    scale = 160.2176634 / volume
    target = np.array([0., reference['diamond_c11']/scale, reference['diamond_c12']/scale])
    solved = np.linalg.solve(matrix, target)
    locked = np.array(chosen['amplitudes'])
    if np.max(abs(solved-locked)) > 2e-11 or np.max(abs(matrix @ locked-target)) > 2e-11:
        raise RuntimeError('saved anchor does not reproduce locked amplitudes')
    rows = []
    for group in groups:
        if np.max(abs(np.array(group['locked_amplitudes'])-locked)) > 2e-12:
            raise RuntimeError('force diagnostic uses a different candidate')
        for kind in ['unconstrained', 'nonnegative']:
            alpha = np.array(group[kind+'_force_only_amplitudes'])
            residual = matrix @ alpha-target
            rows.append(dict(config_type=group['config_type'],declared_xc=group['declared_xc'],
                solution=kind,amplitudes=alpha.tolist(),
                force_RMSE_eV_A=group[kind+'_force_only_RMSE_eV_A'],
                anchor_residual_eV=residual.tolist(),
                normal_stress_at_reference_geometry_GPa=float(residual[0]*scale),
                unrelaxed_normal_curvature_residual_GPa=(residual[1:]*scale).tolist(),
                exact_normal_anchors_preserved=bool(np.max(abs(residual)) <= 2e-11)))
    sources = [chosen_name,groups_name,
               'primary_sources/test-results/model-CASTEP_ASE-test-bulk_diamond-properties.json',
               'diagnose_anchor_tradeoff.py']
    result = dict(complete=True,selected_run=selection['selected_run'],selected_profile=chosen['profile'],
        shape_fixed=True,anchor_matrix=matrix.tolist(),anchor_target_eV=target.tolist(),
        anchor_row_order=['dE/dexx','d2E/dexx2','d2E/dexxdeyy'],
        anchor_singular_values_eV=singular.tolist(),anchor_condition_number=float(singular[0]/singular[-1]),
        solved_locked_amplitudes=solved.tolist(),locked_anchor_residual_eV=(matrix@locked-target).tolist(),
        rows=rows,inputs=[dict(path=n,sha256=hashlib.sha256((root/n).read_bytes()).hexdigest()) for n in sources],
        new_geometry_evaluations=0,new_DFT=0,new_MD=0,selected_parameters_changed=False,
        material_approved=False,family_impossibility_proved=False,
        interpretation='fixed-shape exact normal anchors uniquely fix three amplitudes; force-only alternatives break them; unrelaxed curvature changes are not re-equilibrated elastic constants')
    args.output.write_bytes((json.dumps(result,indent=2,allow_nan=False)+'\n').encode())
    print(json.dumps(dict(complete=True,rows=len(rows),anchor_condition_number=result['anchor_condition_number'],
                         all_force_only_alternatives_break_anchors=all(not r['exact_normal_anchors_preserved'] for r in rows))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
