"""Combine actual v22 studies without promoting research data to Al yield.

Static load-return is not a time hold. Reconstruction is not full translation.
Profile spans interpolate a declared sampled registry diagnostic; they are not
fitted continuum core radii, partial separations, or calibration targets.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .current_core_diagnostics import twofold_core_partner,inner_field_difference
from .report_current_material_core import restore_case,compare_states
from .run_current_material_core import ROOT,load_current_material
from .run_source_core_reference import load_source_material,build_source_core
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def read_csv(path):
    with path.open(encoding='utf8',newline='') as stream:
        return list(csv.DictReader(stream))


def sampled_span(rows,b):
    selected=sorted((r for r in rows if r['both_rows_free']=='True'),key=lambda r:float(r['y_over_L0']))
    y=np.array([float(r['y_over_L0']) for r in selected])
    slip=np.array([float(r['affine_subtracted_delta_x'])/b for r in selected])
    crossings=[]; brackets=[]
    for level in (.75,.25):
        ids=np.flatnonzero((slip[:-1]-level)*(slip[1:]-level)<0)
        if len(ids)!=1:
            return dict(unique_quarter_crossings=False,span_over_L0=None)
        i=ids[0]
        crossings.append(float(y[i]+(y[i+1]-y[i])*(level-slip[i])/(slip[i+1]-slip[i])))
        brackets.append((float(y[i]),float(y[i+1])))
    return dict(unique_quarter_crossings=True,span_over_L0=crossings[1]-crossings[0],
        crossing_75_over_L0=crossings[0],crossing_25_over_L0=crossings[1],
        bracket_75_lower=brackets[0][0],bracket_75_upper=brackets[0][1],
        bracket_25_lower=brackets[1][0],bracket_25_upper=brackets[1][1])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT/'results/current_material_core_v22')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--trial',type=Path)
    parser.add_argument('--wider',type=Path)
    args=parser.parse_args(); root=args.root
    if args.out.exists(): raise FileExistsError('fresh combined report required')
    candidate=load_current_material(); source=load_source_material()
    trial=load_current_material(args.trial) if args.trial else None
    wider=load_current_material(args.wider) if args.wider else None
    materials={'current':candidate,'Mishin_target':source,'joint_trial':trial,'wider_trial':wider}
    stress=[]; widths=[]; barriers=[]
    cases=[('current',root/'stress_R8r5',candidate,{}),
           ('Mishin_target',root/'source_reference/stress_R8r5',source,dict(core_builder=build_source_core)),
           ('Mishin_target',root/'source_reference/stress_R10r5',source,dict(core_builder=build_source_core))]
    for size in (16,20):
        if (root/f'source_reference/stress_R{size}r5/completion.json').exists():
            cases.append(('Mishin_target',root/f'source_reference/stress_R{size}r5',source,dict(core_builder=build_source_core)))
    if trial:
        cases.append(('joint_trial',root/'probe_stress_R8r5',trial,{}))
    if wider:
        cases.append(('wider_trial',root/'wider_stress_R8r5',wider,{}))
    for label,path,material,extra in cases:
        protocol=json.loads((path/'protocol.json').read_bytes())
        complete=json.loads((path/'completion.json').read_bytes())
        baseline=ROOT/protocol['baseline']
        base,zero,_,_,_=restore_case(baseline,*material,**extra)
        partner=twofold_core_partner(base,zero)
        for run in complete['records']:
            if not run['executed']: continue
            core,field,meta,result,_=restore_case(path/run['case'],*material,**extra)
            affine=core.far_field.displacement(core.xyz,shear_traction=core.shear_traction,burgers_sign=0)
            defect=field-affine[core.free_ids]
            original_error=inner_field_difference(base,zero,defect,radius=2.)
            partner_error=inner_field_difference(base,partner,defect,radius=2.)
            projected=json.loads((path/run['case']/'projected_circulation.json').read_bytes())
            stress.append(dict(model=label,case=run['case'],radius=core.free_radius,ring=core.ring,
                shear_MPa=meta['shear_traction_MPa'],stable=run['stable_fixed_boundary'],
                force_residual_eV_L0=result['maximum_free_force'],
                minimum_H=result['final_stability_probe']['minimum_eigenvalue'],
                defect_difference_to_original=original_error['maximum_vector_change_over_L0'],
                defect_difference_to_twofold=partner_error['maximum_vector_change_over_L0'],
                closer_zero_core='original' if original_error['rms_vector_change_over_L0']<partner_error['rms_vector_change_over_L0'] else 'twofold',
                winding=projected['total_winding'],
                projected_cell_y=projected['cells'][0]['y_reference'] if len(projected['cells'])==1 else None,
                full_lattice_glide_verified=False,physical_time_hold=False,macroscopic_residual_strain=None))
    plot_profiles=[]
    profile_groups=[('current',root,['stable_plus_R8r5','stable_plus_R10r7']),
            ('Mishin_target',root/'source_reference',['escaped_plus_R8r5','stable_plus_R10r5','stable_plus_R16r5'])]
    if (root/'source_reference/stable_plus_R20r5/summary.json').exists():
        profile_groups[1][2].append('stable_plus_R20r5')
    if trial:
        profile_groups.append(('joint_trial',root,['probe_core_plus_R8r5','probe_stationary_R10r7']))
    if wider:
        profile_groups.append(('wider_trial',root,['wider_core_near_R8r5','wider_core_R10r7']))
    for label,parent,names in profile_groups:
        for name in names:
            meta=json.loads((parent/name/'metadata.json').read_bytes())
            profile=read_csv(parent/name/'adjacent_registry.csv')
            widths.append(dict(model=label,case=name,radius=meta['radius_over_L0'],ring=meta['ring'],
                **sampled_span(profile,meta['b_reduced']),
                physical_core_radius=False,used_for_fit=False))
            plot_profiles.append((label,name,profile))
        for path in sorted(parent.glob('reconstruction_*')):
            if not (path/'saddle_summary.json').exists(): continue
            result=json.loads((path/'saddle_summary.json').read_bytes())
            meta=json.loads((path/'metadata.json').read_bytes())
            binding=materials[label][2]['parameter_sha256']
            if meta['parameter_sha256']!=binding:
                continue  # no original-material saddle assigned to the trial
            barriers.append(dict(model=label,case=path.name,radius=meta['radius_over_L0'],ring=meta['ring'],
                barrier_eV_per_straight_row_repeat=result['barrier_eV_per_straight_row_repeat'],
                force_residual=result['maximum_force'],negative_eigenvalues=result['negative_eigenvalues'],
                minimum_H=result['minimum_mode_probe']['minimum_eigenvalue'],
                two_descents_verified=len(result['descents'])==2 and all(x['stable_fixed_boundary'] for x in result['descents']),
                full_translation_Peierls=False,finite_loop_activation=False))
    # Inspect reference tradeoffs, including the unchanged baseline loss.
    rows=read_csv(root/'constrained_force_profile/static_tradeoff.csv')
    fitted=[r for r in rows if r['profile']=='sign_only' and r['previous_role']=='fit']
    baseline_loss=sum(((float(r['baseline'])-float(r['target']))/float(r['scale']))**2 for r in fitted)
    constrained_loss=sum(float(r['normalized_residual'])**2 for r in fitted)
    write_csv(args.out/'signed_static_stress_return.csv',stress)
    write_csv(args.out/'sampled_registry_spans.csv',widths)
    write_csv(args.out/'reconstruction_barriers.csv',barriers)
    comparisons=[]
    checks=[('source_independent_translation',root/'source_reference/translated_plus_R16r5',
             root/'source_reference/stress_R16r5/positive_3_0MPa',source,dict(core_builder=build_source_core))]
    if (root/'source_reference/translated_j2_R20r5/summary.json').exists():
        checks.append(('source_independent_two_row_translation',root/'source_reference/translated_j2_R20r5',
             root/'source_reference/stress_R20r5/positive_3_0MPa',source,dict(core_builder=build_source_core)))
    if trial:
        checks.extend([
            ('trial_independent_stationary_solver',root/'probe_core_plus_R10r7',root/'probe_stationary_R10r7',trial,{}),
            ('trial_combined_domain_and_ring',root/'probe_core_plus_R8r5',root/'probe_stationary_R10r7',trial,{})])
    for name in ('source_seed_current_lbfgs_R8r5','source_seed_current_recovered_R8r5'):
        if (root/name/'summary.json').exists():
            checks.append(('original_source_coordinate_multistart',root/'stable_plus_R8r5',
                           root/name,candidate,{}))
    if trial and (root/'source_seed_trial_R8r5/summary.json').exists():
        checks.append(('trial_source_coordinate_multistart',root/'probe_core_plus_R8r5',
                       root/'source_seed_trial_R8r5',trial,{}))
    if wider:
        checks.append(('wider_combined_domain_and_ring',root/'wider_core_near_R8r5',
                       root/'wider_core_R10r7',wider,{}))
        if (root/'wider_core_R8r7/summary.json').exists():
            checks.extend([
                ('wider_ring_only',root/'wider_core_near_R8r5',root/'wider_core_R8r7',wider,{}),
                ('wider_domain_only',root/'wider_core_R8r7',root/'wider_core_R10r7',wider,{})])
        for name in ('wider_core_lbfgs_R8r5','wider_core_recovered_R8r5'):
            saved=root/name/'summary.json'
            if saved.exists() and json.loads(saved.read_bytes())['final_stability_probe']['stable_on_tested_fixed_boundary']:
                checks.append(('wider_source_coordinate_multistart',root/'wider_core_near_R8r5',
                               saved.parent,wider,{}))
    for kind,first,second,material,extra in checks:
        a,_,ma,sa,da=restore_case(first,*material,**extra)
        _,_,mb,sb,db=restore_case(second,*material,**extra)
        same_operator=all(ma[k]==mb[k] for k in ('radius_over_L0','ring',
            'shear_traction_MPa','burgers_sign','reciprocal_tolerance'))
        for radius in (1.,2.,4.):
            comparisons.append(dict(comparison=kind,first=first.relative_to(ROOT).as_posix(),
                second=second.relative_to(ROOT).as_posix(),energy_difference_eV=sb['energy']-sa['energy'],
                energy_difference_uses_identical_operator=same_operator,
                energy_difference_is_a_transition_barrier=False,
                **compare_states((ma,da),(mb,db),radius=radius,b=a.rows.b,h=a.rows.h),
                infinite_domain_certified=False,macroscopic_residual_plasticity=False))
    write_csv(args.out/'independent_state_checks.csv',comparisons)
    path=root/'source_reference/endpoint_saddle_R16r5/saddle_summary.json'
    endpoint=json.loads(path.read_bytes())
    save_json(args.out/'finite_boundary_translation_connection.json',dict(
        source_only=True,saddle=endpoint,independent_translation_checks=comparisons[:3],
        source_path=path.relative_to(ROOT).as_posix(),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        experimental_yield=False,infinite_domain_Peierls=False,finite_loop_activation=False))
    save_json(args.out/'force_tradeoff.json',dict(previous_other_interface_loss=baseline_loss,
        core_only_constrained_other_interface_loss=constrained_loss,
        interface_error_growth_factor=constrained_loss/baseline_loss,
        statistical_chi_square=False,parameters_adopted=False))
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.2))
    for label,name,rows in plot_profiles:
        p=[r for r in rows if r['both_rows_free']=='True' and abs(float(r['y_over_L0']))<=4.]
        axes[0].plot([float(r['y_over_L0']) for r in p],
            [float(r['affine_subtracted_delta_x']) for r in p],'.-',label=label+' '+name)
    for label in ('current','Mishin_target','wider_trial'):
        rows=sorted((r for r in barriers if r['model']==label and r['ring']==5),key=lambda r:r['radius'])
        if not rows:continue
        axes[1].plot([r['radius'] for r in rows],
            [1e3*r['barrier_eV_per_straight_row_repeat'] for r in rows],'.-',label=label)
    axes[0].set(xlabel='Adjacent-row midpoint y/L0',ylabel='Screw disregistry / L0',
                title='Same physical FCC rows; no fitted core width')
    axes[1].set(xlabel='Free radius R/L0',ylabel='Barrier [meV / straight-row repeat]',
                title='Reconstruction, NOT finite-loop activation')
    for ax in axes: ax.grid(alpha=.2); ax.legend(fontsize=6)
    fig.tight_layout();plt.rcParams['svg.fonttype']='none';plt.rcParams['svg.hashsalt']='current-core-v22-comparison'
    stream=io.StringIO();fig.savefig(stream,format='svg',metadata={'Date':None})
    (args.out/'core_comparison.svg').write_text('\n'.join(x.rstrip() for x in stream.getvalue().splitlines())+'\n',encoding='utf8')
    cache=ROOT/'.cache/current_material_core_v22';cache.mkdir(exist_ok=True,parents=True)
    fig.savefig(cache/'combined_core_comparison.png',dpi=150);plt.close(fig)
    save_json(args.out/'decision.json',dict(completed=True,actual_static_stress_cases=len(stress),
        source_is_validation_target_only=True,full_current_energy_retained=True,
        force_compatibility_failure_fixed_shape_only=True,
        current_parameters_adopted=False,infinite_domain_Peierls_certified=False,
        finite_loop_activation_validated=False,physical_yield_validated=False,
        physical_seconds=False,physical_Hz=False,production_changed=False,
        files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__=='__main__': main()
