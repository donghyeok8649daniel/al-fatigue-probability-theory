"""Join published B tables with completed neutral ML relaxations and Hessians."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main(root):
    run=json.loads((root/'mace_boron/run_summary.json').read_text(encoding='utf-8'))
    source=list(csv.DictReader((root/'sources/boron_local_modes_2025.csv').open(encoding='utf-8')))
    lookup={r['defect']:[float(v) for v in r['top_local_frequencies_cm_minus1'].split(';')] for r in source}
    keys={'b_silicon_64_interstitial':'B_interstitial_neutral','2boron_int_cluster':'2B_interstitial_neutral',
          **{f'3boron_int_cluster_confi_{i}':f'3B_interstitial_conf{i}' for i in [1,2,3]}}
    ordered=['b_silicon_64_interstitial','2boron_int_cluster',*[f'3boron_int_cluster_confi_{i}' for i in [1,2,3]]]
    frequencies=[];states=[]
    for name in ordered:
        h=json.loads((root/'mace_boron_hessian'/(name+'_hessian.json')).read_text(encoding='utf-8'))
        predicted=h['highest_3NB_Gamma_frequencies_cm_minus1'];reference=lookup[keys[name]]
        if len(predicted)!=len(reference):raise ValueError('not enough positive modes; inspect instability first')
        errors=np.array(predicted)-reference
        states.append(dict(structure=name,force_max_eV_A=h['max_force_eV_A'],
            fixed_cell_min_curvature_eV_A2=h['Cartesian_projected_min_eigenvalue_eV_A2'],
            Hessian_status=h['Hessian_status'],max_initial_to_final_displacement_A=h['max_initial_to_final_minimum_image_displacement_A'],
            sorted_top_frequency_RMSE_cm_minus1=float(np.sqrt(np.mean(errors**2))),
            sorted_top_frequency_max_error_cm_minus1=float(np.max(abs(errors))),
            published_max_frequency_cm_minus1=max(reference),MACE_max_frequency_cm_minus1=max(predicted)))
        for rank,(p,t) in enumerate(zip(predicted,reference),1):
            frequencies.append(dict(structure=name,rank_within_top_3NB=rank,
                published_LDA_frequency_cm_minus1=t,MACE_frequency_cm_minus1=p,difference_cm_minus1=p-t,
                correspondence='sorted spectrum only; eigenvectors and final DFT geometries unavailable'))
    write_csv(root/'boron_state_comparison.csv',states)
    write_csv(root/'boron_frequency_comparison.csv',frequencies)
    comparisons=[]
    for i,row in enumerate(run['comparisons'],1):
        state=next(s for s in states if s['structure']==f'3boron_int_cluster_confi_{i}')
        comparisons.append(dict(**row,MACE_Hessian_status=state['Hessian_status'],
            comparison_scope='source-seeded stationary states; not all are local minima'))
    write_csv(root/'boron_relative_energy_comparison.csv',comparisons)
    checks=json.loads((root/'mace_boron/derivative_checks.json').read_text(encoding='utf-8'))
    refinement=[c for c in checks if c.get('step')==1e-4]
    summary=dict(neutral_MACE_input_structures=len(run['cases']),all_force_converged=all(c['converged'] for c in run['cases']),
        source_charged_acceptor_targets_compared_with_neutral_model=False,
        neutral_interstitial_Hessians=len(states),minimum_fixed_cell_curvature=min(r['fixed_cell_min_curvature_eV_A2'] for r in states),
        initial_relaxation_calculator_evaluations=run['total_calculator_evaluations'],
        initial_force_FD_max_abs_error_at_1e_minus4_A=max(c['absolute_error'] for c in refinement if c['kind']=='force'),
        initial_stress_FD_max_abs_error_at_1e_minus4_strain=max(c['absolute_error'] for c in refinement if c['kind']=='stress'),
        relative_energy_comparisons=comparisons,state_comparisons=states,
        unstable_original_states=[r['structure'] for r in states if r['fixed_cell_min_curvature_eV_A2']<0],
        negative_mode_following_results='mace_boron_follow; separate final Hessians and endpoint validation',
        material_status='not adopted: quantitative B energy/spectrum discrepancies and different source theory/final geometry',
        fixed_cell_Gamma_only=True,complete_phonon_dispersion=False,finite_T_free_energy=False,
        new_DFT_runs=0,new_MD_runs=0,kinetic_calibrated=False)
    (root/'boron_validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    x=np.arange(3)
    axes[0].bar(x-.17,[r['published_LDA_relative_formation_energy_eV'] for r in run['comparisons']],.34,label='Published LDA',color='#8a8f99')
    bars=axes[0].bar(x+.17,[r['MACE_relative_energy_eV'] for r in run['comparisons']],.34,label='Neutral MACE stationary states',color='#3274a1')
    bars[1].set_hatch('///')
    axes[0].annotate('MACE saddle',xy=(1.17,run['comparisons'][1]['MACE_relative_energy_eV']),
                     xytext=(.6,2.4),arrowprops=dict(arrowstyle='->'),fontsize=9,color='#a33c32')
    axes[0].set(xticks=x,xticklabels=['3B config 1','3B config 2','3B config 3'],ylabel='Relative energy (eV)',title='Source-seeded states; config 1 reference')
    axes[0].legend()
    x=np.arange(5)
    axes[1].plot(x,[r['published_max_frequency_cm_minus1'] for r in states],'s--',color='#8a8f99',label='Published LDA (reported flat modes)')
    axes[1].plot(x,[r['MACE_max_frequency_cm_minus1'] for r in states],'o-',color='#3274a1',label='Neutral MACE at Gamma')
    axes[1].set(xticks=x,xticklabels=['1B','2B','3B-1','3B-2\n(saddle)','3B-3'],ylabel='Highest frequency (cm$^{-1}$)',title='Stationary interstitial spectra; sorted only')
    axes[1].legend()
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.suptitle('Source-seeded B/Si relaxations: quantitative transfer remains unvalidated')
    fig.savefig(root/'boron_validation.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    main(p.parse_args().root)
