"""Replay the v13 stable states, not the minimization or material fit.

Every matching radius uses the same material and the declared site partition.
Finite-domain static recovery is not a macroscopic residual strain/hold test.
"""
import io
import json
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .core_matching import finite_part_samples
from .isolated_screw_core import IsolatedScrewCore,ScrewFarField
from .report_range_core_repair import restore,common_field_difference
from .run_isolated_screw_core import OUT as OLD,material,EV_J
from .run_low_stress_cyclic_diagnostic import ROOT,write_csv
from .run_vector_registry_audit import save_json


OUT=ROOT/'results/fcc111_active_interface/stable_core_v13'
CORES=OUT/'isolated_core'


def main():
    began=time.perf_counter(); surface,tensor,_=material('historical')
    stored={}; summaries=[]; radial=[]; profiles=[]; checks=[]
    paths=[OLD/name for name in ('historical_R6_r6_stable_plus',
        'historical_R4_r4_mode_minus','historical_R4_r4_stable50MPa',
        'historical_R4_r4_stable_unload')]
    paths+=sorted(p for p in CORES.iterdir() if p.is_dir())
    for path in paths:
        if not (path/'summary.json').exists():
            raise RuntimeError(f'preserve incomplete run: {path.name}')
        meta,state=restore(path); stored[path.name]=meta,state
        saved=json.loads((path/'summary.json').read_bytes())
        far=ScrewFarField(tensor,surface.interface.bulk.geometry.b,tuple(meta['center_over_L0']))
        traction=meta['shear_traction_MPa']*1e6*meta['length_scale_m']**3/EV_J
        core=IsolatedScrewCore(surface,far,free_radius=meta['radius_over_L0'],ring=meta['ring'],
            tolerance=meta['reciprocal_tolerance'],shear_traction=traction)
        field=np.array([state[tuple(i)] for i in core.indices[core.free_ids]])
        value=core.evaluate(field); error=abs(value['energy']-saved['energy'])
        if error>2e-9:
            raise ArithmeticError('stored energy did not replay under its actual model')
        probe=saved['final_stability_probe']
        summaries.append(dict(case=path.name,source=path.relative_to(ROOT).as_posix(),
            parameter_sha256=meta['parameter_sha256'],radius_over_L0=meta['radius_over_L0'],ring=meta['ring'],
            static_resolved_shear_MPa=meta['shear_traction_MPa'],free_rows=len(core.free_ids),
            energy_eV_per_row_repeat=value['energy'],maximum_force_eV_L0=saved['maximum_free_force'],
            minimum_H_eV_L0sq=probe['minimum_eigenvalue'],eigenpair_residual=probe['eigenpair_residual'],
            force_converged=saved['force_converged'],stable_fixed_boundary=probe['stable_on_tested_fixed_boundary'],
            method=saved.get('minimizer','lbfgs'),iterations=saved['iterations'],evaluations=saved['evaluations'],
            elapsed_wall_seconds=saved['elapsed_seconds'],winding=core.boundary_winding(field),
            reciprocal_modes=saved['reciprocal_modes'],last_mode_envelope=saved['last_mode_envelope'],
            energy_replay_error_eV=error,material_accepted=False,physical_yield_validated=False))
        for check in probe['energy_curvature_checks']:
            checks.append(dict(case=path.name,analytic_minimum_H_eV_L0sq=probe['minimum_eigenvalue'],**check))
        if meta['shear_traction_MPa']==0:
            distances=np.linalg.norm(core.xyz[core.energy_ids,1:]-far.center,axis=1)
            remainder=value['quadratic_remainder_site_energy']
            radii=np.arange(2.,meta['radius_over_L0']+.01,.5)
            for sample in finite_part_samples(distances,remainder,radii,
                    log_coefficient=far.log_energy_coefficient,b=far.b,free_radius=meta['radius_over_L0']):
                radial.append(dict(case=path.name,partition='smooth',**sample))
            for R in radii:
                E=float(remainder[distances<R].sum())
                radial.append(dict(case=path.name,partition='hard',radius_over_L0=R,inner_fraction=1.,
                    weighted_remainder_eV=E,continuum_subtraction_eV=far.log_energy_coefficient*np.log(R/far.b),
                    finite_part_at_reference_b_eV=E-far.log_energy_coefficient*np.log(R/far.b),
                    window_constant=0.,window_series_tail=0.,window_series_terms=0,
                    log_coefficient_eV=far.log_energy_coefficient,core_radius_fitted=False,convergence_certified=False))
        full=core.full_field(field); lookup={tuple(p):i for i,p in enumerate(core.indices)}
        for j in range(-int(np.ceil(meta['radius_over_L0']))-2,int(np.ceil(meta['radius_over_L0']))+3):
            if (j,0) not in lookup or (j,1) not in lookup: continue
            lo,hi=lookup[j,0],lookup[j,1]; delta=full[hi]-full[lo]
            profiles.append(dict(case=path.name,y_over_L0=float((core.xyz[lo,1]+core.xyz[hi,1])/2),
                delta_ux_over_L0=delta[0],delta_uy_over_L0=delta[1],delta_uz_over_L0=delta[2],
                both_rows_free=bool(lo in core.free_ids and hi in core.free_ids)))
        print(f'replayed {path.name}: {error:.3g} eV error',flush=True)
    pairs=[('historical_R6_r6_stable_plus','historical_R6_r6_newton_check','algorithm'),
        ('historical_R6_r6_stable_plus','historical_R8_r6_stable_plus','domain'),
        ('historical_R8_r6_stable_plus','historical_R10_r6_stable_plus','domain'),
        ('historical_R6_r8_stable_plus','historical_R8_r8_stable_plus','domain'),
        ('historical_R8_r8_stable_plus','historical_R10_r8_stable_plus','domain'),
        ('historical_R6_r6_stable_plus','historical_R6_r8_stable_plus','ring'),
        ('historical_R8_r6_stable_plus','historical_R8_r8_stable_plus','ring'),
        ('historical_R10_r6_stable_plus','historical_R10_r8_stable_plus','ring'),
        ('historical_R8_r8_stable_plus','historical_R8_r8_50MPa','static_load'),
        ('historical_R8_r8_stable_plus','historical_R8_r8_unload','static_unload'),
        ('historical_R4_r4_mode_minus','historical_R4_r4_stable50MPa','legacy_static_load'),
        ('historical_R4_r4_mode_minus','historical_R4_r4_stable_unload','legacy_static_unload')]
    differences=[]; match_changes=[]
    by_sample={(r['case'],r['partition'],r['radius_over_L0'],r['inner_fraction']):r for r in radial}
    for first,second,kind in pairs:
        if first not in stored or second not in stored:
            raise RuntimeError(f'planned comparison not completed: {first} / {second}')
        if stored[first][0]['parameter_sha256']!=stored[second][0]['parameter_sha256']:
            raise ArithmeticError('do not compare domain changes across materials')
        for radius in (1.5,2.):
            differences.append(dict(first=first,second=second,test=kind,comparison_radius_over_L0=radius,
                **common_field_difference(stored[first],stored[second],radius)))
        if kind in ('domain','ring','algorithm'):
            for (case,partition,radius,alpha),sample in by_sample.items():
                other=by_sample.get((second,partition,radius,alpha))
                if case!=first or other is None: continue
                match_changes.append(dict(first=first,second=second,test=kind,partition=partition,
                    radius_over_L0=radius,inner_fraction=alpha,
                    first_finite_part_eV=sample['finite_part_at_reference_b_eV'],
                    second_finite_part_eV=other['finite_part_at_reference_b_eV'],
                    finite_part_change_eV=other['finite_part_at_reference_b_eV']-sample['finite_part_at_reference_b_eV']))
    write_csv(OUT/'core_summary.csv',summaries); write_csv(OUT/'core_curvature_checks.csv',checks)
    write_csv(OUT/'core_field_comparisons.csv',differences); write_csv(OUT/'core_radial_matching.csv',radial)
    write_csv(OUT/'core_matching_changes.csv',match_changes); write_csv(OUT/'core_disregistry.csv',profiles)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for name in ('historical_R6_r8_stable_plus','historical_R8_r8_stable_plus','historical_R10_r8_stable_plus'):
        for alpha,style in ((.5,'-'),(.75,'--')):
            selected=[r for r in radial if r['case']==name and r['inner_fraction']==alpha]
            axes[0].plot([r['radius_over_L0'] for r in selected],
                [r['finite_part_at_reference_b_eV'] for r in selected],style,
                label=name.split('_')[1]+f', window {alpha}')
    axes[0].set(xlabel='Integration radius R/L0',ylabel='Finite part [eV/row repeat]',
        title='Same-potential logarithm; no fitted core radius')
    for name in ('historical_R8_r8_stable_plus','historical_R8_r8_50MPa','historical_R8_r8_unload'):
        selected=[r for r in profiles if r['case']==name and r['both_rows_free']]
        axes[1].plot([r['y_over_L0'] for r in selected],[r['delta_ux_over_L0'] for r in selected],
            '.-',label=name.removeprefix('historical_R8_r8_'))
    axes[1].set(xlabel='y/L0',ylabel='Adjacent-plane displacement / L0',
        title='0 -> 50 -> 0 MPa static shear; not kinetic hold')
    for ax in axes: ax.legend(fontsize=6); ax.grid(True,alpha=.25)
    fig.tight_layout(); plt.rcParams['svg.fonttype']='none'; plt.rcParams['svg.hashsalt']='stable-core-v13'
    stream=io.StringIO(); fig.savefig(stream,format='svg',metadata={'Date':None},bbox_inches='tight')
    (OUT/'stable_core_matching.svg').write_text('\n'.join(line.rstrip() for line in stream.getvalue().splitlines())+'\n',encoding='utf8')
    fig.savefig(ROOT/'.cache/stable_core_v13_matching.png',dpi=135,bbox_inches='tight'); plt.close(fig)
    save_json(OUT/'core_decision.json',dict(completed=True,core_runs_replayed=len(summaries),
        all_fixed_boundary_stable=all(r['stable_fixed_boundary'] for r in summaries),
        maximum_replay_error=max(r['energy_replay_error_eV'] for r in summaries),
        static_unload_difference=common_field_difference(stored['historical_R8_r8_stable_plus'],
            stored['historical_R8_r8_unload'],2.),finite_part_all_domain_limit_certified=False,
        boundary_is_not_measured_pinning=True,macroscopic_residual_plasticity_validated=False,
        material_accepted=False,physical_yield_validated=False,physical_seconds=False,physical_Hz=False,
        production_PDE_UI_changed=False,elapsed_replay_seconds=time.perf_counter()-began))


if __name__=='__main__':
    main()
