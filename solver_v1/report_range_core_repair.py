"""Re-evaluate saved relaxed cores; report material/core gates independently.

This does NOT rerun calibration/relaxation. Actual runs and their stopping
status are preserved. It evaluates analytic forces/Hessians, site partitions,
same-material elasticity and domain/ring differences at the saved states.
"""
import csv
from functools import lru_cache
import io
import json
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh

from .isolated_screw_core import IsolatedScrewCore,ScrewFarField
from .run_isolated_screw_core import OUT as CORES,material,EV_J
from .run_low_stress_cyclic_diagnostic import ROOT,write_csv
from .run_vector_registry_audit import save_json


OUT=CORES.parent


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def restore(path):
    meta=json.loads((path/'metadata.json').read_bytes())
    data=read_csv(path/'state.csv')
    state={(int(r['j']),int(r['l'])):np.array([float(r[k]) for k in ('ux','uy','uz')]) for r in data}
    return meta,state


def common_field_difference(first,second,radius):
    meta,a=first; _,b=second
    h=float(meta['h_reduced']) if 'h_reduced' in meta else np.sqrt(2/3)
    d=np.sqrt(3)/2; center=np.asarray(meta['center_over_L0'])
    changes=[]
    for (j,l),v in a.items():
        if (j,l) not in b or np.linalg.norm([d*(j+l/3)-center[0],h*l-center[1]])>=radius:
            continue
        diff=b[(j,l)]-v; diff[0]-=np.floor(diff[0]+.5)
        changes.append(diff)
    values=np.array(changes)
    return dict(common_rows=len(values),maximum_vector_change_over_L0=float(np.max(np.linalg.norm(values,axis=1))),
                rms_vector_change_over_L0=float(np.sqrt(np.mean(np.sum(values**2,axis=1)))))


def save_plot(fig,name):
    plt.rcParams['svg.fonttype']='none'; plt.rcParams['svg.hashsalt']='range-core-v12'
    stream=io.StringIO(); fig.savefig(stream,format='svg',metadata={'Date':None},bbox_inches='tight')
    (OUT/(name+'.svg')).write_text('\n'.join(line.rstrip() for line in stream.getvalue().splitlines())+'\n',encoding='utf8')
    cache=ROOT/'.cache/range_core_v12'; cache.mkdir(parents=True,exist_ok=True)
    fig.savefig(cache/(name+'.png'),dpi=125,bbox_inches='tight'); plt.close(fig)


def main():
    began=time.perf_counter(); summaries=[]; radial=[]; profiles=[]; diagnostics=[]; morse=[]; stored={}; cached={}
    models=lru_cache(maxsize=2)(material)
    for path in sorted(CORES.iterdir()):
        if not path.is_dir():
            continue
        if not (path/'summary.json').exists():
            raise RuntimeError(f'incomplete case retained: {path.name}')
        meta,state=restore(path); stored[path.name]=meta,state
        # Normalize ONLY our new output's redundant Hz spelling. Values must
        # agree; this changes neither historical outputs nor kinetic data.
        if 'physical_hz' in meta:
            assert meta['physical_hz']==meta['physical_Hz']==False
            meta.pop('physical_hz'); save_json(path/'metadata.json',meta)
        model_name='range' if path.name.startswith('range_') else 'historical'
        surface,tensor,_=models(model_name)
        far=ScrewFarField(tensor,surface.interface.bulk.geometry.b,tuple(meta['center_over_L0']))
        traction=meta['shear_traction_MPa']*1e6*meta['length_scale_m']**3/EV_J
        core=IsolatedScrewCore(surface,far,free_radius=meta['radius_over_L0'],ring=meta['ring'],
            tolerance=meta['reciprocal_tolerance'],shear_traction=traction)
        field=np.array([state[tuple(i)] for i in core.indices[core.free_ids]])
        saved=json.loads((path/'summary.json').read_bytes()); out=core.evaluate(field)
        if abs(out['energy']-saved['energy'])>2e-9:
            raise ArithmeticError('saved state energy changed unexpectedly')
        eigen=(saved.get('final_stability_probe') or {}).get('minimum_eigenvalue')
        if path.name in ('historical_R4_r4_zero','historical_R6_r4_zero',
                         'historical_R6_r8_zero','historical_R8_r6_zero','range_R4_r4_zero','historical_R4_r4_unload'):
            _,apply=core.linearize(field); n=field.size
            operator=LinearOperator((n,n),matvec=lambda v:apply(v.reshape(field.shape)).ravel(),dtype=float)
            eigen=eigsh(operator,k=1,which='SA',tol=1e-7,maxiter=3000,return_eigenvectors=False,
                         v0=np.cos(np.arange(n)*.7))[0]
        d=np.linalg.norm(core.xyz[core.energy_ids,1:]-far.center,axis=1)
        for radius in np.arange(1.,meta['radius_over_L0']+.01,.5):
            mask=d<radius; raw=float(out['site_energy'][mask].sum())
            linear=float(out['linear_reference_site_work'][mask].sum()); remainder=raw-linear
            radial.append(dict(case=path.name,radius_over_L0=radius,raw_site_energy_eV=raw,
                linear_reference_work_eV=linear,quadratic_remainder_eV=remainder,
                finite_part_at_reference_b_eV=remainder-far.log_energy_coefficient*np.log(radius/far.b),
                elastic_log_coefficient_eV=far.log_energy_coefficient))
        full=core.full_field(field); lookup={tuple(p):i for i,p in enumerate(core.indices)}
        for j in range(-int(np.ceil(meta['radius_over_L0']))-2,int(np.ceil(meta['radius_over_L0']))+3):
            if (j,0) not in lookup or (j,1) not in lookup:
                continue
            lo,hi=lookup[(j,0)],lookup[(j,1)]; slip=full[hi]-full[lo]
            profiles.append(dict(case=path.name,y_over_L0=float((core.xyz[lo,1]+core.xyz[hi,1])/2),
                slip_x_over_L0=slip[0],slip_y_over_L0=slip[1],opening_change_over_L0=slip[2],
                both_rows_free=bool(lo in core.free_ids and hi in core.free_ids)))
        linear_bond=np.broadcast_to(core.reference_bond_gradient,core.destination.shape+(3,))
        linear_gradient=core._scatter(linear_bond)
        linear_gradient[core.energy_ids]-=linear_bond.sum(axis=1)
        check=dict(case=path.name,energy_replay_error=abs(out['energy']-saved['energy']),
            reference_linear_free_gradient_max=float(np.max(abs(linear_gradient[core.free_ids]))),
            total_linear_work_eV=float(out['linear_reference_site_work'].sum()),
            affine_elasticity_error_eV_L0cubed=float(np.max(abs(core.affine_antiplane_hessian()-far.matrix))),
            minimum_tested_core_hessian_eigenvalue=eigen)
        diagnostics.append(check)
        for label,probe in [('initial',meta.get('negative_mode_probe')),('final',saved.get('final_stability_probe'))]:
            if probe is not None:
                for fd in probe['energy_curvature_checks']:
                    morse.append(dict(case=path.name,phase=label,analytic_eigenvalue=probe['minimum_eigenvalue'],
                        eigenpair_residual=probe['eigenpair_residual'],**fd))
        summaries.append(dict(case=path.name,model=model_name,radius=meta['radius_over_L0'],ring=meta['ring'],
            shear_MPa=meta['shear_traction_MPa'],free_sites=len(core.free_ids),
            energy_eV_per_row_period=saved['energy'],force_residual=saved['maximum_free_force'],
            force_converged=saved['force_converged'],winding=core.boundary_winding(field),
            minimum_core_eigenvalue=eigen,elapsed_seconds=saved['elapsed_seconds'],
            stationary_classification=('stability_unchecked' if eigen is None else
                'saddle_not_admissible_initial_core' if eigen<0 else
                'stable_fixed_boundary_candidate' if saved['force_converged'] else 'force_unconverged'),
            physical_yield_validated=False))
        if path.name in ('historical_R8_r6_zero','historical_R6_r6_stable_plus'):
            cached[path.name]=(core.xyz[core.free_ids],field)
        print(f"validated {path.name}; force={saved['maximum_free_force']:.3g}; minH={eigen}",flush=True)
    differences=[]
    pairs=[('historical_R2_r3_zero','historical_R3_r3_zero','domain'),
           ('historical_R3_r3_zero','historical_R4_r3_zero','domain'),
           ('historical_R4_r4_zero','historical_R6_r4_zero','domain'),
           ('historical_R6_r4_zero','historical_R8_r4_zero','domain'),
           ('historical_R6_r6_zero','historical_R8_r6_zero','domain'),
           ('historical_R4_r3_zero','historical_R4_r4_zero','ring'),
           ('historical_R6_r4_zero','historical_R6_r6_zero','ring'),
           ('historical_R6_r6_zero','historical_R6_r8_zero','ring'),
           ('historical_R8_r4_zero','historical_R8_r6_zero','ring'),
           ('historical_R4_r4_zero','historical_R4_r4_25MPa','load'),
           ('historical_R4_r4_zero','historical_R4_r4_50MPa','load'),
           ('historical_R4_r4_zero','historical_R4_r4_unload','static_unload_from_saddle'),
           ('historical_R4_r4_mode_minus','historical_R4_r4_unload','zero_load_control_vs_old_unload'),
           ('historical_R4_r4_mode_plus','historical_R6_r4_mode_plus','stable_domain'),
           ('historical_R6_r4_mode_plus','historical_R6_r6_stable_plus','stable_ring'),
           ('historical_R4_r4_mode_minus','historical_R4_r4_stable25MPa','stable_load'),
           ('historical_R4_r4_mode_minus','historical_R4_r4_stable50MPa','stable_load'),
           ('historical_R4_r4_mode_minus','historical_R4_r4_stable_unload','stable_static_unload')]
    for first,second,kind in pairs:
        for radius in (1.5,2.):
            differences.append(dict(first=first,second=second,test=kind,comparison_radius=radius,
                **common_field_difference(stored[first],stored[second],radius)))
    write_csv(OUT/'core_summary.csv',summaries); write_csv(OUT/'core_analytic_audit.csv',diagnostics)
    write_csv(OUT/'core_morse_audit.csv',morse)
    write_csv(OUT/'core_radial_energy.csv',radial); write_csv(OUT/'core_domain_ring_comparison.csv',differences)
    write_csv(OUT/'core_vector_disregistry.csv',profiles)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for name in ('historical_R6_r4_mode_plus','historical_R6_r6_stable_plus','historical_R8_r6_zero'):
        selected=[r for r in radial if r['case']==name]
        label=name.replace('historical_','')+(' (SADDLE)' if name.endswith('r6_zero') else ' (stable)')
        axes[0].plot([r['radius_over_L0'] for r in selected],[r['finite_part_at_reference_b_eV'] for r in selected],'o-',ms=3,label=label)
    axes[0].set(xlabel='Integration radius R/L0',ylabel='E_remainder(R) - k ln(R/b) [eV/row repeat]',
                title='Declared site partition; not a fitted core radius')
    axes[0].legend(fontsize=6); axes[0].grid(True,alpha=.25)
    for name in ('historical_R4_r4_zero','historical_R4_r4_mode_minus','historical_R4_r4_stable50MPa','historical_R4_r4_stable_unload'):
        selected=[r for r in profiles if r['case']==name]
        axes[1].plot([r['y_over_L0'] for r in selected],[r['slip_x_over_L0'] for r in selected],'.-',label=name.replace('historical_R4_r4_',''))
    axes[1].set(xlabel='y/L0',ylabel='Actual adjacent-plane x displacement / L0',
                title='Stable initial state vs centered saddle; static MPa')
    axes[1].legend(fontsize=6); axes[1].grid(True,alpha=.25)
    fig.tight_layout(); save_plot(fig,'core_matching_and_load')
    fig,axes=plt.subplots(1,2,figsize=(9,4))
    for ax,(name,(xyz,field)) in zip(axes,cached.items()):
        scatter=ax.scatter(xyz[:,1]+field[:,1],xyz[:,2]+field[:,2],c=field[:,0],s=12,cmap='coolwarm')
        ax.set(xlabel='y/L0',ylabel='z/L0',title=name.replace('historical_','')+(' (SADDLE)' if name.endswith('_zero') else ' (stable)')); ax.set_aspect('equal')
        fig.colorbar(scatter,ax=ax,label='Actual u_x/L0')
    fig.tight_layout(); save_plot(fig,'isolated_core_fields')
    save_json(OUT/'decision.json',dict(completed=True,material_accepted=False,physical_yield_validated=False,
        physical_seconds=False,physical_Hz=False,production_PDE_UI_changed=False,
        single_screw_fixed_boundary_force_balance_verified=all(r['force_converged'] for r in summaries),
        all_tested_core_stationary_points_are_stable=all(r['minimum_core_eigenvalue']>0 for r in summaries if r['minimum_core_eigenvalue'] is not None),
        stable_core_cases=[r['case'] for r in summaries if r['force_converged'] and r['minimum_core_eigenvalue'] is not None and r['minimum_core_eigenvalue']>0],
        initial_saddle_unload_matches_zero_load_control=common_field_difference(
            stored['historical_R4_r4_mode_minus'],stored['historical_R4_r4_unload'],2.),
        stable_initial_unload_residual=common_field_difference(
            stored['historical_R4_r4_mode_minus'],stored['historical_R4_r4_stable_unload'],2.),
        finite_source_core_matching_validated=False,all_character_core_energy_available=False,
        explanation='Improved material fit still misses C11 and held-out curvature; isolated straight core is not a calibrated finite source.',
        reporting_elapsed_seconds=time.perf_counter()-began))


if __name__=='__main__':
    main()
