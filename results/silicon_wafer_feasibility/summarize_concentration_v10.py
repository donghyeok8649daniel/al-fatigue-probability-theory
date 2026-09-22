"""Rebuild concentration tables and scientific figures from saved raw results."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_concentration_research import constrained_energy_derivative
from results.silicon_wafer_feasibility.run_concentration_cleavage_v10 import dump,csv_write


def main(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summary=json.loads((args.source/'summary.json').read_text())
    protocol=json.loads((args.source/'protocol.json').read_text())
    if not summary['complete'] or len(summary['cases'])!=len(protocol['plans']):raise ValueError('source incomplete')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    cases=summary['cases'];pure=next(r for r in cases if r['species']=='Si')
    host=np.load(args.source/'Si_00_pure'/'opening_0.000.npz')
    depth=np.minimum(host['positions'][:,2],host['cell'][2,2]-host['positions'][:,2])
    adjacent=depth<=depth.min()+1e-6
    raw_errors=[];rows=[];curves={}
    for result in cases:
        tag=result['case'];folder=args.source/tag
        raw_reference=np.load(folder/'opening_0.000.npz')
        actual_count=int(np.count_nonzero(raw_reference['numbers']!=14))
        if actual_count!=result['dopant_count'] or len(raw_reference['numbers'])!=216:raise AssertionError('chemical count')
        density=actual_count*1e24/np.linalg.det(raw_reference['cell'])
        if not np.isclose(density,result['chemical_concentration_cm3'],rtol=1e-14,atol=0):raise AssertionError('concentration')
        curve=list(csv.DictReader((folder/'curve.csv').open()))
        for row in curve:
            q=float(row['opening_A']);raw=np.load(folder/f'opening_{q:.3f}.npz')
            derivative=constrained_energy_derivative(raw['positions'],raw['forces'],raw['cell'],raw['stress'],
                position_tangent=raw['position_tangent'],cell_tangent=raw['cell_tangent'])
            work=(float(raw['energy'])-float(raw_reference['energy']))/protocol['area_A2']*16.02176634
            raw_errors.append(dict(case=tag,opening_A=q,
                derivative_error_eV_A=derivative-float(row['conjugate_force_eV_A']),
                work_error_J_m2=work-float(row['separation_work_J_m2'])))
        curves[tag]=curve
        rows.append(dict(case=tag,species=result['species'],dopant_count=result['dopant_count'],
            arrangement=result['arrangement'],chemical_concentration_cm3=density,
            atomic_site_fraction=result['dopant_count']/216,
            initial_nearest_plane_dopants=int(adjacent[np.asarray(result['indices'],int)].sum()),
            fixed_material_volume_A3=result['volume_A3'],
            rigid_work_J_m2=result['rigid_separation_work_J_m2'],
            work_change_from_pure_percent=100*(result['rigid_separation_work_J_m2']/pure['rigid_separation_work_J_m2']-1),
            sampled_peak_traction_GPa=result['sampled_peak_traction_GPa'],
            sampled_peak_change_from_pure_percent=100*(result['sampled_peak_traction_GPa']/pure['sampled_peak_traction_GPa']-1),
            bulk_sigma_zz_GPa=result['stress_GPa'][2][2],force_max_eV_A=result['force_max_eV_A'],
            plateau_8_to_10_eV=result['plateau_8_to_10_eV']))
    csv_write(args.output/'concentration_effects.csv',rows)
    csv_write(args.output/'raw_recalculation.csv',raw_errors)
    spreads=[]
    for species in ['B','P']:
        for n in sorted({r['dopant_count'] for r in rows if r['species']==species}):
            group=[r for r in rows if r['species']==species and r['dopant_count']==n]
            values=[r['rigid_work_J_m2'] for r in group]
            spreads.append(dict(species=species,dopant_count=n,chemical_concentration_cm3=group[0]['chemical_concentration_cm3'],
                selected_placement_count=len(group),min_work_J_m2=min(values),max_work_J_m2=max(values),
                work_range_percent_of_pure=100*(max(values)-min(values))/pure['rigid_separation_work_J_m2']))
    csv_write(args.output/'placement_spread.csv',spreads)
    duplicate=[]
    for species in ['B','P']:
        a=next(r for r in cases if r['species']==species and r['dopant_count']==1 and r['arrangement']=='interface_cluster')
        b=next(r for r in cases if r['species']==species and r['dopant_count']==1 and r['arrangement']=='interface_spread')
        duplicate.append(dict(species=species,same_indices=a['indices']==b['indices'],
            work_difference_J_m2=a['rigid_separation_work_J_m2']-b['rigid_separation_work_J_m2']))
    error=max(max(abs(r['derivative_error_eV_A']),abs(r['work_error_J_m2'])) for r in raw_errors)
    if error>1e-10 or any(not r['same_indices'] or abs(r['work_difference_J_m2'])>1e-9 for r in duplicate):
        raise AssertionError('saved raw calculation or identical placement replay mismatch')
    plt.rcParams.update({'font.size':10,'figure.dpi':140,'axes.spines.top':False,'axes.spines.right':False})
    styles={'interface_cluster':('-','o','near-plane cluster'),
            'interface_spread':('--','s','near-plane spread'),
            'bulk_spread':(':','^','volume spread')}
    fig,axes=plt.subplots(2,2,figsize=(11,8),sharex=True)
    for col,species in enumerate(['B','P']):
        for arrangement,(style,marker,label) in styles.items():
            selected=sorted([r for r in rows if r['species']==species and r['arrangement']==arrangement],key=lambda r:r['dopant_count'])
            x=[0]+[r['chemical_concentration_cm3']/1e21 for r in selected]
            for ax,key,baseline in [(axes[0,col],'rigid_work_J_m2',pure['rigid_separation_work_J_m2']),
                                    (axes[1,col],'sampled_peak_traction_GPa',pure['sampled_peak_traction_GPa'])]:
                ax.plot(x,[baseline]+[r[key] for r in selected],linestyle=style,marker=marker,label=label)
                ax.grid(alpha=.2)
        axes[0,col].set_title(species+' substitution / fixed 216 sites')
        axes[1,col].set_xlabel('Chemical concentration (10$^{21}$ cm$^{-3}$)')
        axes[0,col].legend(fontsize=8)
    axes[0,0].set_ylabel('Rigid separation work (J/m$^2$)')
    axes[1,0].set_ylabel('Sampled rigid traction maximum (GPa)')
    fig.suptitle('Neutral MACE concentration and placement diagnostic\nStatic constrained paths; no wafer strength or carrier activation calibration',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(args.output/'concentration_and_placement.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,species in zip(axes,['B','P']):
        for n in [0,1,2,4,8]:
            tag='Si_00_pure' if n==0 else f'{species}_{n:02d}_interface_cluster'
            curve=curves[tag]
            ax.plot([float(r['opening_A']) for r in curve],[float(r['separation_work_J_m2']) for r in curve],marker='.',label=f'{n} dopants')
        ax.set_title(species+' / near-plane cluster');ax.set_xlim(0,6)
        ax.set_xlabel('Rigid interface opening (Angstrom)');ax.set_ylabel('Work / projected area (J/m$^2$)')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Actual potential-energy evaluations; internally relaxed bulk, frozen opening paths')
    fig.tight_layout();fig.savefig(args.output/'rigid_energy_paths.png');plt.close(fig)
    loaded=[];loaded_raw_errors=[]
    for path in args.loaded:
        data=json.loads((path/'summary.json').read_text())
        if not data['complete']:raise ValueError('loaded path incomplete')
        loaded.extend(data['cases'])
        for result in data['cases']:
            folder=path/result['case']
            records=list(csv.DictReader((folder/'evaluations.csv').open()))
            fresh=[r for r in records if r['new_evaluation']=='True']
            files=sorted(folder.glob('evaluation_*.npz'))
            if len(files)!=len(fresh) or len(files)!=result['new_standard_calculator_calls']:
                raise AssertionError('loaded raw evaluation count')
            for file,row in zip(files,fresh):
                raw=np.load(file);dc=np.zeros((3,3));dc[2,2]=1.
                derivative=constrained_energy_derivative(raw['positions'],raw['forces'],raw['cell'],raw['stress'],
                    position_tangent=np.zeros_like(raw['positions']),cell_tangent=dc)
                errors=dict(case=result['case'],opening_A=float(raw['opening_A']),
                    opening_error_A=float(raw['opening_A'])-float(row['opening_A']),
                    energy_error_eV=float(raw['energy'])-float(row['energy_eV']),
                    derivative_error_eV_A=derivative-float(row['energy_derivative_eV_A']))
                if any(abs(errors[k])>1e-10 for k in ['opening_error_A','energy_error_eV','derivative_error_eV_A']):
                    raise AssertionError('loaded raw reconstruction mismatch')
                loaded_raw_errors.append(errors)
    if loaded:
        csv_write(args.output/'loaded_raw_recalculation.csv',loaded_raw_errors)
        keys=[r['case'] for r in loaded]
        if len(set(keys))!=len(keys):raise ValueError('duplicate loaded cases')
        if set(keys)!={r['case'] for r in cases}:raise ValueError('loaded concentration sweep incomplete')
        out=[]
        for r in loaded:
            from scipy.constants import electron_volt
            well=min([v for v in r['stationary_points'] if v['kind']=='minimum'],key=lambda v:abs(v['opening_A']-r['minimum_opening_A']))
            maximum=min([v for v in r['stationary_points'] if v['kind']=='maximum'],key=lambda v:abs(v['opening_A']-r['maximum_opening_A']))
            # Independently convert pressure (Pa)*area(m2)*distance(m) to eV.
            external_work=r['local_traction_MPa']*1e6*r['area_A2']*1e-20*(maximum['opening_A']-well['opening_A'])*1e-10/electron_volt
            barrier=maximum['energy_eV']-well['energy_eV']-external_work
            if abs(barrier-r['constrained_path_barrier_eV'])>1e-10:raise AssertionError('loaded work/energy units')
            out.append({k:r[k] for k in ['case','species','dopant_count','arrangement','chemical_concentration_cm3',
                'local_traction_MPa','minimum_opening_A','maximum_opening_A','constrained_path_barrier_eV',
                'barrier_per_area_J_m2','root_force_max_error_eV_A','root_energy_difference_max_error_eV_A']})
        csv_write(args.output/'loaded_path_barriers.csv',out)
        fig,axes=plt.subplots(1,2,figsize=(11,4.6),sharey=True)
        reference=next(r for r in loaded if r['species']=='Si')['constrained_path_barrier_eV']
        for ax,species in zip(axes,['B','P']):
            for arrangement,(style,marker,label) in styles.items():
                selected=sorted([r for r in loaded if r['species']==species and r['arrangement']==arrangement],key=lambda r:r['dopant_count'])
                ax.plot([0]+[r['chemical_concentration_cm3']/1e21 for r in selected],
                        [0]+[100*(r['constrained_path_barrier_eV']/reference-1) for r in selected],
                        linestyle=style,marker=marker,label=label)
            ax.axhline(0,color='black',linewidth=.7);ax.set_title(species+' substitution')
            ax.set_xlabel('Chemical concentration (10$^{21}$ cm$^{-3}$)');ax.grid(alpha=.2);ax.legend(fontsize=8)
        axes[0].set_ylabel('Constrained path barrier change from pure Si (%)')
        fig.suptitle('Same LOCAL normal traction: 100 MPa\nRigid many-atom paths; not finite-T activation barriers or wafer probabilities',fontsize=12)
        fig.tight_layout(rect=(0,0,1,.93));fig.savefig(args.output/'loaded_concentration_barriers.png');plt.close(fig)
    audit=dict(complete=True,cases=len(cases),raw_states=len(raw_errors),raw_recalculation_max_abs_error=error,
        one_atom_duplicate_placement_controls=duplicate,loaded_cases=len(loaded),loaded_raw_states=len(loaded_raw_errors),
        scope='selected neutral configurations; not a random-alloy average or actual concentration-to-failure law')
    dump(args.output/'summary.json',audit)
    print(json.dumps(audit),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--loaded',nargs='*',type=Path,default=[])
    main(p.parse_args())
