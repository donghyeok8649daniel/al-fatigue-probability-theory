"""Correct prior screened EOS reuse outside the four-neighbor cubic branch.

Prior raw runs are preserved. Reuse their training designs, independently
reevaluate every published perfect EOS point, and report corrected selection.
This audit makes no new DFT/MD evaluations and fits no held-out frames.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from .silicon_sw_density_screening import ScreenedDiamondCell
from .run_silicon_sw_material_fit import dump


def main(args):
    if args.output.exists():raise ValueError('fresh audit output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    manifest=[]
    for path in ['solver_v1/run_silicon_sw_eos_reaudit.py','solver_v1/silicon_sw_density_screening.py',
                 'solver_v1/silicon_environment_research.py']:
        raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    reference=json.loads(args.reference.read_text());ev=np.array(reference['diamond_E_vs_V'])
    anchor=int(np.argmin(ev[:,1]));target=ev[:,1]-ev[anchor,1];rows=[];calls=0;selection_changed=False
    prior_files=[]
    for folder in args.runs:
        profile_file=folder/'profiles.json';old=json.loads(profile_file.read_text())
        prior_files.append(dict(run=folder.name,file='profiles.json',sha256=hashlib.sha256(profile_file.read_bytes()).hexdigest()))
        corrected=[]
        for row in old:
            if row['density_power']==0:
                error=row['EOS_RMSE_eV_atom'];energies=None;delta=None;count=0
            else:
                energies=[]
                for volume,_ in ev:
                    model=ScreenedDiamondCell(row['parameters'],(8*volume)**(1/3),density_power=row['density_power'],
                        density_reference=row['density_reference'],density_gamma=row.get('density_gamma'))
                    energies.append(model.evaluate(np.zeros(9),derivatives=False)/2);calls+=1
                energies=np.array(energies);delta=energies-energies[anchor]-target
                error=float(np.sqrt(np.mean(delta**2)));count=len(ev)
            result=dict(run=folder.name,profile=row['profile'],anchor_profile=row['anchor_profile'],
                density_power=row['density_power'],density_gamma=row.get('density_gamma'),
                original_EOS_RMSE_eV_atom=row['EOS_RMSE_eV_atom'],corrected_EOS_RMSE_eV_atom=error,
                original_objective_squared=row['objective_squared'],
                corrected_objective_squared=row['training_objective_squared']+(error/.02)**2,
                bulk_guard_passed=row['bulk_guard_passed'],new_EOS_energy_evaluations=count,
                EOS_energy_eV_atom=energies.tolist() if energies is not None else None,
                EOS_error_eV_atom=delta.tolist() if delta is not None else None)
            corrected.append(result);rows.append(result);dump(args.output/'profiles.json',rows)
        eligible=[r for r in corrected if r['bulk_guard_passed']]
        best=min(eligible,key=lambda r:r['corrected_objective_squared'])
        prior_selected=json.loads((folder/'selected_fit.json').read_text())
        selection_changed|=best['profile']!=prior_selected['profile']
        dump(args.output/(folder.name+'_selection.json'),best)
    dump(args.output/'prior_raw_bindings.json',prior_files)
    dump(args.output/'summary.json',dict(complete=True,profiles=len(rows),new_EOS_energy_evaluations=calls,
        correction='perfect EOS is unchanged only while four tetrahedral neighbors lie inside cutoff; compressed second neighbors invalidate reuse',
        prior_raw_results_preserved=True,no_new_training_features=True,selection_changed=bool(selection_changed),
        selection_claim='inspect corrected per-run selection files',
        material_approved=False,new_DFT=0,new_MD=0,new_LAMMPS=0,elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
