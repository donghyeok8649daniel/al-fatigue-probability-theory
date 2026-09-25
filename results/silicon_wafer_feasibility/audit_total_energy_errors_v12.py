"""Expose whole-configuration energy errors hidden by per-atom normalization.

All pairs within exact composition use existing predictions. They are neither
independent samples nor transition paths, so these are not barrier uncertainties.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np


def main(args):
    if args.output.exists():raise ValueError('fresh energy audit output required')
    records=[];sources=[]
    for family in ('radial','angular'):
        source=args.results/('oxide_'+family+'_delta')/'held_composition_frames.csv'
        rows=list(csv.DictReader(source.open(encoding='utf-8')))
        for composition in sorted({row['composition'] for row in rows}):
            group=[row for row in rows if row['composition']==composition]
            counts={int(row['atoms']) for row in group}
            if len(counts)!=1:raise ValueError('composition count mismatch')
            atoms=counts.pop();number=len(group);pairs=number*(number-1)//2
            before=np.array([float(row['baseline_energy_eV'])-float(row['reference_energy_eV']) for row in group])
            after=np.array([float(row['candidate_energy_eV'])-float(row['reference_energy_eV']) for row in group])
            # Sum_{i<j}(e_i-e_j)^2 = n*sum_i(e_i-mean(e))^2.
            # Energy constants cancel; no offset is fitted or deployed.
            b=float(np.sqrt(number*np.sum((before-before.mean())**2)/pairs))
            a=float(np.sqrt(number*np.sum((after-after.mean())**2)/pairs))
            reference=json.loads((source.parent/'cross_validation.json').read_text(encoding='utf-8'))
            stored=next(row for row in reference['held_composition_relative_energy'] if row['composition']==composition)
            replay=max(abs(b/atoms-stored['baseline_pair_energy_RMSE_eV_atom']),abs(a/atoms-stored['candidate_pair_energy_RMSE_eV_atom']))
            if replay>1e-10:raise ValueError('whole-configuration and per-atom energy analyses disagree')
            records.append(dict(family=family,composition=composition,frames=number,atoms_per_frame=atoms,pairs=pairs,
                baseline_pair_difference_RMSE_eV=b,candidate_pair_difference_RMSE_eV=a,
                baseline_pair_difference_RMSE_eV_atom=b/atoms,candidate_pair_difference_RMSE_eV_atom=a/atoms,
                baseline_worst_pair_difference_error_eV=float(np.ptp(before)),candidate_worst_pair_difference_error_eV=float(np.ptp(after)),
                per_atom_replay_error_eV_atom=replay))
        sources.append(dict(path=source.relative_to(args.results).as_posix(),sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    args.output.mkdir(parents=True)
    with (args.output/'total_energy_errors.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    report=dict(complete=True,records=records,sources=sources,
        energy_reference='pair differences only within exact chemical composition; no cross-stoichiometry formation energies',
        caveat='many source pairs are correlated and do not lie on validated initiation paths; whole-configuration differences are not barrier-error estimates',
        new_fits=0,new_MACE=0,new_DFT=0,new_MD=0,material_approved=False,physical_clock=None)
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
