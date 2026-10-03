"""Reconstruct Gamma modes for the locked finite shape census, without refits."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def main(args):
    if args.output.exists():
        raise ValueError('fresh diagnostic output required')
    root = args.package
    read = lambda n: json.loads((root/n).read_text(encoding='utf-8'))
    selection = read('final_selection.json')
    gate = {(r['run'],r['profile']):r['eligible'] for r in selection['profiles']}
    ref_name = 'raw_runs/final_selected_phonon_recovery_results/summary.json'
    reference = np.sort(read(ref_name)['reference_gamma_optical_THz'])
    factor = np.sqrt(1.602176634e-19/(1e-20*1.66053906660e-27))/(2*np.pi*1e12)
    rows = []
    inputs = ['final_selection.json',ref_name,'diagnose_optical_census.py']
    for run in sorted({r['run'] for r in selection['profiles']}):
        name = 'raw_runs/'+run+'/profiles.json'
        inputs.append(name)
        for p in read(name):
            h = np.array(p['anchor']['hessian'])[6:,6:]
            k = np.linalg.eigvalsh(h)
            if min(k) <= 0:
                raise RuntimeError('nonpositive internal curvature in census')
            f = np.sqrt(2*k/28.0855)*factor
            rows.append(dict(run=run,profile=p['profile'],eligible=gate[run,p['profile']],
                reconstructed_gamma_optical_THz=f.tolist(),
                optical_frequency_RMSE_THz=float(np.sqrt(np.mean((f-reference)**2))),
                declared_training_objective_squared=p.get('objective_squared')))
    if len(rows) != selection['examined_profiles']:
        raise RuntimeError('locked shape census changed')
    eligible = [r for r in rows if r['eligible']]
    if len(eligible) != selection['eligible_profiles']:
        raise RuntimeError('locked eligibility changed')
    result = dict(complete=True,profiles=len(rows),eligible_profiles=len(eligible),
        reference_gamma_optical_THz=reference.tolist(),atomic_mass_amu=28.0855,
        eligible_optical_frequency_range_THz=[min(min(r['reconstructed_gamma_optical_THz']) for r in eligible),
                                            max(max(r['reconstructed_gamma_optical_THz']) for r in eligible)],
        rows=rows,inputs=[dict(path=n,sha256=hashlib.sha256((root/n).read_bytes()).hexdigest()) for n in inputs],
        new_geometry_evaluations=0,new_DFT=0,new_MD=0,selection_changed=False,
        full_bands_recomputed=False,global_shape_optimum_certified=False,family_impossibility_proved=False,
        material_approved=False,physical_clock_calibrated=False,
        interpretation='post-selection relative-basis curvature census of the finite 30 shapes only; no alternative selected by validation error')
    args.output.write_bytes((json.dumps(result,indent=2,allow_nan=False)+'\n').encode())
    print(json.dumps({k:result[k] for k in ['complete','profiles','eligible_profiles',
                                          'eligible_optical_frequency_range_THz']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
