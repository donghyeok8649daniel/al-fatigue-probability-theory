"""Independent stationary/Morse/descent audit of sampled reconstruction peaks."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from .run_local_crack_audit import save_json,polish,atomic_hessian,spectrum
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--path',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--labels',nargs='+',default=['forward_007','forward_017','reverse_005','reverse_016'])
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('use an empty output')
    start=time.perf_counter()
    reference=np.load(args.references/'initial.npz',allow_pickle=False)
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    path=np.load(args.path,allow_pickle=False)
    coordinates=RelaxedCoordinates(reference['reference_positions'],reference['fixed'],bond=reference['bond'])
    p,_=source_parameters();model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    gap=float(reference['gap_A']);sites=reference['sites'];rows=[]
    for label in args.labels:
        row={'source_label':label};saved={'source_positions':path[label]}
        try:
            stationary,correction=polish(model,coordinates,path[label],reference['fixed'],values=[gap])
            h,_,_=atomic_hessian(model,coordinates,stationary,reference['fixed'])
            spec,unstable=spectrum(h,8)
            row.update(energy_from_original_eV=float(np.sum(model.evaluate(stationary).site_energy-sites)),
                correction=correction,conditional_spectrum=spec,stationary_index=spec['negative_index'])
            saved.update(stationary=stationary,unstable_mode=unstable)
            if spec['negative_index']==1:
                descents=[]
                for sign in [-1.,1.]:
                    seed=coordinates.positions(coordinates.encode(stationary)+sign*.03*unstable,[gap])
                    r,info=relax_atoms(model,coordinates,initial=seed,values=[gap],tolerance=2e-6)
                    r,correction_min=polish(model,coordinates,r,reference['fixed'],values=[gap])
                    hmin,_,_=atomic_hessian(model,coordinates,r,reference['fixed'])
                    specmin,_=spectrum(hmin,5)
                    saved[f'descent_{int(sign):+d}']=r
                    descents.append(dict(sign=sign,energy_from_original_eV=float(np.sum(model.evaluate(r).site_energy-sites)),
                        distance_from_original_A=float(np.linalg.norm(r-reference['positions'])),
                        distance_from_reconstructed_A=float(np.linalg.norm(r-path['forward_032'])),
                        conditional_spectrum=specmin,relaxation=info,correction=correction_min))
                row['descents']=descents
        except (RuntimeError,np.linalg.LinAlgError) as error:
            row['failure']=str(error)
        rows.append(row);np.savez_compressed(args.output/(label+'.npz'),**saved)
        save_json(args.output/'running.json',rows)
        print(label,row.get('energy_from_original_eV'),row.get('stationary_index'),row.get('failure'),flush=True)
    save_json(args.output/'summary.json',dict(rows=rows,elapsed_seconds=time.perf_counter()-start,
        scope='same fixed q; each saddle and its descents audited independently; not a full minimum-energy-path certificate',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_path_sha256=hashlib.sha256(args.path.read_bytes()).hexdigest(),physical_transition_rate=None))


if __name__=='__main__':main()
