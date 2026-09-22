"""Recompute metrics from raw arrays and cross-check the force-only backend."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from solver_v1.silicon_atomistic_reference import periodic_basis_audit


def main(args):
    import torch
    from mace.calculators import MACECalculator
    from ase.build import bulk
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model mismatch')
    args.output.mkdir(parents=True,exist_ok=True)
    new_file=np.load(args.results/'predictions.npz')
    old_file=np.load(args.baseline/'dft_predictions.npz')
    new={k:new_file[k] for k in new_file.files};old={k:old_file[k] for k in old_file.files}
    for k in ['offsets','reference_force','reference_energy']:
        if not np.array_equal(new[k],old[k]):raise AssertionError('v5/v9 source arrays differ: '+k)
    rows=list(csv.DictReader((args.results/'frame_metrics.csv').open(encoding='utf-8')))
    report=json.loads((args.results/'summary.json').read_text(encoding='utf-8'))
    baseline=json.loads((args.baseline/'dft_material_audit.json').read_text(encoding='utf-8'))
    offsets=new['offsets'];groups=defaultdict(list);residuals=[]
    if len(rows)!=2475 or len(new['predicted_energy'])!=2475:raise AssertionError('incomplete frame count')
    for i,row in enumerate(rows):
        left,right=offsets[i:i+2]
        difference=new['predicted_force'][left:right]-new['reference_force'][left:right]
        rmse=float(np.sqrt(np.mean(difference*difference)))
        residuals.append(abs(rmse-float(row['component_rmse_eV_A'])))
        if int(row['frame'])!=i:raise AssertionError('frame ordering')
        groups[(row['config_type'],row['declared_xc'])].append(i)
    if max(residuals)>1e-12:raise AssertionError('raw force metric mismatch')
    comparisons=[]
    for (kind,xc),indices in sorted(groups.items()):
        indices_atoms=np.concatenate([np.arange(offsets[i],offsets[i+1]) for i in indices])
        reference=new['reference_force'][indices_atoms]
        d=dict(config_type=kind,declared_xc=xc,frames=len(indices),atoms=len(indices_atoms),
               reference_component_rms_eV_A=float(np.sqrt(np.mean(reference**2))))
        for name,prediction in [('MACE',new['predicted_force']),('SW',old['original_sw_1985_force']),('Tersoff',old['tersoff_1989_force'])]:
            d[name+'_component_rmse_eV_A']=float(np.sqrt(np.mean((prediction[indices_atoms]-reference)**2)))
        saved=next(s for s in report['group_results'] if s['config_type']==kind and s['declared_xc']==xc)
        if abs(d['MACE_component_rmse_eV_A']-saved['component_rmse_eV_A'])>1e-12:
            raise AssertionError('group metrics mismatch')
        comparisons.append(d)
    with (args.output/'same_geometry_comparison.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
    frames,_=read_frames(args.archive)
    crack_indices=groups[('crack_111_1-10','PW91')]
    left_handed=[i for i,a in enumerate(frames) if np.linalg.det(a.cell.array)<0]
    controls=sorted(set([0,548,crack_indices[0],crack_indices[-1],left_handed[0],left_handed[-1],
                        groups[('surface_111_pandey','PBE')][0]]))
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    checks=[]
    for i in controls:
        a=frames[i].copy();a.set_cell(periodic_basis_audit(a.cell.array,a.pbc)[0],scale_atoms=False)
        ep,fp=force_only(calc,a);left,right=offsets[i:i+2]
        de=abs(ep-new['predicted_energy'][i]);df=float(np.max(abs(fp-new['predicted_force'][left:right])))
        if de>1e-9 or df>1e-9:raise AssertionError('force-only/ASE mismatch')
        checks.append(dict(frame=i,atoms=len(a),energy_difference_eV=de,force_difference_eV_A=df))
    # A deterministic distorted cell; rotate both positions and cell together.
    a=bulk('Si','diamond',a=5.47,cubic=True)
    a.positions+=.03*np.sin(np.arange(24).reshape(8,3)+.17)
    ep,fp=force_only(calc,a)
    theta=.437;rotation=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    b=a.copy();b.positions=a.positions@rotation.T;b.set_cell(a.cell.array@rotation.T,scale_atoms=False)
    er,fr=force_only(calc,b)
    repeated=a.repeat((2,1,1));et,ft=force_only(calc,repeated)
    invariance=dict(rotation_energy_error_eV=abs(er-ep),
        rotation_force_error_eV_A=float(np.max(abs(fr-fp@rotation.T))),
        replication_energy_error_eV=abs(et-2*ep),
        replication_force_error_eV_A=float(np.max(abs(ft-np.tile(fp,(2,1))))))
    if max(invariance.values())>1e-9:raise AssertionError('symmetry/extensivity failure')
    final=dict(raw_frames_checked=len(rows),source_arrays_exactly_equal_to_v5=True,
        max_frame_metric_recompute_error=max(residuals),group_count=len(comparisons),
        force_only_vs_standard_ASE=checks,symmetry_and_replication=invariance,
        model_SHA256=MODEL_SHA256,material_calibrated=False,kinetic_calibrated=False,
        scope='independent aggregation and same-model execution/derivative controls; not independent DFT')
    (args.output/'raw_validation.json').write_text(json.dumps(final,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(final,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected=[('dia','PW91'),('crack_111_1-10','PW91'),('crack_110_1-10','PW91'),('surface_111','PW91'),('surface_111_pandey','PBE'),('decohesion','UNSPECIFIED')]
    labels=['Diamond\nPW91','(111) crack\nPW91','(110) crack\nPW91','(111) surface\nPW91','Pandey surface\nPBE','Decohesion\nXC unspecified']
    fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
    x=np.arange(len(selected))
    for j,(name,color) in enumerate([('SW','#8a8f99'),('Tersoff','#dc8f3a'),('MACE','#3274a1')]):
        values=[next(r for r in comparisons if (r['config_type'],r['declared_xc'])==key)[name+'_component_rmse_eV_A'] for key in selected]
        ax.bar(x+(j-1)*.25,values,width=.24,label=name,color=color)
    ax.set(xticks=x,xticklabels=labels,ylabel='Force component RMSE (eV / Angstrom)',
           title='Same published atomic geometries: no fit, no new DFT')
    ax.legend();ax.grid(axis='y',alpha=.2)
    fig.savefig(args.output/'same_geometry_forces.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
