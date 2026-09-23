"""Plot saved static loading/unloading structures; no new model evaluation."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    root = args.results
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    used = []
    def use(relative):
        path = root / relative
        used.append(dict(path=relative, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        return path
    def read(relative):
        return json.loads(use(relative).read_text(encoding='utf-8'))
    audit = read('prism_unload_replay/summary.json')
    if not audit['loading_converged'] or not audit['unloading_converged']:
        raise ValueError('converged comparison required for this figure')
    paths = {
        'baseline_raw': 'intact_prism_360_linesearch_v2/state_004/raw.npz',
        'unload_raw': 'prism_unload_10_to_8/raw.npz',
        'geometry': 'intact_prism_360_linesearch_v2/geometry.npz',
    }
    arrays = {}
    for role, relative in paths.items():
        path = use(relative)
        expected = next(s['sha256'] for s in audit['sources'] if s['role'] == role)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError('audit source changed: ' + role)
        with np.load(path) as data:
            arrays[role] = {key: data[key].copy() for key in data.files}
    r0, r1 = arrays['baseline_raw']['positions'], arrays['unload_raw']['positions']
    free = arrays['geometry']['free']
    distance = np.linalg.norm(r1[free]-r0[free], axis=1)
    rms = float(np.sqrt(np.mean(distance**2)))
    if abs(rms-audit['labelled_free_atom_RMS_difference_A']) > 1e-12:
        raise ValueError('geometry descriptor changed')
    curve = []
    excluded = 0
    for path in sorted((root/'intact_prism_360_linesearch_v2').glob('state_*/result.json')):
        row = read(path.relative_to(root).as_posix())
        if row['converged']:
            curve.append((100*row['strain'], row['nominal_stress_GPa']))
        else:
            excluded += 1
    continued = read('prism_10pct_continuation/summary.json')
    if not continued['converged']:
        raise ValueError('converged continuation required')
    curve.append((100*continued['strain'], continued['nominal_stress_GPa']))
    curve = np.asarray(curve)
    fig, axes = plt.subplots(1, 3, figsize=(11, 6), gridspec_kw={'width_ratios':[.9,.9,1.3]})
    for axis, coordinate in zip(axes[:2], (0,1)):
        axis.scatter(r0[~free,coordinate],r0[~free,2],s=12,c='#9aa8b4',label='Fixed grips')
        axis.scatter(r0[free,coordinate],r0[free,2],s=22,facecolors='none',edgecolors='#a0a0a0',
                     linewidths=.6,label='First loading to 8%')
        sc=axis.scatter(r1[free,coordinate],r1[free,2],s=15,c=distance,cmap='viridis',
                        vmin=0,vmax=distance.max(),label='10% then back to 8%')
        axis.set_aspect('equal');axis.set_xlabel(('x','y')[coordinate]+' (Angstrom)')
        axis.set_ylabel('z (Angstrom)')
        axis.set_title('Same fixed grips; '+('x-z','y-z')[coordinate]+' view',fontsize=10)
    axis=axes[2]
    axis.plot(curve[:,0],curve[:,1],'o-',c='#267493',label='Loading: converged states')
    axis.plot([curve[-1,0],100*audit['target_strain']],
              [curve[-1,1],audit['unloading_stress_GPa']],'s--',c='#c96e29',label='Static return to 8%')
    axis.set_xlabel('Imposed strain from initial reference (%)')
    axis.set_ylabel('Nominal axial stress (GPa)')
    axis.grid(alpha=.2);axis.legend(fontsize=8,loc='upper left')
    axis.set_title('Saved endpoints; lines only guide the eye',fontsize=10)
    axes[0].legend(fontsize=7,loc='lower left')
    fig.subplots_adjust(left=.055,right=.98,bottom=.30,top=.86,wspace=.35)
    caxis=fig.add_axes([.09,.12,.43,.025])
    fig.colorbar(sc,cax=caxis,orientation='horizontal',label='Labelled free-atom position difference (Angstrom)')
    fig.suptitle('Static preparation dependence at identical 8% boundary displacement',fontsize=13,y=.97)
    fig.text(.69,.13,f'Energy difference at same grips: {audit["energy_difference_at_same_grips_eV"]:.5f} eV\n'
             f'RMS position difference: {rms:.5f} Angstrom',ha='center',fontsize=9)
    fig.text(.5,.025,'Geometric and static-force comparison only; no physical time, dissipated-work, plasticity, or first-crack claim.',
             ha='center',fontsize=8)
    args.output.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output/'static_return.png',dpi=180)
    plt.close(fig)
    result=dict(source_files=used,loading_points=len(curve),excluded_unconverged_loading_states=excluded,
                matched_atom_permutation_count=audit['minimum_matching_permuted_atoms'],
                labelled_RMS_A=rms,energy_difference_at_same_grips_eV=audit['energy_difference_at_same_grips_eV'],
                new_potential_calls=0,new_DFT=0,new_MD=0,
                scope='saved static endpoints and raw geometry; no first-crack or physical-cycle inference')
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
