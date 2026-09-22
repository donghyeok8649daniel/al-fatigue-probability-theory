"""Plot measured numerical controls, not a fabricated transition path."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(root):
    endpoint=json.loads((root/'boron_endpoint_validation.json').read_text())
    mode=json.loads((root/'mace_boron_largest_mode/summary.json').read_text())
    saddle=json.loads((root/'mace_boron_hessian/3boron_int_cluster_confi_2_hessian.json').read_text())
    fig,axes=plt.subplots(1,3,figsize=(13,4.5),layout='constrained')
    values=[saddle['Cartesian_projected_min_eigenvalue_eV_A2']]+[r['min_curvature_eV_A2'] for r in endpoint['corrected_endpoints']]
    axes[0].bar(np.arange(3),values,color=['#ac413b','#3274a1','#3274a1'])
    axes[0].axhline(0,color='#777',lw=.8)
    axes[0].set(xticks=np.arange(3),xticklabels=['Original\n3B-2','Minus\nendpoint','Plus\nendpoint'],
        ylabel='Lowest fixed-cell curvature (eV / Angstrom²)',title='A force-converged saddle was detected')
    for i,v in enumerate(values):axes[0].annotate(f'{v:+.5f}',(i,v),xytext=(0,5 if v>0 else -15),textcoords='offset points',ha='center',fontsize=9)
    axes[0].set_ylim(-.2,.5)
    association=endpoint['association'];x=np.arange(len(association))
    y=[r['association_energy_eV'] for r in association]
    axes[1].plot(x,y,'o-',color='#3274a1')
    axes[1].set(xticks=x,xticklabels=[str(r['host_atoms_per_cell']) for r in association],
        xlabel='Si atoms per periodic cell',ylabel='Neutral association energy (eV)',
        xlim=(-.25,2.25),ylim=(0,.65),title='Endpoint energy; not a kinetic barrier')
    for i,v in enumerate(y):axes[1].annotate(f'{v:.5f}',(i,v),xytext=(0,7),textcoords='offset points',ha='center',fontsize=9)
    frequencies=[r['highest_Gamma_frequency_cm_minus1'] for r in mode['cases']]
    axes[2].plot(x,frequencies,'o-',color='#3274a1',label='Neutral MACE, Gamma')
    axes[2].axhline(840,color='#8a8f99',ls='--',label='Published LDA, 1B maximum')
    axes[2].set(xticks=x,xticklabels=['64','216','512'],xlabel='Si atoms per periodic cell',
        ylabel='Highest frequency (cm$^{-1}$)',xlim=(-.25,2.25),ylim=(745,853),title='Cell size does not close this gap')
    for i,v in enumerate(frequencies):axes[2].annotate(f'{v:.2f}',(i,v),xytext=(0,8),textcoords='offset points',ha='center',fontsize=9)
    axes[2].legend(fontsize=8,loc='lower right')
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.suptitle('Neutral B / Si controls: fixed volume, no new DFT and no physical clock')
    fig.savefig(root/'boron_controls.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    main(p.parse_args().root)
