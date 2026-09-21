"""Plot saved observations and explicit failure diagnostics; no new simulations."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('results/silicon_thermal_v4'))
    args=parser.parse_args();root=args.root
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(15,8.8),layout='constrained')
    reference=read(root/'references/front4/summary.json')['rows']
    old=sorted(reference,key=lambda r:r['gap_A'])
    new=sorted(read(root/'reconstructed_branch/summary.json')['rows'],key=lambda r:r['gap_A'])
    ax=axes[0,0]
    ax.plot([r['gap_A'] for r in old],[r['energy_difference_eV'] for r in old],'-o',ms=3,label='Original conditional minimum')
    ax.plot([r['gap_A'] for r in new],[r['energy_from_original_initial_eV'] for r in new],'-o',ms=3,label='Rearranged conditional minimum')
    ax.set(xlabel='Selected gap q (Angstrom)',ylabel='Energy from original state (eV)',title='A. Two stable bath structures at the same q')
    ax.legend(fontsize=8)
    ax=axes[0,1]
    first=read(root/'first_local_reconstruction15/summary.json')['saddle_audit']['energy_from_initial_eV']
    second=read(root/'second_local_continuation/summary.json')['stationary_audits'][0]['energy_from_original_eV']
    third=next(r['energy_from_original_eV'] for r in read(root/'localized_reconstruction_second/summary.json')['rows'] if r.get('index')==1)
    values=[0,first,.11637498949198566,second,-.3167982195313912,third,-2.213422554312857]
    ax.plot(range(7),values,'-',color='.65',lw=1)
    ax.scatter([0,2,4,6],np.array(values)[[0,2,4,6]],color='#286b8e',label='Conditional minima')
    ax.scatter([1,3,5],np.array(values)[[1,3,5]],marker='^',color='#c85d2d',label='Index-one stationary saddles')
    ax.set(xlabel='Stationary sequence (schematic connections)',ylabel='Energy from original state (eV)',title='B. Sequential rearrangement at fixed q')
    ax.text(.03,.42,'Forward local barriers (eV)\n0.299645 / 0.049806 / 0.006078',transform=ax.transAxes,fontsize=9)
    ax.legend(fontsize=8,loc='lower left')
    ax=axes[0,2]
    releases=read(root/'reconstruction_dynamics/summary.json')['records']
    for r in sorted(releases,key=lambda r:(r['temperature_K'],r['seed'])):
        data=np.load(root/'reconstruction_dynamics'/(r['label']+'.npz'),allow_pickle=False)
        order=data['reconstruction_fraction_every_step'];t=np.arange(len(order))*.001
        ax.plot(t[::100],order[::100],lw=.9,color='#21718a' if r['temperature_K']==300 else '#d4762c',
                alpha=.75,label=f"{r['temperature_K']:g} K / seed {r['seed']}")
    ax.axhline(1,color='.5',ls=':',lw=.8)
    ax.set(xlabel='Actual reference MD time (ps)',ylabel='Rearrangement projection / path length',title='C. New metastable structural releases')
    ax.legend(fontsize=7)
    data=np.load(root/'thermal_dynamics/correlations.npz',allow_pickle=False)
    ax=axes[1,0]
    for temperature,color in [(100,'#4f8f7e'),(300,'#ba653c')]:
        label=f'initial_T{temperature}_dt0.00100';t=data[label+'_time_ps'];c=data[label+'_correlations']
        mask=t<=.5
        ax.plot(t[mask],c.mean(0)[mask],color=color,label=f'{temperature} K, four preparations')
    harmonic=np.load(root/'memory_quadrature/initial_histories.npz',allow_pickle=False)
    mask=harmonic['time_ps']<=.5
    ax.plot(harmonic['time_ps'][mask],harmonic['exact_memory_eV_A2'][mask],'--',color='.25',lw=1,label='0 K harmonic kernel')
    ax.axhline(0,color='.7',lw=.6)
    ax.set(xlabel='Lag (ps)',ylabel='Force correlation / kBT (eV/Angstrom²)',title='D. Nonlinear frozen-gap force correlations')
    ax.legend(fontsize=8)
    ax=axes[1,1]
    label='initial_T300_dt0.00100';t=data[label+'_time_ps'];curves=data[label+'_integrals']
    ax.fill_between(t[::10],curves.min(0)[::10],curves.max(0)[::10],color='#c5dce2',label='Observed replica range, not a CI')
    ax.plot(t[::10],curves.mean(0)[::10],color='#20677c',label='300 K mean, front 4')
    ax.axhline(0,color='.3',ls=':',lw=1)
    ax.set(xlabel='Integration cutoff (ps)',ylabel='Finite integral (eV ps/Angstrom²)',title='E. No constant drag selected')
    ax.legend(fontsize=8)
    ax=axes[1,2]
    records=[r for r in read(root/'memory_quadrature/summary.json')['records'] if r['state']=='initial']
    for window in ['1.0','5.0','20.0']:
        ax.loglog([r['rank'] for r in records],[max(r['response_max_errors'][window],1e-16) for r in records],'-o',ms=3,label=window+' ps window')
    ax.set(xlabel='Retained harmonic bath oscillators',ylabel='Max normalized gap-response error',title='F. Short-time compression is not long-time closure')
    ax.legend(fontsize=8)
    fig.suptitle('Si v4: pure original SW / fixed grips / research reference\nNo wafer-strength, fatigue-life or production-Hz calibration',fontsize=14)
    fig.savefig(root/'thermal_research.png',dpi=170);plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(11.5,8),layout='constrained')
    records=read(root/'front8_dynamics/summary.json')['records']
    for row in sorted(records,key=lambda r:(r['temperature_K'],r['seed'])):
        data=np.load(root/'front8_dynamics'/(row['label']+'.npz'),allow_pickle=False)
        col=0 if row['temperature_K']==100 else 1
        color='#b7302e' if row['reflection_count'] else None
        label=f"seed {row['seed']} / {row['reflection_count']} reflections"
        axes[0,col].plot(data['time_ps'][::20],data['observations'][::20,1],lw=.9,color=color,label=label)
        axes[1,col].plot(data['time_ps'][::20],data['energy_residual_eV'][::20]/row['initial_kinetic_energy_eV'],lw=.9,color=color)
    for col,temperature in enumerate([100,300]):
        axes[0,col].axhline(1,color='.4',ls=':',lw=1)
        axes[0,col].set(title=f'Front 8, {temperature} K: all four records retained',ylabel='Max absolute bath coordinate (Angstrom)',xlabel='MD time (ps)')
        axes[0,col].legend(fontsize=8)
        axes[1,col].set(xlabel='MD time (ps)',ylabel='Energy residual / initial kinetic energy')
        axes[1,col].axhline(0,color='.7',lw=.5)
    fig.suptitle('Finite-size control: boundary and energy-error failure remains visible\nThe reflected 300 K record is not used to certify unconfined dynamics',fontsize=13)
    fig.savefig(root/'boundary_diagnostic.png',dpi=170);plt.close(fig)
    print('saved thermal_research.png and boundary_diagnostic.png',flush=True)


if __name__=='__main__':main()
