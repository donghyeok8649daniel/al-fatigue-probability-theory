"""Compare finite-band and DHO mobility estimates without certifying either."""
import argparse
import json
from pathlib import Path
import numpy as np
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def summarize(root, output):
    root,output=Path(root),Path(output)
    if output.exists(): raise FileExistsError('fresh output required')
    controls=json.loads((root/'control_comparison/summary.json').read_bytes())
    rows=[];comparison=[]
    for control in controls['controls']:
        tag=control['study']
        spectral=json.loads((root/('spectrum_'+tag)/'summary.json').read_bytes())
        for r in spectral['records']:
            if not r['valid_finite_band']:continue
            M=np.asarray(r['response']['mobility_m2_per_J_second'])
            K=np.asarray(r['band']['integral_proxy_m2_seconds'])
            for axis,label in enumerate(('normal','direct110','transverse')):
                rows.append(dict(study=tag,axis=label,parts=r['parts'],block=r['block'],
                    nw=r['nw'],lower_cycle_THz=r['lower_cycle_THz'],
                    upper_cycle_THz=r['upper_cycle_THz'],mobility_m2_per_J_s=M[axis,axis],
                    integral_proxy_m2_s=K[axis,axis],
                    cell_extensivity_hypothesis_mobility_m2_per_J_s=4*control['planes']**2*M[axis,axis],
                    zero_frequency_certified=False))
        for label in ('normal','direct110','transverse'):
            # Full records and all bands/tapers: no selected favorable band.
            spec=[r['mobility_m2_per_J_s'] for r in rows
                  if r['study']==tag and r['axis']==label and r['parts']==1]
            dho=[r['mobility_m2_per_J_s'] for r in controls['mobility_records']
                 if r['study']==tag and r['axis']==label]
            comparison.append(dict(study=tag,axis=label,
                spectral_full_min=min(spec),spectral_full_max=max(spec),
                dho_min=min(dho),dho_max=max(dho),
                spectral_to_dho_min=min(spec)/max(dho),
                spectral_to_dho_max=max(spec)/min(dho),
                ranges_overlap=bool(min(spec)<=max(dho) and min(dho)<=max(spec))))
    output.mkdir(parents=True)
    write_csv(output/'finite_band_mobility.csv',rows)
    write_csv(output/'spectral_vs_dho.csv',comparison)
    save_json(output/'summary.json',dict(completed=True,comparisons=comparison,
        production_clock_calibrated=False,
        ranges_are_sensitivity_not_confidence_intervals=True))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=[c['study'] for c in controls['controls']]
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for ax,label in zip(axes,('normal','direct110')):
        for index,name in enumerate(names):
            r=next(r for r in comparison if r['axis']==label and r['study']==name)
            for offset,prefix,color in ((-.1,'spectral_full','C0'),(.1,'dho','C1')):
                low,high=r[prefix+'_min'],r[prefix+'_max']
                ax.plot([index+offset]*2,[low,high],color=color,lw=4)
                ax.plot(index+offset,np.sqrt(low*high),'o',color=color)
        ax.set_yscale('log');ax.set_title(label+' plane-gap PMF')
        ax.set_xticks(range(len(names)),[n.replace('1000ps_','') for n in names],rotation=35,ha='right',fontsize=8)
        ax.set_ylabel('Conditional mobility [m²/(J s)]')
        ax.grid(axis='y',alpha=.25)
    axes[0].plot([],[],color='C0',label='Finite-band spectrum: all full-record bands')
    axes[0].plot([],[],color='C1',label='DHO: all windows and temporal halves')
    axes[0].legend(fontsize=8)
    fig.suptitle('Reference-MD estimates; not production mobility certification')
    fig.savefig(output/'mobility_comparison.png',dpi=160)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();summarize(a.root,a.out)
