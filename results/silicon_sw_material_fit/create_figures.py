"""Scientific figures from the exported raw calculations, without model fits."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    root=args.package;out=root/'figures';out.mkdir(exist_ok=True)
    read=lambda name:json.loads((root/name).read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160})
    colors=['#758497','#146580','#b95e38']
    metrics=read('reaggregated_group_metrics.json')
    def metric(case,kind,field):
        return next(r[field] for r in metrics if r['case']==case and r['config_type']==kind and r['declared_xc']=='PW91')
    cases=['original_SW','angular_shape','joint_static_candidate'];labels=['Original SW','Earlier candidate (C44 low)','Joint elastic candidate']
    fig,ax=plt.subplots(1,2,figsize=(8.4,3.1));x=np.arange(2);width=.24
    for k,(case,label) in enumerate(zip(cases,labels)):
        for axis,field in zip(ax,['force_RMSE_eV_A','relative_energy_RMSE_eV_atom']):
            axis.bar(x+(k-1)*width,[metric(case,t,field) for t in ['dia','surface_111']],width,label=label,color=colors[k])
    for axis in ax:axis.set_xticks(x,['Bulk (training)','Si(111) (excluded from loss)']);axis.grid(axis='y',alpha=.17);axis.set_axisbelow(True)
    ax[0].set_ylabel('Force RMSE (eV / A)');ax[1].set_ylabel('Relative energy RMSE (eV / atom)')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,0,1,.89));fig.savefig(out/'material_errors.png',bbox_inches='tight');plt.close(fig)
    external=read('raw_runs/final_selected_external_results/group_metrics.json');fig,ax=plt.subplots(figsize=(8.2,2.8));x=np.arange(4)
    chosen=next(r['case'] for r in external if r['case'] not in ['original_SW','selected_plain_control'])
    for k,case in enumerate(['original_SW',chosen]):
        values=[next(r['force_RMSE_eV_A'] for r in external if r['case']==case and r['config_type']==t) for t in ['GB','diinter','aSi','SF']]
        ax.bar(x+(k-.5)*.34,values,.34,color=colors[k],label=['Original SW','Joint elastic candidate'][k])
    ax.set_xticks(x,['GB / PW91 (1)','Di-interstitial / PW91 (6)','aSi / PW91 (1)','SF / unspecified XC (20)'])
    ax.set_ylabel('Force RMSE (eV / A)');ax.legend(frameon=False);ax.grid(axis='y',alpha=.17);ax.set_axisbelow(True)
    fig.tight_layout();fig.savefig(out/'external_errors.png',bbox_inches='tight');plt.close(fig)
    band=np.load(root/'raw_runs/final_selected_phonon_recovery_results/final_candidate.npz')
    print('phonon arrays:',{k:band[k].shape for k in band.files})
    fig,ax=plt.subplots(figsize=(8.2,2.9))
    # Raw bands are stored as either (3,50,6) or the concatenated path.
    predicted=band['frequencies_THz'].reshape(-1,6);target=band['reference_THz'].reshape(-1,6)
    for b in range(6):
        ax.plot(target[:,b],color=colors[0],lw=1.1,label='Published DFT' if b==0 else None)
        ax.plot(predicted[:,b],color=colors[1],lw=1.1,label='Final angular candidate' if b==0 else None)
    ax.set_xticks([0,49,99,149],['Gamma','X','Gamma (equivalent)','L']);ax.set_ylabel('Harmonic frequency (THz)')
    ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.15);fig.tight_layout();fig.savefig(out/'phonon_bands.png',bbox_inches='tight');plt.close(fig)
    with (root/'raw_runs/interface_recovery_results/normal_branch.csv').open(newline='') as f:rows=list(csv.DictReader(f))
    with (root/'raw_runs/final_selected_interface_results/normal_branch.csv').open(newline='') as f:rows+=list(csv.DictReader(f))
    fig,ax=plt.subplots(1,2,figsize=(8.2,3.0))
    for k,case in enumerate(['original_SW','selected_joint_candidate']):
        for cut,style in [('shuffle','-'),('glide','--')]:
            group=[r for r in rows if r['case']==case and r['cut_kind']==cut];opening=[float(r['opening_A']) for r in group]
            for a,field in zip(ax,['work_J_m2','normal_traction_GPa']):
                a.plot(opening,[float(r[field]) for r in group],style,color=colors[k],lw=1.2,label=('SW ' if k==0 else 'Candidate ')+cut)
    ax[0].set_ylabel('Rigid work (J / m2)');ax[1].set_ylabel('Rigid normal traction (GPa)')
    for a in ax:a.set_xlabel('Prescribed interface opening (A)');a.grid(alpha=.15)
    ax[0].legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(out/'rigid_interface.png',bbox_inches='tight');plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--package',type=Path,default=Path(__file__).resolve().parent)
    main(parser.parse_args())
