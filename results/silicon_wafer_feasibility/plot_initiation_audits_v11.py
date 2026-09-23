"""Figures for completed paper/material/mathematical audits, not Si lifetimes."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def rows(path):
    with path.open(encoding='utf-8') as f:return list(csv.DictReader(f))


def main(root,out):
    if out.exists() and any(out.iterdir()):raise ValueError('fresh figure folder required')
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120})
    source=json.loads((root/'thermal_oxide_source.json').read_text(encoding='utf-8'))
    data=source['rows'];geometry=rows(root/'thermal_oxide_audit/geometry_conventions.csv')
    stress=rows(root/'thermal_oxide_audit/axial_stress_partition.csv')
    fig,axes=plt.subplots(1,3,figsize=(14,4.4),layout='constrained')
    labels=['Bare','50 nm','100 nm','200 nm','Oxide\nremoved']
    strength=[r['nominal_strength_GPa'] for r in data]
    bars=axes[0].bar(labels,strength,color=['#64748b','#3c85ad','#3c85ad','#3c85ad','#a76d56'])
    axes[0].bar_label(bars,fmt='%.2f',padding=3);axes[0].set_ylim(0,4.7)
    axes[0].set_ylabel('Reported nominal fracture strength (GPa)');axes[0].set_title('Measured final fracture, not onset timing')
    for convention,label,style in [('listed_outer_per_face_t','Listed outer; t per face','-o'),
        ('listed_outer_total_t','Listed outer; t across both faces','--s'),
        ('planar_44percent_consumption_per_face_t','Planar 44% consumption scenario',':^')]:
        selected=[r for r in geometry if r['convention']==convention]
        axes[1].plot([50,100,200],[float(r['modulus_GPa']) for r in selected],style,label=label)
    axes[1].scatter([50,100,200],[r['E_theory_GPa'] for r in data[1:4]],marker='x',s=60,color='black',label='Paper theoretical values',zorder=5)
    axes[1].set_xlabel('Nominal oxide thickness (nm)');axes[1].set_ylabel('Axial composite modulus (GPa)');axes[1].set_title('Thickness convention changes the core');axes[1].legend(fontsize=7,loc='best')
    selected=[r for r in stress if r['convention']=='listed_outer_per_face_t' and r['residual_model']=='oxide_eigenstress_at_zero_common_strain']
    axes[2].plot([50,100,200],[r['FEM_Si_axial_GPa'] for r in data[1:4]],'ok-',label='Paper FEM Si point stress')
    axes[2].plot([50,100,200],[float(r['Si_axial_GPa']) for r in selected],'s--',label='1D phase stress audit')
    axes[2].set_xlabel('Nominal oxide thickness (nm)');axes[2].set_ylabel('Si axial stress (GPa)');axes[2].set_title('Measured fracture force supplied as input');axes[2].legend(fontsize=8)
    fig.suptitle('Thermal-oxide paper audit: geometry and stress bookkeeping',fontsize=14)
    fig.savefig(out/'thermal_oxide_audit.png',dpi=180);plt.close(fig)

    source_rows=rows(root/'oxide_mace_mtpu_1159_v2/frames.csv');groups=rows(root/'oxide_mace_mtpu_1159_v2/groups.csv')
    order=['Si_only','SiO2_stoichiometry','mixed_other_stoichiometry','O_only']
    colors=dict(zip(order,['#4d7192','#55a0a0','#bc7658','#8b7ca9']))
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained');lookup={r['group']:r for r in groups}
    x=np.arange(4);names=['Si only','SiO2 ratio','Other Si/O','O only']
    axes[0].bar(x-.16,[float(lookup[k]['force_component_RMSE_eV_A']) for k in order],.32,label='Component RMSE',color='#446b89')
    axes[0].bar(x+.16,[float(lookup[k]['force_component_MAE_eV_A']) for k in order],.32,label='Component MAE',color='#a3becf')
    axes[0].set_xticks(x,names);axes[0].set_ylabel('Force error (eV / Angstrom)');axes[0].legend();axes[0].set_title('All 1,159 source configurations retained')
    zero=0
    for group in order:
        subset=[r for r in source_rows if r['group']==group]
        a=np.array([float(r['source_force_component_RMS_eV_A']) for r in subset]);b=np.array([float(r['force_component_RMSE_eV_A']) for r in subset])
        mask=(a>0)&(b>0);zero+=int((~mask).sum())
        axes[1].scatter(a[mask],b[mask],s=12,alpha=.6,label=names[order.index(group)],color=colors[group])
    axes[1].set_xscale('log');axes[1].set_yscale('log');axes[1].set_xlabel('DFT force component RMS (eV / Angstrom)')
    axes[1].set_ylabel('MACE force component RMSE (eV / Angstrom)');axes[1].legend(fontsize=8)
    axes[1].set_title(f'Per-frame error; {zero} exact-zero points outside log axes')
    fig.suptitle('Neutral MACE versus public PBE Si/O data: material audit only',fontsize=13)
    fig.savefig(out/'oxide_force_audit.png',dpi=180);plt.close(fig)

    probability=rows(root/'first_passage_math/synthetic_first_passage.csv');refinement=rows(root/'first_passage_math/diffusion_refinement.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    t=np.array([float(r['time_model_units']) for r in probability]);a=np.array([float(r['surface_demo_CDF']) for r in probability]);b=np.array([float(r['interface_demo_CDF']) for r in probability])
    axes[0].plot(t,a+b,label='First entry to defined B',color='#234f70',lw=2)
    axes[0].plot(t,a,'--',label='Cause 1');axes[0].plot(t,b,'--',label='Cause 2')
    axes[0].plot(t,[float(r['premature_intermediate_CDF']) for r in probability],':',label='Premature middle-state entry',color='#ad654a',lw=2)
    axes[0].set_xlim(0,10);axes[0].set_ylim(0,1.03);axes[0].set_xlabel('Model time');axes[0].set_ylabel('Synthetic first-passage probability');axes[0].legend(fontsize=8)
    h=np.array([float(r['spacing']) for r in refinement]);error=np.array([float(r['error']) for r in refinement])
    axes[1].loglog(h,error,'o-',label='SG first-passage error');axes[1].loglog(h,error[-1]*(h/h[-1])**2,':',label='Second-order reference')
    axes[1].set_xlabel('Grid spacing');axes[1].set_ylabel('Survival error versus analytic diffusion');axes[1].legend(fontsize=8)
    fig.suptitle('Probability-accounting validation with synthetic rates; no Si kinetics',fontsize=13)
    fig.savefig(out/'first_passage_validation.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.results,a.output)
