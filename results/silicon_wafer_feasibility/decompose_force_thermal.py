"""Cached Cartesian component audit of the boundary-safe canonical identity.

The three axes are the source Cartesian frame, not assigned cleavage/slip
directions. Correlated block diagnostics are not confidence intervals.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann,electron_volt
from solver_v1.silicon_thermal_identity import canonical_temperature_terms
from solver_v1.silicon_thermal_research import block_statistics

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def arrays(path):
    with np.load(path) as z:return {k:z[k].copy() for k in z.files}
def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def main(a):
    if a.output.exists():raise ValueError('fresh output required')
    a.output.mkdir(parents=True)
    source_paths=[a.root/f'results/silicon_initiation_v13/dense_{name}/raw_hessian.npz' for name in ('loading8','return8')]
    sources=[arrays(p) for p in source_paths]
    common=(sources[0]['positions']+sources[1]['positions'])/2
    kT=Boltzmann/electron_volt*300.
    results=[];saved={};inputs={};maximum_sum_error=0.
    for method,pilot,analysis in [
        ('pcn',a.root/'results/silicon_crack_thermal_validation/nonlinear_pilot',a.root/'results/silicon_crack_thermal_validation/analysis'),
        ('hmc',a.root/'results/silicon_force_thermal_comparison/atomistic_pilot',a.root/'results/silicon_force_thermal_comparison/analysis')]:
        summary=json.loads((analysis/'summary.json').read_text());inputs[str((analysis/'summary.json').relative_to(a.root)).replace('\\','/')]=sha(analysis/'summary.json')
        for chain,source in enumerate(sources):
            folder=pilot/f'chain_{chain}';by_index={int(p.stem.split('_')[1]):p for p in folder.glob('evaluation_*.npz')}
            free=np.any(source['basis'].reshape(len(source['positions']),3,-1)!=0,axis=(1,2))
            residual=[];num=[];div=[]
            for index in summary['chains'][chain]['retained_evaluation_indices']:
                path=by_index[index];value=arrays(path)
                inputs[str(path.relative_to(a.root)).replace('\\','/')]=sha(path)
                term=canonical_temperature_terms(value['positions'],value['forces'],free,common,source['positions'],
                    temperature_K=300.,halfwidth_A=3.,minimum_pair_A=1.,taper_width_A=.4)
                component=term['numerator_components_eV']-kT*term['divergence_components']
                maximum_sum_error=max(maximum_sum_error,abs(float(component.sum())-term['residual_eV']))
                residual.append(component);num.append(term['numerator_components_eV']);div.append(term['divergence_components'])
            residual=np.asarray(residual);axes=residual.sum(axis=1);block=block_statistics(axes,blocks=8)
            mean=axes.mean(axis=0);per_atom=residual.mean(axis=0)
            if abs(mean.sum()-summary['chains'][chain]['mean_identity_residual_eV'])>1e-11:
                raise ValueError('component sum does not reproduce audited global residual')
            results.append(dict(method=method,chain=chain,correlated_draws=len(residual),
                Cartesian_xyz_mean_residual_eV=mean.tolist(),
                Cartesian_xyz_block_SE_diagnostic_eV=block['block_standard_error'].tolist(),
                Cartesian_xyz_first_half_residual_eV=block['first_half_mean'].tolist(),
                Cartesian_xyz_second_half_residual_eV=block['second_half_mean'].tolist(),
                global_mean_residual_eV=float(mean.sum()),
                sum_absolute_atom_axis_mean_residuals_eV=float(abs(per_atom).sum()),
                physical_axis_assignment_verified=False,confidence_interval_claim=False,
                per_component_equilibrium_certified=False))
            saved[f'{method}_chain{chain}_component_residuals_eV']=residual
            saved[f'{method}_chain{chain}_numerator_components_eV']=np.asarray(num)
            saved[f'{method}_chain{chain}_divergence_components']=np.asarray(div)
    if maximum_sum_error>1e-11:raise ValueError('thermal decomposition mismatch')
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),temperature_K=300.,
        Cartesian_axis_names=['x','y','z'],results=results,maximum_sum_error_eV=maximum_sum_error,
        identity='each g_i vanishes at its Cartesian face and every pair wall; integrate each d_i(g_i exp(-U/kBT))',
        aggregate_cancellation_is_not_per_component_validation=True,new_model_calls=0,new_MD=0,new_DFT=0,
        kBT_fitted=False,equilibrium_certified=False,production_changed=False,
        source_sha256={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in source_paths+
            [a.root/'solver_v1/silicon_thermal_identity.py']},runner_sha256=sha(Path(__file__)))
    np.savez_compressed(a.output/'component_terms.npz',**saved)
    dump(a.output/'summary.json',result);dump(a.output/'input_manifest.json',inputs)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(8.6,3.5),layout='constrained')
    for chain,axis in enumerate(axes):
        for offset,method in [(-.17,'pcn'),(.17,'hmc')]:
            row=next(r for r in results if r['method']==method and r['chain']==chain)
            axis.bar(np.arange(3)+offset,row['Cartesian_xyz_mean_residual_eV'],width=.32,
                     yerr=row['Cartesian_xyz_block_SE_diagnostic_eV'],label=method.upper(),capsize=3)
        axis.set(xticks=np.arange(3),xticklabels=['x','y','z'],ylabel='Mean residual (eV)',
                 title=['Initial 8% start','Returned 8% start'][chain])
        axis.axhline(0,color='black',lw=.7);axis.grid(axis='y',alpha=.2);axis.legend()
    fig.suptitle('Cartesian thermal-identity components; error bars are block diagnostics, NOT CI',fontsize=9)
    fig.savefig(a.output/'component_thermal_diagnostics.png',dpi=150);plt.close(fig)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
