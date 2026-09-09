"""Fixed-range necessary stability cuts; unchanged static target loss."""
import json
import time

import numpy as np

from .spectral_material_constraints import BulkHessianCoefficientBasis,stability_halfspace,fit_spectral_faces
from .profiled_range_calibration import fit_profiled_coefficients
from .range_resolved_material import RangeObservationCache,build_range_surface,RangeBulkValidation
from .run_low_stress_cyclic_diagnostic import ROOT,write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .yield_elastic_metric import cubic_metric_problem
from .vector_registry_audit import stationary_state
from .vector_material_calibration import UNITS,IDEAL_H


OUT=ROOT/'results/fcc111_active_interface/stable_core_v13/material_spectral'


def wavepoints():
    return [(path,float(fraction),np.array(end)*fraction)
        for path,end in [('Gamma_X',[1,0,0]),('Gamma_L',[.5,.5,.5]),('Gamma_K',[.75,.75,0])]
        for fraction in (.02,.1,.25,.5,.75,1.)]


def main():
    if (OUT/'profiles.json').exists():
        raise FileExistsError('preserve completed spectral study')
    began=time.perf_counter(); source,observations,states=source_and_targets()
    old=json.loads((ROOT/'results/fcc111_active_interface/range_core_v12/material/calibration.json').read_bytes())
    new=json.loads((OUT.parent/'material_reprofile/calibration.json').read_bytes())
    starts=dict(v12_same_ranges=old['best']['cubic_5pct'],v13_other_basin=new['best'])
    cache=RangeObservationCache(observations); histories={}; summaries=[]; waves=[]; residuals=[]; roots=[]
    for name,start in starts.items():
        decays=[start[k] for k in ('scalar_decay','odd_decay','quadrupole_decay')]
        matrix,obs=cubic_metric_problem(cache.matrix(decays),observations)
        basis=BulkHessianCoefficientBasis(decays,radius=8.)
        points=wavepoints(); columns=[basis.matrices(q) for _,_,q in points]
        cuts=[]; history=[]; fit=fit_profiled_coefficients(matrix,obs)
        for iteration in range(7):
            probes=[stability_halfspace(C,fit['coefficients']) for C in columns]
            index=min(range(len(probes)),key=lambda k:probes[k]['minimum_eigenvalue'])
            lowest=probes[index]; path,fraction,q=points[index]
            # This is a linear-algebra roundoff stopping floor, not physical
            # stability certification. Independent radius errors are tested below.
            numerical_floor=256*np.finfo(float).eps*np.sum(abs(fit['coefficients']))*max(
                np.max(abs(C)) for C in columns)
            history.append(dict(iteration=iteration,fit=fit,worst_path=path,worst_fraction=fraction,
                minimum_eigenvalue=lowest['minimum_eigenvalue'],roundoff_floor=numerical_floor,
                separating_polarization=lowest['polarization'],halfspace=lowest['halfspace']))
            print(f"{name} cut {iteration}: loss={fit['squared_loss']:.7g}, minH={lowest['minimum_eigenvalue']:.7g}",flush=True)
            save_json(OUT/f'{name}_progress.json',dict(completed=False,history=history,decays=decays))
            if lowest['minimum_eigenvalue']>=-numerical_floor or iteration==6:
                break
            cuts.append(lowest['halfspace'])
            fit=fit_spectral_faces(matrix,obs,cuts)
        histories[name]=dict(decays=decays,history=history,cut_count=len(cuts),final_fit=fit)
        final=build_range_surface(*decays,fit['coefficients'])
        previous={}; minimums={}; maximum_change=0.
        for radius in (6.,8.,10.):
            check=RangeBulkValidation(final,cutoff=radius); current={}
            for path,fraction,q in points:
                key=path,fraction; value=check.evaluate(check.crystallographic_wavevector(q))
                current[key]=value['matrix']; change=(float(np.linalg.norm(current[key]-previous[key],2))
                    if key in previous else None)
                waves.append(dict(candidate=name,radius_over_L0=radius,path=path,fraction=fraction,
                    minimum_eigenvalue_eV_L0sq=value['eigenvalues'][0],matrix_spectral_change=change))
            previous=current; minimums[radius]=min(np.linalg.eigvalsh(v)[0] for v in current.values())
            if radius==10.:
                maximum_change=max(r['matrix_spectral_change'] for r in waves if r['candidate']==name and r['radius_over_L0']==10.)
        hold=[i for i,o in enumerate(obs) if o.role=='heldout']
        res=(matrix@fit['coefficients']-np.array([o.target for o in obs]))/np.array([o.scale for o in obs])
        summaries.append(dict(candidate=name,fit_loss=fit['squared_loss'],cuts=len(cuts),
            **dict(zip(['C11_GPa','C12_GPa','C44_GPa'],matrix[2:5]@fit['coefficients'])),
            heldout_normalized_rms=float(np.sqrt(np.mean(res[hold]**2))),
            minimum_radius8_eigenvalue=minimums[8.],minimum_radius10_eigenvalue=minimums[10.],
            maximum_radius8_to10_matrix_change=maximum_change,
            positive_sampled_spectrum_margin=minimums[10.]>maximum_change,
            material_accepted=False))
        residuals.extend(dict(candidate=name,observable=o.name,role=o.role,units=o.units,
            target=o.target,prediction=value,normalized_residual=r)
            for o,value,r in zip(obs,matrix@fit['coefficients'],res))
        for label,guess,index in [('perfect',[IDEAL_H,0.,0.],0),('fault',states['fault'],0),('saddle',states['saddle'],1)]:
            root=stationary_state(final,guess,expected_index=index)
            roots.append(dict(candidate=name,state=label,valid=root['valid'],
                energy_J_m2=float(UNITS.energy_to_surface(root['evaluation'].energy)),
                a=root['q'][0],x=root['q'][1],y=root['q'][2],
                force_residual=root['force_residual'],eigenvalues=root['eigenvalues'].tolist()))
        save_json(OUT/f'{name}_progress.json',dict(completed=True,history=history,decays=decays))
    write_csv(OUT/'summary.csv',summaries); write_csv(OUT/'finite_q_refinement.csv',waves)
    write_csv(OUT/'residuals.csv',residuals); write_csv(OUT/'stationary_states.csv',roots)
    save_json(OUT/'profiles.json',dict(completed=True,histories=histories,
        elapsed_seconds=time.perf_counter()-began,source_sha256=source.reference.sha256,
        energy_family_changed=False,targets_scales_changed=False,
        radial_optimization_under_spectral_constraints=False,
        sampled_psd_is_not_global_stability=True,material_accepted=False,physical_yield_validated=False))


if __name__=='__main__':
    main()
