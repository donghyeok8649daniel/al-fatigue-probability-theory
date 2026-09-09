"""Actual radial reoptimization with a necessary, not sufficient, stability cut.

The negative mode found in the completed unconstrained study supplies one
fixed q/polarization. Its coefficient plane is recomputed for EVERY range.
This is not the fixed-range comparison and not whole-zone PSD certification.
"""
from functools import lru_cache
import json
import time

import numpy as np
from scipy.optimize import least_squares

from .spectral_material_constraints import BulkHessianCoefficientBasis,fit_spectral_faces
from .range_resolved_material import RangeObservationCache,RangeBulkValidation,build_range_surface
from .run_spectral_material_profiles import wavepoints
from .run_low_stress_cyclic_diagnostic import ROOT,write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .yield_elastic_metric import cubic_metric_problem


OUT=ROOT/'results/fcc111_active_interface/stable_core_v13/material_stable_ranges'


def main():
    if (OUT/'calibration.json').exists(): raise FileExistsError('preserve completed research calibration')
    began=time.perf_counter(); source,observations,_=source_and_targets()
    prior=json.loads((OUT.parent/'material_spectral/profiles.json').read_bytes())
    negative=prior['histories']['v13_other_basin']['history'][0]
    v=np.asarray(negative['separating_polarization'])
    points=wavepoints()
    q=next(q for path,fraction,q in points if path==negative['worst_path'] and fraction==negative['worst_fraction'])
    cache=RangeObservationCache(observations); profiles=[]; optimizers=[]
    @lru_cache(maxsize=256)
    def evaluate(logs):
        decays=np.exp(logs); matrix,obs=cubic_metric_problem(cache.matrix(decays),observations)
        columns=BulkHessianCoefficientBasis(decays,radius=8.).matrices(q)
        cut=np.einsum('i,cij,j->c',v,columns,v)
        fit=fit_spectral_faces(matrix,obs,[cut])
        residual=fit['residuals']; hold=[i for i,o in enumerate(obs) if o.role=='heldout']
        row=dict(decays=decays,**fit,cubic_GPa=matrix[2:5]@fit['coefficients'],
            necessary_stability_polarization=float(cut@fit['coefficients']),
            heldout_normalized_rms=float(np.sqrt(np.mean(residual[hold]**2))))
        profiles.append(row)
        if len(profiles)==1 or len(profiles)%10==0:
            print(f"constrained radial {len(profiles)}: loss={fit['squared_loss']:.7g}, decays={decays}",flush=True)
            save_json(OUT/'progress.json',dict(completed=False,profiles=len(profiles),
                best_loss=min(r['squared_loss'] for r in profiles),elapsed_seconds=time.perf_counter()-began))
        return row
    starts=[prior['histories'][key]['decays'] for key in ('v12_same_ranges','v13_other_basin')]
    for start in starts:
        def fun(logs):
            fit=evaluate(tuple(logs)); return fit['residuals'][fit['selected_rows']]
        run=least_squares(fun,np.log(start),jac='2-point',diff_step=2e-4,
            bounds=np.log([[.7,2.,2.],[8.,12.,12.]]),max_nfev=30,ftol=2e-8,xtol=2e-8,gtol=2e-7)
        optimizers.append(dict(start=start,success=bool(run.success),message=str(run.message),
            nfev=run.nfev,njev=run.njev,optimality=float(run.optimality),
            final_decays=np.exp(run.x),final_loss=float(run.fun@run.fun)))
        save_json(ROOT/'.cache/stable_core_v13/material_stable_ranges/optimizer_progress.json',
            dict(completed=False,optimizers=optimizers,profiles=profiles))
    best=min((r for r in profiles if r['strictly_positive_LJ_resolved']),key=lambda r:r['squared_loss'])
    decays=best['decays']; final=build_range_surface(*decays,best['coefficients'])
    waves=[]; previous={}; minima=[]; changes=[]
    for radius in (6.,8.,10.):
        check=RangeBulkValidation(final,cutoff=radius); current={}; local=[]
        for path,fraction,qpoint in points:
            key=path,fraction; value=check.evaluate(check.crystallographic_wavevector(qpoint))
            current[key]=value['matrix']; local.append(value['eigenvalues'][0])
            change=float(np.linalg.norm(current[key]-previous[key],2)) if key in previous else None
            waves.append(dict(radius_over_L0=radius,path=path,fraction=fraction,
                minimum_H_eV_L0sq=value['eigenvalues'][0],matrix_radius_change=change))
            if radius==10.: changes.append(change)
        minima.append(min(local)); previous=current
    matrix,obs=cubic_metric_problem(cache.matrix(decays),observations)
    write_csv(OUT/'residuals.csv',[dict(observable=o.name,role=o.role,units=o.units,
        target=o.target,prediction=p,scale=o.scale,normalized_residual=r)
        for o,p,r in zip(obs,matrix@best['coefficients'],best['residuals'])])
    write_csv(OUT/'finite_q_refinement.csv',waves)
    save_json(OUT/'calibration.json',dict(completed=True,best=best,profiles=profiles,optimizers=optimizers,
        necessary_constraint_q_cubic=q,necessary_constraint_polarization=v,constraint_radius_over_L0=8.,
        minimum_radius6_8_10_H=minima,maximum_radius8_to10_matrix_change=max(changes),
        sampled_stability_above_radius_change=bool(minima[-1]>max(changes)),
        constraint_is_necessary_not_sufficient=True,energy_family_changed=False,targets_scales_changed=False,
        source_sha256=source.reference.sha256,elapsed_seconds=time.perf_counter()-began,
        material_accepted=False,physical_yield_validated=False,physical_Hz=False))
    save_json(OUT/'progress.json',dict(completed=True,profiles=len(profiles),best_loss=best['squared_loss']))
    print(f"completed constrained ranges: loss={best['squared_loss']:.7g}; C={best['cubic_GPa']}; minH={minima[-1]:.6g}",flush=True)


if __name__=='__main__':
    main()
