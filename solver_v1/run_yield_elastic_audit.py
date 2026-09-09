"""Run an explicit metric/elastic-compatibility study before source strength.

No yield data enter fitting; no new energy terms or production-model change.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from .run_low_stress_cyclic_diagnostic import FIT, write_csv
from .run_vector_registry_audit import save_json
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import VectorCoefficientBasis, observation_matrix, fit_coefficients, coefficient_identifiability
from .yield_elastic_metric import cubic_metric_problem, fit_exact_bulk, transformed_mode_covariance

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fcc111_active_interface/yield_bridge_v11/material_metric_refined'


def execute(out, local_evaluations):
    started=time.perf_counter(); source,obs,states=source_and_targets()
    old=json.loads(FIT.read_bytes())['angular_monotone_opening']
    prior=json.loads((ROOT/'results/fcc111_active_interface/material_strength_v10/opening_exchange/calibration.json').read_text())['best']
    profiles=[]; cache={}; logs=[]
    def evaluate(d,k,label,which=None):
        key=float(d),float(k)
        if key not in cache:
            base=observation_matrix(VectorCoefficientBasis(*key),obs)
            cm,co=cubic_metric_problem(base,obs)
            cache[key]=base,cm,co
        base,cm,co=cache[key]
        for mode in (['cubic_5pct','bulk_exact'] if which is None else [which]):
            fit=fit_coefficients(cm,co) if mode=='cubic_5pct' else fit_exact_bulk(cm,co)
            c=fit['coefficients']
            row=dict(scalar_decay=key[0],angular_decay=key[1],metric=mode,stage=label,**fit)
            if c is not None:
                row['cubic_GPa']=cm[2:5]@c
                row['heldout_normalized_rms']=float(np.sqrt(np.mean([
                    ((m@c-o.target)/o.scale)**2 for m,o in zip(base,obs) if o.role=='heldout'])))
            profiles.append(row)
            print(f"{len(profiles)} {mode} {key}: loss={row.get('squared_loss')} feasible={row['admissible']}",flush=True)
        save_json(out/'progress.json',dict(completed=False,profiles=profiles,elapsed_seconds=time.perf_counter()-started))
        return row
    keys=[(old['scalar_decay'],old['angular_decay']),(prior['scalar_decay'],prior['angular_decay'])]
    keys+= [(d,k) for d in (.7,1.287051801314886,2.4,4.5,8.) for k in (2.,3.2,4.898979485566356,7.8,12.)]
    keys=list(dict.fromkeys(keys))
    for d,k in keys:
        evaluate(d,k,'predeclared_grid')
    for metric in ('cubic_5pct','bulk_exact'):
        eligible=[r for r in profiles if r['metric']==metric and r['admissible'] and r['strictly_positive_LJ_resolved']]
        if not eligible:
            logs.append(dict(metric=metric,no_admissible_start=True)); continue
        best=min(eligible,key=lambda r:r['squared_loss'])
        if local_evaluations:
            start=[best['scalar_decay'],best['angular_decay']]
            def objective(x):
                r=evaluate(*x,'local_Powell',metric)
                return r.get('squared_loss',np.inf)
            result=minimize(objective,start,method='Powell',bounds=((.7,8.),(2.,12.)),
                options=dict(maxfev=local_evaluations,ftol=1e-6,xtol=2e-4))
            logs.append(dict(metric=metric,start=start,success=bool(result.success),
                message=str(result.message),nfev=int(result.nfev),objective=float(result.fun)))
    bests={}; table=[]; identities={}
    for mode in ('cubic_5pct','bulk_exact'):
        eligible=[r for r in profiles if r['metric']==mode and r['admissible'] and r['strictly_positive_LJ_resolved']]
        if not eligible: continue
        bests[mode]=min(eligible,key=lambda r:r['squared_loss'])
    for label,fit in [('historical',old),('v10_exchange',prior),*bests.items()]:
        key=fit['scalar_decay'],fit['angular_decay']
        if key not in cache:
            base=observation_matrix(VectorCoefficientBasis(*key),obs)
            cache[key]=(base,*cubic_metric_problem(base,obs))
        base,cm,co=cache[key]; c=np.asarray(fit['coefficients'])
        for matrix,observations,metric in [(base,obs,'legacy_modes'),(cm,co,'independent_cubic')]:
            for row,o in zip(matrix,observations):
                pred=float(row@c)
                table.append(dict(candidate=label,metric=metric,target=o.name,reference=o.target,
                    prediction=pred,units=o.units,role=o.role,scale=o.scale,
                    error=pred-o.target,normalized_residual=(pred-o.target)/o.scale))
        identities[label]=coefficient_identifiability(cm,co,c)
    save_json(out/'observations.json',dict(legacy=[asdict(o) for o in obs],
        independent_cubic=[asdict(o) for o in co],source_sha256=source.reference.sha256,
        candidate_sha256=hashlib.sha256(FIT.read_bytes()).hexdigest(),temperature_K=0,
        zero_traction_shape_prior=False,experimental_strength_in_loss=False,
        legacy_elastic_covariance_in_cubic_coordinates=transformed_mode_covariance([o.scale for o in obs[2:5]])))
    write_csv(out/'residuals.csv',table)
    save_json(out/'identifiability.json',identities)
    save_json(out/'calibration.json',dict(completed=True,best=bests,profiles=profiles,
        deterministic_grid=keys,local_optimizations=logs,elapsed_seconds=time.perf_counter()-started,
        material_accepted=False,production_promoted=False,physical_yield_validated=False,
        notes='Metric sensitivity and bulk/interface compatibility. No global optimum claim.'))
    save_json(out/'progress.json',dict(completed=True,profiles=profiles,elapsed_seconds=time.perf_counter()-started))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--local-evaluations',type=int,default=100)
    args=parser.parse_args(); execute(args.output,args.local_evaluations)


if __name__=='__main__': main()
