"""Actual deterministic radial continuation of the unchanged analytic family."""
from dataclasses import asdict
from functools import lru_cache
import json
import time

import numpy as np
from scipy.optimize import least_squares

from .profiled_range_calibration import fit_profiled_coefficients
from .range_resolved_material import RangeObservationCache,radial_identifiability
from .run_low_stress_cyclic_diagnostic import ROOT,write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .yield_elastic_metric import cubic_metric_problem


OUT=ROOT/'results/fcc111_active_interface/stable_core_v13/material_reprofile'
CHECKPOINT=ROOT/'.cache/stable_core_v13/material_reprofile'


def main():
    if (OUT/'calibration.json').exists():
        raise FileExistsError('preserve completed study; use a new explicitly named run')
    began=time.perf_counter(); source,observations,states=source_and_targets()
    old=json.loads((ROOT/'results/fcc111_active_interface/range_core_v12/material/calibration.json').read_bytes())
    cache=RangeObservationCache(observations); profiles=[]; optimizers=[]
    @lru_cache(maxsize=512)
    def evaluate(point):
        decays=np.exp(point)
        matrix,obs=cubic_metric_problem(cache.matrix(decays),observations)
        fit=fit_profiled_coefficients(matrix,obs)
        row=dict(scalar_decay=decays[0],odd_decay=decays[1],quadrupole_decay=decays[2],**fit)
        row['cubic_GPa']=matrix[2:5]@fit['coefficients']
        hold=[i for i,o in enumerate(obs) if o.role=='heldout']
        row['heldout_normalized_rms']=float(np.sqrt(np.mean(fit['residuals'][hold]**2)))
        profiles.append(row)
        if len(profiles)==1 or len(profiles)%10==0:
            save_json(OUT/'progress.json',dict(completed=False,profiles=len(profiles),
                best_loss=min(r['squared_loss'] for r in profiles),elapsed_seconds=time.perf_counter()-began))
            print(f"profile {len(profiles)}: loss={row['squared_loss']:.8g}, ranges={decays}",flush=True)
        return row
    # Fixed declared starts: the two actual v12 candidates and a distinct
    # already evaluated interior range. No random search and no new target.
    starts=[[old['best'][key][k] for k in ('scalar_decay','odd_decay','quadrupole_decay')]
            for key in ('cubic_5pct','bulk_exact')]
    starts.append([2.7,6.4,10.])
    for index,start in enumerate(starts):
        def residual(logs):
            fit=evaluate(tuple(logs)); return fit['residuals'][fit['selected_rows']]
        run=least_squares(residual,np.log(start),jac='3-point',diff_step=2e-4,
            bounds=np.log([[.7,2.,2.],[8.,12.,12.]]),max_nfev=65,
            ftol=2e-9,xtol=2e-9,gtol=2e-7)
        optimizers.append(dict(start=start,success=bool(run.success),message=str(run.message),
            nfev=run.nfev,njev=run.njev,optimality=float(run.optimality),
            final_ranges=np.exp(run.x),final_loss=float(run.fun@run.fun)))
        save_json(CHECKPOINT/f'after_start_{index}.json',dict(optimizers=optimizers,profiles=profiles))
    eligible=[r for r in profiles if r['strictly_positive_LJ_resolved']]
    best=min(eligible,key=lambda row:row['squared_loss'])
    decays=[best[k] for k in ('scalar_decay','odd_decay','quadrupole_decay')]
    matrix,obs=cubic_metric_problem(cache.matrix(decays),observations)
    rows=[dict(observable=o.name,role=o.role,units=o.units,target=o.target,prediction=value,
        scale=o.scale,normalized_residual=(value-o.target)/o.scale)
        for value,o in zip(matrix@best['coefficients'],obs)]
    write_csv(OUT/'residuals.csv',rows)
    modes=np.array([[1/3,2/3,0.],[.5,-.5,0.],[0.,0.,1.]])
    pieces=modes@matrix[2:5]*best['coefficients']
    write_csv(OUT/'elastic_channels.csv',[dict(channel=channel,target=(modes@np.array([114.,62.,32.]))[i],
        prediction=float(sum(pieces[i])),**dict(zip(best['coefficient_order'],pieces[i])))
        for i,channel in enumerate(('bulk_modulus','tetragonal_Cprime','C44'))])
    save_json(OUT/'identifiability.json',radial_identifiability(decays,best['coefficients'],observations))
    save_json(OUT/'calibration.json',dict(completed=True,best=best,profiles=profiles,
        starts=starts,optimizers=optimizers,elapsed_seconds=time.perf_counter()-began,
        source_sha256=source.reference.sha256,observations=[asdict(o) for o in obs],states=states,
        energy_family_changed=False,targets_weights_changed=False,material_accepted=False,
        physical_yield_validated=False,kinetics_calibrated=False))
    save_json(OUT/'progress.json',dict(completed=True,profiles=len(profiles),best_loss=best['squared_loss']))
    print(f"completed: best loss={best['squared_loss']}, C={best['cubic_GPa']}",flush=True)


if __name__=='__main__':
    main()
