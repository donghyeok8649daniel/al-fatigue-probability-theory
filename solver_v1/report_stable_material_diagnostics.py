"""Independent replay of completed v13 fits: admissibility before adoption."""
import hashlib
import csv
import json

import numpy as np

from .range_resolved_material import RangeObservationCache,build_range_surface,RangeBulkValidation,radial_identifiability
from .run_spectral_material_profiles import wavepoints
from .run_low_stress_cyclic_diagnostic import ROOT,write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .vector_registry_audit import stationary_state
from .vector_material_calibration import UNITS,IDEAL_H
from .yield_elastic_metric import cubic_metric_problem


OUT=ROOT/'results/fcc111_active_interface/stable_core_v13'


def validate_targeted_tails():
    """Resolve the marginal sign at the actually observed worst wavevector.

    One resolved negative mode suffices for rejection; positive values at one
    wavevector can never certify the rest of the Brillouin zone.
    """
    spectral=json.loads((OUT/'material_spectral/profiles.json').read_bytes())['histories']['v13_other_basin']
    stable=json.loads((OUT/'material_stable_ranges/calibration.json').read_bytes())['best']
    fits=[('v13_fixed_range_spectral',spectral['decays'],spectral['final_fit']['coefficients']),
          ('v13_reoptimized_necessary_spectral',stable['decays'],stable['coefficients'])]
    with (OUT/'material_finite_q_refinement.csv').open(newline='') as stream:
        previous_rows=list(csv.DictReader(stream))
    rows=[]; summaries=[]
    for name,decays,c in fits:
        candidates=[r for r in previous_rows if r['candidate']==name and float(r['radius_over_L0'])==10.]
        worst=min(candidates,key=lambda r:float(r['minimum_H_eV_L0sq']))
        path=worst['path']; fraction=float(worst['fraction'])
        q=next(q for p,f,q in wavepoints() if p==path and f==fraction)
        model=build_range_surface(*decays,c); previous=None
        for radius in (10.,12.,16.,20.):
            check=RangeBulkValidation(model,cutoff=radius)
            value=check.evaluate(check.crystallographic_wavevector(q))
            change=float(np.linalg.norm(value['matrix']-previous,2)) if previous is not None else None
            rows.append(dict(candidate=name,path=path,fraction=fraction,radius_over_L0=radius,
                minimum_H_eV_L0sq=value['eigenvalues'][0],matrix_radius_change=change))
            previous=value['matrix']
        last=rows[-1]
        summaries.append(dict(candidate=name,path=path,fraction=fraction,
            minimum_radius20_H=last['minimum_H_eV_L0sq'],radius16_to20_matrix_change=last['matrix_radius_change'],
            negative_sign_above_last_radius_change=bool(last['minimum_H_eV_L0sq'] < -last['matrix_radius_change']),
            radius_change_is_not_rigorous_tail_bound=True,material_accepted=False))
        print(f"large-radius sign {name}: minH={last['minimum_H_eV_L0sq']:.7g}, change={last['matrix_radius_change']:.4g}",flush=True)
    write_csv(OUT/'material_spectral_tail_refinement.csv',rows)
    write_csv(OUT/'material_spectral_tail_summary.csv',summaries)
    with (OUT/'material_comparison.csv').open(newline='') as stream:
        comparisons=list(csv.DictReader(stream))
    tails={r['candidate']:r for r in summaries}
    for row in comparisons:
        initial=row.pop('spectral_classification',row.get('radius10_spectral_classification'))
        row['radius10_spectral_classification']=initial
        extra=tails.get(row['candidate'])
        row['final_spectral_classification']=(
            'negative_targeted_mode_after_radius_refinement' if extra and extra['negative_sign_above_last_radius_change']
            else 'targeted_tail_sign_unresolved' if extra else initial)
        row['final_checked_radius_over_L0']=20. if extra else 10.
        row['final_minimum_H_eV_L0sq']=extra['minimum_radius20_H'] if extra else row['minimum_radius10_H']
    write_csv(OUT/'material_comparison.csv',comparisons)
    decision=json.loads((OUT/'material_decision.json').read_bytes())
    decision['targeted_tail_classification']=summaries
    save_json(OUT/'material_decision.json',decision)


def main():
    source,observations,states=source_and_targets(); cache=RangeObservationCache(observations)
    paths=[OUT.parent/'range_core_v12/material/calibration.json',
           OUT/'material_reprofile/calibration.json',OUT/'material_spectral/profiles.json',
           OUT/'material_stable_ranges/calibration.json']
    old,new,spectral,stable=[json.loads(p.read_bytes()) for p in paths]
    base=old['best']['cubic_5pct']; unconstrained=new['best']
    fixed=spectral['histories']['v13_other_basin']
    fits=[('v12_reference_ranges',[base[k] for k in ('scalar_decay','odd_decay','quadrupole_decay')],base['coefficients']),
          ('v13_unconstrained_ranges',[unconstrained[k] for k in ('scalar_decay','odd_decay','quadrupole_decay')],unconstrained['coefficients']),
          ('v13_fixed_range_spectral',fixed['decays'],fixed['final_fit']['coefficients']),
          ('v13_reoptimized_necessary_spectral',stable['best']['decays'],stable['best']['coefficients'])]
    summaries=[]; residuals=[]; waves=[]; roots=[]; channels=[]
    for name,decays,coefficients in fits:
        c=np.asarray(coefficients); matrix,obs=cubic_metric_problem(cache.matrix(decays),observations)
        prediction=matrix@c; r=(prediction-np.array([o.target for o in obs]))/np.array([o.scale for o in obs])
        fit=[i for i,o in enumerate(obs) if o.role=='fit']; hold=[i for i,o in enumerate(obs) if o.role=='heldout']
        model=build_range_surface(*decays,c); previous={}; minima=[]; changes=[]
        for radius in (6.,8.,10.):
            check=RangeBulkValidation(model,cutoff=radius); current={}; local=[]
            for path,fraction,q in wavepoints():
                value=check.evaluate(check.crystallographic_wavevector(q)); local.append(value['eigenvalues'][0])
                key=path,fraction; current[key]=value['matrix']
                change=float(np.linalg.norm(current[key]-previous[key],2)) if key in previous else None
                waves.append(dict(candidate=name,radius_over_L0=radius,path=path,fraction=fraction,
                    minimum_H_eV_L0sq=value['eigenvalues'][0],matrix_radius_change=change))
                if radius==10.: changes.append(change)
            minima.append(min(local)); previous=current
        floor=max(changes)
        classification=('resolved_sampled_instability' if minima[-1]<-floor else
            'positive_sampled_margin_not_global_certificate' if minima[-1]>floor else 'marginal_below_tail_resolution')
        C=prediction[2:5]
        summaries.append(dict(candidate=name,scalar_decay=decays[0],odd_decay=decays[1],quadrupole_decay=decays[2],
            fit_loss=float(r[fit]@r[fit]),heldout_normalized_rms=float(np.sqrt(np.mean(r[hold]**2))),
            maximum_exact_residual=float(np.max(abs(r[:2]))),C11_GPa=C[0],C12_GPa=C[1],C44_GPa=C[2],
            B_GPa=(C[0]+2*C[1])/3,Cprime_GPa=(C[0]-C[1])/2,
            minimum_radius6_H=minima[0],minimum_radius8_H=minima[1],minimum_radius10_H=minima[2],
            maximum_radius8_to10_matrix_change=floor,spectral_classification=classification,
            material_accepted=False,experimental_yield_validated=False))
        modes=np.array([[1/3,2/3,0],[.5,-.5,0],[0,0,1.]])
        pieces=(modes@matrix[2:5])*c
        for channel,terms in zip(('B','Cprime','C44'),pieces):
            channels.append(dict(candidate=name,channel=channel,units='GPa',prediction=float(sum(terms)),
                **dict(zip(('u','v','A','B','C','D3','D1','D2'),terms))))
        for o,p,res in zip(obs,prediction,r):
            residuals.append(dict(candidate=name,observable=o.name,role=o.role,units=o.units,
                target=o.target,prediction=p,scale=o.scale,normalized_residual=res))
        for label,guess,index in [('perfect',[IDEAL_H,0.,0.],0),('fault',states['fault'],0),('saddle',states['saddle'],1)]:
            result=stationary_state(model,guess,expected_index=index)
            roots.append(dict(candidate=name,state=label,valid=result['valid'],
                energy_J_m2=float(UNITS.energy_to_surface(result['evaluation'].energy)),
                a=result['q'][0],x=result['q'][1],y=result['q'][2],
                force_residual=result['force_residual'],eigenvalues=result['eigenvalues'].tolist()))
        print(f"replayed material {name}: loss={summaries[-1]['fit_loss']:.7g}, {classification}",flush=True)
    write_csv(OUT/'material_comparison.csv',summaries); write_csv(OUT/'material_residuals.csv',residuals)
    write_csv(OUT/'material_finite_q_refinement.csv',waves); write_csv(OUT/'material_stationary_states.csv',roots)
    write_csv(OUT/'material_elastic_channels.csv',channels)
    save_json(OUT/'material_stable_ranges/identifiability.json',radial_identifiability(
        stable['best']['decays'],stable['best']['coefficients'],observations))
    save_json(OUT/'material_decision.json',dict(completed=True,
        inputs=[dict(path=p.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
        source_sha256=source.reference.sha256,material_accepted=False,physical_yield_validated=False,
        whole_zone_stability_certified=False,new_energy_terms_added=False,physical_seconds=False,physical_Hz=False,
        conclusion='Fit loss alone is not adoption; finite-q, held-out, and core/source gates remain separate.'))
    validate_targeted_tails()


if __name__=='__main__':
    main()
