"""Lock a positive-boundary static candidate using training and declared gates."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np


RUNS=['angular_boundary_fit_results','refined_boundary_fit_results','local_boundary_fit_results','guard_boundary_fit_results']


def main(args):
    output=args.cache/'final_selection.json'
    if output.exists():raise ValueError('final selection already locked')
    records=[];eligible=[]
    for name in RUNS:
        folder=args.cache/name;summary=json.loads((folder/'summary.json').read_text())
        if not summary['complete']:raise RuntimeError('shape run incomplete')
        profiles=json.loads((folder/'profiles.json').read_text())
        for r in profiles:
            valid=bool(r.get('bulk_guard_passed') and r.get('training_evaluated'))
            if valid:
                a=r['anchor'];target=np.array([153.28991,56.25009,72.17693])
                actual=np.array([a['C11_GPa'],a['C12_GPa'],a['C44_GPa']])
                valid=bool(np.max(abs(actual-target))<1e-7 and np.max(abs(a['gradient']))<1e-9 and
                    min(a['internal_eigenvalues'])>0 and min(a['elastic_eigenvalues'])>0 and
                    r['sampled_minimum_eigenvalue_eV_A2']>=-r['numerical_stability_floor_eV_A2'])
            records.append(dict(run=name,profile=r['profile'],eligible=valid,training_plus_EOS_objective_squared=r.get('objective_squared'),
                bulk_guard_passed=r.get('bulk_guard_passed'),excluded_metrics_used=False))
            if valid:eligible.append((float(r['objective_squared']),name,r['profile']))
    if not eligible:raise RuntimeError('no eligible candidate')
    score,name,profile=min(eligible)
    selection=args.cache/name/'selected_fit.json';chosen=json.loads(selection.read_text())
    if chosen['profile']!=profile or chosen['objective_squared']!=score:raise RuntimeError('per-run selection does not match global training selection')
    report=dict(selected_run=name,selected_profile=profile,selected_training_plus_EOS_objective_squared=score,
        selected_fit_sha256=hashlib.sha256(selection.read_bytes()).hexdigest(),
        selection_rule='minimum same declared training+EOS objective among positive-boundary roots passing exact a0/C11/C12/C44, positive local curvature, sampled-q numerical gate and original bulk energy guard',
        examined_profiles=len(records),eligible_profiles=len(eligible),profiles=records,
        global_shape_optimality_certified=False,excluded_metrics_used=False,material_approved=False,
        first_initiation_validated=False,physical_clock_validated=False,
        selector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output.write_bytes((json.dumps(report,indent=2,allow_nan=False)+'\n').encode())
    print(json.dumps({k:v for k,v in report.items() if k!='profiles'}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True);main(p.parse_args())
