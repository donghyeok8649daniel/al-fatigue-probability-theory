"""Closed-form replay of the synthetic charge-dependent first-entry example."""
import argparse,csv,hashlib,json
from pathlib import Path


def main(args):
    output=args.result/'closed_form_replay.json'
    if output.exists():raise ValueError('do not overwrite an existing replay')
    rows=list(csv.DictReader((args.result/'fast_charge_hard_basin.csv').open(encoding='utf-8')))
    maximum=0.;cases={}
    for row in rows:
        p=float(row['occupation']);nu=float(row['charge_rate']);r=float(row['spatial_rate'])
        exact=3/r+3*(1-p)*(nu*nu+3*nu*r+r*r)/(nu*p*(nu+r)*(nu+3*r))
        difference=abs(exact-float(row['MFPT_model']));maximum=max(maximum,difference)
        cases[(p,nu,r)]=dict(occupation=p,charge_rate=nu,spatial_rate=r,exact_MFPT=exact,error=difference)
    if maximum>1e-7:raise ValueError('closed-form MFPT does not reproduce the stored linear-system result')
    result=dict(complete=True,cases=len(cases),maximum_absolute_error_model_time=maximum,
        formula='3/r + 3*(1-p)*(nu^2+3*nu*r+r^2)/(nu*p*(nu+r)*(nu+3*r))',
        derivation='solve five backward first-entry equations with T(q2,sector1)=0 and initial q0 charge weights (1-p,p)',
        large_nu_excess='T-3/r = 3*(1-p)/(p*nu) + O(nu^-2)',
        small_nu_leading='T = (1-p)/(p*nu) + O(1), for fixed 0<p<1',
        p_zero='sector0 becomes a closed nonabsorbing class; no finite mean first-entry time',
        occupation_is_not_dopant_atomic_fraction=True,physical_clock=None,Si_initiation_calibrated=False,
        source_csv_sha256=hashlib.sha256((args.result/'fast_charge_hard_basin.csv').read_bytes()).hexdigest(),
        values=list(cases.values()))
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='values'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);main(p.parse_args())
