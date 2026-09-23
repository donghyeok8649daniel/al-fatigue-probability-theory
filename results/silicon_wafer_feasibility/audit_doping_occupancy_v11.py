"""Exact independent-site occupancy arithmetic for the actual nominal prism.

This is a declared statistical example, not a measured dopant placement law,
activation model, specimen crack probability or material calibration.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import binom


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    with np.load(args.geometry) as d:
        sites=len(d['reference']);volume=float(np.prod(d['nominal_lengths_A']))*1e-24
        free=int(d['free'].sum());area=float(d['area_A2'])
    density=sites/volume;rows=[];max_formula_error=0.
    for concentration in [1e15,1e16,1e17,1e18,1e19,1e20,1e21]:
        p=concentration/density
        if not 0<=p<=1:raise ValueError('requested concentration exceeds site density')
        at_least_one=float(-np.expm1(sites*np.log1p(-p)))
        one=float(sites*p*np.exp((sites-1)*np.log1p(-p)))
        scipy_at_least_one=float(binom.sf(0,sites,p))
        max_formula_error=max(max_formula_error,abs(at_least_one-scipy_at_least_one),abs(one-float(binom.pmf(1,sites,p))))
        rows.append(dict(dopant_concentration_cm3=concentration,independent_site_probability=p,
            expected_dopants=concentration*volume,probability_no_dopants=float(np.exp(sites*np.log1p(-p))),
            probability_exactly_one=one,probability_at_least_one=at_least_one,
            probability_two_or_more=float(binom.sf(1,sites,p))))
    if max_formula_error>1e-13:raise ValueError('independent binomial formula mismatch')
    with (args.output/'independent_site_occupancy.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    result=dict(geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),
        lattice_sites=sites,free_atoms=free,nominal_volume_cm3=volume,nominal_area_A2=area,
        nominal_site_density_cm3=density,one_dopant_per_prism_concentration_cm3=1/volume,
        maximum_binomial_formula_difference=max_formula_error,
        assumptions='uniform independent substitutional sites over nominal prism including grips; no segregation, charge or correlations',
        sample_distribution_measured=False,stress_active_site_ensemble_validated=False,
        crack_initiation_probability=None,new_potential_calls=0,new_DFT=0,new_MD=0)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--geometry',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
