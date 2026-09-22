"""Chemical preparation coverage, not a substitute for carrier/rupture data."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_concentration_research import independent_site_configuration_coverage
from results.silicon_wafer_feasibility.run_concentration_cleavage_v10 import dump,csv_write


def main(args):
    protocol=json.loads(args.protocol.read_text())
    n=protocol['atoms'];volume=protocol['reference_material_volume_A3'];site_density=n*1e24/volume
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for species in ['B','P']:
        patterns=[p['substitution_indices'] for p in protocol['plans'] if p['species'] in ['Si',species]]
        for c in sorted(set([0.,1e15,1e16,1e17,1e18,1e19,1e20,1e21]+[i*1e24/volume for i in [1,2,4,8]])):
            x=c/site_density;r=independent_site_configuration_coverage(n,x,patterns)
            log_none=n*np.log1p(-x)
            rows.append(dict(species=species,chemical_concentration_cm3=c,site_fraction=x,
                expected_initial_dopant_count=n*x,probability_no_local_chemical_dopant=float(np.exp(log_none)),
                selected_initial_patterns=len(r['unique_initial_patterns']),duplicates_removed=r['duplicates_removed'],
                literal_selected_mass=r['covered_mass'],unrepresented_mass=r['missing_mass'],
                host_sites_for_one_periodic_dopant=None if c==0 else site_density/c,
                independent_site_assumption=True,symmetry_augmented=False,
                physical_failure_probability=None))
    csv_write(args.output/'chemical_coverage.csv',rows)
    dump(args.output/'summary.json',dict(rows=len(rows),site_count=n,material_volume_A3=volume,site_density_cm3=site_density,
        scope='bookkeeping for literal initial patterns under an explicitly hypothetical IID substitution law',
        coverage_is_not_accuracy_or_failure_probability=True,
        local_zero_dopant_does_not_mean_unchanged_carrier_environment=True,
        no_symmetry_augmentation=True,no_renormalization=True,monte_carlo_samples=0))
    print(json.dumps(rows[-4:]),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--protocol',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
