"""Persist independent v8 numerical checks; no new DFT or atomistic energy.

Run with python -m results.silicon_wafer_feasibility.audit_doping_v8 --output DIR.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
from scipy.optimize import brentq

from solver_v1.silicon_doping_research import (
    directional_young_GPa, equilibrated_charge_branches, load_elastic_samples,
)
from solver_v1.silicon_doping_validation import (
    fourier_surface_compliance, two_orbital_canonical, two_orbital_product,
)
from solver_v1.silicon_specimen_research import AnisotropicModeI

ROOT = Path(__file__).resolve().parents[2]


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(output):
    output = Path(output)
    if (output/'independent_validation.json').exists():
        raise FileExistsError('preserve existing audit; use a new output')
    output.mkdir(parents=True, exist_ok=True)
    elastic, electronic, support = [], [], []
    for sample in load_elastic_samples().values():
        support.append(dict(sample_id=sample.sample_id, min_C=sample.temperature_min_C,
                            max_C=sample.temperature_max_C))
        scenarios = [(t, np.zeros(3), 'measured_fit') for t in
                     [sample.temperature_min_C, 25., sample.temperature_max_C]]
        scenarios += [(25., np.array(sign)*[.3, .1, .2], 'c0_sensitivity_corner')
                      for sign in itertools.product([-1, 1], repeat=3)]
        for t, dc, label in scenarios:
            c = sample.constants_GPa(t)+dc
            independent = fourier_surface_compliance(*c)
            airy = AnisotropicModeI(*c).H_GPa_inv
            error = abs(independent['H_GPa_inv']/airy-1)
            assert error < 3e-12
            elastic.append(dict(sample_id=sample.sample_id, temperature_C=t, scenario=label,
                                delta_c11_GPa=dc[0], delta_c12_GPa=dc[1], delta_c44_GPa=dc[2],
                                H_Airy_GPa_inv=airy, H_Fourier_GPa_inv=independent['H_GPa_inv'],
                                relative_error=error, subspace_residual=independent['subspace_residual']))
    for q, mu, kt in itertools.product([[-.1, .05], [0., 0.], [.12, -.08]], [-.1, 0., .1], [.01, .025, .08]):
        sectors = two_orbital_canonical(q, kt)
        r = equilibrated_charge_branches(*sectors, excess_electron_counts=[0, 1, 2], mu_eV=mu, kBT_eV=kt)
        exact = two_orbital_product(q, mu, kt)
        errors = {key: float(np.max(abs(np.asarray(r[key])-exact[key])))
                  for key in ['grand_potential_eV', 'gradient', 'hessian', 'mean_excess_electrons']}
        for key in errors:
            assert_allclose(r[key], exact[key], rtol=1e-12, atol=1e-14)
        electronic.append(dict(q0=q[0], q1=q[1], mu_eV=mu, kBT_eV=kt, **errors))
    q, mu, kt = np.array([.03, -.02]), .01, .025
    exact = two_orbital_product(q, mu, kt)
    dn = exact['dmean_dq']
    fixed_mean_h = exact['hessian']+np.outer(dn, dn)/exact['dmean_dmu']
    derivatives = []
    for step in [1e-3, 5e-4, 2.5e-4, 1e-4, 5e-5, 1e-5]:
        columns = []
        for axis in np.eye(2):
            gradients = []
            for sign in [1, -1]:
                point = q+sign*step*axis
                chosen = brentq(lambda m: two_orbital_product(point, m, kt)['mean_excess_electrons']
                               -exact['mean_excess_electrons'], -.8, .8, xtol=1e-14)
                gradients.append(two_orbital_product(point, chosen, kt)['gradient'])
            columns.append((gradients[0]-gradients[1])/(2*step))
        derivatives.append(dict(step=step, max_hessian_error=float(np.max(abs(np.array(columns).T-fixed_mean_h)))))
    assert derivatives[-1]['max_hessian_error'] < 1e-8
    # The published E's are rounded to 0.1 GPa. Check their rounding boxes
    # intersect the exact cubic identity, without treating it as a Cij fit.
    with (ROOT/'results/silicon_doping_v7/sources/noda_2023_ideal_strength.csv').open() as stream:
        source = list(csv.DictReader(stream))
    rounded = []
    for row in source:
        a, b, c = [float(row[f'E{d}_LDA_GPa']) for d in ['100', '110', '111']]
        lower = 4/(b+.05)-1/(a-.05)-3/(c-.05)
        upper = 4/(b-.05)-1/(a+.05)-3/(c+.05)
        assert lower <= 0 <= upper
        rounded.append(dict(carrier_type=row['carrier_type'], carrier_cm3=row['excess_carrier_cm3'],
                            residual_GPa_inv=4/b-1/a-3/c, rounding_lower=lower, rounding_upper=upper))
    # Nonidentifiability witness: S11 and Q=S11-S12-S44/2 fixed; vary S12.
    witness = []
    base = np.array([[165.5, 63.7, 63.7], [63.7, 165.5, 63.7], [63.7, 63.7, 165.5]])
    s11, s12 = np.linalg.inv(base)[0, :2]
    invariant = s11-s12-1/(2*79.6)
    for delta in [-.0004, 0., .0004]:
        off = s12+delta
        compliance = np.full((3, 3), off)
        np.fill_diagonal(compliance, s11)
        stiffness = np.linalg.inv(compliance)
        cs = (stiffness[0, 0], stiffness[0, 1], 1/(2*(s11-off-invariant)))
        witness.append(dict(c11_GPa=cs[0], c12_GPa=cs[1], c44_GPa=cs[2],
                            E100_GPa=directional_young_GPa(*cs, [1, 0, 0]),
                            E110_GPa=directional_young_GPa(*cs, [1, 1, 0]),
                            E111_GPa=directional_young_GPa(*cs, [1, 1, 1]),
                            H_GPa_inv=fourier_surface_compliance(*cs)['H_GPa_inv']))
    for key in ['E100_GPa', 'E110_GPa', 'E111_GPa']:
        assert_allclose([r[key] for r in witness], witness[1][key], rtol=1e-13)
    assert np.ptp([r['H_GPa_inv'] for r in witness]) > 1e-5
    files = ['solver_v1/silicon_doping_research.py', 'solver_v1/silicon_doping_validation.py',
             'solver_v1/silicon_specimen_research.py', 'solver_v1/test_silicon_doping_validation.py',
             'results/silicon_wafer_feasibility/audit_doping_v8.py',
             'results/silicon_doping_v8/sources/jaakkola_2014_elastic.csv',
             'results/silicon_doping_v7/sources/noda_2023_ideal_strength.csv']
    summary = dict(status='independent mathematical checks passed; actual doped fracture uncalibrated',
                   corrected_source_support=support,
                   elastic_cases=len(elastic), max_relative_H_error=max(r['relative_error'] for r in elastic),
                   electronic_cases=len(electronic), electronic_max_absolute_errors={key:max(r[key] for r in electronic) for key in errors},
                   fixed_mean_charge=dict(gradient_constraint='mean N, not exact canonical N',
                                          fixed_mu_hessian=exact['hessian'].tolist(),
                                          fixed_mean_hessian=fixed_mean_h.tolist(),
                                          refinement=derivatives),
                   source_young_identity_checks=rounded, nonidentifiability_witness=witness,
                   truncated_sector_counterexample=dict(complete_mean_N=1., truncated_mean_N=2/3,
                       grand_potential_bias_eV=.025*np.log(4/3), normalization_is_not_completeness=True),
                   new_DFT=0, new_MD=0, new_atomistic_energy=0, physical_clock_available=False,
                   source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    write_csv(output/'independent_elasticity.csv', elastic)
    write_csv(output/'independent_electronic.csv', electronic)
    write_csv(output/'young_nonidentifiability.csv', witness)
    (output/'independent_validation.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ['elastic_cases', 'max_relative_H_error', 'electronic_cases', 'electronic_max_absolute_errors']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output)
