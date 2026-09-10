"""Post-fit analytic tail/source-support audit; no optimization or cutoff fit."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from scipy.special import zeta

from .validate_coordination_screening import load_candidate
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import UNITS, LENGTH_M
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh tail-audit output required')
    source, _, _ = source_and_targets()
    support = float(source.reference.r[-1])
    rows, errors, relative_errors, terms, bindings = [], [], [], [], []
    for path in args.models:
        model, shape, c, law, definition = load_candidate(path)
        if definition['source_sha256'] != source.reference.sha256:
            raise ValueError('source binding mismatch')
        bindings.append(dict(model=path.parent.name+'/'+path.name, shape=shape,
            coefficients=c, screening_law=law, calibration_sha256=hashlib.sha256(
                (path/'calibration.json').read_bytes()).hexdigest()))
        asymptotic_errors = []
        for ratio in (2.5, 3., 4., 6., 10., 20., 40., 80.):
            a = ratio*model.h
            value = model.evaluate((a, 0., 0.))
            target = source.evaluate((a, 0., 0.))
            pair = value.components['pair_3'][1]+value.components['pair_6'][1]
            attractive = value.components['pair_3'][1]
            np.testing.assert_allclose(sum(jet[1] for jet in value.components.values()),
                                       value.gradient[0], rtol=1e-13, atol=1e-16)
            for component, jet in value.components.items():
                terms.append(dict(model=path.parent.name+'/'+path.name, a_over_h=ratio,
                    component=component, force_MPa=float(UNITS.force_to_traction_mpa(jet[1])),
                    Haa_eV_L0sq=jet[4], energy_eV_cell=jet[0]))
            area = np.sqrt(3)/2
            mean = 2*np.pi*c[1]/(area*model.h**5)*(zeta(4, ratio)-(ratio-1)*zeta(5, ratio))
            leading = np.pi*c[1]/(6*area*model.h**2*a**3)
            if ratio >= 10:
                # Nonzero reciprocal terms are exponentially negligible here.
                np.testing.assert_allclose(attractive, mean, rtol=2e-10, atol=1e-16)
                asymptotic_errors.append(abs(attractive/leading-1))
                errors.append(abs(attractive-mean))
                relative_errors.append(abs(attractive/mean-1))
            if a*LENGTH_M/1e-10 > support:
                np.testing.assert_array_equal(target.gradient, np.zeros(3))
                np.testing.assert_array_equal(target.hessian, np.zeros((3, 3)))
            rows.append(dict(model=path.parent.name+'/'+path.name, a_over_h=ratio,
                a_angstrom=a*LENGTH_M/1e-10, source_force_MPa=float(UNITS.force_to_traction_mpa(target.gradient[0])),
                total_force_MPa=float(UNITS.force_to_traction_mpa(value.gradient[0])),
                pair_force_MPa=float(UNITS.force_to_traction_mpa(pair)),
                environment_force_MPa=float(UNITS.force_to_traction_mpa(value.gradient[0]-pair)),
                attractive_pair_eV_L0=attractive, exact_reciprocal_zero_eV_L0=mean,
                leading_asymptotic_eV_L0=leading, leading_relative_change=attractive/leading-1,
                target_support_angstrom=support, no_parameter_refit=True))
        if not np.all(np.diff(asymptotic_errors) < 0):
            raise ArithmeticError('large-opening attractive tail did not converge to the derived asymptote')
    write_csv(args.out/'tail_comparison.csv', rows)
    write_csv(args.out/'per_term_tail.csv', terms)
    save_json(args.out/'completion.json', dict(completed=True,
        max_large_a_pair_mean_error=max(errors), max_large_a_pair_mean_relative_error=max(relative_errors),
        source_declared_cutoff_angstrom=source.reference.cutoff,
        source_effective_support_angstrom=support, source_sha256=source.reference.sha256,
        model_bindings=bindings,
        candidate_cutoff_introduced=False, source_targets_changed=False,
        tail_difference_is_not_whole_family_impossibility=True, material_accepted=False))


if __name__ == '__main__':
    main()
