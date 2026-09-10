"""Audit SOURCE interpolation/cutoff normal jets; never alter calibration data."""
import argparse
from pathlib import Path

import numpy as np

from .interface_normal_development_targets import normal_development_observations
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .vector_material_calibration import UNITS


def crossing_gaps(reference, *, low_ratio=1., high_ratio=4.):
    """Exact geometric loss of a cross-interface neighbor at the tabulated edge.

    This finite cutoff belongs to the external EAM SOURCE. The canonical
    LJ/Bessel model is not truncated, altered, or recalibrated by this audit.
    """
    cutoff = reference.r[-1]
    result = []
    for k in range(1, int(np.ceil(cutoff/reference.h))+2):
        xy = reference.R+reference.geometry.abc_shift(k)
        radii_squared = np.sum(xy*xy, axis=1)
        valid = radii_squared < cutoff**2
        gaps = np.sqrt(cutoff**2-radii_squared[valid])-(k-1)*reference.h
        ratios = gaps/reference.h
        for value in np.unique(np.round(ratios, 12)):
            if low_ratio <= value <= high_ratio:
                count = int(np.sum(abs(ratios-value) < 2e-12))
                result.append(dict(layer=k, a_over_h=float(value), inplane_degeneracy=count,
                    cross_pair_multiplicity=k*count))
    return sorted(result, key=lambda r: (r['a_over_h'], r['layer']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('preserve source diagnostic')
    source, previous, states = source_and_targets()
    obs, _ = normal_development_observations(source, previous, states)
    ref = source.reference; crossings = crossing_gaps(ref)
    rows, neighborhoods = [], []
    selected = list(dict.fromkeys(o.state for o in obs if o.state is not None
                                  and o.state[1] == 0. and o.state[2] == 0.))
    for state in selected:
        q = np.asarray(state); value = source.evaluate(q)
        ratio = q[0]/source.h
        distance = min(abs(ratio-r['a_over_h']) for r in crossings)
        for step in (4e-5, 2e-5, 1e-5):
            plus, minus = source.evaluate(q+[step, 0., 0.]), source.evaluate(q-[step, 0., 0.])
            rows.append(dict(a_over_h=ratio, step_over_L0=step,
                source_traction_MPa=float(UNITS.force_to_traction_mpa(value.gradient[0])),
                source_Haa=value.hessian[0, 0],
                gradient_FD_error=(plus.energy-minus.energy)/(2*step)-value.gradient[0],
                Haa_FD_error=(plus.gradient[0]-minus.gradient[0])/(2*step)-value.hessian[0, 0],
                nearest_cutoff_crossing_distance_over_h=distance,
                stencil_avoids_cutoff_crossing=bool(step/source.h < distance)))
    for ratio in sorted({r['a_over_h'] for r in crossings}):
        for epsilon in (4e-5, 2e-5, 1e-5):
            left = source.evaluate(((ratio-epsilon)*source.h, 0., 0.))
            right = source.evaluate(((ratio+epsilon)*source.h, 0., 0.))
            neighborhoods.append(dict(a_over_h=ratio, epsilon_over_h=epsilon,
                force_left=left.gradient[0], force_right=right.gradient[0],
                Haa_left=left.hessian[0, 0], Haa_right=right.hessian[0, 0],
                energy_difference=right.energy-left.energy,
                Haa_difference=right.hessian[0, 0]-left.hessian[0, 0]))
    r = ref.r[-1]; z = float(ref._rphi(r)); zp = float(ref._rphi(r, 1)); zpp = float(ref._rphi(r, 2))
    save_json(args.out/'scope.json', dict(completed=True, source_sha256=ref.sha256,
        nominal_cutoff_angstrom=ref.cutoff, actual_last_tabulated_radius_angstrom=r,
        radial_table_spacing_angstrom=float(ref.r[1]-ref.r[0]),
        last_pair_eV=z/r, last_pair_first_eV_angstrom=zp/r-z/r**2,
        last_pair_second_eV_angstrom_sq=zpp/r-2*zp/r**2+2*z/r**3,
        last_density=[float(ref._rho(r, i)) for i in range(3)],
        outside_source_cutoff_jets=[0., 0., 0.],
        source_modified=False, fitting_observations_modified=False,
        candidate_error_explained_by_cutoff=False,
        interpretation='quantitative source smoothness audit; do not excuse candidate errors without comparison'))
    write_csv(args.out/'crossing_geometry.csv', crossings)
    write_csv(args.out/'source_derivative_refinement.csv', rows)
    write_csv(args.out/'cutoff_neighborhoods.csv', neighborhoods)
    print('source normal/cutoff diagnostic completed', len(crossings), 'layer-shell crossings', flush=True)


if __name__ == '__main__':
    main()
