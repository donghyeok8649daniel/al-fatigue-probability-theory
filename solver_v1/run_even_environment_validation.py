"""Independent v16/v17 interface/curvature checks, not an optimizer or PDE run."""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from .isotropic_bulk_validation import IsotropicBulkBasis
from .quadrupole_saturation import SaturatedQuadrupoleInterface
from .quartic_angular_material import QuarticSymmetryInterface
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .validate_tail_calibration import load_material
from .vector_material_calibration import UNITS, LENGTH_M


def resolved_opening_extrema(model, samples):
    """Refine near h as well as the remote branch; no unimodal assumption.

    Large alpha can create narrow extrema near the perfect interface that
    TWO uniform coarse grids both miss. The geometric part is numerical
    resolution, not a new material length or an empirical stress threshold.
    """
    if int(samples) != samples or samples < 16:
        raise ValueError('at least 16 declared bracket samples required')
    grid = model.h*np.unique(np.r_[np.linspace(1., 5., samples),
                                   1+np.geomspace(1e-8, .1, samples)])
    at = lambda a: model.evaluate([float(a), 0., 0.])
    values = [at(a) for a in grid]
    rows = []
    for a, b, va, vb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if va.hessian[0, 0]*vb.hessian[0, 0] < 0:
            root = brentq(lambda z: at(z).hessian[0, 0], a, b, xtol=2e-12, rtol=2e-12)
            v = at(root)
            rows.append(dict(bracket_samples=samples, actual_grid_points=len(grid),
                grid_kind='near_h_geometric_plus_uniform', a_over_h=root/model.h,
                energy_J_m2=float(UNITS.energy_to_surface(v.energy)),
                traction_MPa=float(UNITS.force_to_traction_mpa(v.gradient[0])),
                curvature_eV_L0sq=v.hessian[0, 0],
                local_type='maximum_traction' if va.hessian[0, 0] > 0 else 'minimum_traction'))
    return rows


def curvature_refinement(basis, coefficients, *, alpha, polarization):
    """Report all FD levels, including roundoff; no best-looking level selection.

    Harmonic analytic columns remain valid at cubic dilation. Direct energies
    independently use the *nonlinear* per-site law at every nonzero amplitude.
    Both quartic truncation and eventual roundoff must remain visible.
    """
    v = np.asarray(polarization, float)
    columns, tails = basis.evaluate((1/3, 1/3, 1/3))
    H = np.einsum('c,cij->ij', coefficients, columns)
    exact = float(v@H@v)
    options = dict(planes=6, mode=2, polarization_plane=v,
                   validation_radius=basis.radius, quadrupole_saturation=alpha)
    zero = basis.direct_sinusoidal_energy(coefficients, amplitude=0., **options)
    rows, previous, previous_fourth = [], None, None
    for h in (4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5, 1.25e-5):
        plus = basis.direct_sinusoidal_energy(coefficients, amplitude=h, **options)
        minus = basis.direct_sinusoidal_energy(coefficients, amplitude=-h, **options)
        central = 2*(plus+minus-2*zero)/h**2
        fourth = None if previous is None else (4*central-previous)/3
        sixth = None if previous_fourth is None else (16*fourth-previous_fourth)/15
        rows.append(dict(amplitude_over_L0=h, analytic_H=exact, central=central,
                         fourth_order=fourth, sixth_order=sixth,
                         central_error=central-exact,
                         fourth_error=None if fourth is None else fourth-exact,
                         sixth_error=None if sixth is None else sixth-exact,
                         harmonic_tail_bound=float(tails@abs(coefficients))))
        previous, previous_fourth = central, fourth
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('preserve prior independent validation')
    began = time.perf_counter()
    source, _, _ = source_and_targets()
    budget, derivatives, truncation, direct, extrema, waves, bindings, widths = [], [], [], [], [], [], [], []
    for directory in args.models:
        data, definition, model, _, _ = load_material(directory)
        if not (definition.get('even_development') or definition.get('normal_development')):
            raise ValueError('explicit v16/v17 dataset required')
        if definition['source_sha256'] != source.reference.sha256:
            raise ValueError('target-source mismatch')
        name = directory.name
        shape = np.asarray(data['best']['decays'])
        c = np.asarray(data['best']['coefficients'])
        alpha = float(shape[3]) if definition.get('quadrupole_saturation_extension') else 0.
        if definition.get('rank_one_range_extension'):
            from .rank_one_range_material import RankOneRangeInterface
            tighter = RankOneRangeInterface(shape, c, tolerance=2e-13)
        else:
            tighter = (SaturatedQuadrupoleInterface(shape, c, tolerance=2e-13) if alpha else
                       QuarticSymmetryInterface(shape[:3], c, tolerance=2e-13))
        bindings.append(dict(model=name,
            calibration_sha256=hashlib.sha256((directory/'calibration.json').read_bytes()).hexdigest(),
            shape=shape, strictly_positive_LJ=bool(np.all(c[:2] > 0)),
            quadrupole_saturation=alpha,
            rank1_decay=shape[4] if definition.get('rank_one_range_extension') else shape[1]))
        if alpha:
            for label, sites, normalization in (
                ('rank2', model.even_sites, model.even_reference_curvature),
                ('Eg', model.eg_sites, 1.)):
                jets, _ = sites.site_jets((model.h, 0., 0.))
                for axis, coordinate in enumerate(('a', 's1', 's2'), start=1):
                    norm2 = float(np.max(np.sum(jets[:, axis]**2, axis=1)))
                    width = np.sqrt(normalization/(alpha*norm2)) if norm2 > 0 else None
                    widths.append(dict(model=name, channel=label, coordinate=coordinate,
                        alpha=alpha, reference_norm_squared=normalization,
                        maximum_site_derivative_norm_squared=norm2,
                        crossover_over_L0=width,
                        crossover_angstrom=None if width is None else width*LENGTH_M/1e-10,
                        interpretation='linear-jet saturation scale, NOT a yield/source/activation length'))
        for label, q in [('inspected039', (model.h, .39, 0.)),
                         ('heldout031', (model.h, .31, 0.)),
                         ('heldout_off', (1.16*model.h, .32, -.12)),
                         ('heldout_opening', (1.93*model.h, 0., 0.))]:
            q = np.asarray(q); value = model.evaluate(q); reference = source.evaluate(q)
            for component, jet in value.components.items():
                budget.append(dict(model=name, state=label, component=component,
                    energy_eV_cell=jet[0], normal_force_eV_L0=jet[1], Haa_eV_L0sq=jet[4],
                    source_Haa_eV_L0sq=reference.hessian[0, 0]))
            tight = tighter.evaluate(q)
            truncation.append(dict(model=name, state=label, tolerance=2e-11, finer_tolerance=2e-13,
                energy_change=abs(tight.energy-value.energy),
                gradient_max_change=float(np.max(abs(tight.gradient-value.gradient))),
                hessian_max_change=float(np.max(abs(tight.hessian-value.hessian)))))
            for step in (4e-5, 2e-5, 1e-5):
                g, H = [], []
                for axis in range(3):
                    d = np.eye(3)[axis]*step
                    plus, minus = model.evaluate(q+d), model.evaluate(q-d)
                    g.append((plus.energy-minus.energy)/(2*step))
                    H.append((plus.gradient-minus.gradient)/(2*step))
                derivatives.append(dict(model=name, state=label, step=step,
                    gradient_max_error=float(np.max(abs(np.array(g)-value.gradient))),
                    hessian_max_error=float(np.max(abs(np.array(H).T-value.hessian)))))
            if alpha:
                exact = np.array([v[0] for v in model.even_jets(q)[:2]])
                for radius in (6, 12, 18):
                    finite = model.direct_even_energy(q, radius=radius, layers=radius)
                    for j, channel in enumerate(('rank2', 'Eg')):
                        direct.append(dict(model=name, state=label, channel=channel,
                            radius=radius, layers=radius, reciprocal=exact[j], direct=finite[j],
                            absolute_error=abs(finite[j]-exact[j]), validation_cutoff_only=True))
            if definition.get('rank_one_range_extension'):
                moment = next(m for _, m in model.base.moments if m.invariant.rank == 1)
                rank1_exact = model.base._angular(q, moment)[0][0]
                for radius in (6, 12, 18):
                    finite = model.direct_rank1_energy(q, radius=radius, layers=radius)
                    direct.append(dict(model=name, state=label, channel='rank1_independent_range',
                        radius=radius, layers=radius, reciprocal=rank1_exact, direct=finite,
                        absolute_error=abs(finite-rank1_exact), validation_cutoff_only=True))
        for samples in (129, 257):
            extrema.extend(dict(model=name, **row) for row in resolved_opening_extrema(model, samples))
        basis = IsotropicBulkBasis(shape[:3], stretch=1.007, radius=12.,
            rank1_decay=shape[4] if definition.get('rank_one_range_extension') else None)
        for j, v in enumerate(np.eye(3)):
            waves.extend(dict(model=name, polarization=j, alpha=alpha, **row)
                         for row in curvature_refinement(basis, c, alpha=alpha, polarization=v))
        print('independent nonlinear/tail/opening checks complete', name, flush=True)
    for samples in (129, 257):
        extrema.extend(dict(model='source_Al99', **row) for row in resolved_opening_extrema(source, samples))
    for name, rows in [('term_budget', budget), ('derivative_refinement', derivatives),
                       ('reciprocal_refinement', truncation), ('direct_even_sums', direct),
                       ('opening_extrema', extrema), ('nonlinear_finite_q', waves),
                       ('saturation_crossover', widths)]:
        write_csv(args.out/(name+'.csv'), rows)
    save_json(args.out/'scope.json', dict(completed=True, bindings=bindings,
        source_sha256=source.reference.sha256, per_site_nonlinearity=True,
        canonical_finite_neighbor_cutoff=False, whole_zone_proof=False,
        validation_state_labels_are_historical=True,
        blind_fit_validation=False,
        full_material_accepted=False, physical_yield_validated=False,
        physical_seconds=False, physical_Hz=False, production_changed=False,
        elapsed_seconds=time.perf_counter()-began))


if __name__ == '__main__':
    main()
