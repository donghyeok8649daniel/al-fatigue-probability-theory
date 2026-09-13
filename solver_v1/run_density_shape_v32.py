"""Reference auxiliary-density shape audit, NOT a new Al material calibration."""
import argparse
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares
from .reference_eam_targets import MishinRigidFCCReference, SOURCE_URL
from .polynomial_exponential_density import QuadraticEnvelopeDensity
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def run(out):
    out = Path(out)
    if out.exists():
        raise FileExistsError('fresh reference-shape audit required')
    source = MishinRigidFCCReference()
    fit_r = np.arange(2., 6.001, .25)
    held_r = np.arange(2.125, 6., .25)
    starts = [(1., 2., 1.8, .3), (.1, 1., 1., 1.), (10., 3., 2., .7)]
    out.mkdir(parents=True)
    save_json(out/'definition.json', dict(source_sha256=source.sha256, source_url=SOURCE_URL,
        source_role='auxiliary density of target-only Al99, not observable electron density',
        radius_unit='Angstrom', density_unit='original setfl auxiliary density gauge',
        fit_radii=fit_r, excluded_radii=held_r, quadratic_starts=starts,
        discrepancy='max(.1*target_density,.01*maximum_training_density)',
        fitting='values only; derivatives and midpoint values are validation',
        material_accepted=False, production_changed=False))
    target = source._rho(fit_r)
    if np.any(target <= 0) or np.any(source._rho(held_r) <= 0):
        raise ValueError('positive source density required for logarithmic shape audit')
    scale = np.maximum(.1*target, .01*target.max())
    began = time.perf_counter()
    fits, rows = [], []
    base = least_squares(lambda p: (np.exp(p[0]-np.exp(p[1])*fit_r)-target)/scale,
        np.log([.5, 1.]), bounds=(np.log([1e-10, .01]), np.log([1e8, 30.])),
        jac='3-point', max_nfev=1000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    candidates = [('single_exponential', base, None)]
    for index, start in enumerate(starts):
        def decode(p):
            return QuadraticEnvelopeDensity(np.exp(p[0]), np.exp(p[1]), p[2], np.exp(p[3]))
        fit = least_squares(lambda p: (decode(p).radial(fit_r)-target)/scale,
            [np.log(start[0]), np.log(start[1]), start[2], np.log(start[3])],
            bounds=([np.log(1e-10), np.log(.01), 0., np.log(.01)],
                    [np.log(1e8), np.log(30.), 6., np.log(5.)]),
            jac='3-point', max_nfev=1000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        candidates.append((f'quadratic_start_{index}', fit, decode(fit.x)))
    for name, fit, model in candidates:
        singular = np.linalg.svd(fit.jac, compute_uv=False)
        fits.append(dict(name=name, optimizer_success=bool(fit.success), message=str(fit.message),
            parameters=dict(C=np.exp(fit.x[0]), k=np.exp(fit.x[1])) if model is None else vars(model),
            squared_loss=float(fit.fun@fit.fun), optimality=float(fit.optimality), nfev=fit.nfev,
            singular_values=singular, jacobian=fit.jac, active_mask=fit.active_mask,
            parameter_coordinates='log(C),log(k),center_A,log(width_A)' if model else 'log(C),log(k)',
            material_accepted=False))
        for role, radii in [('fit', fit_r), ('heldout', held_r)]:
            for r in radii:
                value, first, second = [float(source._rho(r, order)) for order in (0, 1, 2)]
                if model is None:
                    C, k = np.exp(fit.x)
                    pred = C*np.exp(-k*r); fp, fpp = -k*pred, k*k*pred
                else:
                    pred, fp, fpp = [float(model.radial(r, order)) for order in (0, 1, 2)]
                rows.append(dict(model=name, role=role, radius_A=r, target_density=value,
                    predicted_density=pred, normalized_value_residual=(pred-value)/max(.1*value,.01*target.max()),
                    reference_log_slope=first/value, predicted_log_slope=fp/pred,
                    reference_log_curvature=second/value-(first/value)**2,
                    predicted_log_curvature=fpp/pred-(fp/pred)**2))
    write_csv(out/'radial_comparison.csv', rows)
    save_json(out/'fits.json', fits)
    nearest = source.geometry.b
    f, fp, fpp = [float(source._rho(nearest, n)) for n in (0, 1, 2)]
    summary = dict(completed=True, elapsed_seconds=time.perf_counter()-began,
        source_sha256=source.sha256, nearest_neighbor_A=nearest,
        nearest_log_curvature_per_A2=fpp/f-(fp/f)**2,
        profiles=[dict(name=item['name'], squared_loss=item['squared_loss'],
                       optimizer_success=item['optimizer_success']) for item in fits],
        minimum_density_only_fit_loss=min(item['squared_loss'] for item in fits),
        positive_exponential_mixture_cannot_match_negative_log_curvature=True,
        whole_energy_family_impossibility_proved=False, material_accepted=False,
        production_changed=False, physical_time_calibrated=False)
    save_json(out/'summary.json', summary)
    print(summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--out', type=Path, required=True)
    run(parser.parse_args().out)
