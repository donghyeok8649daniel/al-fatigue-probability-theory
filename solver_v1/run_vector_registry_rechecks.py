"""Independent constraint/stencil checks after the vector_registry_v9 run.

No parameter fit. Fixed-path maxima are labeled by their holding tractions;
they must not be confused with an unrestrained three-coordinate spinodal.
"""
from functools import lru_cache
import time

import numpy as np
from scipy.optimize import brentq, least_squares

from .reference_eam_targets import MishinRigidFCCReference
from .run_low_stress_cyclic_diagnostic import build_surface, write_csv
from .run_vector_registry_audit import OUT, point_row, save_json
from .vector_interface_reference import FullRegistryInterface, MishinVectorInterfaceReference
from .vector_registry_audit import relax_at_x


def constrained_point(model, x, mode, guess):
    if mode == 'a_y_fixed':
        q = np.array([model.h, x, 0.]); out = model.evaluate(q)
        return q, out, out.hessian[1, 1], 0.
    if mode == 'a_free_y_fixed':
        @lru_cache(maxsize=1)
        def value(a):
            return model.evaluate([a, x, 0.])
        solved = least_squares(lambda a: [value(float(a[0])).gradient[0]], [guess[0]],
            jac=lambda a: [[value(float(a[0])).hessian[0, 0]]],
            bounds=([.7*model.h], [1.7*model.h]), xtol=2e-13, ftol=2e-13, gtol=2e-13)
        q = np.array([solved.x[0], x, 0.]); out = value(float(q[0]))
        if abs(out.gradient[0]) > 2e-9 or out.hessian[0, 0] <= 0:
            raise RuntimeError('normal-only branch unresolved')
        return q, out, out.hessian[1, 1]-out.hessian[0, 1]**2/out.hessian[0, 0], abs(out.gradient[0])
    result = relax_at_x(model, x, guess[[0, 2]])
    return result['q'], result['evaluation'], result['schur_curvature'], result['transverse_force_residual']


def first_fold(model, mode, n):
    guess = np.array([model.h, 0., 0.]); previous = None
    for x in np.linspace(0., .45, n):
        q, out, curvature, residual = constrained_point(model, x, mode, guess)
        if previous is not None and previous[2] > 0 >= curvature:
            left = previous[0][1]; right = x; reference_q = (previous[0]+q)/2
            def value(t):
                return constrained_point(model, t, mode, reference_q)
            root = brentq(lambda t: value(t)[2], left, right, xtol=2e-11)
            return value(root)
        guess = q; previous = (q, out, curvature, residual)
    raise RuntimeError('first force maximum not bracketed; not called a spinodal')


def main():
    started = time.perf_counter()
    surface, units, metadata = build_surface(tolerance=2e-11)
    models = [('analytic_candidate', FullRegistryInterface(surface)),
              ('Mishin_source_only', MishinVectorInterfaceReference(MishinRigidFCCReference(), units.length_scale_m/1e-10))]
    rows = []; stencil_rows = []
    for name, model in models:
        for mode in ('a_y_fixed', 'a_free_y_fixed', 'a_y_free'):
            for n in (31, 61):
                q, out, curvature, residual = first_fold(model, mode, n)
                row = point_row(name, mode, q, out, units)
                row.update(samples=n, constrained_curvature_eV_L0sq=curvature,
                    solved_force_residual_eV_L0=residual,
                    required_holding_normal_MPa=float(units.force_to_traction_mpa(out.gradient[0])) if mode == 'a_y_fixed' else 0.,
                    required_holding_transverse_MPa=float(units.force_to_traction_mpa(out.gradient[2])) if mode != 'a_y_free' else 0.,
                    full_coordinate_spinodal=(mode == 'a_y_free'), experimental_yield=False)
                rows.append(row)
            print(name, mode, row['traction_x_MPa'], 'MPa, static ideal only', flush=True)
            write_csv(OUT/'constraint_ablation.csv', rows)
        # Source spline knot crossings affect finite-difference stencils. Do
        # not hide the coarser errors in derivative_validation.csv.
        q = np.array([1.51*model.h, .307, -.071]); out = model.evaluate(q)
        for step in (1e-5, 5e-6, 2.5e-6, 1.25e-6, .625e-6):
            H = np.column_stack([(model.evaluate(q+step*e).gradient-
                                  model.evaluate(q-step*e).gradient)/(2*step) for e in np.eye(3)])
            error = float(np.max(abs(out.hessian-H)))
            stencil_rows.append(dict(model=name, step_L0=step, a_L0=q[0], ux_L0=q[1], uy_L0=q[2],
                hessian_max_abs_error_eV_L0sq=error,
                hessian_relative_norm_error=float(np.linalg.norm(out.hessian-H)/np.linalg.norm(out.hessian))))
        if stencil_rows[-1]['hessian_relative_norm_error'] > 2e-7:
            raise RuntimeError('fine-stencil SOURCE/analytic derivative validation failed')
    write_csv(OUT/'source_stencil_refinement.csv', stencil_rows)
    save_json(OUT/'rechecks_status.json', dict(completed=True, wall_seconds=time.perf_counter()-started,
        parameter_sha256=metadata['parameter_sha256'], no_refit=True, all_first_folds_bracketed=True,
        ablations=rows, source_stencil_refinement=stencil_rows,
        interpretation='removing artificial path constraints lowers IDEAL strength, not a calibration to real yield'))


if __name__ == '__main__':
    main()
