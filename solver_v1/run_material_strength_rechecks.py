"""Constraint exchange and independent checks; no strength tuning or promotion.

The original finite-grid fit remains in joint_fit. At its fixed radial decays,
add actual negative-traction extremum locations as coefficient constraints.
Never clip the force. A repaired opening path is not a material certificate.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from .run_low_stress_cyclic_diagnostic import FIT, build_surface, write_csv
from .run_vector_registry_audit import save_json, point_row, direct_rows, derivative_rows
from .run_vector_registry_rechecks import first_fold
from .run_vector_material_calibration import source_and_targets, OPENING_RATIOS
from .vector_interface_reference import FullRegistryInterface
from .vector_registry_audit import stationary_state
from .vector_material_calibration import (
    IDEAL_H, UNITS, COEFFICIENTS, VectorCoefficientBasis, build_coefficient_surface,
    observation_matrix, fit_coefficients, coefficient_identifiability, opening_traction_extrema,
)


ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/fcc111_active_interface/material_strength_v10'


def repair(out):
    started=time.perf_counter()
    data=json.loads((BASE/'joint_fit/calibration.json').read_text())
    previous=data['best']; decays=previous['scalar_decay'],previous['angular_decay']
    basis=VectorCoefficientBasis(*decays)
    source,obs,states=source_and_targets(); matrix=observation_matrix(basis,obs)
    ratios=list(OPENING_RATIOS); result=previous; history=[]; passed=False
    for iteration in range(12):
        c=np.array(result['coefficients'])
        model=build_coefficient_surface(*decays,c)
        tight=build_coefficient_surface(*decays,c,tolerance=2e-13)
        extrema=opening_traction_extrema(model,121)
        negative=[]
        for p in extrema:
            q=(p['a_over_h']*IDEAL_H,0.,0.)
            g=tight.evaluate(q).gradient[0]
            error=abs(g-p['force_eV_L0'])
            roundoff=128*np.finfo(float).eps*np.sum(abs(basis.jet(q)[1]*c))
            p['force_resolution_eV_L0']=max(error,roundoff)
            if p['force_eV_L0'] < -p['force_resolution_eV_L0']:
                negative.append(p['a_over_h'])
        history.append(dict(iteration=iteration,coefficients=c,loss=result['squared_loss'],
                            constraint_ratios=ratios.copy(),extrema=extrema))
        save_json(out/'exchange_progress.json',dict(completed=False,history=history))
        print('exchange',iteration,'loss',result['squared_loss'],'negative extrema',negative,flush=True)
        if not negative:
            passed=True; break
        ratios.extend(negative)
        inequalities=np.array([basis.jet((r*IDEAL_H,0.,0.))[1] for r in ratios])
        result=fit_coefficients(matrix,obs,opening_inequalities=inequalities)
        if not result['admissible'] or not result['strictly_positive_LJ_resolved']:
            raise RuntimeError('constraint-exchange QP failed; cannot accept a repaired candidate')
    best=dict(scalar_decay=decays[0],angular_decay=decays[1],**result)
    rows=[]
    for o,p in zip(obs,matrix@best['coefficients']):
        rows.append(dict(target=o.name,reference=o.target,prediction=p,units=o.units,
                         normalized_residual=(p-o.target)/o.scale,role=o.role))
    write_csv(out/'fit_and_heldout_residuals.csv',rows)
    save_json(out/'identifiability.json',coefficient_identifiability(matrix,obs,best['coefficients']))
    save_json(out/'calibration.json',dict(completed=True,best=best,
        parent='joint_fit/calibration.json',radial_decays_held_fixed=True,history=history,
        sign_bracketed_opening_extrema_resolved=passed,coefficient_order=COEFFICIENTS,
        no_force_clipping=True,material_accepted=False,production_promoted=False,
        elapsed_seconds=time.perf_counter()-started))
    save_json(out/'exchange_progress.json',dict(completed=True,history=history,
        sign_bracketed_opening_extrema_resolved=passed))


def verify(out):
    started=time.perf_counter(); fit=json.loads((out/'calibration.json').read_text())['best']
    decays=fit['scalar_decay'],fit['angular_decay']; c=fit['coefficients']
    model=build_coefficient_surface(*decays,c)
    tight=build_coefficient_surface(*decays,c,tolerance=2e-13)
    direct=direct_rows(model); write_csv(out/'independent_direct_lattice.csv',direct)
    derivatives=derivative_rows(model,'opening_exchange')
    write_csv(out/'independent_derivatives.csv',derivatives)
    tolerance=[]
    for q in [(IDEAL_H,0.,0.),(1.09*IDEAL_H,.237,.087),(1.028*IDEAL_H,.5,np.sqrt(3)/6),
              (2.8*IDEAL_H,0.,0.)]:
        a,b=model.evaluate(q),tight.evaluate(q)
        tolerance.append(dict(a_L0=q[0],x_L0=q[1],y_L0=q[2],
            energy_change_eV=abs(a.energy-b.energy),gradient_change_eV_L0=float(max(abs(a.gradient-b.gradient))),
            Hessian_change_eV_L0sq=float(np.max(abs(a.hessian-b.hessian)))))
    write_csv(out/'reciprocal_depth_refinement.csv',tolerance)
    extrema=[]
    for n,m in [(121,model),(241,tight)]:
        extrema.extend(dict(samples=n,**p) for p in opening_traction_extrema(m,n))
    write_csv(out/'extrema_independent_refinement.csv',extrema)
    folds=[]
    for n in (31,61):
        q,v,curvature,residual=first_fold(model,'a_y_free',n)
        folds.append(dict(samples=n,**point_row('opening_exchange','uniform_ideal_NOT_yield',q,v,UNITS),
                          Schur_curvature=curvature,eliminated_force_residual=residual))
    write_csv(out/'fold_independent_refinement.csv',folds)
    source,_,_=source_and_targets()
    old=FullRegistryInterface(build_surface(tolerance=2e-11)[0])
    responses=[]
    # The same independently prescribed tractions in MPa, NOT strain fitting.
    # .1,.5,.8 experimental plastic strains are NOT imposed on this interface.
    loads=(0.,3.,4.595238095238095,6.928571428571429,9.880952380952381,25.,50.,100.,0.,-50.,0.)
    for name,m in [('old_candidate',old),('opening_exchange',model),('Mishin_source_only',source)]:
        q=np.array([m.h,0.,0.]); origin=None
        for step,traction in enumerate(loads):
            force=UNITS.traction_mpa_to_force([0.,traction,0.])
            r=stationary_state(m,q,force=force,expected_index=0)
            if not r['valid']:
                raise RuntimeError('low-stress intact branch not validated')
            q=r['q']; origin=q.copy() if origin is None else origin
            responses.append(dict(model=name,step=step,shear_x_MPa=traction,**{
                'a_L0':q[0],'ux_L0':q[1],'uy_L0':q[2]},
                displacement_from_initial_L0=float(np.linalg.norm(q-origin)),
                force_residual_eV_L0=r['force_residual'],lambda_min=r['eigenvalues'][0],
                static_unload=bool(step and not traction),
                experimental_plastic_strain_prediction_available=False,
                physical_time_hold=False,interpretation='local zero-K rigid-half branch, not defect-containing wire'))
    write_csv(out/'experimental_stress_range_static_check.csv',responses)
    save_json(out/'independent_verification.json',dict(completed=True,elapsed_seconds=time.perf_counter()-started,
        max_direct_fine_energy_error_eV=max(r['abs_error_eV_cell'] for r in direct if r['direct_radial_index']==72),
        max_gradient_FD_error=max(r['gradient_max_abs_error_eV_L0'] for r in derivatives),
        max_Hessian_FD_error=max(r['hessian_max_abs_error_eV_L0sq'] for r in derivatives),
        extrema=extrema,material_accepted=False,finite_source_implemented=False,
        physical_Hz=False,production_promoted=False))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['repair','verify'])
    p.add_argument('--output',type=Path,default=BASE/'opening_exchange')
    args=p.parse_args()
    (repair if args.stage=='repair' else verify)(args.output)


if __name__=='__main__':
    main()
