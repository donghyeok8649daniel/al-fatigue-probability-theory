"""Independent static checks of completed v14 fits, never a refit/PDE run."""
import argparse
from pathlib import Path
import json
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit

from .range_resolved_material import build_range_surface,RangeObservationCache
from .symmetry_resolved_material import (SymmetryResolvedInterface,SymmetryObservationCache,
    CubicSymmetryInterface,CubicSymmetryObservationCache,SymmetryTailBulkBasis,CubicSymmetryTailBulkBasis)
from .quartic_angular_material import QuarticSymmetryInterface,QuarticSymmetryObservationCache,QuarticSymmetryTailBulkBasis
from .tail_constrained_material import TailBulkCoefficientBasis,declared_wavepoints
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json,point_row
from .run_low_stress_cyclic_diagnostic import write_csv
from .vector_registry_audit import stationary_state
from .vector_material_calibration import UNITS,IDEAL_H,LENGTH_M
from .yield_elastic_metric import cubic_metric_problem,cubic_to_mode_matrix
from .material_calibration_controls import append_fixed_pair


def load_material(directory):
    data=json.loads((directory/'calibration.json').read_bytes())
    if not data['completed']: raise ValueError('completed actual calibration required')
    definition=json.loads((directory/'definition.json').read_bytes())
    best=data['best'];decays=best['decays'];c=np.asarray(best['coefficients'])
    if not best['strictly_positive_LJ'] or not data.get('admissible_solution_found',True):
        raise ValueError('diagnostic zero-pair closure is not an admissible LJ material')
    if definition.get('quadrupole_saturation_extension'):
        from .quadrupole_saturation import SaturatedQuadrupoleInterface,SaturatedQuadrupoleCache,SaturatedQuadrupoleBulk
        model=SaturatedQuadrupoleInterface(decays,c)
        cache_type,bulk_type=SaturatedQuadrupoleCache,SaturatedQuadrupoleBulk
    elif definition.get('density_angular_cross_extension'):
        from .density_angular_cross_material import CrossDensityAngularInterface,CrossDensityAngularCache,CrossDensityAngularBulk
        model=CrossDensityAngularInterface(decays,c)
        cache_type,bulk_type=CrossDensityAngularCache,CrossDensityAngularBulk
    elif definition.get('density_mixture_extension'):
        from .mixed_density_material import MixedDensityQuarticInterface, MixedDensityObservationCache, MixedDensityTailBulkBasis
        model=MixedDensityQuarticInterface(decays,c)
        cache_type,bulk_type=MixedDensityObservationCache,MixedDensityTailBulkBasis
    elif definition.get('quartic_angular_extension'):
        model=QuarticSymmetryInterface(decays[:3],c,saturation=decays[3] if len(decays)==4 else 0.)
        cache_type,bulk_type=QuarticSymmetryObservationCache,QuarticSymmetryTailBulkBasis
    elif definition.get('cubic_density_extension'):
        model=CubicSymmetryInterface(decays,c);cache_type,bulk_type=CubicSymmetryObservationCache,CubicSymmetryTailBulkBasis
    elif definition.get('extra_amplitude'):
        model=SymmetryResolvedInterface(decays,c);cache_type,bulk_type=SymmetryObservationCache,SymmetryTailBulkBasis
    else:
        model=build_range_surface(*decays,c);cache_type,bulk_type=RangeObservationCache,TailBulkCoefficientBasis
    return data,definition,model,cache_type,bulk_type


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--grid-step',type=float,default=.1)
    parser.add_argument('--out',type=Path,help='fresh independent validation directory')
    args=parser.parse_args();directory=args.directory;out=args.out or directory/'validation'
    if out.exists(): raise FileExistsError('preserve existing validation')
    began=time.perf_counter();data,definition,model,cache_type,bulk_type=load_material(directory)
    source,observations,states=source_and_targets();best=data['best'];c=np.asarray(best['coefficients']);decays=best['decays']
    if definition.get('interface_development'):
        from .interface_development_targets import development_observations
        observations,_=development_observations(source,observations,states)
    if definition.get('even_development'):
        from .interface_even_development_targets import even_development_observations
        observations,_=even_development_observations(source,observations,states)
    if source.reference.sha256!=definition['source_sha256']: raise ValueError('source hash mismatch')
    cache=cache_type(observations);raw=cache.matrix(decays)
    matrix,obs=cubic_metric_problem(raw[:,:8],observations)
    if raw.shape[1]>8:
        extra=raw[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5])
        matrix=np.column_stack([matrix,extra])
    bloch_targets=()
    if definition.get('static_bloch_targets'):
        from .static_bloch_targets import source_bloch_targets,append_bloch_problem
        bloch_targets=source_bloch_targets(source.reference)
        basis=bulk_type(decays if definition.get('density_mixture_extension') else decays[:3],radius=definition['radius_over_L0'])
        matrix,obs,_=append_bloch_problem(matrix,obs,basis,bloch_targets)
    matrix,obs=append_fixed_pair(matrix,obs,definition.get('fixed_pair_control'))
    predictions=matrix@c
    if np.max(abs(predictions-np.asarray(best['predictions'])))>2e-8:
        raise ArithmeticError('saved fit fails exact analytic replay')
    roots=[];root_data={};curve=[]
    for label,index in (('perfect',0),('fault',0),('saddle',1)):
        root=stationary_state(model,states[label],expected_index=index)
        root_data[label]=dict(q=root['q'],valid=root['valid'],morse_index=root['morse_index'] if 'morse_index' in root else None)
        roots.append(dict(state=label,valid=root['valid'],**{k:v for k,v in point_row('candidate',label,root['q'],root['evaluation'],UNITS).items() if k!='state'}))
    write_csv(out/'stationary_states.csv',roots)
    # Explicit registry and opening development/independent curve points.
    tau=np.array([.5,np.sqrt(3)/6])
    for path,direction in (('shockley',tau),('direct110',np.array([1.,0.]))):
        for fraction in np.linspace(0.,1.,41):
            q=np.r_[IDEAL_H,fraction*direction];v=model.evaluate(q);s=source.evaluate(q)
            curve.append(dict(path=path,fraction=fraction,a_over_h=1.,
                candidate_J_m2=float(UNITS.energy_to_surface(v.energy)),source_J_m2=float(UNITS.energy_to_surface(s.energy)),
                normal_traction_MPa=float(UNITS.force_to_traction_mpa(v.gradient[0]))))
    for ratio in np.r_[np.linspace(1.,2.,41),np.linspace(2.1,4.,20),8.,20.,40.]:
        q=np.array([ratio*IDEAL_H,0.,0.]);v=model.evaluate(q);s=source.evaluate(q)
        curve.append(dict(path='opening',fraction=0.,a_over_h=ratio,
            candidate_J_m2=float(UNITS.energy_to_surface(v.energy)),source_J_m2=float(UNITS.energy_to_surface(s.energy)),
            normal_traction_MPa=float(UNITS.force_to_traction_mpa(v.gradient[0]))))
    write_csv(out/'interface_curves.csv',curve)
    print('stationary and curve checks complete',flush=True)
    points=declared_wavepoints(args.grid_step);waves=[];previous=None;worst=None
    for radius in (12.,16.):
        basis=bulk_type(decays if definition.get('density_mixture_extension') else decays[:3],radius=radius);current=[]
        for i,q in enumerate(points):
            columns,tails=basis.evaluate(q);H=np.einsum('c,cij->ij',c,columns)
            values=np.linalg.eigvalsh(H);bound=float(tails@abs(c));current.append(H)
            waves.append(dict(radius_over_L0=radius,qx=q[0],qy=q[1],qz=q[2],minimum_H=values[0],
                tail_bound=bound,robust_margin=values[0]-bound,
                radius_change=None if previous is None else float(np.linalg.norm(H-previous[i],2))))
        previous=current
        worst=min(waves[-len(points):],key=lambda r:r['robust_margin'])
        print(f"radius {radius}: {len(points)} wavevectors, worst margin={worst['robust_margin']:.7g}",flush=True)
    write_csv(out/'finite_q_validation.csv',waves)
    # Continuous local minimizations over the irreducible wedge. Normalizing
    # by |q|^2 distinguishes acoustic zero from a true soft finite-q mode.
    basis=bulk_type(decays if definition.get('density_mixture_extension') else decays[:3],radius=16.)
    def objective(q):
        cols,_=basis.evaluate(q);H=np.einsum('c,cij->ij',c,cols)
        return np.linalg.eigvalsh(H)[0]/(q@q)
    starts=[np.array([r['qx'],r['qy'],r['qz']]) for r in sorted(waves[-len(points):],
        key=lambda r:r['minimum_H']/(r['qx']**2+r['qy']**2+r['qz']**2))[:3]]
    minima=[]
    for start in starts:
        run=minimize(objective,start,method='SLSQP',bounds=[(1e-3,1.),(0.,1.),(0.,1.)],
            constraints=[dict(type='ineq',fun=lambda q:np.r_[q[0]-q[1],q[1]-q[2],1.5-q.sum()])],
            options=dict(maxiter=70,ftol=1e-9))
        cols,tail=basis.evaluate(run.x);ev=float(np.linalg.eigvalsh(np.einsum('c,cij->ij',c,cols))[0])
        minima.append(dict(start=start,q=run.x,minimum_H=ev,tail_bound=float(tail@abs(c)),
            success=bool(run.success),message=str(run.message),nfev=run.nfev))
    save_json(out/'continuous_q_search.json',minima)
    norm=(predictions-np.array([o.target for o in obs]))/np.array([o.scale for o in obs])
    write_csv(out/'fit_replay.csv',[dict(observable=o.name,role=('exact' if definition.get('exact_bulk_stage') and i<5 else o.role),
        units=o.units,target=o.target,prediction=p,scale=o.scale,normalized_residual=r)
        for i,(o,p,r) in enumerate(zip(obs,predictions,norm))])
    fit_selected=[i for i,o in enumerate(obs) if o.role=='fit'];hold=[i for i,o in enumerate(obs) if o.role=='heldout']
    # Local sensitivity at actual coefficients, including radial/shape parameters.
    equal=matrix[:2]
    tangent=np.vstack([-np.linalg.solve(equal[:,:2],equal[:,2:]),np.eye(len(c)-2)])
    scale=np.array([o.scale for o in obs]);normalization=np.maximum(1.,abs(c[2:]))
    jac=(matrix[fit_selected]@tangent)*normalization/scale[fit_selected,None]
    radial_refinements=[];matrix_derivative_refinements=[]
    for step in (2e-4,1e-4):
        radial=[];matrix_derivatives=[]
        for axis in range(len(decays)):
            ps=[];matrices=[]
            for sign in (-1,1):
                d=np.array(decays)
                if definition.get('density_mixture_extension') and axis==3:
                    d[axis]=expit(logit(d[axis])+sign*step)
                else:d[axis]*=np.exp(sign*step)
                r=cache.matrix(d)
                m,_=cubic_metric_problem(r[:,:8],observations)
                if r.shape[1]>8:
                    extra=r[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5]);m=np.column_stack([m,extra])
                if bloch_targets:
                    basis=bulk_type(d if definition.get('density_mixture_extension') else d[:3],radius=definition['radius_over_L0'])
                    m,_,_=append_bloch_problem(m,observations,basis,bloch_targets)
                m,_=append_fixed_pair(m,observations,definition.get('fixed_pair_control'))
                matrices.append(m)
                z=c.copy();z[:2]=np.linalg.solve(m[:2,:2],np.array([o.target for o in obs[:2]])-m[:2,2:]@c[2:])
                ps.append(m@z)
            radial.append((ps[1]-ps[0])[fit_selected]/(2*step*scale[fit_selected]))
            matrix_derivatives.append((matrices[1]-matrices[0])/(2*step))
        radial_refinements.append(np.column_stack(radial))
        matrix_derivative_refinements.append(np.array(matrix_derivatives))
    full=np.column_stack([jac,radial_refinements[-1]])
    _,singular,vt=np.linalg.svd(full,full_matrices=False)
    norms=np.linalg.norm(full,axis=0)
    save_json(out/'identifiability.json',dict(jacobian=full,singular_values=singular,
        radial_derivative_refinement_change=float(np.max(abs(radial_refinements[0]-radial_refinements[1]))),
        radial_derivative_steps=[2e-4,1e-4],right_singular_vectors=vt,
        column_cosines=full.T@full/(norms[:,None]*norms[None,:]),
        coefficient_normalization=normalization,condition_number=float(singular[0]/singular[-1]),
        rank=int(np.linalg.matrix_rank(full)),statistical_confidence_claimed=False,
        coordinate_system=('force/cohesion eliminated; normalized energy coefficients + log microscopic ranges + logit density mixture'
            if definition.get('density_mixture_extension') else
            'force/cohesion eliminated; normalized energy coefficients + log microscopic ranges/shape')))
    # The historical diagnostic above allows the three elastic constants to
    # vary. It is NOT the tangent of an exact-bulk stage. Record the actual
    # staged equality manifold separately rather than interpreting numerical
    # full rank in a larger parameter space as identifiability of this fit.
    from .profiled_identifiability import equality_tangent_sensitivity
    exact_rows=sorted(set(([0,1,2,3,4] if definition.get('exact_bulk_stage') else [])
                         +[i for i,o in enumerate(obs) if o.role=='exact']))
    tangent_checks=[equality_tangent_sensitivity(matrix,obs,c,d,exact_rows=exact_rows)
                    for d in matrix_derivative_refinements]
    exact_tangent=tangent_checks[-1]
    exact_tangent['shape_derivative_steps']=[2e-4,1e-4]
    exact_tangent['jacobian_step_refinement_change']=float(np.max(abs(
        tangent_checks[1]['jacobian']-tangent_checks[0]['jacobian'])))
    save_json(out/'exact_stage_identifiability.json',exact_tangent)
    sampled_stable=bool(worst['robust_margin']>0)
    continuous_stable=all(r['success'] and r['minimum_H']>r['tail_bound'] for r in minima)
    save_json(out/'decision.json',dict(completed=True,source_sha256=source.reference.sha256,
        length_scale_m=LENGTH_M,energy_scale_J=1.602176634e-19,cubic_GPa=predictions[2:5],
        cohesive_eV_atom=float(predictions[1]),normal_force_eV_strain=float(predictions[0]),
        root_data=root_data,finite_q_grid_positive_above_tail=sampled_stable,
        continuous_searches_positive_above_tail=continuous_stable,
        whole_zone_proof=False,heldout_normalized_rms=float(np.sqrt(np.mean(norm[hold]**2))),
        bulk_calibration_only=bool(np.max(abs(norm[:5]))<1e-7 and sampled_stable
            and continuous_stable and root_data['perfect']['valid']),
        full_Al_material_accepted=False,experimental_yield_validated=False,
        physical_seconds=False,physical_Hz=False,elapsed_seconds=time.perf_counter()-began))
    print('completed independent validation',flush=True)


if __name__=='__main__': main()
