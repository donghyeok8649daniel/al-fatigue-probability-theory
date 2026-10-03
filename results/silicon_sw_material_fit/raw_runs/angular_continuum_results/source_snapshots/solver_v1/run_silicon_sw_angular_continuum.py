"""Exhaustive one-dimensional polynomial candidates on completed angle profiles.

Energy, force and virial are exactly quadratic in beta at fixed shape/anchors.
Verify this identity on every completed grid design before forming the quartic
loss and bulk guard. Enumerate real stationary/guard roots and both endpoints.
No new geometries, no extra angular coefficients, no held-out selection.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from .run_silicon_sw_material_fit import dump


def real_roots(coefficients,low,high):
    c=np.trim_zeros(np.array(coefficients,float),'b')
    if len(c)<2:return []
    roots=np.polynomial.polynomial.polyroots(c)
    result=[]
    for r in roots:
        if abs(r.imag)<1e-9*max(1.,abs(r.real)) and low-1e-12<=r.real<=high+1e-12:
            result.append(float(r.real))
    return result


def squared_polynomial(d):
    return np.array([d[0]@d[0],2*d[0]@d[1],d[1]@d[1]+2*d[0]@d[2],2*d[1]@d[2],d[2]@d[2]])


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    raw=Path(__file__).read_bytes();target=args.output/'source_snapshots/solver_v1/run_silicon_sw_angular_continuum.py'
    target.parent.mkdir(parents=True);target.write_bytes(raw)
    dump(args.output/'source_manifest.json',[dict(path='solver_v1/run_silicon_sw_angular_continuum.py',sha256=hashlib.sha256(raw).hexdigest())])
    profiles=json.loads((args.profiles/'profiles.json').read_text());groups={}
    for r in profiles:groups.setdefault(r['anchor_profile'],[]).append(r)
    rows=[];allbest=None;checks=[];width=.375
    for anchor,group in sorted(groups.items()):
        controls=[next(r for r in group if r['angle_beta']==v) for v in [-width,0.,width]]
        arrays=[]
        for r in controls:
            d=np.load(args.profiles/f"profile_{r['profile']}_design.npz")
            arrays.append({k:d[k].copy() for k in ['design','target','bulk_design','bulk_target','maximum_bulk_sq']})
        minus,zero,plus=arrays;alpha=np.array(controls[1]['amplitudes'])
        expanded={}
        for key in ['design','bulk_design']:
            expanded[key]=np.stack([zero[key],(plus[key]-minus[key])/(2*width),
                                   (plus[key]+minus[key]-2*zero[key])/(2*width**2)])
        eos=np.array([r['EOS_error_eV_atom'] for r in controls])
        ep=np.stack([eos[1],(eos[2]-eos[0])/(2*width),(eos[2]+eos[0]-2*eos[1])/(2*width**2)])
        maximum_error=0.
        for r in group:
            beta=r['angle_beta'];d=np.load(args.profiles/f"profile_{r['profile']}_design.npz")
            for key in ['design','bulk_design']:
                reconstructed=np.einsum('i,ijk->jk',[1.,beta,beta**2],expanded[key])
                error=float(np.max(abs(reconstructed-d[key])));maximum_error=max(maximum_error,error)
                if error>2e-9:raise RuntimeError('quadratic beta design identity failed')
            if max(abs(ep[0]+beta*ep[1]+beta**2*ep[2]-r['EOS_error_eV_atom']))>2e-10:
                raise RuntimeError('quadratic beta EOS identity failed')
        checks.append(dict(anchor_profile=anchor,maximum_design_reconstruction_error=maximum_error,grid_profiles_verified=len(group)))
        residual=np.einsum('ijk,k->ij',expanded['design'],alpha);residual[0]-=zero['target']
        loss_vectors=np.concatenate([residual,ep/.02],axis=1);loss=squared_polynomial(loss_vectors)
        br=np.einsum('ijk,k->ij',expanded['bulk_design'],alpha);br[0]-=zero['bulk_target']
        guard=squared_polynomial(br);bound=float(zero['maximum_bulk_sq']);gc=guard.copy();gc[0]-=bound
        derivative=np.polynomial.polynomial.polyder(loss)
        roots=real_roots(derivative,-width,width);boundaries=real_roots(gc,-width,width)
        candidates=[(-width,'left endpoint'),(width,'right endpoint')]+[(v,'stationary root') for v in roots]+[(v,'bulk guard root') for v in boundaries]
        accepted=[];evaluated=[]
        for beta,kind in candidates:
            if not -width<=beta<=width:continue
            vector=loss_vectors[0]+beta*loss_vectors[1]+beta**2*loss_vectors[2]
            bv=br[0]+beta*br[1]+beta**2*br[2];actual_guard=float(bv@bv)
            objective=float(vector@vector)
            item=dict(beta=beta,kind=kind,objective_squared=objective,bulk_energy_squared=actual_guard,
                bulk_guard_passed=actual_guard<=bound+2e-9,loss_derivative=float(np.polynomial.polynomial.polyval(beta,derivative)),
                bulk_boundary_residual=actual_guard-bound)
            evaluated.append(item)
            if item['bulk_guard_passed']:accepted.append(item)
        row=dict(anchor_profile=anchor,loss_polynomial_ascending=loss.tolist(),bulk_polynomial_ascending=guard.tolist(),
            maximum_bulk_energy_squared=bound,stationary_roots=roots,guard_roots=boundaries,candidates=evaluated,
            positive_coefficients_and_elastic_anchors_from_prior_profile=True)
        if accepted:
            best=min(accepted,key=lambda r:r['objective_squared']);row['selected']=best
            winner=dict(best,anchor_profile=anchor,parameters=controls[1]['parameters'],amplitudes=alpha.tolist(),shape=controls[1]['shape'])
            if allbest is None or best['objective_squared']<allbest['objective_squared']:allbest=winner
        rows.append(row);dump(args.output/'profiles.json',rows)
    dump(args.output/'quadratic_design_identity.json',checks)
    dump(args.output/'selected_fit.json',allbest)
    dump(args.output/'summary.json',dict(complete=True,shape_profiles=len(rows),source_grid_profiles=len(profiles),
        all_grid_designs_verified=True,selected_anchor=allbest['anchor_profile'] if allbest else None,
        selected_beta=allbest['beta'] if allbest else None,new_geometry_evaluations=0,
        mathematical_search='all floating-point real roots of quartic guard and cubic derivative plus endpoints; no interval-arithmetic global certificate',
        no_global_shape_optimum_claim=True,material_approved=False,first_initiation_validated=False,
        physical_clock_validated=False,elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profiles',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
