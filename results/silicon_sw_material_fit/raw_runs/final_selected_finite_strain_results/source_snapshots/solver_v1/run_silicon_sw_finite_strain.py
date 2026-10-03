"""Conservative stress fits with the historical cubic strain pattern.

Reads a primary-source protocol but does not execute the old optimizer or
quippy engine. Internal displacements are resolved by independent atomic-force
roots. Thus this isolates finite strain/pattern effects, not the whole original
benchmark environment. No universal correction is applied to DFT targets.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import root
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_sw_material_fit import force_features
from .silicon_sw_anharmonic import AngularDiamondCell
from .run_silicon_sw_material_fit import dump


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();base=Path(__file__).resolve().parents[1]
    selected=json.loads(args.selected.read_text());published=json.loads(args.published.read_text());ref=json.loads(args.reference.read_text())
    p,_=source_parameters();p=dict(p,epsilon=2.1675)
    beta=selected['angle_beta'] if 'angle_beta' in selected else selected['selected']['beta']
    configs=[('selected_candidate',selected['parameters'],beta,ref['diamond_a0'])]
    if not args.candidate_only:configs.insert(0,('published_SW_parameters',p,0.,published['diamond_a0']))
    manifests=[]
    for relative in ['solver_v1/run_silicon_sw_finite_strain.py','solver_v1/silicon_sw_material_fit.py',
        'solver_v1/silicon_sw_anharmonic.py','solver_v1/silicon_environment_research.py']:
        raw=(base/relative).read_bytes();target=args.output/'source_snapshots'/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        manifests.append(dict(path=relative,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifests)
    for source in [args.selected,args.published,args.reference,args.historical_manifest]:
        target=args.output/'input_snapshots'/source.name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(source.read_bytes())
    dump(args.output/'protocol.json',dict(patterns=['exx=gamma_yz=e','pure gamma_yz=e control'],
        samples_per_fit=5,strain_increments=[.01,.003,.001],optimizer='atomic-force root, residual <=2e-8 eV/A',
        stress='conservative all-atom current-configuration strain derivative / actual volume',
        old_matscipy_reproduced=False,old_quippy_reproduced=False,DFT_corrected_by_SW=False,
        historical_source=json.loads(args.historical_manifest.read_text()),new_DFT=0,new_MD=0))
    states=[];fits=[];calls=0
    for name,parameters,beta,lattice in configs:
        cell=lattice/2*np.array([[0,1,1],[1,0,1],[1,1,0.]])
        basis=np.array([[0.,0.,0.],[lattice/4]*3]);pbc=np.ones(3,bool)
        model=AngularDiamondCell(parameters,lattice,angle_beta=beta);j=model.evaluate(np.zeros(9))
        hh=j.hessian;rel=hh[:6,:6]-hh[:6,6:]@np.linalg.solve(hh[6:,6:],hh[6:,:6]);scale=160.2176634/model.volume
        exact=rel[[0,0,3],[0,1,3]]*scale
        for pattern in ['cubic_exx_plus_yz','pure_yz']:
            for increment in [.01,.003,.001]:
                group=[]
                for e in increment*np.arange(-2,3):
                    f=np.eye(3);f[1,2]=f[2,1]=e/2
                    if pattern=='cubic_exx_plus_yz':f[0,0]+=e
                    xx=basis@f.T;cc=cell@f.T
                    def compute(u,stress=False):
                        nonlocal calls
                        calls+=1;x=xx.copy();x[1]+=u
                        return force_features(x,cc,pbc,parameters,strain_derivative=stress,angle_beta=beta)
                    rr=root(lambda u:-compute(u)[1][1].sum(axis=1),np.zeros(3),options={'xtol':1e-10})
                    energy,force,work=compute(rr.x,True);maximum=float(np.max(abs(force.sum(axis=2))))
                    if maximum>2e-8:raise RuntimeError('internal force root unresolved: '+rr.message)
                    q=np.zeros(9);q[3]=e;q[6:]=rr.x
                    if pattern=='cubic_exx_plus_yz':q[0]=e
                    # Independent real-space Jet energy validates cell/shift counting.
                    ej=model.evaluate(q,derivatives=False);difference=abs(ej-energy.sum())
                    if difference>3e-10:raise RuntimeError('independent periodic energy differs')
                    stress=work.sum(axis=2)/abs(np.linalg.det(cc))*160.2176634
                    row=dict(case=name,pattern=pattern,increment=increment,strain=float(e),
                        shift_A=rr.x.tolist(),maximum_force_eV_A=maximum,stress_GPa=stress.tolist(),
                        energy_eV=float(energy.sum()),independent_Jet_energy_difference_eV=float(difference))
                    states.append(row);group.append(row)
                x=np.array([r['strain'] for r in group]);s=np.array([r['stress_GPa'] for r in group])
                slopes=np.einsum('i,iab->ab',x,s)/float(x@x)
                fits.append(dict(case=name,pattern=pattern,increment=increment,C11_GPa=float(slopes[0,0]) if pattern.startswith('cubic') else None,
                    C12_GPa=float((slopes[1,1]+slopes[2,2])/2) if pattern.startswith('cubic') else None,
                    C44_GPa=float(slopes[1,2]),exact_local_constants_GPa=exact.tolist()))
    dump(args.output/'states.json',states);dump(args.output/'fits.json',fits)
    dump(args.output/'summary.json',dict(complete=True,new_conservative_feature_calls=calls,states=len(states),fits=len(fits),
        published_SW_C44_GPa=published['diamond_c44'],old_environment_fully_reproduced=False,
        conclusion='separate finite-pattern diagnostic; differences cannot supply a universal correction of DFT constants',
        elapsed_seconds=time.perf_counter()-start,material_approved=False,first_initiation_validated=False,
        physical_clock_validated=False,new_DFT=0,new_MD=0))
    print(json.dumps(fits),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['selected','published','reference','historical_manifest','output']:parser.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    parser.add_argument('--candidate-only',action='store_true',help='evaluate only new selected candidate, without repeating prior SW control')
    main(parser.parse_args())
