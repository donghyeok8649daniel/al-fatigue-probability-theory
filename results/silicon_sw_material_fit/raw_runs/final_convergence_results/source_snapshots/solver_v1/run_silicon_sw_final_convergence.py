"""Resolution audit of the selected static Bessel candidate, not dynamics."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from .silicon_sw_anharmonic import AngularSW111BesselInterface
from .run_silicon_sw_material_fit import dump


def encode(j):
    return dict(energy=float(j.value),gradient=j.gradient.tolist(),hessian=j.hessian.tolist())


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter()
    selected=json.loads(args.selected.read_text());reference=json.loads(args.reference.read_text())
    roots=json.loads(args.branches.read_text());base=Path(__file__).resolve().parents[1]
    source=[]
    for relative in ['solver_v1/run_silicon_sw_final_convergence.py','solver_v1/silicon_sw_anharmonic.py',
                     'solver_v1/silicon_bessel_reference.py','solver_v1/silicon_environment_research.py',
                     'solver_v1/run_silicon_sw_material_fit.py']:
        raw=(base/relative).read_bytes();target=args.output/'source_snapshots'/relative
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        source.append(dict(path=relative,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',source)
    for role,path in [('selection',args.selected),('reference',args.reference),('branches',args.branches)]:
        target=args.output/'input_snapshots'/role/(path.name);target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(path.read_bytes())
    records=[];direct_calls=0
    # Change one resolution control at a time and include a simultaneous refinement.
    levels=[(48,512),(96,256),(96,512),(128,512),(128,768)]
    for cut in ['shuffle','glide']:
        opening=next(r['opening_A'] for r in roots if r['cut_kind']==cut and r['mesh_points']==161)
        states=[[opening,0.,0.],[.25,.15,-.07]]
        model=AngularSW111BesselInterface(selected['parameters'],reference['diamond_a0'],
                                         angle_beta=selected['angle_beta'],cut_kind=cut)
        exact=[model.evaluate(q,method='direct') for q in states];direct_calls+=len(states)
        for shells,nodes in levels:
            model=AngularSW111BesselInterface(selected['parameters'],reference['diamond_a0'],
                angle_beta=selected['angle_beta'],cut_kind=cut,shell_index=shells,nodes=nodes)
            for q,d in zip(states,exact):
                b=model.evaluate(q,method='bessel')
                records.append(dict(cut=cut,state=q,shell_index=shells,nodes=nodes,bessel=encode(b),direct=encode(d),
                    energy_error_eV=abs(float(b.value-d.value)),gradient_error_eV_A=float(np.max(abs(b.gradient-d.gradient))),
                    hessian_error_eV_A2=float(np.max(abs(b.hessian-d.hessian)))))
                dump(args.output/'resolution_checks.json',records)
            print(json.dumps(dict(cut=cut,shell_index=shells,nodes=nodes,completed_states=len(records))),flush=True)
    fine=[r for r in records if r['shell_index']==128 and r['nodes']==768]
    production=[r for r in records if r['shell_index']==96 and r['nodes']==512]
    dump(args.output/'summary.json',dict(complete=True,new_bessel_jets=len(records),new_direct_jets=direct_calls,
        final_refinement_maximum_errors={key:max(r[key] for r in fine) for key in ['energy_error_eV','gradient_error_eV_A','hessian_error_eV_A2']},
        nominal_resolution_maximum_errors={key:max(r[key] for r in production) for key in ['energy_error_eV','gradient_error_eV_A','hessian_error_eV_A2']},
        elapsed_seconds=time.perf_counter()-start,certified_uniform_tail_bound=False,
        statement='sampled rigid static states only; no first-crack or kinetic validation',
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,new_DFT=0,new_MD=0))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['selected','reference','branches','output']:parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
