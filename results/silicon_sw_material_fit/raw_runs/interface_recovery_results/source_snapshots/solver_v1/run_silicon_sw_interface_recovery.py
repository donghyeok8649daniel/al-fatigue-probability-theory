"""Correct stress display units from recorded dimensional jets, without reruns.

The first interface audit used the eV/A^2 surface-energy factor for traction
in eV/A^3. Its GPa display was therefore ten times too small. Raw energy,
gradient, Hessian, roots, separation work and FD comparisons remain valid.
Preserve that audit; create a new authoritative report with SI-derived units.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from .run_silicon_sw_material_fit import dump


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);read=lambda name:json.loads((args.run/name).read_text())
    summary=read('summary.json')
    if not summary['complete']:raise RuntimeError('raw branch audit incomplete')
    bindings=read('source_manifest.json')
    for record in bindings:
        raw=(args.run/'source_snapshots'/record['path']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=record['sha256']:raise RuntimeError('raw calculation source differs')
    roots=read('branch_roots.json');converted=[];area={};conversion=1.602176634e-19/(1e-10)**3/1e9
    for row in roots:
        raw_traction=row['jet']['gradient_eV_A'][0]/row['atomic_cell_area_A2']
        if abs(row['normal_traction_GPa']-raw_traction*16.02176634)>1e-10:
            raise RuntimeError('raw report does not have the identified display conversion')
        row=dict(row,normal_traction_GPa=raw_traction*conversion)
        converted.append(row);area[(row['case'],row['cut_kind'])]=row['atomic_cell_area_A2']
    with (args.run/'normal_branch.csv').open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
    for row in rows:
        row['normal_traction_GPa']=float(row['normal_traction_GPa'])*10.
    with (args.output/'normal_branch.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fd=read('finite_difference_checks.json');maximum_fd=max(r['maximum_hessian_error_eV_A2'] for s in fd for r in s['checks'])
    if maximum_fd>3e-6:raise RuntimeError('independent force-Hessian finite difference failed')
    bs=read('bessel_audit.json')
    for row in bs:
        b,d=row['bessel'],row['direct']
        if abs(b['energy_eV']-d['energy_eV'])>2e-8 or np.max(abs(np.array(b['gradient_eV_A'])-d['gradient_eV_A']))>2e-7:
            raise RuntimeError('raw Bessel force-energy audit failed')
        if np.max(abs(np.array(b['hessian_eV_A2'])-d['hessian_eV_A2']))>3e-6:raise RuntimeError('raw Bessel Hessian audit failed')
    root_differences=[]
    for key in area:
        pair=[r for r in converted if (r['case'],r['cut_kind'])==key]
        if len(pair)!=2 or sorted(r['mesh_points'] for r in pair)!=[81,161]:raise RuntimeError('root mesh audit incomplete')
        root_differences.append(abs(pair[0]['opening_A']-pair[1]['opening_A']))
    dump(args.output/'branch_roots.json',converted)
    dump(args.output/'summary.json',dict(complete=True,corrected_quantity='normal_traction_GPa only',
        unit_conversion_eV_A3_to_GPa=conversion,previous_unit_error='used surface-energy eV/A2-to-J/m2 factor for traction',
        raw_run_complete=True,raw_sources_verified=True,maximum_Hessian_FD_error_eV_A2=maximum_fd,
        maximum_root_mesh_difference_A=max(root_differences),new_geometry_evaluations=0,
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False))
    path='solver_v1/run_silicon_sw_interface_recovery.py';raw=Path(__file__).read_bytes();target=args.output/'source_snapshots'/path
    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    dump(args.output/'source_manifest.json',[dict(path=path,sha256=hashlib.sha256(raw).hexdigest())])
    print(json.dumps(dict(maximum_FD=maximum_fd,max_root_difference_A=max(root_differences),conversion=conversion)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['run','output']:parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
