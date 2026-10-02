"""Reproduce the static Fourier--Bessel SW audit in a NEW output directory."""
from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np
from scipy.optimize import brentq
from results.silicon_wafer_feasibility.run_static_probe import source_parameters, pair
from .silicon_bessel_reference import (SW111BesselInterface, SWPlaneBessel, site_energy,
                                       PlaneMoments, _point_moments)
from .silicon_environment_research import DiamondCell, Jet, _site_energy


def array_record(jet):
    return dict(energy_eV=float(jet.value), gradient_eV_per_A=jet.gradient.tolist(),
                hessian_eV_per_A2=jet.hessian.tolist())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists():
        raise ValueError('preserve prior results; choose a new output directory')
    args.output.mkdir(parents=True)
    start=time.perf_counter()
    p,parameter_hash=source_parameters()
    bond=brentq(lambda r: float(pair(r,p,derivative=True)),2.2,2.5,xtol=1e-14)
    lattice=4*bond/np.sqrt(3)
    states=[[0.,0.,0.],[-.1,.12,.07],[0.,.35,.19],[.21,.15,-.12],
            [1.1,.32,.21],[5.,.13,-.07]]
    rows=[]
    for kind in ['shuffle','glide']:
        model=SW111BesselInterface(p,lattice,cut_kind=kind)
        for q in states:
            a,b=[model.evaluate(q,method=m) for m in ['bessel','direct']]
            error=dict(energy_eV=abs(float(a.value-b.value)),
                       gradient_eV_per_A=float(np.max(abs(a.gradient-b.gradient))),
                       hessian_eV_per_A2=float(np.max(abs(a.hessian-b.hessian))))
            row=dict(cut=kind,state_opening_x_y_A=q,bessel=array_record(a),
                     direct=array_record(b),absolute_error=error)
            rows.append(row)
            print(json.dumps(dict(cut=kind,state=q,error=error)),flush=True)
        planes=model._planes(0.)
        bulk=0.
        for basis in [0,1]:
            center=next(x for x in planes if x[:2] == (0,basis))
            bulk += site_energy(model._center_moments(center,planes,np.zeros(3),False),p).value
        reference_bulk=DiamondCell(p,lattice).evaluate(np.zeros(9),'direct').value
        assert abs(bulk-reference_bulk) < 1e-8
    (args.output/'states.json').write_text(json.dumps(rows,indent=2)+'\n')
    # A witness that independently squaring each plane misses cross-plane angles.
    vectors=np.array([[.2,-.3,2.4],[1.8,.4,1.5]])
    fields=[_point_moments(v,p) for v in vectors]
    def energy(value):
        return site_energy(PlaneMoments(value,np.zeros((15,3)),np.zeros((15,3,3))),p).value
    combined=energy(sum(fields))
    explicit=_site_energy(vectors,p,'direct')
    incorrectly_separate=sum(energy(v) for v in fields)
    assert abs(combined-explicit) < 3e-12
    assert abs(incorrectly_separate-explicit) > 1e-4
    (args.output/'angular_cross_check.json').write_text(json.dumps(dict(
        combined_moments_energy_eV=float(combined),explicit_three_body_energy_eV=float(explicit),
        incorrect_plane_by_plane_energy_eV=float(incorrectly_separate),
        missing_cross_term_eV=float(explicit-incorrectly_separate)),indent=2)+'\n')
    # Separate reciprocal and radial-quadrature refinement at a nonsymmetric point.
    refinement=[];q=np.array([.79,1.7,.9])
    for shell,nodes in [(16,1024),(32,1024),(48,1024),(64,1024),(96,1024),(96,256),(96,512)]:
        plane=SWPlaneBessel(p,lattice,shell_index=shell,nodes=nodes)
        a=site_energy(plane.evaluate(q[0],q[1:]),p)
        x=[Jet.variable(v,3,i) for i,v in enumerate(q)]
        vectors=[]
        for v in plane.direct_vectors(q[0],q[1:]):
            base=v-q[[1,2,0]]
            vectors.append([base[0]+x[1],base[1]+x[2],base[2]+x[0]])
        b=_site_energy(vectors,p,'direct')
        row=dict(shell_index=shell,quadrature_nodes=nodes,
                 bessel=array_record(a),direct=array_record(b),
                 energy_error_eV=abs(float(a.value-b.value)),
                 gradient_error_eV_per_A=float(np.max(abs(a.gradient-b.gradient))),
                 hessian_error_eV_per_A2=float(np.max(abs(a.hessian-b.hessian))))
        refinement.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ('bessel','direct')}),flush=True)
    (args.output/'refinement.json').write_text(json.dumps(refinement,indent=2)+'\n')
    summary=dict(scope='unchanged pure-Si original SW, static rigid periodic (111) interface',
                 lattice_A=lattice,bond_A=bond,parameters=p,parameter_sha256=parameter_hash,
                 reciprocal_shell_index=96,quadrature_nodes=512,interface_states=len(rows),
                 max_energy_error_eV=max(r['absolute_error']['energy_eV'] for r in rows),
                 max_gradient_error_eV_per_A=max(r['absolute_error']['gradient_eV_per_A'] for r in rows),
                 max_hessian_error_eV_per_A2=max(r['absolute_error']['hessian_eV_per_A2'] for r in rows),
                 elapsed_seconds=time.perf_counter()-start,
                 new_DFT=0,new_MD=0,material_fit=False,first_initiation_validated=False,
                 probability_evaluations=0,physical_clock_validated=False,
                 refinement_is_certified_tail_bound=False)
    assert summary['max_energy_error_eV'] < 2e-8
    assert summary['max_gradient_error_eV_per_A'] < 2e-7
    assert summary['max_hessian_error_eV_per_A2'] < 3e-6
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    root=Path(__file__).resolve().parent.parent
    sources=['solver_v1/silicon_bessel_reference.py','solver_v1/run_silicon_bessel_reference.py',
             'solver_v1/test_silicon_bessel_reference.py','solver_v1/silicon_environment_research.py',
             'solver_v1/fcc111_geometry.py','solver_v1/fcc111_lattice_sum.py',
             'results/silicon_wafer_feasibility/run_static_probe.py',
             'results/silicon_wafer_feasibility/source_Si.sw']
    manifest=[dict(path=f,sha256=hashlib.sha256((root/f).read_bytes()).hexdigest()) for f in sources]
    (args.output/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__ == '__main__':
    main()
