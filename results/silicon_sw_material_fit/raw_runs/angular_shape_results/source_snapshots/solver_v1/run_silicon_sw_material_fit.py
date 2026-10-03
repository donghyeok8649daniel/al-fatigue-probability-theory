"""Fit a minimal SW candidate to source DFT, then test whole held-out families.

No DFT, MD, MACE, probability clock or production parameter update is run.
Requires a fresh output folder and the previously hash-audited Cambridge data.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize_scalar
from results.silicon_wafer_feasibility.run_dft_material_audit import (
    read_frames,reference_fields,ARCHIVE_SHA256,DATA_SHA256)
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_environment_research import DiamondCell
from .silicon_sw_material_fit import force_features,fit_amplitudes,amplitude_parameters
from .silicon_bessel_reference import SW111BesselInterface


TRAIN_TYPES=('dia','surface_001','surface_110')
FORCE_SCALE=.2  # eV/A, declared loss scale, not a source uncertainty or CI.
ENERGY_SCALE=.05  # eV/atom, declared loss scale, not a source uncertainty.


def dump(path,data):
    path.write_bytes((json.dumps(data,indent=2,allow_nan=False)+'\n').encode('utf-8'))


def bulk_diagnostic(p):
    def energy(lattice):
        return DiamondCell(p,lattice).evaluate(np.zeros(9),'direct',derivatives=False)/2
    opt=minimize_scalar(energy,bounds=(4.5,6.5),method='bounded',options={'xatol':1e-11})
    if not opt.success or abs(opt.x-4.5)<1e-4 or abs(opt.x-6.5)<1e-4:
        raise ValueError('diamond equilibrium not resolved within declared lattice interval')
    model=DiamondCell(p,opt.x)
    value=model.evaluate(np.zeros(9),'direct')
    h=value.hessian
    internal=np.linalg.eigvalsh(h[6:,6:])
    if np.min(internal)<=0:
        raise ValueError('diamond internal-displacement Hessian unstable')
    relaxed=h[:6,:6]-h[:6,6:]@np.linalg.solve(h[6:,6:],h[6:,:6])
    scale=160.2176634/model.volume
    return dict(lattice_A=float(opt.x),energy_eV_atom=float(opt.fun),
        gradient=value.gradient.tolist(),full_hessian=h.tolist(),
        internal_eigenvalues=internal.tolist(),relaxed_elastic_hessian=relaxed.tolist(),
        relaxed_elastic_eigenvalues=np.linalg.eigvalsh(relaxed).tolist(),
        C11_GPa=float(relaxed[0,0]*scale),C12_GPa=float(relaxed[0,1]*scale),
        C44_GPa=float(relaxed[3,3]*scale),
        scope='0 K diamond diagnostic; RT experimental elasticity is a different target')


def main(args):
    if args.output.exists():
        raise ValueError('fresh output directory required; preserve all earlier data')
    args.output.mkdir(parents=True)
    start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    p,parameter_hash=source_parameters()
    protocol=dict(train='all explicit-PW91 dia, surface_001 and surface_110 source frames',
        held_out='every other config-type/XC combination; entire Si(111) families excluded from fit',
        loss='force and same-XC diamond-referenced energy; equal total weight per training family',
        force_scale_eV_A=FORCE_SCALE,energy_scale_eV_atom=ENERGY_SCALE,
        scales_are_uncertainty=False,amplitudes=['repulsion','attraction','centered_angle'],
        fixed=['epsilon','sigma','a','gamma','costheta0','p','q'],
        epsilon_gauge='epsilon fixed; candidate maps to A, B, lambda',
        shape_search=False,regularization=False,
        heldout_used_for_optimization=False,source_dataset='GAP fitting database, no GAP fitted here',
        archive_sha256=ARCHIVE_SHA256,member_sha256=DATA_SHA256,
        original_SW_parameter_sha256=parameter_hash,new_DFT=0,new_MD=0,new_MACE=0,
        first_initiation_validated=False,physical_clock_validated=False,production_enabled=False)
    dump(args.output/'protocol.json',protocol)
    frames,_=read_frames(args.archive,download=False)
    refs=[reference_fields(a) for a in frames]
    training=[i for i,a in enumerate(frames) if a.info.get('xc_functional')=='PW91'
              and a.info['config_type'] in TRAIN_TYPES]
    counts=Counter(str(frames[i].info['config_type']) for i in training)
    if counts != Counter(dia=489,surface_001=10,surface_110=12):
        raise ValueError('declared source training selection changed')
    anchor=min([i for i in training if frames[i].info['config_type']=='dia'],
               key=lambda i:refs[i][0]/len(frames[i]))
    # Existing independent LAMMPS predictions are controls, not new engine calls.
    with np.load(root/'results/silicon_atomistic_v5/dft_predictions.npz') as saved:
        old_force=saved['original_sw_1985_force'];old_energy=saved['original_sw_1985_energy']
        old_offsets=saved['offsets'];old_ref=saved['reference_force']
    offsets=np.r_[0,np.cumsum([len(a) for a in frames])]
    assert np.array_equal(offsets,old_offsets)
    reference_force=np.vstack([f for _,f,_ in refs])
    assert np.array_equal(reference_force,old_ref)
    energies=np.empty((len(frames),3));forces=np.empty((offsets[-1],3,3))
    max_e=0.;max_f=0.
    for i,a in enumerate(frames):
        e,f=force_features(a.positions,a.cell.array,a.pbc,p)
        lo,hi=offsets[i:i+2];energies[i]=e;forces[lo:hi]=f
        max_e=max(max_e,abs(float(e.sum()-old_energy[i])))
        max_f=max(max_f,float(np.max(abs(f.sum(axis=2)-old_force[lo:hi]))))
        if max_e > 1e-7 or max_f > 1e-7:
            dump(args.output/'baseline_failure.json',dict(frame=i,max_energy_error_eV=max_e,
                                                         max_force_error_eV_A=max_f))
            raise ValueError('new force basis disagrees with archived independent LAMMPS control')
        if (i+1)%250==0:
            progress=dict(completed=i+1,planned=len(frames),elapsed_seconds=time.perf_counter()-start,
                          max_control_energy_error_eV=max_e,max_control_force_error_eV_A=max_f)
            dump(args.output/'running.json',progress);print(json.dumps(progress),flush=True)
    np.savez_compressed(args.output/'features.npz',offsets=offsets,energy_basis=energies,
        force_basis=forces,reference_energy=np.array([r[0] for r in refs]),reference_force=reference_force)
    design_force=[];target_force=[];design_energy=[];target_energy=[]
    ea=energies[anchor]/len(frames[anchor]);ra=refs[anchor][0]/len(frames[anchor])
    for i in training:
        a=frames[i];n=len(a);lo,hi=offsets[i:i+2];count=counts[str(a.info['config_type'])]
        design_force.append(forces[lo:hi].reshape(-1,3)/(FORCE_SCALE*np.sqrt(3*n*count)))
        target_force.append(refs[i][1].ravel()/(FORCE_SCALE*np.sqrt(3*n*count)))
        design_energy.append((energies[i]/n-ea)/(ENERGY_SCALE*np.sqrt(count)))
        target_energy.append((refs[i][0]/n-ra)/(ENERGY_SCALE*np.sqrt(count)))
    xf,yf=np.vstack(design_force),np.concatenate(target_force)
    xe,ye=np.array(design_energy),np.array(target_energy)
    amplitudes,fit=fit_amplitudes(np.r_[xf,xe],np.r_[yf,ye])
    diagnostic_fits={}
    for name,x,y in [('force_only',xf,yf),('energy_only',xe,ye)]:
        _,diagnostic_fits[name]=fit_amplitudes(x,y)
    dump(args.output/'fit.json',dict(primary=fit,ablation_only=diagnostic_fits,
         training_frames=training,training_type_counts=dict(counts),anchor_frame=anchor,
         anchor_is_relaxed_ground_state=False))
    candidate=amplitude_parameters(p,amplitudes)
    dump(args.output/'candidate_parameters.json',dict(parameters=candidate,material_approved=False,
         parameter_source='new three-amplitude fit; epsilon gauge fixed',amplitudes=amplitudes.tolist()))
    metrics=[]
    for i,a in enumerate(frames):
        lo,hi=offsets[i:i+2];truth=reference_force[lo:hi]
        before=forces[lo:hi].sum(axis=2);after=np.einsum('nck,k->nc',forces[lo:hi],amplitudes)
        xc=str(a.info.get('xc_functional','UNSPECIFIED'))
        row=dict(frame=i,config_type=str(a.info['config_type']),declared_xc=xc,atoms=len(a),
                 training=i in training,baseline_force_RMSE_eV_A=float(np.sqrt(np.mean((before-truth)**2))),
                 candidate_force_RMSE_eV_A=float(np.sqrt(np.mean((after-truth)**2))),
                 baseline_energy_eV=float(energies[i].sum()),candidate_energy_eV=float(energies[i]@amplitudes),
                 reference_energy_eV=refs[i][0],baseline_relative_energy_error_eV_atom=None,
                 candidate_relative_energy_error_eV_atom=None)
        if xc=='PW91':
            delta=energies[i]/len(a)-ea;rd=refs[i][0]/len(a)-ra
            row.update(baseline_relative_energy_error_eV_atom=float(delta.sum()-rd),
                       candidate_relative_energy_error_eV_atom=float(delta@amplitudes-rd))
        metrics.append(row)
    with (args.output/'frame_metrics.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(metrics[0]));writer.writeheader();writer.writerows(metrics)
    groups=defaultdict(list)
    for row in metrics:
        groups[(row['config_type'],row['declared_xc'],row['training'])].append(row)
    group_results=[]
    for (kind,xc,trained),rows in sorted(groups.items()):
        size=sum(r['atoms'] for r in rows)
        result=dict(config_type=kind,declared_xc=xc,training=trained,frames=len(rows),atoms=size)
        for model in ['baseline','candidate']:
            result[model+'_force_RMSE_eV_A']=float(np.sqrt(sum(r['atoms']*r[model+'_force_RMSE_eV_A']**2 for r in rows)/size))
            errors=[r[model+'_relative_energy_error_eV_atom'] for r in rows if r[model+'_relative_energy_error_eV_atom'] is not None]
            result[model+'_relative_energy_RMSE_eV_atom']=float(np.sqrt(np.mean(np.square(errors)))) if errors else None
        group_results.append(result)
    dump(args.output/'group_metrics.json',group_results)
    bulk={name:bulk_diagnostic(par) for name,par in [('original',p),('candidate',candidate)]}
    dump(args.output/'bulk_diagnostic.json',bulk)
    audit=[]
    for kind in ['shuffle','glide']:
        model=SW111BesselInterface(candidate,bulk['candidate']['lattice_A'],cut_kind=kind)
        for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
            a,b=[model.evaluate(state,method=m) for m in ['bessel','direct']]
            row=dict(cut=kind,state=state,
                bessel=dict(energy=float(a.value),gradient=a.gradient.tolist(),hessian=a.hessian.tolist()),
                direct=dict(energy=float(b.value),gradient=b.gradient.tolist(),hessian=b.hessian.tolist()),
                energy_error_eV=abs(float(a.value-b.value)),
                gradient_error_eV_A=float(np.max(abs(a.gradient-b.gradient))),
                hessian_error_eV_A2=float(np.max(abs(a.hessian-b.hessian))))
            audit.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ('bessel','direct')}),flush=True)
    dump(args.output/'bessel_candidate_audit.json',audit)
    bessel_pass=(max(r['energy_error_eV'] for r in audit)<2e-8
                 and max(r['gradient_error_eV_A'] for r in audit)<2e-7
                 and max(r['hessian_error_eV_A2'] for r in audit)<3e-6)
    summary=dict(complete=True,training_frames=len(training),held_out_frames=len(frames)-len(training),
        new_feature_evaluations=len(frames),new_DFT=0,new_MD=0,new_LAMMPS=0,new_MACE=0,
        baseline_control_max_energy_error_eV=max_e,baseline_control_max_force_error_eV_A=max_f,
        fit=fit,candidate_bessel_consistency_passed=bessel_pass,
        diamond_internal_stable=min(bulk['candidate']['internal_eigenvalues'])>0,
        diamond_elastic_stable=min(bulk['candidate']['relaxed_elastic_eigenvalues'])>0,
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        production_enabled=False,elapsed_seconds=time.perf_counter()-start)
    dump(args.output/'summary.json',summary)
    sources=['solver_v1/silicon_sw_material_fit.py','solver_v1/run_silicon_sw_material_fit.py',
             'solver_v1/test_silicon_sw_material_fit.py','solver_v1/silicon_bessel_reference.py',
             'solver_v1/silicon_environment_research.py',
             'results/silicon_wafer_feasibility/run_dft_material_audit.py',
             'results/silicon_wafer_feasibility/run_static_probe.py',
             'results/silicon_wafer_feasibility/source_Si.sw',
             'results/silicon_atomistic_v5/dft_predictions.npz']
    dump(args.output/'source_manifest.json',[dict(path=path,sha256=hashlib.sha256((root/path).read_bytes()).hexdigest(),
        lf_sha256=hashlib.sha256((root/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()) for path in sources])
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
