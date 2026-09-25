"""Bounded, unrelaxed energy/force probes along already computed soft directions.

The directions are Cartesian Hessian eigenvectors, not mass-weighted vibrations.
No path optimization, basin integration, barrier or transition rate is computed.
Thermal amplitudes refer only to a labelled classical quadratic coordinate.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann, electron_volt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if args.output.exists():
        raise ValueError('fresh anharmonic output required')
    if hashlib.sha256(args.model.read_bytes()).hexdigest() != MODEL_SHA256:
        raise ValueError('baseline model SHA differs')
    if args.max_seconds <= 0:
        raise ValueError('positive wall-time budget required')
    sources=[]
    for name in ('loading8', 'return8'):
        folder=args.results/('dense_'+name)
        summary=json.loads((folder/'summary.json').read_text(encoding='utf-8'))
        protocol=json.loads((folder/'protocol.json').read_text(encoding='utf-8'))
        if not summary['complete'] or protocol['ensemble'] != 'displacement':
            raise ValueError('complete fixed-grip Hessian required')
        with np.load(folder/'raw_hessian.npz') as d:
            source={key:d[key].copy() for key in d.files}
        with np.load(folder/'spectrum.npz') as d:
            source.update({key:d[key].copy() for key in d.files})
        if source['eigenvalues'][0] <= 0:
            raise ValueError('positive local curvature required for thermal amplitude')
        source.update(name=name, raw_sha256=hashlib.sha256((folder/'raw_hessian.npz').read_bytes()).hexdigest(),
            spectrum_sha256=hashlib.sha256((folder/'spectrum.npz').read_bytes()).hexdigest())
        sources.append(source)
    args.output.mkdir(parents=True)
    started=time.perf_counter();calls=0;rows=[];skipped=[];records=[]
    kb=Boltzmann/electron_volt;kt=300.*kb
    dump(args.output/'protocol.json', dict(model_sha256=MODEL_SHA256,
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=[dict(state=s['name'], raw_sha256=s['raw_sha256'], spectrum_sha256=s['spectrum_sha256']) for s in sources],
        modes=[0,1,2],fixed_max_atom_displacements_A=[.01,.05,.1,.2],
        quadratic_coordinate_amplitudes_at_300K=[1.,2.],maximum_atom_bound_A=.3,
        bound_policy='record and skip thermal probes above bound; do not clip or substitute',
        force_reference='entire constrained gradient compared with g0 + H d',
        max_seconds=args.max_seconds,new_DFT=0,new_MD=0,relaxation=False,
        interpretation='six straight directions from two approximate stationary states; not paths or barriers'))
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    complete=True;interruption=None
    try:
        for source in sources:
            name=source['name'];r0=source['positions'];basis=source['basis'];g0=source['gradient']
            h=(source['hessian']+source['hessian'].T)/2
            atoms=Atoms(numbers=source['numbers'],positions=r0,cell=source['cell'],pbc=False)
            base_energy,base_forces=force_only(calc,atoms);calls+=1
            replay=dict(state=name,energy_error_eV=base_energy-float(source['energy']),
                force_max_error_eV_A=float(np.max(abs(base_forces-source['forces']))))
            if abs(replay['energy_error_eV'])>1e-8 or replay['force_max_error_eV_A']>1e-8:
                raise ValueError('source energy or forces changed')
            records.append(replay)
            evaluated=[];energies=[];forces=[];directions=[]
            for mode in range(3):
                v=source['eigenvectors'][:,mode].copy()
                if v[np.argmax(abs(v))]<0:v=-v
                direction=(basis@v).reshape(r0.shape);directions.append(direction)
                curvature=float(source['eigenvalues'][mode]);linear=float(g0@v)
                peak=float(np.linalg.norm(direction,axis=1).max())
                amplitudes=[('max_atom_'+str(a),a/peak) for a in (.01,.05,.1,.2)]
                amplitudes += [('quadratic_sigma_'+str(a),a*np.sqrt(kt/curvature)) for a in (1.,2.)]
                for label,q in amplitudes:
                    if q*peak>.3+1e-15:
                        skipped.append(dict(state=name,mode=mode,label=label,coordinate_A=float(q),
                            maximum_atom_displacement_A=float(q*peak),reason='predeclared displacement bound'))
                        continue
                    pair=[]
                    for sign in (-1.,1.):
                        if time.perf_counter()-started>args.max_seconds:
                            raise TimeoutError('bounded nonlinear-probe budget; completed points retained')
                        x=sign*q;atoms.positions[:]=r0+x*direction
                        e,f=force_only(calc,atoms);calls+=1
                        g=-basis.T@f.ravel();delta=e-base_energy
                        predicted=linear*x+.5*curvature*x*x
                        expected_g=g0+x*(h@v)
                        row=dict(state=name,mode=mode,label=label,sign=int(sign),coordinate_A=float(x),
                            maximum_atom_displacement_A=float(abs(x)*peak),curvature_eV_A2=curvature,
                            baseline_directional_gradient_eV_A=linear,actual_delta_energy_eV=float(delta),
                            linear_quadratic_delta_energy_eV=float(predicted),
                            nonlinear_energy_residual_eV=float(delta-predicted),
                            nonlinear_energy_residual_over_kBT300=float((delta-predicted)/kt),
                            actual_directional_gradient_eV_A=float(g@v),
                            nonlinear_full_gradient_norm_eV_A=float(np.linalg.norm(g-expected_g)),
                            nonlinear_directional_gradient_eV_A=float((g-expected_g)@v))
                        rows.append(row);pair.append(row);evaluated.append((mode,x));energies.append(e);forces.append(f)
                        dump(args.output/'points.json', rows)
                        np.savez_compressed(args.output/(name+'_raw.npz'),positions=r0,numbers=source['numbers'],
                            cell=source['cell'],basis=basis,directions=np.asarray(directions),
                            evaluated=np.asarray(evaluated),energies=np.asarray(energies),forces=np.asarray(forces),
                            baseline_energy=base_energy,baseline_forces=base_forces)
                    for row in pair:
                        row['symmetric_nonlinear_energy_eV']=(pair[0]['actual_delta_energy_eV']+pair[1]['actual_delta_energy_eV'])/2-.5*curvature*q*q
                        row['antisymmetric_nonlinear_energy_eV']=(pair[1]['actual_delta_energy_eV']-pair[0]['actual_delta_energy_eV'])/2-linear*q
                    dump(args.output/'points.json', rows)
    except TimeoutError as exc:
        complete=False;interruption=str(exc)
    result=dict(complete=complete,interruption=interruption,actual_points=len(rows),new_model_calls=calls,
        completed_source_replays=records,skipped=skipped,elapsed_s=time.perf_counter()-started,
        maximum_absolute_nonlinear_energy_eV=max((abs(r['nonlinear_energy_residual_eV']) for r in rows),default=None),
        new_DFT=0,new_MD=0,relaxation=False,actual_temperature_sampling=False,
        first_crack_certified=False,physical_clock=None,material_approved=False,
        interpretation='unrelaxed line-probe failures can reject local quadratic extrapolation; agreement cannot certify a basin or finite-T free energy')
    dump(args.output/'summary.json',result)
    if rows:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
        for ax,name in zip(axes,('loading8','return8')):
            for mode in range(3):
                points=sorted((r for r in rows if r['state']==name and r['mode']==mode and r['label'].startswith('max_atom')),key=lambda r:r['coordinate_A'])
                ax.plot([r['sign']*r['maximum_atom_displacement_A'] for r in points],
                    [r['nonlinear_energy_residual_eV'] for r in points],'.-',label='soft direction '+str(mode+1))
            ax.axhline(0,color='#555',lw=.6);ax.set(title=name,xlabel='Signed maximum atom displacement (Angstrom)',ylabel='Energy minus linear/quadratic prediction (eV)');ax.legend(fontsize=8)
        fig.suptitle('Actual unrelaxed energy probes along three soft directions\nNeither transition paths nor crack barriers',fontsize=11)
        fig.savefig(args.output/'anharmonic_line_probes.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','model','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-seconds',type=float,default=600.)
    main(p.parse_args())
