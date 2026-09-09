"""Independent verification and reporting of ACTUALLY saved static row states.

No fitting or replacement of a stopped relaxation by its desired endpoint.
"""
from __future__ import annotations

import argparse
import csv
import json
import time

import numpy as np

from .nonlinear_fcc_screw import NonlinearScrewRows
from .run_low_stress_cyclic_diagnostic import build_surface, save_json, write_csv
from .run_vector_core_study import OUT, read_state
from .vector_fcc_rows import VectorPeriodicCell, VectorRowKernel
from .vector_fcc_validation import direct_row_channels


def validation():
    surface,_,metadata=build_surface(tolerance=2e-11);rows=NonlinearScrewRows(surface,tolerance=2e-12)
    kernel=VectorRowKernel(rows);records=[]
    # Fixed non-random points, both registry signs, small/large row radii.
    vectors=np.array([[.071+.173*i,.67+.071*(i%4),(-1)**i*(.11+.21*(i%3))] for i in range(8)]
                     +[[.271,1.91,.39],[-.339,.15,2.17],[.419,-3.21,.17],[-.113,4.32,-1.87]])
    tick=time.perf_counter();reciprocal=kernel.evaluate(vectors,order=2)
    kernel_seconds=time.perf_counter()-tick
    for images in (96,192,384):
        reference=[direct_row_channels(rows,v,images=images) for v in vectors]
        for field in ('value','gradient','hessian'):
            exact=np.array([r[field] for r in reference]);delta=abs(reciprocal[field]-exact)
            resolved=abs(exact)>np.sqrt(np.finfo(float).eps)*max(np.max(abs(exact)),1.)
            records.append(dict(images_each_direction=images,quantity=field,points=len(vectors),
                max_absolute_error=float(np.max(delta)),max_scaled_error=float(np.max(delta/(1+abs(exact)))),
                max_relative_error_away_from_zero=float(np.max(delta[resolved]/abs(exact[resolved]))),
                relative_display_exclusion='sqrt(machine_epsilon)*max(1,max_reference); not a physical floor'))
    write_csv(OUT/'direct_atomic_validation.csv',records)
    v=np.column_stack([np.sin(np.arange(len(vectors))*.61),np.cos(np.arange(len(vectors))*.39),
                       np.sin(np.arange(len(vectors))*.31+.7)])
    derivatives=[]
    for step in (4e-6,2e-6,1e-6):
        plus=kernel.evaluate(vectors+step*v);minus=kernel.evaluate(vectors-step*v)
        expected=np.einsum('nci,ni->nc',reciprocal['gradient'],v)
        curvature=np.einsum('ncij,nj->nci',reciprocal['hessian'],v)
        derivatives.append(dict(step=step,
            gradient_fd_error=float(np.max(abs((plus['value']-minus['value'])/(2*step)-expected))),
            hessian_fd_error=float(np.max(abs((plus['gradient']-minus['gradient'])/(2*step)-curvature)))))
    write_csv(OUT/'derivative_validation.csv',derivatives)
    save_json(OUT/'verification.json',dict(parameter_sha256=metadata['parameter_sha256'],
        kernel_points=len(vectors),kernel_evaluation_seconds=kernel_seconds,
        reciprocal_modes=reciprocal['modes_used'],last_mode_envelope=reciprocal['maximum_last_mode_envelope'],
        worst_direct_absolute_error=max(r['max_absolute_error'] for r in records if r['images_each_direction']==384),
        direct_reference_is_validation_only=True))
    print('independent vector row validation saved',flush=True)


def tail_study(relative_state):
    state_path=OUT/relative_state
    field,gamma=read_state(state_path)
    surface,_,metadata=build_surface(tolerance=2e-11)
    bound=json.loads((state_path.parent/'metadata.json').read_text(encoding='utf8'))
    if bound['parameter_sha256']!=metadata['parameter_sha256']:
        raise ValueError('state/parameter binding mismatch')
    rows=NonlinearScrewRows(surface,tolerance=2e-12);records=[];results=[]
    for ring in (4,6,8,10,12,16):
        tick=time.perf_counter();cell=VectorPeriodicCell(rows,field.shape[:2],ring=ring)
        out=cell.evaluate(field,gamma)
        records.append(dict(state=relative_state,ring=ring,energy_ev=out['energy'],
            max_force=float(np.max(abs(out['gradient']))),reciprocal_modes=out['reciprocal_modes'],
            last_mode_envelope=out['last_mode_envelope'],
            **cell.pair_zero_mode_tail_bounds(field),wall_seconds=time.perf_counter()-tick))
        results.append(out)
        print(f"tail ring={ring}, E={out['energy']:.10g}, seconds={records[-1]['wall_seconds']:.2f}",flush=True)
    for i,(record,result) in enumerate(zip(records,results)):
        record['energy_error_vs_ring16']=abs(result['energy']-results[-1]['energy'])
        record['force_error_vs_ring16']=float(np.max(abs(result['gradient']-results[-1]['gradient'])))
        record['energy_change_from_previous']=None if not i else abs(result['energy']-results[i-1]['energy'])
        record['force_change_from_previous']=None if not i else float(np.max(abs(result['gradient']-results[i-1]['gradient'])))
    filename=('nonzero_transverse_tail_refinement.csv' if relative_state=='dipole_n32_r8/iteration_0_20.csv.gz'
              else f"{state_path.parent.name}_{state_path.name.removesuffix('.csv.gz')}_tail.csv")
    write_csv(OUT/filename,records)


def source_history(directory):
    """Follow explicit continuation without overwriting the interrupted record."""
    meta=json.loads((directory/'metadata.json').read_text())
    root_initial=json.loads((directory/'initial_state.json').read_text());prefix=[]
    source=meta.get('continuation_from_state')
    if source:
        parent=OUT/source
        old,root_initial=source_history(parent.parent)
        source_meta=json.loads((parent.parent/'metadata.json').read_text())
        if parent.name.startswith('state_') and source_meta['transverse_ring']==meta['transverse_ring']:
            step=int(parent.name.split('_')[1].split('.')[0])
            prefix=[r for r in old if int(r['step'])<=step]
    path=directory/'stress_history.csv'
    if path.exists():
        with path.open() as stream:
            current=list(csv.DictReader(stream))
    else:
        current=[]
    offset=len(prefix)
    for r in current:
        r['state_source']=(directory.relative_to(OUT)/f"state_{r['step']}.csv.gz").as_posix()
        r['step']=str(int(r['step'])+offset)
    history=prefix+current
    if history:
        baseline=history[0]
        for r in history:
            r['affine_change_from_initial']=str(float(r['affine_shear'])-float(baseline['affine_shear']))
            r['registry_change_from_initial']=str(float(r['registry_shear'])-float(baseline['registry_shear']))
    return history,root_initial


def bulk_spectrum_tail():
    """Independent 3x3 Fourier assembly at the PERFECT reference only.

    This is not a nonlinear-core Hessian substitute. The vector real-space
    Hessian product is checked independently before the tail comparison.
    """
    rows=NonlinearScrewRows(build_surface(tolerance=2e-11)[0],tolerance=2e-12)
    def symbols(cell,theta_j,theta_l):
        kernel=cell.kernel.evaluate(cell.reference,order=2)
        phase=theta_j[:,None]*cell.j+theta_l[:,None]*cell.l
        effective=kernel['hessian'][:,0]+2*rows.embedding.first_derivative(rows.rho_bulk)*kernel['hessian'][:,1]
        H=2*np.einsum('qr,rij->qij',1-np.cos(phase),effective)
        exp=np.exp(1j*phase)
        factor=(exp-1)[...,None]-(exp.conj()-1)[...,None]*cell.kernel.parity
        L=np.einsum('qrc,rci->qci',factor,kernel['gradient'])
        weights=np.r_[0.,rows.embedding.second_derivative(rows.rho_bulk),2*rows.angular_weights]
        H+=np.einsum('qci,c,qcj->qij',L.conj(),weights,L).real
        Lg=np.sum(kernel['gradient'][...,0]*cell.reference[:,2,None]*(1+cell.kernel.parity),axis=0)
        affine=(np.sum(effective[:,0,0]*cell.reference[:,2]**2)+np.sum(weights*Lg**2))/rows.h**2
        return H,affine
    small=VectorPeriodicCell(rows,(5,6),ring=4);j,l=np.indices(small.shape)
    tj,tl=2*np.pi/5,4*np.pi/6;polar=np.array([1.,.7,-.3])
    direction=np.cos(tj*j+tl*l)[...,None]*polar
    reference=symbols(small,np.array([tj]),np.array([tl]))[0][0]
    actual=small.evaluate(np.zeros_like(direction),direction=direction)['hessian_vector']
    error=float(np.max(abs(actual-np.cos(tj*j+tl*l)[...,None]*(reference@polar))))
    if error>2e-10:
        raise ArithmeticError('independent vector bulk Fourier Hessian check failed')
    records=[];matrices=[]
    for ring in (4,6,8,10,12,16,20,28):
        tick=time.perf_counter();cell=VectorPeriodicCell(rows,(24,24),ring=ring)
        qj,ql=np.indices(cell.shape);keep=(qj!=0)|(ql!=0)
        H,affine=symbols(cell,2*np.pi*qj[keep]/24,2*np.pi*ql[keep]/24)
        eigen,vectors=np.linalg.eigh(H);index=np.unravel_index(np.argmin(eigen),eigen.shape)
        minimum=float(min(eigen[index],affine));polarization=vectors[index[0],:,index[1]]
        matrices.append(H)
        records.append(dict(size=24,ring=ring,minimum_nongauge_curvature=minimum,
            affine_curvature=affine,independent_hessian_product_error=error,
            minimizing_theta_j=float(2*np.pi*qj[keep][index[0]]/24),
            minimizing_theta_l=float(2*np.pi*ql[keep][index[0]]/24),
            polarization_x=float(polarization[0]),polarization_y=float(polarization[1]),polarization_z=float(polarization[2]),
            elapsed_seconds=time.perf_counter()-tick,reference_is_perfect_bulk_only=True))
        print(f'perfect vector spectrum ring {ring}: lambda_min={minimum:.10g}',flush=True)
    for i,r in enumerate(records):
        r['difference_vs_ring28']=abs(r['minimum_nongauge_curvature']-records[-1]['minimum_nongauge_curvature'])
        r['maximum_full_symbol_error_vs_ring28']=float(np.max(np.linalg.norm(matrices[i]-matrices[-1],ord=2,axis=(-2,-1))))
        r['maximum_full_symbol_change_previous']=None if not i else float(np.max(np.linalg.norm(matrices[i]-matrices[i-1],ord=2,axis=(-2,-1))))
    write_csv(OUT/'perfect_vector_spectrum_tail.csv',records)


def report():
    surface,units,metadata=build_surface(tolerance=2e-11);rows=NonlinearScrewRows(surface,tolerance=2e-12)
    comparison=[];histories=[];states={};stability=[];layer_records=[]
    for directory in sorted(OUT.glob('*_n*_r*')):
        if not directory.is_dir():
            continue
        summary_path=directory/'summary.json'
        if not summary_path.exists() or not json.loads(summary_path.read_text())['completed']:
            print(f'INCOMPLETE case retained, not reported as converged: {directory.name}',flush=True)
            continue
        meta=json.loads((directory/'metadata.json').read_text())
        if meta['parameter_sha256']!=metadata['parameter_sha256']:
            raise ValueError('refusing to combine different material candidates')
        history,start=source_history(directory)
        first=history[0];last=history[-1]
        latest=max(directory.glob('state_*.csv.gz'),key=lambda p:int(p.name.split('_')[1].split('.')[0]))
        state,gamma=read_state(latest)
        states[directory.name]=(state,gamma)
        cell=VectorPeriodicCell(rows,state.shape[:2],ring=meta['transverse_ring'])
        for r in cell.layer_registry(state,gamma):
            layer_records.append(dict(case=directory.name,**r))
        phase=np.exp(2j*np.pi*state[...,0]/rows.b)
        phase*=np.exp(-1j*np.angle(phase.mean()))
        maximum_x_phase=float(np.max(abs(np.angle(phase)))*rows.b/(2*np.pi))
        comparison.append(dict(case=directory.name,size=meta['size_j'],ring=meta['transverse_ring'],
            initial_energy_ev=start['energy_ev'],zero_stress_relaxed_energy_ev=float(first['energy_ev_per_line_repeat']),
            initial_transverse_force=start['max_transverse_force'],
            zero_stress_maximum_force=float(first['maximum_force_residual']),
            initial_x_winding_cores=start['positive_cores']+start['negative_cores'],
            relaxed_x_winding_cores=int(first['positive_cores'])+int(first['negative_cores']),
            final_maximum_x_phase_over_L0=maximum_x_phase,
            final_maximum_transverse_displacement=float(np.max(np.linalg.norm(state[...,1:],axis=-1))),
            initial_affine_shear=start['gamma'],zero_stress_relaxed_shear=float(first['affine_shear']),
            unloaded_affine_change=float(last['affine_change_from_initial']),
            unloaded_x_projected_registry_change=float(last['registry_change_from_initial']),
            maximum_force_all_steps=max(float(r['maximum_force_residual']) for r in history),
            total_wall_seconds=json.loads(summary_path.read_text())['elapsed_seconds'],
            actual_initial_core_separation_nm=None if start['core_separation_over_L0'] is None else
                start['core_separation_over_L0']*units.length_scale_m*1e9))
        if float(first['applied_shear_mpa'])==0. and float(last['applied_shear_mpa'])==0.:
            first_state,first_gamma=read_state(OUT/first['state_source'])
            delta_phase=np.exp(2j*np.pi*(state[...,0]-first_state[...,0])/rows.b)
            delta_phase*=np.exp(-1j*np.angle(delta_phase.mean()))
            comparison[-1].update(
                unloaded_x_phase_change_over_L0=float(np.max(abs(np.angle(delta_phase)))*rows.b/(2*np.pi)),
                unloaded_transverse_change_over_L0=float(np.max(abs(state[...,1:]-first_state[...,1:]))))
        else:
            comparison[-1].update(unloaded_x_phase_change_over_L0=None,unloaded_transverse_change_over_L0=None)
        for r in history:
            histories.append(dict(case=directory.name,**r))
        for path in directory.glob('stability_*.json'):
            stability.append(dict(case=directory.name,step=int(path.stem.split('_')[-1]),
                                  **json.loads(path.read_text())))
    if not comparison:
        raise ValueError('no completed cases to report')
    write_csv(OUT/'domain_and_unload_summary.csv',comparison)
    write_csv(OUT/'all_static_stress_states.csv',[
        {('x_projected_'+key if 'registry_shear' in key or key=='registry_change_from_initial' else key):value
         for key,value in r.items()} for r in histories])
    write_csv(OUT/'vector_stability.csv',stability)
    write_csv(OUT/'final_layer_registry.csv',layer_records)
    save_json(OUT/'scientific_status.json',dict(
        parameter_sha256=metadata['parameter_sha256'],no_refit=True,
        completed_cases=len(comparison),actual_static_states=len(histories),
        maximum_force_residual=max(r['maximum_force_all_steps'] for r in comparison),
        all_initial_x_windings_removed=all(r['relaxed_x_winding_cores']==0 for r in comparison if r['initial_x_winding_cores']),
        x_winding_does_not_classify_full_vector_topology=True,
        scalar_x_registry_partition_does_not_define_vector_plastic_strain=True,
        maximum_unloaded_x_projected_registry_change=max(abs(r['unloaded_x_projected_registry_change']) for r in comparison),
        maximum_unloaded_actual_x_phase_change_over_L0=max(r['unloaded_x_phase_change_over_L0'] or 0. for r in comparison),
        maximum_unloaded_actual_transverse_change_over_L0=max(r['unloaded_transverse_change_over_L0'] or 0. for r in comparison),
        maximum_final_transverse_displacement=max(r['final_maximum_transverse_displacement'] for r in comparison),
        geometry='all 3 local displacements, infinite straight rows, fixed transverse periodic cell shape',
        result_is_not_a_physical_time_trajectory=True,
        supplied_close_opposite_pairs_are_not_Al_defect_population=True,
        residual_plastic_deformation_certified=False,experimental_yield_determined=False,
        material_calibration_accepted=False,finite_loop_or_source_validated=False,
        physical_mobilities_available=False,physical_seconds_available=False,physical_hz_available=False,
        production_pde_registered=False,ui_gate_passed=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for record in comparison:
        directory=OUT/record['case'];path=directory/'optimizer_progress.csv'
        if path.exists():
            with path.open() as stream:
                progress=[r for r in csv.DictReader(stream) if int(r['step'])==0]
            if progress:
                axes[0].semilogy([int(r['iteration']) for r in progress],
                    [max(abs(float(r['energy_ev'])),np.finfo(float).tiny) for r in progress],label=record['case'])
        hist=[r for r in histories if r['case']==record['case']]
        if len(hist)>1:
            axes[1].plot([float(r['applied_shear_mpa']) for r in hist],
                         [float(r['affine_change_from_initial']) for r in hist],'o-',
                         label=f"{record['size']}x{record['size']}, R={record['ring']}")
    # Actual nonzero-transverse checkpoint, not a fabricated core cartoon.
    available=sorted(OUT.glob('dipole_n32_r8/iteration_0_*.csv.gz'),key=lambda p:int(p.name.split('_')[2].split('.')[0]))
    if available:
        u,g=read_state(available[0]);j,l=np.indices(u.shape[:2]);yy=rows.d*(j+l/3);zz=rows.h*l
        p=axes[2].scatter(yy,zz,c=u[...,0],s=6,cmap='coolwarm')
        axes[2].quiver(yy,zz,u[...,1],u[...,2],scale_units='xy',scale=1.,width=.002)
        fig.colorbar(p,ax=axes[2],label='u_x / L0')
        axes[2].set_title('Actual vector checkpoint\nunamplified transverse arrows',fontsize=10)
    axes[0].set(xlabel='Optimizer iteration (NOT time)',ylabel='|excess energy| [eV / row repeat]',title='Zero-stress relaxation')
    axes[1].set(xlabel='Imposed resolved shear [MPa]',ylabel='Affine shear change from relaxed state',title='Static loading / unloading')
    axes[2].set(xlabel='y / L0',ylabel='z / L0',aspect='equal')
    for ax in axes[:2]:
        ax.legend(fontsize=6);ax.grid(alpha=.25)
    fig.suptitle('Unchanged LJ/Bessel candidate; static research, not physical-Hz fatigue',fontsize=11)
    fig.tight_layout();fig.savefig(OUT/'vector_core_audit.png',dpi=160);plt.close(fig)
    print('completed static states reported; physical/material gates remain closed',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['validation','tail','bulk-tail','report'])
    parser.add_argument('--state',default='dipole_n32_r8/iteration_0_20.csv.gz')
    args=parser.parse_args()
    if args.command=='validation':
        validation()
    elif args.command=='tail':
        tail_study(args.state)
    elif args.command=='bulk-tail':
        bulk_spectrum_tail()
    else:
        report()
