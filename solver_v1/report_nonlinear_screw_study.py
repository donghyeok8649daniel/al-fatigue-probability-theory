"""Re-evaluate frozen-coordinate forces and report the actual static study.

This is NOT the relaxation runner: it reads its saved atomic states, checks
their unchanged parameter binding, and independently computes diagnostic
forces/refinement. A plot is not a new time-dependent simulation.
"""
from __future__ import annotations

import csv
import gzip
import json
import time

import numpy as np

from .nonlinear_fcc_screw import NonlinearScrewRows
from .nonlinear_screw_transverse_audit import omitted_transverse_forces, transverse_row_coefficients
from .run_low_stress_cyclic_diagnostic import build_surface, save_json, write_csv
from .run_nonlinear_screw_study import EV_J, OUT


def read_states():
    with gzip.open(OUT/'static_states.csv.gz','rt',encoding='utf8') as stream:
        records=list(csv.DictReader(stream))
    output={}
    for r in records:
        key=(int(r['size']),r['case'],int(r['step']))
        if key not in output:
            output[key]=(np.zeros((key[0],key[0])),float(r['gamma']))
        output[key][0][int(r['j']),int(r['l'])]=float(r['u_over_L0'])
    return output


def independent_transverse_energy_check(surface):
    """Direct real-index nonlinear energy re-evaluation after perturbing y,z.

    This is an independent finite-difference validation, not the canonical
    vector energy implementation. A small deliberately finite validation ring
    is identical on both sides; infinite omission is audited separately.
    """
    rows=NonlinearScrewRows(surface,fixed_ring=2)
    cell=rows.periodic_cell((4,5));j,l=np.indices(cell.shape)
    u=.077*np.cos(.7*j+.9*l);gamma=.013
    v=np.stack([np.sin(.43*j+.61*l),np.cos(.71*j-.3*l)],axis=-1)
    shifts=np.array([np.roll(u,(-int(jr),-int(lr)),axis=(0,1))-u+gamma*rows.h*lr
                     for jr,lr in zip(rows.j,rows.l)])
    dv=np.array([np.roll(v,(-int(jr),-int(lr)),axis=(0,1))-v for jr,lr in zip(rows.j,rows.l)])
    jr=np.repeat(rows.j,cell.size);lr=np.repeat(rows.l,cell.size)
    y0=rows.d*(jr+lr/3);z0=rows.h*lr
    def energy(step):
        values=np.zeros((cell.size,17))
        for mode in range(len(rows.g)+1):
            coef,_=transverse_row_coefficients(rows,jr,lr,mode,
                y=y0+step*dv[...,0].ravel(),z=z0+step*dv[...,1].ravel())
            base,_=transverse_row_coefficients(rows,rows.j,rows.l,mode)
            phase=np.exp(2j*np.pi*mode/rows.b*shifts.reshape(len(rows.j),cell.size))
            values+=np.real(coef.reshape(len(rows.j),cell.size,17)*phase[...,None]-base[:,None,:]).sum(axis=0)
        return (.5*values[:,0].sum()+np.sum(rows.embedding.value(rows.rho_bulk+values[:,1])-
            rows.embedding.value(rows.rho_bulk))+np.sum(rows.angular_weights*values[:,2:]**2))
    analytic=float(np.sum(omitted_transverse_forces(cell,u,gamma,ring=2)['gradient']*v))
    results=[]
    for step in (4e-6,2e-6,1e-6):
        derivative=(energy(step)-energy(-step))/(2*step)
        results.append(dict(step=step,analytic_directional_gradient=analytic,
            finite_difference_directional_gradient=derivative,absolute_error=abs(derivative-analytic)))
    write_csv(OUT/'transverse_energy_derivative_validation.csv',results)
    if results[-1]['absolute_error']>1e-7:
        raise ArithmeticError('independent transverse energy derivative check failed')
    return results


def report():
    started=time.perf_counter()
    surface,units,metadata=build_surface(tolerance=2e-11)
    saved=json.loads((OUT/'metadata.json').read_text(encoding='utf8'))
    if metadata['parameter_sha256']!=saved['parameter_sha256']:
        raise ValueError('saved states and current material parameters differ')
    execution=json.loads((OUT/'execution_summary.json').read_text(encoding='utf8'))
    if not execution['all_force_balanced']:
        raise ValueError('the relaxation study has not completed with force balance')
    rows=NonlinearScrewRows(surface,tolerance=2e-12)
    independent=independent_transverse_energy_check(surface)
    states=read_states();summary=[];audits=[];reference_forces={}
    sizes=execution['domain_sizes']
    with open(OUT/'stress_history.csv',encoding='utf8') as stream:
        history=list(csv.DictReader(stream))
    for size in sizes:
        cell=rows.periodic_cell((size,size))
        initial=[r for r in history if int(r['size_j'])==size and r['case']=='dipole' and int(r['step'])==0][0]
        final=[r for r in history if int(r['size_j'])==size and r['case']=='dipole' and int(r['step'])==10][0]
        peak=[r for r in history if int(r['size_j'])==size and r['case']=='dipole' and int(r['step'])==4][0]
        summary.append(dict(size=size,transverse_period_y_nm=size*rows.d*units.length_scale_m*1e9,
            transverse_period_z_nm=size*rows.h*units.length_scale_m*1e9,
            actual_core_separation_nm=float(initial['core_separation_over_L0'])*units.length_scale_m*1e9,
            zero_stress_energy_ev=float(initial['energy_ev_per_line_repeat']),
            initial_registry_shear=float(initial['registry_shear']),
            initial_affine_shear=float(initial['affine_shear']),
            peak50_affine_change=float(peak['affine_shear_change_from_initial']),
            peak50_new_registry=float(peak['registry_shear_change_from_initial']),
            unloaded_new_registry=float(final['registry_shear_change_from_initial']),
            unloaded_affine_change=float(final['affine_shear_change_from_initial']),
            unloaded_energy_change_ev=float(final['energy_ev_per_line_repeat'])-float(initial['energy_ev_per_line_repeat']),
            same_discrete_core_locations=(initial['core_positions']==final['core_positions'])))
        for step in (0,4):
            field,gamma=states[(size,'dipole',step)]
            previous=None
            for ring in (6,10,14,20):
                tick=time.perf_counter()
                out=omitted_transverse_forces(cell,field,gamma,ring=ring)
                force=out['gradient']
                norm=float(np.max(np.linalg.norm(force,axis=-1)))
                audits.append(dict(size=size,step=step,ring=ring,
                    max_omitted_y_gradient_ev_per_L0=out['maximum_y_gradient'],
                    max_omitted_z_gradient_ev_per_L0=out['maximum_z_gradient'],
                    max_omitted_norm_ev_per_L0=norm,
                    max_omitted_force_nN=norm*EV_J/units.length_scale_m*1e9,
                    maximum_change_from_previous_ring=None if previous is None else float(np.max(abs(force-previous))),
                    net_force_residual=float(np.linalg.norm(out['net_gradient'])),
                    wall_seconds=time.perf_counter()-tick))
                print(f'frozen-coordinate audit size={size}, step={step}, ring={ring}: '
                      f'max transverse gradient={norm:.8g} eV/L0',flush=True)
                previous=force
            reference_forces[(size,step)]=force
            write_csv(OUT/'transverse_force_audit.csv',audits)
    write_csv(OUT/'domain_and_unload_summary.csv',summary)
    reference=summary[-1]
    for r in summary:
        r['energy_change_vs_largest_domain_ev']=r['zero_stress_energy_ev']-reference['zero_stress_energy_ev']
        r['peak_response_change_vs_largest_domain']=r['peak50_affine_change']-reference['peak50_affine_change']
    write_csv(OUT/'domain_and_unload_summary.csv',summary)
    changed=any(float(r['registry_shear_change_from_initial'])!=0. for r in history)
    same_final_topology=all(r['same_discrete_core_locations'] and r['unloaded_new_registry']==0. for r in summary)
    status=dict(static_scalar_force_balance=execution['all_force_balanced'],
        new_net_registry_transfer_observed=changed,
        final_registry_and_core_locations_recovered=same_final_topology,
        maximum_unloaded_affine_change=max(abs(r['unloaded_affine_change']) for r in summary),
        interval_mpa=[min(float(r['applied_shear_mpa']) for r in history),
                      max(float(r['applied_shear_mpa']) for r in history)],experimental_strength_determined=False,
        full_vector_core_validated=False,
        omitted_transverse_force_max_ev_per_L0=max(r['max_omitted_norm_ev_per_L0'] for r in audits if r['ring']==20),
        max_last_ring_force_change=max(r['maximum_change_from_previous_ring'] for r in audits if r['ring']==20),
        independent_transverse_energy_derivative_error=independent[-1]['absolute_error'],
        material_candidate_accepted=False,physical_seconds_hz_available=False,
        nonlinear_infinite_row_identity_verified=True,
        finite_nucleation_energy_available=False,production_ui_gate=False,
        conclusion='restricted scalar force balance only; inspect omitted transverse forces and refinements before a vector-core claim',
        wall_seconds=time.perf_counter()-started)
    save_json(OUT/'physical_and_numerical_status.json',status)
    make_plot(history,states,reference_forces,sizes[0],rows)
    print(json.dumps(status,indent=2),flush=True)


def make_plot(history,states,forces,size,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for case in ('perfect','dipole'):
        data=[r for r in history if int(r['size_j'])==size and r['case']==case]
        axes[0,0].plot([float(r['applied_shear_mpa']) for r in data],
                       [float(r['affine_shear_change_from_initial']) for r in data],'o-',label=case)
        axes[0,1].plot([int(r['step']) for r in data],
                       [float(r['registry_shear_change_from_initial']) for r in data],'o-',label=case)
    axes[0,0].set(xlabel='Applied shear stress [MPa]',ylabel='Affine shear change from initial [-]',
                  title='Athermal static load / unload; not frequency-dependent fatigue')
    axes[0,1].set(xlabel='Static load step (NOT time)',ylabel='New registry shear [-]',
                  title='Existing slip content subtracted; no resolved new motion')
    axes[0,0].legend();axes[0,1].legend()
    u,gamma=states[(size,'dipole',0)]
    im=axes[1,0].imshow(np.mod(u/rows.b,1.).T,origin='lower',cmap='twilight',vmin=0,vmax=1)
    axes[1,0].set(xlabel='Logical atomic row j',ylabel='Atomic layer l',title='Relaxed pre-existing screw pair: row phase u/b')
    fig.colorbar(im,ax=axes[1,0],label='Periodic row phase [-]')
    norm=np.linalg.norm(forces[(size,0)],axis=-1)
    im=axes[1,1].imshow(norm.T,origin='lower',cmap='magma')
    axes[1,1].set(xlabel='Logical atomic row j',ylabel='Atomic layer l',title='Forces in FROZEN transverse coordinates: not equilibrated')
    fig.colorbar(im,ax=axes[1,1],label='Omitted energy-gradient norm [eV/L0]')
    fig.suptitle('Same LJ/Bessel research candidate; scalar anti-plane only; NO material/kinetic calibration')
    fig.savefig(OUT/'static_defect_audit.png',dpi=150)
    plt.close(fig)


if __name__=='__main__':
    report()
