"""Reproduce static common-energy and coordinate diagnostics, no Si calibration.

Run: python -m results.silicon_wafer_feasibility.run_unified_coordinate_audit
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, root

from solver_v1.silicon_environment_research import (
    DiamondCell, RigidInterface, environment_jet, relaxed_hessian,
)
from .run_static_probe import RigidSlab, pair, source_parameters, PARAMETER_URL


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/silicon_unified_v1'
EV_A3_TO_GPA = 160.2176634


def write_csv(name, rows):
    with (OUT/name).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def difference(a, b):
    return {'energy_eV': abs(a.value-b.value),
            'gradient_max_abs': float(np.max(abs(a.gradient-b.gradient))),
            'hessian_max_abs': float(np.max(abs(a.hessian-b.hessian)))}


def bulk_audit(p, lattice):
    cell = DiamondCell(p, lattice)
    zero = cell.evaluate(np.zeros(9))
    effective, response = relaxed_hessian(zero, range(6))
    probes, finite = [], []
    for state in [np.zeros(9), np.array([.01,-.02,.015,.025,0.,.03,.01,-.02,.005]),
                  np.array([-.035,.025,0.,0.,-.04,.02,-.01,.006,.013])]:
        a = cell.evaluate(state)
        b = cell.evaluate(state, 'direct')
        probes.append({'state':state.tolist(), **difference(a,b)})
    for gamma in [4e-3, 2e-3, 1e-3, 5e-4]:
        state = np.zeros(9); state[5] = gamma
        def equation(u):
            q = state.copy(); q[6:] = u
            j = cell.evaluate(q, 'direct')
            return j.gradient[6:], j.hessian[6:,6:]
        solved = root(equation, response[:,5]*gamma, jac=True, options={'xtol':1e-10})
        state[6:] = solved.x
        energy = cell.evaluate(state, 'direct', derivatives=False)
        residual = float(np.max(abs(solved.fun)))
        assert residual < 1e-10
        finite.append({'gamma_xy':gamma, 'ux_A':solved.x[0], 'uy_A':solved.x[1], 'uz_A':solved.x[2],
            'force_residual_eV_A':residual, 'root_success_flag':bool(solved.success),
            'C44_from_energy_GPa':2*(energy-zero.value)/gamma**2/cell.volume*EV_A3_TO_GPA,
            'C44_error_GPa':abs(2*(energy-zero.value)/gamma**2-effective[5,5])/cell.volume*EV_A3_TO_GPA})
    write_csv('bulk_internal_relaxation.csv', finite)
    clamped = zero.hessian[:6,:6]/cell.volume*EV_A3_TO_GPA
    relaxed = effective/cell.volume*EV_A3_TO_GPA
    assert finite[-1]['C44_error_GPa'] < finite[0]['C44_error_GPa']/20
    bigger = DiamondCell(p, lattice, image_shell=3).evaluate(np.zeros(9))
    return {'primitive_cell_atoms':2, 'primitive_cell_volume_A3':cell.volume,
        'equilibrium_energy_eV_cell':zero.value, 'equilibrium_gradient_max_abs':float(np.max(abs(zero.gradient))),
        'clamped_C_GPa':clamped.tolist(), 'relaxed_C_GPa':relaxed.tolist(),
        'C44_reduction_fraction':float(1-relaxed[5,5]/clamped[5,5]),
        'internal_hessian_eV_A2':zero.hessian[6:,6:].tolist(),
        'internal_response_A_per_engineering_strain':response.tolist(),
        'moment_direct_probes':probes, 'image_shell_2_to_3_difference':difference(zero,bigger),
        'finite_strain_checks':finite, 'dynamical_elimination_certified':False}


def continuation(interface, opening, intervals):
    rows, previous = [], 0.
    for fraction in np.linspace(0., .3, intervals+1):
        q = np.array([opening, fraction*interface.period, previous])
        starting = interface.evaluate(q)
        if starting.hessian[2,2] <= 0:
            return rows, 'nonpositive transverse predictor curvature'
        def equation(y):
            state = q.copy(); state[2] = float(y[0])
            j = interface.evaluate(state)
            return j.gradient[2:], j.hessian[2:,2:]
        if abs(starting.gradient[2]) < 1e-12:
            y = previous
        else:
            # A Newton predictor avoids MINPACK's tiny initial trust region when
            # the previous root is a nonzero roundoff-sized symmetry coordinate.
            predictor = previous-starting.gradient[2]/starting.hessian[2,2]
            solved = root(equation, [predictor], jac=True, options={'xtol':1e-10})
            y = float(solved.x[0])
        q[2] = y
        relaxed = interface.evaluate(q)
        if abs(relaxed.gradient[2]) > 1e-10 or relaxed.hessian[2,2] <= 0:
            return rows, 'unresolved transverse stationary minimum'
        if abs(y-previous) > .15*interface.period:
            return rows, 'continuation step too large to identify the same branch'
        previous = y
        clamped_q = q.copy(); clamped_q[2] = 0.
        clamped = interface.evaluate(clamped_q)
        effective, tangent = relaxed_hessian(relaxed, [0,1])
        rows.append({'opening_A':opening, 'slip_fraction':fraction,
            'u1_A':q[1], 'u2_relaxed_A':y,
            'clamped_energy_eV_cell':clamped.value, 'relaxed_energy_eV_cell':relaxed.value,
            'lowering_eV_cell':clamped.value-relaxed.value,
            'omitted_gradient_eV_A_cell':clamped.gradient[2],
            'relaxed_gradient_residual_eV_A_cell':relaxed.gradient[2],
            'relaxed_H22_eV_A2_cell':relaxed.hessian[2,2],
            'retained_H11_clamped_at_relaxed_point':relaxed.hessian[1,1],
            'retained_H11_Schur_eV_A2_cell':effective[1,1],
            'du2_du1':tangent[0,1]})
    return rows, 'completed stationary branch; static only'


def interface_audit(p, lattice, kind):
    slab = RigidSlab(p, lattice, kind, repeats=1, periods=3, image_shell=3)
    interface = RigidInterface(slab, image_shell=5)
    samples, finite = [], []
    for q in [np.zeros(3), np.array([.3,.2*slab.period,0.]),
              np.array([.45,.18*slab.period,.08*slab.period])]:
        j = interface.evaluate(q)
        direct = interface.evaluate(q, 'direct')
        def independent(state):
            positions = slab.positions.copy()
            positions[slab.upper] += state@interface.frame
            return (slab.total_energy(positions)-slab.reference)/slab.cells
        samples.append({'state_A':q.tolist(), 'energy_eV_cell':j.value,
            'gradient_eV_A_cell':j.gradient.tolist(), 'hessian_eV_A2_cell':j.hessian.tolist(),
            'direct_difference':difference(j,direct), 'legacy_numpy_energy_difference_eV_cell':abs(j.value-independent(q))})
        if np.linalg.norm(q) == 0:
            continue
        for h in [1e-3,5e-4,2.5e-4]:
            basis = np.eye(3)*h
            gradient = np.array([(independent(q+d)-independent(q-d))/(2*h) for d in basis])
            hessian = np.empty((3,3))
            for a in range(3):
                for b in range(3):
                    if a == b:
                        hessian[a,a] = (independent(q+basis[a])-2*independent(q)+independent(q-basis[a]))/h**2
                    else:
                        hessian[a,b] = sum(sa*sb*independent(q+sa*basis[a]+sb*basis[b])
                            for sa in [-1,1] for sb in [-1,1])/(4*h**2)
            finite.append({'cut':kind, 'opening_A':q[0], 'u1_A':q[1], 'u2_A':q[2], 'step_A':h,
                'gradient_error_eV_A_cell':float(np.max(abs(gradient-j.gradient))),
                'hessian_error_eV_A2_cell':float(np.max(abs(hessian-j.hessian)))})
    geometry = []
    for repeat, periods, image_shell in [(1,4,5),(2,3,5),(2,4,6)]:
        other = RigidInterface(RigidSlab(p,lattice,kind,repeats=repeat,periods=periods,image_shell=3),image_shell=image_shell)
        q = np.array([.3,.2*slab.period,.07*slab.period])
        geometry.append({'repeat':repeat, 'periods':periods, 'image_shell':image_shell,
                         **difference(interface.evaluate(q),other.evaluate(q))})
    profiles, statuses, refinements = [], [], []
    for opening in [0.,.3,.75]:
        fine, status = continuation(interface, opening, 60)
        coarse, coarse_status = continuation(interface, opening, 30)
        statuses.append({'opening_A':opening, 'fine':status, 'coarse':coarse_status,
                         'fine_points':len(fine), 'coarse_points':len(coarse)})
        if len(fine) == 61 and len(coarse) == 31:
            error = max(abs(a['relaxed_energy_eV_cell']-b['relaxed_energy_eV_cell'])
                        for a,b in zip(fine[::2],coarse))
            y_error = max(abs(a['u2_relaxed_A']-b['u2_relaxed_A']) for a,b in zip(fine[::2],coarse))
            refinements.append({'opening_A':opening,'energy_error_eV_cell':error,'u2_error_A':y_error})
            assert error < 1e-9 and y_error < 1e-8
        profiles.extend(fine)
    write_csv(f'{kind}_transverse_profiles.csv', profiles)
    write_csv(f'{kind}_derivative_refinement.csv', finite)
    return {'cut':kind, 'area_A2':slab.area, 'period_A':slab.period,
        'state_samples':samples, 'geometry_convergence':geometry,
        'finite_difference_refinement':finite, 'branch_status':statuses,
        'continuation_refinement':refinements, 'actual_crack_barrier_certified':False}, profiles


def environment_audit(p, bond):
    tetra = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1.]])*bond/np.sqrt(3)
    rows = []
    for scale in [.88,1.,1.13]:
        f = np.array([[scale,.11,0.],[0.,1.,-.04],[.02,0.,1.]])
        x = tetra@f.T
        a,b = (environment_jet(x,p,m) for m in ['direct','moments'])
        rows.append({'scale':scale, **difference(a,b),
                     'torque_max_abs_eV':float(np.max(abs(np.cross(x,-a.gradient.reshape(-1,3)).sum(axis=0))))})
    return rows


def plot_result(summary, profiles):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    c0 = summary['bulk']['clamped_C_GPa'][5][5]
    c1 = summary['bulk']['relaxed_C_GPa'][5][5]
    axes[0].bar(['Fixed basis','Relaxed basis'],[c0,c1],color=['#9ca3af','#2563eb'])
    axes[0].set(ylabel='C44 (GPa)',title='Same SW energy, different constraints')
    for i,c in enumerate([c0,c1]): axes[0].text(i,c+2,f'{c:.2f}',ha='center')
    rows = [r for r in profiles['glide'] if r['opening_A'] == .3]
    x = np.array([r['slip_fraction'] for r in rows])
    axes[1].plot(x,[r['clamped_energy_eV_cell'] for r in rows],label='u2 fixed at zero',color='#9ca3af')
    axes[1].plot(x,[r['relaxed_energy_eV_cell'] for r in rows],label='u2 stationary branch',color='#2563eb')
    axes[1].set(xlabel='u1 / registry period',ylabel='Energy (eV / interface cell)',title='Si(111) glide, opening = 0.3 A')
    axes[1].legend(fontsize=8)
    axes[2].plot(x,[r['omitted_gradient_eV_A_cell'] for r in rows],color='#dc2626')
    axes[2].axhline(0,color='black',linewidth=.5)
    axes[2].set(xlabel='u1 / registry period',ylabel='dW / du2 (eV / A / cell)',title='Force omitted by the two-coordinate slice')
    for axis in axes: axis.grid(axis='y',alpha=.2)
    fig.suptitle('Original SW, pure Si, 0 K static audit — not a fatigue or lifetime calibration',fontsize=11)
    fig.savefig(OUT/'coordinate_audit.png',dpi=180)
    fig.savefig(OUT/'coordinate_audit.svg')
    plt.close(fig)


def main():
    start = time.perf_counter(); OUT.mkdir(exist_ok=True)
    p, sha = source_parameters()
    bond = brentq(lambda x:float(pair(x,p,derivative=True)),2.2,2.5,xtol=1e-14)
    lattice = 4*bond/np.sqrt(3)
    summary = {'scope':'original SW static representation and coordinate audit; not Si calibration',
        'source_url':PARAMETER_URL,'source_sha256':sha,'parameters':p,'temperature_K':0,
        'lattice_A':lattice,'bond_A':bond,'physical_mobility':None,'physical_time_scale_seconds':None,
        'production_modified':False,'material_calibration_accepted':False,
        'source_files_sha256':{str(path.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [Path(__file__), ROOT/'solver_v1/silicon_environment_research.py']}}
    summary['environment_checks'] = environment_audit(p,bond)
    summary['bulk'] = bulk_audit(p,lattice)
    print('Bulk audit complete',flush=True)
    profiles = {}
    for kind in ['shuffle','glide']:
        summary[kind], profiles[kind] = interface_audit(p,lattice,kind)
        print(kind, summary[kind]['branch_status'],flush=True)
    summary['elapsed_seconds'] = time.perf_counter()-start
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    plot_result(summary,profiles)
    print(json.dumps({'output':'results/silicon_unified_v1','elapsed_seconds':summary['elapsed_seconds']}),flush=True)


if __name__ == '__main__':
    main()
