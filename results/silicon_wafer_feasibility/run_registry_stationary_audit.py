"""Stationary states on a specified rigid glide registry path, not crack kinetics.

Run: python -m results.silicon_wafer_feasibility.run_registry_stationary_audit
"""
from pathlib import Path
import csv
import hashlib
import json
import time

import numpy as np
from scipy.optimize import brentq, minimize

from solver_v1.silicon_environment_research import RigidInterface
from .run_static_probe import RigidSlab, pair, source_parameters


def main():
    start = time.perf_counter()
    out = Path('results/silicon_unified_v1')
    out.mkdir(exist_ok=True)
    p, checksum = source_parameters()
    bond = brentq(lambda x:float(pair(x,p,derivative=True)),2.2,2.5,xtol=1e-14)
    lattice = 4*bond/np.sqrt(3)
    interface = RigidInterface(RigidSlab(p,lattice,'glide',repeats=1,periods=3,image_shell=3),image_shell=6)
    along, transverse = np.array([np.sqrt(3)/2,-.5]), np.array([.5,np.sqrt(3)/2])
    length = interface.period/np.sqrt(3)
    records, paths = [], []
    for opening in [0.,.3,.75]:
        def evaluate(s):
            return interface.evaluate([opening,*(s*along)])
        def derivative(s):
            return float(evaluate(s).gradient[1:]@along)
        roots = []
        for count in [80,160]:
            points = np.linspace(0.,length,count+1)
            gradients = np.array([derivative(s) for s in points])
            found = [0.,length]
            for i in range(1,count):
                if abs(gradients[i]) < 1e-11:
                    found.append(points[i])
                elif gradients[i]*gradients[i+1] < 0:
                    found.append(brentq(derivative,points[i],points[i+1],xtol=1e-13))
            unique = []
            for x in sorted(found):
                if not unique or abs(x-unique[-1]) > 1e-8:
                    unique.append(float(x))
            roots.append(unique)
        assert len(roots[0]) == len(roots[1])
        assert np.max(abs(np.array(roots[0])-roots[1])) < 1e-10
        base = evaluate(0.).value
        stationary = []
        for s in roots[1]:
            jet = evaluate(s)
            eigenvalues, vectors = np.linalg.eigh(jet.hessian[1:,1:])
            residual = float(np.max(abs(jet.gradient[1:])))
            assert residual < 1e-9
            index = int(np.sum(eigenvalues < 0))
            detail = {'s_A':s,'s_over_partial':s/length,'registry_u_A':(s*along).tolist(),
                'energy_eV_cell':jet.value,'energy_above_start_eV_cell':jet.value-base,
                'registry_gradient_max_eV_A_cell':residual,
                'registry_hessian_eigenvalues_eV_A2_cell':eigenvalues.tolist(),
                'registry_morse_index':index,'normal_derivative_eV_A_cell':float(jet.gradient[0])}
            if index == 1:
                destinations = []
                for sign in [-1,1]:
                    def energy(z):
                        j=interface.evaluate([opening,*z])
                        return j.value,j.gradient[1:]
                    initial = s*along+sign*.03*interface.period*vectors[:,0]
                    minimum = minimize(energy,initial,jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':200})
                    final = interface.evaluate([opening,*minimum.x])
                    eig = np.linalg.eigvalsh(final.hessian[1:,1:])
                    assert np.max(abs(final.gradient[1:])) < 2e-8 and eig[0] > 0
                    assert final.value < jet.value
                    destinations.append({'registry_u_A':minimum.x.tolist(),'energy_eV_cell':final.value,
                        'gradient_max_eV_A_cell':float(np.max(abs(final.gradient[1:]))),
                        'hessian_eigenvalues_eV_A2_cell':eig.tolist(),'optimizer_success_flag':bool(minimum.success)})
                assert np.linalg.norm(np.array(destinations[0]['registry_u_A'])-destinations[1]['registry_u_A']) > .05
                detail['downhill_stationary_destinations'] = destinations
            stationary.append(detail)
        path = []
        for s in np.linspace(0.,length,161):
            j=evaluate(s)
            path.append({'opening_A':opening,'s_A':s,'s_over_partial':s/length,
                'energy_eV_cell':j.value,'transverse_gradient_eV_A_cell':float(j.gradient[1:]@transverse),
                'transverse_curvature_eV_A2_cell':float(transverse@j.hessian[1:,1:]@transverse)})
        max_gradient=max(abs(r['transverse_gradient_eV_A_cell']) for r in path)
        min_curvature=min(r['transverse_curvature_eV_A2_cell'] for r in path)
        assert max_gradient < 1e-10 and min_curvature > 0
        periodicity = []
        for fraction in [.0,.15,.55]:
            base_jet, partial_jet, full_jet = (evaluate((fraction+shift)*length) for shift in [0.,1.,3.])
            delta = {'fraction':fraction,
                'partial_shift_energy_difference_eV_cell':partial_jet.value-base_jet.value,
                'full_shift_energy_error_eV_cell':abs(full_jet.value-base_jet.value),
                'full_shift_gradient_error_eV_A_cell':float(np.max(abs(full_jet.gradient-base_jet.gradient))),
                'full_shift_hessian_error_eV_A2_cell':float(np.max(abs(full_jet.hessian-base_jet.hessian)))}
            assert delta['full_shift_energy_error_eV_cell'] < 1e-10
            assert delta['full_shift_gradient_error_eV_A_cell'] < 1e-10
            assert delta['full_shift_hessian_error_eV_A2_cell'] < 1e-9
            periodicity.append(delta)
        records.append({'opening_A':opening,'stationary':stationary,
            'max_transverse_gradient_eV_A_cell':max_gradient,'min_transverse_curvature_eV_A2_cell':min_curvature,
            'stationary_scan_intervals':[80,160], 'endpoint_energy_difference_eV_cell':evaluate(length).value-base,
            'periodicity_checks':periodicity,
            'energy_at_0_L_2L_3L_eV_cell':[evaluate(k*length).value for k in range(4)]})
        paths.extend(path)
        print('opening',opening,'stationary fractions',[round(x/length,9) for x in roots[1]],flush=True)
    result={'scope':'fixed-opening rigid glide registry states; not localized crack activation',
        'source_sha256':checksum,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'module_sha256':hashlib.sha256(Path('solver_v1/silicon_environment_research.py').read_bytes()).hexdigest(),
        'temperature_K':0,'lattice_A':lattice,'partial_translation_A':length,
        'full_registry_period_A':3*length,
        'direction_cubic':(along@interface.frame[1:]).tolist(),
        'endpoint_lattice_coordinates':[2/3,-1/3],
        'kinetic_closure_certified':False,'physical_activation_barrier_eV':None,
        'records':records,'elapsed_seconds':time.perf_counter()-start}
    (out/'registry_stationary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (out/'registry_path.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(paths[0]));writer.writeheader();writer.writerows(paths)
    print('completed',result['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    main()
