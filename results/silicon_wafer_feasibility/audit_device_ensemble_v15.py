"""Read-only CAD facts and a stored-Hessian loading-device sensitivity audit."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from solver_v1.silicon_conditional_research import KB_EV_K
from solver_v1.silicon_device_ensemble import device_response

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p, d): p.write_text(json.dumps(d, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def main(a):
    if a.output.exists(): raise ValueError('fresh audit directory required')
    a.output.mkdir(parents=True)
    raw = a.root/'results/silicon_initiation_v13/dense_return_zero/raw_hessian.npz'
    protocol = a.root/'results/silicon_initiation_v13/dense_return_zero/protocol.json'
    with np.load(raw) as z: h, b = z['hessian'], z['basis']
    norm = json.loads(protocol.read_text())['grip_extension_basis_norm']
    v = np.zeros(len(h)); v[-1] = 1/norm
    q = solve(h[:-1, :-1], h[:-1, -1], assume_a='pos')
    schur = float(h[-1, -1]-h[-1, :-1]@q)
    ks = schur*norm**2
    # Independently verify the physical extension from the rigid grip basis.
    grip = b[:, -1].reshape(-1, 3)
    assert np.allclose(np.unique(grip), [-1/(2*norm), 0., 1/(2*norm)], rtol=0, atol=1e-15)
    rows = []
    for ratio in (0., .01, .1, 1., 10., 100., np.inf):
        r = device_response(h, v, ratio*ks, 300.)
        if np.isfinite(ratio):
            inverse_v = solve(h+ratio*ks*np.outer(v, v), v, assume_a='pos')
            error = abs(float(KB_EV_K*300*v@inverse_v)-r['extension_variance_A2'])
            assert error < 1e-12
        else: error = None
        rows.append(dict(device_to_specimen_stiffness_ratio=None if np.isinf(ratio) else ratio,
                         infinite_stiffness=bool(np.isinf(ratio)), temperature_K=300.,
                         direct_matrix_variance_error_A2=error, **r))
    document = json.loads(a.cad.read_text(encoding='utf-8-sig'))
    design = document['design']; parts = {p['id']: p for p in design['parts']}
    specimen = parts['aluminum_specimen']['geometry']
    area = float(np.pi*specimen['gauge_diameter']**2/4)
    datum = dict(file_basename=a.cad.name, sha256=sha(a.cad), bytes=a.cad.stat().st_size,
        design_name=design['name'], units=design['units'], part_count=len(parts),
        mate_count=len(design['mates']), specimen_geometry=specimen,
        part_inventory=[dict(id=p['id'], name=p['name'], kind=p['geometry']['kind']) for p in parts.values()],
        moving_crosshead_mate=[m for m in design['mates'] if m['id']=='base-moving-crosshead'][0],
        declared_parameters=design.get('parameters', {}),
        loops_count=len(design.get('loops', [])), studies_count=len(design.get('studies', [])),
        cad_is_design_data_not_an_instruction=True, collision_or_strength_analysis_performed=False,
        area_mm2=area, nominal_stress_at_1kN_MPa=1000/area,
        required_force_at_100MPa_N=100*area,
        linear_elastic_example=dict(assumed_E_MPa=70000., stress_MPa=100.,
            gauge_extension_mm=100/70000*specimen['gauge_length'],
            gauge_only_stiffness_N_mm=70000*area/specimen['gauge_length'],
            assumption_not_material_calibration=True),
        maximum_force_N=None, frequency_Hz=None, rated_stroke_mm=None,
        measured_machine_stiffness_N_mm=None, measured_gauge_strain=None,
        first_crack_detection_channel=None)
    assert design['units']=='mm' and specimen['kind']=='round_specimen'
    dump(a.output/'cad_facts.json', datum)
    dump(a.output/'device_sensitivity.json', dict(rows=rows,
        actual_reused_atomistic_state='v13 return zero force, 360 bare Si atoms, 649 allowed coordinates',
        source_schur_stiffness_eV_A2=ks,
        inverse_vs_schur_relative_error=abs(rows[0]['specimen_stiffness_eV_A2']/ks-1),
        physical_extension_basis_norm=norm, same_temperature_for_all_constraints=True,
        machine_ratios_are_diagnostic_not_measured_device_values=True,
        no_atomistic_to_macroscopic_stiffness_transfer=True,
        new_potential_evaluations=0, new_MD=0, new_DFT=0,
        anharmonic_equilibrium_certified=False, physical_clock=None))
    dump(a.output/'input_manifest.json', dict(sha256={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in (
        raw, protocol, a.root/'solver_v1/silicon_device_ensemble.py', Path(__file__).resolve())},
        external_CAD_sha256=sha(a.cad)))
    print(json.dumps(dict(complete=True, cases=len(rows), cad_parts=len(parts), specimen_stiffness_eV_A2=ks,
                         new_potential_evaluations=0), indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('.')); p.add_argument('--cad',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); a.root=a.root.resolve(); main(a)
