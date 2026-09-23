"""Reproduce paper arithmetic and expose alternative geometry conventions."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_oxide_research import (
    axial_composite, rectangular_shell, inverse_equivalent_flaw_nm, precipitate_radius_nm)


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        out = csv.DictWriter(stream, fieldnames=list(rows[0])); out.writeheader(); out.writerows(rows)


def main(args):
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    args.output.mkdir(parents=True, exist_ok=True)
    source = json.loads(args.source.read_text(encoding='utf-8'))
    e_si, e_ox, compression = (source['elastic_inputs'][k]
        for k in ('Si_axial_E_GPa', 'oxide_E_GPa', 'oxide_compressive_input_GPa'))
    arithmetic, geometry_rows, stress_rows = [], [], []
    for row in source['rows']:
        t = row['oxide_nominal_nm']*1e-3
        w,h = row['outer_width_um'],row['outer_height_um']
        arithmetic.append(dict(group=row['group'],
            nominal_strength_change_percent=(row['nominal_strength_GPa']/source['rows'][0]['nominal_strength_GPa']-1)*100,
            inferred_flaw_half_length_nm=inverse_equivalent_flaw_nm(row['FEM_Si_axial_GPa']),
            reported_flaw_half_length_nm=row['reported_flaw_half_length_nm'],
            precipitate_radius_nm=(precipitate_radius_nm(row['oxidation_time_min'])
                if row['reported_precipitate_radius_nm'] is not None else None),
            reported_precipitate_radius_nm=row['reported_precipitate_radius_nm']))
        if t == 0: continue
        # These are competing interpretations/sensitivity cases, not corrections.
        conventions = [
            ('listed_outer_per_face_t', w,h,t),
            ('listed_outer_total_t', w,h,t/2),
            ('planar_44percent_consumption_per_face_t', 4+2*.56*t,5+2*.56*t,t)]
        for name,wo,ho,face_t in conventions:
            areas = rectangular_shell(outer_width=wo,outer_height=ho,thickness_per_face=face_t)
            zero = axial_composite(areas=areas,young_moduli=[e_si,e_ox],eigenstrains=[0,0],nominal_stress=0)
            geometry_rows.append(dict(group=row['group'],convention=name,
                outer_width_um=wo,outer_height_um=ho,thickness_per_face_um=face_t,
                core_width_um=wo-2*face_t,core_height_um=ho-2*face_t,
                Si_area_um2=areas[0],oxide_area_um2=areas[1],oxide_area_fraction=areas[1]/sum(areas),
                modulus_GPa=zero['effective_modulus'],source_modulus_GPa=row['E_theory_GPa'],
                modulus_difference_GPa=zero['effective_modulus']-row['E_theory_GPa']))
            # Hold inferred experimental force fixed when changing cross section.
            force_source = row['nominal_strength_GPa']*w*h
            nominal = force_source/(wo*ho)
            for residual in ('none', 'oxide_eigenstress_at_zero_common_strain', 'self_balanced_actual_oxide_stress'):
                eigen = ([0,0] if residual == 'none' else [0,compression/e_ox])
                if residual == 'self_balanced_actual_oxide_stress':
                    eigen = [-compression*areas[1]/(areas[0]*e_si), compression/e_ox]
                result = axial_composite(areas=areas,young_moduli=[e_si,e_ox],eigenstrains=eigen,nominal_stress=nominal)
                si,ox = result['phase_stress']
                stress_rows.append(dict(group=row['group'],convention=name,residual_model=residual,
                    applied_nominal_GPa=nominal,Si_axial_GPa=si,oxide_axial_GPa=ox,
                    source_FEM_Si_GPa=row['FEM_Si_axial_GPa'],source_FEM_oxide_surface_GPa=row['FEM_oxide_axial_GPa'],
                    Si_difference_percent=(si/row['FEM_Si_axial_GPa']-1)*100,
                    oxide_difference_GPa=ox-row['FEM_oxide_axial_GPa'],
                    axial_strain=result['strain'],force_residual_GPa_um2=result['force_residual']))
    write_csv(args.output/'paper_arithmetic.csv', arithmetic)
    write_csv(args.output/'geometry_conventions.csv', geometry_rows)
    write_csv(args.output/'axial_stress_partition.csv', stress_rows)
    modulus_match = [r for r in geometry_rows if r['convention']=='listed_outer_per_face_t']
    balance = max(abs(r['force_residual_GPa_um2']) for r in stress_rows)
    si_match = [r for r in stress_rows if r['convention']=='listed_outer_per_face_t'
                and r['residual_model']=='oxide_eigenstress_at_zero_common_strain']
    summary = dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        arithmetic_rows=len(arithmetic), geometry_rows=len(geometry_rows), stress_rows=len(stress_rows),
        max_force_balance_residual_GPa_um2=balance,
        max_modulus_difference_listed_per_face_GPa=max(abs(r['modulus_difference_GPa']) for r in modulus_match),
        max_Si_point_stress_difference_percent=max(abs(r['Si_difference_percent']) for r in si_match),
        oxide_removed_width_vs_perface200_core_difference_um=3.91-(4.11-.4),
        oxide_removed_height_vs_perface200_core_difference_um=4.91-(5.11-.4),
        scope='1D common-strain audit; source strength is input; no independent strength prediction',
        limitation='FEM stresses are point samples, while this model supplies uniform phase stresses; no 3D field reproduction',
        actual_new_DFT=0, actual_new_MD=0, calibrated_initiation_probability=None)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--source',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
