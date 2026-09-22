"""Reproduce the Durham B-2025 source audit; no DFT or concentration fitting.

Requires numpy/openpyxl. Original workbook/ZIP bytes are never changed.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
from zipfile import ZipFile
import numpy as np

ZIP_SHA256 = '5ed01b79a17bda9f59cd680bc965d8679bf25151a01d1e8c0ea81da4767fa9dd'
DATASET_DOI = '10.15128/r1mw22v5500'
PAPER_DOI = '10.1088/1361-648X/ae191e'


def parse_cell(text):
    """The archived orthogonal CASTEP input format, with explicit units check."""
    def block(name):
        match = re.search(r'%block\s+'+name+r'\s*\n(.*?)%endblock', text, re.I|re.S)
        if match is None:
            raise ValueError('missing '+name)
        return [line.split() for line in match.group(1).splitlines() if line.strip()]
    abc = np.asarray(block('lattice_abc'), float)
    if abc.shape != (2, 3) or not np.allclose(abc[1], 90):
        raise ValueError('this source reader only supports orthogonal default-Angstrom cells')
    rows = block('positions_frac')
    symbols = [r[0] for r in rows]
    frac = np.asarray([r[1:4] for r in rows], float)
    if any(len(r) != 4 for r in rows) or frac.shape != (len(symbols), 3):
        raise ValueError('unexpected fractional position layout')
    return dict(symbols=symbols, fractional_positions=frac.tolist(),
                cell_A=np.diag(abc[0]).tolist(),
                fix_all_cell=bool(re.search(r'FIX_ALL_CELL\s*:\s*TRUE', text, re.I)))


def save_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def run(archive, output):
    import openpyxl
    raw = archive.read_bytes()
    if hashlib.sha256(raw).hexdigest() != ZIP_SHA256:
        raise ValueError('source ZIP differs from audited Durham download')
    output.mkdir(parents=True, exist_ok=True)
    z = ZipFile(BytesIO(raw))
    structures, geometry = {}, []
    members = []
    for name in sorted(z.namelist()):
        if name.endswith('/'):
            continue
        b = z.read(name)
        members.append(dict(member=name, size_bytes=len(b), sha256=hashlib.sha256(b).hexdigest()))
        if not name.endswith('.cell'):
            continue
        s = parse_cell(b.decode('utf-8-sig'))
        key = Path(name).stem
        structures[key] = s
        count = Counter(s['symbols'])
        cell = np.asarray(s['cell_A']); frac = np.asarray(s['fractional_positions'])
        volume = float(np.linalg.det(cell))
        bfrac = frac[np.array(s['symbols']) == 'B']
        delta = bfrac[:, None, :]-bfrac[None, :, :]
        # Orthogonal cell => componentwise wrapping is the exact minimum image.
        distances = np.linalg.norm((delta-np.rint(delta))@cell, axis=-1)
        pairs = distances[np.triu_indices(len(bfrac), 1)].tolist()
        geometry.append(dict(structure=key, Si_atoms=count['Si'], B_atoms=count['B'],
            total_atoms=len(frac), volume_A3=volume, chemical_B_per_cm3=count['B']/volume*1e24,
            B_atomic_fraction=count['B']/len(frac), B_pair_distances_A=json.dumps(pairs),
            fixed_cell_input=s['fix_all_cell'], status='archived input; not relaxed output'))
    save_csv(output/'structure_audit.csv', geometry)
    (output/'source_structures.json').write_text(json.dumps(structures, indent=2)+'\n', encoding='utf-8')

    hall_name = 'Raw data files/Hall and SIMS data/Hall_data_Liverpool.xlsx'
    wb = openpyxl.load_workbook(BytesIO(z.read(hall_name)), data_only=True)
    wf = openpyxl.load_workbook(BytesIO(z.read(hall_name)), data_only=False)
    e = 1.602176634e-19  # exact SI charge; Rh units cm^3/C, rho Ohm cm
    hall = []
    for rownum, row in enumerate(wb['BDS raw data'].iter_rows(min_row=2, values_only=True), 2):
        comment, sample, temperature, thickness, rho, rh, n, mobility = row
        if rho is None:
            continue
        hall.append(dict(sheet='BDS raw data', row=rownum, sample=sample, comment=comment,
            temperature_K=temperature, thickness_nm=thickness, resistivity_Ohm_cm=rho,
            Hall_Rh_cm3_C=rh, reported_Hall_n_cm3=n, recomputed_Hall_n_cm3=1/(e*rh),
            reported_Hall_mobility_cm2_Vs=mobility, recomputed_Hall_mobility_cm2_Vs=rh/rho,
            concentration_relative_rounding_difference=(1/(e*rh)/n-1),
            mobility_relative_rounding_difference=(rh/rho/mobility-1)))
    save_csv(output/'hall_raw_audit.csv', hall)
    temp = []
    for rownum, row in enumerate(wb['HDS Temperature data'].iter_rows(min_row=2, values_only=True), 2):
        t, rho, voltage, polarity, n, mobility, rh = row
        if t is None:
            continue
        temp.append(dict(sheet='HDS Temperature data', row=rownum, temperature_K=t,
            resistivity_Ohm_cm=rho, Hall_voltage_source_value=voltage, polarity=polarity,
            reported_Hall_n_cm3=n, recomputed_Hall_n_cm3=1/(e*rh),
            reported_Hall_mobility_cm2_Vs=mobility, recomputed_Hall_mobility_cm2_Vs=rh/rho,
            Hall_Rh_cm3_C=rh, concentration_relative_difference=1/(e*rh)/n-1,
            mobility_relative_difference=rh/rho/mobility-1))
    save_csv(output/'hall_temperature_audit.csv', temp)
    aggregates = []
    # Raw IDs carry HDS/LDS names. Do not regard repeats as independent wafers.
    for sample in sorted(set(r['sample'] for r in hall)):
        selected = [r for r in hall if r['sample'] == sample]
        d = dict(sample=sample, measurements=len(selected), independent_wafer_count=None)
        for field in ['reported_Hall_n_cm3', 'resistivity_Ohm_cm', 'reported_Hall_mobility_cm2_Vs']:
            values = [r[field] for r in selected]
            d[field+'_mean'] = float(np.mean(values))
            d[field+'_sample_SD'] = float(np.std(values, ddof=1))
        aggregates.append(d)
    save_csv(output/'hall_recomputed_summaries.csv', aggregates)
    broken = [dict(sheet=sheet.title, cell=c.coordinate, cached_value=c.value,
                   formula=wf[sheet.title][c.coordinate].value)
              for sheet in wb for row in sheet for c in row if c.data_type == 'e']

    sims_name = 'Raw data files/Hall and SIMS data/Boron in Si wafers for DU.xlsx'
    ws = openpyxl.load_workbook(BytesIO(z.read(sims_name)), data_only=True).active
    sims, sims_summary = [], []
    for label, first in [('Agar_nominal_1e15_1e16', 2), ('DU_nominal_1e19_1e20', 7), ('DU_undoped', 12)]:
        sample_rows = []
        for rownum in range(3, ws.max_row+1):
            values = [ws.cell(rownum, first+j).value for j in range(3)]
            if all(v is None for v in values):
                continue
            if not all(isinstance(v, (float, int)) for v in values):
                raise ValueError('unexpected SIMS numeric row')
            sample_rows.append(dict(sample=label, sheet=ws.title, row=rownum,
                time_msec=values[0], B11_counts=values[1], Si29_counts=values[2]))
        sims.extend(sample_rows)
        sims_summary.append(dict(sample=label, rows=len(sample_rows),
            B11_count_mean=float(np.mean([r['B11_counts'] for r in sample_rows])),
            Si29_count_mean=float(np.mean([r['Si29_counts'] for r in sample_rows])),
            calibrated_chemical_concentration_cm3=None,
            ionization_fraction=None, status='counts only; no sensitivity factor or depth calibration'))
    save_csv(output/'sims_raw_counts.csv', sims)
    save_csv(output/'sims_counts_summary.csv', sims_summary)
    report = dict(dataset_DOI=DATASET_DOI, paper_DOI=PAPER_DOI, license='CC BY 4.0',
        zip_sha256=ZIP_SHA256, zip_md5=hashlib.md5(raw).hexdigest(),
        repository_md5='8d0c05222168c004374692040e91bc5f',
        structure_count=len(structures), source_members=members, broken_source_cells=broken,
        Hall_max_raw_n_relative_rounding_difference=max(abs(r['concentration_relative_rounding_difference']) for r in hall),
        Hall_max_raw_mobility_relative_rounding_difference=max(abs(r['mobility_relative_rounding_difference']) for r in hall),
        Hall_temperature_n_range_cm3=[min(r['reported_Hall_n_cm3'] for r in temp), max(r['reported_Hall_n_cm3'] for r in temp)],
        Hall_temperature_n_identity_max_relative_difference=max(abs(r['concentration_relative_difference']) for r in temp),
        Hall_temperature_implied_e_C_median=float(np.median([1/(r['Hall_Rh_cm3_C']*r['reported_Hall_n_cm3']) for r in temp])),
        Hall_temperature_implied_e_C_peak_to_peak=float(np.ptp([1/(r['Hall_Rh_cm3_C']*r['reported_Hall_n_cm3']) for r in temp])),
        Hall_temperature_mobility_identity_max_relative_difference=max(abs(r['mobility_relative_difference']) for r in temp),
        distinctions=[
            'Hall number assumes the source Hall model; Hall factor and compensation are not independently resolved.',
            'Hall mobility is carrier transport, not a/s collective-coordinate mobility.',
            'SIMS counts are not calibrated atomic concentration; p/C_B and ionization are unavailable.',
            'Archived cell files are input geometries, with 64 host sites and 64-67 total atoms.',
            'All archived inputs set FIX_ALL_CELL TRUE; the paper describes variable cells for some calculations.',
            'No relaxed DFT output energy, force, Hessian or charge parameter file is present in this ZIP.',
            'Chemical arrangements are quenched labels unless independently justified to exchange.'])
    if report['zip_md5'] != report['repository_md5']:
        raise ValueError('repository MD5 does not match')
    (output/'source_audit.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.archive, args.output)
    print(json.dumps({k:v for k,v in result.items() if k != 'source_members'}, indent=2))
