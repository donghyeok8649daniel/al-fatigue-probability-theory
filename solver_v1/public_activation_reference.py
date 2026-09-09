"""Published apparent activation data: a kinetic *derivative*, not A_c.

Reads cached values of the source XLSX; never evaluates formulas/macros or
executes the accompanying notebooks. Apparent activation volume is inferred
from rate sensitivity, not a geometric volume invented for the atomistic PDE.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET
import zipfile

import numpy as np

from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
BOOK = 'Data in brief/Final results/Overview results and final selection stress relaxations.xlsx'


def cached_workbook(raw):
    """Small read-only OOXML value reader; formulas remain provenance strings."""
    archive = zipfile.ZipFile(io.BytesIO(raw))
    strings = [''.join(node.itertext()) for node in
               ET.fromstring(archive.read('xl/sharedStrings.xml'))]
    relationships = {r.attrib['Id']: r.attrib['Target'] for r in
                     ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
    result = {}
    for sheet in ET.fromstring(archive.read('xl/workbook.xml')).findall('s:sheets/s:sheet', NS):
        target = relationships[sheet.attrib[REL]]
        target = target.lstrip('/') if target.startswith('/') else 'xl/'+target
        cells = {}
        for cell in ET.fromstring(archive.read(target)).findall('.//s:sheetData/s:row/s:c', NS):
            value = cell.find('s:v', NS)
            if value is None:
                continue
            kind = cell.attrib.get('t')
            data = strings[int(value.text)] if kind == 's' else value.text
            if kind not in ('s', 'str', 'e'):
                data = float(data)
            formula = cell.find('s:f', NS)
            cells[cell.attrib['r']] = dict(value=data,
                cached_formula=None if formula is None else formula.text)
        result[sheet.attrib['name']] = cells
    return result


def parse_selected_activation(raw):
    lines = [s.strip() for s in raw.decode('cp1252').splitlines() if s.strip()]
    header = tuple(lines[0].split('\t'))
    if (header[:5] != ('# Relax', 'Eff stress', 'Eff strain', 'a [b²]', 'V [b^3]')
            or header[6] != 'b²/a' or len(header) != 8):
        raise ValueError('unrecognized apparent-activation source units')
    rows = []
    for line in lines[1:]:
        fields = line.split('\t')
        if len(fields) != 8 or re.fullmatch(r'relax\s*\d+', fields[0], re.I) is None:
            raise ValueError('source relaxation label or column count changed')
        values = np.array(list(map(float, fields[1:])))
        if np.any(~np.isfinite(values)):
            raise ValueError('finite source values required; do not silently repair')
        stress, strain, area, volume, error, inverse, inverse_error = values
        if area <= 0 or volume <= 0 or error < 0:
            raise ValueError('positive source apparent activation and nonnegative area error required')
        rows.append(dict(hold=int(re.search(r'\d+', fields[0])[0]),
            source_effective_stress=stress, source_effective_strain=strain,
            apparent_area_over_b2=area, apparent_volume_over_b3=volume,
            source_area_standard_error_over_b2=error,
            source_inverse_area=inverse, source_signed_inverse_error=inverse_error,
            volume_area_identity_error=float(volume-area)))
    return rows


def activation_unit_conversion(area_over_b2, volume_over_b3, *, burgers_m,
                               temperature_K, boltzmann_J_K=1.380649e-23):
    """a_app=kT/b d ln(rate)/d tau, V_app=b a_app; units are SI.

    This derivative is not a mobility, event count, chosen source size or
    specimen independent-region area. A stress-dependent prefactor or changing
    source population contributes to an experimental *apparent* derivative.
    """
    values = np.array([area_over_b2, volume_over_b3, burgers_m,
                       temperature_K, boltzmann_J_K], float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError('positive finite source values with explicit units required')
    area, volume, b, temperature, kb = values
    return dict(apparent_activation_area_m2=float(area*b*b),
        apparent_activation_volume_m3=float(volume*b**3),
        rate_log_slope_per_MPa=float(volume*b**3/(kb*temperature)*1e6),
        specimen_correlation_area=None, cell_mobility=None, t0_seconds=None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('preserve previous source audits')
    manifest = json.loads((args.source/'manifest.json').read_bytes())
    raw = (args.source/'supplementary.zip').read_bytes()
    if hashlib.sha256(raw).hexdigest() != manifest['supplementary_sha256']:
        raise ValueError('source checksum mismatch')
    outer = zipfile.ZipFile(io.BytesIO(raw))
    archive = zipfile.ZipFile(io.BytesIO(outer.read('mmc2.zip')))
    book_raw = archive.read(BOOK); book = cached_workbook(book_raw)
    numbers = book['Numbers']
    if [numbers[key]['value'] for key in ('B4','D4','B5','D5','B6','D6')] != [
            'b','nm','k','*10^-23 m² kg/s²K','T','K']:
        raise ValueError('source workbook unit declarations changed')
    b = numbers['C4']['value']*1e-9
    kb = numbers['C5']['value']*1e-23
    temperature = numbers['C6']['value']
    tables = {r[0].casefold(): (t['label'], r[1]) for t in manifest['tables']
              if t['label'] in ('Table 1','Table 2') for r in t['rows'][1:]}
    rows = []
    for path in sorted(archive.namelist()):
        if '/Final results/' not in path or not path.endswith('.txt'):
            continue
        data = archive.read(path); sample = Path(path).stem
        condition = 'annealed_500C_2h' if '/Annealed/' in path else 'as_cast'
        for record in parse_selected_activation(data):
            rows.append(dict(sample=sample, condition=condition,
                source_member=path, source_member_sha256=hashlib.sha256(data).hexdigest(),
                source_table_matched=sample.casefold() in tables,
                source_marked_final_selected=True,
                source_analysis_temperature_K=temperature, source_analysis_burgers_m=b,
                **record, **activation_unit_conversion(record['apparent_area_over_b2'],
                    record['apparent_volume_over_b3'], burgers_m=b,
                    temperature_K=temperature, boltzmann_J_K=kb)))
    write_csv(args.out/'source_selected_activation.csv', rows)
    # Explicitly verify the representative specimen's source sheet, without
    # running formulas or assuming its processed stress equals raw tau(t=0).
    checked = []
    for row in (r for r in rows if r['sample']=='Al_20_r_10'):
        cells=book['Al_20_r_10']; index=row['hold']+4
        slope=cells[f'D{index}']['value']; initial=cells[f'B{index}']['value']
        checked.append(dict(sample=row['sample'], hold=row['hold'],
            workbook_tau0_MPa=initial,
            final_effective_stress_literal=row['source_effective_stress'],
            processed_minus_workbook_stress=row['source_effective_stress']-initial,
            workbook_slope_per_MPa=slope,
            replayed_slope_per_MPa=row['rate_log_slope_per_MPa'],
            relative_slope_replay_error=row['rate_log_slope_per_MPa']/slope-1,
            workbook_area_over_b2=cells[f'K{index}']['value'],
            final_area_over_b2=row['apparent_area_over_b2'],
            cached_formula=cells[f'K{index}']['cached_formula']))
    write_csv(args.out/'representative_workbook_replay.csv', checked)
    save_json(args.out/'source_scope.json',dict(completed=True,doi=manifest['doi'],
        source_sha256=manifest['supplementary_sha256'],workbook_sha256=hashlib.sha256(book_raw).hexdigest(),
        source_constants={'burgers_m':b,'temperature_K':temperature,'boltzmann_J_K':kb},
        source_constant_status='published workbook analysis constants, not new thermometry/lattice measurement',
        selected_rows=len(rows),samples=len(set(r['sample'] for r in rows)),
        unmatched_samples=sorted(set(r['sample'] for r in rows if not r['source_table_matched'])),
        minimum_volume_over_b3=min(r['apparent_volume_over_b3'] for r in rows),
        maximum_volume_over_b3=max(r['apparent_volume_over_b3'] for r in rows),
        negative_source_inverse_errors=sum(r['source_signed_inverse_error']<0 for r in rows),
        signed_source_errors_not_silently_repaired=True,
        final_effective_stress_not_identified_with_raw_tau=True,
        formulas_executed=False,source_selection_not_independently_rederived=True,
        model_activation_volume_fitted=False,collective_coordinate_mobility_calibrated=False))
    print(f'{len(rows)} source-selected apparent-activation observations; no clock calibration')


if __name__ == '__main__':
    main()
