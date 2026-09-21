"""Evaluate published Si DFT configurations with SW and ordinary Tersoff.

This runs new classical potential evaluations on archived DFT structures.
It does not run DFT or GAP, fit parameters, or certify material/kinetic gates.
Different declared DFT functionals and missing metadata stay separate.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import time
from urllib.request import urlopen
import zipfile

import numpy as np

from solver_v1.silicon_atomistic_reference import (LammpsSilicon, force_metrics,
                                                 periodic_basis_audit)
from .run_atomistic_controls import sources, checksum, save_json


ARCHIVE_URL = 'https://www.repository.cam.ac.uk/bitstream/handle/1810/317974/Si_PRX_GAP.zip'
ARCHIVE_SHA256 = '1d3efe976c53bcd2889c0556f2d600c26bf4603522b5c92fd8633dff7b7ab5c2'
DATA_MEMBER = 'gp_iter6_sparse9k.xml.xyz'
DATA_SHA256 = '74ac5f387aa0149ec9ea8f3e275d1087d7287db68c15ba8bc4c28a450de67e2a'


def read_frames(archive, *, download=False):
    from ase.io import iread
    archive = Path(archive)
    if not archive.exists():
        if not download:
            raise FileNotFoundError('source archive missing; use --download to obtain the hash-bound source')
        archive.parent.mkdir(parents=True, exist_ok=True)
        partial = archive.with_suffix('.download')
        with urlopen(ARCHIVE_URL, timeout=60) as source, partial.open('xb') as target:
            while chunk := source.read(1024*1024):
                target.write(chunk)
        if checksum(partial) != ARCHIVE_SHA256:
            raise ValueError('downloaded archive checksum differs; partial file preserved')
        partial.replace(archive)
    if checksum(archive) != ARCHIVE_SHA256:
        raise ValueError('archive does not match the audited Cambridge source')
    with zipfile.ZipFile(archive) as z:
        bad = z.testzip()
        if bad:
            raise ValueError(f'archive CRC failed: {bad}')
        raw = z.read(DATA_MEMBER)
        members = [dict(name=m.filename, bytes=m.file_size, crc32=f'{m.CRC:08x}')
                   for m in z.infolist()]
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise ValueError('DFT member checksum differs')
    return list(iread(io.StringIO(raw.decode('utf-8')), format='extxyz')), members


def reference_fields(atoms):
    prefixes = [p for p in ('dft', 'DFT') if p+'_energy' in atoms.info]
    if len(prefixes) != 1:
        raise ValueError('one unambiguous source DFT energy field required')
    p = prefixes[0]
    if p+'_force' not in atoms.arrays or not np.all(atoms.numbers == 14):
        raise ValueError('pure Si with matching DFT force field required')
    energy, force = float(atoms.info[p+'_energy']), np.asarray(atoms.arrays[p+'_force'], float)
    if not np.isfinite(energy) or force.shape != atoms.positions.shape or not np.all(np.isfinite(force)):
        raise ValueError('nonfinite/malformed source DFT data')
    return energy, force, p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('.cache/si-atomistic-v5/Si_PRX_GAP.zip'))
    parser.add_argument('--output', type=Path, default=Path('results/silicon_atomistic_v5'))
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if (args.output/'dft_material_audit.json').exists():
        raise FileExistsError('refusing to replace completed DFT material audit')
    started = time.perf_counter()
    frames, members = read_frames(args.archive, download=args.download)
    refs = [reference_fields(a) for a in frames]
    # A fixed, explicitly NONRELAXED same-geometry anchor. It only removes the
    # per-atom energy zero between engines; it is not a fitted bulk ground state.
    eligible = [i for i, a in enumerate(frames) if a.info['config_type'] == 'dia'
                and a.info.get('xc_functional') == 'PW91']
    anchor = min(eligible, key=lambda i:refs[i][0]/len(frames[i]))
    specs, _ = sources(args.output)
    all_rows, metadata = [], {}
    offsets = np.r_[0, np.cumsum([len(a) for a in frames])]
    saved = dict(offsets=offsets, reference_force=np.vstack([f for _, f, _ in refs]),
                 reference_energy=np.array([e for e, _, _ in refs]))
    transforms = []
    for label, (style, potential, source_hash) in specs.items():
        predicted_forces, predicted_energies = [], []
        with LammpsSilicon(style, potential, sha256=source_hash) as engine:
            metadata[label] = engine.metadata()
            anchor_cell, _ = periodic_basis_audit(frames[anchor].cell.array, frames[anchor].pbc)
            ep0 = engine.evaluate(frames[anchor].positions, anchor_cell).energy/len(frames[anchor])
            ed0 = refs[anchor][0]/len(frames[anchor])
            for i, (atoms, (ed, fd, field)) in enumerate(zip(frames, refs)):
                cell, signs = periodic_basis_audit(atoms.cell.array, atoms.pbc)
                value = engine.evaluate(atoms.positions, cell, periodic=atoms.pbc)
                fp = -value.gradient
                xc = str(atoms.info.get('xc_functional', 'UNSPECIFIED'))
                metric = force_metrics(fp, fd)
                row = dict(model=label, frame=i, atoms=len(atoms), config_type=str(atoms.info['config_type']),
                           declared_xc=xc, dft_field_prefix=field, dft_energy_eV=ed,
                           model_energy_eV=value.energy, basis_signs=' '.join(map(str, signs)),
                           volume_A3=float(abs(np.linalg.det(cell))), **metric)
                # Do not compare PBE or absent per-frame XC metadata with a PW91
                # energy anchor. Their raw energies and force comparisons remain.
                row['diamond_referenced_error_eV_atom'] = (
                    value.energy/len(atoms)-ep0-(ed/len(atoms)-ed0) if xc == 'PW91' else None)
                all_rows.append(row)
                predicted_forces.append(fp)
                predicted_energies.append(value.energy)
                if label == next(iter(specs)) and np.any(signs != 1):
                    transforms.append(dict(frame=i, signs=signs.tolist(), operation='integer periodic basis only; positions unchanged'))
                if (i+1) % 250 == 0:
                    print(json.dumps(dict(model=label, evaluated=i+1, total=len(frames))), flush=True)
        saved[label+'_force'] = np.vstack(predicted_forces)
        saved[label+'_energy'] = np.asarray(predicted_energies)
    groups = defaultdict(list)
    for r in all_rows:
        groups[(r['model'], r['config_type'], r['declared_xc'])].append(r)
    summaries = []
    for (model, kind, xc), records in sorted(groups.items()):
        ncomp = sum(r['components'] for r in records)
        de = [r['diamond_referenced_error_eV_atom'] for r in records
              if r['diamond_referenced_error_eV_atom'] is not None]
        summaries.append(dict(model=model, config_type=kind, declared_xc=xc, frames=len(records),
            atoms=sum(r['atoms'] for r in records),
            component_rmse_eV_A=float(np.sqrt(sum(r['squared_error'] for r in records)/ncomp)),
            frame_equal_weight_component_rmse_eV_A=float(np.sqrt(np.mean([r['component_rmse_eV_A']**2 for r in records]))),
            max_vector_error_eV_A=max(r['max_vector_error_eV_A'] for r in records),
            reference_component_rms_eV_A=float(np.sqrt(sum(r['reference_component_rms_eV_A']**2*r['components'] for r in records)/ncomp)),
            diamond_referenced_energy_rmse_meV_atom=float(1000*np.sqrt(np.mean(np.square(de)))) if de else None,
            diamond_referenced_energy_mean_meV_atom=float(1000*np.mean(de)) if de else None))
    with (args.output/'dft_frame_metrics.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    np.savez_compressed(args.output/'dft_predictions.npz', **saved)
    import ase
    result = dict(scope='new classical run-0 evaluations against published DFT; no new DFT/GAP or fitting',
        source=dict(url=ARCHIVE_URL, doi='10.17863/CAM.65004', archive_sha256=ARCHIVE_SHA256,
                    member=DATA_MEMBER, member_sha256=DATA_SHA256, archive_members=members,
                    paper_doi='10.1103/PhysRevX.8.041048', data_role='GAP training data; not held out for GAP'),
        source_frames=len(frames), source_atoms=int(offsets[-1]), evaluated_pairs=len(all_rows),
        declared_xc_counts=dict(Counter(a.info.get('xc_functional', 'UNSPECIFIED') for a in frames)),
        source_DFT_field_counts=dict(Counter(p for _, _, p in refs)),
        source_stored_GAP_predictions=sum('gap_energy' in a.info for a in frames),
        stored_GAP_predictions_used=False,
        frame_types=dict(Counter(a.info['config_type'] for a in frames)),
        energy_anchor=dict(frame=anchor, source='lowest source E/N among PW91 dia frames; no fitting',
                           atoms=len(frames[anchor]), dft_eV_atom=refs[anchor][0]/len(frames[anchor]),
                           force_max_eV_A=float(np.max(abs(refs[anchor][1]))),
                           is_relaxed_bulk_reference=False,
                           scope='only explicit PW91 frames; energy zero, not a cohesive-energy calibration'),
        cell_basis_changes=transforms, model_metadata=metadata, group_results=summaries,
        ase_version=ase.__version__, elapsed_seconds=time.perf_counter()-started,
        code_sha256={str(Path(p).resolve().relative_to(Path.cwd())).replace('\\', '/'):checksum(p)
                     for p in [__file__, 'solver_v1/silicon_atomistic_reference.py']},
        complete=True, material_calibrated=False, kinetic_calibrated=False)
    save_json(args.output/'dft_material_audit.json', result)
    print(json.dumps(dict(event='complete', frames=len(frames), pairs=len(all_rows),
                         elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
