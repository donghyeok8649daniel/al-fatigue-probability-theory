"""Compact only this study's generated outputs, without numerical rounding.

History CSV compression is byte-for-byte reversible. Complete legacy density
snapshots are preserved under ignored .cache before selecting explicit cycle
indices for the small checked-in result set. Per-cycle statistics never change.
No production inputs, source parameters, or unrelated user files are touched.
"""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from .run_low_stress_cyclic_diagnostic import ROOT, OUT, save_json


def compact():
    records = []
    results_root = (OUT/"cycles").resolve()
    backup_root = (ROOT/".cache/low-stress-full-snapshots").resolve()
    if not backup_root.is_relative_to(ROOT.resolve()):
        raise ValueError("backup must stay inside this worktree")
    for audit_path in sorted(results_root.glob("*/audit.json")):
        directory = audit_path.parent.resolve()
        if not directory.is_relative_to(results_root):
            raise ValueError("refuse unrelated output directory")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        history = directory/"history.csv"
        if history.exists():
            raw = history.read_bytes()
            target = directory/"history.csv.gz"
            compressed = gzip.compress(raw, mtime=0)
            if target.exists() and gzip.decompress(target.read_bytes()) != raw:
                raise ValueError(f"different compressed history already exists: {directory.name}")
            target.write_bytes(compressed)
            if gzip.decompress(target.read_bytes()) != raw:
                raise ValueError("lossless compression verification failed")
            history.unlink()  # Exact bytes remain recoverable from history.csv.gz.
        archive = directory/"snapshots.npz"
        with np.load(archive) as original:
            if "cycle_indices" in original:
                continue
            mass = original["same_phase"].copy()
        with (directory/"per_cycle.csv").open(encoding="utf-8") as stream:
            import csv
            cycle_rows = list(csv.DictReader(stream))
        driven = sum(r["segment"] == "cyclic" for r in cycle_rows)
        total = len(cycle_rows)
        indices = sorted({0, 1, driven-1, driven, min(driven+1, total), total})
        backup_root.mkdir(parents=True, exist_ok=True)
        backup = backup_root/(directory.name+".npz")
        original_bytes = archive.read_bytes()
        if backup.exists() and backup.read_bytes() != original_bytes:
            raise ValueError("refuse to overwrite a different full-snapshot backup")
        if not backup.exists():
            backup.write_bytes(original_bytes)
        np.savez_compressed(archive, cycle_indices=indices, same_phase=mass[indices])
        with np.load(archive) as checked:
            np.testing.assert_array_equal(checked["same_phase"], mass[indices])
        audit["snapshot_cycle_indices"] = indices
        audit["full_snapshot_backup"] = backup.relative_to(ROOT).as_posix()
        audit["full_snapshot_sha256"] = hashlib.sha256(original_bytes).hexdigest()
        save_json(audit_path, audit)
        records.append(dict(run=directory.name, original_bytes=len(original_bytes),
            sampled_bytes=archive.stat().st_size, cycle_indices=indices,
            backup=backup.relative_to(ROOT).as_posix()))
    # Additive record: repeated invocations do not erase prior backup provenance.
    manifest = OUT/"output_compaction.json"
    old = json.loads(manifest.read_text()) if manifest.exists() else []
    save_json(manifest, old+records)
    print(f"Backed up and sampled {len(records)} complete snapshot archives.")


if __name__ == "__main__":
    compact()
