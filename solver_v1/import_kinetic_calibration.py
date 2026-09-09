"""Import actual a/s covariance data into a source-bound clock.

Usage: python -m solver_v1.import_kinetic_calibration metadata.json covariance.csv
       --output user_kinetics.json --diagnostics fit_checks.json

CSV: time_seconds,C_aa_m2,C_as_m2,C_sa_m2,C_ss_m2 (includes lag zero).
Metadata: model_id,length_scale_m,energy_scale_J,temperature_K,source,
relative_tolerance,notes. Output is never the repository's default file.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .kinetic_calibration_workflow import calibration_from_collective_correlations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("covariance", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagnostics", type=Path, required=True)
    args = parser.parse_args()
    if (args.output.resolve() == args.diagnostics.resolve()
            or args.output.exists() or args.diagnostics.exists()):
        parser.error("use distinct NEW output paths; existing calibration files are not overwritten")
    meta = json.loads(args.metadata.read_text(encoding="utf-8"))
    with args.covariance.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    t = [float(row["time_seconds"]) for row in rows]
    C = np.array([[float(row[key]) for key in
                   ("C_aa_m2", "C_as_m2", "C_sa_m2", "C_ss_m2")] for row in rows]).reshape(-1, 2, 2)
    calibration, fit = calibration_from_collective_correlations(t, C, **meta)
    for path, data in ((args.output, calibration.to_dict()), (args.diagnostics, fit)):
        with path.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write("\n")
    print("Correlation consistency passed. Source/coordinate provenance still requires scientific review.")


if __name__ == "__main__":
    main()
