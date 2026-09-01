"""PC-side helper for fatigue-tester telemetry.

The real-time force loop belongs on the MCU. This script is for logging,
visualization, and later coupling measured stress/strain histories into the
research theory solver.

Expected CSV columns:
time_s,cycle,stress_ref_pa,force_ref_n,force_n,displacement_m,strain,
temperature_c,dcpd_v,actuator_command,fault_flags
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "time_s",
    "cycle",
    "stress_ref_pa",
    "force_ref_n",
    "force_n",
    "displacement_m",
    "strain",
    "temperature_c",
    "dcpd_v",
    "actuator_command",
    "fault_flags",
]


def validate_header(fieldnames):
    missing = [name for name in FIELDS if name not in fieldnames]
    if missing:
        raise ValueError(f"missing telemetry fields: {missing}")


def summarize(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        validate_header(reader.fieldnames or [])
        rows = list(reader)

    if not rows:
        print("no telemetry rows")
        return

    last = rows[-1]
    print(f"rows: {len(rows)}")
    print(f"last cycle: {last['cycle']}")
    print(f"last force: {last['force_n']} N")
    print(f"last strain: {last['strain']}")
    print(f"last temperature: {last['temperature_c']} C")
    print(f"fault flags: {last['fault_flags']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", type=Path)
    args = parser.parse_args()
    summarize(args.csv_file)


if __name__ == "__main__":
    main()
