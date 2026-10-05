"""Loss-aware MCU telemetry to CAD measurement handoff; no solver or device I/O.

The source contract is fatigue_tester/pc/fatigue_pc_bridge.py at research
commit 399ffb1f9d711439c4fbed395070dfc40050cc3c. The CAD measurement columns
are the published CAD research handoff at commit 6a81c64c08bce6f74e2723278fc8efecd9b03a37.
Laboratory seconds stay laboratory seconds; this never maps the PDE clock.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, InvalidOperation, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile


CAD_COLUMNS = (
    "time_s", "cycle_index", "force_N", "crosshead_displacement_mm",
    "gauge_strain", "temperature_C", "source", "quality_flag",
)
REQUIRED_COLUMNS = ("time_s", "cycle", "force_n", "displacement_m", "temperature_c")
OPTIONAL_CHANNELS = ("stress_ref_pa", "force_ref_n", "strain", "dcpd_v", "actuator_command", "fault_flags")
FORMAT = "aft.cad-measurements/1"


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _decimal(raw: str, field: str, row: int) -> Decimal | None:
    value = raw.strip()
    if not value:
        return None
    if len(value) > 128:
        raise ValueError(f"Record {row}, {field}: numeric token is too long.")
    try:
        result = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"Record {row}, {field}: expected a number or an empty value.") from error
    if not result.is_finite():
        raise ValueError(f"Record {row}, {field}: expected a finite representable value.")
    numeric = float(result)
    if not math.isfinite(numeric) or (result != 0 and numeric == 0):
        raise ValueError(f"Record {row}, {field}: expected a finite representable value.")
    return result


def _integer(raw: str, field: str, row: int, *, maximum: int) -> int | None:
    result = _decimal(raw, field, row)
    if result is None:
        return None
    if result < 0 or result > maximum or result != result.to_integral_value():
        raise ValueError(f"Record {row}, {field}: expected a nonnegative integer in range.")
    return int(result)


def _faults(raw: str, row: int) -> int | None:
    value = raw.strip()
    if not value:
        return None
    if len(value) > 128:
        raise ValueError(f"Record {row}, fault_flags: mask token is too long.")
    if value.casefold().startswith("0x"):
        try:
            flags = int(value, 16)
        except ValueError as error:
            raise ValueError(f"Record {row}, fault_flags: invalid hexadecimal mask.") from error
        if not 0 <= flags <= (1 << 32) - 1:
            raise ValueError(f"Record {row}, fault_flags: mask exceeds the MCU uint32 range.")
        return flags
    return _integer(value, "fault_flags", row, maximum=(1 << 32) - 1)


def _millimetres(value: Decimal | None, row: int) -> str:
    if value is None:
        return ""
    # Keep decimal input exactly during the only dimensional conversion.
    with localcontext() as context:
        context.prec = max(28, len(value.as_tuple().digits) + 4)
        converted = value * Decimal(1000)
    if not math.isfinite(float(converted)):
        raise ValueError(f"Record {row}, displacement_m: millimetre conversion overflow.")
    return str(converted)


def _value(value: Decimal | int | None) -> str:
    return "" if value is None else str(value)


def convert_csv(input_path: str | Path, output_path: str | Path, *, time_basis: str,
                strain_channel_is_gauge: bool = False, source: str = "mcu-telemetry") -> dict:
    """Validate and produce CAD CSV, exact source copy and lossless record sidecar.

    Empty cells remain empty and are flagged; no missing observation becomes a
    numerical zero. All originals, including unknown added sensor channels,
    remain in the source copy and sidecar. Output files must not already exist.
    """
    if time_basis != "physical_seconds":
        raise ValueError("Use time_basis='physical_seconds' for measured laboratory time. Nondimensional PDE/model time is not supported.")
    if not source.strip() or len(source) > 256 or "\n" in source or "\r" in source:
        raise ValueError("Provide a nonempty source label of at most 256 characters on one line.")
    original = Path(input_path).resolve()
    output = Path(output_path).resolve()
    source_copy = output.with_name(output.stem + ".source.csv")
    sidecar = output.with_name(output.stem + ".metadata.json")
    targets = (output, source_copy, sidecar)
    if output.suffix.casefold() != ".csv":
        raise ValueError("The CAD measurement output must have the .csv extension.")
    if original in targets or len(set(targets)) != 3:
        raise ValueError("Output, source copy and metadata must be distinct from the input file.")
    if any(path.exists() for path in targets):
        raise FileExistsError("A handoff target already exists; choose a new output name to preserve it.")
    if not output.parent.is_dir():
        raise FileNotFoundError("The output directory does not exist.")

    # Stage all files beside the output. Invalid data cannot publish half of a
    # handoff or replace an existing log. No user source is ever rewritten.
    with tempfile.TemporaryDirectory(prefix=".measurement-adapter-", dir=output.parent) as temporary:
        staging = Path(temporary).resolve()
        if staging.parent != output.parent:
            raise RuntimeError("Temporary staging escaped the intended output directory.")
        staged_source = staging / "source.csv"
        staged_csv = staging / "measurements.csv"
        staged_metadata = staging / "metadata.json"
        shutil.copyfile(original, staged_source)
        source_sha = _hash_file(staged_source)
        counters: Counter[str] = Counter()
        row_count = 0
        last_time: Decimal | None = None
        last_cycle: int | None = None
        with staged_source.open("r", encoding="utf-8-sig", newline="") as source_stream, \
                staged_csv.open("w", encoding="utf-8", newline="") as cad_stream, \
                staged_metadata.open("w", encoding="utf-8", newline="") as metadata_stream:
            reader = csv.DictReader(source_stream, restval="")
            columns = reader.fieldnames or []
            if not columns or any(not name for name in columns) or len(set(columns)) != len(columns):
                raise ValueError("Use a nonempty CSV header with unique, named columns.")
            missing = [name for name in REQUIRED_COLUMNS if name not in columns]
            if missing:
                raise ValueError(f"Missing MCU telemetry columns: {', '.join(missing)}")
            if strain_channel_is_gauge and "strain" not in columns:
                raise ValueError("The independent-gauge declaration requires a 'strain' column.")
            writer = csv.DictWriter(cad_stream, fieldnames=CAD_COLUMNS)
            writer.writeheader()
            metadata = {
                "format": FORMAT, "status": "measurement_handoff_not_solver_input",
                "time_basis": "physical_seconds", "model_time_mapping": None,
                "source_file_name": original.name, "source_copy_file_name": source_copy.name,
                "source_sha256": source_sha, "source_bytes": staged_source.stat().st_size,
                "cad_file_name": output.name, "cad_columns": list(CAD_COLUMNS),
                "original_columns": columns, "source_label": source,
                "strain_channel_is_independent_gauge": bool(strain_channel_is_gauge),
                "mapping": {"time_s": "time_s", "cycle": "cycle_index", "force_n": "force_N",
                            "displacement_m": "crosshead_displacement_mm = displacement_m * 1000",
                            "temperature_c": "temperature_C", "strain": "gauge_strain only when explicitly declared independent gauge"},
                "units": {"time_s": "s (laboratory)", "force_N": "N", "crosshead_displacement_mm": "mm",
                          "gauge_strain": "1 (independent gauge only)", "temperature_C": "degC"},
                "unmodified_optional_channels": [name for name in OPTIONAL_CHANNELS if name in columns],
                "missing_value_policy": "Keep empty source cells empty; report missing flags, never impute zero.",
                "notes": ["No probability solver, clock conversion, hardware command or firmware upload occurs.",
                          "No machine force, frequency, stroke, material calibration or safe operating rating is inferred.",
                          "Original records preserve every source channel; command and DCPD channels are not discarded.",
                          "Crosshead travel includes machine/grip compliance and is not specimen gauge strain."],
            }
            metadata_stream.write(json.dumps(metadata, ensure_ascii=False, indent=2)[:-1])
            metadata_stream.write(',\n  "original_records": [')
            for row_number, raw in enumerate(reader, 1):
                if None in raw:
                    raise ValueError(f"Record {row_number}: more values than named columns.")
                quality: list[str] = []
                time = _decimal(raw["time_s"], "time_s", row_number)
                cycle = _integer(raw["cycle"], "cycle", row_number, maximum=(1 << 64) - 1)
                force = _decimal(raw["force_n"], "force_n", row_number)
                displacement = _decimal(raw["displacement_m"], "displacement_m", row_number)
                temperature = _decimal(raw["temperature_c"], "temperature_c", row_number)
                for name, value in (("time_s", time), ("cycle", cycle), ("force_n", force),
                                    ("displacement_m", displacement), ("temperature_c", temperature)):
                    if value is None:
                        quality.append("missing_" + name)
                if time is not None:
                    if last_time is not None and time < last_time:
                        raise ValueError(f"Record {row_number}, time_s: timestamps decrease; split separate runs rather than merging their clocks.")
                    if last_time is not None and time == last_time:
                        quality.append("duplicate_time_s")
                    last_time = time
                if cycle is not None:
                    if last_cycle is not None and cycle < last_cycle:
                        raise ValueError(f"Record {row_number}, cycle: cycle count decreases; split separate runs.")
                    last_cycle = cycle
                if strain_channel_is_gauge:
                    strain = _decimal(raw["strain"], "strain", row_number)
                    if strain is None:
                        quality.append("missing_gauge_strain")
                else:
                    strain = None
                    quality.append("gauge_strain_not_declared")
                if "fault_flags" not in columns:
                    quality.append("fault_flags_unavailable")
                else:
                    flags = _faults(raw["fault_flags"], row_number)
                    if flags is None:
                        quality.append("missing_fault_flags")
                    elif flags:
                        quality.append(f"fault_flags=0x{flags:X}")
                for flag in quality:
                    counters["fault_flags_nonzero" if flag.startswith("fault_flags=") else flag] += 1
                writer.writerow(dict(zip(CAD_COLUMNS, (
                    _value(time), _value(cycle), _value(force), _millimetres(displacement, row_number),
                    _value(strain), _value(temperature), source, ";".join(quality) if quality else "ok",
                ))))
                record = {"record_index": row_number, "ending_source_line": reader.line_num,
                          "channels": raw, "quality_flags": quality}
                metadata_stream.write((",\n" if row_count else "\n") + "    " + json.dumps(record, ensure_ascii=False))
                row_count += 1
            metadata_stream.write('\n  ],\n  "summary": ')
            metadata_stream.write(json.dumps({"row_count": row_count, "quality_counts": dict(counters),
                                             "empty_template": row_count == 0}, ensure_ascii=False))
            # Finish the CSV before binding its exact bytes in the sidecar.
            cad_stream.flush()
            metadata_stream.write(',\n  "cad_sha256": ' + json.dumps(_hash_file(staged_csv)) + "\n}\n")
        published: list[Path] = []
        try:
            # Exclusive creation via same-directory hard links avoids replacing
            # even a file created by another process after initial validation.
            for staged, target in ((staged_source, source_copy), (staged_csv, output), (staged_metadata, sidecar)):
                os.link(staged, target)
                published.append(target)
        except OSError:
            for target in reversed(published):
                if target.resolve().parent != output.parent:
                    raise RuntimeError("Refusing cleanup outside the handoff directory.")
                target.unlink()
            raise
    return {"cad_csv": str(output), "source_csv": str(source_copy), "metadata_json": str(sidecar),
            "row_count": row_count, "quality_counts": dict(counters), "source_sha256": source_sha,
            "time_basis": "physical_seconds", "model_time_mapping": None}


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert physical MCU telemetry to a CAD measurement handoff; no solver or hardware control.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--time-basis", required=True, help="Must be physical_seconds; model/PDE time is refused.")
    parser.add_argument("--strain-channel-is-gauge", action="store_true", help="Declare that source strain is an independently measured specimen gauge channel.")
    parser.add_argument("--source", default="mcu-telemetry", help="Measurement provenance label, not a material or solver calibration.")
    args = parser.parse_args()
    try:
        result = convert_csv(args.input_csv, args.output_csv, time_basis=args.time_basis,
                             strain_channel_is_gauge=args.strain_channel_is_gauge, source=args.source)
    except (OSError, UnicodeError, ValueError, csv.Error) as error:
        parser.exit(2, f"Measurement handoff refused: {error}\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
