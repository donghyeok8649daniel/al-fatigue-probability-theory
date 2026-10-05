"""Experimental handoff must preserve measurements, units and provenance."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


ADAPTER_PATH = Path(__file__).resolve().parents[1] / "measurement_csv_adapter.py"
SPEC = importlib.util.spec_from_file_location("aft_measurement_csv_adapter", ADAPTER_PATH)
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)


def _source(tmp_path, rows, *, header=None, bom=False):
    path = tmp_path / "raw.csv"
    header = header or ["time_s", "cycle", "force_n", "displacement_m", "temperature_c",
                        "strain", "dcpd_v", "actuator_command", "fault_flags", "stress_ref_pa", "force_ref_n", "custom_sensor"]
    with path.open("w", encoding="utf-8-sig" if bom else "utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _row(**overrides):
    data = dict(time_s="0", cycle="0", force_n="-12.5", displacement_m="0.001234",
                temperature_c="24.5", strain="0.00004", dcpd_v="1.2e-6", actuator_command="-0.03",
                fault_flags="0", stress_ref_pa="-25000000", force_ref_n="-12.6", custom_sensor="probe, alpha")
    data.update(overrides)
    return data


def _converted(path):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        return reader.fieldnames, list(reader)


def _metadata(path):
    return json.loads(path.with_name(path.stem + ".metadata.json").read_text(encoding="utf-8"))


def _no_outputs(path):
    assert not path.exists()
    assert not path.with_name(path.stem + ".source.csv").exists()
    assert not path.with_name(path.stem + ".metadata.json").exists()


def test_units_signed_force_source_bytes_and_all_channels_survive(tmp_path):
    source = _source(tmp_path, [_row(time_s="-0.5"), _row(time_s="0.25", cycle="1", fault_flags="0x00000005")], bom=True)
    original = source.read_bytes()
    output = tmp_path / "cad.csv"
    result = adapter.convert_csv(source, output, time_basis="physical_seconds",
                                 strain_channel_is_gauge=True, source="실측 로드셀")
    columns, rows = _converted(output)
    assert columns == list(adapter.CAD_COLUMNS)
    assert rows[0]["force_N"] == "-12.5"
    assert rows[0]["crosshead_displacement_mm"] == "1.234000"
    assert rows[0]["time_s"] == "-0.5"  # A physical pretrigger time is retained.
    assert rows[0]["gauge_strain"] == "0.00004"
    assert rows[0]["source"] == "실측 로드셀" and rows[0]["quality_flag"] == "ok"
    assert rows[1]["quality_flag"] == "fault_flags=0x5"
    assert result["quality_counts"] == {"fault_flags_nonzero": 1}
    assert Path(result["source_csv"]).read_bytes() == original == source.read_bytes()
    metadata = _metadata(output)
    assert metadata["source_sha256"] == hashlib.sha256(original).hexdigest()
    assert metadata["cad_sha256"] == hashlib.sha256(output.read_bytes()).hexdigest()
    raw = metadata["original_records"][1]["channels"]
    assert raw["fault_flags"] == "0x00000005"
    assert raw["dcpd_v"] == "1.2e-6" and raw["actuator_command"] == "-0.03"
    assert raw["stress_ref_pa"] == "-25000000" and raw["force_ref_n"] == "-12.6"
    assert raw["custom_sensor"] == "probe, alpha"
    assert metadata["time_basis"] == "physical_seconds" and metadata["model_time_mapping"] is None


def test_crosshead_or_undeclared_strain_is_never_relabelled_as_gauge(tmp_path):
    source = _source(tmp_path, [_row(strain="0.123456")])
    output = tmp_path / "cad.csv"
    adapter.convert_csv(source, output, time_basis="physical_seconds")
    _, rows = _converted(output)
    assert rows[0]["gauge_strain"] == ""
    assert rows[0]["quality_flag"] == "gauge_strain_not_declared"
    metadata = _metadata(output)
    assert metadata["original_records"][0]["channels"]["strain"] == "0.123456"
    assert metadata["strain_channel_is_independent_gauge"] is False


def test_missing_measurements_and_fault_state_are_not_imputed_as_zero(tmp_path):
    source = _source(tmp_path, [_row(time_s="", cycle="", force_n="", displacement_m="", temperature_c="", strain="", fault_flags="")])
    output = tmp_path / "cad.csv"
    adapter.convert_csv(source, output, time_basis="physical_seconds", strain_channel_is_gauge=True)
    _, rows = _converted(output)
    for name in ("time_s", "cycle_index", "force_N", "crosshead_displacement_mm", "gauge_strain", "temperature_C"):
        assert rows[0][name] == ""
    quality = rows[0]["quality_flag"].split(";")
    assert "missing_force_n" in quality and "missing_fault_flags" in quality and "missing_gauge_strain" in quality
    assert _metadata(output)["original_records"][0]["channels"]["force_n"] == ""


@pytest.mark.parametrize("basis", ["model_time", "nondimensional", "s", "", "physical_seconds_assumed"])
def test_pde_or_undeclared_time_basis_is_refused_before_writing(tmp_path, basis):
    source = _source(tmp_path, [_row()])
    output = tmp_path / "cad.csv"
    with pytest.raises(ValueError, match="Nondimensional PDE/model time"):
        adapter.convert_csv(source, output, time_basis=basis)
    _no_outputs(output)


@pytest.mark.parametrize("field,value", [("force_n", "NaN"), ("time_s", "Infinity"),
    ("temperature_c", "-Inf"), ("displacement_m", "1e309"), ("force_n", "1e-1000"),
    ("cycle", "0.5"), ("cycle", "-1"), ("fault_flags", "invalid"), ("fault_flags", "0x100000000")])
def test_bad_numeric_data_cannot_publish_a_partial_handoff(tmp_path, field, value):
    source = _source(tmp_path, [_row(), _row(time_s="1", **{field: value})] if field != "time_s" else [_row(), _row(time_s=value)])
    original = source.read_bytes()
    output = tmp_path / "cad.csv"
    with pytest.raises(ValueError):
        adapter.convert_csv(source, output, time_basis="physical_seconds")
    _no_outputs(output)
    assert source.read_bytes() == original


@pytest.mark.parametrize("row", [_row(time_s="-1"), _row(time_s="1", cycle="0")])
def test_separate_run_clocks_or_cycles_are_not_silently_stitched(tmp_path, row):
    source = _source(tmp_path, [_row(time_s="0", cycle="1"), row])
    output = tmp_path / "cad.csv"
    with pytest.raises(ValueError, match="decrease"):
        adapter.convert_csv(source, output, time_basis="physical_seconds")
    _no_outputs(output)


def test_duplicate_time_is_preserved_and_flagged(tmp_path):
    source = _source(tmp_path, [_row(), _row(force_n="13.0")])
    output = tmp_path / "cad.csv"
    adapter.convert_csv(source, output, time_basis="physical_seconds", strain_channel_is_gauge=True)
    _, rows = _converted(output)
    assert len(rows) == 2 and rows[1]["force_N"] == "13.0"
    assert rows[1]["quality_flag"] == "duplicate_time_s"


def test_precision_is_preserved_in_the_metre_to_mm_conversion(tmp_path):
    source = _source(tmp_path, [_row(displacement_m="1.234567891234567891234567891234567891e-6")])
    output = tmp_path / "cad.csv"
    adapter.convert_csv(source, output, time_basis="physical_seconds")
    _, rows = _converted(output)
    assert rows[0]["crosshead_displacement_mm"] == "0.001234567891234567891234567891234567891000"


def test_unknown_fault_channel_and_empty_template_are_explicit(tmp_path):
    header = ["time_s", "cycle", "force_n", "displacement_m", "temperature_c"]
    source = _source(tmp_path, [dict(zip(header, ("0", "0", "1", "0.001", "20")))], header=header)
    output = tmp_path / "cad.csv"
    adapter.convert_csv(source, output, time_basis="physical_seconds")
    assert "fault_flags_unavailable" in _converted(output)[1][0]["quality_flag"]
    empty = _source(tmp_path, [], header=header)
    template = tmp_path / "template.csv"
    result = adapter.convert_csv(empty, template, time_basis="physical_seconds")
    assert result["row_count"] == 0
    assert _metadata(template)["summary"]["empty_template"] is True


def test_existing_output_and_source_copy_collisions_do_not_replace_user_data(tmp_path):
    source = _source(tmp_path, [_row()])
    output = tmp_path / "cad.csv"
    output.write_bytes(b"existing user's measurements")
    with pytest.raises(FileExistsError):
        adapter.convert_csv(source, output, time_basis="physical_seconds")
    assert output.read_bytes() == b"existing user's measurements"
    assert not (tmp_path / "cad.source.csv").exists()
    with pytest.raises(ValueError, match="distinct"):
        adapter.convert_csv(source, source, time_basis="physical_seconds")
    dangerous_output = tmp_path / "raw.csv"
    dangerous_input = tmp_path / "raw.source.csv"
    dangerous_input.write_bytes(source.read_bytes())
    with pytest.raises(ValueError, match="distinct"):
        adapter.convert_csv(dangerous_input, dangerous_output, time_basis="physical_seconds")


@pytest.mark.parametrize("header", ["time_s,cycle,force_n,force_n,displacement_m,temperature_c", "time_s,cycle,displacement_m,temperature_c"])
def test_ambiguous_or_missing_named_measurement_columns_are_refused(tmp_path, header):
    source = tmp_path / "raw.csv"
    source.write_text(header + "\n", encoding="utf-8")
    output = tmp_path / "cad.csv"
    with pytest.raises(ValueError):
        adapter.convert_csv(source, output, time_basis="physical_seconds")
    _no_outputs(output)


def test_extra_csv_values_cannot_disappear_from_named_channel_mapping(tmp_path):
    source = tmp_path / "raw.csv"
    source.write_text("time_s,cycle,force_n,displacement_m,temperature_c\n0,0,1,0.001,20,unexpected\n", encoding="utf-8")
    output = tmp_path / "cad.csv"
    with pytest.raises(ValueError, match="more values"):
        adapter.convert_csv(source, output, time_basis="physical_seconds")
    _no_outputs(output)


def test_publish_collision_rolls_back_only_our_new_files(tmp_path, monkeypatch):
    source = _source(tmp_path, [_row()])
    output = tmp_path / "cad.csv"
    real_link = adapter.os.link
    def collision(staged, target):
        if Path(target) == output:
            output.write_bytes(b"concurrent user's log")
        return real_link(staged, target)
    monkeypatch.setattr(adapter.os, "link", collision)
    with pytest.raises(FileExistsError):
        adapter.convert_csv(source, output, time_basis="physical_seconds")
    assert output.read_bytes() == b"concurrent user's log"
    assert not (tmp_path / "cad.source.csv").exists()
    assert not (tmp_path / "cad.metadata.json").exists()


def test_cli_requires_the_physical_time_declaration_and_runs_without_solver(tmp_path):
    source = _source(tmp_path, [_row()])
    output = tmp_path / "cad.csv"
    missing = subprocess.run([sys.executable, str(ADAPTER_PATH), str(source), str(output)], capture_output=True, text=True)
    assert missing.returncode == 2 and "--time-basis" in missing.stderr
    result = subprocess.run([sys.executable, str(ADAPTER_PATH), str(source), str(output), "--time-basis", "physical_seconds", "--strain-channel-is-gauge"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["row_count"] == 1
    assert _converted(output)[1][0]["quality_flag"] == "ok"
