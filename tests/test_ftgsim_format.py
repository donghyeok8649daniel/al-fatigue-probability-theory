import json
from pathlib import Path
import zipfile

import pytest

from simulations.fem_tension_app import (
    TensionRunConfig,
    config_from_ftgsim,
    save_tension_ftgsim,
)
from simulations.ftgsim_format import (MANIFEST_NAME, create_ftgsim, extract_geometry,
                                       extract_results, open_ftgsim)


def test_tension_project_round_trip_and_result_extraction(tmp_path: Path):
    result_dir = tmp_path / "run"
    result_dir.mkdir()
    (result_dir / "nodes.csv").write_text("step,node\n0,0\n", encoding="utf-8")
    (result_dir / "elements.csv").write_text("step,element\n0,0\n", encoding="utf-8")
    config = TensionRunConfig(elements=17, frequency_hz=8.0, deformation_scale=3.0)
    geometry_source = tmp_path / "part.obj"
    geometry_source.write_text("v 0 0 0\nv 1 0 0\nl 1 2\n", encoding="utf-8")
    project = save_tension_ftgsim(tmp_path / "sample", config, result_dir, view="3D",
                                  geometry_source=geometry_source)
    assert project.suffix == ".ftgsim"

    loaded, geometry, display = config_from_ftgsim(project)
    assert loaded == config
    assert geometry["mesh_dimension"] == 1
    assert geometry["loading_axis"] == [1.0, 0.0, 0.0]
    assert display["view"] == "3D"
    assert geometry["source_member"] == "geometry/source.obj"
    bundle = open_ftgsim(project)
    extracted = extract_results(bundle, tmp_path / "opened")
    assert {item.name for item in extracted} == {"nodes.csv", "elements.csv"}
    assert (tmp_path / "opened" / "nodes.csv").read_text(encoding="utf-8") == "step,node\n0,0\n"
    geometry_files = extract_geometry(bundle, tmp_path / "opened_geometry")
    assert len(geometry_files) == 1
    assert geometry_files[0].read_text(encoding="utf-8").startswith("v 0 0 0")
    (tmp_path / "opened" / "nodes.csv").write_text("user data\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        extract_results(bundle, tmp_path / "opened")


def test_checksum_tampering_is_rejected(tmp_path: Path):
    path = create_ftgsim(tmp_path / "clean.ftgsim", setup={}, geometry={}, display={})
    rewritten = tmp_path / "tampered.ftgsim"
    with zipfile.ZipFile(path) as source, zipfile.ZipFile(rewritten, "w") as target:
        for info in source.infolist():
            data = source.read(info.filename)
            if info.filename == "setup.json":
                data = b"{}\n "
            target.writestr(info, data)
    with pytest.raises(ValueError, match="checksum mismatch"):
        open_ftgsim(rewritten)


def test_path_traversal_and_executable_members_are_rejected(tmp_path: Path):
    for member in ("../escape.csv", "payload.exe"):
        path = tmp_path / (member.replace(".", "x").replace("/", "_") + ".ftgsim")
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(MANIFEST_NAME, json.dumps({"format": "ftgsim"}))
            archive.writestr(member, b"bad")
        with pytest.raises(ValueError):
            open_ftgsim(path)


def test_wrong_extension_and_physics_model_are_rejected(tmp_path: Path):
    wrong = create_ftgsim(tmp_path / "other.ftgsim", setup={"physics_model": "multiaxial"},
                          geometry={}, display={})
    with pytest.raises(ValueError, match="not a 1D normal-tensile"):
        config_from_ftgsim(wrong)
    renamed = tmp_path / "other.zip"
    renamed.write_bytes(wrong.read_bytes())
    with pytest.raises(ValueError, match="expected a .ftgsim"):
        open_ftgsim(renamed)
