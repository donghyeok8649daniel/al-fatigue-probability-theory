"""Portable, inspectable and non-executable ``.ftgsim`` bundle format.

The container is ZIP/ZIP64 with UTF-8 JSON/CSV/PNG members.  Loading never
executes bundle content and extraction is limited to explicitly allowed data
members after path and size validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import tempfile
from typing import Mapping
import uuid
import zipfile


EXTENSION = ".ftgsim"
FORMAT_NAME = "ftgsim"
SCHEMA_VERSION = "1.0.0"
MANIFEST_NAME = "ftgsim-manifest.json"
MAX_FILES = 256
MAX_MEMBER_BYTES = 128 * 1024 * 1024
MAX_TOTAL_BYTES = 512 * 1024 * 1024
ALLOWED_SUFFIXES = {".json", ".csv", ".png", ".svg", ".txt", ".md",
                    ".obj", ".stl", ".ply", ".vtk"}


@dataclass(frozen=True)
class FTGSimBundle:
    path: Path
    manifest: dict
    setup: dict
    geometry: dict
    display: dict
    members: tuple[str, ...]


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n").encode("utf-8")


def _safe_member(name: str) -> PurePosixPath:
    if not name or "\\" in name:
        raise ValueError(f"invalid archive member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or ":" in path.parts[0]:
        raise ValueError(f"unsafe archive member path: {name!r}")
    if path.suffix.lower() not in ALLOWED_SUFFIXES:
        raise ValueError(f"unsupported bundle member type: {name!r}")
    return path


def _validate_manifest(manifest: dict) -> None:
    if manifest.get("format") != FORMAT_NAME:
        raise ValueError("not an ftgsim bundle")
    version = str(manifest.get("schema_version", ""))
    try:
        major = int(version.split(".")[0])
    except (ValueError, IndexError):
        raise ValueError("invalid ftgsim schema_version") from None
    if major != 1:
        raise ValueError(f"unsupported ftgsim schema major version: {major}")
    if manifest.get("bundle_kind") not in {"project", "result"}:
        raise ValueError("bundle_kind must be project or result")
    paths = manifest.get("paths")
    if not isinstance(paths, dict):
        raise ValueError("manifest paths must be an object")
    for key in ("setup", "geometry", "display"):
        if key not in paths:
            raise ValueError(f"manifest missing paths.{key}")
        _safe_member(str(paths[key]))


def create_ftgsim(
    path: Path,
    *,
    setup: Mapping,
    geometry: Mapping,
    display: Mapping | None = None,
    files: Mapping[str, Path | bytes] | None = None,
    bundle_kind: str = "project",
    generator: Mapping | None = None,
) -> Path:
    """Atomically write a versioned `.ftgsim` bundle."""
    target = Path(path)
    if target.suffix.lower() != EXTENSION:
        target = target.with_suffix(EXTENSION)
    target.parent.mkdir(parents=True, exist_ok=True)
    payloads: dict[str, bytes] = {
        "setup.json": _json_bytes(dict(setup)),
        "geometry.json": _json_bytes(dict(geometry)),
        "display.json": _json_bytes(dict(display or {})),
    }
    for name, source in (files or {}).items():
        safe = _safe_member(name).as_posix()
        if safe in payloads or safe == MANIFEST_NAME:
            raise ValueError(f"duplicate or reserved bundle member: {safe}")
        data = bytes(source) if isinstance(source, bytes) else Path(source).read_bytes()
        if len(data) > MAX_MEMBER_BYTES:
            raise ValueError(f"bundle member exceeds size limit: {safe}")
        payloads[safe] = data
    if len(payloads) + 1 > MAX_FILES or sum(map(len, payloads.values())) > MAX_TOTAL_BYTES:
        raise ValueError("bundle exceeds file-count or total-size limit")
    checksums = {name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()}
    manifest = {
        "format": FORMAT_NAME,
        "schema_version": SCHEMA_VERSION,
        "bundle_kind": bundle_kind,
        "bundle_id": str(uuid.uuid4()),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "generator": dict(generator or {}),
        "scope": {
            "material": "pure single-crystal aluminum",
            "theory": "1D normal tensile crack-initiation model",
            "shear_or_multiaxial": False,
        },
        "paths": {"setup": "setup.json", "geometry": "geometry.json", "display": "display.json"},
        "checksums_sha256": checksums,
    }
    _validate_manifest(manifest)
    handle, temporary_name = tempfile.mkstemp(prefix=target.stem + "-", suffix=".tmp", dir=target.parent)
    os.close(handle)
    temporary = Path(temporary_name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
            archive.writestr(MANIFEST_NAME, _json_bytes(manifest))
            for name, data in payloads.items():
                archive.writestr(name, data)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()
    return target


def open_ftgsim(path: Path) -> FTGSimBundle:
    """Validate metadata, paths, sizes and SHA-256 checksums before reading."""
    source = Path(path)
    if source.suffix.lower() != EXTENSION:
        raise ValueError(f"expected a {EXTENSION} file")
    with zipfile.ZipFile(source, "r") as archive:
        infos = archive.infolist()
        if len(infos) > MAX_FILES:
            raise ValueError("bundle contains too many members")
        names: list[str] = []
        total = 0
        for info in infos:
            name = _safe_member(info.filename).as_posix()
            if name in names:
                raise ValueError(f"duplicate archive member: {name}")
            # Unix file-type bits in external_attr: reject symlinks.
            if ((info.external_attr >> 16) & 0o170000) == 0o120000:
                raise ValueError(f"symbolic links are forbidden: {name}")
            if info.file_size > MAX_MEMBER_BYTES:
                raise ValueError(f"bundle member exceeds size limit: {name}")
            total += info.file_size
            names.append(name)
        if total > MAX_TOTAL_BYTES or MANIFEST_NAME not in names:
            raise ValueError("invalid bundle size or missing manifest")
        manifest = json.loads(archive.read(MANIFEST_NAME))
        _validate_manifest(manifest)
        checksums = manifest.get("checksums_sha256", {})
        for name, expected in checksums.items():
            safe = _safe_member(name).as_posix()
            if safe not in names:
                raise ValueError(f"manifest references missing member: {safe}")
            actual = hashlib.sha256(archive.read(safe)).hexdigest()
            if actual != expected:
                raise ValueError(f"checksum mismatch: {safe}")
        paths = manifest["paths"]
        setup = json.loads(archive.read(paths["setup"]))
        geometry = json.loads(archive.read(paths["geometry"]))
        display = json.loads(archive.read(paths["display"]))
    return FTGSimBundle(source, manifest, setup, geometry, display, tuple(names))


def extract_results(bundle: FTGSimBundle, destination: Path) -> tuple[Path, ...]:
    """Extract only checksummed `results/*.csv|json` members safely."""
    target = Path(destination).resolve()
    target.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    checksums = bundle.manifest.get("checksums_sha256", {})
    with zipfile.ZipFile(bundle.path, "r") as archive:
        for name in bundle.members:
            pure = _safe_member(name)
            if not pure.parts or pure.parts[0] != "results" or pure.suffix.lower() not in {".csv", ".json"}:
                continue
            if name not in checksums:
                raise ValueError(f"result member is not checksummed: {name}")
            relative = Path(*pure.parts[1:])
            output = (target / relative).resolve()
            if target != output and target not in output.parents:
                raise ValueError(f"unsafe extraction target: {name}")
            output.parent.mkdir(parents=True, exist_ok=True)
            data = archive.read(name)
            if output.exists() and output.read_bytes() != data:
                raise FileExistsError(f"refusing to overwrite an existing result: {output}")
            if not output.exists():
                output.write_bytes(data)
            extracted.append(output)
    return tuple(extracted)


def extract_geometry(bundle: FTGSimBundle, destination: Path) -> tuple[Path, ...]:
    """Extract only a checksummed, non-executable mesh below `geometry/`."""
    target = Path(destination).resolve()
    target.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    checksums = bundle.manifest.get("checksums_sha256", {})
    with zipfile.ZipFile(bundle.path, "r") as archive:
        for name in bundle.members:
            pure = _safe_member(name)
            if not pure.parts or pure.parts[0] != "geometry" or pure.suffix.lower() not in {".obj", ".stl", ".ply", ".vtk"}:
                continue
            if name not in checksums:
                raise ValueError(f"geometry member is not checksummed: {name}")
            output = (target / Path(*pure.parts[1:])).resolve()
            if target != output and target not in output.parents:
                raise ValueError(f"unsafe extraction target: {name}")
            data = archive.read(name)
            if output.exists() and output.read_bytes() != data:
                raise FileExistsError(f"refusing to overwrite existing geometry: {output}")
            if not output.exists():
                output.write_bytes(data)
            extracted.append(output)
    return tuple(extracted)
