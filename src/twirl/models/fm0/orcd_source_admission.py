"""Fail-closed source admission for ORCD-resident FM light curves.

This module verifies that a later sector is ready for FM-specific HDF5 and
six-view validation.  It does not promote the sector to accepted A2v1, freeze
a temporal panel, or claim that the light curves themselves have been opened.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .registry import FM0ContractError

SOURCE_RECEIPT_SCHEMA_VERSION = "twirl_fm0_orcd_source_receipt_v1"
SOURCE_READY_STATE = "FM_ORCD_SOURCE_READY"
LC_VERIFIED_STATE = "FM_ORCD_LC_VERIFIED"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_CELL = re.compile(r"^s(?P<sector>[0-9]{4})_o(?P<orbit>[0-9]+)_cam(?P<camera>[1-4])_ccd(?P<ccd>[1-4])$")


def _read_json(path: Path, *, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} must be a JSON object: {path}")
    return value


def _digest(value: Any, *, label: str) -> str:
    result = str(value).strip().lower()
    if _SHA256.fullmatch(result) is None:
        raise FM0ContractError(f"{label} is not a lowercase SHA-256 digest")
    return result


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _expected_cells(sector: int, orbits: Sequence[int]) -> set[str]:
    if len(orbits) != 2 or len(set(orbits)) != 2:
        raise FM0ContractError("exactly two distinct expected orbits are required")
    return {
        f"s{sector:04d}_o{orbit}_cam{camera}_ccd{ccd}"
        for orbit in orbits
        for camera in range(1, 5)
        for ccd in range(1, 5)
    }


def _retention_cell_label(value: Any) -> str:
    """Normalize the two historical retention-receipt cell encodings."""

    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        label = str(value.get("label", "")).strip()
        match = _CELL.fullmatch(label)
        if match is None:
            return ""
        try:
            identity = (
                int(value.get("sector")),
                int(value.get("orbit")),
                int(value.get("camera")),
                int(value.get("ccd")),
            )
        except (TypeError, ValueError):
            return ""
        parsed = (
            int(match.group("sector")),
            int(match.group("orbit")),
            int(match.group("camera")),
            int(match.group("ccd")),
        )
        return label if identity == parsed else ""
    return ""


def verify_archive_source(
    *, archive_dir: str | Path, sector: int, expected_orbits: Sequence[int]
) -> dict[str, Any]:
    """Verify one already-ingressed checksum-bound sector tar authority."""

    root = Path(archive_dir).expanduser().resolve()
    tag = f"s{sector:04d}"
    archive = root / f"{tag}_A2v1_raw_hdf5.tar"
    sidecar = root / f"{archive.name}.sha256"
    summary_path = root / "summary.json"
    ready = root / "READY"
    ready_orcd = root / "READY_ORCD"
    for path in (archive, sidecar, summary_path, ready, ready_orcd):
        if not path.is_file() or path.stat().st_size <= 0:
            raise FM0ContractError(f"missing ORCD archive authority: {path}")

    summary = _read_json(summary_path, label="archive summary")
    if (
        summary.get("schema_version") != "twirl_fm0_1_a2v1_sector_tar_v1"
        or int(summary.get("sector", -1)) != sector
        or tuple(int(v) for v in summary.get("orbits", ())) != tuple(expected_orbits)
        or summary.get("archive") != archive.name
        or summary.get("accepted_sources_modified") is not False
    ):
        raise FM0ContractError("ORCD archive summary identity drifted")
    if int(summary.get("n_hdf5_files", 0)) <= 0:
        raise FM0ContractError("ORCD archive contains no declared HDF5 products")
    if int(summary.get("archive_bytes", -1)) != archive.stat().st_size:
        raise FM0ContractError("ORCD archive byte count drifted")

    expected = _digest(summary.get("archive_sha256"), label="archive SHA-256")
    fields = sidecar.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2 or fields[1] != archive.name or fields[0] != expected:
        raise FM0ContractError("ORCD archive checksum sidecar drifted")
    if ready.read_text(encoding="utf-8").strip() != expected:
        raise FM0ContractError("ORCD archive READY digest drifted")
    ingress = _read_json(ready_orcd, label="ORCD ingress receipt")
    if ingress.get("tag") != tag or not str(ingress.get("task_id", "")).strip():
        raise FM0ContractError("ORCD ingress receipt identity drifted")

    return {
        "schema_version": SOURCE_RECEIPT_SCHEMA_VERSION,
        "sector": sector,
        "expected_orbits": list(expected_orbits),
        "source_form": "checksum_bound_sector_tar",
        "source_state": SOURCE_READY_STATE,
        "product_state": "A2V1_ACCEPTED",
        "pdo_sector_accepted": True,
        "n_cells": 32,
        "n_hdf5_products_declared": int(summary["n_hdf5_files"]),
        "source_root": str(root),
        "source_bindings": {
            "archive_path": str(archive),
            "archive_sha256": expected,
            "summary_sha256": _sha256(summary_path),
            "ingress_receipt_sha256": _sha256(ready_orcd),
        },
        "hdf5_openability_verified": False,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "claim_limit": "ORCD FM source readiness only; HDF5/shard gates remain required",
    }


def verify_retained_sector_source(
    *, sector_root: str | Path, sector: int, expected_orbits: Sequence[int]
) -> dict[str, Any]:
    """Verify all 32 retained ORCD cell receipts for one later sector."""

    root = Path(sector_root).expanduser().resolve()
    outputs = root / "outputs"
    if not outputs.is_dir():
        raise FM0ContractError(f"missing retained sector outputs: {outputs}")
    expected_cells = _expected_cells(sector, expected_orbits)
    observed = {path.name for path in outputs.iterdir() if path.is_dir()}
    if observed != expected_cells:
        missing = sorted(expected_cells - observed)
        extra = sorted(observed - expected_cells)
        raise FM0ContractError(
            f"S{sector} retained cells differ; missing={missing}, extra={extra}"
        )

    bindings: list[dict[str, Any]] = []
    total_hdf5 = 0
    for cell in sorted(expected_cells):
        match = _CELL.fullmatch(cell)
        assert match is not None
        cell_root = outputs / cell
        complete_paths = sorted((cell_root / "complete").glob("*.json"))
        modern_retention = sorted(
            (cell_root / "retained").glob("*.retention.json")
        )
        retention_paths = modern_retention or sorted(
            (cell_root / "retained").glob("*/retention.json")
        )
        if len(complete_paths) != 1 or len(retention_paths) != 1:
            raise FM0ContractError(f"{cell} lacks one completion/retention receipt")
        complete_path, retention_path = complete_paths[0], retention_paths[0]
        complete = _read_json(complete_path, label=f"{cell} completion receipt")
        retained = _read_json(retention_path, label=f"{cell} retention receipt")
        if (
            complete.get("schema") != "twirl-a2v1-orcd-cell-complete-v1"
            or complete.get("state") != "complete"
            or complete.get("cell") != cell
            or retained.get("schema") != "twirl-a2v1-orcd-retention-v1"
            or retained.get("ok") is not True
            or _retention_cell_label(retained.get("cell")) != cell
            or retained.get("attempt_id") != complete.get("attempt_id")
            or retained.get("pdo_return_deferred") is not True
            or retained.get("pdo_sector_accepted") is not False
        ):
            raise FM0ContractError(f"{cell} completion/retention identity drifted")
        validation = retained.get("retained_root_validation")
        if not isinstance(validation, Mapping) or validation.get("ok") is not True:
            raise FM0ContractError(f"{cell} retained-root validation did not pass")
        checks = validation.get("checks")
        if not isinstance(checks, Mapping) or not checks or not all(
            value is True for value in checks.values()
        ):
            raise FM0ContractError(f"{cell} retained-root checks did not all pass")
        validation_hashes = validation.get("hashes", validation)
        if not isinstance(validation_hashes, Mapping):
            raise FM0ContractError(f"{cell} retained-root hashes are invalid")
        for key in (
            "input_manifest_sha256",
            "output_manifest_sha256",
            "environment_manifest_sha256",
        ):
            if _digest(complete.get(key), label=f"{cell} {key}") != _digest(
                validation_hashes.get(key), label=f"{cell} retained {key}"
            ):
                raise FM0ContractError(f"{cell} {key} differs after retention")
        if _digest(
            validation_hashes.get("completion_json_sha256"),
            label=f"{cell} completion receipt digest",
        ) != _sha256(complete_path):
            raise FM0ContractError(f"{cell} completion receipt digest drifted")
        hdf5 = validation.get("outputs", {}).get("hdf5", {})
        if not isinstance(hdf5, Mapping):
            raise FM0ContractError(f"{cell} lacks retained HDF5 inventory")
        validated = int(hdf5.get("validated_tics", 0))
        if validated <= 0 or list(hdf5.get("missing_expected", ())) or list(
            hdf5.get("extra_unrequested", ())
        ):
            raise FM0ContractError(f"{cell} retained HDF5 authority did not pass")
        retained_root = Path(str(retained.get("retained_root", ""))).resolve()
        if cell_root not in retained_root.parents or not retained_root.is_dir():
            raise FM0ContractError(f"{cell} retained root is missing or escapes cell")
        total_hdf5 += validated
        bindings.append(
            {
                "cell": cell,
                "completion_receipt": str(complete_path),
                "completion_receipt_sha256": _sha256(complete_path),
                "retention_receipt": str(retention_path),
                "retention_receipt_sha256": _sha256(retention_path),
                "retained_root": str(retained_root),
                "n_hdf5_products_declared": validated,
            }
        )

    return {
        "schema_version": SOURCE_RECEIPT_SCHEMA_VERSION,
        "sector": sector,
        "expected_orbits": list(expected_orbits),
        "source_form": "retained_orcd_cells",
        "source_state": SOURCE_READY_STATE,
        "product_state": "ORCD_COMPLETE_DEFERRED",
        "pdo_sector_accepted": False,
        "n_cells": len(bindings),
        "n_hdf5_products_declared": total_hdf5,
        "source_root": str(root),
        "cell_bindings": bindings,
        "hdf5_openability_verified": False,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "claim_limit": "ORCD FM source readiness only; HDF5/shard gates remain required",
    }


def write_source_receipt(receipt: Mapping[str, Any], output: str | Path) -> Path:
    """Write one immutable source receipt without replacing prior evidence."""

    target = Path(output).expanduser().resolve()
    if target.exists():
        raise FM0ContractError(f"refusing to overwrite source receipt: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return target
