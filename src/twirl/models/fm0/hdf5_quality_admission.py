"""Read-only HDF5 openability and cadence-quality admission on ORCD."""

from __future__ import annotations

import json
import multiprocessing
import re
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .mission_quality_admission import (
    MISSION_QUALITY_READY_STATE,
    MISSION_QUALITY_RECEIPT_SCHEMA_VERSION,
    _read_flag_file,
    _read_quat,
    verify_mission_quality_sources,
)
from .mission_quality_transfer import MISSION_QUALITY_TRANSFER_SCHEMA_VERSION
from .orcd_source_admission import SOURCE_READY_STATE, SOURCE_RECEIPT_SCHEMA_VERSION
from .registry import FM0ContractError, sha256_file

HDF5_QUALITY_RECEIPT_SCHEMA_VERSION = "twirl_fm0_hdf5_quality_admission_v1"
HDF5_QUALITY_READY_STATE = "FM_HDF5_QUALITY_READY"
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_CELL = re.compile(
    r"^s(?P<sector>[0-9]{4})_o(?P<orbit>[0-9]+)_cam(?P<camera>[1-4])_ccd(?P<ccd>[1-4])$"
)


@dataclass(frozen=True)
class Hdf5QualityAuthority:
    sector: int
    quaternion: Mapping[tuple[int, int], frozenset[int]]
    qflag: Mapping[tuple[int, int, int], Mapping[int, int]]
    mission: Mapping[tuple[int, int], Mapping[int, int]]
    exclusions: Mapping[tuple[int, int], frozenset[int]]


def _read_json(path: Path, *, label: str) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not isinstance(payload, Mapping):
        raise FM0ContractError(f"{label} must be a JSON object: {path}")
    return payload


def _verify_transfer_sector(
    *, transfer_root: Path, sector: int, expected_orbits: Sequence[int]
) -> tuple[Hdf5QualityAuthority, Mapping[str, Any], Path, Mapping[str, Any]]:
    manifest_path = transfer_root / "manifest.json"
    manifest = _read_json(manifest_path, label="mission-quality transfer manifest")
    if (
        manifest.get("schema_version") != MISSION_QUALITY_TRANSFER_SCHEMA_VERSION
        or sector not in [int(value) for value in manifest.get("sectors", ())]
    ):
        raise FM0ContractError(f"S{sector} is absent from the quality transfer")
    rows = [
        row
        for row in list(manifest.get("receipts", ()))
        + list(manifest.get("sources", ()))
        if isinstance(row, Mapping) and int(row.get("sector", -1)) == sector
    ]
    if not rows:
        raise FM0ContractError(f"S{sector} quality transfer inventory is empty")
    for row in rows:
        staged = transfer_root / str(row.get("staged_path", ""))
        if not staged.is_file() or sha256_file(staged) != str(row.get("sha256", "")):
            raise FM0ContractError(f"S{sector} transferred quality source drifted: {staged}")

    receipt_rows = [
        row
        for row in manifest.get("receipts", ())
        if isinstance(row, Mapping) and int(row.get("sector", -1)) == sector
    ]
    if len(receipt_rows) != 1:
        raise FM0ContractError(f"S{sector} transfer lacks one quality receipt")
    quality_receipt_path = transfer_root / str(receipt_rows[0]["staged_path"])
    quality_receipt = _read_json(quality_receipt_path, label="mission-quality receipt")
    if (
        quality_receipt.get("schema_version")
        != MISSION_QUALITY_RECEIPT_SCHEMA_VERSION
        or quality_receipt.get("quality_state") != MISSION_QUALITY_READY_STATE
        or quality_receipt.get("passed") is not True
        or int(quality_receipt.get("sector", -1)) != sector
    ):
        raise FM0ContractError(f"S{sector} transferred quality receipt is invalid")

    qlp_root = transfer_root / f"sources/s{sector:04d}/qlp"
    mission_root = transfer_root / f"sources/s{sector:04d}/mission"
    relocated_receipt = verify_mission_quality_sources(
        sector=sector,
        expected_orbits=expected_orbits,
        qlp_root=qlp_root,
        mission_flag_root=mission_root,
    )
    for key in (
        "quality_policy",
        "mission_quality_type",
        "n_qflag_rows",
        "n_mission_quality_rows",
        "n_mission_quality_rows_retained",
        "n_mission_quality_rows_excluded_by_quat",
        "mission_quality_authority_exclusions_sha256",
    ):
        if relocated_receipt.get(key) != quality_receipt.get(key):
            raise FM0ContractError(
                f"S{sector} relocated quality receipt disagrees on {key}"
            )

    provider = str(quality_receipt["mission_quality_type"])
    prefix = "spocffiflag" if provider == "spoc" else "ticaffiflag"
    quaternion: dict[tuple[int, int], frozenset[int]] = {}
    qflag: dict[tuple[int, int, int], Mapping[int, int]] = {}
    mission: dict[tuple[int, int], Mapping[int, int]] = {}
    for orbit in expected_orbits:
        for camera in range(1, 5):
            run = qlp_root / f"orbit-{orbit}/ffi/run"
            quaternion[(int(orbit), camera)] = frozenset(
                _read_quat(run / f"cam{camera}_quat.txt")
            )
            for ccd in range(1, 5):
                qflag[(int(orbit), camera, ccd)] = _read_flag_file(
                    run / f"cam{camera}ccd{ccd}_qflag.txt",
                    label="transferred QLP qflag authority",
                )
    for camera in range(1, 5):
        for ccd in range(1, 5):
            mission[(camera, ccd)] = _read_flag_file(
                mission_root / f"{prefix}_s{sector}_cam{camera}_ccd{ccd}.txt",
                label=f"transferred {provider.upper()} mission authority",
            )

    exclusions_payload = quality_receipt.get("mission_quality_authority_exclusions")
    if not isinstance(exclusions_payload, Mapping) or not isinstance(
        exclusions_payload.get("by_detector"), Mapping
    ):
        raise FM0ContractError(f"S{sector} quality exclusions are missing")
    exclusions: dict[tuple[int, int], frozenset[int]] = {}
    for camera in range(1, 5):
        for ccd in range(1, 5):
            entry = exclusions_payload["by_detector"].get(f"cam{camera}_ccd{ccd}")
            if not isinstance(entry, Mapping) or not isinstance(entry.get("rows"), list):
                raise FM0ContractError(f"S{sector} quality exclusion entry is invalid")
            values = frozenset(int(row["cadence"]) for row in entry["rows"])
            if len(values) != int(entry.get("n_rows", -1)):
                raise FM0ContractError(f"S{sector} quality exclusions are duplicated")
            exclusions[(camera, ccd)] = values
    return (
        Hdf5QualityAuthority(
            sector=sector,
            quaternion=quaternion,
            qflag=qflag,
            mission=mission,
            exclusions=exclusions,
        ),
        manifest,
        manifest_path,
        quality_receipt,
    )


def check_hdf5_product(
    *,
    path: str | Path,
    sector: int,
    orbit: int,
    camera: int,
    ccd: int,
    authority: Hdf5QualityAuthority,
) -> dict[str, int]:
    """Open and validate one retained HDF5 without loading photometry arrays."""

    import h5py

    resolved = Path(path).resolve()
    with h5py.File(resolved, "r") as handle:
        expected_attrs = {
            "Sector": int(sector),
            "Orbit": int(orbit),
            "Camera": int(camera),
            "CCD": int(ccd),
        }
        for name, expected in expected_attrs.items():
            if int(handle.attrs.get(name, -1)) != expected:
                raise FM0ContractError(f"{resolved} has wrong {name} attribute")
        tic_id = int(handle.attrs.get("TIC ID", -1))
        if tic_id <= 0 or resolved.stem != str(tic_id):
            raise FM0ContractError(f"{resolved} TIC identity disagrees with filename")
        if "LightCurve/Cadence" not in handle or "LightCurve/QualityFlag" not in handle:
            raise FM0ContractError(f"{resolved} lacks cadence or internal quality")
        cadence = np.asarray(handle["LightCurve/Cadence"][()])
        internal = np.asarray(handle["LightCurve/QualityFlag"][()])
    if (
        cadence.ndim != 1
        or internal.ndim != 1
        or cadence.size == 0
        or cadence.shape != internal.shape
        or not np.issubdtype(cadence.dtype, np.integer)
        or not np.issubdtype(internal.dtype, np.integer)
        or np.any(cadence <= 0)
        or np.any(internal < 0)
        or np.unique(cadence).size != cadence.size
    ):
        raise FM0ContractError(f"{resolved} cadence/internal-quality arrays are invalid")

    quat = authority.quaternion[(int(orbit), int(camera))]
    qflag = authority.qflag[(int(orbit), int(camera), int(ccd))]
    mission = authority.mission[(int(camera), int(ccd))]
    exclusions = authority.exclusions[(int(camera), int(ccd))]
    n_external_bad = 0
    n_excluded = 0
    for value in cadence:
        cad = int(value)
        if cad in exclusions:
            n_excluded += 1
            n_external_bad += 1
            continue
        if cad not in quat or cad not in qflag or cad not in mission:
            raise FM0ContractError(f"{resolved} has unexplained cadence {cad}")
        n_external_bad += int(qflag[cad] != 0 or mission[cad] != 0)
    return {
        "n_cadences": int(cadence.size),
        "n_internal_bad": int(np.count_nonzero(internal)),
        "n_external_bad": int(n_external_bad),
        "n_authority_excluded": int(n_excluded),
    }


_WORKER_AUTHORITY: Hdf5QualityAuthority | None = None


def _initialize_worker(authority: Hdf5QualityAuthority) -> None:
    global _WORKER_AUTHORITY
    _WORKER_AUTHORITY = authority


def _check_task(task: tuple[str, int, int, int, int]) -> dict[str, Any]:
    path, sector, orbit, camera, ccd = task
    try:
        if _WORKER_AUTHORITY is None:  # pragma: no cover - worker invariant
            raise RuntimeError("HDF5 quality worker lacks authority")
        counts = check_hdf5_product(
            path=path,
            sector=sector,
            orbit=orbit,
            camera=camera,
            ccd=ccd,
            authority=_WORKER_AUTHORITY,
        )
        return {"path": path, "passed": True, **counts}
    except (FM0ContractError, KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
        return {
            "path": path,
            "passed": False,
            "error": f"{type(exc).__name__}: {exc}",
        }


def audit_retained_sector_hdf5_quality(
    *,
    source_receipt_path: str | Path,
    quality_transfer_root: str | Path,
    sector: int,
    expected_orbits: Sequence[int],
    workers: int,
    producer_git_sha: str,
) -> dict[str, Any]:
    """Open every declared retained HDF5 and reconcile cadence quality."""

    if workers not in range(1, 17):
        raise FM0ContractError("workers must be in 1..16")
    if _GIT_SHA.fullmatch(producer_git_sha) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    orbits = tuple(int(value) for value in expected_orbits)
    if orbits != (2 * sector + 7, 2 * sector + 8):
        raise FM0ContractError(f"S{sector} expected orbit identity drifted")
    source_path = Path(source_receipt_path).expanduser().resolve()
    source = _read_json(source_path, label="ORCD source receipt")
    if (
        source.get("schema_version") != SOURCE_RECEIPT_SCHEMA_VERSION
        or source.get("source_state") != SOURCE_READY_STATE
        or source.get("source_form") != "retained_orcd_cells"
        or int(source.get("sector", -1)) != sector
        or tuple(int(value) for value in source.get("expected_orbits", ())) != orbits
        or int(source.get("n_cells", -1)) != 32
    ):
        raise FM0ContractError(f"S{sector} ORCD source receipt is invalid")

    transfer_root = Path(quality_transfer_root).expanduser().resolve()
    authority, transfer_manifest, transfer_manifest_path, quality_receipt = (
        _verify_transfer_sector(
            transfer_root=transfer_root,
            sector=sector,
            expected_orbits=orbits,
        )
    )
    tasks: list[tuple[str, int, int, int, int]] = []
    observed_paths: set[Path] = set()
    bindings = source.get("cell_bindings")
    if not isinstance(bindings, list) or len(bindings) != 32:
        raise FM0ContractError(f"S{sector} source receipt lacks 32 cell bindings")
    for binding in bindings:
        if not isinstance(binding, Mapping):
            raise FM0ContractError(f"S{sector} source cell binding is invalid")
        match = _CELL.fullmatch(str(binding.get("cell", "")))
        if match is None or int(match.group("sector")) != sector:
            raise FM0ContractError(f"S{sector} source cell identity is invalid")
        orbit = int(match.group("orbit"))
        camera = int(match.group("camera"))
        ccd = int(match.group("ccd"))
        root = Path(str(binding.get("retained_root", ""))).resolve()
        files = sorted((root / "LC").glob("*.h5"))
        if len(files) != int(binding.get("n_hdf5_products_declared", -1)):
            raise FM0ContractError(
                f"{binding.get('cell')} actual HDF5 count differs from receipt"
            )
        for path in files:
            if path in observed_paths:
                raise FM0ContractError(f"S{sector} HDF5 path appears in multiple cells")
            observed_paths.add(path)
            tasks.append((str(path), sector, orbit, camera, ccd))
    if len(tasks) != int(source.get("n_hdf5_products_declared", -1)):
        raise FM0ContractError(f"S{sector} sector HDF5 count differs from receipt")

    totals = {
        "n_hdf5_products": len(tasks),
        "n_hdf5_opened": 0,
        "n_unreadable_hdf5": 0,
        "n_cadences_checked": 0,
        "n_internal_bad_cadences": 0,
        "n_external_bad_cadences": 0,
        "n_authority_excluded_cadences": 0,
    }
    failures: list[dict[str, str]] = []
    if workers == 1:
        _initialize_worker(authority)
        results = map(_check_task, tasks)
        pool = None
    else:
        context = multiprocessing.get_context("fork")
        pool = ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
            initializer=_initialize_worker,
            initargs=(authority,),
        )
        results = pool.map(_check_task, tasks, chunksize=16)
    try:
        for index, result in enumerate(results, start=1):
            if result["passed"] is not True:
                totals["n_unreadable_hdf5"] += 1
                if len(failures) < 20:
                    failures.append(
                        {"path": str(result["path"]), "error": str(result["error"])}
                    )
            else:
                totals["n_hdf5_opened"] += 1
                totals["n_cadences_checked"] += int(result["n_cadences"])
                totals["n_internal_bad_cadences"] += int(result["n_internal_bad"])
                totals["n_external_bad_cadences"] += int(result["n_external_bad"])
                totals["n_authority_excluded_cadences"] += int(
                    result["n_authority_excluded"]
                )
            if index % 1000 == 0 or index == len(tasks):
                print(
                    f"[hdf5-quality] S{sector}: {index:,}/{len(tasks):,} "
                    f"opened={totals['n_hdf5_opened']:,} "
                    f"failed={totals['n_unreadable_hdf5']:,}",
                    flush=True,
                )
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    if failures or totals["n_hdf5_opened"] != totals["n_hdf5_products"]:
        raise FM0ContractError(
            f"S{sector} HDF5 quality admission failed: "
            f"n_failures={totals['n_unreadable_hdf5']}, examples={failures[:5]}"
        )
    return {
        "schema_version": HDF5_QUALITY_RECEIPT_SCHEMA_VERSION,
        "sector": sector,
        "expected_orbits": list(orbits),
        "quality_state": HDF5_QUALITY_READY_STATE,
        "passed": True,
        "producer_git_sha": producer_git_sha,
        "producer_code_path": str(Path(__file__).resolve()),
        "producer_code_sha256": sha256_file(Path(__file__).resolve()),
        "workers": workers,
        **totals,
        "source_receipt_path": str(source_path),
        "source_receipt_sha256": sha256_file(source_path),
        "quality_transfer_manifest_path": str(transfer_manifest_path),
        "quality_transfer_manifest_sha256": sha256_file(transfer_manifest_path),
        "quality_transfer_producer_git_sha": transfer_manifest["producer_git_sha"],
        "mission_quality_receipt_schema_version": quality_receipt["schema_version"],
        "mission_quality_type": quality_receipt["mission_quality_type"],
        "mission_quality_authority_exclusions_sha256": quality_receipt[
            "mission_quality_authority_exclusions_sha256"
        ],
        "hdf5_openability_verified": True,
        "internal_cadence_quality_verified": True,
        "external_cadence_quality_verified": True,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "claim_limit": (
            "HDF5 openability and cadence-quality join only; checksum-bound "
            "HDF5 provenance, six-view shards, identity, and panel gates remain"
        ),
    }


def write_hdf5_quality_receipt(receipt: Mapping[str, Any], output: str | Path) -> Path:
    if receipt.get("passed") is not True:
        raise FM0ContractError("refusing to publish a failing HDF5 quality receipt")
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise FM0ContractError(f"refusing to overwrite HDF5 quality receipt: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return target


__all__ = [
    "HDF5_QUALITY_READY_STATE",
    "HDF5_QUALITY_RECEIPT_SCHEMA_VERSION",
    "Hdf5QualityAuthority",
    "audit_retained_sector_hdf5_quality",
    "check_hdf5_product",
    "write_hdf5_quality_receipt",
]
