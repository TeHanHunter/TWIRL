"""Build a hash-bound cadence/quality authority table for A2v1 Tier-1 QA.

The QLP ``camC_quat.txt`` files are authoritative for cadence membership and
orbit assignment.  For Sector 56, the final external-quality mask combines
the detector-level SPOC bitmask with the binary QLP detector qflag.  This
module joins those inputs without inferring or filling cadences: every
detector in both quality sources must have exactly the cadence set declared by
its camera's quaternion files.

The precomputed SPOC-quality input must contain one row per detector cadence
with columns ``sector``, ``camera``, ``ccd``, ``cadenceno``, and ``quality``.
An optional ``orbitid`` column is checked against the quaternion-derived orbit
rather than trusted as an authority.  A mandatory provenance sidecar binds
that table to all original detector-level SPOC flag files.  The builder
rehashes and rereads those originals, so the sidecar is not accepted solely
on its declaration.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd


CADENCE_REFERENCE_BUILDER_VERSION = "a2v1_cadence_reference_builder_v2"
SPOC_QUALITY_TABLE_BUILDER_VERSION = "spoc_quality_table_builder_v1"
SPOC_QUALITY_DERIVATION_CONTRACT = "spocffiflag_detector_cadence_quality_v1"
SPOC_QUALITY_TABLE_COLUMNS = (
    "sector",
    "orbitid",
    "camera",
    "ccd",
    "cadenceno",
    "quality",
)
CADENCE_REFERENCE_COLUMNS = (
    "sector",
    "orbitid",
    "camera",
    "ccd",
    "cadenceno",
    "spoc_quality",
    "qlp_quality",
    "external_quality",
)
QLP_QUALITY_EXTERNAL_BIT = 30
S56_EXPECTED_ORBITS = (119, 120)
S56_EXPECTED_DETECTORS = tuple(
    (camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)
)


@dataclass(frozen=True)
class QuatSource:
    """One QLP camera quaternion file and its externally known orbit."""

    orbitid: int
    camera: int
    path: Path

    def __post_init__(self) -> None:
        if int(self.orbitid) <= 0:
            raise ValueError("quat orbitid must be positive")
        if int(self.camera) not in range(1, 5):
            raise ValueError("quat camera must be in 1..4")
        object.__setattr__(self, "orbitid", int(self.orbitid))
        object.__setattr__(self, "camera", int(self.camera))
        object.__setattr__(self, "path", Path(self.path))


@dataclass(frozen=True)
class SpocFlagSource:
    """One original detector-level SPOC flag file declared by a sidecar."""

    camera: int
    ccd: int
    path: Path
    sha256: str

    def __post_init__(self) -> None:
        if int(self.camera) not in range(1, 5) or int(self.ccd) not in range(1, 5):
            raise ValueError("SPOC flag source camera/ccd must be in 1..4")
        digest = str(self.sha256).lower()
        if len(digest) != 64 or any(
            value not in "0123456789abcdef" for value in digest
        ):
            raise ValueError(
                "SPOC flag source sha256 must be 64 hexadecimal characters"
            )
        object.__setattr__(self, "camera", int(self.camera))
        object.__setattr__(self, "ccd", int(self.ccd))
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(self, "sha256", digest)


@dataclass(frozen=True)
class SpocFlagInput:
    """One explicitly identified detector-level SPOC flag input."""

    camera: int
    ccd: int
    path: Path

    def __post_init__(self) -> None:
        if int(self.camera) not in range(1, 5) or int(self.ccd) not in range(1, 5):
            raise ValueError("SPOC flag input camera/ccd must be in 1..4")
        object.__setattr__(self, "camera", int(self.camera))
        object.__setattr__(self, "ccd", int(self.ccd))
        object.__setattr__(self, "path", Path(self.path))


@dataclass(frozen=True)
class QlpQflagSource:
    """One orbit/detector-level QLP binary qflag authority file."""

    orbitid: int
    camera: int
    ccd: int
    path: Path

    def __post_init__(self) -> None:
        if int(self.orbitid) <= 0:
            raise ValueError("QLP qflag orbitid must be positive")
        if int(self.camera) not in range(1, 5) or int(self.ccd) not in range(1, 5):
            raise ValueError("QLP qflag camera/ccd must be in 1..4")
        object.__setattr__(self, "orbitid", int(self.orbitid))
        object.__setattr__(self, "camera", int(self.camera))
        object.__setattr__(self, "ccd", int(self.ccd))
        object.__setattr__(self, "path", Path(self.path))


def file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest for one file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def parse_quat_spec(value: str) -> QuatSource:
    """Parse ``ORBIT,CAMERA,PATH`` for the stage-1 command-line interface."""

    parts = value.split(",", 2)
    if len(parts) != 3:
        raise ValueError("quat specification must be ORBIT,CAMERA,PATH")
    orbit, camera, path = parts
    if not path.strip():
        raise ValueError("quat specification has an empty path")
    return QuatSource(orbitid=int(orbit), camera=int(camera), path=Path(path))


def parse_spoc_flag_spec(value: str) -> SpocFlagInput:
    """Parse ``CAMERA,CCD,PATH`` for the S56 SPOC-quality builder CLI."""

    parts = value.split(",", 2)
    if len(parts) != 3:
        raise ValueError("SPOC flag specification must be CAMERA,CCD,PATH")
    camera, ccd, path = parts
    if not path.strip():
        raise ValueError("SPOC flag specification has an empty path")
    return SpocFlagInput(camera=int(camera), ccd=int(ccd), path=Path(path))


def parse_qlp_qflag_spec(value: str) -> QlpQflagSource:
    """Parse ``ORBIT,CAMERA,CCD,PATH`` for the final-reference CLI."""

    parts = value.split(",", 3)
    if len(parts) != 4:
        raise ValueError(
            "QLP qflag specification must be ORBIT,CAMERA,CCD,PATH"
        )
    orbit, camera, ccd, path = parts
    if not path.strip():
        raise ValueError("QLP qflag specification has an empty path")
    return QlpQflagSource(
        orbitid=int(orbit), camera=int(camera), ccd=int(ccd), path=Path(path)
    )


def _integer_value(value: object, *, context: str) -> int:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}: expected an integer, got {value!r}") from exc
    if not np.isfinite(number) or not number.is_integer():
        raise ValueError(f"{context}: expected an integer, got {value!r}")
    return int(number)


def read_qlp_quat_cadences(path: Path) -> np.ndarray:
    """Read the required ``cadence`` column from a QLP ``camC_quat.txt``.

    Unlike the older cadence-repair helper, this authority builder fails on a
    malformed nonblank row rather than silently skipping it.
    """

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    cadences: list[int] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{path}: empty quaternion file") from exc
        normalized = [str(value).strip().lower() for value in header]
        if "cadence" not in normalized:
            raise ValueError(f"{path}: no 'cadence' column in header {header!r}")
        cadence_column = normalized.index("cadence")
        for line_number, row in enumerate(reader, start=2):
            if not row or not any(str(value).strip() for value in row):
                continue
            if cadence_column >= len(row) or not str(row[cadence_column]).strip():
                raise ValueError(f"{path}:{line_number}: missing cadence value")
            cadence = _integer_value(
                row[cadence_column], context=f"{path}:{line_number} cadence"
            )
            if cadence < 0:
                raise ValueError(f"{path}:{line_number}: cadence must be nonnegative")
            cadences.append(cadence)
    if not cadences:
        raise ValueError(f"{path}: no cadences parsed")
    if len(cadences) != len(set(cadences)):
        duplicates = pd.Series(cadences)[pd.Series(cadences).duplicated()].unique()
        raise ValueError(
            f"{path}: duplicate cadence IDs: {duplicates[:10].astype(int).tolist()}"
        )
    return np.asarray(cadences, dtype=np.int64)


def read_spoc_flag_file(path: Path) -> dict[int, int]:
    """Read one original headerless ``cadenceno,quality`` SPOC flag file."""

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    values: dict[int, int] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle), start=1):
            if not row or not any(str(value).strip() for value in row):
                continue
            if len(row) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected cadenceno,quality; got {row!r}"
                )
            cadence = _integer_value(
                row[0], context=f"{path}:{line_number} cadenceno"
            )
            quality = _integer_value(
                row[1], context=f"{path}:{line_number} quality"
            )
            if cadence < 0 or quality < 0:
                raise ValueError(
                    f"{path}:{line_number}: cadence and quality must be nonnegative"
                )
            if cadence in values:
                raise ValueError(f"{path}:{line_number}: duplicate cadence {cadence}")
            values[cadence] = quality
    if not values:
        raise ValueError(f"{path}: no SPOC quality rows parsed")
    return values


def read_qlp_qflag_file(path: Path) -> dict[int, int]:
    """Read one headerless whitespace-delimited ``cadence flag`` QLP file."""

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    values: dict[int, int] = {}
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected whitespace-delimited "
                    f"cadence flag; got {line.rstrip()!r}"
                )
            cadence = _integer_value(
                fields[0], context=f"{path}:{line_number} cadence"
            )
            flag = _integer_value(fields[1], context=f"{path}:{line_number} flag")
            if cadence < 0:
                raise ValueError(f"{path}:{line_number}: cadence must be nonnegative")
            if flag not in {0, 1}:
                raise ValueError(
                    f"{path}:{line_number}: QLP qflag must be raw 0 or 1"
                )
            if cadence in values:
                raise ValueError(f"{path}:{line_number}: duplicate cadence {cadence}")
            values[cadence] = flag
    if not values:
        raise ValueError(f"{path}: no QLP qflag rows parsed")
    return values


def load_spoc_quality_provenance(
    path: Path,
    *,
    spoc_quality_table: Path,
    sector: int,
    expected_detectors: Sequence[tuple[int, int]],
) -> tuple[dict[str, Any], tuple[SpocFlagSource, ...]]:
    """Validate and rehash the mandatory SPOC-quality provenance sidecar."""

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("SPOC-quality provenance sidecar must contain a JSON object")
    required = {
        "contract_version",
        "derivation_contract",
        "sector",
        "quality_authority",
        "table_sha256",
        "n_rows",
        "detectors",
        "source_flag_files",
    }
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError(f"SPOC-quality provenance is missing fields: {missing}")
    expected_contract = f"s{int(sector)}_spoc_quality_table_provenance_v1"
    if str(payload["contract_version"]) != expected_contract:
        raise ValueError("SPOC-quality provenance contract mismatch")
    if str(payload["derivation_contract"]) != SPOC_QUALITY_DERIVATION_CONTRACT:
        raise ValueError("SPOC-quality derivation contract mismatch")
    if str(payload["quality_authority"]) != "spoc_quality_flags":
        raise ValueError("SPOC-quality provenance does not declare SPOC authority")
    try:
        manifest_sector = int(payload["sector"])
        manifest_rows = int(payload["n_rows"])
    except (TypeError, ValueError) as exc:
        raise ValueError("SPOC-quality sector/n_rows must be integers") from exc
    if manifest_sector != int(sector):
        raise ValueError("SPOC-quality provenance sector mismatch")
    table = _read_table(Path(spoc_quality_table))
    if manifest_rows != len(table):
        raise ValueError("SPOC-quality provenance row-count mismatch")
    table_digest = file_sha256(Path(spoc_quality_table))
    if str(payload["table_sha256"]).lower() != table_digest:
        raise ValueError("SPOC-quality provenance table hash mismatch")

    detector_pairs = _normalize_pairs(expected_detectors)
    expected_names = [f"cam{camera}_ccd{ccd}" for camera, ccd in detector_pairs]
    if not isinstance(payload["detectors"], list):
        raise ValueError("SPOC-quality provenance detectors must be a list")
    observed_names = sorted(str(value) for value in payload["detectors"])
    if observed_names != sorted(expected_names):
        raise ValueError(
            "SPOC-quality provenance detector inventory mismatch: "
            f"expected={expected_names}, observed={observed_names}"
        )
    if not isinstance(payload["source_flag_files"], list):
        raise ValueError("SPOC-quality provenance source_flag_files must be a list")
    sources: list[SpocFlagSource] = []
    for index, item in enumerate(payload["source_flag_files"]):
        if not isinstance(item, dict):
            raise ValueError(f"SPOC source entry {index} must be an object")
        item_missing = sorted({"camera", "ccd", "path", "sha256"} - set(item))
        if item_missing:
            raise ValueError(
                f"SPOC source entry {index} is missing fields: {item_missing}"
            )
        source_path = Path(str(item["path"]))
        if not source_path.is_absolute():
            raise ValueError(
                f"SPOC source entry {index} must record an absolute original path"
            )
        expected_name = (
            f"spocffiflag_s{int(sector)}_cam{int(item['camera'])}_"
            f"ccd{int(item['ccd'])}.txt"
        )
        if source_path.name != expected_name:
            raise ValueError(
                f"SPOC source entry {index} filename mismatch: "
                f"expected {expected_name}, observed {source_path.name}"
            )
        sources.append(
            SpocFlagSource(
                camera=int(item["camera"]),
                ccd=int(item["ccd"]),
                path=source_path.resolve(),
                sha256=str(item["sha256"]),
            )
        )
    observed_pairs = [(source.camera, source.ccd) for source in sources]
    if len(observed_pairs) != len(set(observed_pairs)):
        raise ValueError("SPOC-quality provenance has duplicate detector sources")
    if set(observed_pairs) != set(detector_pairs):
        missing_pairs = sorted(set(detector_pairs) - set(observed_pairs))
        extra_pairs = sorted(set(observed_pairs) - set(detector_pairs))
        raise ValueError(
            "SPOC-quality provenance source inventory is incomplete: "
            f"missing={missing_pairs}, extra={extra_pairs}"
        )
    resolved_paths = [source.path for source in sources]
    if len(resolved_paths) != len(set(resolved_paths)):
        raise ValueError(
            "SPOC-quality provenance reuses one file for multiple detectors"
        )
    for source in sources:
        if not source.path.is_file():
            raise FileNotFoundError(
                f"SPOC authority file is not accessible: {source.path}"
            )
        observed_digest = file_sha256(source.path)
        if observed_digest != source.sha256:
            raise ValueError(
                "SPOC authority file hash mismatch: "
                f"cam{source.camera}/ccd{source.ccd} {source.path}"
            )
    return payload, tuple(sorted(sources, key=lambda item: (item.camera, item.ccd)))


def _read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv" or path.name.endswith(".csv.gz"):
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"unsupported SPOC-quality table format: {path}")


def _integer_column(frame: pd.DataFrame, column: str) -> pd.Series:
    numeric = pd.to_numeric(frame[column], errors="raise")
    values = np.asarray(numeric, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"SPOC-quality column {column!r} contains non-finite values")
    if not np.all(values == np.floor(values)):
        raise ValueError(f"SPOC-quality column {column!r} must be integer-valued")
    return pd.Series(values.astype(np.int64), index=frame.index, name=column)


def _normalize_pairs(values: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    pairs = tuple(sorted({(int(camera), int(ccd)) for camera, ccd in values}))
    if not pairs:
        raise ValueError("expected_detectors cannot be empty")
    invalid = [
        pair
        for pair in pairs
        if pair[0] not in range(1, 5) or pair[1] not in range(1, 5)
    ]
    if invalid:
        raise ValueError(f"invalid expected detector pairs: {invalid}")
    return pairs


def build_spoc_quality_table(
    *,
    sector: int,
    quat_sources: Sequence[QuatSource],
    spoc_flag_sources: Sequence[SpocFlagInput],
    expected_orbits: Sequence[int],
    expected_detectors: Sequence[tuple[int, int]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Derive a detector quality table directly from QLP and SPOC authorities.

    QLP camera quaternion files define the exact accepted cadence set and its
    orbit assignment.  Each detector-level SPOC file supplies quality values.
    SPOC-only cadences are retained in the returned derivation metadata as
    explicit exclusions; missing quaternion cadences fail closed.
    """

    sector = int(sector)
    if sector <= 0:
        raise ValueError("sector must be positive")
    orbit_set = {int(value) for value in expected_orbits}
    if not orbit_set or any(value <= 0 for value in orbit_set):
        raise ValueError("expected_orbits must contain positive orbit IDs")
    detector_pairs = _normalize_pairs(expected_detectors)
    expected_cameras = {camera for camera, _ in detector_pairs}

    quats = tuple(quat_sources)
    quat_keys = [(source.orbitid, source.camera) for source in quats]
    if len(quat_keys) != len(set(quat_keys)):
        raise ValueError("duplicate quaternion source for an orbit/camera pair")
    expected_quat_keys = {
        (orbitid, camera) for orbitid in orbit_set for camera in expected_cameras
    }
    if set(quat_keys) != expected_quat_keys:
        missing = sorted(expected_quat_keys - set(quat_keys))
        extra = sorted(set(quat_keys) - expected_quat_keys)
        raise ValueError(
            f"quaternion source inventory mismatch: missing={missing}, extra={extra}"
        )
    quat_paths = [source.path.resolve() for source in quats]
    if len(quat_paths) != len(set(quat_paths)):
        raise ValueError("one quaternion file cannot represent multiple inputs")

    cadence_orbits: dict[tuple[int, int], int] = {}
    quat_summaries: list[dict[str, Any]] = []
    for source in sorted(quats, key=lambda item: (item.orbitid, item.camera)):
        resolved = source.path.resolve()
        cadences = read_qlp_quat_cadences(resolved)
        for cadence_value in cadences:
            cadence = int(cadence_value)
            key = (source.camera, cadence)
            if key in cadence_orbits:
                raise ValueError(
                    "quaternion files assign one camera cadence to multiple orbits: "
                    f"camera={source.camera}, cadenceno={cadence}, "
                    f"orbits={cadence_orbits[key]},{source.orbitid}"
                )
            cadence_orbits[key] = source.orbitid
        quat_summaries.append(
            {
                "orbitid": source.orbitid,
                "camera": source.camera,
                "path": str(resolved),
                "n_cadences": int(len(cadences)),
            }
        )

    flags = tuple(spoc_flag_sources)
    flag_keys = [(source.camera, source.ccd) for source in flags]
    if len(flag_keys) != len(set(flag_keys)):
        raise ValueError("duplicate SPOC flag source for a detector")
    if set(flag_keys) != set(detector_pairs):
        missing = sorted(set(detector_pairs) - set(flag_keys))
        extra = sorted(set(flag_keys) - set(detector_pairs))
        raise ValueError(
            f"SPOC flag source inventory mismatch: missing={missing}, extra={extra}"
        )
    flag_paths = [source.path.resolve() for source in flags]
    if len(flag_paths) != len(set(flag_paths)):
        raise ValueError("one SPOC flag file cannot represent multiple detectors")
    if set(flag_paths) & set(quat_paths):
        raise ValueError("one authority file cannot be both quaternion and SPOC input")

    rows: list[dict[str, int]] = []
    flag_summaries: list[dict[str, Any]] = []
    for source in sorted(flags, key=lambda item: (item.camera, item.ccd)):
        resolved = source.path.resolve()
        expected_name = (
            f"spocffiflag_s{sector}_cam{source.camera}_ccd{source.ccd}.txt"
        )
        if resolved.name != expected_name:
            raise ValueError(
                "SPOC flag filename mismatch: "
                f"expected {expected_name}, observed {resolved.name}"
            )
        raw_quality = read_spoc_flag_file(resolved)
        expected_cadences = {
            cadence: orbitid
            for (camera, cadence), orbitid in cadence_orbits.items()
            if camera == source.camera
        }
        if not expected_cadences:
            raise ValueError(f"no quaternion cadences for camera {source.camera}")
        missing_cadences = sorted(set(expected_cadences) - set(raw_quality))
        if missing_cadences:
            raise ValueError(
                "SPOC authority is missing quaternion cadences: "
                f"cam{source.camera}/ccd{source.ccd}, "
                f"n_missing={len(missing_cadences)}, "
                f"examples={missing_cadences[:20]}"
            )
        excluded = sorted(set(raw_quality) - set(expected_cadences))
        for cadence, orbitid in sorted(
            expected_cadences.items(), key=lambda item: (item[1], item[0])
        ):
            rows.append(
                {
                    "sector": sector,
                    "orbitid": int(orbitid),
                    "camera": source.camera,
                    "ccd": source.ccd,
                    "cadenceno": int(cadence),
                    "quality": int(raw_quality[cadence]),
                }
            )
        flag_summaries.append(
            {
                "camera": source.camera,
                "ccd": source.ccd,
                "path": str(resolved),
                "n_rows_raw": int(len(raw_quality)),
                "n_rows_retained": int(len(expected_cadences)),
                "n_rows_excluded_by_quat": int(len(excluded)),
                "excluded_rows": [
                    {"cadenceno": int(cadence), "quality": int(raw_quality[cadence])}
                    for cadence in excluded
                ],
            }
        )

    frame = pd.DataFrame(rows, columns=SPOC_QUALITY_TABLE_COLUMNS)
    if frame.empty:
        raise ValueError("derived SPOC quality table is empty")
    for column in SPOC_QUALITY_TABLE_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype(np.int64)
    frame = frame.sort_values(
        ["orbitid", "camera", "ccd", "cadenceno"], kind="stable"
    ).reset_index(drop=True)
    if frame.duplicated(["camera", "ccd", "cadenceno"]).any():
        raise RuntimeError("internal error: duplicate detector cadence in quality table")
    observed_detectors = set(
        zip(frame["camera"].astype(int), frame["ccd"].astype(int), strict=True)
    )
    if observed_detectors != set(detector_pairs):
        raise RuntimeError("internal error: derived detector inventory is incomplete")
    if set(frame["orbitid"].astype(int)) != orbit_set:
        raise RuntimeError("internal error: derived orbit inventory is incomplete")

    summary = {
        "quat_sources": quat_summaries,
        "spoc_flag_sources": flag_summaries,
        "n_rows_raw": int(sum(item["n_rows_raw"] for item in flag_summaries)),
        "n_rows_retained": int(len(frame)),
        "n_rows_excluded_by_quat": int(
            sum(item["n_rows_excluded_by_quat"] for item in flag_summaries)
        ),
    }
    return frame, summary


def build_cadence_reference(
    *,
    sector: int,
    quat_sources: Sequence[QuatSource],
    qlp_qflag_sources: Sequence[QlpQflagSource],
    spoc_quality_table: Path,
    spoc_quality_provenance: Path,
    expected_orbits: Sequence[int],
    expected_detectors: Sequence[tuple[int, int]],
) -> pd.DataFrame:
    """Construct the exact Tier-1 cadence-reference table.

    The join is fail-closed: missing detector cadences, extra detector
    cadences, duplicate rows, duplicate camera-cadence orbit assignments, any
    supplied ``orbitid`` disagreement, and non-exact QLP qflag coverage all
    raise ``ValueError``.
    """

    sector = int(sector)
    if sector <= 0:
        raise ValueError("sector must be positive")
    expected_orbit_set = {int(value) for value in expected_orbits}
    if not expected_orbit_set or any(value <= 0 for value in expected_orbit_set):
        raise ValueError("expected_orbits must contain positive orbit IDs")
    detector_pairs = _normalize_pairs(expected_detectors)
    expected_cameras = {camera for camera, _ in detector_pairs}
    _, spoc_sources = load_spoc_quality_provenance(
        spoc_quality_provenance,
        spoc_quality_table=spoc_quality_table,
        sector=sector,
        expected_detectors=detector_pairs,
    )

    sources = tuple(quat_sources)
    source_keys = [(source.orbitid, source.camera) for source in sources]
    if len(source_keys) != len(set(source_keys)):
        raise ValueError("duplicate quaternion source for an orbit/camera pair")
    expected_source_keys = {
        (orbitid, camera)
        for orbitid in expected_orbit_set
        for camera in expected_cameras
    }
    observed_source_keys = set(source_keys)
    if observed_source_keys != expected_source_keys:
        missing = sorted(expected_source_keys - observed_source_keys)
        extra = sorted(observed_source_keys - expected_source_keys)
        raise ValueError(
            f"quaternion source inventory mismatch: missing={missing}, extra={extra}"
        )

    cadence_orbits: dict[tuple[int, int], int] = {}
    for source in sorted(sources, key=lambda item: (item.orbitid, item.camera)):
        for cadence in read_qlp_quat_cadences(source.path):
            key = (source.camera, int(cadence))
            if key in cadence_orbits:
                raise ValueError(
                    "quaternion files assign one camera cadence to multiple orbits: "
                    f"camera={source.camera}, cadenceno={int(cadence)}, "
                    f"orbits={cadence_orbits[key]},{source.orbitid}"
                )
            cadence_orbits[key] = source.orbitid

    qflag_sources = tuple(qlp_qflag_sources)
    qflag_keys = [
        (source.orbitid, source.camera, source.ccd) for source in qflag_sources
    ]
    if len(qflag_keys) != len(set(qflag_keys)):
        raise ValueError("duplicate QLP qflag source for an orbit/detector")
    expected_qflag_keys = {
        (orbitid, camera, ccd)
        for orbitid in expected_orbit_set
        for camera, ccd in detector_pairs
    }
    if set(qflag_keys) != expected_qflag_keys:
        missing = sorted(expected_qflag_keys - set(qflag_keys))
        extra = sorted(set(qflag_keys) - expected_qflag_keys)
        raise ValueError(
            f"QLP qflag source inventory mismatch: missing={missing}, extra={extra}"
        )
    qflag_paths = [source.path.resolve() for source in qflag_sources]
    if len(qflag_paths) != len(set(qflag_paths)):
        raise ValueError("one QLP qflag file cannot represent multiple inputs")
    quat_paths = [source.path.resolve() for source in sources]
    spoc_paths = [source.path.resolve() for source in spoc_sources]
    if set(qflag_paths) & (set(quat_paths) | set(spoc_paths)):
        raise ValueError(
            "one authority file cannot represent both a QLP qflag and another input"
        )

    qlp_quality_by_key: dict[tuple[int, int, int, int], int] = {}
    for source in sorted(
        qflag_sources, key=lambda item: (item.orbitid, item.camera, item.ccd)
    ):
        raw_qflag = read_qlp_qflag_file(source.path)
        expected_cadences = {
            cadence
            for (camera, cadence), orbitid in cadence_orbits.items()
            if camera == source.camera and orbitid == source.orbitid
        }
        missing_cadences = sorted(expected_cadences - set(raw_qflag))
        extra_cadences = sorted(set(raw_qflag) - expected_cadences)
        if missing_cadences or extra_cadences:
            raise ValueError(
                "QLP qflag cadence coverage is not exact: "
                f"orbit={source.orbitid}, cam{source.camera}/ccd{source.ccd}, "
                f"missing={missing_cadences[:20]}, extra={extra_cadences[:20]}"
            )
        for cadence in sorted(expected_cadences):
            qlp_quality_by_key[
                (source.orbitid, source.camera, source.ccd, cadence)
            ] = raw_qflag[cadence]

    quality = _read_table(Path(spoc_quality_table)).copy()
    required = {"sector", "camera", "ccd", "cadenceno", "quality"}
    missing_columns = sorted(required - set(quality.columns))
    if missing_columns:
        raise KeyError(f"SPOC-quality table is missing columns: {missing_columns}")
    for column in sorted(required | ({"orbitid"} & set(quality.columns))):
        quality[column] = _integer_column(quality, column)
    if not len(quality):
        raise ValueError("SPOC-quality table is empty")
    if set(quality["sector"].astype(int)) != {sector}:
        raise ValueError(
            f"SPOC-quality sector mismatch: expected {sector}, "
            f"observed {sorted(set(quality['sector'].astype(int)))}"
        )
    if (quality["cadenceno"] < 0).any():
        raise ValueError("SPOC-quality cadenceno values must be nonnegative")
    if (quality["quality"] < 0).any():
        raise ValueError("SPOC quality bitmasks must be nonnegative")
    quality = quality.rename(columns={"quality": "spoc_quality"})
    key_columns = ["sector", "camera", "ccd", "cadenceno"]
    duplicate_mask = quality.duplicated(key_columns, keep=False)
    if duplicate_mask.any():
        examples = quality.loc[duplicate_mask, key_columns].head(10).to_dict("records")
        raise ValueError(
            f"SPOC-quality table has duplicate detector cadences: {examples}"
        )
    observed_detectors = set(
        zip(quality["camera"].astype(int), quality["ccd"].astype(int), strict=True)
    )
    if observed_detectors != set(detector_pairs):
        missing = sorted(set(detector_pairs) - observed_detectors)
        extra = sorted(observed_detectors - set(detector_pairs))
        raise ValueError(
            "SPOC-quality detector inventory mismatch: "
            f"missing={missing}, extra={extra}"
        )

    expected_keys: list[dict[str, int]] = []
    for camera, ccd in detector_pairs:
        camera_cadences = sorted(
            (cadence, orbitid)
            for (candidate_camera, cadence), orbitid in cadence_orbits.items()
            if candidate_camera == camera
        )
        if not camera_cadences:
            raise ValueError(f"no quaternion cadences for camera {camera}")
        expected_keys.extend(
            {
                "sector": sector,
                "orbitid": orbitid,
                "camera": camera,
                "ccd": ccd,
                "cadenceno": cadence,
                "quat_orbitid": orbitid,
                "qlp_quality": qlp_quality_by_key[
                    (orbitid, camera, ccd, cadence)
                ],
            }
            for cadence, orbitid in camera_cadences
        )
    expected = pd.DataFrame(expected_keys)
    joined = expected.merge(
        quality,
        how="outer",
        on=key_columns,
        validate="one_to_one",
        indicator=True,
        suffixes=("", "_spoc"),
    )
    missing_quality = joined.loc[joined["_merge"] == "left_only", key_columns]
    extra_quality = joined.loc[joined["_merge"] == "right_only", key_columns]
    if len(missing_quality) or len(extra_quality):
        raise ValueError(
            "quaternion/SPOC cadence join is not exact: "
            f"missing_quality={len(missing_quality)} "
            f"examples={missing_quality.head(5).to_dict('records')}; "
            f"extra_quality={len(extra_quality)} "
            f"examples={extra_quality.head(5).to_dict('records')}"
        )
    if "orbitid_spoc" in joined and not np.array_equal(
        joined["orbitid_spoc"].to_numpy(dtype=np.int64),
        joined["quat_orbitid"].to_numpy(dtype=np.int64),
    ):
        mismatches = joined.loc[
            joined["orbitid_spoc"] != joined["quat_orbitid"],
            [*key_columns, "orbitid_spoc", "quat_orbitid"],
        ]
        raise ValueError(
            "SPOC-quality orbitid disagrees with QLP quaternion authority: "
            f"{mismatches.head(10).to_dict('records')}"
        )
    if "orbitid_spoc" in joined:
        joined = joined.drop(columns="orbitid_spoc")

    result = joined.drop(columns=["quat_orbitid", "_merge"])
    result["external_quality"] = np.bitwise_or(
        result["spoc_quality"].to_numpy(dtype=np.int64),
        np.left_shift(
            result["qlp_quality"].to_numpy(dtype=np.int64),
            QLP_QUALITY_EXTERNAL_BIT,
        ),
    )
    result = result.loc[:, list(CADENCE_REFERENCE_COLUMNS)].copy()
    for column in CADENCE_REFERENCE_COLUMNS:
        result[column] = pd.to_numeric(result[column], errors="raise").astype(np.int64)
    result = result.sort_values(
        ["orbitid", "camera", "ccd", "cadenceno"], kind="stable"
    ).reset_index(drop=True)
    if result.duplicated(["camera", "ccd", "cadenceno"]).any():
        raise RuntimeError("internal error: duplicate detector cadence in output")
    for source in spoc_sources:
        raw_quality = read_spoc_flag_file(source.path)
        detector = result.loc[
            (result["camera"] == source.camera) & (result["ccd"] == source.ccd),
            ["cadenceno", "spoc_quality"],
        ]
        expected_quality = {
            int(row.cadenceno): int(row.spoc_quality)
            for row in detector.itertuples(index=False)
        }
        missing_raw = sorted(set(expected_quality) - set(raw_quality))
        mismatched = sorted(
            cadence
            for cadence in set(expected_quality) & set(raw_quality)
            if expected_quality[cadence] != raw_quality[cadence]
        )
        if missing_raw or mismatched:
            raise ValueError(
                "precomputed SPOC quality does not match its original authority: "
                f"cam{source.camera}/ccd{source.ccd}, "
                f"missing_raw={missing_raw[:10]}, "
                f"quality_mismatch={mismatched[:10]}"
            )
    return result


def _temporary_path(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.stem}.", suffix=path.suffix
    )
    os.close(descriptor)
    os.chmod(temporary, 0o640)
    return Path(temporary)


def _write_table(frame: pd.DataFrame, path: Path) -> None:
    if path.suffix.lower() == ".csv":
        frame.to_csv(path, index=False, lineterminator="\n")
    elif path.suffix.lower() == ".parquet":
        frame.to_parquet(path, index=False)
    else:
        raise ValueError("output table must end in .csv or .parquet")
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def _write_json(payload: dict[str, Any], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(Path(path), os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publish_pair_with_rollback(
    *,
    table_temporary: Path,
    manifest_temporary: Path,
    output_table: Path,
    output_manifest: Path,
) -> None:
    """Publish two prepared files and restore the prior pair on an exception."""

    outputs = (Path(output_table), Path(output_manifest))
    backups: dict[Path, Path] = {}
    existed = {path: path.exists() for path in outputs}
    try:
        for path in outputs:
            if not existed[path]:
                continue
            backup = _temporary_path(path)
            try:
                shutil.copyfile(path, backup)
            except Exception:
                backup.unlink(missing_ok=True)
                raise
            with backup.open("rb") as handle:
                os.fsync(handle.fileno())
            backups[path] = backup
        os.replace(table_temporary, output_table)
        os.replace(manifest_temporary, output_manifest)
        for parent in {path.parent for path in outputs}:
            _fsync_directory(parent)
    except Exception:
        for path in outputs:
            backup = backups.get(path)
            if backup is not None and backup.exists():
                os.replace(backup, path)
            elif not existed[path]:
                path.unlink(missing_ok=True)
        for parent in {path.parent for path in outputs}:
            _fsync_directory(parent)
        raise
    finally:
        for backup in backups.values():
            backup.unlink(missing_ok=True)


def write_spoc_quality_table(
    *,
    sector: int,
    quat_sources: Sequence[QuatSource],
    spoc_flag_sources: Sequence[SpocFlagInput],
    output_table: Path,
    output_provenance: Path,
    expected_orbits: Sequence[int],
    expected_detectors: Sequence[tuple[int, int]],
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write a hash-bound SPOC pair with rollback on publication errors."""

    output_table = Path(output_table)
    output_provenance = Path(output_provenance)
    if output_table.suffix.lower() != ".csv":
        raise ValueError("SPOC quality output table must end in .csv")
    if output_provenance.suffix.lower() != ".json":
        raise ValueError("SPOC quality provenance must end in .json")
    if output_table.resolve() == output_provenance.resolve():
        raise ValueError("SPOC quality table and provenance outputs must differ")

    quats = tuple(quat_sources)
    flags = tuple(spoc_flag_sources)
    inputs = [*(source.path for source in quats), *(source.path for source in flags)]
    resolved_inputs = [path.resolve() for path in inputs]
    if len(resolved_inputs) != len(set(resolved_inputs)):
        raise ValueError("authority input paths must be unique")
    output_paths = {output_table.resolve(), output_provenance.resolve()}
    collisions = sorted(str(path) for path in set(resolved_inputs) & output_paths)
    if collisions:
        raise ValueError(f"output path collides with an authority input: {collisions}")
    existing = [path for path in (output_table, output_provenance) if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"output already exists: {existing}")

    initial_source_hashes = {
        str(path): file_sha256(path) for path in sorted(set(resolved_inputs))
    }
    frame, derivation = build_spoc_quality_table(
        sector=sector,
        quat_sources=quats,
        spoc_flag_sources=flags,
        expected_orbits=expected_orbits,
        expected_detectors=expected_detectors,
    )
    table_temporary = _temporary_path(output_table)
    provenance_temporary = _temporary_path(output_provenance)
    try:
        _write_table(frame, table_temporary)
        table_sha256 = file_sha256(table_temporary)

        source_quat_files: list[dict[str, Any]] = []
        for item in derivation["quat_sources"]:
            source_path = Path(str(item["path"])).resolve()
            source_quat_files.append(
                {
                    **item,
                    "path": str(source_path),
                    "sha256": initial_source_hashes[str(source_path)],
                }
            )
        source_flag_files: list[dict[str, Any]] = []
        for item in derivation["spoc_flag_sources"]:
            source_path = Path(str(item["path"])).resolve()
            source_flag_files.append(
                {
                    **item,
                    "path": str(source_path),
                    "sha256": initial_source_hashes[str(source_path)],
                }
            )
        detector_pairs = _normalize_pairs(expected_detectors)
        excluded_by_detector = {
            f"cam{int(item['camera'])}_ccd{int(item['ccd'])}": {
                "n_rows": int(item["n_rows_excluded_by_quat"]),
                "rows": item["excluded_rows"],
            }
            for item in source_flag_files
        }
        manifest: dict[str, Any] = {
            "contract_version": f"s{int(sector)}_spoc_quality_table_provenance_v1",
            "builder_version": SPOC_QUALITY_TABLE_BUILDER_VERSION,
            "derivation_contract": SPOC_QUALITY_DERIVATION_CONTRACT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "sector": int(sector),
            "cadence_authority": "qlp_cam_quat",
            "quality_authority": "spoc_quality_flags",
            "table_path": str(output_table.resolve()),
            "table_sha256": table_sha256,
            "table_columns": list(SPOC_QUALITY_TABLE_COLUMNS),
            "n_rows": int(len(frame)),
            "n_rows_raw": int(derivation["n_rows_raw"]),
            "n_rows_retained": int(derivation["n_rows_retained"]),
            "n_rows_excluded_by_quat": int(
                derivation["n_rows_excluded_by_quat"]
            ),
            "n_quaternion_files": len(source_quat_files),
            "n_spoc_flag_files": len(source_flag_files),
            "detectors": [
                f"cam{camera}_ccd{ccd}" for camera, ccd in detector_pairs
            ],
            "orbits": sorted({int(value) for value in expected_orbits}),
            "n_rows_by_detector": {
                f"cam{int(camera)}_ccd{int(ccd)}": int(len(group))
                for (camera, ccd), group in frame.groupby(
                    ["camera", "ccd"], sort=True
                )
            },
            "n_rows_by_orbit": {
                str(int(orbitid)): int(len(group))
                for orbitid, group in frame.groupby("orbitid", sort=True)
            },
            "n_nonzero_quality": int(
                np.count_nonzero(frame["quality"].to_numpy(dtype=np.int64))
            ),
            "exclusions": {
                "policy": "exclude SPOC cadences absent from the explicit QLP quaternion authorities",
                "n_rows": int(derivation["n_rows_excluded_by_quat"]),
                "by_detector": excluded_by_detector,
            },
            "source_quat_files": source_quat_files,
            "source_flag_files": source_flag_files,
            "source_file_sha256": dict(sorted(initial_source_hashes.items())),
        }
        _write_json(manifest, provenance_temporary)

        # Validate the exact temporary bytes against the existing downstream
        # loader before either output is published.
        load_spoc_quality_provenance(
            provenance_temporary,
            spoc_quality_table=table_temporary,
            sector=sector,
            expected_detectors=detector_pairs,
        )
        final_source_hashes = {
            str(path): file_sha256(path) for path in sorted(set(resolved_inputs))
        }
        if final_source_hashes != initial_source_hashes:
            raise RuntimeError("an authority input changed while the table was built")

        _publish_pair_with_rollback(
            table_temporary=table_temporary,
            manifest_temporary=provenance_temporary,
            output_table=output_table,
            output_manifest=output_provenance,
        )
        if file_sha256(output_table) != table_sha256:
            raise RuntimeError("published SPOC quality table hash changed")
        return manifest
    finally:
        table_temporary.unlink(missing_ok=True)
        provenance_temporary.unlink(missing_ok=True)


def write_cadence_reference(
    *,
    sector: int,
    quat_sources: Sequence[QuatSource],
    qlp_qflag_sources: Sequence[QlpQflagSource],
    spoc_quality_table: Path,
    spoc_quality_provenance: Path,
    output_table: Path,
    output_manifest: Path,
    expected_orbits: Sequence[int],
    expected_detectors: Sequence[tuple[int, int]],
    overwrite: bool = False,
) -> dict[str, Any]:
    """Build and publish the table/manifest pair with exception rollback."""

    output_table = Path(output_table)
    output_manifest = Path(output_manifest)
    spoc_quality_table = Path(spoc_quality_table)
    spoc_quality_provenance = Path(spoc_quality_provenance)
    if output_manifest.suffix.lower() != ".json":
        raise ValueError("output manifest must end in .json")
    provenance_hash_before_load = file_sha256(spoc_quality_provenance)
    provenance_payload, spoc_sources = load_spoc_quality_provenance(
        spoc_quality_provenance,
        spoc_quality_table=spoc_quality_table,
        sector=sector,
        expected_detectors=expected_detectors,
    )
    inputs = [
        spoc_quality_table,
        spoc_quality_provenance,
        *(source.path for source in quat_sources),
        *(source.path for source in qlp_qflag_sources),
        *(source.path for source in spoc_sources),
    ]
    resolved_inputs = [path.resolve() for path in inputs]
    input_paths = set(resolved_inputs)
    if len(resolved_inputs) != len(input_paths):
        raise ValueError("cadence-reference authority input paths must be unique")
    quat_paths = [source.path.resolve() for source in quat_sources]
    if len(quat_paths) != len(set(quat_paths)):
        raise ValueError(
            "one quaternion file cannot represent multiple orbit/camera inputs"
        )
    qflag_paths = [source.path.resolve() for source in qlp_qflag_sources]
    if len(qflag_paths) != len(set(qflag_paths)):
        raise ValueError("one QLP qflag file cannot represent multiple inputs")
    output_paths = {output_table.resolve(), output_manifest.resolve()}
    if len(output_paths) != 2:
        raise ValueError("output table and output manifest must be different files")
    collisions = sorted(str(path) for path in input_paths & output_paths)
    if collisions:
        raise ValueError(f"output path collides with an input: {collisions}")
    existing = [path for path in (output_table, output_manifest) if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"output already exists: {existing}")

    initial_source_hashes = {
        str(path.resolve()): file_sha256(path.resolve()) for path in inputs
    }
    if (
        initial_source_hashes[str(spoc_quality_provenance.resolve())]
        != provenance_hash_before_load
    ):
        raise RuntimeError("the SPOC provenance sidecar changed while it was loaded")
    frame = build_cadence_reference(
        sector=sector,
        quat_sources=quat_sources,
        qlp_qflag_sources=qlp_qflag_sources,
        spoc_quality_table=spoc_quality_table,
        spoc_quality_provenance=spoc_quality_provenance,
        expected_orbits=expected_orbits,
        expected_detectors=expected_detectors,
    )
    table_temporary = _temporary_path(output_table)
    manifest_temporary = _temporary_path(output_manifest)
    try:
        _write_table(frame, table_temporary)
        table_hash = file_sha256(table_temporary)
        final_source_hashes = {
            str(path.resolve()): file_sha256(path.resolve()) for path in inputs
        }
        if final_source_hashes != initial_source_hashes:
            raise RuntimeError(
                "an authority input changed while the reference was built"
            )
        sources: list[dict[str, Any]] = []
        ordered_sources = sorted(
            quat_sources, key=lambda item: (item.orbitid, item.camera)
        )
        for source in ordered_sources:
            resolved = source.path.resolve()
            digest = initial_source_hashes[str(resolved)]
            sources.append(
                {
                    "role": "qlp_cam_quat",
                    "orbitid": source.orbitid,
                    "camera": source.camera,
                    "path": str(resolved),
                    "sha256": digest,
                }
            )
        ordered_qflag_sources = sorted(
            qlp_qflag_sources,
            key=lambda item: (item.orbitid, item.camera, item.ccd),
        )
        n_qflag_rows = 0
        for source in ordered_qflag_sources:
            resolved = source.path.resolve()
            digest = initial_source_hashes[str(resolved)]
            n_rows = len(read_qlp_qflag_file(resolved))
            n_qflag_rows += n_rows
            sources.append(
                {
                    "role": "qlp_detector_qflag",
                    "orbitid": source.orbitid,
                    "camera": source.camera,
                    "ccd": source.ccd,
                    "path": str(resolved),
                    "sha256": digest,
                    "n_rows": n_rows,
                }
            )
        quality_resolved = spoc_quality_table.resolve()
        quality_hash = initial_source_hashes[str(quality_resolved)]
        sources.append(
            {
                "role": "spoc_quality_table",
                "path": str(quality_resolved),
                "sha256": quality_hash,
            }
        )
        provenance_resolved = spoc_quality_provenance.resolve()
        provenance_hash = initial_source_hashes[str(provenance_resolved)]
        sources.append(
            {
                "role": "spoc_quality_provenance",
                "path": str(provenance_resolved),
                "sha256": provenance_hash,
            }
        )
        n_raw_rows = 0
        for source in spoc_sources:
            digest = initial_source_hashes[str(source.path)]
            raw_rows = len(read_spoc_flag_file(source.path))
            n_raw_rows += raw_rows
            sources.append(
                {
                    "role": "spoc_flag_file",
                    "camera": source.camera,
                    "ccd": source.ccd,
                    "path": str(source.path),
                    "sha256": digest,
                    "n_rows": raw_rows,
                }
            )
        detector_pairs = _normalize_pairs(expected_detectors)
        manifest: dict[str, Any] = {
            "contract_version": f"s{int(sector)}_a2v1_cadence_reference_v1",
            "builder_version": CADENCE_REFERENCE_BUILDER_VERSION,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "sector": int(sector),
            "cadence_authority": "qlp_cam_quat",
            "quality_authority": "spoc_and_qlp_quality_flags",
            "quality_composition": {
                "external_quality": "spoc_quality | (qlp_quality << 30)",
                "qlp_quality_raw_values": [0, 1],
                "qlp_quality_external_bit": QLP_QUALITY_EXTERNAL_BIT,
            },
            "spoc_quality_derivation_contract": provenance_payload[
                "derivation_contract"
            ],
            "spoc_quality_provenance_sha256": provenance_hash,
            "n_spoc_authority_files_verified": len(spoc_sources),
            "n_qlp_qflag_files_verified": len(ordered_qflag_sources),
            "n_qlp_qflag_rows": n_qflag_rows,
            "n_spoc_raw_rows": n_raw_rows,
            "n_spoc_rows_excluded_by_quat": n_raw_rows - int(len(frame)),
            "table_sha256": table_hash,
            "table_columns": list(CADENCE_REFERENCE_COLUMNS),
            "n_rows": int(len(frame)),
            "detectors": [f"cam{camera}_ccd{ccd}" for camera, ccd in detector_pairs],
            "orbits": sorted({int(value) for value in expected_orbits}),
            "n_rows_by_detector": {
                f"cam{int(camera)}_ccd{int(ccd)}": int(len(group))
                for (camera, ccd), group in frame.groupby(["camera", "ccd"], sort=True)
            },
            "n_nonzero_spoc_quality": int(
                np.count_nonzero(frame["spoc_quality"].to_numpy())
            ),
            "n_nonzero_qlp_quality": int(
                np.count_nonzero(frame["qlp_quality"].to_numpy())
            ),
            "n_nonzero_external_quality": int(
                np.count_nonzero(frame["external_quality"].to_numpy())
            ),
            "source_file_sha256": dict(sorted(initial_source_hashes.items())),
            "sources": sources,
        }
        _write_json(manifest, manifest_temporary)
        final_source_hashes = {
            str(path.resolve()): file_sha256(path.resolve()) for path in inputs
        }
        if final_source_hashes != initial_source_hashes:
            raise RuntimeError(
                "an authority input changed while the manifest was built"
            )
        _publish_pair_with_rollback(
            table_temporary=table_temporary,
            manifest_temporary=manifest_temporary,
            output_table=output_table,
            output_manifest=output_manifest,
        )
        return manifest
    finally:
        table_temporary.unlink(missing_ok=True)
        manifest_temporary.unlink(missing_ok=True)
