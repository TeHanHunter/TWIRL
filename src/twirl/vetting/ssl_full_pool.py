"""Freeze the checksum-bound real-light-curve pool for Teacher SSL.

The pool is observation keyed: a TIC observed in more than one sector keeps
one row per ``(sector, tic)``.  Exclusions are deliberately host keyed.  Every
observation of a TIC is removed when that TIC belongs to either the immutable
Teacher-v3 fixed test split or the prospective Sector 63 reservation.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import time
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.teacher_split_registry import (
    REGISTRY_COLUMNS,
    validate_tic_split_registry,
)


FULL_POOL_CONTRACT_VERSION = "twirl_teacher_ssl_full_pool_v1"
FULL_POOL_SUMMARY_SCHEMA_VERSION = "twirl_teacher_ssl_full_pool_summary_v1"
FULL_POOL_SELECTION_POLICY_VERSION = (
    "whole_tic_exclude_teacher_v3_fixed_test_and_s63_reservation_v1"
)
S63_RESERVATION_SUMMARY_SCHEMA_VERSION = (
    "twirl_teacher_ssl_reserved_sector_tics_v1"
)
EXPECTED_SECTORS: tuple[int, ...] = tuple(range(56, 63))
EXPECTED_TIME_UNIT = "BJD - 2457000"
EXPECTED_FLUX_COLUMNS: tuple[str, str] = (
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP",
)
MANIFEST_FIELDS: frozenset[str] = frozenset(
    {
        "created_utc",
        "sector",
        "hlsp_root",
        "out_h5",
        "time_unit",
        "requested_columns",
        "n_discovered_files",
        "n_exported_targets",
        "skipped",
        "records",
    }
)
MANIFEST_SKIPPED_FIELDS: frozenset[str] = frozenset(
    {
        "read_failed",
        "tic_filter",
        "duplicate_tic",
        "no_flux_columns",
    }
)
MANIFEST_RECORD_FIELDS: frozenset[str] = frozenset(
    {
        "tic",
        "sector",
        "camera",
        "ccd",
        "tessmag",
        "n_cadences",
        "flux_columns",
        "source_fits",
    }
)
POOL_COLUMNS: tuple[str, ...] = (
    "pool_contract_version",
    "observation_id",
    "sector",
    "tic",
    "camera",
    "ccd",
    "tessmag",
    "n_cadences",
    "flux_columns_json",
    "compact_h5_path",
    "compact_h5_sha256",
    "compact_h5_size_bytes",
    "compact_group_path",
    "compact_manifest_path",
    "compact_manifest_sha256",
    "source_fits",
)
MAX_INT64 = np.iinfo(np.int64).max
ProgressCallback = Callable[[str], None]


@dataclass(frozen=True)
class _FileSnapshot:
    path: Path
    size_bytes: int
    mtime_ns: int
    device: int
    inode: int
    sha256: str

    def public_metadata(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }


def _stat_identity(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return (
        int(stat.st_size),
        int(stat.st_mtime_ns),
        int(stat.st_dev),
        int(stat.st_ino),
    )


def _snapshot_file(
    path: Path,
    *,
    progress: ProgressCallback | None = None,
    label: str | None = None,
) -> _FileSnapshot:
    resolved = Path(path).expanduser().resolve(strict=True)
    before = _stat_identity(resolved)
    digest = hashlib.sha256()
    bytes_read = 0
    last_update = time.monotonic()
    with resolved.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
            bytes_read += len(block)
            now = time.monotonic()
            if progress is not None and now - last_update >= 30.0:
                description = label or resolved.name
                progress(
                    f"[ssl-full-pool] hashing {description}: "
                    f"{bytes_read / (1024**3):.1f} GiB"
                )
                last_update = now
    after = _stat_identity(resolved)
    if before != after or bytes_read != after[0]:
        raise RuntimeError(f"input changed while it was hashed: {resolved}")
    return _FileSnapshot(
        path=resolved,
        size_bytes=after[0],
        mtime_ns=after[1],
        device=after[2],
        inode=after[3],
        sha256=digest.hexdigest(),
    )


def _verify_snapshot(snapshot: _FileSnapshot) -> None:
    try:
        observed = _stat_identity(snapshot.path)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"input disappeared before pool publication: {snapshot.path}"
        ) from exc
    expected = (
        snapshot.size_bytes,
        snapshot.mtime_ns,
        snapshot.device,
        snapshot.inode,
    )
    if observed != expected:
        raise RuntimeError(
            f"input changed before pool publication: {snapshot.path}"
        )


def _read_snapshot_bytes(snapshot: _FileSnapshot) -> bytes:
    payload = snapshot.path.read_bytes()
    if len(payload) != snapshot.size_bytes:
        raise RuntimeError(f"input changed while it was read: {snapshot.path}")
    if hashlib.sha256(payload).hexdigest() != snapshot.sha256:
        raise RuntimeError(f"input hash changed while it was read: {snapshot.path}")
    _verify_snapshot(snapshot)
    return payload


def _no_duplicate_json_keys(
    pairs: list[tuple[str, Any]],
) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for key, value in pairs:
        if key in output:
            raise ValueError(f"JSON object contains duplicate key {key!r}")
        output[key] = value
    return output


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON contains non-finite numeric constant {value!r}")


def _load_json_bytes(payload: bytes, *, context: str) -> Any:
    try:
        return json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=_no_duplicate_json_keys,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{context} is not strict UTF-8 JSON") from exc


def _require_fields(
    value: Any,
    expected: frozenset[str],
    *,
    context: str,
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{context} must be a JSON object")
    observed = set(value)
    if observed != set(expected):
        missing = sorted(set(expected) - observed)
        unexpected = sorted(observed - set(expected))
        raise ValueError(
            f"{context} has wrong fields; missing={missing}, "
            f"unexpected={unexpected}"
        )
    return value


def _strict_integer(
    value: Any,
    *,
    context: str,
    minimum: int = 0,
    maximum: int = MAX_INT64,
) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, np.integer)
    ):
        raise ValueError(f"{context} must be an integer")
    normalized = int(value)
    if normalized < minimum or normalized > maximum:
        raise ValueError(
            f"{context} must be in the closed interval "
            f"[{minimum}, {maximum}]"
        )
    return normalized


def _strict_float(value: Any, *, context: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{context} must be a finite number")
    normalized = float(value)
    if not np.isfinite(normalized):
        raise ValueError(f"{context} must be a finite number")
    return normalized


def _nonblank_text(value: Any, *, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} must be a nonblank string")
    if value != value.strip():
        raise ValueError(f"{context} must not have surrounding whitespace")
    return value


def _validate_flux_columns(value: Any, *, context: str) -> tuple[str, ...]:
    if not isinstance(value, list) or any(
        not isinstance(column, str) for column in value
    ):
        raise ValueError(f"{context} must be a JSON string array")
    normalized = tuple(value)
    if normalized != EXPECTED_FLUX_COLUMNS:
        raise ValueError(
            f"{context} must equal {list(EXPECTED_FLUX_COLUMNS)!r}; "
            f"observed={list(normalized)!r}"
        )
    return normalized


def _validate_export_manifest(
    value: Any,
    *,
    path: Path,
) -> tuple[int, list[dict[str, Any]], dict[str, Any]]:
    manifest = _require_fields(
        value,
        MANIFEST_FIELDS,
        context=f"compact export manifest {path}",
    )
    sector = _strict_integer(
        manifest["sector"],
        context=f"{path}: sector",
        minimum=1,
    )
    _nonblank_text(manifest["created_utc"], context=f"{path}: created_utc")
    _nonblank_text(manifest["hlsp_root"], context=f"{path}: hlsp_root")
    declared_out_h5 = _nonblank_text(
        manifest["out_h5"],
        context=f"{path}: out_h5",
    )
    if manifest["time_unit"] != EXPECTED_TIME_UNIT:
        raise ValueError(
            f"{path}: time_unit must equal {EXPECTED_TIME_UNIT!r}"
        )
    _validate_flux_columns(
        manifest["requested_columns"],
        context=f"{path}: requested_columns",
    )

    discovered = _strict_integer(
        manifest["n_discovered_files"],
        context=f"{path}: n_discovered_files",
        minimum=1,
    )
    exported = _strict_integer(
        manifest["n_exported_targets"],
        context=f"{path}: n_exported_targets",
        minimum=1,
    )
    skipped = _require_fields(
        manifest["skipped"],
        MANIFEST_SKIPPED_FIELDS,
        context=f"{path}: skipped",
    )
    normalized_skips = {
        key: _strict_integer(
            skipped[key],
            context=f"{path}: skipped.{key}",
        )
        for key in sorted(MANIFEST_SKIPPED_FIELDS)
    }
    if any(normalized_skips.values()):
        raise ValueError(
            f"{path}: compact export is not lossless; skipped={normalized_skips}"
        )

    records_value = manifest["records"]
    if not isinstance(records_value, list) or not records_value:
        raise ValueError(f"{path}: records must be a nonempty JSON array")
    if exported != len(records_value) or discovered != len(records_value):
        raise ValueError(
            f"{path}: manifest counts disagree with records; "
            f"discovered={discovered}, exported={exported}, "
            f"records={len(records_value)}"
        )

    records: list[dict[str, Any]] = []
    for index, raw_record in enumerate(records_value):
        context = f"{path}: records[{index}]"
        record = _require_fields(
            raw_record,
            MANIFEST_RECORD_FIELDS,
            context=context,
        )
        tic = _strict_integer(
            record["tic"],
            context=f"{context}.tic",
            minimum=1,
        )
        record_sector = _strict_integer(
            record["sector"],
            context=f"{context}.sector",
            minimum=1,
        )
        if record_sector != sector:
            raise ValueError(
                f"{context}.sector={record_sector} disagrees with "
                f"manifest sector={sector}"
            )
        camera = _strict_integer(
            record["camera"],
            context=f"{context}.camera",
            minimum=1,
            maximum=4,
        )
        ccd = _strict_integer(
            record["ccd"],
            context=f"{context}.ccd",
            minimum=1,
            maximum=4,
        )
        tessmag = _strict_float(
            record["tessmag"],
            context=f"{context}.tessmag",
        )
        n_cadences = _strict_integer(
            record["n_cadences"],
            context=f"{context}.n_cadences",
            minimum=1,
        )
        flux_columns = _validate_flux_columns(
            record["flux_columns"],
            context=f"{context}.flux_columns",
        )
        source_fits = _nonblank_text(
            record["source_fits"],
            context=f"{context}.source_fits",
        )
        records.append(
            {
                "tic": tic,
                "sector": record_sector,
                "camera": camera,
                "ccd": ccd,
                "tessmag": tessmag,
                "n_cadences": n_cadences,
                "flux_columns": flux_columns,
                "source_fits": source_fits,
            }
        )

    record_frame = pd.DataFrame.from_records(records)
    duplicate_tic = record_frame["tic"].duplicated(keep=False)
    if duplicate_tic.any():
        examples = sorted(
            record_frame.loc[duplicate_tic, "tic"].astype(int).unique()
        )[:5]
        raise ValueError(
            f"{path}: records contain duplicate sector/TIC observations; "
            f"first TICs={examples}"
        )
    duplicate_source = record_frame["source_fits"].duplicated(keep=False)
    if duplicate_source.any():
        examples = (
            record_frame.loc[duplicate_source, "source_fits"]
            .drop_duplicates()
            .head(5)
            .tolist()
        )
        raise ValueError(
            f"{path}: records reuse source_fits paths; first={examples}"
        )
    metadata = {
        "created_utc": str(manifest["created_utc"]),
        "declared_hlsp_root": str(manifest["hlsp_root"]),
        "declared_out_h5": declared_out_h5,
        "n_discovered_files": discovered,
        "n_exported_targets": exported,
    }
    return sector, records, metadata


def _decode_h5_text(value: Any, *, context: str) -> str:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError(f"{context} must be a scalar string")
    scalar = array.reshape(-1)[0]
    if isinstance(scalar, (bytes, np.bytes_)):
        try:
            scalar = bytes(scalar).decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"{context} must be UTF-8") from exc
    return _nonblank_text(scalar, context=context)


def _h5_integer_attribute(
    attrs: Any,
    name: str,
    *,
    path: Path,
    minimum: int,
) -> int:
    if name not in attrs:
        raise ValueError(f"{path}: missing HDF5 root attribute {name!r}")
    value = np.asarray(attrs[name])
    if value.size != 1:
        raise ValueError(f"{path}: HDF5 root attribute {name!r} is not scalar")
    scalar = value.reshape(-1)[0]
    return _strict_integer(
        scalar,
        context=f"{path}: HDF5 root attribute {name}",
        minimum=minimum,
    )


def _validate_compact_h5(
    *,
    path: Path,
    sector: int,
    records: Sequence[Mapping[str, Any]],
) -> None:
    import h5py

    with h5py.File(path, "r") as h5:
        observed_sector = _h5_integer_attribute(
            h5.attrs,
            "sector",
            path=path,
            minimum=1,
        )
        if observed_sector != sector:
            raise ValueError(
                f"{path}: HDF5 sector={observed_sector} disagrees with "
                f"manifest sector={sector}"
            )
        observed_count = _h5_integer_attribute(
            h5.attrs,
            "n_targets",
            path=path,
            minimum=1,
        )
        if observed_count != len(records):
            raise ValueError(
                f"{path}: HDF5 n_targets={observed_count} disagrees with "
                f"manifest records={len(records)}"
            )
        if "time_unit" not in h5.attrs:
            raise ValueError(f"{path}: missing HDF5 root attribute 'time_unit'")
        time_unit = _decode_h5_text(
            h5.attrs["time_unit"],
            context=f"{path}: HDF5 time_unit",
        )
        if time_unit != EXPECTED_TIME_UNIT:
            raise ValueError(
                f"{path}: HDF5 time_unit must equal {EXPECTED_TIME_UNIT!r}"
            )
        if "flux_columns" not in h5.attrs:
            raise ValueError(
                f"{path}: missing HDF5 root attribute 'flux_columns'"
            )
        flux_columns_json = _decode_h5_text(
            h5.attrs["flux_columns"],
            context=f"{path}: HDF5 flux_columns",
        )
        try:
            flux_columns = json.loads(
                flux_columns_json,
                parse_constant=_reject_json_constant,
            )
        except (json.JSONDecodeError, ValueError) as exc:
            raise ValueError(
                f"{path}: HDF5 flux_columns is not strict JSON"
            ) from exc
        _validate_flux_columns(
            flux_columns,
            context=f"{path}: HDF5 flux_columns",
        )
        if "targets" not in h5 or not isinstance(h5["targets"], h5py.Group):
            raise ValueError(f"{path}: missing HDF5 /targets group")
        observed_groups = sorted(str(value) for value in h5["targets"].keys())
        expected_groups = sorted(
            f"{int(record['tic']):016d}" for record in records
        )
        if observed_groups != expected_groups:
            missing = sorted(set(expected_groups) - set(observed_groups))[:5]
            unexpected = sorted(set(observed_groups) - set(expected_groups))[:5]
            raise ValueError(
                f"{path}: HDF5 /targets inventory disagrees with manifest; "
                f"missing={missing}, unexpected={unexpected}"
            )


def _read_split_registry(
    snapshot: _FileSnapshot,
) -> tuple[pd.DataFrame, set[int]]:
    if snapshot.path.suffix.lower() != ".csv":
        raise ValueError("Teacher-v3 TIC split registry must be a CSV file")
    payload = _read_snapshot_bytes(snapshot)
    try:
        registry = pd.read_csv(
            io.BytesIO(payload),
            dtype=str,
            keep_default_na=False,
        )
    except Exception as exc:
        raise ValueError(
            f"unable to read TIC split registry CSV: {snapshot.path}"
        ) from exc
    if tuple(registry.columns) != REGISTRY_COLUMNS:
        raise ValueError(
            "TIC split registry columns/order disagree with the frozen "
            f"schema; observed={registry.columns.tolist()}"
        )
    checked = validate_tic_split_registry(registry)
    fixed_test = set(
        checked.loc[checked["fixed_split"].eq("test"), "tic"]
        .astype(np.int64)
        .astype(int)
    )
    if not fixed_test:
        raise ValueError("TIC split registry contains no fixed-test TICs")
    return checked, fixed_test


def _tic_inventory_sha256(tics: Sequence[int] | set[int]) -> str:
    normalized = sorted(int(value) for value in tics)
    if (
        not normalized
        or any(value <= 0 or value > MAX_INT64 for value in normalized)
        or len(normalized) != len(set(normalized))
    ):
        raise ValueError("TIC inventory must contain unique positive int64 IDs")
    payload = "".join(f"{value}\n" for value in normalized).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _read_s63_reserved_tics(snapshot: _FileSnapshot) -> set[int]:
    if snapshot.path.suffix.lower() not in {".txt", ".list"}:
        raise ValueError(
            "S63 reserved TIC inventory must be one-TIC-per-line text"
        )
    payload = _read_snapshot_bytes(snapshot)
    if not payload.endswith(b"\n"):
        raise ValueError("S63 reserved TIC inventory must end with a newline")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ValueError("S63 reserved TIC inventory must be ASCII") from exc
    lines = text.splitlines()
    if not lines or any(not line for line in lines):
        raise ValueError(
            "S63 reserved TIC inventory must contain one TIC on every line"
        )
    if any(not line.isdigit() for line in lines):
        raise ValueError(
            "S63 reserved TIC inventory contains a non-decimal TIC value"
        )
    tics = [int(line) for line in lines]
    if any(
        value <= 0
        or value > MAX_INT64
        or str(value) != line
        for value, line in zip(tics, lines, strict=True)
    ):
        raise ValueError(
            "S63 reserved TIC inventory must use canonical positive int64 IDs"
        )
    if tics != sorted(set(tics)):
        raise ValueError(
            "S63 reserved TIC inventory must be strictly increasing and unique"
        )
    return set(tics)


def _validate_s63_summary(
    snapshot: _FileSnapshot,
    *,
    reservation: _FileSnapshot,
    reserved_tics: set[int],
) -> dict[str, Any]:
    raw = _load_json_bytes(
        _read_snapshot_bytes(snapshot),
        context=f"S63 reservation summary {snapshot.path}",
    )
    if not isinstance(raw, Mapping):
        raise ValueError("S63 reservation summary must be a JSON object")
    required = {
        "schema_version",
        "identity_only",
        "light_curves_opened",
        "sector",
        "orbits",
        "n_reserved_tics",
        "reserved_tics_sha256",
    }
    missing = sorted(required - set(raw))
    if missing:
        raise ValueError(
            f"S63 reservation summary is missing fields: {missing}"
        )
    if raw["schema_version"] != S63_RESERVATION_SUMMARY_SCHEMA_VERSION:
        raise ValueError("S63 reservation summary schema_version mismatch")
    if raw["identity_only"] is not True or raw["light_curves_opened"] is not False:
        raise ValueError("S63 reservation summary is not identity-only")
    if _strict_integer(raw["sector"], context="S63 summary sector") != 63:
        raise ValueError("S63 reservation summary must declare sector 63")
    if raw["orbits"] != [133, 134]:
        raise ValueError("S63 reservation summary must declare orbits [133, 134]")
    count = _strict_integer(
        raw["n_reserved_tics"],
        context="S63 summary n_reserved_tics",
        minimum=1,
    )
    if count != len(reserved_tics):
        raise ValueError("S63 reservation summary TIC count mismatch")
    if str(raw["reserved_tics_sha256"]).lower() != reservation.sha256:
        raise ValueError("S63 reservation summary SHA-256 mismatch")
    return dict(raw)


def _identity_sha256(rows: pd.DataFrame) -> str:
    identities = (
        rows.loc[:, ["sector", "tic"]]
        .sort_values(["sector", "tic"], kind="stable")
        .to_dict(orient="records")
    )
    digest = hashlib.sha256()
    for identity in identities:
        digest.update(
            json.dumps(
                {
                    "sector": int(identity["sector"]),
                    "tic": int(identity["tic"]),
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _table_payload(table: pd.DataFrame, *, parquet: bool) -> bytes:
    if not parquet:
        return table.to_csv(index=False, lineterminator="\n").encode("utf-8")
    buffer = io.BytesIO()
    table.to_parquet(
        buffer,
        engine="pyarrow",
        compression="zstd",
        index=False,
    )
    return buffer.getvalue()


def _publish_immutable(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                "refusing to replace immutable output with different bytes: "
                f"{path}"
            )
        return
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            if path.read_bytes() != payload:
                raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _artifact_metadata(path: Path, payload: bytes, **extra: Any) -> dict[str, Any]:
    return {
        "path": str(path),
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        **extra,
    }


def _resolve_default_h5(manifest_path: Path) -> Path:
    suffix = ".manifest.json"
    if not manifest_path.name.endswith(suffix):
        raise ValueError(
            "compact export manifest filename must end in '.manifest.json': "
            f"{manifest_path}"
        )
    return manifest_path.with_name(manifest_path.name[: -len(suffix)] + ".h5")


def _sector_counts(
    rows: pd.DataFrame,
    *,
    fixed_test_tics: set[int],
    s63_reserved_tics: set[int],
) -> dict[str, Any]:
    tics = rows["tic"].astype(np.int64)
    fixed = tics.isin(fixed_test_tics)
    s63 = tics.isin(s63_reserved_tics)
    excluded = fixed | s63
    retained = ~excluded

    def summarize(mask: pd.Series) -> dict[str, int]:
        selected = rows.loc[mask]
        return {
            "n_observations": int(len(selected)),
            "n_unique_tics": int(selected["tic"].nunique()),
        }

    return {
        "input": summarize(pd.Series(True, index=rows.index)),
        "excluded": summarize(excluded),
        "retained": summarize(retained),
        "exclusion_categories": {
            "fixed_test_only": summarize(fixed & ~s63),
            "s63_reserved_only": summarize(s63 & ~fixed),
            "fixed_test_and_s63_reserved": summarize(fixed & s63),
        },
    }


def write_ssl_full_pool(
    *,
    compact_manifest_paths: Sequence[Path],
    split_registry_path: Path,
    s63_reserved_tics_path: Path,
    out_dir: Path,
    compact_h5_by_sector: Mapping[int, Path] | None = None,
    s63_reserved_summary_path: Path | None = None,
    progress: ProgressCallback | None = None,
) -> dict[str, Any]:
    """Validate, filter, and immutably publish the S56--S62 SSL pool."""

    if len(compact_manifest_paths) != len(EXPECTED_SECTORS):
        raise ValueError(
            f"exactly {len(EXPECTED_SECTORS)} compact export manifests are "
            f"required for sectors {EXPECTED_SECTORS}"
        )
    manifest_paths = [Path(value).expanduser() for value in compact_manifest_paths]
    resolved_manifests = [path.resolve(strict=True) for path in manifest_paths]
    if len(set(resolved_manifests)) != len(resolved_manifests):
        raise ValueError("compact export manifest paths must be unique")

    overrides = {
        int(sector): Path(path).expanduser()
        for sector, path in (compact_h5_by_sector or {}).items()
    }
    invalid_override_sectors = sorted(set(overrides) - set(EXPECTED_SECTORS))
    if invalid_override_sectors:
        raise ValueError(
            f"compact HDF5 overrides contain invalid sectors: "
            f"{invalid_override_sectors}"
        )

    manifest_inputs: dict[int, dict[str, Any]] = {}
    snapshots: list[_FileSnapshot] = []
    for path in resolved_manifests:
        snapshot = _snapshot_file(path)
        snapshots.append(snapshot)
        raw = _load_json_bytes(
            _read_snapshot_bytes(snapshot),
            context=f"compact export manifest {path}",
        )
        sector, records, metadata = _validate_export_manifest(raw, path=path)
        if sector in manifest_inputs:
            raise ValueError(
                f"duplicate compact export manifest for sector {sector}"
            )
        manifest_inputs[sector] = {
            "snapshot": snapshot,
            "records": records,
            "metadata": metadata,
        }
    observed_sectors = tuple(sorted(manifest_inputs))
    if observed_sectors != EXPECTED_SECTORS:
        missing = sorted(set(EXPECTED_SECTORS) - set(observed_sectors))
        unexpected = sorted(set(observed_sectors) - set(EXPECTED_SECTORS))
        raise ValueError(
            "compact export sector coverage must be exactly S56--S62; "
            f"missing={missing}, unexpected={unexpected}"
        )
    unused_overrides = sorted(set(overrides) - set(manifest_inputs))
    if unused_overrides:
        raise ValueError(f"unused compact HDF5 overrides: {unused_overrides}")

    for sector in EXPECTED_SECTORS:
        manifest_snapshot = manifest_inputs[sector]["snapshot"]
        h5_path = overrides.get(
            sector,
            _resolve_default_h5(manifest_snapshot.path),
        )
        if progress is not None:
            progress(
                f"[ssl-full-pool] hashing S{sector} compact HDF5: "
                f"{Path(h5_path).expanduser()}"
            )
        h5_snapshot = _snapshot_file(
            h5_path,
            progress=progress,
            label=f"S{sector} compact HDF5",
        )
        snapshots.append(h5_snapshot)
        _validate_compact_h5(
            path=h5_snapshot.path,
            sector=sector,
            records=manifest_inputs[sector]["records"],
        )
        _verify_snapshot(h5_snapshot)
        manifest_inputs[sector]["h5_snapshot"] = h5_snapshot
        if progress is not None:
            progress(
                f"[ssl-full-pool] validated S{sector}: "
                f"{len(manifest_inputs[sector]['records']):,} targets, "
                f"sha256={h5_snapshot.sha256}"
            )

    registry_snapshot = _snapshot_file(split_registry_path)
    snapshots.append(registry_snapshot)
    registry, fixed_test_tics = _read_split_registry(registry_snapshot)

    s63_snapshot = _snapshot_file(s63_reserved_tics_path)
    snapshots.append(s63_snapshot)
    s63_reserved_tics = _read_s63_reserved_tics(s63_snapshot)

    reservation_summary_snapshot: _FileSnapshot | None = None
    reservation_summary: dict[str, Any] | None = None
    if s63_reserved_summary_path is None:
        adjacent = s63_snapshot.path.with_suffix(".summary.json")
        if adjacent.is_file():
            s63_reserved_summary_path = adjacent
    if s63_reserved_summary_path is not None:
        reservation_summary_snapshot = _snapshot_file(
            s63_reserved_summary_path
        )
        snapshots.append(reservation_summary_snapshot)
        reservation_summary = _validate_s63_summary(
            reservation_summary_snapshot,
            reservation=s63_snapshot,
            reserved_tics=s63_reserved_tics,
        )

    full_records: list[dict[str, Any]] = []
    input_export_metadata: list[dict[str, Any]] = []
    for sector in EXPECTED_SECTORS:
        item = manifest_inputs[sector]
        manifest_snapshot = item["snapshot"]
        h5_snapshot = item["h5_snapshot"]
        records = item["records"]
        for record in records:
            full_records.append(
                {
                    **record,
                    "compact_h5_path": str(h5_snapshot.path),
                    "compact_h5_sha256": h5_snapshot.sha256,
                    "compact_h5_size_bytes": h5_snapshot.size_bytes,
                    "compact_manifest_path": str(manifest_snapshot.path),
                    "compact_manifest_sha256": manifest_snapshot.sha256,
                }
            )
        sector_identity = pd.DataFrame(
            {
                "sector": [sector] * len(records),
                "tic": [int(record["tic"]) for record in records],
            }
        )
        input_export_metadata.append(
            {
                "sector": sector,
                "manifest": manifest_snapshot.public_metadata(),
                "compact_h5": h5_snapshot.public_metadata(),
                **item["metadata"],
                "observation_identity_sha256": _identity_sha256(
                    sector_identity
                ),
            }
        )
    full = pd.DataFrame.from_records(full_records)
    full = full.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )
    duplicate_observation = full.duplicated(["sector", "tic"], keep=False)
    if duplicate_observation.any():
        examples = (
            full.loc[duplicate_observation, ["sector", "tic"]]
            .head(5)
            .to_dict(orient="records")
        )
        raise ValueError(
            "compact manifests contain duplicate (sector, tic) observations "
            f"across inputs; first={examples}"
        )

    fixed_mask = full["tic"].isin(fixed_test_tics)
    s63_mask = full["tic"].isin(s63_reserved_tics)
    excluded_mask = fixed_mask | s63_mask
    retained = full.loc[~excluded_mask].copy()
    if retained.empty:
        raise ValueError("whole-host exclusions leave an empty SSL pool")
    retained_sectors = tuple(sorted(retained["sector"].astype(int).unique()))
    if retained_sectors != EXPECTED_SECTORS:
        raise ValueError(
            "whole-host exclusions leave one or more sectors empty; "
            f"retained={retained_sectors}"
        )
    if retained["tic"].isin(fixed_test_tics | s63_reserved_tics).any():
        raise AssertionError("whole-host exclusion leakage survived filtering")

    retained["pool_contract_version"] = FULL_POOL_CONTRACT_VERSION
    retained["observation_id"] = [
        f"s{int(sector):04d}-tic{int(tic):016d}"
        for sector, tic in zip(
            retained["sector"],
            retained["tic"],
            strict=True,
        )
    ]
    retained["flux_columns_json"] = json.dumps(
        list(EXPECTED_FLUX_COLUMNS),
        separators=(",", ":"),
    )
    retained["compact_group_path"] = [
        f"targets/{int(tic):016d}" for tic in retained["tic"]
    ]
    pool = retained.loc[:, list(POOL_COLUMNS)].copy()
    for column in (
        "sector",
        "tic",
        "camera",
        "ccd",
        "n_cadences",
        "compact_h5_size_bytes",
    ):
        pool[column] = pool[column].astype(np.int64)
    pool["tessmag"] = pool["tessmag"].astype(np.float64)
    pool = pool.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )
    if pool["observation_id"].duplicated().any():
        raise AssertionError("constructed duplicate observation IDs")

    out_dir = Path(out_dir).expanduser().resolve()
    csv_path = out_dir / "teacher_ssl_full_pool_observations.csv"
    parquet_path = out_dir / "teacher_ssl_full_pool_observations.parquet"
    summary_path = out_dir / "teacher_ssl_full_pool_manifest.summary.json"
    allowlist_paths = {
        sector: out_dir / "allowlists" / f"s{sector}_tics.csv"
        for sector in EXPECTED_SECTORS
    }
    output_paths = {
        csv_path,
        parquet_path,
        summary_path,
        *allowlist_paths.values(),
    }
    input_paths = {snapshot.path for snapshot in snapshots}
    if output_paths & input_paths:
        raise ValueError("input and output paths must be distinct")

    csv_payload = _table_payload(pool, parquet=False)
    parquet_payload = _table_payload(pool, parquet=True)
    artifact_payloads: dict[Path, bytes] = {
        csv_path: csv_payload,
        parquet_path: parquet_payload,
    }
    allowlist_metadata: dict[str, dict[str, Any]] = {}
    for sector in EXPECTED_SECTORS:
        sector_tics = (
            pool.loc[pool["sector"].eq(sector), "tic"]
            .astype(np.int64)
            .astype(int)
            .tolist()
        )
        if sector_tics != sorted(set(sector_tics)) or not sector_tics:
            raise AssertionError(f"S{sector} allowlist is not sorted and unique")
        allowlist_payload = (
            "tic\n" + "".join(f"{tic}\n" for tic in sector_tics)
        ).encode("ascii")
        path = allowlist_paths[sector]
        artifact_payloads[path] = allowlist_payload
        allowlist_metadata[str(sector)] = _artifact_metadata(
            path,
            allowlist_payload,
            n_tics=len(sector_tics),
            tic_inventory_sha256=_tic_inventory_sha256(sector_tics),
        )

    fixed_only = fixed_mask & ~s63_mask
    s63_only = s63_mask & ~fixed_mask
    both = fixed_mask & s63_mask
    all_input_tics = set(full["tic"].astype(np.int64).astype(int))
    retained_tics = set(pool["tic"].astype(np.int64).astype(int))
    per_sector = {
        str(sector): _sector_counts(
            full.loc[full["sector"].eq(sector)],
            fixed_test_tics=fixed_test_tics,
            s63_reserved_tics=s63_reserved_tics,
        )
        for sector in EXPECTED_SECTORS
    }
    summary: dict[str, Any] = {
        "passed": True,
        "schema_version": FULL_POOL_SUMMARY_SCHEMA_VERSION,
        "pool_contract_version": FULL_POOL_CONTRACT_VERSION,
        "selection_policy_version": FULL_POOL_SELECTION_POLICY_VERSION,
        "sectors": list(EXPECTED_SECTORS),
        "observation_identity_columns": ["sector", "tic"],
        "exclusion_identity_column": "tic",
        "exclusion_scope": (
            "whole host: exclude every S56--S62 observation of any TIC in "
            "the frozen fixed-test registry or prospective S63 reservation"
        ),
        "required_flux_columns": list(EXPECTED_FLUX_COLUMNS),
        "required_time_unit": EXPECTED_TIME_UNIT,
        "inputs": {
            "compact_exports": input_export_metadata,
            "tic_split_registry": {
                **registry_snapshot.public_metadata(),
                "n_registry_tics": int(len(registry)),
                "n_fixed_test_tics": len(fixed_test_tics),
                "fixed_test_tic_inventory_sha256": _tic_inventory_sha256(
                    fixed_test_tics
                ),
            },
            "s63_reserved_tics": {
                **s63_snapshot.public_metadata(),
                "n_reserved_tics": len(s63_reserved_tics),
                "tic_inventory_sha256": _tic_inventory_sha256(
                    s63_reserved_tics
                ),
                "summary": (
                    reservation_summary_snapshot.public_metadata()
                    if reservation_summary_snapshot is not None
                    else None
                ),
                "summary_validated": reservation_summary is not None,
            },
        },
        "counts": {
            "input": {
                "n_observations": int(len(full)),
                "n_unique_tics": len(all_input_tics),
            },
            "declared_exclusion_tics": {
                "n_fixed_test": len(fixed_test_tics),
                "n_s63_reserved": len(s63_reserved_tics),
                "n_intersection": len(
                    fixed_test_tics & s63_reserved_tics
                ),
                "n_union": len(fixed_test_tics | s63_reserved_tics),
            },
            "excluded_input": {
                "n_observations": int(excluded_mask.sum()),
                "n_unique_tics": len(
                    all_input_tics
                    & (fixed_test_tics | s63_reserved_tics)
                ),
                "categories": {
                    "fixed_test_only": {
                        "n_observations": int(fixed_only.sum()),
                        "n_unique_tics": int(
                            full.loc[fixed_only, "tic"].nunique()
                        ),
                    },
                    "s63_reserved_only": {
                        "n_observations": int(s63_only.sum()),
                        "n_unique_tics": int(
                            full.loc[s63_only, "tic"].nunique()
                        ),
                    },
                    "fixed_test_and_s63_reserved": {
                        "n_observations": int(both.sum()),
                        "n_unique_tics": int(
                            full.loc[both, "tic"].nunique()
                        ),
                    },
                },
            },
            "retained": {
                "n_observations": int(len(pool)),
                "n_unique_tics": len(retained_tics),
                "n_multisector_tics": int(
                    pool.groupby("tic")["sector"].nunique().gt(1).sum()
                ),
            },
        },
        "per_sector": per_sector,
        "identity_hashes": {
            "input_observations_sha256": _identity_sha256(full),
            "retained_observations_sha256": _identity_sha256(pool),
        },
        "leakage_audit": {
            "fixed_test_observations_retained": int(
                pool["tic"].isin(fixed_test_tics).sum()
            ),
            "s63_reserved_observations_retained": int(
                pool["tic"].isin(s63_reserved_tics).sum()
            ),
        },
        "outputs": {
            "csv": _artifact_metadata(
                csv_path,
                csv_payload,
                n_rows=int(len(pool)),
            ),
            "parquet": _artifact_metadata(
                parquet_path,
                parquet_payload,
                n_rows=int(len(pool)),
            ),
            "sector_allowlists": allowlist_metadata,
        },
    }
    summary_payload = (
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    artifact_payloads[summary_path] = summary_payload

    for path, payload in artifact_payloads.items():
        if path.exists() and path.read_bytes() != payload:
            raise FileExistsError(
                "refusing to replace immutable output with different bytes: "
                f"{path}"
            )
    for snapshot in snapshots:
        _verify_snapshot(snapshot)
    for path, payload in artifact_payloads.items():
        _publish_immutable(path, payload)
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        expected = hashlib.sha256(payload).hexdigest()
        if observed != expected:
            raise RuntimeError(f"published artifact hash mismatch: {path}")
    return summary


__all__ = [
    "EXPECTED_FLUX_COLUMNS",
    "EXPECTED_SECTORS",
    "FULL_POOL_CONTRACT_VERSION",
    "FULL_POOL_SELECTION_POLICY_VERSION",
    "FULL_POOL_SUMMARY_SCHEMA_VERSION",
    "POOL_COLUMNS",
    "write_ssl_full_pool",
]
