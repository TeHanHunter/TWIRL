"""Global checksum boundary for the S56--S62 full-pool BLS products."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.search.a2v1_bls_contract import (
    approved_a2v1_teacher_bls_config,
    bls_config_sha256,
)
from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
)
from twirl.vetting.ssl_full_pool import (
    EXPECTED_SECTORS,
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SELECTION_POLICY_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
    POOL_COLUMNS,
)


GLOBAL_BLS_CONTRACT_VERSION = "twirl_teacher_ssl_full_pool_bls_global_v1"
GLOBAL_BLS_SUMMARY_SCHEMA_VERSION = (
    "twirl_teacher_ssl_full_pool_bls_global_summary_v1"
)
FULL_POOL_BLS_SOURCE_PRODUCT_TAG = "A2v1-fullpool-v1"
TARGET_SELECTION_CONTRACT_VERSION = "a2v1_bls_target_allowlist_v1"
ORBITID_RECONCILIATION_CONTRACT_VERSION = (
    "a2v1_compact_orbitid_reconciliation_v1"
)
ORBITID_POLICY_BY_SECTOR: dict[int, str] = {
    **{sector: "strict" for sector in range(56, 62)},
    62: "reference_by_cadence",
}
MAX_INT64 = np.iinfo(np.int64).max
_BLS_SUMMARY_REQUIRED_FIELDS: frozenset[str] = frozenset(
    {
        "passed",
        "sector",
        "contract_version",
        "bls_search_contract_version",
        "bls_config_sha256",
        "compact_lc",
        "compact_lc_sha256",
        "cadence_reference",
        "cadence_reference_sha256",
        "cadence_reference_manifest",
        "cadence_reference_manifest_sha256",
        "cadence_reference_contract_version",
        "cadence_reference_cadence_authority",
        "cadence_reference_quality_authority",
        "cadence_reference_source_hashes_sha256",
        "authority_exclusion_policy_contract",
        "authority_exclusion_external_bit",
        "authority_exclusions_sha256",
        "n_authority_exclusions",
        "target_selection_contract_version",
        "target_allowlist",
        "target_allowlist_sha256",
        "target_allowlist_count",
        "target_allowlist_tics_sha256",
        "orbitid_policy",
        "orbitid_reconciliation_contract_version",
        "n_cad_orbitid_reference_matched",
        "n_cad_orbitid_mismatch",
        "n_cad_orbitid_corrected",
        "n_targets_orbitid_mismatch",
        "orbitid_corrections_sha256",
        "apertures",
        "n_targets",
        "n_targets_total",
        "n_unique_tics",
        "n_rows",
        "n_periods",
        "n_peaks",
        "n_shards",
        "shard_index",
        "n_source_shards",
        "source_product_tag",
        "config",
        "peak_table_sha256",
        "outputs",
    }
)
_BLS_ROW_REQUIRED_COLUMNS: frozenset[str] = frozenset(
    {
        "sector",
        "tic",
        "aperture",
        "peak_rank",
        "status",
        "bls_search_branch",
        "adp_only_contract_version",
        "source_product_tag",
        "bls_n_periods",
        "bls_n_peaks",
        "bls_p_min_d",
        "bls_p_max_cap_d",
        "bls_max_period_fraction",
        "bls_sigma_clip",
        "bls_orbit_edge_trim_d",
        "bls_search_contract_version",
        "bls_config_sha256",
        "orbitid_policy",
        "orbitid_reconciliation_contract_version",
        "n_cad_orbitid_reference_matched",
        "n_cad_orbitid_mismatch",
        "n_cad_orbitid_corrected",
        "orbitid_correction_signature_sha256",
    }
)


@dataclass(frozen=True)
class _FileSnapshot:
    path: Path
    size_bytes: int
    mtime_ns: int
    device: int
    inode: int
    sha256: str

    def metadata(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }


@dataclass(frozen=True)
class _PoolAuthority:
    summary: Mapping[str, Any]
    summary_snapshot: _FileSnapshot
    identities: pd.DataFrame
    compact_exports: Mapping[int, Mapping[str, Any]]
    allowlist_metadata: Mapping[int, Mapping[str, Any]]
    allowlist_tics: Mapping[int, tuple[int, ...]]
    artifact_snapshots: tuple[_FileSnapshot, ...]


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
    cache: dict[Path, _FileSnapshot],
) -> _FileSnapshot:
    resolved = Path(path).expanduser().resolve(strict=True)
    if resolved in cache:
        return cache[resolved]
    before = _stat_identity(resolved)
    digest = hashlib.sha256()
    bytes_read = 0
    with resolved.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
            bytes_read += len(block)
    after = _stat_identity(resolved)
    if before != after or bytes_read != after[0]:
        raise RuntimeError(f"input changed while it was hashed: {resolved}")
    snapshot = _FileSnapshot(
        path=resolved,
        size_bytes=after[0],
        mtime_ns=after[1],
        device=after[2],
        inode=after[3],
        sha256=digest.hexdigest(),
    )
    cache[resolved] = snapshot
    return snapshot


def _verify_snapshot(snapshot: _FileSnapshot) -> None:
    expected = (
        snapshot.size_bytes,
        snapshot.mtime_ns,
        snapshot.device,
        snapshot.inode,
    )
    try:
        observed = _stat_identity(snapshot.path)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"input disappeared before global BLS publication: {snapshot.path}"
        ) from exc
    if observed != expected:
        raise RuntimeError(
            f"input changed before global BLS publication: {snapshot.path}"
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
    raise ValueError(f"JSON contains non-finite constant {value!r}")


def _load_json(snapshot: _FileSnapshot, *, context: str) -> Mapping[str, Any]:
    try:
        value = json.loads(
            _read_snapshot_bytes(snapshot).decode("utf-8"),
            object_pairs_hook=_no_duplicate_json_keys,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{context} is not strict UTF-8 JSON") from exc
    if not isinstance(value, Mapping):
        raise ValueError(f"{context} must be a JSON object")
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
            f"{context} must be in [{minimum}, {maximum}]"
        )
    return normalized


def _digest(value: Any, *, context: str) -> str:
    normalized = str(value).lower()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise ValueError(f"{context} must be a lowercase SHA-256 digest")
    return normalized


def _canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _tic_inventory_sha256(tics: Sequence[int] | set[int]) -> str:
    normalized = sorted(int(value) for value in tics)
    if (
        not normalized
        or any(value <= 0 or value > MAX_INT64 for value in normalized)
        or len(normalized) != len(set(normalized))
    ):
        raise ValueError("TIC inventory must contain unique positive int64 IDs")
    payload = "".join(f"{tic}\n" for tic in normalized).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _normalize_identity_frame(
    frame: pd.DataFrame,
    *,
    context: str,
) -> pd.DataFrame:
    missing = sorted({"sector", "tic"} - set(frame.columns))
    if missing:
        raise ValueError(f"{context} is missing identity columns: {missing}")
    output = frame.loc[:, ["sector", "tic"]].copy()
    for column in ("sector", "tic"):
        numeric = pd.to_numeric(output[column], errors="coerce")
        array = numeric.to_numpy(dtype=float)
        invalid = (
            ~np.isfinite(array)
            | (array <= 0)
            | (array > MAX_INT64)
            | ~np.equal(array, np.rint(array))
        )
        if invalid.any():
            raise ValueError(f"{context} contains invalid {column} values")
        output[column] = np.rint(array).astype(np.int64)
    duplicate = output.duplicated(["sector", "tic"], keep=False)
    if duplicate.any():
        examples = output.loc[duplicate].head(5).to_dict(orient="records")
        raise ValueError(
            f"{context} contains duplicate (sector, tic) observations; "
            f"first={examples}"
        )
    return output.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )


def _identity_sha256(identities: pd.DataFrame) -> str:
    normalized = _normalize_identity_frame(
        identities,
        context="observation identity table",
    )
    digest = hashlib.sha256()
    for row in normalized.itertuples(index=False):
        digest.update(
            json.dumps(
                {"sector": int(row.sector), "tic": int(row.tic)},
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _validated_artifact(
    metadata: Any,
    *,
    context: str,
    cache: dict[Path, _FileSnapshot],
) -> _FileSnapshot:
    if not isinstance(metadata, Mapping):
        raise ValueError(f"{context} metadata must be an object")
    missing = sorted({"path", "size_bytes", "sha256"} - set(metadata))
    if missing:
        raise ValueError(f"{context} metadata is missing fields: {missing}")
    snapshot = _snapshot_file(Path(str(metadata["path"])), cache=cache)
    declared_size = _strict_integer(
        metadata["size_bytes"],
        context=f"{context} size_bytes",
        minimum=1,
    )
    if snapshot.size_bytes != declared_size:
        raise ValueError(f"{context} size disagrees with its frozen summary")
    if snapshot.sha256 != _digest(
        metadata["sha256"],
        context=f"{context} sha256",
    ):
        raise ValueError(f"{context} SHA-256 disagrees with its frozen summary")
    return snapshot


def _read_canonical_allowlist(
    snapshot: _FileSnapshot,
) -> tuple[int, ...]:
    if snapshot.path.suffix.lower() != ".csv":
        raise ValueError("frozen sector allowlists must use CSV")
    try:
        frame = pd.read_csv(
            snapshot.path,
            dtype=str,
            keep_default_na=False,
        )
    except Exception as exc:
        raise ValueError(
            f"unable to read frozen sector allowlist: {snapshot.path}"
        ) from exc
    if frame.columns.tolist() != ["tic"] or frame.empty:
        raise ValueError(
            "frozen sector allowlist must have exactly one nonempty TIC column"
        )
    values = frame["tic"].tolist()
    if any(
        not value
        or not value.isascii()
        or not value.isdigit()
        or str(int(value)) != value
        for value in values
    ):
        raise ValueError("frozen sector allowlist contains noncanonical TICs")
    tics = tuple(int(value) for value in values)
    if (
        any(tic <= 0 or tic > MAX_INT64 for tic in tics)
        or list(tics) != sorted(set(tics))
    ):
        raise ValueError(
            "frozen sector allowlist must be sorted, unique, positive int64"
        )
    return tics


def _validate_frozen_pool(
    *,
    summary_path: Path,
    cache: dict[Path, _FileSnapshot],
) -> _PoolAuthority:
    summary_snapshot = _snapshot_file(summary_path, cache=cache)
    summary = _load_json(
        summary_snapshot,
        context="frozen full-pool summary",
    )
    if summary.get("passed") is not True:
        raise ValueError("frozen full-pool summary did not pass")
    if summary.get("schema_version") != FULL_POOL_SUMMARY_SCHEMA_VERSION:
        raise ValueError("frozen full-pool summary schema mismatch")
    if summary.get("pool_contract_version") != FULL_POOL_CONTRACT_VERSION:
        raise ValueError("frozen full-pool contract mismatch")
    if (
        summary.get("selection_policy_version")
        != FULL_POOL_SELECTION_POLICY_VERSION
    ):
        raise ValueError("frozen full-pool selection policy mismatch")
    if summary.get("sectors") != list(EXPECTED_SECTORS):
        raise ValueError("frozen full-pool sector coverage is not S56--S62")
    leakage = summary.get("leakage_audit")
    if not isinstance(leakage, Mapping) or any(
        _strict_integer(value, context=f"pool leakage {key}") != 0
        for key, value in leakage.items()
    ):
        raise ValueError("frozen full-pool summary reports host leakage")

    outputs = summary.get("outputs")
    if not isinstance(outputs, Mapping):
        raise ValueError("frozen full-pool summary lacks outputs")
    csv_snapshot = _validated_artifact(
        outputs.get("csv"),
        context="frozen pool CSV",
        cache=cache,
    )
    parquet_snapshot = _validated_artifact(
        outputs.get("parquet"),
        context="frozen pool Parquet",
        cache=cache,
    )
    pool = pd.read_parquet(parquet_snapshot.path)
    if tuple(pool.columns) != POOL_COLUMNS:
        raise ValueError(
            "frozen pool Parquet columns disagree with the pool contract"
        )
    identities = _normalize_identity_frame(
        pool,
        context="frozen pool Parquet",
    )
    csv_identities = _normalize_identity_frame(
        pd.read_csv(csv_snapshot.path, usecols=["sector", "tic"]),
        context="frozen pool CSV",
    )
    try:
        pd.testing.assert_frame_equal(identities, csv_identities)
    except AssertionError as exc:
        raise ValueError(
            "frozen pool CSV and Parquet observation identities disagree"
        ) from exc
    if not identities["sector"].isin(EXPECTED_SECTORS).all():
        raise ValueError("frozen pool contains sectors outside S56--S62")
    observed_sectors = tuple(sorted(identities["sector"].astype(int).unique()))
    if observed_sectors != EXPECTED_SECTORS:
        raise ValueError("frozen pool does not contain every S56--S62 sector")

    retained_counts = summary.get("counts", {}).get("retained", {})
    expected_rows = _strict_integer(
        retained_counts.get("n_observations"),
        context="frozen pool retained n_observations",
        minimum=1,
    )
    expected_tics = _strict_integer(
        retained_counts.get("n_unique_tics"),
        context="frozen pool retained n_unique_tics",
        minimum=1,
    )
    if expected_rows != len(identities):
        raise ValueError("frozen pool retained observation count mismatch")
    if expected_tics != identities["tic"].nunique():
        raise ValueError("frozen pool retained TIC count mismatch")
    for name in ("csv", "parquet"):
        declared_rows = _strict_integer(
            outputs[name].get("n_rows"),
            context=f"frozen pool {name} n_rows",
            minimum=1,
        )
        if declared_rows != len(identities):
            raise ValueError(f"frozen pool {name} row count mismatch")
    declared_identity_hash = (
        summary.get("identity_hashes", {}).get(
            "retained_observations_sha256"
        )
    )
    if _digest(
        declared_identity_hash,
        context="frozen pool retained identity hash",
    ) != _identity_sha256(identities):
        raise ValueError("frozen pool retained identity hash mismatch")

    inputs = summary.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ValueError("frozen full-pool summary lacks inputs")
    exports_value = inputs.get("compact_exports")
    if not isinstance(exports_value, list) or len(exports_value) != len(
        EXPECTED_SECTORS
    ):
        raise ValueError("frozen pool must declare seven compact exports")
    compact_exports: dict[int, Mapping[str, Any]] = {}
    for item in exports_value:
        if not isinstance(item, Mapping):
            raise ValueError("frozen pool compact export must be an object")
        sector = _strict_integer(
            item.get("sector"),
            context="frozen pool compact-export sector",
            minimum=1,
        )
        if sector in compact_exports:
            raise ValueError(
                f"frozen pool repeats compact export for sector {sector}"
            )
        compact_h5 = item.get("compact_h5")
        if not isinstance(compact_h5, Mapping):
            raise ValueError(f"frozen S{sector} compact HDF5 metadata is absent")
        missing = sorted(
            {"path", "size_bytes", "sha256"} - set(compact_h5)
        )
        if missing:
            raise ValueError(
                f"frozen S{sector} compact HDF5 metadata lacks {missing}"
            )
        compact_path = Path(str(compact_h5["path"])).expanduser().resolve(
            strict=True
        )
        if compact_path.stat().st_size != _strict_integer(
            compact_h5["size_bytes"],
            context=f"frozen S{sector} compact HDF5 size",
            minimum=1,
        ):
            raise ValueError(
                f"frozen S{sector} compact HDF5 size has changed"
            )
        _digest(
            compact_h5["sha256"],
            context=f"frozen S{sector} compact HDF5 sha256",
        )
        compact_exports[sector] = item
    if tuple(sorted(compact_exports)) != EXPECTED_SECTORS:
        raise ValueError("frozen pool compact-export sectors are not S56--S62")

    allowlists_value = outputs.get("sector_allowlists")
    if not isinstance(allowlists_value, Mapping) or set(allowlists_value) != {
        str(sector) for sector in EXPECTED_SECTORS
    }:
        raise ValueError("frozen pool sector-allowlist coverage is not S56--S62")
    allowlist_metadata: dict[int, Mapping[str, Any]] = {}
    allowlist_tics: dict[int, tuple[int, ...]] = {}
    artifact_snapshots: list[_FileSnapshot] = [
        csv_snapshot,
        parquet_snapshot,
    ]
    for sector in EXPECTED_SECTORS:
        metadata = allowlists_value[str(sector)]
        snapshot = _validated_artifact(
            metadata,
            context=f"frozen S{sector} allowlist",
            cache=cache,
        )
        artifact_snapshots.append(snapshot)
        tics = _read_canonical_allowlist(snapshot)
        declared_count = _strict_integer(
            metadata.get("n_tics"),
            context=f"frozen S{sector} allowlist n_tics",
            minimum=1,
        )
        if declared_count != len(tics):
            raise ValueError(f"frozen S{sector} allowlist count mismatch")
        semantic_hash = _tic_inventory_sha256(tics)
        if semantic_hash != _digest(
            metadata.get("tic_inventory_sha256"),
            context=f"frozen S{sector} allowlist TIC hash",
        ):
            raise ValueError(
                f"frozen S{sector} allowlist TIC-inventory hash mismatch"
            )
        expected_sector_tics = tuple(
            identities.loc[identities["sector"].eq(sector), "tic"]
            .astype(np.int64)
            .astype(int)
            .tolist()
        )
        if tics != expected_sector_tics:
            raise ValueError(
                f"frozen S{sector} allowlist does not exactly cover pool rows"
            )
        per_sector = summary.get("per_sector", {}).get(str(sector), {})
        retained = per_sector.get("retained", {})
        if _strict_integer(
            retained.get("n_observations"),
            context=f"frozen S{sector} retained observation count",
            minimum=1,
        ) != len(tics):
            raise ValueError(f"frozen S{sector} retained count mismatch")
        allowlist_metadata[sector] = metadata
        allowlist_tics[sector] = tics
    return _PoolAuthority(
        summary=summary,
        summary_snapshot=summary_snapshot,
        identities=identities,
        compact_exports=compact_exports,
        allowlist_metadata=allowlist_metadata,
        allowlist_tics=allowlist_tics,
        artifact_snapshots=tuple(artifact_snapshots),
    )


def _expected_bls_config() -> dict[str, Any]:
    config = approved_a2v1_teacher_bls_config()
    config["source_product_tag"] = FULL_POOL_BLS_SOURCE_PRODUCT_TAG
    return config


def _same_path(left: Any, right: Any) -> bool:
    try:
        return Path(str(left)).expanduser().resolve() == Path(
            str(right)
        ).expanduser().resolve()
    except (OSError, RuntimeError):
        return False


def _validate_orbit_rows(
    frame: pd.DataFrame,
    *,
    policy: str,
    summary: Mapping[str, Any],
    sector: int,
) -> dict[str, Any]:
    count_columns = (
        "n_cad_orbitid_reference_matched",
        "n_cad_orbitid_mismatch",
        "n_cad_orbitid_corrected",
    )
    if set(frame["orbitid_policy"].astype(str)) != {policy}:
        raise ValueError(f"S{sector} BLS rows have the wrong orbit-ID policy")
    if set(
        frame["orbitid_reconciliation_contract_version"].astype(str)
    ) != {ORBITID_RECONCILIATION_CONTRACT_VERSION}:
        raise ValueError(
            f"S{sector} BLS rows have the wrong orbit-ID contract"
        )
    consistency_columns = [
        *count_columns,
        "orbitid_correction_signature_sha256",
    ]
    if (
        frame.groupby("tic", sort=False)[consistency_columns]
        .nunique(dropna=False)
        .gt(1)
        .any()
        .any()
    ):
        raise ValueError(
            f"S{sector} orbit-ID fields disagree within a TIC"
        )
    targets = frame.drop_duplicates("tic", keep="first").copy()
    for column in count_columns:
        values = pd.to_numeric(targets[column], errors="coerce")
        if (
            values.isna().any()
            or (values < 0).any()
            or (values != np.floor(values)).any()
        ):
            raise ValueError(f"S{sector} BLS rows have invalid {column}")
        targets[column] = values.astype(np.int64)
    if (
        targets["n_cad_orbitid_mismatch"]
        > targets["n_cad_orbitid_reference_matched"]
    ).any():
        raise ValueError(
            f"S{sector} orbit-ID mismatches exceed matched cadences"
        )
    if policy == "strict":
        if (
            targets[
                ["n_cad_orbitid_mismatch", "n_cad_orbitid_corrected"]
            ].to_numpy(dtype=np.int64)
            != 0
        ).any():
            raise ValueError(
                f"S{sector} strict orbit-ID policy contains corrections"
            )
    elif policy == "reference_by_cadence":
        if not targets["n_cad_orbitid_corrected"].eq(
            targets["n_cad_orbitid_mismatch"]
        ).all():
            raise ValueError(
                f"S{sector} reference_by_cadence did not correct all mismatches"
            )
    else:
        raise ValueError(f"unsupported orbit-ID policy {policy!r}")
    signatures = targets["orbitid_correction_signature_sha256"].astype(str)
    if not signatures.str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError(
            f"S{sector} BLS rows contain invalid orbit-ID signatures"
        )
    records = [
        {
            "tic": int(row.tic),
            "n_cad_orbitid_reference_matched": int(
                row.n_cad_orbitid_reference_matched
            ),
            "n_cad_orbitid_mismatch": int(row.n_cad_orbitid_mismatch),
            "n_cad_orbitid_corrected": int(row.n_cad_orbitid_corrected),
            "orbitid_correction_signature_sha256": str(
                row.orbitid_correction_signature_sha256
            ),
        }
        for row in targets.sort_values("tic", kind="stable").itertuples(
            index=False
        )
    ]
    observed = {
        "n_cad_orbitid_reference_matched": int(
            targets["n_cad_orbitid_reference_matched"].sum()
        ),
        "n_cad_orbitid_mismatch": int(
            targets["n_cad_orbitid_mismatch"].sum()
        ),
        "n_cad_orbitid_corrected": int(
            targets["n_cad_orbitid_corrected"].sum()
        ),
        "n_targets_orbitid_mismatch": int(
            targets["n_cad_orbitid_mismatch"].gt(0).sum()
        ),
        "orbitid_corrections_sha256": _canonical_json_sha256(
            {
                "contract_version": (
                    ORBITID_RECONCILIATION_CONTRACT_VERSION
                ),
                "policy": policy,
                "targets": records,
            }
        ),
    }
    for field, value in observed.items():
        if summary.get(field) != value:
            raise ValueError(
                f"S{sector} BLS orbit-ID summary mismatch for {field}"
            )
    return observed


def _validate_constant_column(
    frame: pd.DataFrame,
    column: str,
    expected: Any,
    *,
    sector: int,
) -> None:
    values = frame[column]
    if isinstance(expected, float):
        numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
        if (
            not np.isfinite(numeric).all()
            or not np.equal(numeric, expected).all()
        ):
            raise ValueError(
                f"S{sector} BLS row column {column} disagrees with config"
            )
    elif set(values.astype(str)) != {str(expected)}:
        raise ValueError(
            f"S{sector} BLS row column {column} disagrees with config"
        )


def _validate_sector_bls(
    *,
    sector: int,
    table_path: Path,
    summary_path: Path,
    pool: _PoolAuthority,
    expected_columns: tuple[str, ...] | None,
    expected_dtypes: tuple[str, ...] | None,
    cache: dict[Path, _FileSnapshot],
) -> tuple[
    pd.DataFrame,
    Mapping[str, Any],
    _FileSnapshot,
    _FileSnapshot,
    dict[str, Any],
]:
    table_snapshot = _snapshot_file(table_path, cache=cache)
    summary_snapshot = _snapshot_file(summary_path, cache=cache)
    summary = _load_json(
        summary_snapshot,
        context=f"S{sector} merged BLS summary",
    )
    missing_fields = sorted(_BLS_SUMMARY_REQUIRED_FIELDS - set(summary))
    if missing_fields:
        raise ValueError(
            f"S{sector} merged BLS summary is missing fields: {missing_fields}"
        )
    if summary["passed"] is not True:
        raise ValueError(f"S{sector} merged BLS summary did not pass")
    if _strict_integer(
        summary["sector"],
        context=f"S{sector} BLS summary sector",
        minimum=1,
    ) != sector:
        raise ValueError(f"S{sector} BLS summary declares the wrong sector")
    if summary["contract_version"] != ADP_ONLY_CONTRACT_VERSION:
        raise ValueError(f"S{sector} BLS ADP-only contract mismatch")
    if _strict_integer(
        summary["n_shards"],
        context=f"S{sector} BLS n_shards",
        minimum=1,
    ) != 1 or _strict_integer(
        summary["shard_index"],
        context=f"S{sector} BLS shard_index",
    ) != 0:
        raise ValueError(f"S{sector} input is not a merged BLS product")
    source_shards = _strict_integer(
        summary["n_source_shards"],
        context=f"S{sector} BLS n_source_shards",
        minimum=1,
    )
    if source_shards < 1:
        raise ValueError(f"S{sector} BLS has no source shards")
    if _digest(
        summary["peak_table_sha256"],
        context=f"S{sector} BLS peak-table sha256",
    ) != table_snapshot.sha256:
        raise ValueError(f"S{sector} BLS peak-table SHA-256 mismatch")
    outputs = summary["outputs"]
    if not isinstance(outputs, Mapping):
        raise ValueError(f"S{sector} BLS outputs must be an object")
    if not _same_path(outputs.get("peak_table"), table_snapshot.path):
        raise ValueError(f"S{sector} BLS output peak-table path mismatch")
    if not _same_path(outputs.get("summary"), summary_snapshot.path):
        raise ValueError(f"S{sector} BLS output summary path mismatch")

    config = _expected_bls_config()
    expected_config_hash = bls_config_sha256(config)
    if summary["config"] != config:
        raise ValueError(f"S{sector} BLS config is not the locked full-pool config")
    if summary["bls_config_sha256"] != expected_config_hash:
        raise ValueError(f"S{sector} BLS config SHA-256 mismatch")
    if summary["source_product_tag"] != FULL_POOL_BLS_SOURCE_PRODUCT_TAG:
        raise ValueError(f"S{sector} BLS source-product tag mismatch")
    if _strict_integer(
        summary["n_periods"],
        context=f"S{sector} BLS n_periods",
        minimum=1,
    ) != int(config["n_periods"]):
        raise ValueError(f"S{sector} BLS n_periods mismatch")
    if _strict_integer(
        summary["n_peaks"],
        context=f"S{sector} BLS n_peaks",
        minimum=1,
    ) != int(config["n_peaks"]):
        raise ValueError(f"S{sector} BLS n_peaks mismatch")
    if summary["apertures"] != list(ADP_ONLY_APERTURES):
        raise ValueError(f"S{sector} BLS aperture contract mismatch")
    expected_search_contract = (
        "custom"
        if config != approved_a2v1_teacher_bls_config()
        else str(summary["bls_search_contract_version"])
    )
    if summary["bls_search_contract_version"] != expected_search_contract:
        raise ValueError(f"S{sector} BLS search-contract declaration mismatch")

    compact = pool.compact_exports[sector]["compact_h5"]
    if not _same_path(summary["compact_lc"], compact["path"]):
        raise ValueError(
            f"S{sector} BLS compact-HDF5 path disagrees with frozen pool"
        )
    if _digest(
        summary["compact_lc_sha256"],
        context=f"S{sector} BLS compact sha256",
    ) != _digest(
        compact["sha256"],
        context=f"frozen S{sector} compact sha256",
    ):
        raise ValueError(
            f"S{sector} BLS compact-HDF5 SHA-256 disagrees with frozen pool"
        )
    allowlist = pool.allowlist_metadata[sector]
    if (
        summary["target_selection_contract_version"]
        != TARGET_SELECTION_CONTRACT_VERSION
    ):
        raise ValueError(f"S{sector} BLS target-selection contract mismatch")
    if not _same_path(summary["target_allowlist"], allowlist["path"]):
        raise ValueError(
            f"S{sector} BLS allowlist path disagrees with frozen pool"
        )
    if _digest(
        summary["target_allowlist_sha256"],
        context=f"S{sector} BLS allowlist sha256",
    ) != _digest(
        allowlist["sha256"],
        context=f"frozen S{sector} allowlist sha256",
    ):
        raise ValueError(
            f"S{sector} BLS allowlist SHA-256 disagrees with frozen pool"
        )
    expected_tics = pool.allowlist_tics[sector]
    if _strict_integer(
        summary["target_allowlist_count"],
        context=f"S{sector} BLS target_allowlist_count",
        minimum=1,
    ) != len(expected_tics):
        raise ValueError(
            f"S{sector} BLS allowlist count disagrees with frozen pool"
        )
    if _digest(
        summary["target_allowlist_tics_sha256"],
        context=f"S{sector} BLS allowlist TIC hash",
    ) != _digest(
        allowlist["tic_inventory_sha256"],
        context=f"frozen S{sector} allowlist TIC hash",
    ):
        raise ValueError(
            f"S{sector} BLS allowlist TIC hash disagrees with frozen pool"
        )

    policy = ORBITID_POLICY_BY_SECTOR[sector]
    if summary["orbitid_policy"] != policy:
        raise ValueError(
            f"S{sector} BLS orbit-ID policy must equal {policy!r}"
        )
    if (
        summary["orbitid_reconciliation_contract_version"]
        != ORBITID_RECONCILIATION_CONTRACT_VERSION
    ):
        raise ValueError(f"S{sector} BLS orbit-ID contract mismatch")
    if summary["cadence_reference_contract_version"] != (
        f"s{sector}_a2v1_cadence_reference_v1"
    ):
        raise ValueError(f"S{sector} BLS cadence-reference contract mismatch")
    for field in (
        "cadence_reference_sha256",
        "cadence_reference_manifest_sha256",
        "cadence_reference_source_hashes_sha256",
        "authority_exclusions_sha256",
    ):
        _digest(summary[field], context=f"S{sector} BLS {field}")

    frame = pd.read_parquet(table_snapshot.path)
    columns = tuple(str(column) for column in frame.columns)
    dtypes = tuple(str(dtype) for dtype in frame.dtypes)
    if expected_columns is not None and columns != expected_columns:
        raise ValueError("sector BLS tables have different column schemas")
    if expected_dtypes is not None and dtypes != expected_dtypes:
        raise ValueError("sector BLS tables have different column dtypes")
    missing_columns = sorted(_BLS_ROW_REQUIRED_COLUMNS - set(frame.columns))
    if missing_columns:
        raise ValueError(
            f"S{sector} BLS table lacks required columns: {missing_columns}"
        )
    if _strict_integer(
        summary["n_rows"],
        context=f"S{sector} BLS n_rows",
        minimum=1,
    ) != len(frame):
        raise ValueError(f"S{sector} BLS table row-count mismatch")
    sector_values = pd.to_numeric(frame["sector"], errors="coerce")
    if (
        sector_values.isna().any()
        or not sector_values.eq(sector).all()
    ):
        raise ValueError(f"S{sector} BLS table contains wrong-sector rows")
    tic_values = pd.to_numeric(frame["tic"], errors="coerce")
    tic_array = tic_values.to_numpy(dtype=float)
    if (
        (~np.isfinite(tic_array)).any()
        or (tic_array <= 0).any()
        or (tic_array > MAX_INT64).any()
        or (~np.equal(tic_array, np.rint(tic_array))).any()
    ):
        raise ValueError(f"S{sector} BLS table contains invalid TICs")
    frame["sector"] = np.full(len(frame), sector, dtype=np.int64)
    frame["tic"] = np.rint(tic_array).astype(np.int64)
    observed_tics = tuple(sorted(frame["tic"].astype(int).unique()))
    if observed_tics != expected_tics:
        missing = sorted(set(expected_tics) - set(observed_tics))[:5]
        unexpected = sorted(set(observed_tics) - set(expected_tics))[:5]
        raise ValueError(
            f"S{sector} BLS coverage disagrees with frozen pool; "
            f"missing={missing}, unexpected={unexpected}"
        )
    for field in ("n_targets", "n_targets_total", "n_unique_tics"):
        if _strict_integer(
            summary[field],
            context=f"S{sector} BLS {field}",
            minimum=1,
        ) != len(expected_tics):
            raise ValueError(
                f"S{sector} BLS {field} disagrees with frozen pool"
            )

    apertures = frame["aperture"].fillna("").astype(str)
    if set(apertures) != set(ADP_ONLY_APERTURES):
        raise ValueError(f"S{sector} BLS rows do not use the exact ADP pair")
    aperture_coverage = (
        frame.assign(_aperture=apertures)
        .groupby("tic", sort=False)["_aperture"]
        .nunique()
    )
    if not aperture_coverage.eq(len(ADP_ONLY_APERTURES)).all():
        raise ValueError(
            f"S{sector} BLS targets do not all cover both ADP apertures"
        )
    ranks = pd.to_numeric(frame["peak_rank"], errors="coerce")
    if (
        ranks.isna().any()
        or (ranks < 0).any()
        or (ranks != np.floor(ranks)).any()
    ):
        raise ValueError(f"S{sector} BLS table has invalid peak ranks")
    frame["peak_rank"] = ranks.astype(np.int64)
    if frame.duplicated(["tic", "aperture", "peak_rank"]).any():
        raise ValueError(
            f"S{sector} BLS table has duplicate TIC/aperture/peak-rank rows"
        )

    row_constants = {
        "bls_search_branch": "current_adp",
        "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
        "source_product_tag": FULL_POOL_BLS_SOURCE_PRODUCT_TAG,
        "bls_n_periods": int(config["n_periods"]),
        "bls_n_peaks": int(config["n_peaks"]),
        "bls_p_min_d": float(config["p_min_d"]),
        "bls_p_max_cap_d": float(config["p_max_cap_d"]),
        "bls_max_period_fraction": float(config["max_period_fraction"]),
        "bls_sigma_clip": float(config["sigma_clip"]),
        "bls_orbit_edge_trim_d": float(config["orbit_edge_trim_d"]),
        "bls_search_contract_version": expected_search_contract,
        "bls_config_sha256": expected_config_hash,
    }
    for column, expected in row_constants.items():
        _validate_constant_column(
            frame,
            column,
            expected,
            sector=sector,
        )
    orbit_summary = _validate_orbit_rows(
        frame,
        policy=policy,
        summary=summary,
        sector=sector,
    )
    _verify_snapshot(table_snapshot)
    _verify_snapshot(summary_snapshot)
    frame = frame.sort_values(
        ["sector", "tic", "aperture", "peak_rank"],
        kind="stable",
    ).reset_index(drop=True)
    product = {
        "sector": sector,
        "peak_table": table_snapshot.metadata(),
        "summary": summary_snapshot.metadata(),
        "compact_h5": {
            "path": str(
                Path(str(compact["path"])).expanduser().resolve()
            ),
            "size_bytes": int(compact["size_bytes"]),
            "sha256": str(compact["sha256"]),
        },
        "allowlist": {
            "path": str(
                Path(str(allowlist["path"])).expanduser().resolve()
            ),
            "size_bytes": int(allowlist["size_bytes"]),
            "sha256": str(allowlist["sha256"]),
            "n_tics": len(expected_tics),
            "tic_inventory_sha256": str(
                allowlist["tic_inventory_sha256"]
            ),
        },
        "orbitid_policy": policy,
        "n_source_shards": source_shards,
        "n_observations": len(expected_tics),
        "n_rows": int(len(frame)),
        **orbit_summary,
    }
    return (
        frame,
        summary,
        table_snapshot,
        summary_snapshot,
        product,
    )


def _preflight_immutable_file(
    path: Path,
    *,
    expected_size: int,
    expected_sha256: str,
) -> None:
    if not path.exists():
        return
    stat = path.stat()
    if stat.st_size != expected_size:
        raise FileExistsError(
            f"refusing to replace immutable output with different bytes: {path}"
        )
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != expected_sha256:
        raise FileExistsError(
            f"refusing to replace immutable output with different bytes: {path}"
        )


def _publish_temp_immutable(
    temporary: Path,
    output: Path,
    *,
    expected_size: int,
    expected_sha256: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    _preflight_immutable_file(
        output,
        expected_size=expected_size,
        expected_sha256=expected_sha256,
    )
    if output.exists():
        return
    try:
        os.link(temporary, output)
    except FileExistsError:
        _preflight_immutable_file(
            output,
            expected_size=expected_size,
            expected_sha256=expected_sha256,
        )


def _publish_bytes_immutable(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    expected_sha256 = hashlib.sha256(payload).hexdigest()
    _preflight_immutable_file(
        path,
        expected_size=len(payload),
        expected_sha256=expected_sha256,
    )
    if path.exists():
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
            _preflight_immutable_file(
                path,
                expected_size=len(payload),
                expected_sha256=expected_sha256,
            )
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_global_full_pool_bls(
    *,
    pool_summary_path: Path,
    sector_bls_paths: Mapping[int, Path],
    out_parquet_path: Path,
    sector_summary_paths: Mapping[int, Path] | None = None,
    out_summary_path: Path | None = None,
) -> dict[str, Any]:
    """Validate and immutably consolidate seven merged sector BLS products."""

    normalized_bls = {
        int(sector): Path(path).expanduser()
        for sector, path in sector_bls_paths.items()
    }
    if tuple(sorted(normalized_bls)) != EXPECTED_SECTORS:
        missing = sorted(set(EXPECTED_SECTORS) - set(normalized_bls))
        unexpected = sorted(set(normalized_bls) - set(EXPECTED_SECTORS))
        raise ValueError(
            "sector BLS products must cover exactly S56--S62; "
            f"missing={missing}, unexpected={unexpected}"
        )
    normalized_summaries = {
        int(sector): Path(path).expanduser()
        for sector, path in (sector_summary_paths or {}).items()
    }
    invalid_summary_sectors = sorted(
        set(normalized_summaries) - set(EXPECTED_SECTORS)
    )
    if invalid_summary_sectors:
        raise ValueError(
            "sector BLS summary overrides contain invalid sectors: "
            f"{invalid_summary_sectors}"
        )

    out_parquet_path = Path(out_parquet_path).expanduser().resolve()
    if out_parquet_path.suffix.lower() not in {".parquet", ".pq"}:
        raise ValueError("global BLS output must use Parquet")
    out_summary_path = (
        out_parquet_path.with_suffix(".summary.json")
        if out_summary_path is None
        else Path(out_summary_path).expanduser().resolve()
    )
    if out_summary_path.suffix.lower() != ".json":
        raise ValueError("global BLS summary output must use JSON")
    if out_summary_path == out_parquet_path:
        raise ValueError("global BLS Parquet and summary paths must differ")

    cache: dict[Path, _FileSnapshot] = {}
    pool = _validate_frozen_pool(
        summary_path=pool_summary_path,
        cache=cache,
    )
    input_paths = {
        pool.summary_snapshot.path,
        *(snapshot.path for snapshot in pool.artifact_snapshots),
    }
    out_parquet_path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    writer = None
    expected_columns: tuple[str, ...] | None = None
    expected_dtypes: tuple[str, ...] | None = None
    arrow_schema = None
    sector_products: list[dict[str, Any]] = []
    observed_identities: list[pd.DataFrame] = []
    n_rows = 0
    global_tics: set[int] = set()
    status_counts: Counter[str] = Counter()
    aperture_counts: Counter[str] = Counter()
    search_contract: str | None = None
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq

        with tempfile.NamedTemporaryFile(
            dir=out_parquet_path.parent,
            prefix=f".{out_parquet_path.name}.",
            suffix=".tmp.parquet",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
        for sector in EXPECTED_SECTORS:
            table_path = normalized_bls[sector]
            summary_path = normalized_summaries.get(
                sector,
                table_path.with_suffix(".summary.json"),
            )
            (
                frame,
                summary,
                table_snapshot,
                summary_snapshot,
                product,
            ) = _validate_sector_bls(
                sector=sector,
                table_path=table_path,
                summary_path=summary_path,
                pool=pool,
                expected_columns=expected_columns,
                expected_dtypes=expected_dtypes,
                cache=cache,
            )
            input_paths.update(
                {table_snapshot.path, summary_snapshot.path}
            )
            columns = tuple(str(column) for column in frame.columns)
            dtypes = tuple(str(dtype) for dtype in frame.dtypes)
            if expected_columns is None:
                expected_columns = columns
                expected_dtypes = dtypes
            current_search_contract = str(
                summary["bls_search_contract_version"]
            )
            if search_contract is None:
                search_contract = current_search_contract
            elif search_contract != current_search_contract:
                raise ValueError(
                    "sector BLS summaries disagree on search contract"
                )
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                arrow_schema = table.schema
                writer = pq.ParquetWriter(
                    temporary,
                    arrow_schema,
                    compression="zstd",
                )
            elif table.schema != arrow_schema:
                raise ValueError(
                    "sector BLS tables have different Arrow schemas"
                )
            writer.write_table(table)

            identities = frame.loc[:, ["sector", "tic"]].drop_duplicates()
            observed_identities.append(identities)
            n_rows += int(len(frame))
            global_tics.update(frame["tic"].astype(int))
            status_counts.update(
                frame["status"].fillna("").astype(str).tolist()
            )
            aperture_counts.update(
                frame["aperture"].fillna("").astype(str).tolist()
            )
            sector_products.append(product)
            del frame, table
        if writer is None:
            raise ValueError("global BLS consolidation received no rows")
        writer.close()
        writer = None

        global_identities = _normalize_identity_frame(
            pd.concat(observed_identities, ignore_index=True),
            context="global BLS observations",
        )
        try:
            pd.testing.assert_frame_equal(
                global_identities,
                pool.identities,
            )
        except AssertionError as exc:
            raise ValueError(
                "global BLS (sector, tic) coverage does not exactly equal "
                "the frozen full pool"
            ) from exc
        identity_hash = _identity_sha256(global_identities)
        declared_pool_hash = _digest(
            pool.summary["identity_hashes"][
                "retained_observations_sha256"
            ],
            context="frozen pool identity hash",
        )
        if identity_hash != declared_pool_hash:
            raise ValueError("global BLS identity hash disagrees with frozen pool")

        for snapshot in cache.values():
            _verify_snapshot(snapshot)
        output_snapshot_cache: dict[Path, _FileSnapshot] = {}
        temporary_snapshot = _snapshot_file(
            temporary,
            cache=output_snapshot_cache,
        )
        if temporary_snapshot.size_bytes <= 0:
            raise RuntimeError("global BLS temporary Parquet is empty")
        expected_config = _expected_bls_config()
        summary: dict[str, Any] = {
            "passed": True,
            "schema_version": GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
            "contract_version": GLOBAL_BLS_CONTRACT_VERSION,
            "sectors": list(EXPECTED_SECTORS),
            "observation_identity_columns": ["sector", "tic"],
            "frozen_pool": {
                "summary": pool.summary_snapshot.metadata(),
                "contract_version": FULL_POOL_CONTRACT_VERSION,
                "selection_policy_version": (
                    FULL_POOL_SELECTION_POLICY_VERSION
                ),
                "observation_identity_sha256": identity_hash,
                "n_observations": int(len(pool.identities)),
                "n_unique_tics": int(pool.identities["tic"].nunique()),
            },
            "bls_contract": {
                "config": expected_config,
                "config_sha256": bls_config_sha256(expected_config),
                "search_contract_version": search_contract,
                "source_product_tag": FULL_POOL_BLS_SOURCE_PRODUCT_TAG,
                "apertures": list(ADP_ONLY_APERTURES),
                "orbitid_reconciliation_contract_version": (
                    ORBITID_RECONCILIATION_CONTRACT_VERSION
                ),
                "orbitid_policy_by_sector": {
                    str(sector): ORBITID_POLICY_BY_SECTOR[sector]
                    for sector in EXPECTED_SECTORS
                },
            },
            "counts": {
                "n_rows": n_rows,
                "n_observations": int(len(global_identities)),
                "n_unique_tics": len(global_tics),
                "n_multisector_tics": int(
                    global_identities.groupby("tic")["sector"]
                    .nunique()
                    .gt(1)
                    .sum()
                ),
                "status_rows": {
                    key: int(value)
                    for key, value in sorted(status_counts.items())
                },
                "aperture_rows": {
                    key: int(value)
                    for key, value in sorted(aperture_counts.items())
                },
            },
            "coverage_audit": {
                "missing_frozen_observations": 0,
                "unexpected_bls_observations": 0,
                "observation_identity_sha256": identity_hash,
            },
            "sector_products": sector_products,
            "output": {
                "path": str(out_parquet_path),
                "size_bytes": temporary_snapshot.size_bytes,
                "sha256": temporary_snapshot.sha256,
                "n_rows": n_rows,
                "n_observations": int(len(global_identities)),
            },
        }
        summary_payload = (
            json.dumps(
                summary,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
        if out_parquet_path in input_paths or out_summary_path in input_paths:
            raise ValueError("global BLS input and output paths must be distinct")
        _preflight_immutable_file(
            out_parquet_path,
            expected_size=temporary_snapshot.size_bytes,
            expected_sha256=temporary_snapshot.sha256,
        )
        _preflight_immutable_file(
            out_summary_path,
            expected_size=len(summary_payload),
            expected_sha256=hashlib.sha256(summary_payload).hexdigest(),
        )
        for snapshot in cache.values():
            _verify_snapshot(snapshot)
        _publish_temp_immutable(
            temporary,
            out_parquet_path,
            expected_size=temporary_snapshot.size_bytes,
            expected_sha256=temporary_snapshot.sha256,
        )
        _publish_bytes_immutable(out_summary_path, summary_payload)
        return summary
    finally:
        if writer is not None:
            writer.close()
        if temporary is not None:
            temporary.unlink(missing_ok=True)


__all__ = [
    "FULL_POOL_BLS_SOURCE_PRODUCT_TAG",
    "GLOBAL_BLS_CONTRACT_VERSION",
    "GLOBAL_BLS_SUMMARY_SCHEMA_VERSION",
    "ORBITID_POLICY_BY_SECTOR",
    "write_global_full_pool_bls",
]
