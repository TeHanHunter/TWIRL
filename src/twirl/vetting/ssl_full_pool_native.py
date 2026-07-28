"""Native raw/error preparation for the label-free S56--S62 SSL pool.

This module deliberately does not call ``build_raw_pair_export``.  That
function retains the immutable Teacher-v3 S62 release binding (997 rows and
exact release hashes).  The broad SSL pool needs the same model-facing
channels, but a different observation-keyed, real-only contract whose S62
orbit-ID reconciliation is validated per cadence rather than by a frozen
row-count signature.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.lightcurves.a2v1_cadence_reference import (
    S56_EXPECTED_DETECTORS,
    S56_EXPECTED_ORBITS,
)
from twirl.lightcurves.external_quality import (
    EXPECTED_CADENCE_AUTHORITY,
    ORBITID_POLICY_REFERENCE,
    ORBITID_POLICY_STRICT,
    ORBITID_RECONCILIATION_CONTRACT_VERSION,
    load_external_quality_reference,
)
from twirl.vetting import harmonic_export
from twirl.vetting.harmonic_inputs import (
    CHRONOLOGY_SMALL_CHANNELS,
    CHRONOLOGY_SUPPLEMENTAL_CHANNELS,
    HARMONIC_VIEW_CHANNELS,
    NATIVE_DATASETS,
    ORBITID_RECONCILIATION_MASK_DATASET,
    PERIODOGRAM_CHANNELS,
    PERIODOGRAM_DATASETS,
    RAW_PAIR_EXTERNAL_QUALITY_ATTRS,
    RAW_PAIR_ORBITID_RECONCILIATION_ATTRS,
    RAW_PAIR_QUALITY_COUNT_NAMES,
    read_native_light_curve_from_h5,
)
from twirl.vetting.ssl_full_pool import (
    EXPECTED_SECTORS,
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
)
from twirl.vetting.teacher_native_registry import (
    file_sha256,
    read_table,
    write_native_input_registry,
)


FULL_POOL_RAW_SOURCE_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_tglc_raw_source_v1"
)
FULL_POOL_NATIVE_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_real_native_v1"
)
FULL_POOL_NATIVE_RELEASE_BINDING = (
    "teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp_real_only"
)
FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_native_shard_summary_v1"
)
FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_native_registry_source_v1"
)
FULL_POOL_NATIVE_SECTORS = EXPECTED_SECTORS
S62_CORRECTION_CADENCE_MIN = 766_048
S62_CORRECTION_CADENCE_MAX = 766_136
S62_REFERENCE_ORBIT_131_MIN = 760_742
S62_REFERENCE_ORBIT_131_MAX = 766_136
S62_REFERENCE_ORBIT_132_MIN = 766_228
S62_REFERENCE_ORBIT_132_MAX = 771_851
S62_CORRECTION_INPUT_ORBIT = 132
S62_CORRECTION_REFERENCE_ORBIT = 131


@dataclass(frozen=True)
class FileBinding:
    """Immutable file identity captured before a long-running build."""

    path: Path
    size_bytes: int
    mtime_ns: int
    device: int
    inode: int
    sha256: str

    def assert_unchanged(self) -> None:
        stat = self.path.stat()
        observed = (
            int(stat.st_size),
            int(stat.st_mtime_ns),
            int(stat.st_dev),
            int(stat.st_ino),
        )
        expected = (
            self.size_bytes,
            self.mtime_ns,
            self.device,
            self.inode,
        )
        if observed != expected:
            raise RuntimeError(f"input changed while in use: {self.path}")


@dataclass(frozen=True)
class SectorPoolAuthority:
    """One sector's exact frozen-pool rows and transitive input bindings."""

    sector: int
    rows: pd.DataFrame
    pool: FileBinding
    pool_summary: FileBinding
    allowlist: FileBinding
    compact_h5_sha256: str
    compact_h5_size_bytes: int
    observation_identity_sha256: str


def _bind_file(
    path: Path,
    *,
    expected_sha256: str | None = None,
    expected_size_bytes: int | None = None,
) -> FileBinding:
    resolved = Path(path).expanduser().resolve(strict=True)
    before = resolved.stat()
    digest = file_sha256(resolved)
    after = resolved.stat()
    before_identity = (
        int(before.st_size),
        int(before.st_mtime_ns),
        int(before.st_dev),
        int(before.st_ino),
    )
    after_identity = (
        int(after.st_size),
        int(after.st_mtime_ns),
        int(after.st_dev),
        int(after.st_ino),
    )
    if before_identity != after_identity:
        raise RuntimeError(f"input changed while it was hashed: {resolved}")
    if expected_sha256 is not None and digest != str(expected_sha256):
        raise ValueError(
            f"SHA-256 mismatch for {resolved}: {digest} != {expected_sha256}"
        )
    if (
        expected_size_bytes is not None
        and int(after.st_size) != int(expected_size_bytes)
    ):
        raise ValueError(
            f"size mismatch for {resolved}: {after.st_size} != "
            f"{expected_size_bytes}"
        )
    return FileBinding(
        path=resolved,
        size_bytes=int(after.st_size),
        mtime_ns=int(after.st_mtime_ns),
        device=int(after.st_dev),
        inode=int(after.st_ino),
        sha256=digest,
    )


def _observation_identity_sha256(rows: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for row in (
        rows.loc[:, ["sector", "tic"]]
        .sort_values(["sector", "tic"], kind="stable")
        .itertuples(index=False)
    ):
        digest.update(
            json.dumps(
                {"sector": int(row.sector), "tic": int(row.tic)},
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _tic_inventory_sha256(tics: Sequence[int]) -> str:
    digest = hashlib.sha256()
    for tic in sorted(int(value) for value in tics):
        digest.update(f"{tic}\n".encode("ascii"))
    return digest.hexdigest()


def _positive_tics(values: pd.Series, *, context: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    invalid = (
        ~np.isfinite(numeric)
        | (numeric <= 0)
        | ~np.equal(numeric, np.rint(numeric))
    )
    if invalid.any():
        examples = values.loc[invalid].head(5).tolist()
        raise ValueError(f"{context} contains invalid TIC values; first={examples}")
    return pd.Series(
        np.rint(numeric).astype(np.int64),
        index=values.index,
        name="tic",
    )


def _read_summary(path: Path) -> tuple[dict[str, Any], FileBinding]:
    binding = _bind_file(path)
    try:
        value = json.loads(binding.path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid frozen-pool summary: {binding.path}") from exc
    if not isinstance(value, dict):
        raise ValueError("frozen-pool summary must be a JSON object")
    if value.get("passed") is not True:
        raise ValueError("frozen-pool summary did not pass")
    if value.get("schema_version") != FULL_POOL_SUMMARY_SCHEMA_VERSION:
        raise ValueError("frozen-pool summary has the wrong schema")
    if value.get("pool_contract_version") != FULL_POOL_CONTRACT_VERSION:
        raise ValueError("frozen-pool summary has the wrong pool contract")
    if value.get("sectors") != list(FULL_POOL_NATIVE_SECTORS):
        raise ValueError("frozen-pool summary is not the S56--S62 release")
    return value, binding


def load_sector_pool_authority(
    *,
    sector: int,
    pool_path: Path,
    pool_summary_path: Path,
    allowlist_path: Path,
) -> SectorPoolAuthority:
    """Validate one staged sector against the frozen full-pool release.

    The pool and allowlist may be byte-identical staged copies on PDO.  Their
    paths therefore need not equal the ORCD paths recorded in the summary;
    their sizes, SHA-256 digests, row counts, and TIC identities must.
    """

    sector = int(sector)
    if sector not in FULL_POOL_NATIVE_SECTORS:
        raise ValueError(
            f"full-pool native preparation is bounded to {FULL_POOL_NATIVE_SECTORS}"
        )
    summary, summary_binding = _read_summary(pool_summary_path)
    pool_binding = _bind_file(pool_path)
    suffix = pool_binding.path.suffix.lower()
    output_name = "parquet" if suffix in {".parquet", ".pq"} else "csv"
    declared_pool = summary.get("outputs", {}).get(output_name)
    if not isinstance(declared_pool, dict):
        raise ValueError(
            f"frozen-pool summary lacks {output_name} output metadata"
        )
    if (
        declared_pool.get("sha256") != pool_binding.sha256
        or int(declared_pool.get("size_bytes", -1)) != pool_binding.size_bytes
    ):
        raise ValueError("staged frozen-pool table does not match its summary")

    pool = read_table(pool_binding.path)
    required = {
        "pool_contract_version",
        "sector",
        "tic",
        "camera",
        "ccd",
        "compact_h5_sha256",
        "compact_h5_size_bytes",
        "compact_group_path",
    }
    missing = sorted(required - set(pool.columns))
    if missing:
        raise KeyError(f"frozen-pool table lacks columns: {missing}")
    pool["sector"] = pd.to_numeric(pool["sector"], errors="raise").astype(
        np.int64
    )
    pool["tic"] = _positive_tics(pool["tic"], context="frozen-pool table")
    if pool.duplicated(["sector", "tic"]).any():
        raise ValueError("frozen-pool table contains duplicate observation keys")
    if not pool["pool_contract_version"].astype(str).eq(
        FULL_POOL_CONTRACT_VERSION
    ).all():
        raise ValueError("frozen-pool table has the wrong contract version")
    if int(declared_pool.get("n_rows", -1)) != len(pool):
        raise ValueError("frozen-pool table row count differs from summary")
    observed_sectors = tuple(sorted(pool["sector"].astype(int).unique()))
    if observed_sectors != FULL_POOL_NATIVE_SECTORS:
        raise ValueError("frozen-pool table is not exact S56--S62 coverage")
    retained = summary.get("counts", {}).get("retained", {})
    expected_retained = {
        "n_observations": int(len(pool)),
        "n_unique_tics": int(pool["tic"].nunique()),
        "n_multisector_tics": int(
            pool.groupby("tic")["sector"].nunique().gt(1).sum()
        ),
    }
    if retained != expected_retained:
        raise ValueError("frozen-pool retained counts differ from table")
    if summary.get("leakage_audit") != {
        "fixed_test_observations_retained": 0,
        "s63_reserved_observations_retained": 0,
    }:
        raise ValueError("frozen-pool leakage audit did not pass")
    if summary.get("identity_hashes", {}).get(
        "retained_observations_sha256"
    ) != _observation_identity_sha256(pool):
        raise ValueError("frozen-pool observation identity differs from summary")

    rows = (
        pool.loc[pool["sector"].eq(sector)]
        .sort_values("tic", kind="stable")
        .reset_index(drop=True)
    )
    if rows.empty:
        raise ValueError(f"frozen pool contains no S{sector} observations")
    for name in ("camera", "ccd"):
        numeric = pd.to_numeric(rows[name], errors="coerce")
        if numeric.isna().any() or not numeric.isin(range(1, 5)).all():
            raise ValueError(f"S{sector} frozen-pool {name} mapping is invalid")
        rows[name] = numeric.astype(np.int16)
    expected_group = rows["tic"].map(lambda tic: f"targets/{int(tic):016d}")
    if not rows["compact_group_path"].astype(str).eq(expected_group).all():
        raise ValueError(f"S{sector} compact group paths are not TIC keyed")
    compact_hashes = sorted(
        rows["compact_h5_sha256"].fillna("").astype(str).unique()
    )
    compact_sizes = sorted(
        pd.to_numeric(
            rows["compact_h5_size_bytes"], errors="coerce"
        ).dropna().astype(np.int64).unique()
    )
    if len(compact_hashes) != 1 or not pd.Series(compact_hashes).str.fullmatch(
        r"[0-9a-f]{64}"
    ).all():
        raise ValueError(f"S{sector} has ambiguous compact HDF5 hashes")
    if len(compact_sizes) != 1 or compact_sizes[0] <= 0:
        raise ValueError(f"S{sector} has ambiguous compact HDF5 sizes")

    compact_exports = summary.get("inputs", {}).get("compact_exports")
    if not isinstance(compact_exports, list):
        raise ValueError("frozen-pool summary lacks compact-export bindings")
    matches = [
        item
        for item in compact_exports
        if isinstance(item, dict) and int(item.get("sector", -1)) == sector
    ]
    if len(matches) != 1:
        raise ValueError(f"frozen-pool summary has ambiguous S{sector} export")
    declared_compact = matches[0].get("compact_h5")
    if not isinstance(declared_compact, dict):
        raise ValueError(f"frozen-pool summary lacks S{sector} compact HDF5")
    if (
        declared_compact.get("sha256") != compact_hashes[0]
        or int(declared_compact.get("size_bytes", -1)) != compact_sizes[0]
    ):
        raise ValueError(f"S{sector} compact binding disagrees with pool rows")

    declared_allowlist = summary.get("outputs", {}).get(
        "sector_allowlists", {}
    ).get(str(sector))
    if not isinstance(declared_allowlist, dict):
        raise ValueError(f"frozen-pool summary lacks S{sector} allowlist")
    allowlist_binding = _bind_file(
        allowlist_path,
        expected_sha256=str(declared_allowlist.get("sha256", "")),
        expected_size_bytes=int(declared_allowlist.get("size_bytes", -1)),
    )
    allowlist = pd.read_csv(allowlist_binding.path)
    if list(allowlist.columns) != ["tic"]:
        raise ValueError("sector allowlist must contain exactly one tic column")
    allowlist_tics = _positive_tics(
        allowlist["tic"], context=f"S{sector} allowlist"
    )
    if allowlist_tics.tolist() != sorted(set(allowlist_tics.tolist())):
        raise ValueError(f"S{sector} allowlist must be sorted and unique")
    expected_tics = rows["tic"].astype(np.int64).tolist()
    if allowlist_tics.tolist() != expected_tics:
        raise ValueError(f"S{sector} allowlist differs from frozen-pool rows")
    if (
        int(declared_allowlist.get("n_tics", -1)) != len(expected_tics)
        or declared_allowlist.get("tic_inventory_sha256")
        != _tic_inventory_sha256(expected_tics)
    ):
        raise ValueError(f"S{sector} allowlist metadata is inconsistent")

    identity = _observation_identity_sha256(rows)
    return SectorPoolAuthority(
        sector=sector,
        rows=rows,
        pool=pool_binding,
        pool_summary=summary_binding,
        allowlist=allowlist_binding,
        compact_h5_sha256=compact_hashes[0],
        compact_h5_size_bytes=int(compact_sizes[0]),
        observation_identity_sha256=identity,
    )


def shard_for_tic(*, sector: int, tic: int, n_shards: int) -> int:
    """Return the stable full-pool native shard for one observation."""

    if int(n_shards) < 1:
        raise ValueError("n_shards must be positive")
    digest = hashlib.sha256(
        f"s{int(sector):04d}-tic{int(tic):016d}".encode("ascii")
    ).digest()
    return int.from_bytes(digest[:8], "big") % int(n_shards)


def select_sector_shard(
    authority: SectorPoolAuthority,
    *,
    shard_index: int,
    n_shards: int,
) -> pd.DataFrame:
    """Select an exact deterministic shard from one sector authority."""

    shard_index = int(shard_index)
    n_shards = int(n_shards)
    if n_shards < 1 or shard_index < 0 or shard_index >= n_shards:
        raise ValueError("shard_index must satisfy 0 <= shard_index < n_shards")
    assignments = authority.rows["tic"].map(
        lambda tic: shard_for_tic(
            sector=authority.sector,
            tic=int(tic),
            n_shards=n_shards,
        )
    )
    selected = authority.rows.loc[assignments.eq(shard_index)].copy()
    if selected.empty:
        raise ValueError(
            f"S{authority.sector} shard {shard_index}/{n_shards} is empty"
        )
    return selected.sort_values("tic", kind="stable").reset_index(drop=True)


def _atomic_h5_path(path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite immutable output: {path}")
    return path.with_name(f".{path.name}.{os.getpid()}.tmp")


def _publish_immutable_bytes(path: Path, payload: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace immutable output: {path}"
            )
        return
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            if path.read_bytes() != payload:
                raise
    finally:
        temporary.unlink(missing_ok=True)


def export_full_pool_raw_source_shard(
    *,
    sector: int,
    pool_path: Path,
    pool_summary_path: Path,
    allowlist_path: Path,
    raw_root: Path,
    orbits: Sequence[int],
    out_h5: Path,
    shard_index: int = 0,
    n_shards: int = 1,
) -> dict[str, Any]:
    """Export one exact real-only raw/error shard on PDO."""

    import h5py

    normalized_orbits = tuple(int(value) for value in orbits)
    if (
        not normalized_orbits
        or len(set(normalized_orbits)) != len(normalized_orbits)
        or any(value <= 0 for value in normalized_orbits)
    ):
        raise ValueError("orbits must contain distinct positive integers")
    authority = load_sector_pool_authority(
        sector=sector,
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        allowlist_path=allowlist_path,
    )
    rows = select_sector_shard(
        authority,
        shard_index=shard_index,
        n_shards=n_shards,
    )
    out_h5 = Path(out_h5).expanduser().resolve()
    temporary = _atomic_h5_path(out_h5)
    failures: list[dict[str, Any]] = []
    try:
        with h5py.File(temporary, "w") as output:
            output.attrs["contract_version"] = (
                harmonic_export.RAW_SOURCE_CONTRACT_VERSION
            )
            output.attrs["full_pool_raw_source_contract_version"] = (
                FULL_POOL_RAW_SOURCE_CONTRACT_VERSION
            )
            output.attrs["created_utc"] = datetime.now(
                timezone.utc
            ).isoformat()
            output.attrs["sector"] = int(sector)
            output.attrs["shard_index"] = int(shard_index)
            output.attrs["n_shards"] = int(n_shards)
            output.attrs["orbits"] = json.dumps(list(normalized_orbits))
            output.attrs["time_system"] = "BJD"
            output.attrs["real_only"] = 1
            output.attrs["full_pool_path"] = str(authority.pool.path)
            output.attrs["full_pool_sha256"] = authority.pool.sha256
            output.attrs["full_pool_summary_path"] = str(
                authority.pool_summary.path
            )
            output.attrs["full_pool_summary_sha256"] = (
                authority.pool_summary.sha256
            )
            output.attrs["sector_allowlist_path"] = str(
                authority.allowlist.path
            )
            output.attrs["sector_allowlist_sha256"] = (
                authority.allowlist.sha256
            )
            output.attrs["sector_observation_identity_sha256"] = (
                authority.observation_identity_sha256
            )
            output.attrs["shard_observation_identity_sha256"] = (
                _observation_identity_sha256(rows)
            )
            output.attrs["n_sector_observations"] = int(len(authority.rows))
            output.attrs["n_shard_observations"] = int(len(rows))
            targets = output.create_group("targets")
            for count, row in enumerate(rows.itertuples(index=False), start=1):
                tic = int(row.tic)
                camera = int(row.camera)
                ccd = int(row.ccd)
                paths = harmonic_export.discover_tglc_paths(
                    raw_root=Path(raw_root),
                    tic=tic,
                    camera=camera,
                    ccd=ccd,
                    orbits=normalized_orbits,
                )
                try:
                    if not paths:
                        raise FileNotFoundError(
                            f"no raw TGLC HDF5 paths found for TIC {tic}"
                        )
                    payload = harmonic_export.merge_tglc_raw_paths(paths)
                    group = targets.create_group(f"{tic:016d}")
                    group.attrs["tic"] = tic
                    group.attrs["sector"] = int(sector)
                    group.attrs["camera"] = camera
                    group.attrs["ccd"] = ccd
                    group.attrs["source_paths"] = json.dumps(
                        [str(path) for path in paths]
                    )
                    for name, values in payload.items():
                        harmonic_export._write_dataset(group, name, values)
                except Exception as exc:
                    failures.append(
                        {
                            "tic": tic,
                            "error": str(exc),
                            "paths": "|".join(str(path) for path in paths),
                        }
                    )
                if count % 100 == 0:
                    print(
                        f"[fullpool-raw] S{sector} shard "
                        f"{shard_index}/{n_shards}: {count}/{len(rows)}",
                        flush=True,
                    )
        if failures:
            pd.DataFrame(failures).to_csv(
                out_h5.with_suffix(".failures.csv"), index=False
            )
            raise RuntimeError(
                f"raw export failed for {len(failures)} of {len(rows)} TICs"
            )
        for binding in (
            authority.pool,
            authority.pool_summary,
            authority.allowlist,
        ):
            binding.assert_unchanged()
        temporary.replace(out_h5)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    failures_path = out_h5.with_suffix(".failures.csv")
    failures_path.unlink(missing_ok=True)
    return {
        "schema_version": FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION,
        "stage": "pdo_raw_source_export",
        "sector": int(sector),
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
        "n_sector_observations": int(len(authority.rows)),
        "n_shard_observations": int(len(rows)),
        "sector_observation_identity_sha256": (
            authority.observation_identity_sha256
        ),
        "shard_observation_identity_sha256": (
            _observation_identity_sha256(rows)
        ),
        "full_pool_sha256": authority.pool.sha256,
        "full_pool_summary_sha256": authority.pool_summary.sha256,
        "sector_allowlist_sha256": authority.allowlist.sha256,
        "out_h5": str(out_h5),
        "out_h5_size_bytes": int(out_h5.stat().st_size),
        "out_h5_sha256": file_sha256(out_h5),
        "real_only": True,
    }


def _full_pool_reconciliation_summary(
    *,
    policy: str,
    stats: Mapping[str, Any],
) -> dict[str, Any]:
    summary = harmonic_export._orbitid_reconciliation_summary(
        policy=policy,
        stats=stats,
    )
    if policy == ORBITID_POLICY_REFERENCE:
        summary["scope"] = "teacher_ssl_fullpool_s62_real_only"
        summary["release_binding"] = FULL_POOL_NATIVE_RELEASE_BINDING
    return summary


def _validate_s62_group_reconciliation(
    *,
    tic: int,
    cadenceno: np.ndarray,
    compact_orbitid: np.ndarray,
    raw_orbitid: np.ndarray,
    resolved_orbitid: np.ndarray,
    correction_mask: np.ndarray,
) -> None:
    """Validate generic S62 reconciliation without a release row-count bind."""

    cadences = np.asarray(cadenceno, dtype=np.int64)
    compact = np.asarray(compact_orbitid, dtype=np.int64)
    raw = np.asarray(raw_orbitid, dtype=np.int64)
    resolved = np.asarray(resolved_orbitid, dtype=np.int64)
    corrected = np.asarray(correction_mask, dtype=bool)
    if len({len(cadences), len(compact), len(raw), len(resolved), len(corrected)}) != 1:
        raise ValueError(f"TIC {tic}: orbit-ID reconciliation lengths differ")
    if not np.array_equal(raw, compact):
        raise ValueError(
            f"TIC {tic}: raw-source and compact orbitid arrays disagree "
            "before reconciliation"
        )
    mismatch = compact != resolved
    if not np.array_equal(mismatch, corrected):
        raise ValueError(
            f"TIC {tic}: not all and only authority mismatches were corrected"
        )
    if np.any(compact[corrected] != S62_CORRECTION_INPUT_ORBIT) or np.any(
        resolved[corrected] != S62_CORRECTION_REFERENCE_ORBIT
    ):
        raise ValueError(f"TIC {tic}: correction is not the bounded 132->131 map")
    corrected_cadences = cadences[corrected]
    if corrected_cadences.size and (
        int(np.min(corrected_cadences)) < S62_CORRECTION_CADENCE_MIN
        or int(np.max(corrected_cadences)) > S62_CORRECTION_CADENCE_MAX
    ):
        raise ValueError(f"TIC {tic}: correction exceeds the S62 cadence bound")
    reference_131 = (
        (cadences >= S62_REFERENCE_ORBIT_131_MIN)
        & (cadences <= S62_REFERENCE_ORBIT_131_MAX)
    )
    reference_132 = (
        (cadences >= S62_REFERENCE_ORBIT_132_MIN)
        & (cadences <= S62_REFERENCE_ORBIT_132_MAX)
    )
    if np.any(~(reference_131 | reference_132)):
        raise ValueError(
            f"TIC {tic}: cadence lies outside the bounded S62 authority"
        )
    if np.any(resolved[reference_131] != 131) or np.any(
        resolved[reference_132] != 132
    ):
        raise ValueError(
            f"TIC {tic}: resolved orbitid disagrees with S62 cadence authority"
        )


def build_full_pool_native_shard(
    *,
    sector: int,
    pool_path: Path,
    pool_summary_path: Path,
    allowlist_path: Path,
    raw_source_h5: Path,
    compact_adp_h5: Path,
    cadence_reference_table: Path,
    cadence_reference_manifest: Path,
    out_h5: Path,
    shard_index: int = 0,
    n_shards: int = 1,
    n_periods: int = harmonic_export.DEFAULT_BLS_PERIODS,
    orbitid_policy: str = ORBITID_POLICY_STRICT,
) -> dict[str, Any]:
    """Build one real-only native HDF5 shard on ORCD."""

    import h5py

    sector = int(sector)
    if orbitid_policy not in {
        ORBITID_POLICY_STRICT,
        ORBITID_POLICY_REFERENCE,
    }:
        raise ValueError("orbitid_policy must be strict or reference_by_cadence")
    if orbitid_policy == ORBITID_POLICY_REFERENCE and sector != 62:
        raise ValueError("reference_by_cadence is bounded to S62")
    if sector == 62 and orbitid_policy != ORBITID_POLICY_REFERENCE:
        raise ValueError("S62 full-pool native preparation requires reference_by_cadence")
    if sector != 62 and orbitid_policy != ORBITID_POLICY_STRICT:
        raise ValueError("S56--S61 full-pool native preparation requires strict")
    if int(n_periods) < 1:
        raise ValueError("n_periods must be positive")

    authority = load_sector_pool_authority(
        sector=sector,
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        allowlist_path=allowlist_path,
    )
    rows = select_sector_shard(
        authority,
        shard_index=shard_index,
        n_shards=n_shards,
    )
    raw_binding = _bind_file(raw_source_h5)
    compact_binding = _bind_file(
        compact_adp_h5,
        expected_sha256=authority.compact_h5_sha256,
        expected_size_bytes=authority.compact_h5_size_bytes,
    )
    reference_table_binding = _bind_file(cadence_reference_table)
    reference_manifest_binding = _bind_file(cadence_reference_manifest)
    quality_reference = load_external_quality_reference(
        table_path=reference_table_binding.path,
        manifest_path=reference_manifest_binding.path,
        sector=sector,
        expected_orbits=S56_EXPECTED_ORBITS if sector == 56 else None,
        expected_detectors=S56_EXPECTED_DETECTORS,
    )
    if (
        quality_reference.table_sha256 != reference_table_binding.sha256
        or quality_reference.manifest_sha256
        != reference_manifest_binding.sha256
    ):
        raise RuntimeError("cadence-reference loader changed input identity")

    out_h5 = Path(out_h5).expanduser().resolve()
    temporary = _atomic_h5_path(out_h5)
    failures: list[dict[str, Any]] = []
    quality_totals = {name: 0 for name in RAW_PAIR_QUALITY_COUNT_NAMES}
    orbitid_stats = harmonic_export._new_orbitid_reconciliation_stats()
    orbitid_summary: dict[str, Any] = {}
    try:
        with (
            h5py.File(raw_binding.path, "r") as raw_file,
            h5py.File(compact_binding.path, "r") as compact_file,
            h5py.File(temporary, "w") as output,
        ):
            if str(raw_file.attrs.get("contract_version", "")) != (
                harmonic_export.RAW_SOURCE_CONTRACT_VERSION
            ):
                raise ValueError("raw source has the wrong base contract")
            if str(
                raw_file.attrs.get(
                    "full_pool_raw_source_contract_version", ""
                )
            ) != FULL_POOL_RAW_SOURCE_CONTRACT_VERSION:
                raise ValueError("raw source lacks the full-pool contract")
            expected_raw_attrs: dict[str, Any] = {
                "sector": sector,
                "shard_index": int(shard_index),
                "n_shards": int(n_shards),
                "n_sector_observations": int(len(authority.rows)),
                "n_shard_observations": int(len(rows)),
                "full_pool_sha256": authority.pool.sha256,
                "full_pool_summary_sha256": authority.pool_summary.sha256,
                "sector_allowlist_sha256": authority.allowlist.sha256,
                "sector_observation_identity_sha256": (
                    authority.observation_identity_sha256
                ),
                "shard_observation_identity_sha256": (
                    _observation_identity_sha256(rows)
                ),
                "real_only": 1,
            }
            mismatches = {
                name: {
                    "expected": expected,
                    "observed": raw_file.attrs.get(name),
                }
                for name, expected in expected_raw_attrs.items()
                if raw_file.attrs.get(name) != expected
            }
            if mismatches:
                raise ValueError(
                    "raw-source full-pool binding mismatch: "
                    + json.dumps(mismatches, sort_keys=True)
                )
            raw_tics = sorted(int(value) for value in raw_file["targets"])
            expected_tics = rows["tic"].astype(int).tolist()
            if raw_tics != expected_tics:
                raise ValueError("raw-source TIC inventory differs from shard")

            output.attrs["contract_version"] = (
                FULL_POOL_NATIVE_CONTRACT_VERSION
            )
            output.attrs["created_utc"] = datetime.now(
                timezone.utc
            ).isoformat()
            output.attrs["time_system"] = "BJD"
            output.attrs["real_only"] = 1
            output.attrs["sector"] = sector
            output.attrs["shard_index"] = int(shard_index)
            output.attrs["n_shards"] = int(n_shards)
            output.attrs["n_sector_observations"] = int(len(authority.rows))
            output.attrs["n_shard_observations"] = int(len(rows))
            output.attrs["full_pool_path"] = str(authority.pool.path)
            output.attrs["full_pool_sha256"] = authority.pool.sha256
            output.attrs["full_pool_summary_path"] = str(
                authority.pool_summary.path
            )
            output.attrs["full_pool_summary_sha256"] = (
                authority.pool_summary.sha256
            )
            output.attrs["sector_allowlist_path"] = str(
                authority.allowlist.path
            )
            output.attrs["sector_allowlist_sha256"] = (
                authority.allowlist.sha256
            )
            output.attrs["sector_observation_identity_sha256"] = (
                authority.observation_identity_sha256
            )
            output.attrs["shard_observation_identity_sha256"] = (
                _observation_identity_sha256(rows)
            )
            output.attrs["raw_source_h5"] = str(raw_binding.path)
            output.attrs["raw_source_h5_sha256"] = raw_binding.sha256
            output.attrs["compact_adp_h5"] = str(compact_binding.path)
            output.attrs["compact_adp_h5_sha256"] = compact_binding.sha256
            output.attrs["periodogram_grid"] = "log10_period_d"
            output.attrs["periodogram_n"] = int(n_periods)
            output.attrs["chronology_small_channels"] = json.dumps(
                CHRONOLOGY_SMALL_CHANNELS
            )
            output.attrs["chronology_supplemental_channels"] = json.dumps(
                CHRONOLOGY_SUPPLEMENTAL_CHANNELS
            )
            output.attrs["harmonic_view_channels"] = json.dumps(
                HARMONIC_VIEW_CHANNELS
            )
            output.attrs["periodogram_channels"] = json.dumps(
                PERIODOGRAM_CHANNELS
            )
            harmonic_export._write_external_quality_attrs(
                output, quality_reference
            )
            targets = output.create_group("targets")
            for count, row in enumerate(rows.itertuples(index=False), start=1):
                tic = int(row.tic)
                raw_path = f"targets/{tic:016d}"
                compact_path = str(row.compact_group_path)
                try:
                    raw_group = raw_file[raw_path]
                    compact_group = compact_file[compact_path]
                    for name, expected in (
                        ("sector", sector),
                        ("tic", tic),
                        ("camera", int(row.camera)),
                        ("ccd", int(row.ccd)),
                    ):
                        if int(compact_group.attrs.get(name, -1)) != expected:
                            raise ValueError(
                                f"compact {name} mapping differs from pool"
                            )
                    raw = harmonic_export._raw_source_payload(raw_group)
                    aligned = harmonic_export.align_raw_by_cadence(
                        raw,
                        cadenceno=np.asarray(compact_group["cadenceno"]),
                        time=np.asarray(compact_group["time"]),
                    )
                    compact_orbitid = np.asarray(
                        compact_group["orbitid"], dtype=np.int64
                    )
                    raw_orbitid = np.asarray(
                        aligned["orbitid"], dtype=np.int64
                    )
                    if not np.array_equal(raw_orbitid, compact_orbitid):
                        raise ValueError(
                            "raw-source and compact orbitid arrays disagree "
                            "before reconciliation"
                        )
                    quality_overlay = quality_reference.apply(
                        sector=sector,
                        camera=int(row.camera),
                        ccd=int(row.ccd),
                        cadenceno=np.asarray(compact_group["cadenceno"]),
                        orbitid=compact_orbitid,
                        internal_quality=np.asarray(
                            compact_group["quality"]
                        ),
                        context=f"S{sector} TIC {tic}",
                        orbitid_policy=orbitid_policy,
                    )
                    correction_mask = np.asarray(
                        quality_overlay.orbitid_reference_correction_mask,
                        dtype=bool,
                    )
                    resolved_orbitid = np.asarray(
                        quality_overlay.resolved_orbitid, dtype=np.int64
                    )
                    if orbitid_policy == ORBITID_POLICY_REFERENCE:
                        _validate_s62_group_reconciliation(
                            tic=tic,
                            cadenceno=np.asarray(
                                compact_group["cadenceno"]
                            ),
                            compact_orbitid=compact_orbitid,
                            raw_orbitid=raw_orbitid,
                            resolved_orbitid=resolved_orbitid,
                            correction_mask=correction_mask,
                        )
                    elif correction_mask.any() or not np.array_equal(
                        resolved_orbitid, compact_orbitid
                    ):
                        raise ValueError(
                            "strict orbit-ID policy changed one or more cadences"
                        )

                    corrected = harmonic_export._record_orbitid_reconciliation(
                        orbitid_stats,
                        camera=int(row.camera),
                        ccd=int(row.ccd),
                        cadenceno=np.asarray(compact_group["cadenceno"]),
                        input_orbitid=compact_orbitid,
                        resolved_orbitid=resolved_orbitid,
                        correction_mask=correction_mask,
                    )
                    group = targets.create_group(f"{tic:016d}")
                    harmonic_export._copy_attrs(compact_group, group)
                    harmonic_export._write_quality_counts(
                        group, quality_overlay.counts
                    )
                    group.attrs["orbitid_reconciliation_policy"] = (
                        orbitid_policy
                    )
                    group.attrs[
                        "n_cad_orbitid_reference_corrected"
                    ] = int(corrected)
                    group.attrs[
                        "orbitid_reconciliation_source_agreement"
                    ] = 1
                    group.attrs["raw_source_paths"] = raw_group.attrs.get(
                        "source_paths", ""
                    )
                    group.attrs["raw_adp_time_delta_max_s"] = float(
                        np.nanmax(aligned["_time_delta_s"])
                    )
                    payload = {
                        "time": harmonic_export._absolute_bjd(
                            np.asarray(compact_group["time"])
                        ),
                        "cadenceno": np.asarray(
                            compact_group["cadenceno"], dtype=np.int64
                        ),
                        "orbitid": resolved_orbitid.astype(np.int32),
                        "orbitid_pre_reference": compact_orbitid.astype(
                            np.int32
                        ),
                        "raw_orbitid_pre_reference": raw_orbitid.astype(
                            np.int32
                        ),
                        ORBITID_RECONCILIATION_MASK_DATASET: correction_mask.astype(
                            np.uint8
                        ),
                        "quality": np.asarray(
                            quality_overlay.quality, dtype=np.int32
                        ),
                        "raw_flux_small": aligned["raw_flux_small"],
                        "raw_flux_err_small": aligned[
                            "raw_flux_err_small"
                        ],
                        "raw_flux_primary": aligned["raw_flux_primary"],
                        "raw_flux_err_primary": aligned[
                            "raw_flux_err_primary"
                        ],
                        "det_flux_adp_sml": np.asarray(
                            compact_group["DET_FLUX_ADP_SML"]
                        ),
                        "det_flux_adp": np.asarray(
                            compact_group["DET_FLUX_ADP"]
                        ),
                    }
                    payload.update(
                        harmonic_export._periodogram_payload(
                            time=payload["time"],
                            quality=payload["quality"],
                            small=payload["det_flux_adp_sml"],
                            primary=payload["det_flux_adp"],
                            n_periods=int(n_periods),
                        )
                    )
                    for name, values in payload.items():
                        harmonic_export._write_dataset(group, name, values)
                    harmonic_export._add_quality_counts(
                        quality_totals, quality_overlay.counts
                    )
                except Exception as exc:
                    failures.append(
                        {"sector": sector, "tic": tic, "error": str(exc)}
                    )
                if count % 100 == 0:
                    print(
                        f"[fullpool-native] S{sector} shard "
                        f"{shard_index}/{n_shards}: {count}/{len(rows)}",
                        flush=True,
                    )
            for name, value in quality_totals.items():
                output.attrs[f"quality_overlay_{name}"] = int(value)
            orbitid_summary = _full_pool_reconciliation_summary(
                policy=orbitid_policy,
                stats=orbitid_stats,
            )
            harmonic_export._write_orbitid_reconciliation_attrs(
                output, orbitid_summary
            )
        if failures:
            pd.DataFrame(failures).to_csv(
                out_h5.with_suffix(".failures.csv"), index=False
            )
            raise RuntimeError(
                f"native export failed for {len(failures)} of {len(rows)} TICs"
            )
        for binding in (
            authority.pool,
            authority.pool_summary,
            authority.allowlist,
            raw_binding,
            compact_binding,
            reference_table_binding,
            reference_manifest_binding,
        ):
            binding.assert_unchanged()
        quality_reference.assert_unchanged()
        temporary.replace(out_h5)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    failures_path = out_h5.with_suffix(".failures.csv")
    failures_path.unlink(missing_ok=True)
    verification = verify_full_pool_native_shard(
        out_h5,
        expected_sector=sector,
        expected_shard_index=int(shard_index),
        expected_n_shards=int(n_shards),
        expected_observations=len(rows),
        require_periodograms=True,
    )
    if not verification["passed"]:
        raise RuntimeError(
            "published full-pool native shard failed verification: "
            + "; ".join(verification["failures"][:10])
        )
    return {
        "schema_version": FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION,
        "stage": "orcd_native_shard_build",
        "sector": sector,
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
        "n_sector_observations": int(len(authority.rows)),
        "n_shard_observations": int(len(rows)),
        "sector_observation_identity_sha256": (
            authority.observation_identity_sha256
        ),
        "shard_observation_identity_sha256": (
            _observation_identity_sha256(rows)
        ),
        "full_pool_sha256": authority.pool.sha256,
        "full_pool_summary_sha256": authority.pool_summary.sha256,
        "sector_allowlist_sha256": authority.allowlist.sha256,
        "raw_source_h5_sha256": raw_binding.sha256,
        "compact_adp_h5_sha256": compact_binding.sha256,
        "cadence_reference_table_sha256": reference_table_binding.sha256,
        "cadence_reference_manifest_sha256": (
            reference_manifest_binding.sha256
        ),
        "orbitid_reconciliation": orbitid_summary,
        "out_h5": str(out_h5),
        "out_h5_size_bytes": int(out_h5.stat().st_size),
        "out_h5_sha256": file_sha256(out_h5),
        "native_contract_version": FULL_POOL_NATIVE_CONTRACT_VERSION,
        "real_only": True,
        "verification": verification,
    }


def full_pool_native_root_failures(attrs: Mapping[str, Any]) -> list[str]:
    """Validate the generic full-pool orbit/provenance root envelope."""

    failures: list[str] = []
    required = {
        "real_only",
        "sector",
        "shard_index",
        "n_shards",
        "n_sector_observations",
        "n_shard_observations",
        "full_pool_sha256",
        "full_pool_summary_sha256",
        "sector_allowlist_sha256",
        "sector_observation_identity_sha256",
        "shard_observation_identity_sha256",
        "raw_source_h5_sha256",
        "compact_adp_h5_sha256",
        *RAW_PAIR_ORBITID_RECONCILIATION_ATTRS,
    }
    missing = sorted(name for name in required if name not in attrs)
    if missing:
        return [f"full-pool native root lacks attrs: {','.join(missing)}"]
    try:
        sector = int(attrs["sector"])
        shard_index = int(attrs["shard_index"])
        n_shards = int(attrs["n_shards"])
        n_sector = int(attrs["n_sector_observations"])
        n_shard = int(attrs["n_shard_observations"])
        real_only = int(attrs["real_only"])
    except (TypeError, ValueError):
        return ["full-pool native root has invalid integer metadata"]
    if (
        sector not in FULL_POOL_NATIVE_SECTORS
        or n_shards < 1
        or shard_index < 0
        or shard_index >= n_shards
        or n_sector < 1
        or n_shard < 1
        or n_shard > n_sector
        or real_only != 1
    ):
        failures.append("full-pool native root metadata is outside bounds")
    for name in (
        "full_pool_sha256",
        "full_pool_summary_sha256",
        "sector_allowlist_sha256",
        "sector_observation_identity_sha256",
        "shard_observation_identity_sha256",
        "raw_source_h5_sha256",
        "compact_adp_h5_sha256",
    ):
        if not isinstance(attrs[name], (str, bytes)) or not str(
            attrs[name]
        ).strip("b'").strip("'").lower().isalnum():
            failures.append(f"{name} is invalid")
        text = str(attrs[name])
        if len(text) != 64 or any(value not in "0123456789abcdef" for value in text):
            failures.append(f"{name} is not a lowercase SHA-256")

    if str(attrs["orbitid_reconciliation_contract_version"]) != (
        ORBITID_RECONCILIATION_CONTRACT_VERSION
    ):
        failures.append("orbitid reconciliation contract mismatch")
    if str(attrs["orbitid_reconciliation_authority"]) != (
        EXPECTED_CADENCE_AUTHORITY
    ):
        failures.append("orbitid reconciliation authority mismatch")
    policy = str(attrs["orbitid_reconciliation_policy"])
    integer_names = (
        "orbitid_reconciliation_bounded_sector",
        "orbitid_reconciliation_source_agreement_required",
        "orbitid_reconciliation_n_groups_examined",
        "orbitid_reconciliation_n_groups_corrected",
        "orbitid_reconciliation_n_groups_unmodified",
        "orbitid_reconciliation_n_cad_corrected",
        "orbitid_reconciliation_min_cadenceno_corrected",
        "orbitid_reconciliation_max_cadenceno_corrected",
        "orbitid_reconciliation_input_orbitid",
        "orbitid_reconciliation_reference_orbitid",
    )
    try:
        integers = {name: int(attrs[name]) for name in integer_names}
    except (TypeError, ValueError):
        failures.append("orbitid reconciliation root integers are invalid")
        return failures
    if integers["orbitid_reconciliation_n_groups_examined"] != n_shard:
        failures.append("orbitid reconciliation examined count differs from shard")
    if (
        integers["orbitid_reconciliation_n_groups_corrected"]
        + integers["orbitid_reconciliation_n_groups_unmodified"]
        != integers["orbitid_reconciliation_n_groups_examined"]
    ):
        failures.append("orbitid reconciliation group arithmetic is inconsistent")
    mapping_names = (
        "orbitid_reconciliation_n_cad_corrected_by_camera",
        "orbitid_reconciliation_n_groups_corrected_by_camera",
        "orbitid_reconciliation_n_cad_corrected_by_detector",
    )
    mappings: dict[str, dict[str, int]] = {}
    for name in mapping_names:
        try:
            value = json.loads(str(attrs[name]))
            if not isinstance(value, dict):
                raise TypeError
            mappings[name] = {str(key): int(count) for key, count in value.items()}
            if any(count < 0 for count in mappings[name].values()):
                raise ValueError
        except (TypeError, ValueError, json.JSONDecodeError):
            failures.append(f"{name} is invalid")
    if mappings:
        cadence_total = integers["orbitid_reconciliation_n_cad_corrected"]
        group_total = integers["orbitid_reconciliation_n_groups_corrected"]
        for name in (
            "orbitid_reconciliation_n_cad_corrected_by_camera",
            "orbitid_reconciliation_n_cad_corrected_by_detector",
        ):
            if name in mappings and sum(mappings[name].values()) != cadence_total:
                failures.append(f"{name} does not sum to cadence total")
        group_name = "orbitid_reconciliation_n_groups_corrected_by_camera"
        if group_name in mappings and sum(mappings[group_name].values()) != group_total:
            failures.append(f"{group_name} does not sum to group total")

    if policy == ORBITID_POLICY_STRICT:
        if sector == 62:
            failures.append("S62 may not use strict orbit-ID policy")
        expected = {
            "orbitid_reconciliation_scope": "strict_input_match",
            "orbitid_reconciliation_release_binding": "",
            "orbitid_reconciliation_bounded_sector": -1,
            "orbitid_reconciliation_source_agreement_required": 0,
            "orbitid_reconciliation_n_groups_corrected": 0,
            "orbitid_reconciliation_n_cad_corrected": 0,
            "orbitid_reconciliation_min_cadenceno_corrected": -1,
            "orbitid_reconciliation_max_cadenceno_corrected": -1,
            "orbitid_reconciliation_input_orbitid": -1,
            "orbitid_reconciliation_reference_orbitid": -1,
        }
    elif policy == ORBITID_POLICY_REFERENCE:
        if sector != 62:
            failures.append("reference orbit-ID policy is outside S62")
        expected = {
            "orbitid_reconciliation_scope": (
                "teacher_ssl_fullpool_s62_real_only"
            ),
            "orbitid_reconciliation_release_binding": (
                FULL_POOL_NATIVE_RELEASE_BINDING
            ),
            "orbitid_reconciliation_bounded_sector": 62,
            "orbitid_reconciliation_source_agreement_required": 1,
        }
        corrected = integers["orbitid_reconciliation_n_cad_corrected"]
        if corrected:
            if (
                integers["orbitid_reconciliation_input_orbitid"]
                != S62_CORRECTION_INPUT_ORBIT
                or integers["orbitid_reconciliation_reference_orbitid"]
                != S62_CORRECTION_REFERENCE_ORBIT
                or integers["orbitid_reconciliation_min_cadenceno_corrected"]
                < S62_CORRECTION_CADENCE_MIN
                or integers["orbitid_reconciliation_max_cadenceno_corrected"]
                > S62_CORRECTION_CADENCE_MAX
            ):
                failures.append("S62 correction signature exceeds generic bounds")
        else:
            for name in (
                "orbitid_reconciliation_input_orbitid",
                "orbitid_reconciliation_reference_orbitid",
                "orbitid_reconciliation_min_cadenceno_corrected",
                "orbitid_reconciliation_max_cadenceno_corrected",
            ):
                if integers[name] != -1:
                    failures.append(f"{name} must be -1 without corrections")
    else:
        return failures + ["full-pool orbit-ID policy is invalid"]
    for name, expected_value in expected.items():
        if name in integers:
            observed: Any = integers[name]
        else:
            observed = str(attrs[name])
        if observed != expected_value:
            failures.append(f"{name} differs from full-pool policy")
    return failures


def full_pool_native_group_failures(
    group: Any,
    *,
    context: str,
    root_policy: str,
) -> tuple[list[str], int]:
    """Validate one group's auditable pre/post-reconciliation arrays."""

    failures: list[str] = []
    required_attrs = (
        "orbitid_reconciliation_policy",
        "n_cad_orbitid_reference_corrected",
        "orbitid_reconciliation_source_agreement",
        "sector",
        "tic",
    )
    missing_attrs = [name for name in required_attrs if name not in group.attrs]
    if missing_attrs:
        return [f"{context}:missing attrs {','.join(missing_attrs)}"], 0
    try:
        corrected = int(group.attrs["n_cad_orbitid_reference_corrected"])
        agreement = int(group.attrs["orbitid_reconciliation_source_agreement"])
        sector = int(group.attrs["sector"])
        tic = int(group.attrs["tic"])
    except (TypeError, ValueError):
        return [f"{context}:invalid orbit-ID audit attrs"], 0
    if agreement != 1:
        failures.append(f"{context}:raw/compact source agreement failed")
    if str(group.attrs["orbitid_reconciliation_policy"]) != root_policy:
        failures.append(f"{context}:orbit-ID policy differs from root")
    required_datasets = (
        "cadenceno",
        "orbitid",
        "orbitid_pre_reference",
        "raw_orbitid_pre_reference",
        ORBITID_RECONCILIATION_MASK_DATASET,
    )
    missing = [name for name in required_datasets if name not in group]
    if missing:
        return failures + [f"{context}:missing={','.join(missing)}"], corrected
    cadences = np.asarray(group["cadenceno"], dtype=np.int64)
    stored = np.asarray(group["orbitid"], dtype=np.int64)
    compact_pre = np.asarray(group["orbitid_pre_reference"], dtype=np.int64)
    raw_pre = np.asarray(group["raw_orbitid_pre_reference"], dtype=np.int64)
    raw_mask = np.asarray(group[ORBITID_RECONCILIATION_MASK_DATASET])
    if (
        len({len(cadences), len(stored), len(compact_pre), len(raw_pre), len(raw_mask)})
        != 1
    ):
        return failures + [f"{context}:orbit-ID arrays differ in length"], corrected
    if not np.isin(raw_mask, (0, 1, False, True)).all():
        return failures + [f"{context}:correction mask is not binary"], corrected
    mask = raw_mask.astype(bool)
    if int(np.count_nonzero(mask)) != corrected:
        failures.append(f"{context}:correction count differs from mask")
    if not np.array_equal(raw_pre, compact_pre):
        failures.append(f"{context}:raw/compact pre-correction orbitids differ")
    if root_policy == ORBITID_POLICY_STRICT:
        if corrected or mask.any() or not np.array_equal(compact_pre, stored):
            failures.append(f"{context}:strict policy changed orbit IDs")
    elif root_policy == ORBITID_POLICY_REFERENCE:
        if sector != 62:
            failures.append(f"{context}:reference correction is outside S62")
        try:
            _validate_s62_group_reconciliation(
                tic=tic,
                cadenceno=cadences,
                compact_orbitid=compact_pre,
                raw_orbitid=raw_pre,
                resolved_orbitid=stored,
                correction_mask=mask,
            )
        except ValueError as exc:
            failures.append(f"{context}:{exc}")
    else:
        failures.append(f"{context}:invalid root orbit-ID policy")
    return failures, corrected


def verify_full_pool_native_shard(
    path: Path,
    *,
    expected_sector: int | None = None,
    expected_shard_index: int | None = None,
    expected_n_shards: int | None = None,
    expected_observations: int | None = None,
    require_periodograms: bool = True,
) -> dict[str, Any]:
    """Read and validate every real group in one full-pool native shard."""

    import h5py

    path = Path(path)
    failures: list[str] = []
    counts = {"targets": 0, "injections": 0}
    if not path.is_file():
        return {
            "passed": False,
            "failures": [f"missing file: {path}"],
            "counts": counts,
        }
    with h5py.File(path, "r") as h5:
        if str(h5.attrs.get("contract_version", "")) != (
            FULL_POOL_NATIVE_CONTRACT_VERSION
        ):
            failures.append("wrong full-pool native contract")
        failures.extend(full_pool_native_root_failures(h5.attrs))
        for name, expected in (
            ("sector", expected_sector),
            ("shard_index", expected_shard_index),
            ("n_shards", expected_n_shards),
            ("n_shard_observations", expected_observations),
        ):
            if expected is not None and int(h5.attrs.get(name, -1)) != int(expected):
                failures.append(f"{name} differs from expected value")
        if "injections" in h5 and len(h5["injections"]):
            counts["injections"] = len(h5["injections"])
            failures.append("full-pool native shard contains injections")
        if "targets" not in h5:
            failures.append("full-pool native shard lacks targets")
        else:
            counts["targets"] = len(h5["targets"])
            for key, group in h5["targets"].items():
                group_path = f"targets/{key}"
                try:
                    light_curve = read_native_light_curve_from_h5(
                        h5,
                        group_path=group_path,
                        require_errors=True,
                    )
                    if require_periodograms and any(
                        len(getattr(light_curve, name)) == 0
                        for name in (
                            "bls_power_small",
                            "bls_sde_small",
                            "bls_power_primary",
                            "bls_sde_primary",
                        )
                    ):
                        failures.append(f"/{group_path}:missing periodograms")
                    missing_native = [
                        name for name in NATIVE_DATASETS if name not in group
                    ]
                    missing_periodogram = [
                        name for name in PERIODOGRAM_DATASETS if name not in group
                    ]
                    if missing_native:
                        failures.append(
                            f"/{group_path}:missing={','.join(missing_native)}"
                        )
                    if require_periodograms and missing_periodogram:
                        failures.append(
                            f"/{group_path}:missing={','.join(missing_periodogram)}"
                        )
                except Exception as exc:
                    failures.append(f"/{group_path}:{exc}")
        if counts["targets"] != int(
            h5.attrs.get("n_shard_observations", -1)
        ):
            failures.append("target count differs from root declaration")
    return {
        "passed": not failures,
        "failures": failures,
        "counts": counts,
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "size_bytes": int(path.stat().st_size),
    }


def write_full_pool_native_registry(
    *,
    pool_path: Path,
    pool_summary_path: Path,
    native_shard_paths: Sequence[Path],
    source_path: Path,
    registry_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    """Freeze exact shard coverage into an observation-keyed native registry."""

    pool_binding = _bind_file(pool_path)
    summary, summary_binding = _read_summary(pool_summary_path)
    suffix = pool_binding.path.suffix.lower()
    output_name = "parquet" if suffix in {".parquet", ".pq"} else "csv"
    declared = summary.get("outputs", {}).get(output_name, {})
    if (
        declared.get("sha256") != pool_binding.sha256
        or int(declared.get("size_bytes", -1)) != pool_binding.size_bytes
    ):
        raise ValueError("frozen pool does not match its summary")
    pool = read_table(pool_binding.path)
    if "pool_contract_version" not in pool:
        raise KeyError("frozen pool lacks pool_contract_version")
    pool["sector"] = pd.to_numeric(pool["sector"], errors="raise").astype(
        np.int64
    )
    pool["tic"] = _positive_tics(pool["tic"], context="frozen pool")
    if not pool["pool_contract_version"].astype(str).eq(
        FULL_POOL_CONTRACT_VERSION
    ).all():
        raise ValueError("frozen pool has the wrong contract")
    if tuple(sorted(pool["sector"].astype(int).unique())) != (
        FULL_POOL_NATIVE_SECTORS
    ):
        raise ValueError("frozen pool is not exact S56--S62 coverage")
    retained = summary.get("counts", {}).get("retained", {})
    if retained != {
        "n_observations": int(len(pool)),
        "n_unique_tics": int(pool["tic"].nunique()),
        "n_multisector_tics": int(
            pool.groupby("tic")["sector"].nunique().gt(1).sum()
        ),
    }:
        raise ValueError("frozen pool counts differ from summary")
    if summary.get("identity_hashes", {}).get(
        "retained_observations_sha256"
    ) != _observation_identity_sha256(pool):
        raise ValueError("frozen pool identity differs from summary")
    if summary.get("leakage_audit") != {
        "fixed_test_observations_retained": 0,
        "s63_reserved_observations_retained": 0,
    }:
        raise ValueError("frozen pool leakage audit did not pass")
    expected_keys = set(
        zip(pool["sector"].astype(int), pool["tic"].astype(int))
    )
    if not expected_keys or len(expected_keys) != len(pool):
        raise ValueError("frozen pool has empty or duplicate observation keys")

    records: list[dict[str, Any]] = []
    shard_metadata: list[dict[str, Any]] = []
    paths = [Path(value).expanduser().resolve(strict=True) for value in native_shard_paths]
    if not paths or len(paths) != len(set(paths)):
        raise ValueError("native shard paths must be nonempty and unique")
    import h5py

    for path in sorted(paths):
        audit = verify_full_pool_native_shard(path)
        if not audit["passed"]:
            raise ValueError(
                f"native shard failed verification {path}: "
                + "; ".join(audit["failures"][:10])
            )
        with h5py.File(path, "r") as h5:
            if h5.attrs.get("full_pool_sha256") != pool_binding.sha256:
                raise ValueError(f"native shard binds a different pool: {path}")
            if (
                h5.attrs.get("full_pool_summary_sha256")
                != summary_binding.sha256
            ):
                raise ValueError(
                    f"native shard binds a different pool summary: {path}"
                )
            sector = int(h5.attrs["sector"])
            for key, group in h5["targets"].items():
                tic = int(key)
                if (
                    int(group.attrs.get("sector", -1)) != sector
                    or int(group.attrs.get("tic", -1)) != tic
                ):
                    raise ValueError(f"native group identity mismatch in {path}")
                records.append(
                    {
                        "source_schema_version": (
                            FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION
                        ),
                        "sector": sector,
                        "tic": tic,
                        "native_h5_path": str(path),
                        "native_group_path": f"targets/{key}",
                        "full_pool_sha256": pool_binding.sha256,
                        "full_pool_summary_sha256": summary_binding.sha256,
                    }
                )
        shard_metadata.append(
            {
                "path": str(path),
                "sha256": audit["sha256"],
                "size_bytes": audit["size_bytes"],
                "sector": int(audit.get("sector", sector)),
                "n_observations": int(audit["counts"]["targets"]),
            }
        )
    source = pd.DataFrame.from_records(records)
    observed_keys = set(
        zip(source["sector"].astype(int), source["tic"].astype(int))
    )
    if source.duplicated(["sector", "tic"]).any() or observed_keys != expected_keys:
        missing = sorted(expected_keys - observed_keys)[:10]
        unexpected = sorted(observed_keys - expected_keys)[:10]
        raise ValueError(
            "native shards do not exactly cover the frozen pool; "
            f"missing={missing}, unexpected={unexpected}"
        )
    source = source.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )
    source_path = Path(source_path)
    source_payload = source.to_csv(
        index=False, lineterminator="\n"
    ).encode("utf-8")
    _publish_immutable_bytes(source_path, source_payload)
    registry_summary = write_native_input_registry(
        source_path=source_path,
        registry_path=registry_path,
        summary_path=summary_path,
        expected_contract_version=FULL_POOL_NATIVE_CONTRACT_VERSION,
    )
    pool_binding.assert_unchanged()
    summary_binding.assert_unchanged()
    return {
        **registry_summary,
        "full_pool_sha256": pool_binding.sha256,
        "full_pool_summary_sha256": summary_binding.sha256,
        "native_shards": shard_metadata,
        "coverage_identity_sha256": _observation_identity_sha256(source),
    }


__all__ = [
    "FULL_POOL_NATIVE_CONTRACT_VERSION",
    "FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION",
    "FULL_POOL_NATIVE_RELEASE_BINDING",
    "FULL_POOL_NATIVE_SECTORS",
    "FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION",
    "FULL_POOL_RAW_SOURCE_CONTRACT_VERSION",
    "SectorPoolAuthority",
    "build_full_pool_native_shard",
    "export_full_pool_raw_source_shard",
    "full_pool_native_group_failures",
    "full_pool_native_root_failures",
    "load_sector_pool_authority",
    "select_sector_shard",
    "shard_for_tic",
    "verify_full_pool_native_shard",
    "write_full_pool_native_registry",
]
