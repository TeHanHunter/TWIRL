"""Native raw/error preparation for the label-free S56--S62 SSL pool.

This module deliberately does not call ``build_raw_pair_export``.  That
function retains the immutable Teacher-v3 S62 release binding (997 rows and
exact release hashes).  The broad SSL pool needs the same model-facing
channels, but a different observation-keyed, real-only contract whose S62
orbit-ID reconciliation is validated per cadence rather than by a frozen
row-count signature.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings

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
from twirl.lightcurves.detrend_presets import adp03q_config
from twirl.lightcurves.flux_detrend import flux_space_detrend_result
from twirl.vetting import harmonic_export
from twirl.vetting.harmonic_inputs import (
    CHRONOLOGY_SMALL_CHANNELS,
    CHRONOLOGY_SUPPLEMENTAL_CHANNELS,
    HARMONIC_VIEW_CHANNELS,
    NATIVE_DATASETS,
    ORBITID_RECONCILIATION_MASK_DATASET,
    PERIODOGRAM_CHANNELS,
    PERIODOGRAM_DATASETS,
    RAW_PAIR_ORBITID_RECONCILIATION_ATTRS,
    RAW_PAIR_QUALITY_COUNT_NAMES,
    read_native_light_curve_from_h5,
)
from twirl.vetting.ssl_full_pool import (
    EXPECTED_SECTORS,
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
)
from twirl.vetting.ssl_full_pool_eligibility import (
    EligibilityAuthority,
    NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
    PRODUCTION_FROZEN_POOL_CSV_SHA256,
    load_native_model_eligibility,
    observation_identity_sha256,
)
from twirl.vetting.teacher_native_registry import (
    file_sha256,
    read_table,
    validate_native_input_registry_path,
    write_native_input_registry,
)


FULL_POOL_RAW_SOURCE_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_tglc_raw_source_v1"
)
FULL_POOL_NATIVE_CONTRACT_VERSION_V1 = (
    "twirl_teacher_ssl_fullpool_real_native_v1"
)
FULL_POOL_NATIVE_CONTRACT_VERSION_V2 = (
    "twirl_teacher_ssl_fullpool_real_native_v2"
)
FULL_POOL_NATIVE_CONTRACT_VERSION_V3 = (
    "twirl_teacher_ssl_fullpool_real_native_v3_effective_quality_adp"
)
FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1 = (
    "twirl_teacher_ssl_fullpool_real_native_v3r1_effective_quality_adp_btjd"
)
FULL_POOL_NATIVE_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_real_native_v4_detector_consistent_raw_v1_"
    "effective_quality_adp_btjd"
)
FULL_POOL_NATIVE_RELEASE_BINDING_V1 = (
    "teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp_real_only"
)
FULL_POOL_NATIVE_RELEASE_BINDING_V2 = (
    "teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible"
)
FULL_POOL_NATIVE_RELEASE_BINDING_V3 = (
    "teacher_ssl_fullpool_v3_s56_s62_a2v1_effective_quality_adp_bls_eligible"
)
FULL_POOL_NATIVE_RELEASE_BINDING_V3R1 = (
    "teacher_ssl_fullpool_v3r1_s56_s62_a2v1_effective_quality_adp_btjd_bls_eligible"
)
FULL_POOL_NATIVE_RELEASE_BINDING = (
    "teacher_ssl_fullpool_v4_s56_s62_a2v1_raw_v1_detector_consistent_"
    "effective_quality_adp_btjd_bls_eligible"
)
FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION_V2 = (
    "twirl_teacher_ssl_fullpool_native_shard_summary_v2"
)
FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION_V3 = (
    "twirl_teacher_ssl_fullpool_native_shard_summary_v3"
)
FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION_V4 = (
    "twirl_teacher_ssl_fullpool_native_shard_summary_v4"
)
FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_native_shard_summary_v5"
)
FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION_V2 = (
    "twirl_teacher_ssl_fullpool_native_registry_source_v2"
)
FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION_V3 = (
    "twirl_teacher_ssl_fullpool_native_registry_source_v3"
)
FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION_V4 = (
    "twirl_teacher_ssl_fullpool_native_registry_source_v4"
)
FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_native_registry_source_v5"
)
FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION_V2 = (
    "twirl_teacher_ssl_fullpool_native_registry_release_v2"
)
FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION_V3 = (
    "twirl_teacher_ssl_fullpool_native_registry_release_v3"
)
FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION_V4 = (
    "twirl_teacher_ssl_fullpool_native_registry_release_v4"
)
FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_native_registry_release_v5"
)
FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION_V1 = (
    "twirl_teacher_ssl_fullpool_effective_quality_adp_builder_v1"
)
FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION_V2 = (
    "twirl_teacher_ssl_fullpool_effective_quality_adp_btjd_builder_v2"
)
FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_detector_consistent_raw_v1_"
    "effective_quality_adp_btjd_builder_v3"
)
FULL_POOL_NATIVE_PERIODOGRAM_N = int(harmonic_export.DEFAULT_BLS_PERIODS)
FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION_V1 = (
    "twirl_fs_adp03q_effective_quality_v1"
)
FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION = (
    "twirl_fs_adp03q_effective_quality_btjd_v2"
)
FULL_POOL_NATIVE_DETREND_CONFIG_JSON = json.dumps(
    asdict(adp03q_config()),
    sort_keys=True,
    separators=(",", ":"),
    allow_nan=False,
)
FULL_POOL_NATIVE_DETREND_CONFIG_SHA256 = hashlib.sha256(
    FULL_POOL_NATIVE_DETREND_CONFIG_JSON.encode("ascii")
).hexdigest()
FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_adp_detrend_time_btjd_v1"
)
FULL_POOL_NATIVE_DETREND_TIME_DATASET = "adp_detrend_time_btjd"
FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D = 2_457_000.0
FULL_POOL_NATIVE_BTJD_MIN = 0.0
FULL_POOL_NATIVE_BTJD_MAX = 100_000.0
FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY = (
    "twirl_teacher_ssl_fullpool_exact_numpy_rankwarning_capture_v1"
)
FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY = "require_zero"
FULL_POOL_NATIVE_RANK_WARNING_MESSAGE = "Polyfit may be poorly conditioned"
_NUMPY_RANK_WARNING = getattr(np, "RankWarning", None)
if _NUMPY_RANK_WARNING is None:  # NumPy >= 2.0
    _NUMPY_RANK_WARNING = np.exceptions.RankWarning
FULL_POOL_NATIVE_RANK_WARNING_CATEGORY = (
    f"{_NUMPY_RANK_WARNING.__module__}.{_NUMPY_RANK_WARNING.__qualname__}"
)
FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON = "[]"
FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256 = hashlib.sha256(
    FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON.encode("ascii")
).hexdigest()
FULL_POOL_NATIVE_SECTORS = EXPECTED_SECTORS
FULL_POOL_NATIVE_SHARDS_PER_SECTOR = 16
FULL_POOL_RAW_EXPORT_CONTROLLER_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_pdo_raw_export_controller_v1"
)
FULL_POOL_RAW_TRANSFER_SCHEMA_VERSION = (
    "twirl_teacher_ssl_raw_transfer_validation_v1"
)
FULL_POOL_RAW_SHARD_SUMMARY_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_native_shard_summary_v1"
)
PRODUCTION_RAW_CODE_REVISION = (
    "991a94b9b36cdb07a4543b87602748c9bb2267f9"
)
PRODUCTION_RAW_EXPORT_COMPLETE_SHA256 = (
    "e61e94b623059e12081d3b5f6b76d8d7fe31bca035d85ca12f4360cae6112749"
)
PRODUCTION_RAW_EXPORT_LAUNCHER_SHA256 = (
    "204f4ad953495ae548786e554b9f2a00671b2aa966dd35bec797ca7936ddb4cf"
)
PRODUCTION_FULL_POOL_SUMMARY_SHA256 = (
    "f9bb5afa6672db0d74f4a7de6c0c6064295e4e0bfae420bf78777ec101b303c4"
)
PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR: Mapping[int, str] = {
    56: "f2d40348cc2cc5350370d103b703300c9ea7eda6753110f7e2b08b4707bbd07c",
    57: "25c38fc3613a4309aa663bac8058f5575eed74bff6ad33ff0dbefe24f42407a0",
    58: "7902be769b6b3b535545f621899dcee8b0d912964e142dd4a3badd7cfecdb8c5",
    59: "2c6b495ba25637a1a3718a8e91d39d1784cceb3b92f97193215a84f408a507a1",
    60: "3698b2942ffcff4a614a0ff82c29ed9673ac2663d800cd485b4301a3339a02c5",
    61: "505fdf42cb74efcdc343023797c374d58a42ee7a4e6a9728915107a1e3d10f89",
    62: "76491700efff698c9f7a7b5583eed27762d3bece35048edb0919ab37a25d77eb",
}
PRODUCTION_RAW_SECTOR_OBSERVATIONS: Mapping[int, int] = {
    56: 30_814,
    57: 26_777,
    58: 22_701,
    59: 21_646,
    60: 26_581,
    61: 27_002,
    62: 19_845,
}
PRODUCTION_RAW_SECTOR_IDENTITY_SHA256: Mapping[int, str] = {
    56: "8c80982649a0befb1a9b4e697008873fcaa4f3e37f62e0ec60e48f3fc1b93933",
    57: "ed5cace2aafcd87d53942937e5b49a9f92bc786f67ea26137d373c63d539bc4d",
    58: "00b9a11b1a8e328ea7d3b58f34d5e1344456801aaecf9e0a118999d0f98f57b4",
    59: "45385a3a1f6b3ffa0e292e8979107211b84d3ac4394cb283075961004215a838",
    60: "53adf7edd4ca1731683dfd604b7e2744fd08dd1a2302cd7b38ccf4b69625d6e7",
    61: "8df73671ea4f26a91d8f111dcd90143e8247c77e29156e0b8b532e6c2bcb6834",
    62: "ac856f3f194fb58448d503c8768a7e695e9be3a888411bf97b4d81f8453ca8ff",
}
PRODUCTION_RAW_SECTOR_ALLOWLIST_SHA256: Mapping[int, str] = {
    56: "dd95621e7a6ce09d145423821497127c03e1eace3fe4fba77bcb69060ae866eb",
    57: "2913dde75162152a5bce6fae4265ff93e733467dde0db3615dc2d315d69a3238",
    58: "c263751edd9f6df0b60fb094b72a4fee7849f6fd7e64c85f836030db73561b53",
    59: "1f0ae5e0a5541efde37ecff2887365fa935f973522069255ab5ec7e4824b0160",
    60: "e3943c1caba35ba68ea196854b0ef691053f9d1b3ccb5fbde5be80b6ac70e6e6",
    61: "6ccb5e83f484ba7327072af02a07cbe39cd18ddcbff1e0d73a3fdf45c6be31dd",
    62: "e9809a0c9150a7a48bef0319133011fa7dc3f360161b00ce5fa77d4be6e4760f",
}
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


@dataclass(frozen=True)
class RawSourceReleaseAuthority:
    """Exact validated v1 raw-source release consumed by one v2 shard."""

    raw_source: FileBinding
    raw_source_summary: FileBinding
    export_complete: FileBinding
    transfer_validation: FileBinding
    code_revision: str


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
    return observation_identity_sha256(rows)


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


def _read_json_binding(
    path: Path,
    *,
    context: str,
) -> tuple[dict[str, Any], FileBinding]:
    binding = _bind_file(path)
    try:
        value = json.loads(binding.path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid {context}: {binding.path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a JSON object")
    return value, binding


def load_production_raw_source_release(
    *,
    authority: SectorPoolAuthority,
    source_rows: pd.DataFrame,
    sector: int,
    shard_index: int,
    n_shards: int,
    raw_source_h5: Path,
    raw_source_summary_path: Path,
    raw_export_complete_path: Path,
    raw_transfer_validation_path: Path,
) -> RawSourceReleaseAuthority:
    """Authenticate one reused raw-v1 shard against the exact published release."""

    sector = int(sector)
    shard_index = int(shard_index)
    n_shards = int(n_shards)
    if (
        sector not in FULL_POOL_NATIVE_SECTORS
        or n_shards != FULL_POOL_NATIVE_SHARDS_PER_SECTOR
        or shard_index not in range(FULL_POOL_NATIVE_SHARDS_PER_SECTOR)
    ):
        raise ValueError("raw-source production release requires exact 7x16 layout")
    if (
        authority.sector != sector
        or authority.pool.sha256 != PRODUCTION_FROZEN_POOL_CSV_SHA256
        or authority.pool_summary.sha256
        != PRODUCTION_FULL_POOL_SUMMARY_SHA256
        or authority.allowlist.sha256
        != PRODUCTION_RAW_SECTOR_ALLOWLIST_SHA256[sector]
        or authority.observation_identity_sha256
        != PRODUCTION_RAW_SECTOR_IDENTITY_SHA256[sector]
        or len(authority.rows) != PRODUCTION_RAW_SECTOR_OBSERVATIONS[sector]
    ):
        raise ValueError("raw-source sector authority differs from production lock")

    aggregate, aggregate_binding = _read_json_binding(
        raw_export_complete_path,
        context="raw export completion manifest",
    )
    if aggregate_binding.sha256 != PRODUCTION_RAW_EXPORT_COMPLETE_SHA256:
        raise ValueError("raw export completion manifest differs from production")
    expected_aggregate = {
        "passed": True,
        "schema_version": FULL_POOL_RAW_EXPORT_CONTROLLER_SCHEMA_VERSION,
        "code_release": PRODUCTION_RAW_CODE_REVISION,
        "launcher_sha256": PRODUCTION_RAW_EXPORT_LAUNCHER_SHA256,
        "n_observations": sum(PRODUCTION_RAW_SECTOR_OBSERVATIONS.values()),
        "n_sectors": len(FULL_POOL_NATIVE_SECTORS),
        "n_shards": (
            len(FULL_POOL_NATIVE_SECTORS)
            * FULL_POOL_NATIVE_SHARDS_PER_SECTOR
        ),
        "n_shards_per_sector": FULL_POOL_NATIVE_SHARDS_PER_SECTOR,
        "full_pool_sha256": [PRODUCTION_FROZEN_POOL_CSV_SHA256],
        "full_pool_summary_sha256": [PRODUCTION_FULL_POOL_SUMMARY_SHA256],
        "sector_observation_identity_sha256": {
            str(key): value
            for key, value in PRODUCTION_RAW_SECTOR_IDENTITY_SHA256.items()
        },
    }
    for name, expected in expected_aggregate.items():
        if aggregate.get(name) != expected:
            raise ValueError(
                f"raw export completion manifest {name} differs from production"
            )

    transfer, transfer_binding = _read_json_binding(
        raw_transfer_validation_path,
        context=f"S{sector} raw transfer validation",
    )
    if (
        transfer_binding.sha256
        != PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR[sector]
    ):
        raise ValueError(
            f"S{sector} raw transfer validation differs from production"
        )
    expected_transfer = {
        "passed": True,
        "schema_version": FULL_POOL_RAW_TRANSFER_SCHEMA_VERSION,
        "sector": sector,
        "n_observations": PRODUCTION_RAW_SECTOR_OBSERVATIONS[sector],
        "n_shards": FULL_POOL_NATIVE_SHARDS_PER_SECTOR,
        "full_pool_sha256": PRODUCTION_FROZEN_POOL_CSV_SHA256,
        "full_pool_summary_sha256": PRODUCTION_FULL_POOL_SUMMARY_SHA256,
        "sector_allowlist_sha256": (
            PRODUCTION_RAW_SECTOR_ALLOWLIST_SHA256[sector]
        ),
        "sector_observation_identity_sha256": (
            PRODUCTION_RAW_SECTOR_IDENTITY_SHA256[sector]
        ),
    }
    for name, expected in expected_transfer.items():
        if transfer.get(name) != expected:
            raise ValueError(
                f"S{sector} raw transfer validation {name} differs"
            )
    shard_inventory = transfer.get("shards")
    if not isinstance(shard_inventory, list):
        raise ValueError("raw transfer validation lacks shard inventory")
    by_index: dict[int, Mapping[str, Any]] = {}
    for item in shard_inventory:
        if not isinstance(item, Mapping):
            raise ValueError("raw transfer shard inventory is malformed")
        index = int(item.get("shard_index", -1))
        if index in by_index:
            raise ValueError("raw transfer shard inventory has duplicate indices")
        by_index[index] = item
    if set(by_index) != set(range(FULL_POOL_NATIVE_SHARDS_PER_SECTOR)):
        raise ValueError("raw transfer shard inventory is incomplete")
    shard_meta = by_index[shard_index]

    raw_binding = _bind_file(
        raw_source_h5,
        expected_sha256=str(shard_meta.get("sha256", "")),
        expected_size_bytes=int(shard_meta.get("size_bytes", -1)),
    )
    expected_source_identity = _observation_identity_sha256(source_rows)
    if (
        int(shard_meta.get("n_observations", -1)) != len(source_rows)
        or shard_meta.get("shard_observation_identity_sha256")
        != expected_source_identity
    ):
        raise ValueError("raw transfer shard metadata differs from frozen pool")

    shard_summary, shard_summary_binding = _read_json_binding(
        raw_source_summary_path,
        context=f"S{sector} raw source shard {shard_index} summary",
    )
    expected_summary = {
        "schema_version": FULL_POOL_RAW_SHARD_SUMMARY_SCHEMA_VERSION,
        "stage": "pdo_raw_source_export",
        "sector": sector,
        "shard_index": shard_index,
        "n_shards": FULL_POOL_NATIVE_SHARDS_PER_SECTOR,
        "n_sector_observations": PRODUCTION_RAW_SECTOR_OBSERVATIONS[sector],
        "n_shard_observations": len(source_rows),
        "sector_observation_identity_sha256": (
            PRODUCTION_RAW_SECTOR_IDENTITY_SHA256[sector]
        ),
        "shard_observation_identity_sha256": expected_source_identity,
        "full_pool_sha256": PRODUCTION_FROZEN_POOL_CSV_SHA256,
        "full_pool_summary_sha256": PRODUCTION_FULL_POOL_SUMMARY_SHA256,
        "sector_allowlist_sha256": (
            PRODUCTION_RAW_SECTOR_ALLOWLIST_SHA256[sector]
        ),
        "out_h5_size_bytes": raw_binding.size_bytes,
        "out_h5_sha256": raw_binding.sha256,
        "real_only": True,
    }
    for name, expected in expected_summary.items():
        if shard_summary.get(name) != expected:
            raise ValueError(f"raw source shard summary {name} differs")

    for binding in (
        authority.pool,
        authority.pool_summary,
        authority.allowlist,
        raw_binding,
        shard_summary_binding,
        aggregate_binding,
        transfer_binding,
    ):
        binding.assert_unchanged()
    return RawSourceReleaseAuthority(
        raw_source=raw_binding,
        raw_source_summary=shard_summary_binding,
        export_complete=aggregate_binding,
        transfer_validation=transfer_binding,
        code_revision=PRODUCTION_RAW_CODE_REVISION,
    )


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
    raw_orbitid: np.ndarray,
    resolved_orbitid: np.ndarray,
    correction_mask: np.ndarray,
) -> None:
    """Validate S62 reconciliation from the raw-v1 orbit authority."""

    cadences = np.asarray(cadenceno, dtype=np.int64)
    raw = np.asarray(raw_orbitid, dtype=np.int64)
    resolved = np.asarray(resolved_orbitid, dtype=np.int64)
    corrected = np.asarray(correction_mask, dtype=bool)
    if len({len(cadences), len(raw), len(resolved), len(corrected)}) != 1:
        raise ValueError(f"TIC {tic}: orbit-ID reconciliation lengths differ")
    mismatch = raw != resolved
    if not np.array_equal(mismatch, corrected):
        raise ValueError(
            f"TIC {tic}: not all and only authority mismatches were corrected"
        )
    if np.any(raw[corrected] != S62_CORRECTION_INPUT_ORBIT) or np.any(
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


def _validated_btjd_time(
    values: np.ndarray,
    *,
    context: str,
) -> np.ndarray:
    """Return a bounded finite BTJD vector suitable for flux detrending."""

    time_btjd = np.asarray(values, dtype=np.float64)
    if time_btjd.ndim != 1 or time_btjd.size == 0:
        raise ValueError(f"{context}: detrend time must be a nonempty 1-D array")
    if not np.isfinite(time_btjd).all():
        raise ValueError(f"{context}: detrend BTJD contains nonfinite values")
    if (
        float(np.min(time_btjd)) < FULL_POOL_NATIVE_BTJD_MIN
        or float(np.max(time_btjd)) >= FULL_POOL_NATIVE_BTJD_MAX
    ):
        raise ValueError(f"{context}: detrend time is outside the BTJD bound")
    if np.any(np.diff(time_btjd) <= 0.0):
        raise ValueError(f"{context}: detrend BTJD is not strictly increasing")
    return time_btjd


def _published_bjd_from_btjd(
    time_btjd: np.ndarray,
    *,
    context: str,
) -> np.ndarray:
    """Convert validated BTJD to the exact absolute-BJD payload relation."""

    validated = _validated_btjd_time(time_btjd, context=context)
    published = validated + FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
    if (
        not np.isfinite(published).all()
        or not np.array_equal(
            published,
            validated + FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
        )
    ):
        raise ValueError(f"{context}: published BJD does not exactly encode BTJD")
    return published


def _rank_warning_ledger_json(
    entries: Sequence[Mapping[str, Any]],
) -> str:
    """Canonicalize the exact captured RankWarning audit ledger."""

    normalized = [
        {
            "aperture": str(entry["aperture"]),
            "category": str(entry["category"]),
            "message": str(entry["message"]),
            "sector": int(entry["sector"]),
            "tic": int(entry["tic"]),
        }
        for entry in entries
    ]
    normalized.sort(
        key=lambda entry: (
            entry["sector"],
            entry["tic"],
            entry["aperture"],
            entry["category"],
            entry["message"],
        )
    )
    return json.dumps(
        normalized,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _effective_quality_adp03q(
    *,
    time_btjd: np.ndarray,
    raw_flux: np.ndarray,
    raw_error: np.ndarray,
    quality: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Rebuild and independently recenter one ADP03q aperture.

    The compact input is intentionally absent from this interface: v3 uses it
    only to select the established cadence inventory.  Every photometric value
    is derived from the immutable raw-v1 flux/error arrays with the final
    effective-quality overlay applied to the spline fit.
    """

    time_btjd = _validated_btjd_time(
        time_btjd,
        context="effective-quality ADP",
    )
    raw_flux = np.asarray(raw_flux, dtype=np.float64)
    raw_error = np.asarray(raw_error, dtype=np.float64)
    quality = np.asarray(quality, dtype=np.int32)
    lengths = {len(time_btjd), len(raw_flux), len(raw_error), len(quality)}
    if len(lengths) != 1:
        raise ValueError("effective-quality ADP inputs differ in length")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = flux_space_detrend_result(
            time_btjd,
            raw_flux,
            quality=quality,
            flux_err=raw_error,
            cfg=adp03q_config(),
        )
    rank_warning_categories: list[str] = []
    rank_warning_messages: list[str] = []
    for warning in caught:
        category = warning.category
        message = str(warning.message)
        if (
            category is not _NUMPY_RANK_WARNING
            or message != FULL_POOL_NATIVE_RANK_WARNING_MESSAGE
        ):
            category_name = (
                f"{category.__module__}.{category.__qualname__}"
            )
            raise RuntimeError(
                "effective-quality ADP emitted an unauthorized warning: "
                f"{category_name}: {message}"
            )
        rank_warning_categories.append(FULL_POOL_NATIVE_RANK_WARNING_CATEGORY)
        rank_warning_messages.append(message)
    detrended = np.asarray(result.det_flux, dtype=np.float64)
    effective_good = (quality == 0) & np.isfinite(detrended)
    if not np.any(effective_good):
        raise ValueError("effective-quality ADP has no finite quality-zero values")
    center = float(np.median(detrended[effective_good]))
    if not np.isfinite(center):
        raise ValueError("effective-quality ADP recenter median is nonfinite")
    detrended = detrended - center + 1.0
    if not np.isfinite(detrended[effective_good]).all():
        raise ValueError("effective-quality ADP produced nonfinite quality-zero values")
    diagnostics = {
        "fit_count": int(result.fit_count),
        "n_segments": int(result.n_segments),
        "scale": float(result.scale),
        "scale_source": str(result.scale_source),
        "cotrend_status": str(result.cotrend_status),
        "recenter_median_before": center,
        "n_effective_good": int(np.count_nonzero(effective_good)),
        "n_finite_output": int(np.count_nonzero(np.isfinite(detrended))),
        "detrend_time_contract_version": (
            FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
        ),
        "detrend_time_system": "BTJD",
        "detrend_time_min": float(np.min(time_btjd)),
        "detrend_time_max": float(np.max(time_btjd)),
        "warning_capture_policy": FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY,
        "rank_warning_publication_policy": (
            FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
        ),
        "rank_warning_count": len(rank_warning_messages),
        "rank_warning_categories_json": json.dumps(
            rank_warning_categories,
            separators=(",", ":"),
        ),
        "rank_warning_messages_json": json.dumps(
            rank_warning_messages,
            separators=(",", ":"),
        ),
    }
    return detrended, diagnostics


def _write_effective_quality_adp_diagnostics(
    group: Any,
    *,
    aperture: str,
    diagnostics: Mapping[str, Any],
) -> None:
    prefix = f"effective_quality_adp_{aperture}"
    for name, value in diagnostics.items():
        group.attrs[f"{prefix}_{name}"] = value


def build_full_pool_native_shard(
    *,
    sector: int,
    pool_path: Path,
    pool_summary_path: Path,
    eligibility_exclusions_path: Path,
    eligibility_summary_path: Path,
    allowlist_path: Path,
    raw_source_h5: Path,
    raw_source_summary_path: Path | None = None,
    raw_export_complete_path: Path | None = None,
    raw_transfer_validation_path: Path | None = None,
    compact_adp_h5: Path,
    cadence_reference_table: Path,
    cadence_reference_manifest: Path,
    out_h5: Path,
    expected_git_sha: str,
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
    if int(n_periods) != FULL_POOL_NATIVE_PERIODOGRAM_N:
        raise ValueError(
            "native-v3 periodogram resolution must equal "
            f"{FULL_POOL_NATIVE_PERIODOGRAM_N}"
        )
    expected_git_sha = str(expected_git_sha).strip().lower()
    if len(expected_git_sha) != 40 or any(
        value not in "0123456789abcdef" for value in expected_git_sha
    ):
        raise ValueError("expected_git_sha must be an exact 40-character Git SHA")

    authority = load_sector_pool_authority(
        sector=sector,
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        allowlist_path=allowlist_path,
    )
    eligibility = load_native_model_eligibility(
        eligibility_exclusions_path,
        eligibility_summary_path,
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        # The preceding v4 eligibility stage independently rederives the
        # corrected BLS complement, then publishes this final authority only
        # when it satisfies the preregistered production partition.
        production_lock=True,
        rederive_from_bls=False,
    )
    source_rows = select_sector_shard(
        authority,
        shard_index=shard_index,
        n_shards=n_shards,
    )
    source_keys = set(
        zip(
            source_rows["sector"].astype(int),
            source_rows["tic"].astype(int),
        )
    )
    if not source_keys.issubset(eligibility.full_keys):
        raise ValueError("native source shard lies outside eligibility full pool")
    eligible_keys = source_keys & eligibility.eligible_keys
    excluded_keys = source_keys & eligibility.excluded_keys
    if eligible_keys | excluded_keys != source_keys or eligible_keys & excluded_keys:
        raise ValueError("native source shard is not exactly eligibility partitioned")
    rows = source_rows.loc[
        [
            (int(row.sector), int(row.tic)) in eligible_keys
            for row in source_rows.itertuples(index=False)
        ]
    ].reset_index(drop=True)
    excluded_rows = source_rows.loc[
        [
            (int(row.sector), int(row.tic)) in excluded_keys
            for row in source_rows.itertuples(index=False)
        ]
    ].reset_index(drop=True)
    if rows.empty:
        raise ValueError("eligibility-authorized native shard is empty")
    eligibility_exclusions_binding = eligibility.bindings["exclusions"]
    eligibility_summary_binding = eligibility.bindings["eligibility_summary"]
    raw_release = load_production_raw_source_release(
        authority=authority,
        source_rows=source_rows,
        sector=sector,
        shard_index=shard_index,
        n_shards=n_shards,
        raw_source_h5=raw_source_h5,
        raw_source_summary_path=raw_source_summary_path,
        raw_export_complete_path=raw_export_complete_path,
        raw_transfer_validation_path=raw_transfer_validation_path,
    )
    raw_binding = raw_release.raw_source
    compact_binding = _bind_file(
        compact_adp_h5,
        expected_sha256=authority.compact_h5_sha256,
        expected_size_bytes=authority.compact_h5_size_bytes,
    )
    reference_table_binding = _bind_file(cadence_reference_table)
    reference_manifest_binding = _bind_file(cadence_reference_manifest)
    builder_binding = _bind_file(Path(__file__))
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
    rank_warning_ledger: list[dict[str, Any]] = []
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
                "n_shard_observations": int(len(source_rows)),
                "full_pool_sha256": authority.pool.sha256,
                "full_pool_summary_sha256": authority.pool_summary.sha256,
                "sector_allowlist_sha256": authority.allowlist.sha256,
                "sector_observation_identity_sha256": (
                    authority.observation_identity_sha256
                ),
                "shard_observation_identity_sha256": (
                    _observation_identity_sha256(source_rows)
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
            expected_tics = source_rows["tic"].astype(int).tolist()
            if raw_tics != expected_tics:
                raise ValueError("raw-source TIC inventory differs from shard")

            output.attrs["contract_version"] = (
                FULL_POOL_NATIVE_CONTRACT_VERSION
            )
            output.attrs["release_binding"] = FULL_POOL_NATIVE_RELEASE_BINDING
            output.attrs["builder_contract_version"] = (
                FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
            )
            output.attrs["builder_code_sha256"] = builder_binding.sha256
            output.attrs["expected_git_sha"] = expected_git_sha
            output.attrs["detrend_contract_version"] = (
                FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
            )
            output.attrs["detrend_config_json"] = (
                FULL_POOL_NATIVE_DETREND_CONFIG_JSON
            )
            output.attrs["detrend_config_sha256"] = (
                FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
            )
            output.attrs["detrend_quality_source"] = (
                "final_effective_quality"
            )
            output.attrs["detrend_time_contract_version"] = (
                FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
            )
            output.attrs["detrend_time_dataset"] = (
                FULL_POOL_NATIVE_DETREND_TIME_DATASET
            )
            output.attrs["detrend_time_system"] = "BTJD"
            output.attrs["published_time_system"] = "BJD"
            output.attrs["btjd_to_bjd_offset_d"] = (
                FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
            )
            output.attrs["detrend_btjd_min_inclusive"] = (
                FULL_POOL_NATIVE_BTJD_MIN
            )
            output.attrs["detrend_btjd_max_exclusive"] = (
                FULL_POOL_NATIVE_BTJD_MAX
            )
            output.attrs["warning_capture_policy"] = (
                FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY
            )
            output.attrs["rank_warning_publication_policy"] = (
                FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
            )
            output.attrs["raw_photometry_only"] = 1
            output.attrs["compact_adp_photometry_reused"] = 0
            output.attrs["compact_adp_flux_reused"] = 0
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
            output.attrs["n_source_shard_observations"] = int(
                len(source_rows)
            )
            output.attrs["n_shard_excluded_observations"] = int(
                len(excluded_rows)
            )
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
            output.attrs["source_shard_observation_identity_sha256"] = (
                _observation_identity_sha256(source_rows)
            )
            output.attrs["shard_excluded_observation_identity_sha256"] = (
                _observation_identity_sha256(excluded_rows)
            )
            output.attrs["eligibility_contract_version"] = (
                eligibility.contract_version
            )
            output.attrs["eligibility_exclusions_sha256"] = (
                eligibility_exclusions_binding.sha256
            )
            output.attrs["eligibility_summary_sha256"] = (
                eligibility_summary_binding.sha256
            )
            output.attrs["native_model_full_identity_sha256"] = (
                eligibility.full_observation_identity_sha256
            )
            output.attrs["native_model_eligible_identity_sha256"] = (
                eligibility.eligible_observation_identity_sha256
            )
            output.attrs["native_model_excluded_identity_sha256"] = (
                eligibility.excluded_observation_identity_sha256
            )
            output.attrs["raw_source_h5"] = str(raw_binding.path)
            output.attrs["raw_source_h5_sha256"] = raw_binding.sha256
            output.attrs["raw_source_release_code_revision"] = (
                raw_release.code_revision
            )
            output.attrs["raw_source_summary_sha256"] = (
                raw_release.raw_source_summary.sha256
            )
            output.attrs["raw_export_complete_sha256"] = (
                raw_release.export_complete.sha256
            )
            output.attrs["raw_transfer_validation_sha256"] = (
                raw_release.transfer_validation.sha256
            )
            output.attrs["compact_adp_h5"] = str(compact_binding.path)
            output.attrs["compact_adp_h5_sha256"] = compact_binding.sha256
            output.attrs["compact_adp_role"] = (
                "cadence_time_inventory_reference_only"
            )
            output.attrs["raw_detector_mapping_authority"] = 1
            output.attrs["raw_orbitid_authority"] = 1
            output.attrs["raw_internal_quality_authority"] = 1
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
            excluded_structurally_validated = 0
            for count, row in enumerate(
                source_rows.itertuples(index=False), start=1
            ):
                tic = int(row.tic)
                include_native = (sector, tic) in eligible_keys
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
                        if int(raw_group.attrs.get(name, -1)) != expected:
                            raise ValueError(
                                f"raw-source {name} mapping differs from pool"
                            )
                    raw = harmonic_export._raw_source_payload(raw_group)
                    required_raw = {
                        "time",
                        "cadenceno",
                        "orbitid",
                        "quality",
                        "raw_flux_small",
                        "raw_flux_err_small",
                        "raw_flux_primary",
                        "raw_flux_err_primary",
                    }
                    missing_raw = sorted(required_raw - set(raw))
                    if missing_raw:
                        raise ValueError(
                            f"raw-source group lacks arrays: {missing_raw}"
                        )
                    lengths = {name: len(raw[name]) for name in required_raw}
                    if len(set(lengths.values())) != 1 or next(
                        iter(lengths.values()), 0
                    ) < 1:
                        raise ValueError(
                            f"raw-source arrays differ in length: {lengths}"
                        )
                    raw_cadenceno = np.asarray(raw["cadenceno"], dtype=np.int64)
                    if (
                        len(np.unique(raw_cadenceno)) != len(raw_cadenceno)
                        or np.any(np.diff(raw_cadenceno) <= 0)
                    ):
                        raise ValueError(
                            "raw-source cadences are not unique and strictly sorted"
                        )
                    compact_cadenceno = np.asarray(
                        compact_group["cadenceno"], dtype=np.int64
                    )
                    if not np.array_equal(raw_cadenceno, compact_cadenceno):
                        raise ValueError(
                            "raw-source cadence inventory differs from frozen compact"
                        )
                    raw_time_bjd = np.asarray(raw["time"], dtype=np.float64)
                    if (
                        not np.isfinite(raw_time_bjd).all()
                        or np.any(raw_time_bjd < 2_000_000.0)
                    ):
                        raise ValueError(
                            "raw-source time is not finite absolute BJD"
                        )
                    raw_time_btjd = _validated_btjd_time(
                        raw_time_bjd - FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
                        context=f"S{sector} TIC {tic}",
                    )
                    compact_time = np.asarray(
                        compact_group["time"], dtype=np.float64
                    )
                    if len(compact_time) != len(raw_time_bjd):
                        raise ValueError("frozen compact time length differs")
                    compact_time_bjd = np.where(
                        compact_time < 100_000.0,
                        compact_time + FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
                        compact_time,
                    )
                    raw_compact_time_delta_s = np.abs(
                        raw_time_bjd - compact_time_bjd
                    ) * 86_400.0
                    if (
                        not np.isfinite(raw_compact_time_delta_s).all()
                        or float(np.max(raw_compact_time_delta_s)) > 2.0
                    ):
                        raise ValueError(
                            "raw-source/frozen-compact times differ by more than 2 s"
                        )
                    raw_orbitid = np.asarray(raw["orbitid"], dtype=np.int64)
                    raw_internal_quality = np.asarray(
                        raw["quality"], dtype=np.int64
                    )
                    quality_overlay = quality_reference.apply(
                        sector=sector,
                        camera=int(row.camera),
                        ccd=int(row.ccd),
                        cadenceno=raw_cadenceno,
                        orbitid=raw_orbitid,
                        internal_quality=raw_internal_quality,
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
                            cadenceno=raw_cadenceno,
                            raw_orbitid=raw_orbitid,
                            resolved_orbitid=resolved_orbitid,
                            correction_mask=correction_mask,
                        )
                    elif correction_mask.any() or not np.array_equal(
                        resolved_orbitid, raw_orbitid
                    ):
                        raise ValueError(
                            "strict orbit-ID policy changed one or more cadences"
                        )

                    if not include_native:
                        excluded_structurally_validated += 1
                        continue

                    corrected = harmonic_export._record_orbitid_reconciliation(
                        orbitid_stats,
                        camera=int(row.camera),
                        ccd=int(row.ccd),
                        cadenceno=raw_cadenceno,
                        input_orbitid=raw_orbitid,
                        resolved_orbitid=resolved_orbitid,
                        correction_mask=correction_mask,
                    )
                    group = targets.create_group(f"{tic:016d}")
                    harmonic_export._copy_attrs(raw_group, group)
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
                    group.attrs["raw_detector_mapping_authority"] = 1
                    group.attrs["raw_orbitid_authority"] = 1
                    group.attrs["raw_internal_quality_authority"] = 1
                    group.attrs["compact_cadence_time_inventory_match"] = 1
                    group.attrs["compact_cadence_time_inventory_role"] = (
                        "reference_only"
                    )
                    group.attrs["photometry_source"] = (
                        "immutable_raw_v1_effective_quality_adp03q"
                    )
                    group.attrs["compact_adp_photometry_reused"] = 0
                    group.attrs["compact_adp_flux_reused"] = 0
                    group.attrs["detrend_contract_version"] = (
                        FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
                    )
                    group.attrs["detrend_config_sha256"] = (
                        FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
                    )
                    group.attrs["detrend_quality_source"] = (
                        "final_effective_quality"
                    )
                    group.attrs["detrend_time_contract_version"] = (
                        FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
                    )
                    group.attrs["detrend_time_dataset"] = (
                        FULL_POOL_NATIVE_DETREND_TIME_DATASET
                    )
                    group.attrs["detrend_time_system"] = "BTJD"
                    group.attrs["published_time_system"] = "BJD"
                    group.attrs["btjd_to_bjd_offset_d"] = (
                        FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
                    )
                    group.attrs["raw_source_paths"] = raw_group.attrs.get(
                        "source_paths", ""
                    )
                    group.attrs["raw_compact_time_delta_max_s"] = float(
                        np.max(raw_compact_time_delta_s)
                    )
                    output_time = _published_bjd_from_btjd(
                        raw_time_btjd,
                        context=f"S{sector} TIC {tic}",
                    )
                    group.attrs["detrend_time_min"] = float(
                        np.min(raw_time_btjd)
                    )
                    group.attrs["detrend_time_max"] = float(
                        np.max(raw_time_btjd)
                    )
                    effective_quality = np.asarray(
                        quality_overlay.quality, dtype=np.int32
                    )
                    det_small, small_diagnostics = (
                        _effective_quality_adp03q(
                            time_btjd=raw_time_btjd,
                            raw_flux=raw["raw_flux_small"],
                            raw_error=raw["raw_flux_err_small"],
                            quality=effective_quality,
                        )
                    )
                    det_primary, primary_diagnostics = (
                        _effective_quality_adp03q(
                            time_btjd=raw_time_btjd,
                            raw_flux=raw["raw_flux_primary"],
                            raw_error=raw["raw_flux_err_primary"],
                            quality=effective_quality,
                        )
                    )
                    _write_effective_quality_adp_diagnostics(
                        group,
                        aperture="small",
                        diagnostics=small_diagnostics,
                    )
                    _write_effective_quality_adp_diagnostics(
                        group,
                        aperture="primary",
                        diagnostics=primary_diagnostics,
                    )
                    for aperture, diagnostics in (
                        ("small", small_diagnostics),
                        ("primary", primary_diagnostics),
                    ):
                        categories = json.loads(
                            str(
                                diagnostics[
                                    "rank_warning_categories_json"
                                ]
                            )
                        )
                        messages = json.loads(
                            str(
                                diagnostics[
                                    "rank_warning_messages_json"
                                ]
                            )
                        )
                        if len(categories) != len(messages) or len(
                            messages
                        ) != int(diagnostics["rank_warning_count"]):
                            raise RuntimeError(
                                "effective-quality ADP warning diagnostics "
                                "are internally inconsistent"
                            )
                        rank_warning_ledger.extend(
                            {
                                "sector": sector,
                                "tic": tic,
                                "aperture": aperture,
                                "category": category,
                                "message": message,
                            }
                            for category, message in zip(
                                categories,
                                messages,
                                strict=True,
                            )
                        )
                    payload = {
                        "time": output_time,
                        FULL_POOL_NATIVE_DETREND_TIME_DATASET: (
                            raw_time_btjd
                        ),
                        "cadenceno": raw_cadenceno,
                        "orbitid": resolved_orbitid.astype(np.int32),
                        "orbitid_pre_reference": raw_orbitid.astype(
                            np.int32
                        ),
                        "raw_orbitid_pre_reference": raw_orbitid.astype(
                            np.int32
                        ),
                        ORBITID_RECONCILIATION_MASK_DATASET: correction_mask.astype(
                            np.uint8
                        ),
                        "quality": effective_quality,
                        "raw_flux_small": raw["raw_flux_small"],
                        "raw_flux_err_small": raw["raw_flux_err_small"],
                        "raw_flux_primary": raw["raw_flux_primary"],
                        "raw_flux_err_primary": raw[
                            "raw_flux_err_primary"
                        ],
                        "det_flux_adp_sml": det_small,
                        "det_flux_adp": det_primary,
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
                        f"{shard_index}/{n_shards}: "
                        f"{count}/{len(source_rows)}",
                        flush=True,
                    )
            if excluded_structurally_validated != len(excluded_rows):
                raise RuntimeError(
                    "excluded-row structural validation count differs from "
                    "eligibility partition"
                )
            output.attrs[
                "n_shard_excluded_structurally_validated"
            ] = int(excluded_structurally_validated)
            rank_warning_ledger_json = _rank_warning_ledger_json(
                rank_warning_ledger
            )
            output.attrs["rank_warning_count"] = int(
                len(rank_warning_ledger)
            )
            output.attrs["rank_warning_ledger_json"] = (
                rank_warning_ledger_json
            )
            output.attrs["rank_warning_ledger_sha256"] = hashlib.sha256(
                rank_warning_ledger_json.encode("ascii")
            ).hexdigest()
            if rank_warning_ledger:
                failures.append(
                    {
                        "sector": sector,
                        "tic": -1,
                        "error": (
                            "native publication requires zero exact "
                            f"RankWarnings; observed={len(rank_warning_ledger)}"
                        ),
                    }
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
                "native export failed for "
                f"{len(failures)} of {len(source_rows)} TICs"
            )
        for binding in (
            authority.pool,
            authority.pool_summary,
            authority.allowlist,
            raw_binding,
            compact_binding,
            reference_table_binding,
            reference_manifest_binding,
            builder_binding,
            raw_release.raw_source_summary,
            raw_release.export_complete,
            raw_release.transfer_validation,
        ):
            binding.assert_unchanged()
        for binding in (
            eligibility_exclusions_binding,
            eligibility_summary_binding,
        ):
            _bind_file(
                binding.path,
                expected_sha256=binding.sha256,
                expected_size_bytes=binding.size_bytes,
            )
        quality_reference.assert_unchanged()
        verification = verify_full_pool_native_shard(
            temporary,
            expected_git_sha=expected_git_sha,
            expected_sector=sector,
            expected_shard_index=int(shard_index),
            expected_n_shards=int(n_shards),
            expected_observations=len(rows),
            expected_source_observations=len(source_rows),
            expected_excluded_observations=len(excluded_rows),
            expected_source_identity_sha256=_observation_identity_sha256(
                source_rows
            ),
            expected_output_identity_sha256=_observation_identity_sha256(
                rows
            ),
            expected_excluded_identity_sha256=_observation_identity_sha256(
                excluded_rows
            ),
            require_periodograms=True,
        )
        if not verification["passed"]:
            raise RuntimeError(
                "temporary full-pool native shard failed verification: "
                + "; ".join(verification["failures"][:10])
            )
        # Publish without a clobber window.  ``Path.replace`` could overwrite
        # an output created after the initial existence check; the hard link
        # is atomic and fails closed if another writer won that race.
        os.link(temporary, out_h5)
        temporary.unlink()
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    failures_path = out_h5.with_suffix(".failures.csv")
    failures_path.unlink(missing_ok=True)
    verification = {
        **verification,
        "path": str(out_h5),
    }
    return {
        "schema_version": FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION,
        "stage": "orcd_native_shard_build",
        "sector": sector,
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
        "n_sector_observations": int(len(authority.rows)),
        "n_source_shard_observations": int(len(source_rows)),
        "n_shard_observations": int(len(rows)),
        "n_shard_excluded_observations": int(len(excluded_rows)),
        "n_shard_excluded_structurally_validated": int(
            len(excluded_rows)
        ),
        "sector_observation_identity_sha256": (
            authority.observation_identity_sha256
        ),
        "source_shard_observation_identity_sha256": (
            _observation_identity_sha256(source_rows)
        ),
        "shard_observation_identity_sha256": (
            _observation_identity_sha256(rows)
        ),
        "shard_excluded_observation_identity_sha256": (
            _observation_identity_sha256(excluded_rows)
        ),
        "full_pool_sha256": authority.pool.sha256,
        "full_pool_summary_sha256": authority.pool_summary.sha256,
        "sector_allowlist_sha256": authority.allowlist.sha256,
        "raw_source_h5_sha256": raw_binding.sha256,
        "raw_source_release_code_revision": raw_release.code_revision,
        "raw_source_summary_sha256": raw_release.raw_source_summary.sha256,
        "raw_export_complete_sha256": raw_release.export_complete.sha256,
        "raw_transfer_validation_sha256": (
            raw_release.transfer_validation.sha256
        ),
        "compact_adp_h5_sha256": compact_binding.sha256,
        "compact_adp_role": "cadence_time_inventory_reference_only",
        "raw_detector_mapping_authority": True,
        "raw_orbitid_authority": True,
        "raw_internal_quality_authority": True,
        "cadence_reference_table_sha256": reference_table_binding.sha256,
        "cadence_reference_manifest_sha256": (
            reference_manifest_binding.sha256
        ),
        "expected_git_sha": expected_git_sha,
        "builder_contract_version": (
            FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
        ),
        "builder_code_sha256": builder_binding.sha256,
        "detrend_contract_version": (
            FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
        ),
        "detrend_config_sha256": (
            FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
        ),
        "detrend_config_json": FULL_POOL_NATIVE_DETREND_CONFIG_JSON,
        "detrend_quality_source": "final_effective_quality",
        "detrend_time_contract_version": (
            FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
        ),
        "detrend_time_dataset": FULL_POOL_NATIVE_DETREND_TIME_DATASET,
        "detrend_time_system": "BTJD",
        "published_time_system": "BJD",
        "btjd_to_bjd_offset_d": FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
        "warning_capture_policy": FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY,
        "rank_warning_publication_policy": (
            FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
        ),
        "rank_warning_count": len(rank_warning_ledger),
        "rank_warning_ledger_json": _rank_warning_ledger_json(
            rank_warning_ledger
        ),
        "rank_warning_ledger_sha256": hashlib.sha256(
            _rank_warning_ledger_json(rank_warning_ledger).encode("ascii")
        ).hexdigest(),
        "raw_photometry_only": True,
        "compact_adp_photometry_reused": False,
        "compact_adp_flux_reused": False,
        "periodogram_n": FULL_POOL_NATIVE_PERIODOGRAM_N,
        "eligibility_contract_version": eligibility.contract_version,
        "eligibility_exclusions_sha256": (
            eligibility_exclusions_binding.sha256
        ),
        "eligibility_summary_sha256": eligibility_summary_binding.sha256,
        "native_model_full_identity_sha256": (
            eligibility.full_observation_identity_sha256
        ),
        "native_model_eligible_identity_sha256": (
            eligibility.eligible_observation_identity_sha256
        ),
        "native_model_excluded_identity_sha256": (
            eligibility.excluded_observation_identity_sha256
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
    contract = str(attrs.get("contract_version", ""))
    if contract not in {
        FULL_POOL_NATIVE_CONTRACT_VERSION_V1,
        FULL_POOL_NATIVE_CONTRACT_VERSION_V2,
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3,
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1,
        FULL_POOL_NATIVE_CONTRACT_VERSION,
    }:
        return ["wrong full-pool native contract"]
    is_eligibility_partitioned = contract in {
        FULL_POOL_NATIVE_CONTRACT_VERSION_V2,
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3,
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1,
        FULL_POOL_NATIVE_CONTRACT_VERSION,
    }
    is_effective_quality_adp = contract in {
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3,
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1,
        FULL_POOL_NATIVE_CONTRACT_VERSION,
    }
    is_btjd = contract in {
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1,
        FULL_POOL_NATIVE_CONTRACT_VERSION,
    }
    is_detector_consistent_v4 = contract == FULL_POOL_NATIVE_CONTRACT_VERSION
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
    if is_eligibility_partitioned:
        required.update(
            {
                "eligibility_contract_version",
                "eligibility_exclusions_sha256",
                "eligibility_summary_sha256",
                "native_model_full_identity_sha256",
                "native_model_eligible_identity_sha256",
                "native_model_excluded_identity_sha256",
                "n_source_shard_observations",
                "n_shard_excluded_observations",
                "source_shard_observation_identity_sha256",
                "shard_excluded_observation_identity_sha256",
                "raw_source_release_code_revision",
                "raw_source_summary_sha256",
                "raw_export_complete_sha256",
                "raw_transfer_validation_sha256",
            }
        )
    if is_effective_quality_adp:
        required.update(
            {
                "release_binding",
                "builder_contract_version",
                "builder_code_sha256",
                "expected_git_sha",
                "detrend_contract_version",
                "detrend_config_json",
                "detrend_config_sha256",
                "periodogram_n",
                "detrend_quality_source",
                "raw_photometry_only",
                "compact_adp_photometry_reused",
                "compact_adp_flux_reused",
                "compact_adp_role",
                "n_shard_excluded_structurally_validated",
            }
        )
        try:
            if int(attrs["periodogram_n"]) != FULL_POOL_NATIVE_PERIODOGRAM_N:
                failures.append("periodogram_n differs from native-v3 contract")
        except (KeyError, TypeError, ValueError):
            failures.append("periodogram_n is invalid")
    if is_btjd:
        required.update(
            {
                "detrend_time_contract_version",
                "detrend_time_dataset",
                "detrend_time_system",
                "published_time_system",
                "btjd_to_bjd_offset_d",
                "detrend_btjd_min_inclusive",
                "detrend_btjd_max_exclusive",
                "warning_capture_policy",
                "rank_warning_publication_policy",
                "rank_warning_count",
                "rank_warning_ledger_json",
                "rank_warning_ledger_sha256",
            }
        )
    if is_detector_consistent_v4:
        required.update(
            {
                "raw_detector_mapping_authority",
                "raw_orbitid_authority",
                "raw_internal_quality_authority",
            }
        )
    missing = sorted(name for name in required if name not in attrs)
    if missing:
        return [f"full-pool native root lacks attrs: {','.join(missing)}"]
    try:
        sector = int(attrs["sector"])
        shard_index = int(attrs["shard_index"])
        n_shards = int(attrs["n_shards"])
        n_sector = int(attrs["n_sector_observations"])
        n_shard = int(attrs["n_shard_observations"])
        n_source_shard = (
            int(attrs["n_source_shard_observations"])
            if is_eligibility_partitioned
            else n_shard
        )
        n_excluded_shard = (
            int(attrs["n_shard_excluded_observations"])
            if is_eligibility_partitioned
            else 0
        )
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
        or n_source_shard < 1
        or n_excluded_shard < 0
        or n_source_shard > n_sector
        or n_shard + n_excluded_shard != n_source_shard
        or real_only != 1
    ):
        failures.append("full-pool native root metadata is outside bounds")
    sha_names = [
        "full_pool_sha256",
        "full_pool_summary_sha256",
        "sector_allowlist_sha256",
        "sector_observation_identity_sha256",
        "shard_observation_identity_sha256",
        "raw_source_h5_sha256",
        "compact_adp_h5_sha256",
    ]
    if is_eligibility_partitioned:
        sha_names.extend(
            [
                "eligibility_exclusions_sha256",
                "eligibility_summary_sha256",
                "native_model_full_identity_sha256",
                "native_model_eligible_identity_sha256",
                "native_model_excluded_identity_sha256",
                "source_shard_observation_identity_sha256",
                "shard_excluded_observation_identity_sha256",
                "raw_source_summary_sha256",
                "raw_export_complete_sha256",
                "raw_transfer_validation_sha256",
            ]
        )
        if str(attrs["eligibility_contract_version"]) != (
            NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION
        ):
            failures.append("native/model eligibility contract mismatch")
        if str(attrs["raw_source_release_code_revision"]) != (
            PRODUCTION_RAW_CODE_REVISION
        ):
            failures.append("raw-source code revision differs from production")
        if str(attrs["raw_export_complete_sha256"]) != (
            PRODUCTION_RAW_EXPORT_COMPLETE_SHA256
        ):
            failures.append("raw export release differs from production")
        if str(attrs["raw_transfer_validation_sha256"]) != (
            PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR.get(sector, "")
        ):
            failures.append("raw transfer release differs from production")
    if is_effective_quality_adp:
        sha_names.extend(
            [
                "builder_code_sha256",
                "detrend_config_sha256",
            ]
        )
        expected_release_binding = (
            FULL_POOL_NATIVE_RELEASE_BINDING
            if is_detector_consistent_v4
            else FULL_POOL_NATIVE_RELEASE_BINDING_V3R1
            if contract == FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1
            else FULL_POOL_NATIVE_RELEASE_BINDING_V3
        )
        expected_builder_contract = (
            FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
            if is_detector_consistent_v4
            else FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION_V2
            if contract == FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1
            else FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION_V1
        )
        expected_detrend_contract = (
            FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
            if is_btjd
            else FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION_V1
        )
        if str(attrs["release_binding"]) != expected_release_binding:
            failures.append("v3 release binding mismatch")
        if str(attrs["builder_contract_version"]) != (
            expected_builder_contract
        ):
            failures.append("v3 builder contract mismatch")
        if str(attrs["detrend_contract_version"]) != (
            expected_detrend_contract
        ):
            failures.append("v3 detrend contract mismatch")
        if str(attrs["detrend_config_json"]) != (
            FULL_POOL_NATIVE_DETREND_CONFIG_JSON
        ):
            failures.append("v3 detrend config JSON mismatch")
        if str(attrs["detrend_config_sha256"]) != (
            FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
        ):
            failures.append("v3 detrend config SHA-256 mismatch")
        if str(attrs["detrend_quality_source"]) != "final_effective_quality":
            failures.append("v3 detrend did not use final effective quality")
        expected_compact_role = (
            "cadence_time_inventory_reference_only"
            if is_detector_consistent_v4
            else "cadence_time_orbit_detector_internal_quality_mapping_only"
        )
        if str(attrs["compact_adp_role"]) != expected_compact_role:
            failures.append("v3 compact input role mismatch")
        if is_detector_consistent_v4:
            try:
                authority_flags = {
                    name: int(attrs[name])
                    for name in (
                        "raw_detector_mapping_authority",
                        "raw_orbitid_authority",
                        "raw_internal_quality_authority",
                    )
                }
            except (KeyError, TypeError, ValueError):
                failures.append("native-v4 raw authority flags are invalid")
            else:
                if any(value != 1 for value in authority_flags.values()):
                    failures.append("native-v4 raw authority flags differ")
        try:
            v3_integers = {
                name: int(attrs[name])
                for name in (
                    "raw_photometry_only",
                    "compact_adp_photometry_reused",
                    "compact_adp_flux_reused",
                    "n_shard_excluded_structurally_validated",
                )
            }
        except (TypeError, ValueError):
            failures.append("v3 photometry provenance integers are invalid")
        else:
            if (
                v3_integers["raw_photometry_only"] != 1
                or v3_integers["compact_adp_photometry_reused"] != 0
                or v3_integers["compact_adp_flux_reused"] != 0
            ):
                failures.append("v3 photometry provenance permits compact reuse")
            if (
                v3_integers["n_shard_excluded_structurally_validated"]
                != n_excluded_shard
            ):
                failures.append(
                    "v3 excluded structural validation count mismatch"
                )
        expected_git_sha = str(attrs["expected_git_sha"])
        if len(expected_git_sha) != 40 or any(
            value not in "0123456789abcdef" for value in expected_git_sha
        ):
            failures.append("expected_git_sha is not an exact Git SHA")
    if is_btjd:
        sha_names.append("rank_warning_ledger_sha256")
        expected_time_attrs: dict[str, Any] = {
            "detrend_time_contract_version": (
                FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
            ),
            "detrend_time_dataset": FULL_POOL_NATIVE_DETREND_TIME_DATASET,
            "detrend_time_system": "BTJD",
            "published_time_system": "BJD",
            "warning_capture_policy": (
                FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY
            ),
            "rank_warning_publication_policy": (
                FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
            ),
            "rank_warning_ledger_json": (
                FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON
            ),
            "rank_warning_ledger_sha256": (
                FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256
            ),
        }
        for name, expected_value in expected_time_attrs.items():
            if str(attrs[name]) != expected_value:
                failures.append(f"{name} differs from native-v3r1 contract")
        try:
            if (
                float(attrs["btjd_to_bjd_offset_d"])
                != FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
                or float(attrs["detrend_btjd_min_inclusive"])
                != FULL_POOL_NATIVE_BTJD_MIN
                or float(attrs["detrend_btjd_max_exclusive"])
                != FULL_POOL_NATIVE_BTJD_MAX
            ):
                failures.append("native-v3r1 time bounds differ")
            if int(attrs["rank_warning_count"]) != 0:
                failures.append("native-v3r1 publication contains RankWarnings")
        except (TypeError, ValueError):
            failures.append("native-v3r1 time/warning metadata is invalid")
    for name in sha_names:
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
                if is_detector_consistent_v4
                else FULL_POOL_NATIVE_RELEASE_BINDING_V3R1
                if contract == FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1
                else FULL_POOL_NATIVE_RELEASE_BINDING_V3
                if contract == FULL_POOL_NATIVE_CONTRACT_VERSION_V3
                else FULL_POOL_NATIVE_RELEASE_BINDING_V2
                if contract == FULL_POOL_NATIVE_CONTRACT_VERSION_V2
                else FULL_POOL_NATIVE_RELEASE_BINDING_V1
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
    root_contract: str,
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
    if root_contract in {
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3,
        FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1,
        FULL_POOL_NATIVE_CONTRACT_VERSION,
    }:
        is_btjd = root_contract in {
            FULL_POOL_NATIVE_CONTRACT_VERSION_V3R1,
            FULL_POOL_NATIVE_CONTRACT_VERSION,
        }
        is_detector_consistent_v4 = (
            root_contract == FULL_POOL_NATIVE_CONTRACT_VERSION
        )
        v3_required_attrs = {
            "photometry_source",
            "compact_adp_photometry_reused",
            "compact_adp_flux_reused",
            "detrend_contract_version",
            "detrend_config_sha256",
            "detrend_quality_source",
        }
        if is_detector_consistent_v4:
            v3_required_attrs.update(
                {
                    "raw_detector_mapping_authority",
                    "raw_orbitid_authority",
                    "raw_internal_quality_authority",
                    "compact_cadence_time_inventory_match",
                    "compact_cadence_time_inventory_role",
                    "raw_compact_time_delta_max_s",
                }
            )
        else:
            v3_required_attrs.add(
                "raw_compact_internal_quality_agreement"
            )
        for aperture in ("small", "primary"):
            v3_required_attrs.update(
                {
                    f"effective_quality_adp_{aperture}_fit_count",
                    f"effective_quality_adp_{aperture}_n_segments",
                    f"effective_quality_adp_{aperture}_scale",
                    f"effective_quality_adp_{aperture}_scale_source",
                    f"effective_quality_adp_{aperture}_cotrend_status",
                    (
                        f"effective_quality_adp_{aperture}_"
                        "recenter_median_before"
                    ),
                    f"effective_quality_adp_{aperture}_n_effective_good",
                    f"effective_quality_adp_{aperture}_n_finite_output",
                }
            )
            if is_btjd:
                v3_required_attrs.update(
                    {
                        (
                            f"effective_quality_adp_{aperture}_"
                            "detrend_time_contract_version"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "detrend_time_system"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "detrend_time_min"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "detrend_time_max"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "warning_capture_policy"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "rank_warning_publication_policy"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "rank_warning_count"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "rank_warning_categories_json"
                        ),
                        (
                            f"effective_quality_adp_{aperture}_"
                            "rank_warning_messages_json"
                        ),
                    }
                )
        if is_btjd:
            v3_required_attrs.update(
                {
                    "detrend_time_contract_version",
                    "detrend_time_dataset",
                    "detrend_time_system",
                    "published_time_system",
                    "btjd_to_bjd_offset_d",
                    "detrend_time_min",
                    "detrend_time_max",
                }
            )
        missing_v3 = sorted(
            name for name in v3_required_attrs if name not in group.attrs
        )
        if missing_v3:
            return [
                f"{context}:missing v3 attrs {','.join(missing_v3)}"
            ], 0
        authority_mismatch = False
        if is_detector_consistent_v4:
            try:
                authority_mismatch = (
                    any(
                        int(group.attrs[name]) != 1
                        for name in (
                            "raw_detector_mapping_authority",
                            "raw_orbitid_authority",
                            "raw_internal_quality_authority",
                            "compact_cadence_time_inventory_match",
                        )
                    )
                    or str(
                        group.attrs["compact_cadence_time_inventory_role"]
                    )
                    != "reference_only"
                    or not np.isfinite(
                        float(group.attrs["raw_compact_time_delta_max_s"])
                    )
                    or float(group.attrs["raw_compact_time_delta_max_s"]) > 2.0
                )
            except (KeyError, TypeError, ValueError):
                authority_mismatch = True
        else:
            authority_mismatch = (
                int(group.attrs["raw_compact_internal_quality_agreement"])
                != 1
            )
        if (
            authority_mismatch
            or int(group.attrs["compact_adp_photometry_reused"]) != 0
            or int(group.attrs["compact_adp_flux_reused"]) != 0
            or str(group.attrs["photometry_source"])
            != "immutable_raw_v1_effective_quality_adp03q"
            or str(group.attrs["detrend_contract_version"])
            != (
                FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
                if is_btjd
                else FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION_V1
            )
            or str(group.attrs["detrend_config_sha256"])
            != FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
            or str(group.attrs["detrend_quality_source"])
            != "final_effective_quality"
        ):
            failures.append(f"{context}:v3 photometry provenance mismatch")
        for aperture in ("small", "primary"):
            prefix = f"effective_quality_adp_{aperture}"
            try:
                fit_count = int(group.attrs[f"{prefix}_fit_count"])
                n_segments = int(group.attrs[f"{prefix}_n_segments"])
                n_effective_good = int(
                    group.attrs[f"{prefix}_n_effective_good"]
                )
                n_finite_output = int(
                    group.attrs[f"{prefix}_n_finite_output"]
                )
                scale = float(group.attrs[f"{prefix}_scale"])
                center = float(
                    group.attrs[f"{prefix}_recenter_median_before"]
                )
            except (TypeError, ValueError):
                failures.append(
                    f"{context}:{aperture} ADP diagnostics are invalid"
                )
                continue
            if (
                fit_count < 1
                or n_segments < 1
                or n_effective_good < 1
                or n_finite_output < n_effective_good
                or not np.isfinite(scale)
                or scale <= 0
                or not np.isfinite(center)
                or not str(group.attrs[f"{prefix}_scale_source"])
                or not str(group.attrs[f"{prefix}_cotrend_status"])
            ):
                failures.append(
                    f"{context}:{aperture} ADP diagnostics are outside bounds"
                )
            if is_btjd:
                warning_attrs = {
                    "detrend_time_contract_version": (
                        FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
                    ),
                    "detrend_time_system": "BTJD",
                    "warning_capture_policy": (
                        FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY
                    ),
                    "rank_warning_publication_policy": (
                        FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
                    ),
                    "rank_warning_categories_json": "[]",
                    "rank_warning_messages_json": "[]",
                }
                for suffix, expected in warning_attrs.items():
                    if str(group.attrs[f"{prefix}_{suffix}"]) != expected:
                        failures.append(
                            f"{context}:{aperture} {suffix} mismatch"
                        )
                try:
                    if (
                        int(group.attrs[f"{prefix}_rank_warning_count"])
                        != 0
                        or float(group.attrs[f"{prefix}_detrend_time_min"])
                        != float(group.attrs["detrend_time_min"])
                        or float(group.attrs[f"{prefix}_detrend_time_max"])
                        != float(group.attrs["detrend_time_max"])
                    ):
                        failures.append(
                            f"{context}:{aperture} BTJD/warning audit mismatch"
                        )
                except (TypeError, ValueError):
                    failures.append(
                        f"{context}:{aperture} BTJD/warning audit is invalid"
                    )
        if is_btjd:
            if (
                str(group.attrs["detrend_time_contract_version"])
                != FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
                or str(group.attrs["detrend_time_dataset"])
                != FULL_POOL_NATIVE_DETREND_TIME_DATASET
                or str(group.attrs["detrend_time_system"]) != "BTJD"
                or str(group.attrs["published_time_system"]) != "BJD"
            ):
                failures.append(f"{context}:v3r1 time provenance mismatch")
            try:
                if float(group.attrs["btjd_to_bjd_offset_d"]) != (
                    FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
                ):
                    failures.append(f"{context}:v3r1 BJD offset mismatch")
            except (TypeError, ValueError):
                failures.append(f"{context}:v3r1 BJD offset is invalid")
            if (
                FULL_POOL_NATIVE_DETREND_TIME_DATASET not in group
                or "time" not in group
            ):
                failures.append(f"{context}:missing v3r1 BTJD/BJD time arrays")
            else:
                time_btjd = np.asarray(
                    group[FULL_POOL_NATIVE_DETREND_TIME_DATASET],
                    dtype=np.float64,
                )
                published_bjd = np.asarray(group["time"], dtype=np.float64)
                try:
                    time_btjd = _validated_btjd_time(
                        time_btjd,
                        context=context,
                    )
                except ValueError as exc:
                    failures.append(str(exc))
                else:
                    if (
                        len(published_bjd) != len(time_btjd)
                        or not np.isfinite(published_bjd).all()
                        or not np.array_equal(
                            published_bjd,
                            time_btjd
                            + FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
                        )
                    ):
                        failures.append(
                            f"{context}:published BJD differs from BTJD offset"
                        )
                    if (
                        float(group.attrs["detrend_time_min"])
                        != float(np.min(time_btjd))
                        or float(group.attrs["detrend_time_max"])
                        != float(np.max(time_btjd))
                    ):
                        failures.append(
                            f"{context}:BTJD extrema metadata mismatch"
                        )
        for name in PERIODOGRAM_DATASETS:
            if name not in group:
                failures.append(f"{context}:missing v3 periodogram {name}")
            elif len(group[name]) != FULL_POOL_NATIVE_PERIODOGRAM_N:
                failures.append(
                    f"{context}:{name} length differs from native-v3 contract"
                )
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
    input_pre = np.asarray(group["orbitid_pre_reference"], dtype=np.int64)
    raw_pre = np.asarray(group["raw_orbitid_pre_reference"], dtype=np.int64)
    raw_mask = np.asarray(group[ORBITID_RECONCILIATION_MASK_DATASET])
    if (
        len({len(cadences), len(stored), len(input_pre), len(raw_pre), len(raw_mask)})
        != 1
    ):
        return failures + [f"{context}:orbit-ID arrays differ in length"], corrected
    if not np.isin(raw_mask, (0, 1, False, True)).all():
        return failures + [f"{context}:correction mask is not binary"], corrected
    mask = raw_mask.astype(bool)
    if int(np.count_nonzero(mask)) != corrected:
        failures.append(f"{context}:correction count differs from mask")
    if not np.array_equal(raw_pre, input_pre):
        failures.append(f"{context}:raw pre-correction orbitids differ")
    if root_policy == ORBITID_POLICY_STRICT:
        if corrected or mask.any() or not np.array_equal(input_pre, stored):
            failures.append(f"{context}:strict policy changed orbit IDs")
    elif root_policy == ORBITID_POLICY_REFERENCE:
        if sector != 62:
            failures.append(f"{context}:reference correction is outside S62")
        try:
            _validate_s62_group_reconciliation(
                tic=tic,
                cadenceno=cadences,
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
    expected_git_sha: str | None = None,
    expected_sector: int | None = None,
    expected_shard_index: int | None = None,
    expected_n_shards: int | None = None,
    expected_observations: int | None = None,
    expected_source_observations: int | None = None,
    expected_excluded_observations: int | None = None,
    expected_source_identity_sha256: str | None = None,
    expected_output_identity_sha256: str | None = None,
    expected_excluded_identity_sha256: str | None = None,
    require_periodograms: bool = True,
) -> dict[str, Any]:
    """Read and validate every real group in one full-pool native shard."""

    import h5py

    path = Path(path)
    failures: list[str] = []
    counts = {"targets": 0, "injections": 0}
    group_rank_warning_count = 0
    if not path.is_file():
        return {
            "passed": False,
            "failures": [f"missing file: {path}"],
            "counts": counts,
        }
    metadata: dict[str, Any] = {}
    with h5py.File(path, "r") as h5:
        if str(h5.attrs.get("contract_version", "")) != (
            FULL_POOL_NATIVE_CONTRACT_VERSION
        ):
            failures.append("wrong full-pool native contract")
        failures.extend(full_pool_native_root_failures(h5.attrs))
        if str(h5.attrs.get("builder_code_sha256", "")) != file_sha256(
            Path(__file__)
        ):
            failures.append("builder_code_sha256 differs from verifier source")
        if expected_git_sha is not None and str(
            h5.attrs.get("expected_git_sha", "")
        ) != str(expected_git_sha).strip().lower():
            failures.append("expected_git_sha differs from expected value")
        for name, expected in (
            ("sector", expected_sector),
            ("shard_index", expected_shard_index),
            ("n_shards", expected_n_shards),
            ("n_shard_observations", expected_observations),
            ("n_source_shard_observations", expected_source_observations),
            (
                "n_shard_excluded_observations",
                expected_excluded_observations,
            ),
            (
                "n_shard_excluded_structurally_validated",
                expected_excluded_observations,
            ),
        ):
            if expected is not None and int(h5.attrs.get(name, -1)) != int(expected):
                failures.append(f"{name} differs from expected value")
        for name, expected in (
            (
                "source_shard_observation_identity_sha256",
                expected_source_identity_sha256,
            ),
            (
                "shard_observation_identity_sha256",
                expected_output_identity_sha256,
            ),
            (
                "shard_excluded_observation_identity_sha256",
                expected_excluded_identity_sha256,
            ),
        ):
            if expected is not None and str(h5.attrs.get(name, "")) != str(
                expected
            ):
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
                    for aperture in ("small", "primary"):
                        group_rank_warning_count += int(
                            group.attrs.get(
                                (
                                    "effective_quality_adp_"
                                    f"{aperture}_rank_warning_count"
                                ),
                                -1,
                            )
                        )
                except Exception as exc:
                    failures.append(f"/{group_path}:{exc}")
        if counts["targets"] != int(
            h5.attrs.get("n_shard_observations", -1)
        ):
            failures.append("target count differs from root declaration")
        if (
            group_rank_warning_count
            != int(h5.attrs.get("rank_warning_count", -1))
            or group_rank_warning_count != 0
        ):
            failures.append(
                "group/root RankWarning count differs or is nonzero"
            )
        metadata = {
            "contract_version": str(h5.attrs.get("contract_version", "")),
            "sector": int(h5.attrs.get("sector", -1)),
            "shard_index": int(h5.attrs.get("shard_index", -1)),
            "n_shards": int(h5.attrs.get("n_shards", -1)),
            "n_source_shard_observations": int(
                h5.attrs.get("n_source_shard_observations", -1)
            ),
            "n_shard_observations": int(
                h5.attrs.get("n_shard_observations", -1)
            ),
            "n_shard_excluded_observations": int(
                h5.attrs.get("n_shard_excluded_observations", -1)
            ),
            "source_shard_observation_identity_sha256": str(
                h5.attrs.get(
                    "source_shard_observation_identity_sha256", ""
                )
            ),
            "shard_observation_identity_sha256": str(
                h5.attrs.get("shard_observation_identity_sha256", "")
            ),
            "shard_excluded_observation_identity_sha256": str(
                h5.attrs.get(
                    "shard_excluded_observation_identity_sha256", ""
                )
            ),
            "eligibility_exclusions_sha256": str(
                h5.attrs.get("eligibility_exclusions_sha256", "")
            ),
            "eligibility_summary_sha256": str(
                h5.attrs.get("eligibility_summary_sha256", "")
            ),
            "raw_source_release_code_revision": str(
                h5.attrs.get("raw_source_release_code_revision", "")
            ),
            "raw_source_summary_sha256": str(
                h5.attrs.get("raw_source_summary_sha256", "")
            ),
            "raw_export_complete_sha256": str(
                h5.attrs.get("raw_export_complete_sha256", "")
            ),
            "raw_transfer_validation_sha256": str(
                h5.attrs.get("raw_transfer_validation_sha256", "")
            ),
            "expected_git_sha": str(
                h5.attrs.get("expected_git_sha", "")
            ),
            "builder_contract_version": str(
                h5.attrs.get("builder_contract_version", "")
            ),
            "builder_code_sha256": str(
                h5.attrs.get("builder_code_sha256", "")
            ),
            "detrend_contract_version": str(
                h5.attrs.get("detrend_contract_version", "")
            ),
            "detrend_config_sha256": str(
                h5.attrs.get("detrend_config_sha256", "")
            ),
            "periodogram_n": int(h5.attrs.get("periodogram_n", -1)),
            "detrend_time_contract_version": str(
                h5.attrs.get("detrend_time_contract_version", "")
            ),
            "detrend_time_dataset": str(
                h5.attrs.get("detrend_time_dataset", "")
            ),
            "detrend_time_system": str(
                h5.attrs.get("detrend_time_system", "")
            ),
            "published_time_system": str(
                h5.attrs.get("published_time_system", "")
            ),
            "btjd_to_bjd_offset_d": float(
                h5.attrs.get("btjd_to_bjd_offset_d", np.nan)
            ),
            "warning_capture_policy": str(
                h5.attrs.get("warning_capture_policy", "")
            ),
            "rank_warning_publication_policy": str(
                h5.attrs.get("rank_warning_publication_policy", "")
            ),
            "rank_warning_count": int(
                h5.attrs.get("rank_warning_count", -1)
            ),
            "rank_warning_ledger_json": str(
                h5.attrs.get("rank_warning_ledger_json", "")
            ),
            "rank_warning_ledger_sha256": str(
                h5.attrs.get("rank_warning_ledger_sha256", "")
            ),
        }
    return {
        "passed": not failures,
        "failures": failures,
        "counts": counts,
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "size_bytes": int(path.stat().st_size),
        **metadata,
    }


def write_full_pool_native_registry(
    *,
    pool_path: Path,
    pool_summary_path: Path,
    eligibility_exclusions_path: Path,
    eligibility_summary_path: Path,
    native_shard_paths: Sequence[Path],
    native_shard_summary_paths: Sequence[Path],
    source_path: Path,
    registry_path: Path,
    summary_path: Path,
    release_summary_path: Path,
    expected_shards_per_sector: int = FULL_POOL_NATIVE_SHARDS_PER_SECTOR,
    eligibility_production_lock: bool = True,
) -> dict[str, Any]:
    """Freeze the exact eligible native coverage and its full-pool audit."""

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
    full_keys = set(
        zip(pool["sector"].astype(int), pool["tic"].astype(int))
    )
    if not full_keys or len(full_keys) != len(pool):
        raise ValueError("frozen pool has empty or duplicate observation keys")
    eligibility = load_native_model_eligibility(
        eligibility_exclusions_path,
        eligibility_summary_path,
        pool_path=pool_binding.path,
        pool_summary_path=summary_binding.path,
        production_lock=bool(eligibility_production_lock),
        rederive_from_bls=False,
    )
    expected_keys = set(eligibility.eligible_keys)
    excluded_keys = set(eligibility.excluded_keys)
    if (
        set(eligibility.full_keys) != full_keys
        or expected_keys | excluded_keys != full_keys
        or expected_keys & excluded_keys
    ):
        raise ValueError("eligibility authority does not partition frozen pool")
    expected_shards_per_sector = int(expected_shards_per_sector)
    if expected_shards_per_sector < 1:
        raise ValueError("expected_shards_per_sector must be positive")

    records: list[dict[str, Any]] = []
    shard_metadata: list[dict[str, Any]] = []
    paths = [
        Path(value).expanduser().resolve(strict=True)
        for value in native_shard_paths
    ]
    shard_summaries = [
        Path(value).expanduser().resolve(strict=True)
        for value in native_shard_summary_paths
    ]
    if not paths or len(paths) != len(set(paths)):
        raise ValueError("native shard paths must be nonempty and unique")
    if (
        len(shard_summaries) != len(paths)
        or len(shard_summaries) != len(set(shard_summaries))
    ):
        raise ValueError(
            "native shard summary paths must be unique and pair every shard"
        )
    if len(paths) != (
        len(FULL_POOL_NATIVE_SECTORS) * expected_shards_per_sector
    ):
        raise ValueError("native shard inventory count is not exact")
    import h5py

    observed_shards: dict[int, set[int]] = {
        sector: set() for sector in FULL_POOL_NATIVE_SECTORS
    }
    expected_git_shas: set[str] = set()
    builder_code_shas: set[str] = set()
    for path, shard_summary_path in sorted(
        zip(paths, shard_summaries, strict=True),
        key=lambda pair: str(pair[0]),
    ):
        try:
            shard_summary = json.loads(
                shard_summary_path.read_text(encoding="utf-8")
            )
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ValueError(
                f"invalid native shard summary: {shard_summary_path}"
            ) from exc
        if not isinstance(shard_summary, dict):
            raise ValueError("native shard summary must be a JSON object")
        if (
            shard_summary.get("schema_version")
            != FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION
            or shard_summary.get("stage") != "orcd_native_shard_build"
            or shard_summary.get("native_contract_version")
            != FULL_POOL_NATIVE_CONTRACT_VERSION
            or shard_summary.get("real_only") is not True
            or shard_summary.get("verification", {}).get("passed") is not True
            or shard_summary.get("builder_contract_version")
            != FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
            or shard_summary.get("detrend_contract_version")
            != FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
            or shard_summary.get("detrend_config_sha256")
            != FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
            or shard_summary.get("detrend_config_json")
            != FULL_POOL_NATIVE_DETREND_CONFIG_JSON
            or shard_summary.get("detrend_quality_source")
            != "final_effective_quality"
            or shard_summary.get("detrend_time_contract_version")
            != FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
            or shard_summary.get("detrend_time_dataset")
            != FULL_POOL_NATIVE_DETREND_TIME_DATASET
            or shard_summary.get("detrend_time_system") != "BTJD"
            or shard_summary.get("published_time_system") != "BJD"
            or float(shard_summary.get("btjd_to_bjd_offset_d", np.nan))
            != FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
            or shard_summary.get("warning_capture_policy")
            != FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY
            or shard_summary.get("rank_warning_publication_policy")
            != FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
            or int(shard_summary.get("rank_warning_count", -1)) != 0
            or shard_summary.get("rank_warning_ledger_json")
            != FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON
            or shard_summary.get("rank_warning_ledger_sha256")
            != FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256
            or shard_summary.get("raw_photometry_only") is not True
            or shard_summary.get("compact_adp_photometry_reused") is not False
            or shard_summary.get("compact_adp_flux_reused") is not False
            or int(shard_summary.get("periodogram_n", -1))
            != FULL_POOL_NATIVE_PERIODOGRAM_N
        ):
            raise ValueError(
                f"native shard summary did not pass v3 contract: {path}"
            )
        expected_git_sha = str(shard_summary.get("expected_git_sha", ""))
        builder_code_sha = str(
            shard_summary.get("builder_code_sha256", "")
        )
        if (
            len(expected_git_sha) != 40
            or any(
                value not in "0123456789abcdef"
                for value in expected_git_sha
            )
            or len(builder_code_sha) != 64
            or any(
                value not in "0123456789abcdef"
                for value in builder_code_sha
            )
        ):
            raise ValueError("native shard summary code identity is invalid")
        expected_git_shas.add(expected_git_sha)
        builder_code_shas.add(builder_code_sha)
        if Path(str(shard_summary.get("out_h5", ""))).resolve() != path:
            raise ValueError(f"native shard summary binds a different HDF5: {path}")
        if (
            shard_summary.get("out_h5_sha256") != file_sha256(path)
            or int(shard_summary.get("out_h5_size_bytes", -1))
            != path.stat().st_size
        ):
            raise ValueError(f"native shard differs from its summary: {path}")
        sector = int(shard_summary.get("sector", -1))
        shard_index = int(shard_summary.get("shard_index", -1))
        n_shards = int(shard_summary.get("n_shards", -1))
        if (
            sector not in FULL_POOL_NATIVE_SECTORS
            or n_shards != expected_shards_per_sector
            or shard_index not in range(expected_shards_per_sector)
            or shard_index in observed_shards[sector]
        ):
            raise ValueError("native shard sector/index inventory is invalid")
        observed_shards[sector].add(shard_index)

        sector_rows = pool.loc[pool["sector"].eq(sector)].copy()
        assignments = sector_rows["tic"].map(
            lambda tic: shard_for_tic(
                sector=sector,
                tic=int(tic),
                n_shards=expected_shards_per_sector,
            )
        )
        source_rows = (
            sector_rows.loc[assignments.eq(shard_index)]
            .sort_values("tic", kind="stable")
            .reset_index(drop=True)
        )
        source_keys = set(
            zip(
                source_rows["sector"].astype(int),
                source_rows["tic"].astype(int),
            )
        )
        eligible_rows = source_rows.loc[
            [
                (int(row.sector), int(row.tic)) in expected_keys
                for row in source_rows.itertuples(index=False)
            ]
        ].reset_index(drop=True)
        excluded_rows = source_rows.loc[
            [
                (int(row.sector), int(row.tic)) in excluded_keys
                for row in source_rows.itertuples(index=False)
            ]
        ].reset_index(drop=True)
        if (
            source_keys
            != set(
                zip(
                    eligible_rows["sector"].astype(int),
                    eligible_rows["tic"].astype(int),
                )
            )
            | set(
                zip(
                    excluded_rows["sector"].astype(int),
                    excluded_rows["tic"].astype(int),
                )
            )
        ):
            raise ValueError("native shard source partition is incomplete")
        audit = verify_full_pool_native_shard(
            path,
            expected_sector=sector,
            expected_shard_index=shard_index,
            expected_n_shards=expected_shards_per_sector,
            expected_observations=len(eligible_rows),
            expected_source_observations=len(source_rows),
            expected_excluded_observations=len(excluded_rows),
            expected_source_identity_sha256=_observation_identity_sha256(
                source_rows
            ),
            expected_output_identity_sha256=_observation_identity_sha256(
                eligible_rows
            ),
            expected_excluded_identity_sha256=_observation_identity_sha256(
                excluded_rows
            ),
        )
        if not audit["passed"]:
            raise ValueError(
                f"native shard failed verification {path}: "
                + "; ".join(audit["failures"][:10])
            )
        for name, expected in (
            ("full_pool_sha256", pool_binding.sha256),
            ("full_pool_summary_sha256", summary_binding.sha256),
            (
                "eligibility_exclusions_sha256",
                eligibility.bindings["exclusions"].sha256,
            ),
            (
                "eligibility_summary_sha256",
                eligibility.bindings["eligibility_summary"].sha256,
            ),
            (
                "native_model_full_identity_sha256",
                eligibility.full_observation_identity_sha256,
            ),
            (
                "native_model_eligible_identity_sha256",
                eligibility.eligible_observation_identity_sha256,
            ),
            (
                "native_model_excluded_identity_sha256",
                eligibility.excluded_observation_identity_sha256,
            ),
            (
                "raw_source_release_code_revision",
                PRODUCTION_RAW_CODE_REVISION,
            ),
            (
                "raw_export_complete_sha256",
                PRODUCTION_RAW_EXPORT_COMPLETE_SHA256,
            ),
            (
                "raw_transfer_validation_sha256",
                PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR[sector],
            ),
        ):
            if shard_summary.get(name) != expected:
                raise ValueError(
                    f"native shard summary {name} binding differs: {path}"
                )
        for name, expected in (
            ("n_source_shard_observations", len(source_rows)),
            ("n_shard_observations", len(eligible_rows)),
            ("n_shard_excluded_observations", len(excluded_rows)),
            (
                "n_shard_excluded_structurally_validated",
                len(excluded_rows),
            ),
        ):
            if int(shard_summary.get(name, -1)) != expected:
                raise ValueError(
                    f"native shard summary {name} differs: {path}"
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
            if h5.attrs.get("eligibility_exclusions_sha256") != (
                eligibility.bindings["exclusions"].sha256
            ) or h5.attrs.get("eligibility_summary_sha256") != (
                eligibility.bindings["eligibility_summary"].sha256
            ):
                raise ValueError(
                    f"native shard binds a different eligibility release: {path}"
                )
            for name, expected in (
                (
                    "release_binding",
                    FULL_POOL_NATIVE_RELEASE_BINDING,
                ),
                (
                    "builder_contract_version",
                    FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION,
                ),
                ("builder_code_sha256", builder_code_sha),
                ("expected_git_sha", expected_git_sha),
                (
                    "detrend_contract_version",
                    FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION,
                ),
                (
                    "detrend_config_sha256",
                    FULL_POOL_NATIVE_DETREND_CONFIG_SHA256,
                ),
                (
                    "detrend_config_json",
                    FULL_POOL_NATIVE_DETREND_CONFIG_JSON,
                ),
                (
                    "detrend_quality_source",
                    "final_effective_quality",
                ),
                (
                    "detrend_time_contract_version",
                    FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION,
                ),
                (
                    "detrend_time_dataset",
                    FULL_POOL_NATIVE_DETREND_TIME_DATASET,
                ),
                ("detrend_time_system", "BTJD"),
                ("published_time_system", "BJD"),
                (
                    "btjd_to_bjd_offset_d",
                    FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
                ),
                (
                    "warning_capture_policy",
                    FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY,
                ),
                (
                    "rank_warning_publication_policy",
                    FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY,
                ),
                ("rank_warning_count", 0),
                (
                    "rank_warning_ledger_json",
                    FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON,
                ),
                (
                    "rank_warning_ledger_sha256",
                    FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256,
                ),
                ("raw_photometry_only", 1),
                ("compact_adp_photometry_reused", 0),
                ("compact_adp_flux_reused", 0),
                (
                    "raw_source_release_code_revision",
                    PRODUCTION_RAW_CODE_REVISION,
                ),
                (
                    "raw_export_complete_sha256",
                    PRODUCTION_RAW_EXPORT_COMPLETE_SHA256,
                ),
                (
                    "raw_transfer_validation_sha256",
                    PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR[sector],
                ),
                (
                    "raw_source_summary_sha256",
                    shard_summary.get("raw_source_summary_sha256"),
                ),
            ):
                if h5.attrs.get(name) != expected:
                    raise ValueError(
                        f"native shard raw-release binding differs: {path}"
                    )
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
                        "eligibility_exclusions_sha256": (
                            eligibility.bindings["exclusions"].sha256
                        ),
                        "eligibility_summary_sha256": (
                            eligibility.bindings["eligibility_summary"].sha256
                        ),
                    }
                )
        shard_metadata.append(
            {
                "path": str(path),
                "sha256": audit["sha256"],
                "size_bytes": audit["size_bytes"],
                "summary_path": str(shard_summary_path),
                "summary_sha256": file_sha256(shard_summary_path),
                "sector": sector,
                "shard_index": shard_index,
                "n_shards": n_shards,
                "n_source_observations": len(source_rows),
                "n_observations": len(eligible_rows),
                "n_excluded_observations": len(excluded_rows),
                "source_observation_identity_sha256": (
                    _observation_identity_sha256(source_rows)
                ),
                "observation_identity_sha256": (
                    _observation_identity_sha256(eligible_rows)
                ),
                "excluded_observation_identity_sha256": (
                    _observation_identity_sha256(excluded_rows)
                ),
                "raw_source_h5_sha256": shard_summary.get(
                    "raw_source_h5_sha256"
                ),
                "raw_source_summary_sha256": shard_summary.get(
                    "raw_source_summary_sha256"
                ),
                "raw_export_complete_sha256": shard_summary.get(
                    "raw_export_complete_sha256"
                ),
                "raw_transfer_validation_sha256": shard_summary.get(
                    "raw_transfer_validation_sha256"
                ),
                "raw_source_release_code_revision": shard_summary.get(
                    "raw_source_release_code_revision"
                ),
                "expected_git_sha": expected_git_sha,
                "builder_contract_version": (
                    FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
                ),
                "builder_code_sha256": builder_code_sha,
                "detrend_contract_version": (
                    FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
                ),
                "detrend_config_sha256": (
                    FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
                ),
                "detrend_quality_source": "final_effective_quality",
                "detrend_time_contract_version": (
                    FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
                ),
                "detrend_time_dataset": (
                    FULL_POOL_NATIVE_DETREND_TIME_DATASET
                ),
                "detrend_time_system": "BTJD",
                "published_time_system": "BJD",
                "btjd_to_bjd_offset_d": (
                    FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
                ),
                "warning_capture_policy": (
                    FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY
                ),
                "rank_warning_publication_policy": (
                    FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
                ),
                "rank_warning_count": 0,
                "rank_warning_ledger_json": (
                    FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON
                ),
                "rank_warning_ledger_sha256": (
                    FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256
                ),
                "raw_photometry_only": True,
                "compact_adp_photometry_reused": False,
                "compact_adp_flux_reused": False,
                "periodogram_n": FULL_POOL_NATIVE_PERIODOGRAM_N,
                "n_excluded_structurally_validated": len(excluded_rows),
                "verification_passed": True,
            }
        )
    for sector, indices in observed_shards.items():
        if indices != set(range(expected_shards_per_sector)):
            raise ValueError(f"S{sector} native shard index coverage is incomplete")
    if len(expected_git_shas) != 1 or len(builder_code_shas) != 1:
        raise ValueError("native shards do not share one exact code identity")
    if builder_code_shas != {file_sha256(Path(__file__))}:
        raise ValueError("native shard builder SHA differs from registry source")
    expected_git_sha = next(iter(expected_git_shas))
    builder_code_sha = next(iter(builder_code_shas))
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
    source_binding = _bind_file(source_path)
    registry_binding = _bind_file(registry_path)
    registry_summary_binding = _bind_file(summary_path)
    release_summary = {
        "passed": True,
        "schema_version": FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION,
        "release_binding": FULL_POOL_NATIVE_RELEASE_BINDING,
        "native_contract_version": FULL_POOL_NATIVE_CONTRACT_VERSION,
        "expected_git_sha": expected_git_sha,
        "builder_contract_version": (
            FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
        ),
        "builder_code_sha256": builder_code_sha,
        "detrend_contract_version": (
            FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
        ),
        "detrend_config_json": FULL_POOL_NATIVE_DETREND_CONFIG_JSON,
        "detrend_config_sha256": (
            FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
        ),
        "detrend_quality_source": "final_effective_quality",
        "detrend_time_contract_version": (
            FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
        ),
        "detrend_time_dataset": FULL_POOL_NATIVE_DETREND_TIME_DATASET,
        "detrend_time_system": "BTJD",
        "published_time_system": "BJD",
        "btjd_to_bjd_offset_d": FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
        "warning_capture_policy": FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY,
        "rank_warning_publication_policy": (
            FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
        ),
        "rank_warning_count": 0,
        "rank_warning_ledger_json": (
            FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON
        ),
        "rank_warning_ledger_sha256": (
            FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256
        ),
        "raw_photometry_only": True,
        "compact_adp_photometry_reused": False,
        "compact_adp_flux_reused": False,
        "periodogram_n": FULL_POOL_NATIVE_PERIODOGRAM_N,
        "eligibility_contract_version": eligibility.contract_version,
        "sectors": list(FULL_POOL_NATIVE_SECTORS),
        "shards_per_sector": expected_shards_per_sector,
        "counts": {
            "full_observations": len(full_keys),
            "eligible_observations": len(expected_keys),
            "excluded_observations": len(excluded_keys),
            "native_registry_observations": len(source),
            "native_shards": len(shard_metadata),
        },
        "identity_hashes": {
            "full_observations_sha256": (
                eligibility.full_observation_identity_sha256
            ),
            "eligible_observations_sha256": (
                eligibility.eligible_observation_identity_sha256
            ),
            "excluded_observations_sha256": (
                eligibility.excluded_observation_identity_sha256
            ),
            "native_registry_observations_sha256": (
                _observation_identity_sha256(source)
            ),
        },
        "partition_audit": {
            "eligible_excluded_disjoint": not bool(
                expected_keys & excluded_keys
            ),
            "eligible_excluded_union_equals_full": (
                expected_keys | excluded_keys == full_keys
            ),
            "native_registry_equals_eligible": observed_keys == expected_keys,
            "excluded_present_in_native_registry": bool(
                excluded_keys & observed_keys
            ),
        },
        "source_authorities": {
            "frozen_pool": {
                "path": str(pool_binding.path),
                "sha256": pool_binding.sha256,
                "size_bytes": pool_binding.size_bytes,
            },
            "frozen_pool_summary": {
                "path": str(summary_binding.path),
                "sha256": summary_binding.sha256,
                "size_bytes": summary_binding.size_bytes,
            },
            "eligibility_exclusions": (
                eligibility.bindings["exclusions"].metadata()
            ),
            "eligibility_summary": (
                eligibility.bindings["eligibility_summary"].metadata()
            ),
            "raw_v1_release": {
                "code_revision": PRODUCTION_RAW_CODE_REVISION,
                "raw_export_complete_sha256": (
                    PRODUCTION_RAW_EXPORT_COMPLETE_SHA256
                ),
                "raw_transfer_validation_sha256_by_sector": {
                    str(key): value
                    for key, value in (
                        PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR.items()
                    )
                },
            },
        },
        "outputs": {
            "source_table": {
                "path": str(source_binding.path),
                "sha256": source_binding.sha256,
                "size_bytes": source_binding.size_bytes,
            },
            "native_registry": {
                "path": str(registry_binding.path),
                "sha256": registry_binding.sha256,
                "size_bytes": registry_binding.size_bytes,
            },
            "native_registry_summary": {
                "path": str(registry_summary_binding.path),
                "sha256": registry_summary_binding.sha256,
                "size_bytes": registry_summary_binding.size_bytes,
            },
        },
        "native_shards": shard_metadata,
    }
    release_summary_path = Path(release_summary_path).expanduser().resolve()
    _publish_immutable_bytes(
        release_summary_path,
        (
            json.dumps(
                release_summary,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8"),
    )
    pool_binding.assert_unchanged()
    summary_binding.assert_unchanged()
    return {
        **registry_summary,
        "passed": True,
        "full_pool_sha256": pool_binding.sha256,
        "full_pool_summary_sha256": summary_binding.sha256,
        "eligibility_exclusions_sha256": (
            eligibility.bindings["exclusions"].sha256
        ),
        "eligibility_summary_sha256": (
            eligibility.bindings["eligibility_summary"].sha256
        ),
        "native_shards": shard_metadata,
        "coverage_identity_sha256": _observation_identity_sha256(source),
        "release_summary": str(release_summary_path),
        "release_summary_sha256": file_sha256(release_summary_path),
    }


def load_full_pool_native_registry_release(
    *,
    registry_path: Path,
    registry_summary_path: Path,
    release_summary_path: Path,
    eligibility: EligibilityAuthority | None = None,
    verify_shard_files: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate the generic registry and the stronger full-pool v2 release.

    ``verify_shard_files=False`` retains all release, inventory, and declared
    binding checks but defers the expensive HDF5 hashes.  The training gate
    uses that mode immediately before hashing every selected native file
    itself, avoiding a duplicate full read of the same 112 shards.
    """

    registry_path = Path(registry_path).expanduser().resolve(strict=True)
    registry_summary_path = (
        Path(registry_summary_path).expanduser().resolve(strict=True)
    )
    release_summary_path = (
        Path(release_summary_path).expanduser().resolve(strict=True)
    )
    validate_native_input_registry_path(
        registry_path=registry_path,
        summary_path=registry_summary_path,
        expected_contract_version=FULL_POOL_NATIVE_CONTRACT_VERSION,
        verify_files=verify_shard_files,
    )
    registry = read_table(registry_path)
    try:
        release = json.loads(
            release_summary_path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(
            f"invalid full-pool native release summary: {release_summary_path}"
        ) from exc
    if not isinstance(release, dict):
        raise ValueError("full-pool native release summary must be a JSON object")
    if (
        release.get("passed") is not True
        or release.get("schema_version")
        != FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION
        or release.get("release_binding") != FULL_POOL_NATIVE_RELEASE_BINDING
        or release.get("native_contract_version")
        != FULL_POOL_NATIVE_CONTRACT_VERSION
        or release.get("builder_contract_version")
        != FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
        or release.get("builder_code_sha256")
        != file_sha256(Path(__file__))
        or release.get("detrend_contract_version")
        != FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
        or release.get("detrend_config_json")
        != FULL_POOL_NATIVE_DETREND_CONFIG_JSON
        or release.get("detrend_config_sha256")
        != FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
        or release.get("detrend_quality_source")
        != "final_effective_quality"
        or release.get("detrend_time_contract_version")
        != FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION
        or release.get("detrend_time_dataset")
        != FULL_POOL_NATIVE_DETREND_TIME_DATASET
        or release.get("detrend_time_system") != "BTJD"
        or release.get("published_time_system") != "BJD"
        or float(release.get("btjd_to_bjd_offset_d", np.nan))
        != FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D
        or release.get("warning_capture_policy")
        != FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY
        or release.get("rank_warning_publication_policy")
        != FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY
        or int(release.get("rank_warning_count", -1)) != 0
        or release.get("rank_warning_ledger_json")
        != FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON
        or release.get("rank_warning_ledger_sha256")
        != FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256
        or release.get("raw_photometry_only") is not True
        or release.get("compact_adp_photometry_reused") is not False
        or release.get("compact_adp_flux_reused") is not False
        or int(release.get("periodogram_n", -1))
        != FULL_POOL_NATIVE_PERIODOGRAM_N
        or release.get("eligibility_contract_version")
        != NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION
        or release.get("sectors") != list(FULL_POOL_NATIVE_SECTORS)
    ):
        raise ValueError("full-pool native release has the wrong v3 contract")
    expected_git_sha = str(release.get("expected_git_sha", ""))
    if len(expected_git_sha) != 40 or any(
        value not in "0123456789abcdef" for value in expected_git_sha
    ):
        raise ValueError("full-pool native release has invalid Git identity")
    outputs = release.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError("full-pool native release lacks outputs")
    for name, path in (
        ("native_registry", registry_path),
        ("native_registry_summary", registry_summary_path),
    ):
        metadata = outputs.get(name)
        if not isinstance(metadata, dict) or (
            metadata.get("sha256") != file_sha256(path)
            or int(metadata.get("size_bytes", -1)) != path.stat().st_size
        ):
            raise ValueError(f"full-pool native release {name} binding differs")
    keys = set(
        zip(
            pd.to_numeric(registry["sector"], errors="raise").astype(int),
            pd.to_numeric(registry["tic"], errors="raise").astype(int),
        )
    )
    if len(keys) != len(registry):
        raise ValueError("full-pool native registry has duplicate keys")
    counts = release.get("counts")
    identities = release.get("identity_hashes")
    partition = release.get("partition_audit")
    if (
        not isinstance(counts, dict)
        or not isinstance(identities, dict)
        or partition
        != {
            "eligible_excluded_disjoint": True,
            "eligible_excluded_union_equals_full": True,
            "native_registry_equals_eligible": True,
            "excluded_present_in_native_registry": False,
        }
        or int(counts.get("native_registry_observations", -1))
        != len(registry)
        or identities.get("native_registry_observations_sha256")
        != _observation_identity_sha256(registry)
    ):
        raise ValueError("full-pool native release coverage audit failed")
    if eligibility is not None:
        if (
            keys != set(eligibility.eligible_keys)
            or int(counts.get("full_observations", -1))
            != eligibility.n_full
            or int(counts.get("eligible_observations", -1))
            != eligibility.n_eligible
            or int(counts.get("excluded_observations", -1))
            != eligibility.n_excluded
            or identities.get("full_observations_sha256")
            != eligibility.full_observation_identity_sha256
            or identities.get("eligible_observations_sha256")
            != eligibility.eligible_observation_identity_sha256
            or identities.get("excluded_observations_sha256")
            != eligibility.excluded_observation_identity_sha256
        ):
            raise ValueError(
                "full-pool native release differs from eligibility authority"
            )
    source_authorities = release.get("source_authorities")
    raw_release = (
        source_authorities.get("raw_v1_release")
        if isinstance(source_authorities, dict)
        else None
    )
    if raw_release != {
        "code_revision": PRODUCTION_RAW_CODE_REVISION,
        "raw_export_complete_sha256": PRODUCTION_RAW_EXPORT_COMPLETE_SHA256,
        "raw_transfer_validation_sha256_by_sector": {
            str(key): value
            for key, value in PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR.items()
        },
    }:
        raise ValueError("full-pool native release raw-v1 authority differs")
    shards = release.get("native_shards")
    if not isinstance(shards, list) or int(
        counts.get("native_shards", -1)
    ) != len(shards):
        raise ValueError("full-pool native release shard inventory is invalid")
    expected_pairs = {
        (sector, shard_index)
        for sector in FULL_POOL_NATIVE_SECTORS
        for shard_index in range(FULL_POOL_NATIVE_SHARDS_PER_SECTOR)
    }
    observed_pairs: set[tuple[int, int]] = set()
    observed_paths: set[str] = set()
    observed_summary_paths: set[str] = set()
    for item in shards:
        if not isinstance(item, dict) or item.get("verification_passed") is not True:
            raise ValueError("full-pool native release contains an unverified shard")
        sector = int(item.get("sector", -1))
        shard_index = int(item.get("shard_index", -1))
        pair = (sector, shard_index)
        if (
            pair in observed_pairs
            or pair not in expected_pairs
            or int(item.get("n_shards", -1))
            != FULL_POOL_NATIVE_SHARDS_PER_SECTOR
        ):
            raise ValueError("full-pool native release shard layout is invalid")
        observed_pairs.add(pair)
        shard_path = Path(str(item.get("path", ""))).resolve()
        shard_summary_path = Path(
            str(item.get("summary_path", ""))
        ).resolve()
        if (
            str(shard_path) in observed_paths
            or str(shard_summary_path) in observed_summary_paths
        ):
            raise ValueError("full-pool native release repeats a shard artifact")
        observed_paths.add(str(shard_path))
        observed_summary_paths.add(str(shard_summary_path))
        if (
            not shard_path.is_file()
            or not shard_summary_path.is_file()
            or int(item.get("size_bytes", -1)) != shard_path.stat().st_size
            or item.get("summary_sha256")
            != file_sha256(shard_summary_path)
        ):
            raise ValueError("full-pool native release shard binding changed")
        if verify_shard_files and item.get("sha256") != file_sha256(
            shard_path
        ):
            raise ValueError("full-pool native release shard hash changed")
        shard_rows = registry.loc[
            registry["native_h5_path"].astype(str).eq(str(shard_path))
        ]
        if (
            len(shard_rows) != int(item.get("n_observations", -1))
            or shard_rows.empty
            or not shard_rows["sector"].astype(int).eq(sector).all()
            or _observation_identity_sha256(shard_rows)
            != item.get("observation_identity_sha256")
        ):
            raise ValueError(
                "full-pool native release shard rows differ from registry"
            )
        expected_groups = shard_rows["tic"].astype(int).map(
            lambda tic: f"targets/{tic:016d}"
        )
        if not shard_rows["native_group_path"].astype(str).eq(
            expected_groups
        ).all():
            raise ValueError(
                "full-pool native release shard group paths are invalid"
            )
        for name, expected in (
            ("expected_git_sha", expected_git_sha),
            (
                "builder_contract_version",
                FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION,
            ),
            ("builder_code_sha256", file_sha256(Path(__file__))),
            (
                "detrend_contract_version",
                FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION,
            ),
            (
                "detrend_config_sha256",
                FULL_POOL_NATIVE_DETREND_CONFIG_SHA256,
            ),
            ("detrend_quality_source", "final_effective_quality"),
            (
                "detrend_time_contract_version",
                FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION,
            ),
            (
                "detrend_time_dataset",
                FULL_POOL_NATIVE_DETREND_TIME_DATASET,
            ),
            ("detrend_time_system", "BTJD"),
            ("published_time_system", "BJD"),
            (
                "btjd_to_bjd_offset_d",
                FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
            ),
            (
                "warning_capture_policy",
                FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY,
            ),
            (
                "rank_warning_publication_policy",
                FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY,
            ),
            ("rank_warning_count", 0),
            (
                "rank_warning_ledger_json",
                FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON,
            ),
            (
                "rank_warning_ledger_sha256",
                FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256,
            ),
            ("raw_photometry_only", True),
            ("compact_adp_photometry_reused", False),
            ("compact_adp_flux_reused", False),
            ("periodogram_n", FULL_POOL_NATIVE_PERIODOGRAM_N),
            ("raw_source_release_code_revision", PRODUCTION_RAW_CODE_REVISION),
            (
                "raw_export_complete_sha256",
                PRODUCTION_RAW_EXPORT_COMPLETE_SHA256,
            ),
            (
                "raw_transfer_validation_sha256",
                PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR[sector],
            ),
        ):
            if item.get(name) != expected:
                raise ValueError(
                    "full-pool native release raw-source binding differs"
                )
    if observed_pairs != expected_pairs:
        raise ValueError("full-pool native release is not the exact 7x16 inventory")
    registry_paths = set(registry["native_h5_path"].astype(str))
    if observed_paths != registry_paths:
        raise ValueError(
            "full-pool native release shard paths differ from registry"
        )
    return registry, release


__all__ = [
    "FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION",
    "FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION_V1",
    "FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D",
    "FULL_POOL_NATIVE_CONTRACT_VERSION",
    "FULL_POOL_NATIVE_CONTRACT_VERSION_V1",
    "FULL_POOL_NATIVE_CONTRACT_VERSION_V2",
    "FULL_POOL_NATIVE_CONTRACT_VERSION_V3",
    "FULL_POOL_NATIVE_DETREND_CONFIG_JSON",
    "FULL_POOL_NATIVE_DETREND_CONFIG_SHA256",
    "FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION_V1",
    "FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION",
    "FULL_POOL_NATIVE_DETREND_TIME_DATASET",
    "FULL_POOL_NATIVE_PERIODOGRAM_N",
    "FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION",
    "FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON",
    "FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256",
    "FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY",
    "FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION",
    "FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION_V2",
    "FULL_POOL_NATIVE_REGISTRY_RELEASE_SCHEMA_VERSION_V3",
    "FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION",
    "FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION_V2",
    "FULL_POOL_NATIVE_REGISTRY_SOURCE_SCHEMA_VERSION_V3",
    "FULL_POOL_NATIVE_RELEASE_BINDING",
    "FULL_POOL_NATIVE_RELEASE_BINDING_V1",
    "FULL_POOL_NATIVE_RELEASE_BINDING_V2",
    "FULL_POOL_NATIVE_RELEASE_BINDING_V3",
    "FULL_POOL_NATIVE_SECTORS",
    "FULL_POOL_NATIVE_SHARDS_PER_SECTOR",
    "FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION",
    "FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION_V2",
    "FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION_V3",
    "FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY",
    "FULL_POOL_RAW_EXPORT_CONTROLLER_SCHEMA_VERSION",
    "FULL_POOL_RAW_SHARD_SUMMARY_SCHEMA_VERSION",
    "FULL_POOL_RAW_SOURCE_CONTRACT_VERSION",
    "FULL_POOL_RAW_TRANSFER_SCHEMA_VERSION",
    "PRODUCTION_RAW_CODE_REVISION",
    "PRODUCTION_RAW_EXPORT_COMPLETE_SHA256",
    "PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR",
    "RawSourceReleaseAuthority",
    "SectorPoolAuthority",
    "build_full_pool_native_shard",
    "export_full_pool_raw_source_shard",
    "full_pool_native_group_failures",
    "full_pool_native_root_failures",
    "load_production_raw_source_release",
    "load_sector_pool_authority",
    "load_full_pool_native_registry_release",
    "select_sector_shard",
    "shard_for_tic",
    "verify_full_pool_native_shard",
    "write_full_pool_native_registry",
]
