"""Checksum-bound native/model eligibility for the broad SSL light-curve pool.

The frozen full-pool table remains the survey/search and final-audit authority.
This module derives a narrower native/model partition only from the global BLS
ADP-small anchor row for each observation.  It intentionally does not inspect
raw flux uncertainties, human labels, injections, or per-aperture status-row
counts.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from twirl.search.a2v1_bls_contract import (
    approved_a2v1_teacher_bls_config,
    bls_config_sha256,
)
from twirl.vetting.adp_only import ADP_ONLY_APERTURES
from twirl.vetting.ssl_full_pool import (
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
)
from twirl.vetting.ssl_full_pool_bls import (
    FULL_POOL_BLS_SOURCE_PRODUCT_TAG,
    GLOBAL_BLS_CONTRACT_VERSION,
    GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
    ORBITID_POLICY_BY_SECTOR,
    ORBITID_RECONCILIATION_CONTRACT_VERSION,
)
from twirl.vetting.teacher_native_registry import file_sha256, read_table


NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_native_model_eligibility_v2"
)
NATIVE_MODEL_ELIGIBILITY_SUMMARY_SCHEMA_VERSION = (
    "twirl_teacher_ssl_fullpool_native_model_eligibility_summary_v2"
)
NATIVE_MODEL_ELIGIBILITY_RELEASE_BINDING = (
    "teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible"
)
NATIVE_MODEL_ELIGIBILITY_ANCHOR_APERTURE = "DET_FLUX_ADP_SML"
NATIVE_MODEL_ELIGIBILITY_EXCLUSION_REASON = "bls_unsearchable"

PRODUCTION_FULL_OBSERVATIONS = 175_366
PRODUCTION_ELIGIBLE_OBSERVATIONS = 175_347
PRODUCTION_EXCLUDED_OBSERVATIONS = 19
PRODUCTION_FULL_IDENTITY_SHA256 = (
    "8e9e9c12a24d5ebc7be94b03a4e35411cd10066d62a87d921a8443b06cc188d1"
)
PRODUCTION_ELIGIBLE_IDENTITY_SHA256 = (
    "6ddc8e57bb5fb938ce05389c1629221c73e0e73ac3bf40da47a2019e1a5660e6"
)
PRODUCTION_EXCLUDED_IDENTITY_SHA256 = (
    "b9f536144265e54a70bff17c782e0668ddbd96efcdf6c223ebc58f46edb7d976"
)
PRODUCTION_FROZEN_POOL_CSV_SHA256 = (
    "64bf3dc6dc5a34a702c410838b8d792ee89a59a7082719e9c704d11dcae03745"
)
PRODUCTION_FROZEN_POOL_SUMMARY_SHA256 = (
    "f9bb5afa6672db0d74f4a7de6c0c6064295e4e0bfae420bf78777ec101b303c4"
)
PRODUCTION_GLOBAL_BLS_SHA256 = (
    "160bc2b9c9a17b05b3607ef26d77981ec67163565affeb4b5f9cf1f10cd322e1"
)
PRODUCTION_GLOBAL_BLS_SUMMARY_SHA256 = (
    "e31ddd76eaa4946a48901fe76b1e4f28ab79d417e3c684e996d1c9bcb355fe38"
)
# The frozen full-pool BLS job supplied the approved parameters explicitly,
# so its immutable summary records the search-contract label as ``custom``.
# The exact config and config SHA below remain identical to the named preset.
PRODUCTION_GLOBAL_BLS_SEARCH_CONTRACT_VERSION = "custom"
PRODUCTION_EXCLUDED_BY_SECTOR: Mapping[int, int] = {
    56: 4,
    57: 1,
    58: 1,
    59: 4,
    60: 2,
    61: 3,
    62: 4,
}

ELIGIBILITY_EXCLUSION_COLUMNS: tuple[str, ...] = (
    "contract_version",
    "observation_id",
    "sector",
    "tic",
    "anchor_aperture",
    "peak_rank",
    "bls_status",
    "exclusion_reason",
)

_SUMMARY_TOP_LEVEL_FIELDS = {
    "passed",
    "schema_version",
    "contract_version",
    "release_binding",
    "sectors",
    "observation_identity_columns",
    "anchor_policy",
    "source_authorities",
    "counts",
    "by_sector",
    "by_reason",
    "identity_hashes",
    "partition_audit",
    "data_usage_audit",
    "outputs",
}
_GLOBAL_BLS_SUMMARY_FIELDS = {
    "passed",
    "schema_version",
    "contract_version",
    "sectors",
    "observation_identity_columns",
    "frozen_pool",
    "bls_contract",
    "counts",
    "coverage_audit",
    "sector_products",
    "output",
}
_FORBIDDEN_LABEL_COLUMNS = {
    "human_label",
    "label",
    "morphology_label",
    "morphology_target",
    "morphology_target_index",
    "preserve_target",
    "preserve_target_index",
    "harmonic_target",
    "harmonic_target_index",
    "teacher_v3_training_include",
}
_SHA256_PATTERN = r"^[0-9a-f]{64}$"
_MAX_INT64 = np.iinfo(np.int64).max

ObservationKey = tuple[int, int]


@dataclass(frozen=True)
class ArtifactBinding:
    """One immutable artifact identity exposed by an eligibility authority."""

    path: Path
    size_bytes: int
    sha256: str

    def metadata(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "size_bytes": int(self.size_bytes),
            "sha256": str(self.sha256),
        }


@dataclass(frozen=True)
class EligibilityAuthority:
    """Validated full/eligible/excluded partition and transitive bindings."""

    contract_version: str
    release_binding: str
    anchor_aperture: str
    full_keys: frozenset[ObservationKey]
    eligible_keys: frozenset[ObservationKey]
    excluded_keys: frozenset[ObservationKey]
    n_full: int
    n_eligible: int
    n_excluded: int
    full_observation_identity_sha256: str
    eligible_observation_identity_sha256: str
    excluded_observation_identity_sha256: str
    bindings: Mapping[str, ArtifactBinding]
    exclusions: pd.DataFrame
    summary: Mapping[str, Any]


def _no_duplicate_json_keys(
    pairs: list[tuple[str, Any]],
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"JSON object contains duplicate key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON contains non-finite constant {value!r}")


def _load_json(path: Path, *, context: str) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_no_duplicate_json_keys,
            parse_constant=_reject_json_constant,
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{context} is not strict UTF-8 JSON: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a JSON object")
    return value


def _binding(path: Path) -> ArtifactBinding:
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
    return ArtifactBinding(
        path=resolved,
        size_bytes=int(after.st_size),
        sha256=digest,
    )


def _positive_integer_series(values: pd.Series, *, name: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    invalid = (
        ~np.isfinite(numeric)
        | (numeric <= 0)
        | (numeric > _MAX_INT64)
        | ~np.equal(numeric, np.rint(numeric))
    )
    if invalid.any():
        examples = values.loc[invalid].head(5).tolist()
        raise ValueError(f"{name} must contain positive integers; first={examples}")
    return pd.Series(
        np.rint(numeric).astype(np.int64),
        index=values.index,
        name=name,
    )


def _integer_series(values: pd.Series, *, name: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    invalid = (
        ~np.isfinite(numeric)
        | (numeric < np.iinfo(np.int64).min)
        | (numeric > _MAX_INT64)
        | ~np.equal(numeric, np.rint(numeric))
    )
    if invalid.any():
        examples = values.loc[invalid].head(5).tolist()
        raise ValueError(f"{name} must contain finite integers; first={examples}")
    return pd.Series(
        np.rint(numeric).astype(np.int64),
        index=values.index,
        name=name,
    )


def _normalize_keys(
    rows_or_keys: pd.DataFrame | Iterable[ObservationKey],
    *,
    context: str,
) -> list[ObservationKey]:
    if isinstance(rows_or_keys, pd.DataFrame):
        missing = sorted({"sector", "tic"} - set(rows_or_keys.columns))
        if missing:
            raise KeyError(f"{context} lacks columns: {missing}")
        sectors = _positive_integer_series(
            rows_or_keys["sector"],
            name=f"{context} sector",
        )
        tics = _positive_integer_series(
            rows_or_keys["tic"],
            name=f"{context} tic",
        )
        keys = list(zip(sectors.astype(int), tics.astype(int), strict=True))
    else:
        keys = []
        for value in rows_or_keys:
            if (
                not isinstance(value, tuple)
                or len(value) != 2
                or isinstance(value[0], (bool, np.bool_))
                or isinstance(value[1], (bool, np.bool_))
            ):
                raise ValueError(f"{context} keys must be (sector, tic) integer tuples")
            sector, tic = int(value[0]), int(value[1])
            if sector <= 0 or tic <= 0:
                raise ValueError(f"{context} keys must be positive")
            keys.append((sector, tic))
    if len(keys) != len(set(keys)):
        raise ValueError(f"{context} contains duplicate observation keys")
    return sorted(keys)


def observation_identity_sha256(
    rows_or_keys: pd.DataFrame | Iterable[ObservationKey],
) -> str:
    """Return the canonical SHA-256 for a unique ``(sector, tic)`` inventory."""

    keys = _normalize_keys(rows_or_keys, context="observation identity")
    digest = hashlib.sha256()
    for sector, tic in keys:
        digest.update(
            json.dumps(
                {"sector": int(sector), "tic": int(tic)},
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _normalize_frozen_pool(frozen_pool: pd.DataFrame) -> pd.DataFrame:
    required = {"pool_contract_version", "observation_id", "sector", "tic"}
    missing = sorted(required - set(frozen_pool.columns))
    if missing:
        raise KeyError(f"frozen SSL pool lacks columns: {missing}")
    if frozen_pool.empty:
        raise ValueError("frozen SSL pool is empty")
    pool = frozen_pool.copy()
    pool["sector"] = _positive_integer_series(pool["sector"], name="sector")
    pool["tic"] = _positive_integer_series(pool["tic"], name="tic")
    if pool.duplicated(["sector", "tic"]).any():
        raise ValueError("frozen SSL pool has duplicate observation keys")
    if (
        not pool["pool_contract_version"]
        .astype(str)
        .eq(FULL_POOL_CONTRACT_VERSION)
        .all()
    ):
        raise ValueError("frozen SSL pool has the wrong contract version")
    expected_ids = [
        f"s{int(sector):04d}-tic{int(tic):016d}"
        for sector, tic in zip(pool["sector"], pool["tic"], strict=True)
    ]
    if pool["observation_id"].fillna("").astype(str).tolist() != expected_ids:
        raise ValueError("frozen SSL pool observation IDs do not match (sector, tic)")
    return pool.sort_values(["sector", "tic"], kind="stable").reset_index(drop=True)


def derive_anchor_eligibility(
    bls_rows: pd.DataFrame,
    *,
    frozen_pool: pd.DataFrame,
    anchor_aperture: str = NATIVE_MODEL_ELIGIBILITY_ANCHOR_APERTURE,
) -> pd.DataFrame:
    """Derive one native/model eligibility decision per frozen observation.

    The only eligible row is an ADP-small rank-one anchor with ``status=ok``,
    finite positive period, finite epoch, and finite positive duration.  A
    rank-zero status row retains an otherwise unsearchable observation in the
    audit partition.
    """

    forbidden = sorted(_FORBIDDEN_LABEL_COLUMNS & set(bls_rows.columns))
    if forbidden:
        raise ValueError(
            f"global BLS observations must be label-free; found={forbidden}"
        )
    required = {
        "sector",
        "tic",
        "aperture",
        "peak_rank",
        "status",
        "period_d",
        "t0_bjd",
        "duration_min",
    }
    missing = sorted(required - set(bls_rows.columns))
    if missing:
        raise KeyError(f"global BLS table lacks columns: {missing}")
    if bls_rows.empty:
        raise ValueError("global BLS table is empty")
    pool = _normalize_frozen_pool(frozen_pool)

    ranks = _integer_series(
        bls_rows["peak_rank"],
        name="peak_rank",
    )
    anchor_mask = bls_rows["aperture"].fillna("").astype(str).eq(
        str(anchor_aperture)
    ) & ranks.isin((0, 1))
    anchor = bls_rows.loc[anchor_mask, list(required)].copy()
    if anchor.empty:
        raise ValueError(
            f"BLS table contains no rank-one/status rows for {anchor_aperture}"
        )
    anchor["sector"] = _positive_integer_series(
        anchor["sector"],
        name="sector",
    )
    anchor["tic"] = _positive_integer_series(anchor["tic"], name="tic")
    anchor["peak_rank"] = ranks.loc[anchor_mask].astype(np.int64)
    duplicate = anchor.duplicated(["sector", "tic"], keep=False)
    if duplicate.any():
        examples = (
            anchor.loc[duplicate, ["sector", "tic", "peak_rank", "status"]]
            .head(10)
            .to_dict(orient="records")
        )
        raise ValueError(
            "ADP-small rank-one/status anchor is not unique per observation; "
            f"first={examples}"
        )
    pool_keys = set(_normalize_keys(pool, context="frozen SSL pool"))
    anchor_keys = set(_normalize_keys(anchor, context="BLS anchor"))
    if anchor_keys != pool_keys:
        missing_keys = sorted(pool_keys - anchor_keys)[:10]
        extra_keys = sorted(anchor_keys - pool_keys)[:10]
        raise ValueError(
            "BLS anchor coverage differs from the frozen SSL pool; "
            f"missing={missing_keys}, extra={extra_keys}"
        )

    anchor["status"] = anchor["status"].fillna("").astype(str).str.strip().str.lower()
    period = pd.to_numeric(anchor["period_d"], errors="coerce")
    epoch = pd.to_numeric(anchor["t0_bjd"], errors="coerce")
    duration = pd.to_numeric(anchor["duration_min"], errors="coerce")
    eligible = (
        anchor["status"].eq("ok")
        & anchor["peak_rank"].eq(1)
        & np.isfinite(period)
        & period.gt(0)
        & np.isfinite(epoch)
        & np.isfinite(duration)
        & duration.gt(0)
    )
    anchor["native_model_eligible"] = eligible.astype(bool)
    anchor = (
        anchor.loc[
            :,
            [
                "sector",
                "tic",
                "aperture",
                "peak_rank",
                "status",
                "period_d",
                "t0_bjd",
                "duration_min",
                "native_model_eligible",
            ],
        ]
        .merge(
            pool.loc[:, ["sector", "tic", "observation_id"]],
            on=["sector", "tic"],
            how="inner",
            validate="one_to_one",
        )
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    anchor["exclusion_reason"] = np.where(
        anchor["native_model_eligible"],
        "",
        NATIVE_MODEL_ELIGIBILITY_EXCLUSION_REASON,
    )
    return anchor


def _validate_frozen_pool_authority(
    *,
    pool_binding: ArtifactBinding,
    summary_binding: ArtifactBinding,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    pool = _normalize_frozen_pool(read_table(pool_binding.path))
    summary = _load_json(
        summary_binding.path,
        context="frozen full-pool summary",
    )
    if summary.get("passed") is not True:
        raise ValueError("frozen full-pool summary did not pass")
    if summary.get("schema_version") != FULL_POOL_SUMMARY_SCHEMA_VERSION:
        raise ValueError("frozen full-pool summary has the wrong schema")
    if summary.get("pool_contract_version") != FULL_POOL_CONTRACT_VERSION:
        raise ValueError("frozen full-pool summary has the wrong contract")
    sectors = sorted(pool["sector"].astype(int).unique().tolist())
    if summary.get("sectors") != sectors:
        raise ValueError("frozen full-pool summary sectors differ from table")
    if summary.get("observation_identity_columns") != ["sector", "tic"]:
        raise ValueError("frozen full-pool summary has the wrong identity")
    retained = summary.get("counts", {}).get("retained", {})
    expected_retained = {
        "n_observations": int(len(pool)),
        "n_unique_tics": int(pool["tic"].nunique()),
        "n_multisector_tics": int(pool.groupby("tic")["sector"].nunique().gt(1).sum()),
    }
    if retained != expected_retained:
        raise ValueError("frozen full-pool retained counts differ from table")
    identity = observation_identity_sha256(pool)
    if (
        summary.get("identity_hashes", {}).get("retained_observations_sha256")
        != identity
    ):
        raise ValueError("frozen full-pool identity differs from table")
    suffix = pool_binding.path.suffix.lower()
    output_name = "parquet" if suffix in {".parquet", ".pq"} else "csv"
    output = summary.get("outputs", {}).get(output_name)
    if not isinstance(output, Mapping):
        raise ValueError(f"frozen full-pool summary lacks {output_name} output binding")
    if (
        output.get("sha256") != pool_binding.sha256
        or int(output.get("size_bytes", -1)) != pool_binding.size_bytes
        or int(output.get("n_rows", -1)) != len(pool)
    ):
        raise ValueError("frozen full-pool table binding differs from summary")
    return pool, summary


def _validate_global_bls_authority(
    *,
    bls_binding: ArtifactBinding,
    summary_binding: ArtifactBinding,
    frozen_pool: pd.DataFrame,
    frozen_pool_summary_binding: ArtifactBinding,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if bls_binding.path.suffix.lower() not in {".parquet", ".pq"}:
        raise ValueError("authorized global BLS artifact must be Parquet")
    summary = _load_json(
        summary_binding.path,
        context="global full-pool BLS summary",
    )
    if set(summary) != _GLOBAL_BLS_SUMMARY_FIELDS:
        raise ValueError("global full-pool BLS summary has the wrong fields")
    if summary.get("passed") is not True:
        raise ValueError("global full-pool BLS summary did not pass")
    if summary.get("schema_version") != GLOBAL_BLS_SUMMARY_SCHEMA_VERSION:
        raise ValueError("global full-pool BLS summary has the wrong schema")
    if summary.get("contract_version") != GLOBAL_BLS_CONTRACT_VERSION:
        raise ValueError("global full-pool BLS summary has the wrong contract")
    sectors = sorted(frozen_pool["sector"].astype(int).unique().tolist())
    if summary.get("sectors") != sectors:
        raise ValueError("global full-pool BLS sectors differ from frozen pool")
    if summary.get("observation_identity_columns") != ["sector", "tic"]:
        raise ValueError("global full-pool BLS summary has the wrong identity")
    identity = observation_identity_sha256(frozen_pool)
    pool_binding = summary.get("frozen_pool")
    if not isinstance(pool_binding, Mapping):
        raise ValueError("global BLS summary lacks frozen-pool binding")
    pool_summary = pool_binding.get("summary")
    if not isinstance(pool_summary, Mapping):
        raise ValueError("global BLS summary lacks pool-summary binding")
    if (
        pool_summary.get("sha256") != frozen_pool_summary_binding.sha256
        or int(pool_summary.get("size_bytes", -1))
        != frozen_pool_summary_binding.size_bytes
        or pool_binding.get("contract_version") != FULL_POOL_CONTRACT_VERSION
        or pool_binding.get("observation_identity_sha256") != identity
        or int(pool_binding.get("n_observations", -1)) != len(frozen_pool)
        or int(pool_binding.get("n_unique_tics", -1)) != frozen_pool["tic"].nunique()
    ):
        raise ValueError("global BLS summary binds a different frozen pool")
    bls_contract = summary.get("bls_contract")
    if not isinstance(bls_contract, Mapping):
        raise ValueError("global BLS summary lacks BLS contract")
    if bls_contract.get("apertures") != list(ADP_ONLY_APERTURES):
        raise ValueError("global BLS summary does not bind the exact ADP pair")
    coverage = summary.get("coverage_audit")
    if coverage != {
        "missing_frozen_observations": 0,
        "unexpected_bls_observations": 0,
        "observation_identity_sha256": identity,
    }:
        raise ValueError("global full-pool BLS coverage audit failed")
    output = summary.get("output")
    if not isinstance(output, Mapping):
        raise ValueError("global full-pool BLS summary lacks output binding")
    if (
        output.get("sha256") != bls_binding.sha256
        or int(output.get("size_bytes", -1)) != bls_binding.size_bytes
    ):
        raise ValueError("global full-pool BLS artifact differs from summary")
    bls = read_table(bls_binding.path)
    counts = summary.get("counts")
    if not isinstance(counts, Mapping):
        raise ValueError("global full-pool BLS summary lacks counts")
    if (
        int(output.get("n_rows", -1)) != len(bls)
        or int(counts.get("n_rows", -1)) != len(bls)
        or int(output.get("n_observations", -1)) != len(frozen_pool)
        or int(counts.get("n_observations", -1)) != len(frozen_pool)
    ):
        raise ValueError("global full-pool BLS counts differ from artifacts")
    observed_keys = set(
        _normalize_keys(
            bls.loc[:, ["sector", "tic"]].drop_duplicates(),
            context="global BLS table",
        )
    )
    full_keys = set(_normalize_keys(frozen_pool, context="frozen SSL pool"))
    if observed_keys != full_keys:
        raise ValueError("global full-pool BLS coverage differs from frozen pool")
    return bls, summary


def _resolve_source_path(
    *,
    override: Path | None,
    metadata: Mapping[str, Any],
    context: str,
) -> Path:
    if override is not None:
        return Path(override).expanduser().resolve(strict=True)
    value = metadata.get("path")
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"eligibility summary lacks {context} path")
    return Path(value).expanduser().resolve(strict=True)


def _declared_binding(
    metadata: Mapping[str, Any],
    *,
    context: str,
) -> ArtifactBinding:
    value = metadata.get("path")
    digest = str(metadata.get("sha256", "")).lower()
    try:
        size_bytes = int(metadata.get("size_bytes", -1))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"eligibility summary {context} size is invalid") from exc
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"eligibility summary {context} path is invalid")
    if not pd.Series([digest]).str.fullmatch(_SHA256_PATTERN).all():
        raise ValueError(f"eligibility summary {context} SHA-256 is invalid")
    if size_bytes <= 0:
        raise ValueError(f"eligibility summary {context} size is invalid")
    return ArtifactBinding(
        path=Path(value).expanduser().resolve(),
        size_bytes=size_bytes,
        sha256=digest,
    )


def _csv_payload(table: pd.DataFrame) -> bytes:
    return table.to_csv(index=False, lineterminator="\n").encode("utf-8")


def _json_payload(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            value,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _preflight_immutable(path: Path, payload: bytes) -> None:
    if path.exists() and path.read_bytes() != payload:
        raise FileExistsError(
            f"refusing to replace immutable output with different bytes: {path}"
        )


def _publish_immutable(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace immutable output with different bytes: {path}"
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


def _metadata_for_payload(
    path: Path,
    payload: bytes,
    **extra: Any,
) -> dict[str, Any]:
    return {
        "path": str(path),
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        **extra,
    }


def _partition_from_decisions(
    decisions: pd.DataFrame,
) -> tuple[
    frozenset[ObservationKey], frozenset[ObservationKey], frozenset[ObservationKey]
]:
    full_keys = frozenset(_normalize_keys(decisions, context="eligibility decisions"))
    eligible_keys = frozenset(
        _normalize_keys(
            decisions.loc[decisions["native_model_eligible"]],
            context="eligible observations",
        )
    )
    excluded_keys = frozenset(full_keys - eligible_keys)
    return full_keys, eligible_keys, excluded_keys


def _validate_production_partition(
    *,
    full_keys: frozenset[ObservationKey],
    eligible_keys: frozenset[ObservationKey],
    excluded_keys: frozenset[ObservationKey],
) -> None:
    observed = {
        "n_full": len(full_keys),
        "n_eligible": len(eligible_keys),
        "n_excluded": len(excluded_keys),
        "full_sha256": observation_identity_sha256(full_keys),
        "eligible_sha256": observation_identity_sha256(eligible_keys),
        "excluded_sha256": observation_identity_sha256(excluded_keys),
    }
    expected = {
        "n_full": PRODUCTION_FULL_OBSERVATIONS,
        "n_eligible": PRODUCTION_ELIGIBLE_OBSERVATIONS,
        "n_excluded": PRODUCTION_EXCLUDED_OBSERVATIONS,
        "full_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
        "eligible_sha256": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
        "excluded_sha256": PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    }
    if observed != expected:
        raise ValueError(
            "native/model eligibility partition differs from production lock; "
            f"observed={observed}"
        )
    by_sector = Counter(sector for sector, _tic in excluded_keys)
    if dict(sorted(by_sector.items())) != dict(PRODUCTION_EXCLUDED_BY_SECTOR):
        raise ValueError("native/model exclusions differ from production sector counts")


def _production_bls_contract() -> dict[str, Any]:
    config = approved_a2v1_teacher_bls_config()
    config["source_product_tag"] = FULL_POOL_BLS_SOURCE_PRODUCT_TAG
    return {
        "config": config,
        "config_sha256": bls_config_sha256(config),
        "search_contract_version": (
            PRODUCTION_GLOBAL_BLS_SEARCH_CONTRACT_VERSION
        ),
        "source_product_tag": FULL_POOL_BLS_SOURCE_PRODUCT_TAG,
        "apertures": list(ADP_ONLY_APERTURES),
        "orbitid_reconciliation_contract_version": (
            ORBITID_RECONCILIATION_CONTRACT_VERSION
        ),
        "orbitid_policy_by_sector": {
            str(sector): ORBITID_POLICY_BY_SECTOR[sector]
            for sector in sorted(ORBITID_POLICY_BY_SECTOR)
        },
    }


def _validate_production_source_bindings(
    *,
    pool_binding: ArtifactBinding,
    pool_summary_binding: ArtifactBinding,
    bls_binding: ArtifactBinding,
    bls_summary_binding: ArtifactBinding,
    bls_summary: Mapping[str, Any] | None = None,
) -> None:
    observed = {
        "frozen_pool_csv_sha256": pool_binding.sha256,
        "frozen_pool_summary_sha256": pool_summary_binding.sha256,
        "global_bls_sha256": bls_binding.sha256,
        "global_bls_summary_sha256": bls_summary_binding.sha256,
    }
    expected = {
        "frozen_pool_csv_sha256": PRODUCTION_FROZEN_POOL_CSV_SHA256,
        "frozen_pool_summary_sha256": PRODUCTION_FROZEN_POOL_SUMMARY_SHA256,
        "global_bls_sha256": PRODUCTION_GLOBAL_BLS_SHA256,
        "global_bls_summary_sha256": PRODUCTION_GLOBAL_BLS_SUMMARY_SHA256,
    }
    if observed != expected:
        raise ValueError(
            "native/model eligibility source artifacts differ from production "
            f"lock; observed={observed}"
        )
    if bls_summary is not None and bls_summary.get("bls_contract") != (
        _production_bls_contract()
    ):
        raise ValueError(
            "native/model eligibility BLS contract differs from production lock"
        )


def write_native_model_eligibility(
    *,
    pool_path: Path,
    pool_summary_path: Path,
    bls_path: Path,
    bls_summary_path: Path,
    exclusions_path: Path,
    summary_path: Path,
    anchor_aperture: str = NATIVE_MODEL_ELIGIBILITY_ANCHOR_APERTURE,
    production_lock: bool = False,
) -> EligibilityAuthority:
    """Derive and immutably publish the native/model eligibility authority."""

    pool_binding = _binding(pool_path)
    pool_summary_binding = _binding(pool_summary_path)
    bls_binding = _binding(bls_path)
    bls_summary_binding = _binding(bls_summary_path)
    exclusions_path = Path(exclusions_path).expanduser().resolve()
    summary_path = Path(summary_path).expanduser().resolve()
    input_paths = {
        pool_binding.path,
        pool_summary_binding.path,
        bls_binding.path,
        bls_summary_binding.path,
    }
    if (
        exclusions_path == summary_path
        or {
            exclusions_path,
            summary_path,
        }
        & input_paths
    ):
        raise ValueError("eligibility input and output paths must be distinct")

    pool, pool_summary = _validate_frozen_pool_authority(
        pool_binding=pool_binding,
        summary_binding=pool_summary_binding,
    )
    bls, bls_summary = _validate_global_bls_authority(
        bls_binding=bls_binding,
        summary_binding=bls_summary_binding,
        frozen_pool=pool,
        frozen_pool_summary_binding=pool_summary_binding,
    )
    if production_lock:
        _validate_production_source_bindings(
            pool_binding=pool_binding,
            pool_summary_binding=pool_summary_binding,
            bls_binding=bls_binding,
            bls_summary_binding=bls_summary_binding,
            bls_summary=bls_summary,
        )
    decisions = derive_anchor_eligibility(
        bls,
        frozen_pool=pool,
        anchor_aperture=anchor_aperture,
    )
    full_keys, eligible_keys, excluded_keys = _partition_from_decisions(decisions)
    if full_keys != eligible_keys | excluded_keys:
        raise AssertionError("eligibility union does not reconstruct full pool")
    if eligible_keys & excluded_keys:
        raise AssertionError("eligibility partition is not disjoint")
    if production_lock:
        _validate_production_partition(
            full_keys=full_keys,
            eligible_keys=eligible_keys,
            excluded_keys=excluded_keys,
        )

    exclusions = decisions.loc[
        ~decisions["native_model_eligible"],
        [
            "observation_id",
            "sector",
            "tic",
            "aperture",
            "peak_rank",
            "status",
            "exclusion_reason",
        ],
    ].copy()
    exclusions.insert(
        0,
        "contract_version",
        NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
    )
    exclusions = exclusions.rename(
        columns={
            "aperture": "anchor_aperture",
            "status": "bls_status",
        }
    ).loc[:, list(ELIGIBILITY_EXCLUSION_COLUMNS)]
    exclusions = exclusions.sort_values(
        ["sector", "tic"],
        kind="stable",
    ).reset_index(drop=True)
    exclusions_payload = _csv_payload(exclusions)
    full_identity = observation_identity_sha256(full_keys)
    eligible_identity = observation_identity_sha256(eligible_keys)
    excluded_identity = observation_identity_sha256(excluded_keys)

    by_sector: dict[str, dict[str, int]] = {}
    for sector in sorted(pool["sector"].astype(int).unique()):
        sector_full = {(s, tic) for s, tic in full_keys if s == sector}
        sector_eligible = {(s, tic) for s, tic in eligible_keys if s == sector}
        sector_excluded = {(s, tic) for s, tic in excluded_keys if s == sector}
        by_sector[str(sector)] = {
            "full": len(sector_full),
            "eligible": len(sector_eligible),
            "excluded": len(sector_excluded),
        }
    summary: dict[str, Any] = {
        "passed": True,
        "schema_version": (NATIVE_MODEL_ELIGIBILITY_SUMMARY_SCHEMA_VERSION),
        "contract_version": NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
        "release_binding": NATIVE_MODEL_ELIGIBILITY_RELEASE_BINDING,
        "sectors": sorted(pool["sector"].astype(int).unique().tolist()),
        "observation_identity_columns": ["sector", "tic"],
        "anchor_policy": {
            "anchor_aperture": str(anchor_aperture),
            "accepted_peak_ranks": [0, 1],
            "eligible_peak_rank": 1,
            "eligible_status": "ok",
            "requires_finite_positive_period": True,
            "requires_finite_epoch": True,
            "requires_finite_positive_duration": True,
            "exclusion_reason": NATIVE_MODEL_ELIGIBILITY_EXCLUSION_REASON,
            "authority_source": ("authorized global BLS Parquet anchor rows only"),
        },
        "source_authorities": {
            "frozen_pool": {
                "artifact": pool_binding.metadata(),
                "summary": pool_summary_binding.metadata(),
                "contract_version": FULL_POOL_CONTRACT_VERSION,
                "summary_schema_version": FULL_POOL_SUMMARY_SCHEMA_VERSION,
                "summary_declared_output_hash_equal": True,
                "summary_identity_hash_equal": True,
            },
            "global_bls": {
                "artifact": bls_binding.metadata(),
                "summary": bls_summary_binding.metadata(),
                "contract_version": GLOBAL_BLS_CONTRACT_VERSION,
                "summary_schema_version": GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
                "summary_declared_output_hash_equal": True,
                "frozen_pool_binding_equal": True,
            },
        },
        "counts": {
            "full": {
                "n_observations": len(full_keys),
                "n_unique_tics": int(pool["tic"].nunique()),
            },
            "eligible": {
                "n_observations": len(eligible_keys),
                "n_unique_tics": len({tic for _sector, tic in eligible_keys}),
            },
            "excluded": {
                "n_observations": len(excluded_keys),
                "n_unique_tics": len({tic for _sector, tic in excluded_keys}),
            },
        },
        "by_sector": by_sector,
        "by_reason": {
            str(reason): int(count)
            for reason, count in sorted(
                Counter(exclusions["exclusion_reason"].astype(str)).items()
            )
        },
        "identity_hashes": {
            "full_observations_sha256": full_identity,
            "eligible_observations_sha256": eligible_identity,
            "excluded_observations_sha256": excluded_identity,
        },
        "partition_audit": {
            "eligible_excluded_disjoint": not bool(eligible_keys & excluded_keys),
            "eligible_excluded_union_equals_full": (eligible_keys | excluded_keys)
            == full_keys,
            "count_arithmetic_equal": (
                len(eligible_keys) + len(excluded_keys) == len(full_keys)
            ),
            "excluded_rows_equal_exclusions_csv": (
                len(excluded_keys) == len(exclusions)
            ),
            "production_lock_requested": bool(production_lock),
            "production_lock_passed": bool(production_lock),
        },
        "data_usage_audit": {
            "labels_consumed": False,
            "injections_consumed": False,
            "raw_flux_errors_define_eligibility": False,
            "two_aperture_status_row_counts_define_eligibility": False,
        },
        "outputs": {
            "exclusions": _metadata_for_payload(
                exclusions_path,
                exclusions_payload,
                n_rows=len(exclusions),
                columns=list(ELIGIBILITY_EXCLUSION_COLUMNS),
            ),
        },
    }
    summary_payload = _json_payload(summary)
    _preflight_immutable(exclusions_path, exclusions_payload)
    _preflight_immutable(summary_path, summary_payload)
    _publish_immutable(exclusions_path, exclusions_payload)
    _publish_immutable(summary_path, summary_payload)
    # The validating reload below intentionally rederives from the global BLS
    # authority.  Release the first multi-million-row frame before opening it
    # again so publication verification does not double peak memory.
    del bls, bls_summary, decisions, pool, pool_summary
    return load_native_model_eligibility(
        exclusions_path,
        summary_path,
        pool_path=pool_binding.path,
        pool_summary_path=pool_summary_binding.path,
        bls_path=bls_binding.path,
        bls_summary_path=bls_summary_binding.path,
        production_lock=production_lock,
    )


def _validate_exclusions(
    exclusions: pd.DataFrame,
    *,
    anchor_aperture: str,
) -> pd.DataFrame:
    if tuple(exclusions.columns) != ELIGIBILITY_EXCLUSION_COLUMNS:
        raise ValueError("eligibility exclusions columns/order differ from contract")
    checked = exclusions.copy()
    checked["sector"] = _positive_integer_series(
        checked["sector"],
        name="exclusion sector",
    )
    checked["tic"] = _positive_integer_series(
        checked["tic"],
        name="exclusion tic",
    )
    checked["peak_rank"] = _integer_series(
        checked["peak_rank"],
        name="exclusion peak_rank",
    )
    if checked.duplicated(["sector", "tic"]).any():
        raise ValueError("eligibility exclusions contain duplicate keys")
    if (
        not checked["contract_version"]
        .astype(str)
        .eq(NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION)
        .all()
    ):
        raise ValueError("eligibility exclusions have the wrong contract")
    if not checked["anchor_aperture"].astype(str).eq(anchor_aperture).all():
        raise ValueError("eligibility exclusions have the wrong aperture")
    if (
        not checked["exclusion_reason"]
        .astype(str)
        .eq(NATIVE_MODEL_ELIGIBILITY_EXCLUSION_REASON)
        .all()
    ):
        raise ValueError("eligibility exclusions have the wrong reason")
    expected_ids = [
        f"s{int(sector):04d}-tic{int(tic):016d}"
        for sector, tic in zip(
            checked["sector"],
            checked["tic"],
            strict=True,
        )
    ]
    if checked["observation_id"].astype(str).tolist() != expected_ids:
        raise ValueError("eligibility exclusion observation IDs are invalid")
    ordered = checked.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )
    if not checked.reset_index(drop=True).equals(ordered):
        raise ValueError("eligibility exclusions must be sorted by sector/tic")
    return checked


def load_native_model_eligibility(
    exclusions_path: Path,
    summary_path: Path,
    *,
    pool_path: Path | None = None,
    pool_summary_path: Path | None = None,
    bls_path: Path | None = None,
    bls_summary_path: Path | None = None,
    production_lock: bool = False,
    rederive_from_bls: bool = True,
) -> EligibilityAuthority:
    """Load and fully rederive a native/model eligibility authority.

    Source overrides support byte-identical staging at a new location.  When
    omitted, paths recorded in the eligibility summary are used.  Native shard
    jobs may set ``rederive_from_bls=False``: this still verifies the pool,
    eligibility files, exact complement, counts, identities, and the BLS
    metadata frozen inside the eligibility release, but avoids re-reading the
    multi-million-row BLS Parquet in every shard process.
    """

    exclusions_binding = _binding(exclusions_path)
    summary_binding = _binding(summary_path)
    summary = _load_json(
        summary_binding.path,
        context="native/model eligibility summary",
    )
    if set(summary) != _SUMMARY_TOP_LEVEL_FIELDS:
        raise ValueError("native/model eligibility summary has wrong fields")
    if summary.get("passed") is not True:
        raise ValueError("native/model eligibility summary did not pass")
    if (
        summary.get("schema_version") != NATIVE_MODEL_ELIGIBILITY_SUMMARY_SCHEMA_VERSION
        or summary.get("contract_version") != NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION
        or summary.get("release_binding") != NATIVE_MODEL_ELIGIBILITY_RELEASE_BINDING
    ):
        raise ValueError("native/model eligibility summary has wrong contract")
    if summary.get("observation_identity_columns") != ["sector", "tic"]:
        raise ValueError("native/model eligibility summary has wrong identity")
    anchor_policy = summary.get("anchor_policy")
    if not isinstance(anchor_policy, Mapping):
        raise ValueError("eligibility summary lacks anchor policy")
    anchor_aperture = str(anchor_policy.get("anchor_aperture", ""))
    expected_anchor_policy = {
        "anchor_aperture": anchor_aperture,
        "accepted_peak_ranks": [0, 1],
        "eligible_peak_rank": 1,
        "eligible_status": "ok",
        "requires_finite_positive_period": True,
        "requires_finite_epoch": True,
        "requires_finite_positive_duration": True,
        "exclusion_reason": NATIVE_MODEL_ELIGIBILITY_EXCLUSION_REASON,
        "authority_source": "authorized global BLS Parquet anchor rows only",
    }
    if dict(anchor_policy) != expected_anchor_policy or not anchor_aperture:
        raise ValueError("native/model eligibility anchor policy is invalid")
    data_usage = summary.get("data_usage_audit")
    if data_usage != {
        "labels_consumed": False,
        "injections_consumed": False,
        "raw_flux_errors_define_eligibility": False,
        "two_aperture_status_row_counts_define_eligibility": False,
    }:
        raise ValueError("native/model eligibility data-usage audit failed")

    sources = summary.get("source_authorities")
    if not isinstance(sources, Mapping):
        raise ValueError("eligibility summary lacks source authorities")
    pool_meta = sources.get("frozen_pool")
    bls_meta = sources.get("global_bls")
    if not isinstance(pool_meta, Mapping) or not isinstance(bls_meta, Mapping):
        raise ValueError("eligibility summary source authorities are invalid")
    expected_pool_metadata = {
        "contract_version": FULL_POOL_CONTRACT_VERSION,
        "summary_schema_version": FULL_POOL_SUMMARY_SCHEMA_VERSION,
        "summary_declared_output_hash_equal": True,
        "summary_identity_hash_equal": True,
    }
    if any(
        pool_meta.get(key) != value for key, value in expected_pool_metadata.items()
    ):
        raise ValueError("eligibility summary frozen-pool contract is invalid")
    expected_bls_metadata = {
        "contract_version": GLOBAL_BLS_CONTRACT_VERSION,
        "summary_schema_version": GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
        "summary_declared_output_hash_equal": True,
        "frozen_pool_binding_equal": True,
    }
    if any(bls_meta.get(key) != value for key, value in expected_bls_metadata.items()):
        raise ValueError("eligibility summary global-BLS contract is invalid")
    pool_artifact_meta = pool_meta.get("artifact")
    pool_summary_meta = pool_meta.get("summary")
    bls_artifact_meta = bls_meta.get("artifact")
    bls_summary_meta = bls_meta.get("summary")
    if not all(
        isinstance(value, Mapping)
        for value in (
            pool_artifact_meta,
            pool_summary_meta,
            bls_artifact_meta,
            bls_summary_meta,
        )
    ):
        raise ValueError("eligibility summary source bindings are invalid")
    resolved_pool_path = _resolve_source_path(
        override=pool_path,
        metadata=pool_artifact_meta,
        context="frozen-pool artifact",
    )
    resolved_pool_summary_path = _resolve_source_path(
        override=pool_summary_path,
        metadata=pool_summary_meta,
        context="frozen-pool summary",
    )
    pool_binding = _binding(resolved_pool_path)
    pool_summary_binding = _binding(resolved_pool_summary_path)
    for context, observed, declared in (
        ("frozen-pool artifact", pool_binding, pool_artifact_meta),
        ("frozen-pool summary", pool_summary_binding, pool_summary_meta),
    ):
        if observed.sha256 != declared.get("sha256") or observed.size_bytes != int(
            declared.get("size_bytes", -1)
        ):
            raise ValueError(f"{context} differs from native/model eligibility summary")

    pool, _pool_summary = _validate_frozen_pool_authority(
        pool_binding=pool_binding,
        summary_binding=pool_summary_binding,
    )
    loaded_bls_summary: Mapping[str, Any] | None = None
    if rederive_from_bls:
        resolved_bls_path = _resolve_source_path(
            override=bls_path,
            metadata=bls_artifact_meta,
            context="global-BLS artifact",
        )
        resolved_bls_summary_path = _resolve_source_path(
            override=bls_summary_path,
            metadata=bls_summary_meta,
            context="global-BLS summary",
        )
        bls_binding = _binding(resolved_bls_path)
        bls_summary_binding = _binding(resolved_bls_summary_path)
        for context, observed, declared in (
            ("global-BLS artifact", bls_binding, bls_artifact_meta),
            ("global-BLS summary", bls_summary_binding, bls_summary_meta),
        ):
            if observed.sha256 != declared.get("sha256") or observed.size_bytes != int(
                declared.get("size_bytes", -1)
            ):
                raise ValueError(
                    f"{context} differs from native/model eligibility summary"
                )
        bls, loaded_bls_summary = _validate_global_bls_authority(
            bls_binding=bls_binding,
            summary_binding=bls_summary_binding,
            frozen_pool=pool,
            frozen_pool_summary_binding=pool_summary_binding,
        )
        decisions = derive_anchor_eligibility(
            bls,
            frozen_pool=pool,
            anchor_aperture=anchor_aperture,
        )
        full_keys, eligible_keys, excluded_keys = _partition_from_decisions(decisions)
    else:
        bls_binding = _declared_binding(
            bls_artifact_meta,
            context="global-BLS artifact",
        )
        bls_summary_binding = _declared_binding(
            bls_summary_meta,
            context="global-BLS summary",
        )
        decisions = None
        full_keys = frozenset(_normalize_keys(pool, context="frozen SSL pool"))
        eligible_keys = frozenset()
        excluded_keys = frozenset()

    outputs = summary.get("outputs")
    exclusion_meta = outputs.get("exclusions") if isinstance(outputs, Mapping) else None
    if not isinstance(exclusion_meta, Mapping):
        raise ValueError("eligibility summary lacks exclusions output")
    if (
        exclusion_meta.get("sha256") != exclusions_binding.sha256
        or int(exclusion_meta.get("size_bytes", -1)) != exclusions_binding.size_bytes
    ):
        raise ValueError("eligibility exclusions differ from summary")
    try:
        exclusions = pd.read_csv(
            io.BytesIO(exclusions_binding.path.read_bytes()),
            dtype={
                "contract_version": str,
                "observation_id": str,
                "anchor_aperture": str,
                "bls_status": str,
                "exclusion_reason": str,
            },
            keep_default_na=False,
        )
    except Exception as exc:
        raise ValueError("unable to read eligibility exclusions CSV") from exc
    exclusions = _validate_exclusions(
        exclusions,
        anchor_aperture=anchor_aperture,
    )
    excluded_csv_keys = frozenset(
        _normalize_keys(exclusions, context="eligibility exclusions")
    )
    if not excluded_csv_keys <= full_keys:
        extra = sorted(excluded_csv_keys - full_keys)[:10]
        raise ValueError(
            f"eligibility exclusions contain keys outside the full pool; first={extra}"
        )
    if decisions is not None:
        expected_exclusions = decisions.loc[
            ~decisions["native_model_eligible"],
            [
                "observation_id",
                "sector",
                "tic",
                "aperture",
                "peak_rank",
                "status",
                "exclusion_reason",
            ],
        ].copy()
        expected_exclusions.insert(
            0,
            "contract_version",
            NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
        )
        expected_exclusions = expected_exclusions.rename(
            columns={"aperture": "anchor_aperture", "status": "bls_status"}
        ).loc[:, list(ELIGIBILITY_EXCLUSION_COLUMNS)]
        expected_exclusions = expected_exclusions.sort_values(
            ["sector", "tic"],
            kind="stable",
        ).reset_index(drop=True)
        pd.testing.assert_frame_equal(
            exclusions,
            expected_exclusions,
            check_dtype=False,
            check_like=False,
        )
        if excluded_csv_keys != excluded_keys:
            raise ValueError("eligibility exclusions do not equal derived complement")
    else:
        excluded_keys = excluded_csv_keys
        eligible_keys = frozenset(full_keys - excluded_keys)

    partition_audit = summary.get("partition_audit")
    declared_production_lock = (
        partition_audit.get("production_lock_requested")
        if isinstance(partition_audit, Mapping)
        else None
    )
    declared_production_pass = (
        partition_audit.get("production_lock_passed")
        if isinstance(partition_audit, Mapping)
        else None
    )
    if (
        not isinstance(declared_production_lock, bool)
        or not isinstance(declared_production_pass, bool)
        or declared_production_pass != declared_production_lock
    ):
        raise ValueError("eligibility summary production-lock audit is invalid")
    if production_lock and not declared_production_lock:
        raise ValueError("eligibility summary was not production locked")
    if production_lock or declared_production_lock:
        _validate_production_source_bindings(
            pool_binding=pool_binding,
            pool_summary_binding=pool_summary_binding,
            bls_binding=bls_binding,
            bls_summary_binding=bls_summary_binding,
            bls_summary=loaded_bls_summary,
        )
        _validate_production_partition(
            full_keys=full_keys,
            eligible_keys=eligible_keys,
            excluded_keys=excluded_keys,
        )

    full_identity = observation_identity_sha256(full_keys)
    eligible_identity = observation_identity_sha256(eligible_keys)
    excluded_identity = observation_identity_sha256(excluded_keys)
    expected_counts = {
        "full": {
            "n_observations": len(full_keys),
            "n_unique_tics": len({tic for _sector, tic in full_keys}),
        },
        "eligible": {
            "n_observations": len(eligible_keys),
            "n_unique_tics": len({tic for _sector, tic in eligible_keys}),
        },
        "excluded": {
            "n_observations": len(excluded_keys),
            "n_unique_tics": len({tic for _sector, tic in excluded_keys}),
        },
    }
    if summary.get("counts") != expected_counts:
        raise ValueError("eligibility summary counts differ from partition")
    expected_by_sector: dict[str, dict[str, int]] = {}
    for sector in sorted({sector for sector, _tic in full_keys}):
        expected_by_sector[str(sector)] = {
            "full": sum(key[0] == sector for key in full_keys),
            "eligible": sum(key[0] == sector for key in eligible_keys),
            "excluded": sum(key[0] == sector for key in excluded_keys),
        }
    if summary.get("by_sector") != expected_by_sector:
        raise ValueError("eligibility summary sector counts differ")
    expected_by_reason = {
        str(reason): int(count)
        for reason, count in sorted(
            Counter(exclusions["exclusion_reason"].astype(str)).items()
        )
    }
    if summary.get("by_reason") != expected_by_reason:
        raise ValueError("eligibility summary reason counts differ")
    if summary.get("identity_hashes") != {
        "full_observations_sha256": full_identity,
        "eligible_observations_sha256": eligible_identity,
        "excluded_observations_sha256": excluded_identity,
    }:
        raise ValueError("eligibility summary identity hashes differ")
    expected_partition_audit = {
        "eligible_excluded_disjoint": True,
        "eligible_excluded_union_equals_full": True,
        "count_arithmetic_equal": True,
        "excluded_rows_equal_exclusions_csv": True,
        "production_lock_requested": declared_production_lock,
        "production_lock_passed": declared_production_pass,
    }
    if partition_audit != expected_partition_audit:
        raise ValueError("eligibility summary partition audit failed")
    if summary.get("sectors") != sorted({sector for sector, _tic in full_keys}):
        raise ValueError("eligibility summary sectors differ from partition")
    if int(exclusion_meta.get("n_rows", -1)) != len(exclusions):
        raise ValueError("eligibility exclusion row count differs from summary")
    if exclusion_meta.get("columns") != list(ELIGIBILITY_EXCLUSION_COLUMNS):
        raise ValueError("eligibility exclusion schema differs from summary")

    bindings = {
        "frozen_pool": pool_binding,
        "frozen_pool_summary": pool_summary_binding,
        "global_bls": bls_binding,
        "global_bls_summary": bls_summary_binding,
        "exclusions": exclusions_binding,
        "eligibility_summary": summary_binding,
    }
    return EligibilityAuthority(
        contract_version=NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
        release_binding=NATIVE_MODEL_ELIGIBILITY_RELEASE_BINDING,
        anchor_aperture=anchor_aperture,
        full_keys=full_keys,
        eligible_keys=eligible_keys,
        excluded_keys=excluded_keys,
        n_full=len(full_keys),
        n_eligible=len(eligible_keys),
        n_excluded=len(excluded_keys),
        full_observation_identity_sha256=full_identity,
        eligible_observation_identity_sha256=eligible_identity,
        excluded_observation_identity_sha256=excluded_identity,
        bindings=bindings,
        exclusions=exclusions,
        summary=summary,
    )


__all__ = [
    "ArtifactBinding",
    "ELIGIBILITY_EXCLUSION_COLUMNS",
    "EligibilityAuthority",
    "NATIVE_MODEL_ELIGIBILITY_ANCHOR_APERTURE",
    "NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION",
    "NATIVE_MODEL_ELIGIBILITY_EXCLUSION_REASON",
    "NATIVE_MODEL_ELIGIBILITY_RELEASE_BINDING",
    "NATIVE_MODEL_ELIGIBILITY_SUMMARY_SCHEMA_VERSION",
    "PRODUCTION_ELIGIBLE_IDENTITY_SHA256",
    "PRODUCTION_ELIGIBLE_OBSERVATIONS",
    "PRODUCTION_EXCLUDED_BY_SECTOR",
    "PRODUCTION_EXCLUDED_IDENTITY_SHA256",
    "PRODUCTION_EXCLUDED_OBSERVATIONS",
    "PRODUCTION_FROZEN_POOL_CSV_SHA256",
    "PRODUCTION_FULL_IDENTITY_SHA256",
    "PRODUCTION_FULL_OBSERVATIONS",
    "derive_anchor_eligibility",
    "load_native_model_eligibility",
    "observation_identity_sha256",
    "write_native_model_eligibility",
]
