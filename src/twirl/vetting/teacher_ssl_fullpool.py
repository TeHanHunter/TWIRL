"""Fold-local self-supervision over the broad S56--S62 light-curve pool.

This module is deliberately separate from :mod:`teacher_ssl_training`.  The
older module is the frozen, label-table-backed Teacher v4-SSL pilot.  The
full-pool contract instead starts from a rank-one BLS observation table and a
separate observation-keyed native registry.  Human labels are neither required
nor allowed in the resulting SSL registry.

The registry builder applies three host-level rules before any tensor is
constructed:

* every TIC in the frozen Teacher-v3 fixed test is excluded;
* every TIC in the prospective S63 inventory is excluded, including its
  earlier-sector observations;
* development TICs retain their frozen Teacher-v3 fold, while previously
  unseen TICs have no held fold and therefore appear in all five pretraining
  pools.

Each encoder job leaves the matching frozen labeled-development TICs out of
pretraining while retaining all leakage-safe previously unseen TICs.  The
runner writes one atomic resume checkpoint per epoch and never performs the
pilot's quadratic neighbor probe or all-pool embedding export, both of which
are inappropriate for roughly 200,000 observations.
"""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import subprocess
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    MODEL_VERSION,
    HarmonicModelConfig,
    build_harmonic_cnn,
)
from twirl.vetting.harmonic_dataset import HarmonicNativeDataset
from twirl.vetting.harmonic_inputs import MODEL_INPUT_CONTRACT_VERSION
from twirl.vetting.harmonic_ssl import (
    HARMONIC_SSL_CONTRACT_VERSION,
    EventPreservingAugmentationConfig,
    VICRegConfig,
    augment_ssl_batch,
    vicreg_loss,
)
from twirl.vetting.harmonic_training import _loader, _to_device
from twirl.vetting.ssl_full_pool import (
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
    POOL_COLUMNS,
)
from twirl.vetting.ssl_full_pool_bls import (
    GLOBAL_BLS_CONTRACT_VERSION,
    GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
)
from twirl.vetting.ssl_full_pool_eligibility import (
    EligibilityAuthority,
    PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    PRODUCTION_EXCLUDED_OBSERVATIONS,
    PRODUCTION_FULL_IDENTITY_SHA256,
    PRODUCTION_FULL_OBSERVATIONS,
    derive_anchor_eligibility,
    load_native_model_eligibility,
)
from twirl.vetting.ssl_full_pool_native import (
    FULL_POOL_NATIVE_CONTRACT_VERSION,
    load_full_pool_native_registry_release,
)
from twirl.vetting.ssl_full_pool_numeric import (
    MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT,
    validate_numeric_gate_release,
)
from twirl.vetting.teacher_native_registry import file_sha256, read_table
from twirl.vetting.teacher_split_registry import validate_tic_split_assignments


FULLPOOL_SSL_REGISTRY_SCHEMA = "twirl_teacher_ssl_fullpool_registry_v2"
FULLPOOL_SSL_REGISTRY_SUMMARY_SCHEMA = (
    "twirl_teacher_ssl_fullpool_registry_summary_v2"
)
FULLPOOL_SSL_SELECTION_SCHEMA = "twirl_teacher_ssl_fullpool_fold_selection_v3"
FULLPOOL_SSL_RUN_CONTRACT_SCHEMA = "twirl_teacher_ssl_fullpool_fold_run_v3"
FULLPOOL_SSL_RESUME_SCHEMA = "twirl_teacher_ssl_fullpool_resume_v3"
FULLPOOL_SSL_CHECKPOINT_SCHEMA = "twirl_teacher_ssl_fullpool_checkpoint_v3"
FULLPOOL_SSL_SUMMARY_SCHEMA = "twirl_teacher_ssl_fullpool_fold_summary_v3"
FULLPOOL_SSL_RUN_ID = (
    "teacher_ssl_fullpool_v3_s56_s62_a2v1_current_adp_bls_eligible_"
    "effective_quality_mask_v1"
)
FULLPOOL_SSL_ENCODER_NAME = (
    "teacher_ssl_fullpool_v3_effective_quality_mask_v1"
)
FULLPOOL_SSL_MODEL_FACING_NAME = "Teacher v4-SSL full-pool"
FULLPOOL_SSL_PROFILE = "shape_plus_periodogram_bls"
FULLPOOL_SSL_CHECKPOINT_NAMESPACE = (
    "twirl_teacher_ssl_fullpool_v3_s56_s62_a2v1_current_adp_bls_eligible_"
    "effective_quality_mask_v1"
)
FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA = (
    "twirl_teacher_ssl_fullpool_training_authority_v3"
)
FULLPOOL_SSL_SECTORS: tuple[int, ...] = tuple(range(56, 63))
FULLPOOL_SSL_N_FOLDS = 5
FULLPOOL_SSL_DEFAULT_TRAINING_SEED = 560064
FULLPOOL_SSL_ANCHOR_APERTURE = "DET_FLUX_ADP_SML"
_SHA256_PATTERN = r"^[0-9a-f]{64}$"

FULLPOOL_SSL_REGISTRY_COLUMNS: tuple[str, ...] = (
    "registry_schema_version",
    "ssl_observation_id",
    "sector",
    "tic",
    "period_d",
    "t0_bjd",
    "duration_min",
    "bls_status",
    "native_h5_path",
    "native_group_path",
    "native_h5_sha256",
    "native_contract_version",
    "is_injected_row",
    "fixed_test_member",
    "reserved_prospective_member",
    "ssl_pool_include",
    "ssl_pool_exclusion_reason",
    "ssl_held_out_fold",
    "fold_assignment_source",
)

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


def _canonical_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("canonical mappings require string keys")
        return {
            key: _canonical_value(item)
            for key, item in sorted(value.items(), key=lambda pair: pair[0])
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return _canonical_value(value.item())
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return value
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError(f"unsupported canonical payload type: {type(value).__name__}")


def _canonical_sha256(payload: Any) -> str:
    encoded = json.dumps(
        _canonical_value(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _positive_integer_series(values: pd.Series, *, name: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    invalid = (
        ~np.isfinite(numeric)
        | (numeric <= 0)
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
    invalid = ~np.isfinite(numeric) | ~np.equal(numeric, np.rint(numeric))
    if invalid.any():
        examples = values.loc[invalid].head(5).tolist()
        raise ValueError(f"{name} must contain finite integers; first={examples}")
    return pd.Series(
        np.rint(numeric).astype(np.int64),
        index=values.index,
        name=name,
    )


def _strict_bool_series(values: pd.Series, *, name: str) -> pd.Series:
    true_values = {"1", "true", "t", "yes", "y"}
    false_values = {"0", "false", "f", "no", "n"}
    parsed: list[bool] = []
    for value in values.tolist():
        if isinstance(value, (bool, np.bool_)):
            parsed.append(bool(value))
            continue
        if isinstance(value, (int, np.integer)) and int(value) in (0, 1):
            parsed.append(bool(value))
            continue
        if isinstance(value, str):
            text = value.strip().lower()
            if text in true_values:
                parsed.append(True)
                continue
            if text in false_values:
                parsed.append(False)
                continue
        raise ValueError(f"{name} contains an invalid boolean: {value!r}")
    return pd.Series(parsed, index=values.index, dtype=bool, name=name)


def _nonblank_text(values: pd.Series, *, name: str) -> pd.Series:
    normalized = values.fillna("").astype(str).str.strip()
    if normalized.eq("").any():
        examples = values.loc[normalized.eq("")].head(5).tolist()
        raise ValueError(f"{name} contains blank values; first={examples}")
    return normalized.rename(name)


def _frame_sha256(frame: pd.DataFrame, columns: Sequence[str]) -> str:
    records = frame.loc[:, list(columns)].to_dict(orient="records")
    return _canonical_sha256(records)


def _pool_identity_sha256(rows: pd.DataFrame) -> str:
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


def load_frozen_ssl_full_pool(
    *,
    pool_path: Path,
    summary_path: Path,
    validate_allowlists: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate the preregistered LC pool, summary, and sector allowlists."""

    pool_path = Path(pool_path).expanduser().resolve()
    summary_path = Path(summary_path).expanduser().resolve()
    pool = read_table(pool_path)
    missing = sorted(set(POOL_COLUMNS) - set(pool.columns))
    if missing:
        raise KeyError(f"frozen SSL pool lacks columns: {missing}")
    if pool.empty:
        raise ValueError("frozen SSL pool is empty")
    pool = pool.loc[:, list(POOL_COLUMNS)].copy()
    pool["sector"] = _positive_integer_series(pool["sector"], name="sector")
    pool["tic"] = _positive_integer_series(pool["tic"], name="tic")
    if pool.duplicated(["sector", "tic"]).any():
        raise ValueError("frozen SSL pool has duplicate observation keys")
    if pool["observation_id"].fillna("").astype(str).str.strip().eq("").any():
        raise ValueError("frozen SSL pool has blank observation IDs")
    if pool["observation_id"].astype(str).duplicated().any():
        raise ValueError("frozen SSL pool observation IDs are not unique")
    if not pool["pool_contract_version"].astype(str).eq(
        FULL_POOL_CONTRACT_VERSION
    ).all():
        raise ValueError("frozen SSL pool has the wrong contract version")
    observed_sectors = tuple(sorted(pool["sector"].unique().astype(int)))
    if observed_sectors != FULLPOOL_SSL_SECTORS:
        raise ValueError(
            "frozen SSL pool sectors do not match S56--S62: "
            f"{observed_sectors}"
        )

    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid frozen SSL pool summary: {summary_path}") from exc
    if not isinstance(summary, dict):
        raise ValueError("frozen SSL pool summary must be a JSON object")
    if summary.get("schema_version") != FULL_POOL_SUMMARY_SCHEMA_VERSION:
        raise ValueError("frozen SSL pool summary has the wrong schema")
    if summary.get("pool_contract_version") != FULL_POOL_CONTRACT_VERSION:
        raise ValueError("frozen SSL pool summary has the wrong pool contract")
    if summary.get("sectors") != list(FULLPOOL_SSL_SECTORS):
        raise ValueError("frozen SSL pool summary has the wrong sectors")
    if summary.get("leakage_audit") != {
        "fixed_test_observations_retained": 0,
        "s63_reserved_observations_retained": 0,
    }:
        raise ValueError("frozen SSL pool summary did not pass leakage audit")
    retained = summary.get("counts", {}).get("retained", {})
    expected_counts = {
        "n_observations": int(len(pool)),
        "n_unique_tics": int(pool["tic"].nunique()),
        "n_multisector_tics": int(
            pool.groupby("tic")["sector"].nunique().gt(1).sum()
        ),
    }
    if retained != expected_counts:
        raise ValueError("frozen SSL pool retained counts do not match table")
    identity = summary.get("identity_hashes", {}).get(
        "retained_observations_sha256"
    )
    if identity != _pool_identity_sha256(pool):
        raise ValueError("frozen SSL pool identity hash does not match table")

    suffix = pool_path.suffix.lower()
    output_key = "parquet" if suffix in {".parquet", ".pq"} else "csv"
    declared_output = summary.get("outputs", {}).get(output_key)
    if not isinstance(declared_output, dict):
        raise ValueError("frozen SSL pool summary lacks selected table metadata")
    if declared_output.get("sha256") != file_sha256(pool_path):
        raise ValueError("frozen SSL pool file hash does not match summary")
    if int(declared_output.get("size_bytes", -1)) != pool_path.stat().st_size:
        raise ValueError("frozen SSL pool file size does not match summary")
    if int(declared_output.get("n_rows", -1)) != len(pool):
        raise ValueError("frozen SSL pool row count does not match summary")

    if validate_allowlists:
        declared_allowlists = summary.get("outputs", {}).get(
            "sector_allowlists"
        )
        if not isinstance(declared_allowlists, dict):
            raise ValueError("frozen SSL pool summary lacks sector allowlists")
        for sector in FULLPOOL_SSL_SECTORS:
            metadata = declared_allowlists.get(str(sector))
            if not isinstance(metadata, dict):
                raise ValueError(f"frozen SSL pool lacks S{sector} allowlist")
            recorded = Path(str(metadata.get("path", ""))).expanduser()
            adjacent = (
                summary_path.parent
                / "allowlists"
                / f"s{sector}_tics.csv"
            )
            allowlist_path = (
                recorded.resolve()
                if recorded.is_file()
                else adjacent.resolve()
            )
            if not allowlist_path.is_file():
                raise FileNotFoundError(
                    f"frozen SSL pool S{sector} allowlist is missing"
                )
            if metadata.get("sha256") != file_sha256(allowlist_path):
                raise ValueError(
                    f"frozen SSL pool S{sector} allowlist hash mismatch"
                )
            allowlist = pd.read_csv(allowlist_path)
            if list(allowlist.columns) != ["tic"]:
                raise ValueError(f"frozen SSL pool S{sector} allowlist schema")
            allowlist_tics = _positive_integer_series(
                allowlist["tic"], name="tic"
            )
            expected_tics = (
                pool.loc[pool["sector"].eq(sector), "tic"]
                .astype(np.int64)
                .sort_values()
                .reset_index(drop=True)
            )
            if not allowlist_tics.reset_index(drop=True).equals(expected_tics):
                raise ValueError(
                    f"frozen SSL pool S{sector} allowlist does not match pool"
                )
            if int(metadata.get("n_tics", -1)) != len(expected_tics):
                raise ValueError(
                    f"frozen SSL pool S{sector} allowlist count mismatch"
                )
    return pool.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    ), summary


def load_global_full_pool_bls(
    *,
    summary_path: Path,
    frozen_pool: pd.DataFrame,
    frozen_pool_summary_path: Path,
    output_path_override: Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any], Path]:
    """Load only the Parquet artifact authorized by the global BLS summary.

    ``output_path_override`` permits an explicitly staged copy, but never
    changes artifact identity: the copy must match the summary's byte count
    and SHA-256 exactly.
    """

    summary_path = Path(summary_path).expanduser().resolve()
    frozen_pool_summary_path = (
        Path(frozen_pool_summary_path).expanduser().resolve()
    )
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid global full-pool BLS summary: {summary_path}") from exc
    required = {
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
    if not isinstance(summary, dict) or set(summary) != required:
        raise ValueError("global full-pool BLS summary has the wrong fields")
    if summary.get("passed") is not True:
        raise ValueError("global full-pool BLS summary did not pass")
    if summary.get("schema_version") != GLOBAL_BLS_SUMMARY_SCHEMA_VERSION:
        raise ValueError("global full-pool BLS summary has the wrong schema")
    if summary.get("contract_version") != GLOBAL_BLS_CONTRACT_VERSION:
        raise ValueError("global full-pool BLS summary has the wrong contract")
    if summary.get("sectors") != list(FULLPOOL_SSL_SECTORS):
        raise ValueError("global full-pool BLS summary has the wrong sectors")
    if summary.get("observation_identity_columns") != ["sector", "tic"]:
        raise ValueError("global full-pool BLS summary has the wrong identity")

    pool_identity = _pool_identity_sha256(frozen_pool)
    pool_binding = summary.get("frozen_pool")
    if not isinstance(pool_binding, dict):
        raise ValueError("global BLS summary lacks frozen-pool binding")
    declared_pool_summary = pool_binding.get("summary")
    if not isinstance(declared_pool_summary, dict):
        raise ValueError("global BLS summary lacks pool-summary metadata")
    if (
        declared_pool_summary.get("sha256")
        != file_sha256(frozen_pool_summary_path)
        or int(declared_pool_summary.get("size_bytes", -1))
        != frozen_pool_summary_path.stat().st_size
    ):
        raise ValueError("global BLS summary binds a different frozen pool summary")
    if pool_binding.get("contract_version") != FULL_POOL_CONTRACT_VERSION:
        raise ValueError("global BLS summary binds the wrong pool contract")
    if pool_binding.get("observation_identity_sha256") != pool_identity:
        raise ValueError("global BLS summary binds a different pool identity")
    if int(pool_binding.get("n_observations", -1)) != len(frozen_pool):
        raise ValueError("global BLS pool observation count mismatch")
    if int(pool_binding.get("n_unique_tics", -1)) != frozen_pool["tic"].nunique():
        raise ValueError("global BLS pool TIC count mismatch")

    coverage = summary.get("coverage_audit")
    if coverage != {
        "missing_frozen_observations": 0,
        "unexpected_bls_observations": 0,
        "observation_identity_sha256": pool_identity,
    }:
        raise ValueError("global full-pool BLS coverage audit failed")
    output = summary.get("output")
    if not isinstance(output, dict):
        raise ValueError("global full-pool BLS summary lacks output metadata")
    if set(output) != {
        "path",
        "size_bytes",
        "sha256",
        "n_rows",
        "n_observations",
    }:
        raise ValueError("global full-pool BLS output metadata has wrong fields")
    declared_output_path = Path(
        str(output.get("path", ""))
    ).expanduser().resolve()
    output_path = (
        Path(output_path_override).expanduser().resolve()
        if output_path_override is not None
        else declared_output_path
    )
    if output_path.suffix.lower() not in {".parquet", ".pq"}:
        raise ValueError("global full-pool BLS output must be Parquet")
    if not output_path.is_file():
        raise FileNotFoundError(output_path)
    if (
        output.get("sha256") != file_sha256(output_path)
        or int(output.get("size_bytes", -1)) != output_path.stat().st_size
    ):
        raise ValueError("global full-pool BLS output hash/size mismatch")
    bls = read_table(output_path)
    counts = summary.get("counts")
    if not isinstance(counts, dict):
        raise ValueError("global full-pool BLS summary lacks counts")
    declared_rows = int(output.get("n_rows", -1))
    declared_observations = int(output.get("n_observations", -1))
    if declared_rows != int(counts.get("n_rows", -1)) or declared_rows != len(bls):
        raise ValueError("global full-pool BLS output row count mismatch")
    if (
        declared_observations != int(counts.get("n_observations", -1))
        or declared_observations != len(frozen_pool)
    ):
        raise ValueError("global full-pool BLS output observation count mismatch")
    normalized_sector = _positive_integer_series(
        bls["sector"], name="sector"
    )
    normalized_tic = _positive_integer_series(bls["tic"], name="tic")
    identities = (
        pd.DataFrame({"sector": normalized_sector, "tic": normalized_tic})
        .drop_duplicates()
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    expected = (
        frozen_pool.loc[:, ["sector", "tic"]]
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    if not identities.equals(expected):
        raise ValueError("global full-pool BLS identities differ from frozen pool")
    return bls, summary, output_path


def _normalize_reserved_tics(
    reserved_hosts: pd.DataFrame | Sequence[int],
) -> tuple[set[int], str]:
    if isinstance(reserved_hosts, pd.DataFrame):
        if "tic" not in reserved_hosts:
            raise KeyError("reserved-host inventory lacks tic")
        values = _positive_integer_series(reserved_hosts["tic"], name="tic")
    else:
        values = _positive_integer_series(
            pd.Series(list(reserved_hosts), dtype=object),
            name="tic",
        )
    unique = sorted({int(value) for value in values})
    if not unique:
        raise ValueError(
            "reserved-host inventory is empty; an accepted prospective "
            "inventory is required before the broad SSL registry can freeze"
        )
    return set(unique), _canonical_sha256(unique)


def read_tic_inventory(path: Path) -> pd.DataFrame:
    """Read a strict sorted one-TIC-per-line file or a table with ``tic``."""

    path = Path(path).expanduser().resolve()
    if path.suffix.lower() != ".txt":
        table = read_table(path)
        if "tic" not in table:
            raise KeyError(f"TIC inventory lacks tic: {path}")
        values = _positive_integer_series(table["tic"], name="tic")
        if values.duplicated().any():
            raise ValueError("TIC inventory contains duplicate values")
        return pd.DataFrame({"tic": values.astype(np.int64)})

    payload = path.read_text(encoding="ascii")
    lines = payload.splitlines()
    if not lines or any(not line or line != line.strip() for line in lines):
        raise ValueError(
            "text TIC inventory must contain one nonblank integer per line"
        )
    if any(not line.isdecimal() for line in lines):
        raise ValueError("text TIC inventory contains a non-decimal value")
    tics = [int(line) for line in lines]
    if any(tic <= 0 for tic in tics):
        raise ValueError("text TIC inventory contains a non-positive TIC")
    if tics != sorted(set(tics)):
        raise ValueError("text TIC inventory must be sorted and unique")
    return pd.DataFrame({"tic": np.asarray(tics, dtype=np.int64)})


def _rank_one_anchor_rows(
    bls_rows: pd.DataFrame,
    *,
    aperture: str,
) -> pd.DataFrame:
    forbidden = sorted(_FORBIDDEN_LABEL_COLUMNS & set(bls_rows.columns))
    if forbidden:
        raise ValueError(
            "full-pool BLS observations must be label-free; "
            f"found={forbidden}"
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
        raise KeyError(f"full-pool BLS table lacks columns: {missing}")
    if bls_rows.empty:
        raise ValueError("full-pool BLS table is empty")

    work = bls_rows.copy()
    work["sector"] = _positive_integer_series(work["sector"], name="sector")
    work["tic"] = _positive_integer_series(work["tic"], name="tic")
    rank = _integer_series(work["peak_rank"], name="peak_rank")
    anchor = work.loc[
        work["aperture"].fillna("").astype(str).eq(str(aperture))
        & rank.isin((0, 1))
    ].copy()
    if anchor.empty:
        raise ValueError(
            f"BLS table contains no rank-one/status rows for {aperture}"
        )
    duplicate = anchor.duplicated(["sector", "tic"], keep=False)
    if duplicate.any():
        examples = (
            anchor.loc[duplicate, ["sector", "tic", "peak_rank", "status"]]
            .head(10)
            .to_dict(orient="records")
        )
        raise ValueError(
            "rank-one BLS anchor is not unique per (sector, tic); "
            f"first={examples}"
        )
    anchor["peak_rank"] = _integer_series(
        anchor["peak_rank"], name="peak_rank"
    )
    anchor["status"] = (
        anchor["status"].fillna("").astype(str).str.strip().str.lower()
    )
    anchor["ssl_observation_id"] = [
        f"s{int(sector):04d}-tic{int(tic)}"
        for sector, tic in zip(anchor["sector"], anchor["tic"])
    ]
    if anchor["ssl_observation_id"].duplicated().any():
        raise RuntimeError("generated SSL observation identities are not unique")
    return anchor.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )


def build_fullpool_ssl_registry(
    bls_rows: pd.DataFrame,
    native_registry: pd.DataFrame,
    split_registry: pd.DataFrame,
    reserved_hosts: pd.DataFrame | Sequence[int],
    *,
    frozen_pool: pd.DataFrame,
    eligibility: EligibilityAuthority | None = None,
    sectors: Sequence[int] = FULLPOOL_SSL_SECTORS,
    anchor_aperture: str = FULLPOOL_SSL_ANCHOR_APERTURE,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build the label-free, host-safe broad-pool SSL registry.

    ``frozen_pool`` is the preregistered LC allowlist.  ``bls_rows`` must cover
    it exactly: the ADP-small rank-one row is retained for searchable targets,
    while a rank-zero status row carries an unsearchable observation through
    the exclusion audit.  ``native_registry`` is a separate observation-keyed
    storage registry and must not be the Teacher-v3 split-bound training table.
    """

    expected_sectors = tuple(sorted({int(value) for value in sectors}))
    if not expected_sectors or any(value <= 0 for value in expected_sectors):
        raise ValueError("sectors must contain positive integers")
    missing_pool = sorted(set(POOL_COLUMNS) - set(frozen_pool.columns))
    if missing_pool:
        raise KeyError(f"frozen SSL pool lacks columns: {missing_pool}")
    pool = frozen_pool.loc[:, list(POOL_COLUMNS)].copy()
    pool["sector"] = _positive_integer_series(pool["sector"], name="sector")
    pool["tic"] = _positive_integer_series(pool["tic"], name="tic")
    if pool.empty or pool.duplicated(["sector", "tic"]).any():
        raise ValueError(
            "frozen SSL pool must have unique nonempty observation keys"
        )
    if not pool["pool_contract_version"].astype(str).eq(
        FULL_POOL_CONTRACT_VERSION
    ).all():
        raise ValueError("frozen SSL pool has the wrong contract version")
    pool_sectors = tuple(sorted(pool["sector"].unique().astype(int)))
    if pool_sectors != expected_sectors:
        raise ValueError(
            f"frozen SSL pool sectors {pool_sectors} != {expected_sectors}"
        )

    decisions = derive_anchor_eligibility(
        bls_rows,
        frozen_pool=pool,
        anchor_aperture=str(anchor_aperture),
    )
    anchor = decisions.rename(
        columns={"observation_id": "ssl_observation_id"}
    )
    observed_sectors = tuple(sorted(anchor["sector"].unique().astype(int)))
    if observed_sectors != expected_sectors:
        raise ValueError(
            "full-pool BLS sectors do not match the requested release: "
            f"{observed_sectors} != {expected_sectors}"
        )
    pool_keys = set(zip(pool["sector"].astype(int), pool["tic"].astype(int)))
    anchor_keys = set(
        zip(anchor["sector"].astype(int), anchor["tic"].astype(int))
    )
    if anchor_keys != pool_keys:
        missing = sorted(pool_keys - anchor_keys)[:10]
        extra = sorted(anchor_keys - pool_keys)[:10]
        raise ValueError(
            "BLS anchor coverage differs from the frozen SSL pool; "
            f"missing={missing}, extra={extra}"
        )
    anchor = anchor.sort_values(
        ["sector", "tic"], kind="stable"
    ).reset_index(drop=True)
    if eligibility is not None and (
        set(eligibility.full_keys) != pool_keys
        or set(
            zip(
                anchor.loc[
                    anchor["native_model_eligible"], "sector"
                ].astype(int),
                anchor.loc[
                    anchor["native_model_eligible"], "tic"
                ].astype(int),
            )
        )
        != set(eligibility.eligible_keys)
    ):
        raise ValueError(
            "BLS-derived eligibility differs from the frozen eligibility authority"
        )

    native_required = {
        "sector",
        "tic",
        "native_h5_path",
        "native_group_path",
        "native_h5_sha256",
        "native_contract_version",
    }
    missing_native = sorted(native_required - set(native_registry.columns))
    if missing_native:
        raise KeyError(
            f"full-pool native registry lacks columns: {missing_native}"
        )
    forbidden_native = sorted(
        _FORBIDDEN_LABEL_COLUMNS & set(native_registry.columns)
    )
    if forbidden_native:
        raise ValueError(
            "full-pool native registry must be label-free; "
            f"found={forbidden_native}"
        )
    native = native_registry.loc[:, list(native_required)].copy()
    native["sector"] = _positive_integer_series(native["sector"], name="sector")
    native["tic"] = _positive_integer_series(native["tic"], name="tic")
    if native.duplicated(["sector", "tic"]).any():
        raise ValueError(
            "full-pool native registry contains duplicate (sector, tic) keys"
        )
    native["native_h5_path"] = _nonblank_text(
        native["native_h5_path"],
        name="native_h5_path",
    ).map(lambda value: str(Path(value).expanduser().resolve()))
    native["native_group_path"] = _nonblank_text(
        native["native_group_path"],
        name="native_group_path",
    )
    native["native_h5_sha256"] = _nonblank_text(
        native["native_h5_sha256"],
        name="native_h5_sha256",
    ).str.lower()
    if not native["native_h5_sha256"].str.fullmatch(_SHA256_PATTERN).all():
        raise ValueError("full-pool native registry has invalid SHA-256 values")
    native["native_contract_version"] = _nonblank_text(
        native["native_contract_version"],
        name="native_contract_version",
    )
    native_keys = set(
        zip(native["sector"].astype(int), native["tic"].astype(int))
    )

    assignments = validate_tic_split_assignments(
        split_registry,
        require_unique_tics=True,
        require_complete_partitions=False,
    )
    fixed_test_tics = set(
        assignments.loc[
            assignments["fixed_split"].eq("test"), "tic"
        ].astype(int)
    )
    development_folds = {
        int(row.tic): int(row.cv_fold)
        for row in assignments.loc[
            assignments["fixed_split"].eq("development")
        ].itertuples(index=False)
    }
    reserved_tics, reserved_tics_sha256 = _normalize_reserved_tics(
        reserved_hosts
    )

    merged = anchor.merge(
        native,
        on=["sector", "tic"],
        how="left",
        validate="one_to_one",
        indicator="_native_join",
    )
    period = pd.to_numeric(merged["period_d"], errors="coerce")
    epoch = pd.to_numeric(merged["t0_bjd"], errors="coerce")
    duration = pd.to_numeric(merged["duration_min"], errors="coerce")
    rank = _integer_series(merged["peak_rank"], name="peak_rank")
    recomputed_bls_eligible = (
        merged["status"].eq("ok")
        & rank.eq(1)
        & np.isfinite(period)
        & period.gt(0)
        & np.isfinite(epoch)
        & np.isfinite(duration)
        & duration.gt(0)
    )
    bls_eligible = merged["native_model_eligible"].astype(bool)
    if not bls_eligible.equals(recomputed_bls_eligible.astype(bool)):
        raise RuntimeError("shared BLS eligibility decision drifted during merge")
    bls_eligible_keys = set(
        zip(
            merged.loc[bls_eligible, "sector"].astype(int),
            merged.loc[bls_eligible, "tic"].astype(int),
        )
    )
    if native_keys != bls_eligible_keys:
        missing = sorted(bls_eligible_keys - native_keys)[:10]
        unexpected = sorted(native_keys - bls_eligible_keys)[:10]
        if missing:
            raise ValueError(
                "leakage-safe searchable BLS rows lack full-pool native "
                f"mappings; first={missing}"
            )
        raise ValueError(
            "full-pool native registry contains observations outside the "
            f"BLS-eligible partition; first={unexpected}"
        )
    if eligibility is not None:
        if native_keys != set(eligibility.eligible_keys):
            raise ValueError(
                "full-pool native registry differs from eligibility authority"
            )
        if not native["native_contract_version"].eq(
            FULL_POOL_NATIVE_CONTRACT_VERSION
        ).all():
            raise ValueError(
                "production full-pool registry requires native-v2 inputs"
            )
    native_present = merged["_native_join"].eq("both")
    if "is_injected_row" in merged:
        injected = _strict_bool_series(
            merged["is_injected_row"],
            name="is_injected_row",
        )
    else:
        injected = pd.Series(False, index=merged.index, dtype=bool)
    fixed_test = merged["tic"].isin(fixed_test_tics)
    reserved = merged["tic"].isin(reserved_tics)
    if fixed_test.any() or reserved.any():
        raise ValueError(
            "frozen SSL pool retained fixed-test or prospective-reserved hosts"
        )
    missing_for_searchable = (
        bls_eligible
        & ~injected
        & ~fixed_test
        & ~reserved
        & ~native_present
    )
    if missing_for_searchable.any():
        examples = (
            merged.loc[missing_for_searchable, ["sector", "tic"]]
            .head(10)
            .to_dict(orient="records")
        )
        raise ValueError(
            "leakage-safe searchable BLS rows lack full-pool native mappings; "
            f"first={examples}"
        )
    include = bls_eligible & native_present & ~injected & ~fixed_test & ~reserved

    reasons: list[str] = []
    for index in merged.index:
        current: list[str] = []
        if bool(fixed_test.loc[index]):
            current.append("fixed_test_tic")
        if bool(reserved.loc[index]):
            current.append("reserved_prospective_tic")
        if bool(injected.loc[index]):
            current.append("injected_row")
        if not bool(bls_eligible.loc[index]):
            current.append("bls_unsearchable")
        if not bool(native_present.loc[index]):
            current.append("native_missing")
        reasons.append("+".join(current))

    held_out_folds: list[int] = []
    fold_sources: list[str] = []
    for tic, active in zip(merged["tic"].astype(int), include):
        if not bool(active):
            held_out_folds.append(-1)
            fold_sources.append("excluded")
        elif int(tic) in development_folds:
            held_out_folds.append(int(development_folds[int(tic)]))
            fold_sources.append("frozen_development_split")
        else:
            held_out_folds.append(-1)
            fold_sources.append("unlabeled_all_folds")

    registry = pd.DataFrame(
        {
            "registry_schema_version": FULLPOOL_SSL_REGISTRY_SCHEMA,
            "ssl_observation_id": merged["ssl_observation_id"].astype(str),
            "sector": merged["sector"].astype(np.int16),
            "tic": merged["tic"].astype(np.int64),
            "period_d": period.astype(float),
            "t0_bjd": epoch.astype(float),
            "duration_min": duration.astype(float),
            "bls_status": merged["status"].astype(str),
            "native_h5_path": merged["native_h5_path"].fillna("").astype(str),
            "native_group_path": merged["native_group_path"].fillna("").astype(str),
            "native_h5_sha256": merged["native_h5_sha256"].fillna("").astype(str),
            "native_contract_version": (
                merged["native_contract_version"].fillna("").astype(str)
            ),
            "is_injected_row": injected.astype(bool),
            "fixed_test_member": fixed_test.astype(bool),
            "reserved_prospective_member": reserved.astype(bool),
            "ssl_pool_include": include.astype(bool),
            "ssl_pool_exclusion_reason": reasons,
            "ssl_held_out_fold": np.asarray(held_out_folds, dtype=np.int8),
            "fold_assignment_source": fold_sources,
        }
    )
    registry = validate_fullpool_ssl_registry(
        registry,
        expected_sectors=expected_sectors,
    )

    native_records = _native_file_records(registry)
    audit = {
        "summary_schema_version": FULLPOOL_SSL_REGISTRY_SUMMARY_SCHEMA,
        "registry_schema_version": FULLPOOL_SSL_REGISTRY_SCHEMA,
        "sectors": list(expected_sectors),
        "anchor_aperture": str(anchor_aperture),
        "rank_policy": (
            "shared_native_model_eligibility_v2:"
            "peak_rank_1_or_rank_0_status_row"
        ),
        "frozen_pool_contract_version": FULL_POOL_CONTRACT_VERSION,
        "frozen_pool_identity_sha256": _pool_identity_sha256(pool),
        "real_only": True,
        "labels_consumed": False,
        "injections_consumed": False,
        "fixed_test_exclusion_scope": "whole_tic",
        "prospective_exclusion_scope": "whole_tic",
        "fold_assignment_policy": (
            "frozen_development_fold_else_present_in_all_pretraining_folds_v2"
        ),
        "n_bls_rows": int(len(bls_rows)),
        "n_anchor_observations": int(len(registry)),
        "n_anchor_tics": int(registry["tic"].nunique()),
        "n_ssl_pool_observations": int(registry["ssl_pool_include"].sum()),
        "n_ssl_pool_tics": int(
            registry.loc[registry["ssl_pool_include"], "tic"].nunique()
        ),
        "n_fixed_test_observations_excluded": int(fixed_test.sum()),
        "n_fixed_test_tics_in_pool": int(
            registry.loc[fixed_test, "tic"].nunique()
        ),
        "n_reserved_observations_excluded": int(reserved.sum()),
        "n_reserved_tics_in_pool": int(
            registry.loc[reserved, "tic"].nunique()
        ),
        "n_injected_observations_excluded": int(injected.sum()),
        "n_bls_unsearchable_observations": int((~bls_eligible).sum()),
        "native_model_eligibility_contract_version": (
            eligibility.contract_version if eligibility is not None else ""
        ),
        "native_model_full_identity_sha256": (
            eligibility.full_observation_identity_sha256
            if eligibility is not None
            else _pool_identity_sha256(pool)
        ),
        "native_model_eligible_identity_sha256": (
            eligibility.eligible_observation_identity_sha256
            if eligibility is not None
            else _pool_identity_sha256(
                registry.loc[registry["ssl_pool_include"]]
            )
        ),
        "native_model_excluded_identity_sha256": (
            eligibility.excluded_observation_identity_sha256
            if eligibility is not None
            else _pool_identity_sha256(
                registry.loc[~registry["ssl_pool_include"]]
            )
        ),
        "fixed_split_registry_tics": int(assignments["tic"].nunique()),
        "fixed_test_registry_tics": int(len(fixed_test_tics)),
        "reserved_inventory_tics": int(len(reserved_tics)),
        "reserved_inventory_tics_sha256": reserved_tics_sha256,
        "canonical_bls_anchor_sha256": _frame_sha256(
            anchor,
            (
                "ssl_observation_id",
                "sector",
                "tic",
                "period_d",
                "t0_bjd",
                "duration_min",
                "status",
            ),
        ),
        "canonical_native_registry_sha256": _frame_sha256(
            native.sort_values(["sector", "tic"], kind="stable"),
            (
                "sector",
                "tic",
                "native_h5_path",
                "native_group_path",
                "native_h5_sha256",
                "native_contract_version",
            ),
        ),
        "canonical_split_assignments_sha256": _frame_sha256(
            assignments.sort_values("tic", kind="stable"),
            ("tic", "fixed_split", "cv_fold"),
        ),
        "canonical_registry_sha256": _frame_sha256(
            registry,
            FULLPOOL_SSL_REGISTRY_COLUMNS,
        ),
        "held_out_fold_counts": _held_out_fold_counts(registry),
        "native_files": native_records,
    }
    return registry, audit


def validate_fullpool_ssl_registry(
    registry: pd.DataFrame,
    *,
    expected_sectors: Sequence[int] = FULLPOOL_SSL_SECTORS,
) -> pd.DataFrame:
    """Normalize and validate a broad-pool registry without opening HDF5."""

    missing = sorted(set(FULLPOOL_SSL_REGISTRY_COLUMNS) - set(registry.columns))
    if missing:
        raise KeyError(f"full-pool SSL registry lacks columns: {missing}")
    forbidden = sorted(_FORBIDDEN_LABEL_COLUMNS & set(registry.columns))
    if forbidden:
        raise ValueError(
            "full-pool SSL registry must remain label-free; "
            f"found={forbidden}"
        )
    if registry.empty:
        raise ValueError("full-pool SSL registry is empty")
    work = registry.loc[:, list(FULLPOOL_SSL_REGISTRY_COLUMNS)].copy()
    work["registry_schema_version"] = _nonblank_text(
        work["registry_schema_version"],
        name="registry_schema_version",
    )
    if not work["registry_schema_version"].eq(
        FULLPOOL_SSL_REGISTRY_SCHEMA
    ).all():
        raise ValueError("full-pool SSL registry has the wrong schema version")
    work["ssl_observation_id"] = _nonblank_text(
        work["ssl_observation_id"],
        name="ssl_observation_id",
    )
    if work["ssl_observation_id"].duplicated().any():
        raise ValueError("full-pool SSL observation IDs are not unique")
    work["sector"] = _positive_integer_series(work["sector"], name="sector")
    work["tic"] = _positive_integer_series(work["tic"], name="tic")
    if work.duplicated(["sector", "tic"]).any():
        raise ValueError(
            "full-pool SSL registry contains duplicate (sector, tic) keys"
        )
    expected = tuple(sorted({int(value) for value in expected_sectors}))
    observed = tuple(sorted(work["sector"].unique().astype(int)))
    if observed != expected:
        raise ValueError(
            f"full-pool SSL registry sectors {observed} != expected {expected}"
        )
    for column in (
        "is_injected_row",
        "fixed_test_member",
        "reserved_prospective_member",
        "ssl_pool_include",
    ):
        work[column] = _strict_bool_series(work[column], name=column)
    work["ssl_held_out_fold"] = _integer_series(
        work["ssl_held_out_fold"],
        name="ssl_held_out_fold",
    )
    work["bls_status"] = (
        work["bls_status"].fillna("").astype(str).str.strip().str.lower()
    )
    work["ssl_pool_exclusion_reason"] = (
        work["ssl_pool_exclusion_reason"].fillna("").astype(str).str.strip()
    )
    work["fold_assignment_source"] = _nonblank_text(
        work["fold_assignment_source"],
        name="fold_assignment_source",
    )

    include = work["ssl_pool_include"]
    excluded = ~include
    if not include.any():
        raise ValueError("full-pool SSL registry has no included observations")
    forbidden_include = include & (
        work["is_injected_row"]
        | work["fixed_test_member"]
        | work["reserved_prospective_member"]
    )
    if forbidden_include.any():
        raise ValueError(
            "full-pool SSL registry includes injected, fixed-test, or "
            "prospective-reserved observations"
        )
    if not work.loc[include, "ssl_held_out_fold"].isin(
        (-1, *range(FULLPOOL_SSL_N_FOLDS))
    ).all():
        raise ValueError(
            "included full-pool rows require ssl_held_out_fold in {-1,0,...,4}"
        )
    if not work.loc[excluded, "ssl_held_out_fold"].eq(-1).all():
        raise ValueError(
            "excluded full-pool rows require ssl_held_out_fold=-1"
        )
    if not work.loc[excluded, "fold_assignment_source"].eq("excluded").all():
        raise ValueError("excluded full-pool rows require assignment source excluded")
    if work.loc[excluded, "ssl_pool_exclusion_reason"].eq("").any():
        raise ValueError("excluded full-pool rows require an exclusion reason")
    if work.loc[include, "ssl_pool_exclusion_reason"].ne("").any():
        raise ValueError("included full-pool rows may not carry an exclusion reason")
    if not work.loc[include, "fold_assignment_source"].isin(
        {"frozen_development_split", "unlabeled_all_folds"}
    ).all():
        raise ValueError("included full-pool rows have invalid fold assignment sources")
    frozen = include & work["fold_assignment_source"].eq(
        "frozen_development_split"
    )
    unlabeled = include & work["fold_assignment_source"].eq(
        "unlabeled_all_folds"
    )
    if not work.loc[frozen, "ssl_held_out_fold"].isin(
        range(FULLPOOL_SSL_N_FOLDS)
    ).all():
        raise ValueError(
            "frozen development TICs require a held fold in [0,4]"
        )
    if not work.loc[unlabeled, "ssl_held_out_fold"].eq(-1).all():
        raise ValueError(
            "previously unseen TICs must be present in all pretraining folds"
        )

    for column in ("period_d", "t0_bjd", "duration_min"):
        work[column] = pd.to_numeric(work[column], errors="coerce").astype(float)
    valid_ephemeris = (
        np.isfinite(work["period_d"])
        & work["period_d"].gt(0)
        & np.isfinite(work["t0_bjd"])
        & np.isfinite(work["duration_min"])
        & work["duration_min"].gt(0)
    )
    if not valid_ephemeris.loc[include].all():
        raise ValueError("included full-pool rows have invalid BLS ephemerides")
    if not work.loc[include, "bls_status"].eq("ok").all():
        raise ValueError("included full-pool rows require BLS status ok")

    for column in (
        "native_h5_path",
        "native_group_path",
        "native_h5_sha256",
        "native_contract_version",
    ):
        work[column] = work[column].fillna("").astype(str).str.strip()
        if work.loc[include, column].eq("").any():
            raise ValueError(f"included full-pool rows have blank {column}")
    if not work.loc[include, "native_h5_sha256"].str.fullmatch(
        _SHA256_PATTERN
    ).all():
        raise ValueError("included full-pool rows have invalid native hashes")
    if not work.loc[include, "native_contract_version"].eq(
        FULL_POOL_NATIVE_CONTRACT_VERSION
    ).all():
        raise ValueError("included full-pool rows require native-v2 inputs")
    if not work.loc[include, "native_h5_path"].map(
        lambda value: Path(value).is_absolute()
    ).all():
        raise ValueError("included full-pool native paths must be absolute")
    storage_duplicate = work.loc[include].duplicated(
        ["native_h5_path", "native_group_path"],
        keep=False,
    )
    if storage_duplicate.any():
        raise ValueError(
            "full-pool native storage maps to more than one observation"
        )
    for path, group in work.loc[include].groupby("native_h5_path", sort=False):
        if group["native_h5_sha256"].nunique() != 1:
            raise ValueError(f"native file {path} has conflicting hashes")
        if group["native_contract_version"].nunique() != 1:
            raise ValueError(f"native file {path} has conflicting contracts")

    host = work.groupby("tic", sort=False).agg(
        fixed_test_memberships=("fixed_test_member", "nunique"),
        prospective_memberships=("reserved_prospective_member", "nunique"),
        included_folds=(
            "ssl_held_out_fold",
            lambda values: values[values >= 0].nunique(),
        ),
    )
    if host["fixed_test_memberships"].gt(1).any():
        raise ValueError("fixed-test membership is not host-wide")
    if host["prospective_memberships"].gt(1).any():
        raise ValueError("prospective membership is not host-wide")
    if host["included_folds"].gt(1).any():
        raise ValueError("one TIC maps to more than one SSL fold")
    return work.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )


def _held_out_fold_counts(registry: pd.DataFrame) -> list[dict[str, int]]:
    include = registry["ssl_pool_include"]
    records: list[dict[str, int]] = []
    for fold in range(FULLPOOL_SSL_N_FOLDS):
        mask = include & registry["ssl_held_out_fold"].eq(fold)
        records.append(
            {
                "fold": int(fold),
                "n_observations": int(mask.sum()),
                "n_tics": int(registry.loc[mask, "tic"].nunique()),
            }
        )
    return records


def _native_file_records(registry: pd.DataFrame) -> list[dict[str, Any]]:
    include = registry["ssl_pool_include"]
    records: list[dict[str, Any]] = []
    for path, group in registry.loc[include].groupby(
        "native_h5_path", sort=True
    ):
        records.append(
            {
                "native_h5_path": str(path),
                "native_h5_sha256": str(group["native_h5_sha256"].iloc[0]),
                "native_contract_version": str(
                    group["native_contract_version"].iloc[0]
                ),
                "n_observations": int(len(group)),
            }
        )
    return records


def _table_bytes(table: pd.DataFrame, path: Path) -> bytes:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return table.to_csv(
            index=False,
            lineterminator="\n",
            float_format="%.17g",
        ).encode("utf-8")
    if suffix in {".parquet", ".pq"}:
        stream = io.BytesIO()
        table.to_parquet(stream, index=False)
        return stream.getvalue()
    raise ValueError(f"unsupported registry table format: {path}")


def _publish_immutable(path: Path, payload: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace immutable full-pool artifact: {path}"
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


def write_fullpool_ssl_registry(
    registry: pd.DataFrame,
    audit: Mapping[str, Any],
    *,
    registry_path: Path,
    summary_path: Path,
    source_provenance: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Publish the broad-pool registry and its checksum-bound summary once."""

    normalized = validate_fullpool_ssl_registry(registry)
    registry_path = Path(registry_path).expanduser().resolve()
    summary_path = Path(summary_path).expanduser().resolve()
    registry_payload = _table_bytes(normalized, registry_path)
    _publish_immutable(registry_path, registry_payload)
    registry_sha256 = file_sha256(registry_path)
    expected_counts = {
        "n_anchor_observations": int(len(normalized)),
        "n_anchor_tics": int(normalized["tic"].nunique()),
        "n_ssl_pool_observations": int(
            normalized["ssl_pool_include"].sum()
        ),
        "n_ssl_pool_tics": int(
            normalized.loc[normalized["ssl_pool_include"], "tic"].nunique()
        ),
    }
    for name, value in expected_counts.items():
        if int(audit.get(name, -1)) != value:
            raise ValueError(f"registry audit {name} does not match table")
    summary = {
        **dict(audit),
        "summary_schema_version": FULLPOOL_SSL_REGISTRY_SUMMARY_SCHEMA,
        "registry_schema_version": FULLPOOL_SSL_REGISTRY_SCHEMA,
        "registry_path": str(registry_path),
        "registry_sha256": registry_sha256,
        "held_out_fold_counts": _held_out_fold_counts(normalized),
        "native_files": _native_file_records(normalized),
        "source_provenance": dict(source_provenance or {}),
    }
    summary_payload = (
        json.dumps(
            _canonical_value(summary),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    _publish_immutable(summary_path, summary_payload)
    return {
        "registry": str(registry_path),
        "registry_sha256": registry_sha256,
        "summary": str(summary_path),
        "summary_sha256": file_sha256(summary_path),
        **expected_counts,
    }


def load_fullpool_ssl_registry(
    *,
    registry_path: Path,
    summary_path: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load and cross-check one immutable broad-pool registry pair."""

    registry_path = Path(registry_path).expanduser().resolve()
    summary_path = Path(summary_path).expanduser().resolve()
    registry = validate_fullpool_ssl_registry(read_table(registry_path))
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid full-pool summary {summary_path}: {exc}") from exc
    if not isinstance(summary, dict):
        raise ValueError("full-pool registry summary must be a JSON object")
    if summary.get("summary_schema_version") != (
        FULLPOOL_SSL_REGISTRY_SUMMARY_SCHEMA
    ):
        raise ValueError("full-pool registry summary has the wrong schema")
    if summary.get("registry_schema_version") != FULLPOOL_SSL_REGISTRY_SCHEMA:
        raise ValueError("full-pool registry summary has the wrong registry schema")
    if summary.get("registry_sha256") != file_sha256(registry_path):
        raise ValueError("full-pool registry hash does not match its summary")
    counts = {
        "n_anchor_observations": int(len(registry)),
        "n_anchor_tics": int(registry["tic"].nunique()),
        "n_ssl_pool_observations": int(registry["ssl_pool_include"].sum()),
        "n_ssl_pool_tics": int(
            registry.loc[registry["ssl_pool_include"], "tic"].nunique()
        ),
    }
    for name, value in counts.items():
        if int(summary.get(name, -1)) != value:
            raise ValueError(f"full-pool summary {name} does not match registry")
    if summary.get("held_out_fold_counts") != _held_out_fold_counts(registry):
        raise ValueError(
            "full-pool summary held-out fold counts do not match registry"
        )
    if summary.get("native_files") != _native_file_records(registry):
        raise ValueError("full-pool summary native-file records do not match")
    return registry, summary


def select_fullpool_ssl_fold(
    registry: pd.DataFrame,
    *,
    held_out_fold: int,
    max_rows: int | None = None,
    required_observation_ids: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Select the four-fold, real-only corpus for one encoder job."""

    fold = int(held_out_fold)
    if fold not in range(FULLPOOL_SSL_N_FOLDS):
        raise ValueError("held_out_fold must be in [0,4]")
    if isinstance(required_observation_ids, (str, bytes)):
        raise TypeError(
            "required_observation_ids must be a sequence of complete IDs"
        )
    required_values = (
        () if required_observation_ids is None else required_observation_ids
    )
    required: list[str] = []
    for index, value in enumerate(required_values):
        if (
            not isinstance(value, str)
            or not value
            or value != value.strip()
        ):
            raise ValueError(
                "required_observation_ids contains an invalid value at "
                f"index {index}"
            )
        required.append(value)
    if len(required) != len(set(required)):
        raise ValueError("required_observation_ids contains duplicates")
    required = sorted(required)
    if required and max_rows is None:
        raise ValueError(
            "required_observation_ids are only valid with bounded max_rows"
        )
    work = validate_fullpool_ssl_registry(registry)
    include = work["ssl_pool_include"]
    held = include & work["ssl_held_out_fold"].eq(fold)
    selected = work.loc[include & ~held].copy()
    registry_ids = set(work["ssl_observation_id"].astype(str))
    missing_required = sorted(set(required) - registry_ids)
    if missing_required:
        raise ValueError(
            "required observations are absent from the full-pool registry: "
            f"{missing_required}"
        )
    excluded_required = sorted(
        set(required)
        & set(work.loc[~include, "ssl_observation_id"].astype(str))
    )
    if excluded_required:
        raise ValueError(
            "required observations are not model eligible: "
            f"{excluded_required}"
        )
    held_required = sorted(
        set(required)
        & set(work.loc[held, "ssl_observation_id"].astype(str))
    )
    if held_required:
        raise ValueError(
            "required observations belong to the held-out fold: "
            f"{held_required}"
        )
    if max_rows is not None:
        if int(max_rows) < 2:
            raise ValueError("max_rows must be at least 2")
        if len(required) > int(max_rows):
            raise ValueError(
                "required observations exceed the bounded max_rows"
            )
        if len(selected) > int(max_rows):
            order = selected["ssl_observation_id"].map(
                lambda value: hashlib.sha256(
                    f"smoke:{fold}:{value}".encode("utf-8")
                ).hexdigest()
            )
            required_mask = selected["ssl_observation_id"].isin(required)
            required_indices = selected.index[required_mask]
            remaining_indices = order.loc[~required_mask].sort_values(
                kind="stable"
            ).index[: int(max_rows) - len(required_indices)]
            selected = selected.loc[
                list(required_indices) + list(remaining_indices)
            ]
    selected = selected.sort_values(
        ["sector", "tic"], kind="stable"
    ).reset_index(drop=True)
    if len(selected) < 2:
        raise ValueError("fold-local full-pool SSL selection has fewer than two rows")
    selected_tics = set(selected["tic"].astype(int))
    held_tics = set(work.loc[held, "tic"].astype(int))
    fixed_test_tics = set(
        work.loc[work["fixed_test_member"], "tic"].astype(int)
    )
    prospective_tics = set(
        work.loc[work["reserved_prospective_member"], "tic"].astype(int)
    )
    selected_observation_ids = set(selected["ssl_observation_id"].astype(str))
    required_selected = set(required) <= selected_observation_ids
    if not required_selected:
        raise RuntimeError(
            "fold-local selection omitted a required observation"
        )
    disjoint = {
        "held_fold_tics": not bool(selected_tics & held_tics),
        "fixed_test_tics": not bool(selected_tics & fixed_test_tics),
        "reserved_prospective_tics": not bool(
            selected_tics & prospective_tics
        ),
    }
    if not all(disjoint.values()):
        raise RuntimeError(f"full-pool fold selection leaks hosts: {disjoint}")
    selection_columns = (
        "ssl_observation_id",
        "sector",
        "tic",
        "period_d",
        "t0_bjd",
        "duration_min",
        "native_h5_path",
        "native_group_path",
        "native_h5_sha256",
        "native_contract_version",
        "ssl_held_out_fold",
    )
    audit = {
        "selection_schema_version": FULLPOOL_SSL_SELECTION_SCHEMA,
        "held_out_fold": fold,
        "n_registry_observations": int(len(work)),
        "n_eligible_observations": int(include.sum()),
        "n_eligible_tics": int(work.loc[include, "tic"].nunique()),
        "n_held_observations": int(held.sum()),
        "n_held_tics": int(len(held_tics)),
        "n_selected_observations": int(len(selected)),
        "n_selected_tics": int(len(selected_tics)),
        "max_rows": None if max_rows is None else int(max_rows),
        "required_observation_ids": required,
        "n_required_observations": int(len(required)),
        "required_observations_selected": required_selected,
        "selected_rows_sha256": _frame_sha256(selected, selection_columns),
        "selected_tics_sha256": _canonical_sha256(sorted(selected_tics)),
        "tic_disjoint": disjoint,
    }
    return selected, audit


def fullpool_dataset_rows(selected: pd.DataFrame) -> pd.DataFrame:
    """Adapt label-free registry rows to ``HarmonicNativeDataset`` inputs."""

    if selected.empty:
        raise ValueError("cannot build a dataset adapter from empty rows")
    forbidden = sorted(_FORBIDDEN_LABEL_COLUMNS & set(selected.columns))
    if forbidden:
        raise ValueError(f"selected full-pool rows unexpectedly contain {forbidden}")
    rows = selected.copy()
    rows["review_id"] = rows["ssl_observation_id"].astype(str)
    rows["fixed_split"] = "development"
    rows["cv_fold"] = rows["ssl_held_out_fold"].astype(np.int8)
    rows["is_injected_row"] = False
    rows["input_variant"] = "observed"
    for column in (
        "morphology_target_index",
        "preserve_target_index",
        "harmonic_target_index",
    ):
        rows[column] = -1
    for column in (
        "morphology_weight",
        "preserve_weight",
        "harmonic_weight",
        "compact_weight",
    ):
        rows[column] = np.float32(0.0)
    rows["compact_target_index"] = -1
    rows["pretrain_target"] = -1
    return rows.reset_index(drop=True)


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def _json_bytes(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            _json_safe(payload),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _atomic_write(path: Path, payload: bytes) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        digest = file_sha256(temporary)
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    if file_sha256(path) != digest:
        raise RuntimeError(f"atomic output changed while installing: {path}")
    return digest


def _atomic_torch_save(payload: Mapping[str, Any], path: Path) -> str:
    import torch

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        torch.save(dict(payload), temporary)
        digest = file_sha256(temporary)
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    if file_sha256(path) != digest:
        raise RuntimeError(f"checkpoint changed while installing: {path}")
    return digest


def _code_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(__file__).resolve().parents[3],
        check=False,
        capture_output=True,
        text=True,
    )
    revision = completed.stdout.strip()
    if completed.returncode != 0 or len(revision) != 40:
        raise RuntimeError("cannot bind full-pool SSL to a Git revision")
    return revision


def _artifact_metadata(path: Path) -> dict[str, Any]:
    """Return a stable path/size/hash binding for one required authority."""

    resolved = Path(path).expanduser().resolve(strict=True)
    before = resolved.stat()
    if not resolved.is_file() or before.st_size <= 0:
        raise FileNotFoundError(f"training authority is missing or empty: {resolved}")
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
        raise RuntimeError(f"training authority changed while hashing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": int(after.st_size),
        "sha256": digest,
    }


def _metadata_from_record(
    record: Any,
    *,
    context: str,
) -> dict[str, Any]:
    if not isinstance(record, Mapping):
        raise ValueError(f"{context} must be an artifact binding")
    path_text = record.get("path")
    digest = str(record.get("sha256", "")).lower()
    try:
        size_bytes = int(record.get("size_bytes", -1))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} has an invalid size") from exc
    if not isinstance(path_text, str) or not path_text.strip():
        raise ValueError(f"{context} has an invalid path")
    if len(digest) != 64 or any(
        character not in "0123456789abcdef" for character in digest
    ):
        raise ValueError(f"{context} has an invalid SHA-256")
    if size_bytes <= 0:
        raise ValueError(f"{context} has an invalid size")
    return {
        "path": str(Path(path_text).expanduser().resolve()),
        "size_bytes": size_bytes,
        "sha256": digest,
    }


def _observation_keys(rows: pd.DataFrame) -> frozenset[tuple[int, int]]:
    return frozenset(
        zip(
            pd.to_numeric(rows["sector"], errors="raise").astype(int),
            pd.to_numeric(rows["tic"], errors="raise").astype(int),
        )
    )


def _normalized_native_mapping(
    rows: pd.DataFrame,
) -> list[dict[str, Any]]:
    required = {
        "sector",
        "tic",
        "native_h5_path",
        "native_group_path",
        "native_h5_sha256",
        "native_contract_version",
    }
    missing = sorted(required - set(rows))
    if missing:
        raise KeyError(f"native authority mapping lacks columns: {missing}")
    records: list[dict[str, Any]] = []
    for row in rows.sort_values(["sector", "tic"], kind="stable").itertuples(
        index=False
    ):
        path = Path(str(row.native_h5_path)).expanduser()
        if not path.is_absolute():
            raise ValueError("production native mappings must use absolute paths")
        records.append(
            {
                "sector": int(row.sector),
                "tic": int(row.tic),
                "native_h5_path": str(path.resolve()),
                "native_group_path": str(row.native_group_path).strip(),
                "native_h5_sha256": str(row.native_h5_sha256).strip().lower(),
                "native_contract_version": str(
                    row.native_contract_version
                ).strip(),
            }
        )
    return records


def _validate_training_authority_chain(
    *,
    eligibility: EligibilityAuthority,
    native_registry: pd.DataFrame,
    native_release: Mapping[str, Any],
    registry: pd.DataFrame,
    registry_summary: Mapping[str, Any],
    numeric_gate_release: Mapping[str, Any],
    authority_bindings: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Fail closed unless every production training authority agrees exactly."""

    if (
        not isinstance(numeric_gate_release, Mapping)
        or numeric_gate_release.get("passed") is not True
    ):
        raise ValueError("model-input numerical gate did not pass")
    if (
        numeric_gate_release.get("model_input_contract_version")
        != MODEL_INPUT_CONTRACT_VERSION
    ):
        raise ValueError(
            "model-input numerical gate binds the wrong input contract"
        )
    if (
        numeric_gate_release.get("envelope_contract")
        != MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
    ):
        raise ValueError(
            "model-input numerical gate binds the wrong envelope contract"
        )
    expected_partition = {
        "full": {
            "count": PRODUCTION_FULL_OBSERVATIONS,
            "identity_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
        },
        "eligible": {
            "count": PRODUCTION_ELIGIBLE_OBSERVATIONS,
            "identity_sha256": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
        },
        "excluded": {
            "count": PRODUCTION_EXCLUDED_OBSERVATIONS,
            "identity_sha256": PRODUCTION_EXCLUDED_IDENTITY_SHA256,
        },
    }
    observed_partition = {
        "full": {
            "count": int(eligibility.n_full),
            "identity_sha256": (
                eligibility.full_observation_identity_sha256
            ),
        },
        "eligible": {
            "count": int(eligibility.n_eligible),
            "identity_sha256": (
                eligibility.eligible_observation_identity_sha256
            ),
        },
        "excluded": {
            "count": int(eligibility.n_excluded),
            "identity_sha256": (
                eligibility.excluded_observation_identity_sha256
            ),
        },
    }
    if observed_partition != expected_partition:
        raise ValueError(
            "training eligibility differs from the locked production partition"
        )
    if (
        eligibility.eligible_keys & eligibility.excluded_keys
        or eligibility.eligible_keys | eligibility.excluded_keys
        != eligibility.full_keys
    ):
        raise ValueError("training eligibility is not an exact partition")
    for authority_name, eligibility_name in (
        ("eligibility_exclusions", "exclusions"),
        ("eligibility_summary", "eligibility_summary"),
    ):
        observed = _metadata_from_record(
            authority_bindings.get(authority_name),
            context=f"authority_bindings.{authority_name}",
        )
        expected = eligibility.bindings[eligibility_name].metadata()
        if observed != expected:
            raise ValueError(
                f"explicit {authority_name} differs from loaded eligibility"
            )

    full_keys = _observation_keys(registry)
    include = registry["ssl_pool_include"].astype(bool)
    included_keys = _observation_keys(registry.loc[include])
    excluded_keys = _observation_keys(registry.loc[~include])
    if (
        full_keys != eligibility.full_keys
        or included_keys != eligibility.eligible_keys
        or excluded_keys != eligibility.excluded_keys
    ):
        raise ValueError(
            "final SSL registry does not equal the eligibility partition"
        )
    if (
        registry["is_injected_row"].astype(bool).any()
        or registry["fixed_test_member"].astype(bool).any()
        or registry["reserved_prospective_member"].astype(bool).any()
    ):
        raise ValueError(
            "production SSL registry contains forbidden injected or held hosts"
        )
    native_columns = (
        "native_h5_path",
        "native_group_path",
        "native_h5_sha256",
        "native_contract_version",
    )
    excluded_native_present = (
        registry.loc[~include, list(native_columns)]
        .fillna("")
        .astype(str)
        .apply(lambda column: column.str.strip())
        .ne("")
        .to_numpy(dtype=bool)
        .any()
    )
    if excluded_native_present:
        raise ValueError("excluded SSL rows unexpectedly have native mappings")
    if not registry.loc[
        ~include, "ssl_pool_exclusion_reason"
    ].astype(str).eq("bls_unsearchable+native_missing").all():
        raise ValueError(
            "final SSL exclusions do not match the locked BLS/native complement"
        )

    final_mapping = _normalized_native_mapping(registry.loc[include])
    native_mapping = _normalized_native_mapping(native_registry)
    if final_mapping != native_mapping:
        raise ValueError(
            "final SSL registry native mapping differs from native-v2 release"
        )
    mapping_frame = pd.DataFrame.from_records(native_mapping)
    mapping_sha256 = _frame_sha256(
        mapping_frame,
        (
            "sector",
            "tic",
            "native_h5_path",
            "native_group_path",
            "native_h5_sha256",
            "native_contract_version",
        ),
    )

    expected_summary_values = {
        "n_anchor_observations": PRODUCTION_FULL_OBSERVATIONS,
        "n_ssl_pool_observations": PRODUCTION_ELIGIBLE_OBSERVATIONS,
        "n_bls_unsearchable_observations": PRODUCTION_EXCLUDED_OBSERVATIONS,
        "native_model_full_identity_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
        "native_model_eligible_identity_sha256": (
            PRODUCTION_ELIGIBLE_IDENTITY_SHA256
        ),
        "native_model_excluded_identity_sha256": (
            PRODUCTION_EXCLUDED_IDENTITY_SHA256
        ),
        "canonical_native_registry_sha256": mapping_sha256,
        "frozen_pool_identity_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
    }
    for name, expected in expected_summary_values.items():
        if registry_summary.get(name) != expected:
            raise ValueError(
                f"final SSL registry summary {name} differs from production"
            )
    if (
        registry_summary.get("real_only") is not True
        or registry_summary.get("labels_consumed") is not False
        or registry_summary.get("injections_consumed") is not False
        or int(registry_summary.get("n_injected_observations_excluded", -1))
        != 0
    ):
        raise ValueError("final SSL registry data-usage audit failed")
    if registry_summary.get("canonical_registry_sha256") != _frame_sha256(
        registry,
        FULLPOOL_SSL_REGISTRY_COLUMNS,
    ):
        raise ValueError("final SSL registry canonical identity differs")

    release_counts = native_release.get("counts")
    release_identities = native_release.get("identity_hashes")
    if not isinstance(release_counts, Mapping) or not isinstance(
        release_identities, Mapping
    ):
        raise ValueError("native-v2 release lacks partition metadata")
    if {
        "full": int(release_counts.get("full_observations", -1)),
        "eligible": int(release_counts.get("eligible_observations", -1)),
        "excluded": int(release_counts.get("excluded_observations", -1)),
        "native": int(
            release_counts.get("native_registry_observations", -1)
        ),
    } != {
        "full": PRODUCTION_FULL_OBSERVATIONS,
        "eligible": PRODUCTION_ELIGIBLE_OBSERVATIONS,
        "excluded": PRODUCTION_EXCLUDED_OBSERVATIONS,
        "native": PRODUCTION_ELIGIBLE_OBSERVATIONS,
    }:
        raise ValueError("native-v2 release counts differ from production")
    if {
        "full": release_identities.get("full_observations_sha256"),
        "eligible": release_identities.get("eligible_observations_sha256"),
        "excluded": release_identities.get("excluded_observations_sha256"),
        "native": release_identities.get(
            "native_registry_observations_sha256"
        ),
    } != {
        "full": PRODUCTION_FULL_IDENTITY_SHA256,
        "eligible": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
        "excluded": PRODUCTION_EXCLUDED_IDENTITY_SHA256,
        "native": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    }:
        raise ValueError("native-v2 release identities differ from production")

    provenance = registry_summary.get("source_provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("final SSL registry summary lacks source provenance")
    required_provenance = {
        "frozen_pool": eligibility.bindings["frozen_pool"].metadata(),
        "frozen_pool_summary": (
            eligibility.bindings["frozen_pool_summary"].metadata()
        ),
        "bls_summary": eligibility.bindings["global_bls_summary"].metadata(),
        "eligibility_exclusions": authority_bindings[
            "eligibility_exclusions"
        ],
        "eligibility_summary": authority_bindings["eligibility_summary"],
        "native_registry": authority_bindings["native_registry"],
        "native_registry_summary": authority_bindings[
            "native_registry_summary"
        ],
        "native_release_summary": authority_bindings[
            "native_release_summary"
        ],
    }
    for name, expected in required_provenance.items():
        observed = _metadata_from_record(
            provenance.get(name),
            context=f"source_provenance.{name}",
        )
        normalized_expected = _metadata_from_record(
            expected,
            context=f"authority_bindings.{name}",
        )
        if observed != normalized_expected:
            raise ValueError(
                f"final SSL registry provenance {name} differs from authority"
            )
    for name in ("frozen_split_registry", "reserved_hosts"):
        _metadata_from_record(
            provenance.get(name),
            context=f"source_provenance.{name}",
        )
    if "bls_peaks_override" in provenance:
        if _metadata_from_record(
            provenance["bls_peaks_override"],
            context="source_provenance.bls_peaks_override",
        ) != eligibility.bindings["global_bls"].metadata():
            raise ValueError(
                "final SSL registry BLS override differs from eligibility"
            )

    frozen_binding = provenance.get("frozen_pool_authority_bindings")
    split_binding = _metadata_from_record(
        provenance["frozen_split_registry"],
        context="source_provenance.frozen_split_registry",
    )
    reserved_binding = _metadata_from_record(
        provenance["reserved_hosts"],
        context="source_provenance.reserved_hosts",
    )
    expected_frozen_binding = {
        "split_registry_sha256_equal": True,
        "reserved_hosts_sha256_equal": True,
        "split_registry_sha256": split_binding["sha256"],
        "reserved_hosts_sha256": reserved_binding["sha256"],
    }
    if frozen_binding != expected_frozen_binding:
        raise ValueError("final SSL frozen-pool provenance binding failed")

    global_bls_binding = provenance.get("global_bls_authority_binding")
    expected_global_bls_binding = {
        "summary_path": str(
            eligibility.bindings["global_bls_summary"].path
        ),
        "summary_sha256": (
            eligibility.bindings["global_bls_summary"].sha256
        ),
        "artifact_path": str(eligibility.bindings["global_bls"].path),
        "artifact_sha256": eligibility.bindings["global_bls"].sha256,
        "artifact_matches_summary": True,
        "summary_schema_version": GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
        "summary_contract_version": GLOBAL_BLS_CONTRACT_VERSION,
    }
    if global_bls_binding != expected_global_bls_binding:
        raise ValueError("final SSL global-BLS provenance binding failed")

    eligibility_binding = provenance.get(
        "native_model_eligibility_binding"
    )
    expected_eligibility_binding = {
        "contract_version": eligibility.contract_version,
        "release_binding": eligibility.release_binding,
        "full_observations": PRODUCTION_FULL_OBSERVATIONS,
        "eligible_observations": PRODUCTION_ELIGIBLE_OBSERVATIONS,
        "excluded_observations": PRODUCTION_EXCLUDED_OBSERVATIONS,
        "full_observation_identity_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
        "eligible_observation_identity_sha256": (
            PRODUCTION_ELIGIBLE_IDENTITY_SHA256
        ),
        "excluded_observation_identity_sha256": (
            PRODUCTION_EXCLUDED_IDENTITY_SHA256
        ),
    }
    if eligibility_binding != expected_eligibility_binding:
        raise ValueError("final SSL eligibility provenance binding failed")

    native_binding = provenance.get("native_release_binding")
    expected_native_binding = {
        "schema_version": native_release.get("schema_version"),
        "release_binding": native_release.get("release_binding"),
        "native_contract_version": native_release.get(
            "native_contract_version"
        ),
        "release_summary_sha256": authority_bindings[
            "native_release_summary"
        ]["sha256"],
    }
    if native_binding != expected_native_binding:
        raise ValueError("final SSL native-release provenance binding failed")

    return {
        "schema_version": FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA,
        "production_lock_passed": True,
        "partition": expected_partition,
        "native_mapping_sha256": mapping_sha256,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "model_input_numeric_envelope_contract": (
            MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "numeric_gate_release": {
            "binding": dict(authority_bindings["numeric_gate_release"]),
            "schema_version": numeric_gate_release.get("schema_version"),
            "envelope_canonical_sha256": numeric_gate_release.get(
                "envelope_canonical_sha256"
            ),
            "counts": dict(numeric_gate_release.get("counts", {})),
            "identity_hashes": dict(
                numeric_gate_release.get("identity_hashes", {})
            ),
            "code_revision": numeric_gate_release.get("code_revision"),
            "passed": True,
        },
        "source_provenance_verified": True,
        "authority_bindings": {
            name: dict(metadata)
            for name, metadata in sorted(authority_bindings.items())
        },
    }


def _load_training_authority_chain(
    *,
    eligibility_exclusions_path: Path,
    eligibility_summary_path: Path,
    native_registry_path: Path,
    native_registry_summary_path: Path,
    native_release_summary_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    numeric_gate_release_path: Path,
    expected_code_revision: str,
) -> tuple[
    pd.DataFrame,
    dict[str, Any],
    pd.DataFrame,
    dict[str, Any],
    dict[str, Any],
]:
    paths = {
        "eligibility_exclusions": eligibility_exclusions_path,
        "eligibility_summary": eligibility_summary_path,
        "native_registry": native_registry_path,
        "native_registry_summary": native_registry_summary_path,
        "native_release_summary": native_release_summary_path,
        "registry": registry_path,
        "registry_summary": registry_summary_path,
        "numeric_gate_release": numeric_gate_release_path,
    }
    authority_bindings = {
        name: _artifact_metadata(path) for name, path in paths.items()
    }
    eligibility = load_native_model_eligibility(
        authority_bindings["eligibility_exclusions"]["path"],
        authority_bindings["eligibility_summary"]["path"],
        production_lock=True,
        rederive_from_bls=False,
    )
    native_registry, native_release = (
        load_full_pool_native_registry_release(
            registry_path=authority_bindings["native_registry"]["path"],
            registry_summary_path=authority_bindings[
                "native_registry_summary"
            ]["path"],
            release_summary_path=authority_bindings[
                "native_release_summary"
            ]["path"],
            eligibility=eligibility,
            verify_shard_files=False,
        )
    )
    registry, registry_summary = load_fullpool_ssl_registry(
        registry_path=authority_bindings["registry"]["path"],
        summary_path=authority_bindings["registry_summary"]["path"],
    )
    numeric_gate_release = validate_numeric_gate_release(
        Path(authority_bindings["numeric_gate_release"]["path"]),
        expected_code_revision=expected_code_revision,
        expected_authority_bindings={
            "ssl_registry": authority_bindings["registry"],
            "ssl_registry_summary": authority_bindings["registry_summary"],
            "native_registry": authority_bindings["native_registry"],
            "native_registry_summary": authority_bindings[
                "native_registry_summary"
            ],
            "native_release_summary": authority_bindings[
                "native_release_summary"
            ],
        },
    )
    if _artifact_metadata(
        authority_bindings["numeric_gate_release"]["path"]
    ) != authority_bindings["numeric_gate_release"]:
        raise RuntimeError(
            "numeric-gate release changed during authority validation"
        )
    audit = _validate_training_authority_chain(
        eligibility=eligibility,
        native_registry=native_registry,
        native_release=native_release,
        registry=registry,
        registry_summary=registry_summary,
        numeric_gate_release=numeric_gate_release,
        authority_bindings=authority_bindings,
    )
    return (
        registry,
        registry_summary,
        native_registry,
        native_release,
        audit,
    )


def _verify_native_files(
    selected: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Hash selected files and inspect their actual v2 root/group identities."""

    import h5py

    records: list[dict[str, Any]] = []
    for path_text, group in selected.groupby("native_h5_path", sort=True):
        path = Path(str(path_text))
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"missing full-pool native file: {path}")
        declared = str(group["native_h5_sha256"].iloc[0])
        before = path.stat()
        observed = file_sha256(path)
        if observed != declared:
            raise ValueError(f"native file hash changed: {path}")
        with h5py.File(path, "r") as h5:
            actual_contract = str(h5.attrs.get("contract_version", "")).strip()
            if actual_contract != FULL_POOL_NATIVE_CONTRACT_VERSION:
                raise ValueError(
                    f"selected native file is not native-v2: {path}"
                )
            for row in group.itertuples(index=False):
                if str(row.native_contract_version) != actual_contract:
                    raise ValueError(
                        f"selected native contract declaration changed: {path}"
                    )
                group_path = str(row.native_group_path)
                if group_path not in h5 or not isinstance(
                    h5[group_path], h5py.Group
                ):
                    raise KeyError(
                        f"missing selected native group {group_path!r} in {path}"
                    )
                native_group = h5[group_path]
                actual_key = (
                    int(native_group.attrs.get("sector", -1)),
                    int(native_group.attrs.get("tic", -1)),
                )
                expected_key = (int(row.sector), int(row.tic))
                if actual_key != expected_key:
                    raise ValueError(
                        "selected native group identity differs from registry: "
                        f"{path}:{group_path} {actual_key} != {expected_key}"
                    )
        after = path.stat()
        if (
            int(before.st_size),
            int(before.st_mtime_ns),
            int(before.st_dev),
            int(before.st_ino),
        ) != (
            int(after.st_size),
            int(after.st_mtime_ns),
            int(after.st_dev),
            int(after.st_ino),
        ):
            raise RuntimeError(
                f"selected native file changed during verification: {path}"
            )
        records.append(
            {
                "native_h5_path": str(path),
                "native_h5_sha256": declared,
                "native_h5_size_bytes": int(path.stat().st_size),
                "native_contract_version": actual_contract,
                "hash_verified_now": True,
                "root_contract_verified_now": True,
                "group_identities_verified_now": True,
                "n_selected_observations": int(len(group)),
            }
        )
    return records


def _loss_and_components(result: Any) -> tuple[Any, dict[str, float]]:
    if isinstance(result, tuple) and len(result) == 2:
        loss, raw = result
    elif isinstance(result, Mapping) and "loss" in result:
        loss = result["loss"]
        raw = {key: value for key, value in result.items() if key != "loss"}
    else:
        loss = result
        raw = {}
    components: dict[str, float] = {}
    for name, value in raw.items():
        try:
            components[str(name)] = float(value.detach())
        except AttributeError:
            components[str(name)] = float(value)
    return loss, components


def _assert_finite_loss_and_components(
    loss: Any,
    components: Mapping[str, float],
    *,
    fold: int,
    epoch: int,
    batch_index: int,
) -> None:
    """Fail before backpropagation when a training objective is non-finite."""

    import torch

    detached = loss.detach()
    if detached.numel() != 1 or not bool(torch.isfinite(detached).all().item()):
        raise FloatingPointError(
            "non-finite full-pool SSL loss "
            f"(fold={fold}, epoch={epoch}, batch={batch_index})"
        )
    invalid = {
        str(name): float(value)
        for name, value in components.items()
        if not math.isfinite(float(value))
    }
    if invalid:
        raise FloatingPointError(
            "non-finite full-pool SSL loss components "
            f"(fold={fold}, epoch={epoch}, batch={batch_index}): {invalid}"
        )


def _assert_finite_module_state(
    modules: Mapping[str, Any],
    *,
    fold: int,
    epoch: int,
) -> None:
    """Reject a numerically invalid encoder/projector before publication."""

    import torch

    for module_name, module in modules.items():
        for tensor_name, tensor in module.state_dict().items():
            if not torch.is_tensor(tensor) or not bool(
                torch.isfinite(tensor.detach()).all().item()
            ):
                raise FloatingPointError(
                    "non-finite full-pool SSL model state "
                    f"(fold={fold}, epoch={epoch}, "
                    f"tensor={module_name}.{tensor_name})"
                )


def _projection_head(input_dim: int) -> Any:
    from torch import nn

    return nn.Sequential(
        nn.Linear(int(input_dim), 128),
        nn.GELU(),
        nn.Linear(128, 64),
    )


def _prepare_fold_directory(
    *,
    fold_dir: Path,
    run_contract: Mapping[str, Any],
    resume: bool,
) -> tuple[str, bool]:
    """Install or validate the immutable run contract for one fold."""

    fold_dir = Path(fold_dir)
    contract_payload = _json_bytes(run_contract)
    contract_sha256 = hashlib.sha256(contract_payload).hexdigest()
    contract_path = fold_dir / "run_contract.json"
    if fold_dir.exists():
        if not fold_dir.is_dir() or fold_dir.is_symlink():
            raise RuntimeError(f"invalid full-pool fold directory: {fold_dir}")
        prior = next(fold_dir.iterdir(), None)
        if prior is not None and not resume:
            raise FileExistsError(
                "full-pool fold output is not empty; pass --resume only for "
                f"the identical contract: {fold_dir}"
            )
    else:
        fold_dir.mkdir(parents=True, exist_ok=False)
    if contract_path.exists():
        if contract_path.read_bytes() != contract_payload:
            raise RuntimeError(
                "resume contract differs from the existing full-pool fold run"
            )
        return contract_sha256, True
    if any(fold_dir.iterdir()):
        raise RuntimeError(
            "full-pool fold directory has content but no run contract"
        )
    _atomic_write(contract_path, contract_payload)
    return contract_sha256, False


def _load_resume_checkpoint(
    *,
    fold_dir: Path,
    contract_sha256: str,
    device: Any,
) -> dict[str, Any] | None:
    import torch

    state_path = Path(fold_dir) / "resume_state.json"
    if not state_path.exists():
        return None
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if state.get("schema_version") != FULLPOOL_SSL_RESUME_SCHEMA:
        raise ValueError("full-pool resume state has the wrong schema")
    if state.get("run_contract_sha256") != contract_sha256:
        raise ValueError("full-pool resume state belongs to another contract")
    checkpoint_path = Path(str(state.get("checkpoint", "")))
    if not checkpoint_path.is_absolute():
        checkpoint_path = Path(fold_dir) / checkpoint_path
    if not checkpoint_path.is_file():
        raise FileNotFoundError(checkpoint_path)
    if file_sha256(checkpoint_path) != state.get("checkpoint_sha256"):
        raise ValueError("full-pool resume checkpoint hash mismatch")
    checkpoint = torch.load(
        checkpoint_path,
        map_location=device,
        weights_only=False,
    )
    if checkpoint.get("schema_version") != FULLPOOL_SSL_RESUME_SCHEMA:
        raise ValueError("full-pool resume checkpoint has the wrong schema")
    if checkpoint.get("run_contract_sha256") != contract_sha256:
        raise ValueError("full-pool resume checkpoint belongs to another contract")
    return checkpoint


def _publish_resume_checkpoint(
    *,
    fold_dir: Path,
    contract_sha256: str,
    checkpoint: Mapping[str, Any],
    epoch: int,
    global_step: int,
) -> tuple[Path, str]:
    """Publish a versioned checkpoint, then atomically advance its pointer.

    The checkpoint named by ``resume_state.json`` is never overwritten.  A
    crash after the new checkpoint is installed but before the pointer moves
    therefore leaves the prior generation authoritative and loadable.  The
    unreferenced generation can be safely replaced when that epoch is replayed.
    """

    fold_dir = Path(fold_dir)
    epoch = int(epoch)
    global_step = int(global_step)
    if epoch < 1 or global_step < 1:
        raise ValueError("resume generation requires a positive epoch and step")
    checkpoint_path = (
        fold_dir
        / f"resume_epoch_{epoch:04d}_step_{global_step:012d}.pt"
    )
    state_path = fold_dir / "resume_state.json"
    if state_path.exists():
        active = json.loads(state_path.read_text(encoding="utf-8"))
        active_path = Path(str(active.get("checkpoint", "")))
        if not active_path.is_absolute():
            active_path = fold_dir / active_path
        if active_path.resolve() == checkpoint_path.resolve():
            raise RuntimeError(
                "refusing to overwrite the authoritative resume checkpoint"
            )

    checkpoint_sha256 = _atomic_torch_save(checkpoint, checkpoint_path)
    _atomic_write(
        state_path,
        _json_bytes(
            {
                "schema_version": FULLPOOL_SSL_RESUME_SCHEMA,
                "run_contract_sha256": contract_sha256,
                "checkpoint": checkpoint_path.name,
                "checkpoint_sha256": checkpoint_sha256,
                "epoch": epoch,
                "global_step": global_step,
            }
        ),
    )
    return checkpoint_path, checkpoint_sha256


def _capture_rng_state() -> dict[str, Any]:
    import torch

    return {
        "python_random_state": random.getstate(),
        "numpy_random_state": np.random.get_state(),
        "torch_cpu_rng_state": torch.get_rng_state(),
        "torch_cuda_rng_states": (
            torch.cuda.get_rng_state_all() if torch.cuda.is_available() else []
        ),
    }


def _restore_rng_state(state: Mapping[str, Any]) -> None:
    import torch

    required = {
        "python_random_state",
        "numpy_random_state",
        "torch_cpu_rng_state",
        "torch_cuda_rng_states",
    }
    missing = sorted(required - set(state))
    if missing:
        raise KeyError(f"resume checkpoint lacks RNG states: {missing}")
    random.setstate(state["python_random_state"])
    np.random.set_state(state["numpy_random_state"])
    torch.set_rng_state(state["torch_cpu_rng_state"])
    cuda_states = state["torch_cuda_rng_states"]
    if torch.cuda.is_available():
        if len(cuda_states) != torch.cuda.device_count():
            raise ValueError("resume checkpoint CUDA RNG device count changed")
        torch.cuda.set_rng_state_all(cuda_states)


def run_fullpool_ssl_fold(
    *,
    eligibility_exclusions_path: Path,
    eligibility_summary_path: Path,
    native_registry_path: Path,
    native_registry_summary_path: Path,
    native_release_summary_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    numeric_gate_release_path: Path,
    out_root: Path,
    fold: int,
    epochs: int = 20,
    batch_size: int = 64,
    workers: int = 8,
    seed: int = FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
    learning_rate: float = 3.0e-4,
    weight_decay: float = 1.0e-4,
    checkpoint_every: int = 1,
    resume: bool = False,
    require_cuda: bool = True,
    max_rows: int | None = None,
    required_observation_ids: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Train or resume one broad-pool fold-local VICReg encoder."""

    import torch

    fold = int(fold)
    if fold not in range(FULLPOOL_SSL_N_FOLDS):
        raise ValueError("fold must be in [0,4]")
    for name, value in (
        ("epochs", epochs),
        ("batch_size", batch_size),
        ("checkpoint_every", checkpoint_every),
    ):
        if int(value) < 1:
            raise ValueError(f"{name} must be positive")
    if int(workers) < 0:
        raise ValueError("workers must be nonnegative")
    if not np.isfinite(learning_rate) or float(learning_rate) <= 0:
        raise ValueError("learning_rate must be finite and positive")
    if not np.isfinite(weight_decay) or float(weight_decay) < 0:
        raise ValueError("weight_decay must be finite and nonnegative")

    code_revision = _code_revision()
    (
        registry,
        registry_summary,
        _native_registry,
        _native_release,
        authority_audit,
    ) = _load_training_authority_chain(
        eligibility_exclusions_path=eligibility_exclusions_path,
        eligibility_summary_path=eligibility_summary_path,
        native_registry_path=native_registry_path,
        native_registry_summary_path=native_registry_summary_path,
        native_release_summary_path=native_release_summary_path,
        registry_path=registry_path,
        registry_summary_path=registry_summary_path,
        numeric_gate_release_path=numeric_gate_release_path,
        expected_code_revision=code_revision,
    )
    registry_binding = authority_audit["authority_bindings"]["registry"]
    registry_summary_binding = authority_audit["authority_bindings"][
        "registry_summary"
    ]
    registry_path = Path(registry_binding["path"])
    registry_summary_path = Path(registry_summary_binding["path"])
    selected, selection_audit = select_fullpool_ssl_fold(
        registry,
        held_out_fold=fold,
        max_rows=max_rows,
        required_observation_ids=required_observation_ids,
    )
    dataset_rows = fullpool_dataset_rows(selected)
    native_files = _verify_native_files(selected)
    fold_seed = int(seed) + 1000 * fold
    model_config = HarmonicModelConfig(metadata_dim=0)
    augmentation_config = EventPreservingAugmentationConfig()
    vicreg_config = VICRegConfig()
    run_contract = {
        "schema_version": FULLPOOL_SSL_RUN_CONTRACT_SCHEMA,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
        "model_facing_name": FULLPOOL_SSL_MODEL_FACING_NAME,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "model_version": MODEL_VERSION,
        "ssl_contract_version": HARMONIC_SSL_CONTRACT_VERSION,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "model_input_numeric_envelope_contract": (
            MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "numeric_gate_release": authority_audit["numeric_gate_release"],
        "profile": FULLPOOL_SSL_PROFILE,
        "fold": fold,
        "registry_path": str(registry_path),
        "registry_sha256": registry_binding["sha256"],
        "registry_summary_path": str(registry_summary_path),
        "registry_summary_sha256": registry_summary_binding["sha256"],
        "registry_summary_schema": registry_summary.get(
            "summary_schema_version"
        ),
        "training_authority": authority_audit,
        "selection_audit": selection_audit,
        "native_files": native_files,
        "native_hashes_verified_now": True,
        "native_root_contracts_verified_now": True,
        "native_group_identities_verified_now": True,
        "model_config": asdict(model_config),
        "augmentation_config": asdict(augmentation_config),
        "vicreg_config": asdict(vicreg_config),
        "projection_architecture": [2 * model_config.embedding_dim, 128, 64],
        "epochs": int(epochs),
        "batch_size": int(batch_size),
        "workers": int(workers),
        "seed": fold_seed,
        "learning_rate": float(learning_rate),
        "weight_decay": float(weight_decay),
        "checkpoint_every": int(checkpoint_every),
        "require_cuda": bool(require_cuda),
        "max_rows": None if max_rows is None else int(max_rows),
        "required_observation_ids": selection_audit[
            "required_observation_ids"
        ],
        "labels_loaded": False,
        "fixed_test_tensors_constructed": False,
        "prospective_sector_tensors_constructed": False,
        "embedding_export": False,
        "neighbor_probe": False,
        "code_revision": code_revision,
    }
    fold_dir = (
        Path(out_root).expanduser().resolve()
        / "encoder_pretraining"
        / f"fold_{fold}"
    )
    contract_sha256, existed = _prepare_fold_directory(
        fold_dir=fold_dir,
        run_contract=run_contract,
        resume=bool(resume),
    )
    summary_path = fold_dir / "summary.json"
    if summary_path.exists():
        if not resume:
            raise FileExistsError(summary_path)
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("run_contract_sha256") != contract_sha256:
            raise ValueError("completed fold summary belongs to another contract")
        checkpoint_path = Path(str(summary.get("checkpoint", "")))
        if not checkpoint_path.is_file():
            raise FileNotFoundError(checkpoint_path)
        if file_sha256(checkpoint_path) != summary.get("checkpoint_sha256"):
            raise ValueError("completed fold checkpoint hash mismatch")
        return summary
    if existed and not resume:
        raise RuntimeError("existing full-pool contract requires --resume")

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required for full-pool SSL but unavailable")
    torch.manual_seed(fold_seed)
    np.random.seed(fold_seed)
    random.seed(fold_seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(fold_seed)

    model = build_harmonic_cnn(
        model_config,
        profile=FULLPOOL_SSL_PROFILE,
    ).to(device)
    projector = _projection_head(2 * model_config.embedding_dim).to(device)
    optimizer = torch.optim.AdamW(
        list(model.parameters()) + list(projector.parameters()),
        lr=float(learning_rate),
        weight_decay=float(weight_decay),
    )
    history: list[dict[str, Any]] = []
    global_step = 0
    start_epoch = 0
    if resume:
        checkpoint = _load_resume_checkpoint(
            fold_dir=fold_dir,
            contract_sha256=contract_sha256,
            device=device,
        )
        if checkpoint is not None:
            model.load_state_dict(checkpoint["encoder_state_dict"], strict=True)
            projector.load_state_dict(
                checkpoint["projection_state_dict"], strict=True
            )
            optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
            start_epoch = int(checkpoint["epoch"])
            global_step = int(checkpoint["global_step"])
            history = list(checkpoint["history"])
            _restore_rng_state(checkpoint["rng_state"])
            if start_epoch > int(epochs):
                raise ValueError("resume checkpoint is beyond requested epochs")

    dataset = HarmonicNativeDataset(
        dataset_rows,
        native_h5=None,
        metadata=np.empty((len(dataset_rows), 0), dtype=np.float32),
        cache_size=0,
        profile=FULLPOOL_SSL_PROFILE,
    )
    try:
        for epoch in range(start_epoch + 1, int(epochs) + 1):
            epoch_seed = fold_seed + epoch * 1_000_003
            loader = _loader(
                dataset,
                np.arange(len(dataset), dtype=int),
                batch_size=int(batch_size),
                shuffle=True,
                workers=int(workers),
                seed=epoch_seed,
            )
            model.train()
            projector.train()
            totals: dict[str, float] = {"loss": 0.0}
            seen = 0
            for batch_index, raw_batch in enumerate(loader):
                if len(raw_batch["review_id"]) < 2:
                    continue
                batch = _to_device(raw_batch, device)
                view_seed = epoch_seed + batch_index * 17
                first_view = augment_ssl_batch(
                    batch,
                    duration_min=batch["duration_min"],
                    config=augmentation_config,
                    seed=view_seed,
                    view_index=0,
                )
                second_view = augment_ssl_batch(
                    batch,
                    duration_min=batch["duration_min"],
                    config=augmentation_config,
                    seed=view_seed,
                    view_index=1,
                )
                optimizer.zero_grad(set_to_none=True)
                with torch.autocast(
                    device_type=device.type,
                    dtype=torch.bfloat16,
                    enabled=device.type == "cuda",
                ):
                    first = projector(model(first_view)["embedding"])
                    second = projector(model(second_view)["embedding"])
                loss, components = _loss_and_components(
                    vicreg_loss(
                        first.float(),
                        second.float(),
                        config=vicreg_config,
                    )
                )
                _assert_finite_loss_and_components(
                    loss,
                    components,
                    fold=fold,
                    epoch=epoch,
                    batch_index=batch_index,
                )
                loss.backward()
                optimizer.step()
                count = len(raw_batch["review_id"])
                totals["loss"] += float(loss.detach()) * count
                for name, value in components.items():
                    totals[name] = totals.get(name, 0.0) + float(value) * count
                seen += count
                global_step += 1
            if seen < 2:
                raise RuntimeError(
                    f"full-pool SSL fold {fold} epoch {epoch} was empty"
                )
            record = {
                "epoch": int(epoch),
                "n_observations": int(seen),
                **{
                    name: float(value / seen)
                    for name, value in totals.items()
                },
            }
            invalid_record = {
                name: value
                for name, value in record.items()
                if name not in {"epoch", "n_observations"}
                and not math.isfinite(float(value))
            }
            if invalid_record:
                raise FloatingPointError(
                    "non-finite full-pool SSL epoch metrics "
                    f"(fold={fold}, epoch={epoch}): {invalid_record}"
                )
            _assert_finite_module_state(
                {"encoder": model, "projector": projector},
                fold=fold,
                epoch=epoch,
            )
            history.append(record)
            print(
                f"[teacher_ssl_fullpool fold={fold}] epoch={epoch}/{epochs} "
                f"rows={seen} loss={record['loss']:.6f}",
                flush=True,
            )
            if epoch % int(checkpoint_every) == 0 or epoch == int(epochs):
                resume_payload = {
                    "schema_version": FULLPOOL_SSL_RESUME_SCHEMA,
                    "run_contract_sha256": contract_sha256,
                    "fold": fold,
                    "epoch": int(epoch),
                    "global_step": int(global_step),
                    "history": history,
                    "rng_state": _capture_rng_state(),
                    "encoder_state_dict": {
                        name: value.detach().cpu()
                        for name, value in model.state_dict().items()
                    },
                    "projection_state_dict": {
                        name: value.detach().cpu()
                        for name, value in projector.state_dict().items()
                    },
                    "optimizer_state_dict": optimizer.state_dict(),
                }
                _atomic_write(
                    fold_dir / "history.csv",
                    pd.DataFrame(history).to_csv(
                        index=False,
                        lineterminator="\n",
                    ).encode("utf-8"),
                )
                _publish_resume_checkpoint(
                    fold_dir=fold_dir,
                    contract_sha256=contract_sha256,
                    checkpoint=resume_payload,
                    epoch=epoch,
                    global_step=global_step,
                )
    finally:
        dataset.close()
    _atomic_write(
        fold_dir / "history.csv",
        pd.DataFrame(history).to_csv(
            index=False,
            lineterminator="\n",
        ).encode("utf-8"),
    )
    checkpoint_path = fold_dir / "encoder.pt"
    checkpoint = {
        "schema_version": FULLPOOL_SSL_CHECKPOINT_SCHEMA,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
        "model_facing_name": FULLPOOL_SSL_MODEL_FACING_NAME,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "model_version": MODEL_VERSION,
        "ssl_contract_version": HARMONIC_SSL_CONTRACT_VERSION,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "model_input_numeric_envelope_contract": (
            MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "numeric_gate_release": authority_audit["numeric_gate_release"],
        "profile": FULLPOOL_SSL_PROFILE,
        "fold": fold,
        "run_contract_sha256": contract_sha256,
        "selection_audit": selection_audit,
        "model_config": asdict(model_config),
        "augmentation_config": asdict(augmentation_config),
        "vicreg_config": asdict(vicreg_config),
        "projection_architecture": [2 * model_config.embedding_dim, 128, 64],
        "epochs": int(epochs),
        "history": history,
        "encoder_state_dict": {
            name: value.detach().cpu()
            for name, value in model.state_dict().items()
        },
        "projection_state_dict": {
            name: value.detach().cpu()
            for name, value in projector.state_dict().items()
        },
    }
    checkpoint_sha256 = _atomic_torch_save(checkpoint, checkpoint_path)
    summary = {
        "schema_version": FULLPOOL_SSL_SUMMARY_SCHEMA,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
        "model_facing_name": FULLPOOL_SSL_MODEL_FACING_NAME,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "model_input_numeric_envelope_contract": (
            MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "numeric_gate_release": authority_audit["numeric_gate_release"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "fold": fold,
        "run_contract": str(fold_dir / "run_contract.json"),
        "run_contract_sha256": contract_sha256,
        "checkpoint": str(checkpoint_path),
        "checkpoint_sha256": checkpoint_sha256,
        "history": str(fold_dir / "history.csv"),
        "history_sha256": file_sha256(fold_dir / "history.csv"),
        "completed_epochs": int(epochs),
        "global_step": int(global_step),
        "selection_audit": selection_audit,
        "fixed_test_status": "host_excluded_tensors_not_constructed",
        "prospective_status": "host_excluded_tensors_not_constructed",
        "labels_loaded": False,
        "automatic_production_promotion": False,
    }
    _atomic_write(summary_path, _json_bytes(summary))
    return summary


__all__ = [
    "FULLPOOL_SSL_ANCHOR_APERTURE",
    "FULLPOOL_SSL_DEFAULT_TRAINING_SEED",
    "FULLPOOL_SSL_ENCODER_NAME",
    "FULLPOOL_SSL_MODEL_FACING_NAME",
    "FULLPOOL_SSL_REGISTRY_COLUMNS",
    "FULLPOOL_SSL_REGISTRY_SCHEMA",
    "FULLPOOL_SSL_REGISTRY_SUMMARY_SCHEMA",
    "FULLPOOL_SSL_RUN_ID",
    "FULLPOOL_SSL_SECTORS",
    "build_fullpool_ssl_registry",
    "fullpool_dataset_rows",
    "load_fullpool_ssl_registry",
    "load_frozen_ssl_full_pool",
    "load_global_full_pool_bls",
    "read_tic_inventory",
    "run_fullpool_ssl_fold",
    "select_fullpool_ssl_fold",
    "validate_fullpool_ssl_registry",
    "write_fullpool_ssl_registry",
]
