#!/usr/bin/env python3
"""Freeze the fail-closed S63 Teacher-v3 preprocessing handoff.

This manifest is written only after the accepted Stage-1 product, cadence
authority, compact ADP pair, prospective inventories, locked BLS search,
rank-one candidate table, raw-source export, and native 4096-period input all
agree by content hash.  It is the final preprocessing gate before any S63
Teacher-v3 score is computed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.search.a2v1_bls_contract import (  # noqa: E402
    A2V1_TEACHER_BLS_SEARCH_CONTRACT,
    approved_a2v1_teacher_bls_config,
    bls_config_sha256,
)
from twirl.vetting.adp_only import (  # noqa: E402
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
)
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    RAW_PAIR_CONTRACT_VERSION,
)
from twirl.vetting.teacher_v3_prospective import (  # noqa: E402
    load_selection_policy,
    load_prospective_plan,
)
from twirl.vetting.s63_preprocessing import (  # noqa: E402
    require_producer_git_sha,
    validate_producer_git_sha,
    validate_s63_stage1_compact,
    validate_s63_stage1_receipt_attestation,
)


CONTRACT_VERSION = "twirl_teacher_v3_s63_prospective_launch_v1"
AUTHORIZATION_CONTRACT = "s63_teacher_v3_prospective_v1"
SECTOR = 63
ORBITS = (133, 134)
PERIODOGRAM_N = 4096
SHA256_RE = re.compile(r"[0-9a-f]{64}")
GIT_SHA_RE = re.compile(r"[0-9a-f]{40}")
RAW_SOURCE_CONTRACT_VERSION = "s56_tglc_raw_pair_v1"
S63_TARGET_METADATA_CONTRACT = (
    "twirl_teacher_v3_s63_rank1_candidate_metadata_v1"
)
S63_TARGET_METADATA_COLUMNS = (
    "adp_sml_own_even_depth",
    "adp_sml_own_odd_depth",
    "adp_sml_own_even_odd_depth_delta",
    "adp_sml_own_even_odd_sigma_delta",
    "adp_sml_trend_ptp",
    "adp_own_even_depth",
    "adp_own_odd_depth",
    "adp_own_even_odd_depth_delta",
    "adp_own_even_odd_sigma_delta",
    "adp_trend_ptp",
)


def file_sha256(path: Path, *, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def _require_regular_file(path: Path, *, label: str) -> Path:
    resolved = Path(path).expanduser().resolve(strict=True)
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise ValueError(f"{label} must be a nonempty regular file: {resolved}")
    lowered_parts = {part.lower() for part in resolved.parts}
    if any(
        token in part
        for part in lowered_parts
        for token in ("pending", ".tmp", "merge_pending")
    ):
        raise ValueError(f"{label} path is pending/temporary: {resolved}")
    return resolved


def _load_json(path: Path, *, label: str) -> dict[str, Any]:
    path = _require_regular_file(path, label=label)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid {label} JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain one JSON object")
    return payload


def _valid_sha256(value: Any, *, label: str) -> str:
    normalized = str(value).strip().lower()
    if SHA256_RE.fullmatch(normalized) is None:
        raise ValueError(f"{label} must be a resolved SHA-256 digest")
    return normalized


def _artifact(path: Path, *, label: str) -> dict[str, Any]:
    resolved = _require_regular_file(path, label=label)
    digest = _valid_sha256(file_sha256(resolved), label=f"{label} SHA-256")
    return {
        "path": str(resolved),
        "sha256": digest,
        "size_bytes": int(resolved.stat().st_size),
    }


def _expect_hash(
    payload: Mapping[str, Any],
    names: str | Sequence[str],
    expected: str,
    *,
    label: str,
) -> str:
    candidates = (names,) if isinstance(names, str) else tuple(names)
    for name in candidates:
        if name in payload:
            observed = _valid_sha256(payload[name], label=f"{label}.{name}")
            if observed != expected:
                raise ValueError(
                    f"{label}.{name} does not bind the expected artifact"
                )
            return name
    raise ValueError(f"{label} lacks required hash field; expected one of {candidates}")


def _expect_hash_anywhere(
    payload: Mapping[str, Any],
    names: str | Sequence[str],
    expected: str,
    *,
    label: str,
) -> str:
    """Accept a declared hash at top level or in ``artifact_hashes`` only."""

    candidates = (names,) if isinstance(names, str) else tuple(names)
    try:
        return _expect_hash(payload, candidates, expected, label=label)
    except ValueError as top_level_error:
        nested = payload.get("artifact_hashes")
        if not isinstance(nested, Mapping):
            raise top_level_error
        return _expect_hash(nested, candidates, expected, label=f"{label}.artifact_hashes")


def validate_accepted_stage1(
    path: Path,
    *,
    expected_producer_git_sha: str | None = None,
) -> dict[str, Any]:
    """Require the accepted full-schema S63 A2v1 gate, including edge policy."""

    resolved = _require_regular_file(path, label="accepted Stage-1 validation")
    summary = _load_json(resolved, label="accepted Stage-1 validation")
    failures: list[str] = []
    if summary.get("sector") != SECTOR:
        failures.append("sector is not 63")
    if sorted(summary.get("orbits", [])) != list(ORBITS):
        failures.append("orbits are not exactly 133/134")
    for field in ("ok", "ok_h5", "ok_fits"):
        if summary.get(field) is not True:
            failures.append(f"{field} is not true")
    for field in ("n_expected_h5", "n_expected_unique_tics"):
        try:
            positive = int(summary.get(field, -1)) > 0
        except (TypeError, ValueError):
            positive = False
        if not positive:
            failures.append(f"{field} is not positive")
    schema = summary.get("schema")
    if not isinstance(schema, Mapping):
        failures.append("schema declaration is missing")
    else:
        expected_schema = {
            "expected_method": "A2v1",
            "expected_prodtag": "A2v1",
            "schema_only": True,
            "check_h5_open": True,
            "allow_edge_warn_missing": True,
        }
        for field, expected in expected_schema.items():
            if schema.get(field) != expected:
                failures.append(f"schema.{field} != {expected!r}")
        required = set(schema.get("required_a2v1_columns", []))
        for column in ("DET_FLUX_ADP_SML", "DET_FLUX_ADP"):
            if column not in required:
                failures.append(f"schema lacks {column}")
    expected_contract = summary.get("expected_contract")
    if not isinstance(expected_contract, Mapping):
        failures.append("expected_contract is missing")
    else:
        if expected_contract.get("ok") is not True:
            failures.append("expected_contract.ok is not true")
        if sorted(expected_contract.get("requested_orbits", [])) != list(ORBITS):
            failures.append("expected_contract requested orbits are not 133/134")
        if expected_contract.get("missing_requested_orbits") != []:
            failures.append("expected_contract has missing requested orbits")
    h5 = summary.get("h5")
    if not isinstance(h5, Mapping):
        failures.append("HDF5 audit is missing")
    else:
        for field in (
            "n_missing_h5_non_edge",
            "n_zero_byte_h5",
            "n_unreadable_h5",
        ):
            if h5.get(field) != 0:
                failures.append(f"h5.{field} is not zero")
    fits = summary.get("fits")
    if not isinstance(fits, Mapping) or fits.get("skipped") is True:
        failures.append("FITS audit is missing/skipped")
    else:
        for field in ("n_missing_fits_non_edge_tics", "n_bad_checked_fits"):
            if fits.get(field) != 0:
                failures.append(f"fits.{field} is not zero")
    if failures:
        raise ValueError(
            "S63 Stage-1 product is not accepted: " + "; ".join(failures)
        )
    receipt_attestation: dict[str, Any] | None = None
    if expected_producer_git_sha is not None:
        receipt_attestation = validate_s63_stage1_receipt_attestation(
            resolved,
            expected_producer_git_sha=expected_producer_git_sha,
        )
    result = {
        "sector": SECTOR,
        "orbits": list(ORBITS),
        "validation_sha256": file_sha256(resolved),
        "n_expected_h5": int(summary.get("n_expected_h5", -1)),
        "n_expected_unique_tics": int(
            summary.get("n_expected_unique_tics", -1)
        ),
        "n_missing_h5_edge_warn": int(h5.get("n_missing_h5_edge_warn", 0)),
        "n_missing_fits_edge_warn_tics": int(
            fits.get("n_missing_fits_edge_warn_tics", 0)
        ),
        "passed": True,
    }
    if receipt_attestation is not None:
        result["immutable_source"] = receipt_attestation
    return result


def _read_tic_inventory(
    path: Path,
    *,
    label: str,
    allow_duplicate_rows: bool = False,
) -> list[int]:
    path = _require_regular_file(path, label=label)
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "tic" not in frame:
            raise ValueError(f"{label} CSV lacks tic")
        values = frame["tic"].astype(str).str.strip().tolist()
    elif path.suffix.lower() in {".txt", ".list"}:
        values = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if values and values[0].lower() == "tic":
            values = values[1:]
    else:
        raise ValueError(f"{label} must be CSV or one-TIC-per-line text")
    if not values or any(not value.isascii() or not value.isdigit() for value in values):
        raise ValueError(f"{label} contains invalid TIC values")
    tics = [int(value) for value in values]
    if any(value <= 0 for value in tics):
        raise ValueError(f"{label} TICs must be positive integers")
    if not allow_duplicate_rows and len(tics) != len(set(tics)):
        raise ValueError(f"{label} TICs must be unique positive integers")
    return sorted(set(tics))


def _read_table(path: Path, *, label: str) -> pd.DataFrame:
    path = _require_regular_file(path, label=label)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"{label} must be CSV or Parquet")


def _validate_cadence(
    table: Path,
    manifest_path: Path,
    *,
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    table = _require_regular_file(table, label="cadence table")
    manifest = _load_json(manifest_path, label="cadence manifest")
    table_digest = file_sha256(table)
    expected = {
        "contract_version": "s63_a2v1_cadence_reference_v1",
        "sector": SECTOR,
        "cadence_authority": "qlp_cam_quat",
        "quality_authority": "spoc_and_qlp_quality_flags",
        "orbits": list(ORBITS),
        "n_spoc_authority_files_verified": 16,
        "n_qlp_qflag_files_verified": 32,
    }
    for field, value in expected.items():
        if manifest.get(field) != value:
            raise ValueError(f"cadence manifest {field} does not match {value!r}")
    if _valid_sha256(manifest.get("table_sha256"), label="cadence table hash") != table_digest:
        raise ValueError("cadence manifest does not bind the cadence table")
    if int(manifest.get("n_rows", 0)) <= 0:
        raise ValueError("cadence manifest contains no rows")
    require_producer_git_sha(
        manifest, expected_producer_git_sha, label="cadence manifest"
    )
    return {
        "contract_version": str(manifest["contract_version"]),
        "n_rows": int(manifest["n_rows"]),
        "table_sha256": table_digest,
        "manifest_sha256": file_sha256(manifest_path),
        "passed": True,
    }


def _validate_compact(
    compact_h5: Path,
    compact_manifest_path: Path,
    *,
    accepted_validation_path: Path,
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    receipt_audit = validate_s63_stage1_compact(
        accepted_validation_path=accepted_validation_path,
        compact_h5_path=compact_h5,
        compact_manifest_path=compact_manifest_path,
        expected_producer_git_sha=expected_producer_git_sha,
    )
    compact_h5 = _require_regular_file(compact_h5, label="compact ADP HDF5")
    manifest = _load_json(compact_manifest_path, label="compact ADP manifest")
    with h5py.File(compact_h5, "r") as h5:
        if int(h5.attrs.get("sector", -1)) != SECTOR:
            raise ValueError("compact HDF5 sector is not 63")
        try:
            columns = json.loads(str(h5.attrs["flux_columns"]))
        except (KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError("compact HDF5 lacks valid flux_columns") from exc
        if columns != list(ADP_ONLY_APERTURES):
            raise ValueError("compact HDF5 is not exactly the locked ADP pair")
        if "targets" not in h5:
            raise ValueError("compact HDF5 lacks /targets")
        n_targets = len(h5["targets"])
        try:
            target_tics = sorted(int(key) for key in h5["targets"].keys())
        except ValueError as exc:
            raise ValueError("compact HDF5 has a non-TIC target group") from exc
        if len(target_tics) != len(set(target_tics)) or any(
            tic <= 0 for tic in target_tics
        ):
            raise ValueError("compact HDF5 has invalid target identities")
        required_datasets = {
            "time",
            "cadenceno",
            "quality",
            "orbitid",
            *ADP_ONLY_APERTURES,
        }
        for key, group in h5["targets"].items():
            missing = sorted(required_datasets - set(group.keys()))
            if missing:
                raise ValueError(f"compact target {key} lacks datasets: {missing}")
            shapes = {name: group[name].shape for name in required_datasets}
            if any(len(shape) != 1 for shape in shapes.values()):
                raise ValueError(f"compact target {key} has non-vector datasets")
            lengths = {shape[0] for shape in shapes.values()}
            if len(lengths) != 1 or next(iter(lengths), 0) <= 0:
                raise ValueError(
                    f"compact target {key} has empty/misaligned dataset lengths"
                )
            tic = int(key)
            expected_attrs = {
                "tic": tic,
                "sector": SECTOR,
            }
            for field, expected in expected_attrs.items():
                if int(group.attrs.get(field, -1)) != expected:
                    raise ValueError(f"compact target {key} has invalid {field}")
            for field in ("camera", "ccd"):
                if int(group.attrs.get(field, -1)) not in range(1, 5):
                    raise ValueError(f"compact target {key} has invalid {field}")
        if n_targets <= 0 or int(h5.attrs.get("n_targets", -1)) != n_targets:
            raise ValueError("compact HDF5 target count is invalid")
    if manifest.get("sector") != SECTOR:
        raise ValueError("compact manifest sector is not 63")
    if manifest.get("requested_columns") != list(ADP_ONLY_APERTURES):
        raise ValueError("compact manifest is not exactly the locked ADP pair")
    if int(manifest.get("n_exported_targets", -1)) != n_targets:
        raise ValueError("compact manifest target count disagrees with HDF5")
    if target_tics != receipt_audit["target_tics"]:
        raise RuntimeError("compact target identities changed during launch validation")
    return {
        "n_targets": n_targets,
        "_target_tics": target_tics,
        "accepted_stage1_validation_sha256": receipt_audit[
            "accepted_stage1_validation_sha256"
        ],
        "hlsp_root": receipt_audit["hlsp_root"],
        "source_host_out_h5": receipt_audit["source_host_out_h5"],
        "relocated_by_hash": receipt_audit["relocated_by_hash"],
        "passed": True,
    }


def _validate_bls(
    *,
    peaks_path: Path,
    summary_path: Path,
    compact_sha256: str,
    cadence_sha256: str,
    cadence_manifest_sha256: str,
    allowlist_sha256: str,
    allowlist_tics: Sequence[int],
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    peaks_path = _require_regular_file(peaks_path, label="merged BLS peaks")
    summary = _load_json(summary_path, label="merged BLS summary")
    require_producer_git_sha(
        summary, expected_producer_git_sha, label="merged BLS summary"
    )
    expected_config = approved_a2v1_teacher_bls_config()
    expected_config_sha256 = bls_config_sha256(expected_config)
    required = {
        "sector": SECTOR,
        "contract_version": ADP_ONLY_CONTRACT_VERSION,
        "bls_search_contract_version": A2V1_TEACHER_BLS_SEARCH_CONTRACT,
        "bls_config_sha256": expected_config_sha256,
        "config": expected_config,
        "apertures": list(ADP_ONLY_APERTURES),
        "n_periods": 50_000,
        "n_peaks": 10,
        "source_product_tag": "A2v1",
        "orbitid_policy": "strict",
        "n_shards": 1,
        "shard_index": 0,
        "passed": True,
        "target_metadata_contract_version": S63_TARGET_METADATA_CONTRACT,
        "target_metadata_columns": list(S63_TARGET_METADATA_COLUMNS),
    }
    for field, value in required.items():
        if summary.get(field) != value:
            raise ValueError(f"BLS summary {field} does not match {value!r}")
    bindings = {
        "peak_table_sha256": file_sha256(peaks_path),
        "compact_lc_sha256": compact_sha256,
        "cadence_reference_sha256": cadence_sha256,
        "cadence_reference_manifest_sha256": cadence_manifest_sha256,
        "target_allowlist_sha256": allowlist_sha256,
    }
    for field, expected in bindings.items():
        if _valid_sha256(summary.get(field), label=f"BLS {field}") != expected:
            raise ValueError(f"BLS summary {field} binding mismatch")
    target_count = len(allowlist_tics)
    for field in ("n_targets", "n_targets_total", "target_allowlist_count"):
        if int(summary.get(field, -1)) != target_count:
            raise ValueError(f"BLS summary {field} disagrees with model-ready allowlist")
    if int(summary.get("n_source_shards", 1)) < 1:
        raise ValueError("BLS summary has invalid source-shard count")
    finite_counts = summary.get("target_metadata_finite_counts")
    if not isinstance(finite_counts, Mapping) or set(finite_counts) != set(
        S63_TARGET_METADATA_COLUMNS
    ):
        raise ValueError("BLS summary lacks exact target-metadata finite counts")
    if any(int(finite_counts[column]) <= 0 for column in S63_TARGET_METADATA_COLUMNS):
        raise ValueError("BLS summary contains a wholly nonfinite model feature")
    return {
        "n_targets": target_count,
        "bls_search_contract_version": A2V1_TEACHER_BLS_SEARCH_CONTRACT,
        "bls_config_sha256": expected_config_sha256,
        "passed": True,
    }


def _validate_model_ready_summary(
    *,
    summary_path: Path,
    expected_hashes: Mapping[str, str],
    n_reserved: int,
    n_model_ready: int,
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    summary = _load_json(summary_path, label="model-ready summary")
    require_producer_git_sha(
        summary, expected_producer_git_sha, label="model-ready summary"
    )
    if summary.get("schema_version") != (
        "twirl_teacher_v3_s63_model_ready_allowlist_v1"
    ):
        raise ValueError("model-ready summary schema mismatch")
    if summary.get("sector") != SECTOR or summary.get("orbits") != list(ORBITS):
        raise ValueError("model-ready summary is not S63 orbits 133/134")
    if summary.get("passed") is not True or summary.get(
        "light_curve_tensors_read"
    ) is not False:
        raise ValueError("model-ready identity-only derivation did not pass")
    source_hashes = summary.get("source_hashes")
    if not isinstance(source_hashes, Mapping):
        raise ValueError("model-ready summary lacks source hashes")
    for field, expected in expected_hashes.items():
        if field == "model_ready_allowlist_sha256":
            continue
        if _valid_sha256(
            source_hashes.get(field), label=f"model-ready source {field}"
        ) != expected:
            raise ValueError(f"model-ready source hash mismatch for {field}")
    counts = summary.get("counts")
    if not isinstance(counts, Mapping):
        raise ValueError("model-ready summary lacks counts")
    if int(counts.get("n_reserved_tics", -1)) != n_reserved or int(
        counts.get("n_model_ready_tics", -1)
    ) != n_model_ready:
        raise ValueError("model-ready summary count mismatch")
    checks = summary.get("partition_checks")
    if not isinstance(checks, Mapping) or not checks or not all(
        value is True for value in checks.values()
    ):
        raise ValueError("model-ready summary partition checks did not pass")
    outputs = summary.get("outputs")
    declaration = (
        outputs.get("model_ready_allowlist")
        if isinstance(outputs, Mapping)
        else None
    )
    if not isinstance(declaration, Mapping) or _valid_sha256(
        declaration.get("sha256"), label="model-ready output hash"
    ) != expected_hashes["model_ready_allowlist_sha256"]:
        raise ValueError("model-ready summary does not bind its allowlist")
    return {
        "n_reserved": n_reserved,
        "n_model_ready": n_model_ready,
        "passed": True,
    }


def _validate_selection_policy(path: Path) -> dict[str, Any]:
    resolved = _require_regular_file(path, label="selection policy")
    policy, _ = load_selection_policy(
        resolved,
        expected_sha256=file_sha256(resolved),
    )
    return {
        "schema_version": str(policy["schema_version"]),
        "selection_seed": 630056,
        "n_primary": 1000,
        "n_repeated_host": 100,
        "passed": True,
    }


def _validate_candidates(
    *,
    candidates_path: Path,
    summary_path: Path,
    model_ready_tics: set[int],
    expected_hashes: Mapping[str, str],
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    candidates_path = _require_regular_file(candidates_path, label="S63 candidates")
    candidates = _read_table(candidates_path, label="S63 candidates")
    summary = _load_json(summary_path, label="S63 candidate summary")
    require_producer_git_sha(
        summary, expected_producer_git_sha, label="candidate summary"
    )
    required_columns = {
        "tic",
        "sector",
        "rep_aperture",
        "rep_peak_rank",
        "period_d",
        "t0_bjd",
        "duration_min",
        "native_input_include",
        *S63_TARGET_METADATA_COLUMNS,
    }
    missing = sorted(required_columns - set(candidates))
    if missing:
        raise ValueError(f"S63 candidate table lacks columns: {missing}")
    sectors = set(pd.to_numeric(candidates["sector"], errors="raise").astype(int))
    if sectors != {SECTOR}:
        raise ValueError("candidate table is not S63-only")
    if set(candidates["rep_aperture"].astype(str)) != {ADP_ONLY_APERTURES[0]}:
        raise ValueError("candidate table is not rank-one ADP-small")
    if set(pd.to_numeric(candidates["rep_peak_rank"], errors="raise").astype(int)) != {1}:
        raise ValueError("candidate table is not rank-one ADP-small")
    candidate_tics = set(
        pd.to_numeric(candidates["tic"], errors="raise").astype(int).tolist()
    )
    if candidate_tics != model_ready_tics:
        raise ValueError(
            "candidate TICs must exactly equal the model-ready allowlist; "
            f"missing={sorted(model_ready_tics - candidate_tics)[:10]}, "
            f"extra={sorted(candidate_tics - model_ready_tics)[:10]}"
        )
    if candidates["tic"].duplicated().any():
        raise ValueError("candidate table must contain at most one row per TIC")
    if summary.get("schema_version") != "twirl_teacher_v3_s63_rank1_candidates_v1":
        raise ValueError("candidate summary schema mismatch")
    if summary.get("sector") != SECTOR:
        raise ValueError("candidate summary is not S63")
    if int(summary.get("n_rows", -1)) != len(model_ready_tics) or int(
        summary.get("n_unique_tics", -1)
    ) != len(model_ready_tics):
        raise ValueError("candidate summary does not cover every model-ready TIC once")
    if summary.get("science_ready") is not False or summary.get(
        "promotion_enabled"
    ) is not False:
        raise ValueError("candidate summary must disable science claims/promotion")
    if summary.get("teacher_v3_metadata_columns") != list(
        S63_TARGET_METADATA_COLUMNS
    ):
        raise ValueError("candidate summary lacks the exact Teacher-v3 metadata set")
    finite_counts = {
        column: int(
            np.isfinite(
                pd.to_numeric(candidates[column], errors="coerce").to_numpy(
                    dtype=float
                )
            ).sum()
        )
        for column in S63_TARGET_METADATA_COLUMNS
    }
    if any(count <= 0 for count in finite_counts.values()):
        raise ValueError(
            "candidate table contains a wholly nonfinite Teacher-v3 feature: "
            f"{finite_counts}"
        )
    _expect_hash(
        summary,
        "candidate_table_sha256",
        file_sha256(candidates_path),
        label="candidate summary",
    )
    for label, (names, expected) in {
        "prospective plan": ("prospective_plan_sha256", expected_hashes["prospective_plan"]),
        "model-ready allowlist": ("model_ready_allowlist_sha256", expected_hashes["allowlist"]),
        "merged BLS peaks": ("bls_merged_table_sha256", expected_hashes["bls_peaks"]),
    }.items():
        _expect_hash_anywhere(
            summary,
            names,
            expected,
            label=f"candidate summary {label}",
        )
    return {
        "n_candidates": int(len(candidates)),
        "n_candidate_tics": int(len(candidate_tics)),
        "_candidate_tics": sorted(candidate_tics),
        "passed": True,
    }


def _validate_raw(
    *,
    raw_h5: Path,
    raw_summary_path: Path,
    candidate_tics: set[int],
    candidate_sha256: str,
    compact_sha256: str,
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    raw_h5 = _require_regular_file(raw_h5, label="raw-source HDF5")
    summary = _load_json(raw_summary_path, label="raw-source summary")
    require_producer_git_sha(
        summary, expected_producer_git_sha, label="raw-source summary"
    )
    with h5py.File(raw_h5, "r") as h5:
        if validate_producer_git_sha(
            h5.attrs.get("producer_git_sha"),
            label="raw-source HDF5 producer_git_sha",
        ) != validate_producer_git_sha(expected_producer_git_sha):
            raise ValueError("raw-source HDF5 producer Git SHA differs from launch")
        if str(h5.attrs.get("contract_version", "")) != RAW_SOURCE_CONTRACT_VERSION:
            raise ValueError("raw-source HDF5 contract mismatch")
        try:
            orbits = sorted(json.loads(str(h5.attrs["orbits"])))
        except (KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError("raw-source HDF5 lacks valid orbit provenance") from exc
        if orbits != list(ORBITS):
            raise ValueError("raw-source HDF5 is not bound to orbits 133/134")
        if str(h5.attrs.get("training_table_sha256", "")) != candidate_sha256:
            raise ValueError("raw-source HDF5 does not bind the candidate table")
        if str(h5.attrs.get("compact_adp_h5_sha256", "")) != compact_sha256:
            raise ValueError("raw-source HDF5 does not bind the compact ADP input")
        observed_tics = (
            {int(key) for key in h5["targets"]} if "targets" in h5 else set()
        )
        if observed_tics != candidate_tics:
            raise ValueError("raw-source HDF5 target coverage is not exact")
    if int(summary.get("n_requested", -1)) != len(candidate_tics) or int(
        summary.get("n_written", -1)
    ) != len(candidate_tics):
        raise ValueError("raw-source summary target coverage is not exact")
    if _valid_sha256(
        summary.get("out_h5_sha256"), label="raw-source output hash"
    ) != file_sha256(raw_h5):
        raise ValueError("raw-source summary does not bind its HDF5")
    for field, expected in (
        ("training_table_sha256", candidate_sha256),
        ("compact_adp_h5_sha256", compact_sha256),
    ):
        if _valid_sha256(summary.get(field), label=f"raw-source {field}") != expected:
            raise ValueError(f"raw-source summary {field} binding mismatch")
    return {"n_targets": len(candidate_tics), "passed": True}


def _validate_native(
    *,
    native_h5: Path,
    native_summary_path: Path,
    candidate_sha256: str,
    raw_sha256: str,
    compact_sha256: str,
    cadence_sha256: str,
    cadence_manifest_sha256: str,
    candidate_tics: int,
    expected_candidate_tics: set[int],
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    native_h5 = _require_regular_file(native_h5, label="native HDF5")
    summary = _load_json(native_summary_path, label="native summary")
    require_producer_git_sha(
        summary, expected_producer_git_sha, label="native summary"
    )
    if summary.get("authorization_contract") != AUTHORIZATION_CONTRACT:
        raise ValueError("native summary lacks the prospective S63 authorization")
    if summary.get("sector") != SECTOR:
        raise ValueError("native summary is not S63")
    if summary.get("verification", {}).get("passed") is not True:
        raise ValueError("native summary verification did not pass")
    if _valid_sha256(summary.get("out_h5_sha256"), label="native output hash") != file_sha256(native_h5):
        raise ValueError("native summary does not bind the native HDF5")
    if summary.get("exact_count_match") is not True or summary.get(
        "exact_group_identity_match"
    ) is not True:
        raise ValueError("native merge lacks exact count/group identity")
    with h5py.File(native_h5, "r") as h5:
        if validate_producer_git_sha(
            h5.attrs.get("producer_git_sha"),
            label="native HDF5 producer_git_sha",
        ) != validate_producer_git_sha(expected_producer_git_sha):
            raise ValueError("native HDF5 producer Git SHA differs from launch")
        expected_attrs = {
            "contract_version": RAW_PAIR_CONTRACT_VERSION,
            "training_table_sha256": candidate_sha256,
            "raw_source_h5_sha256": raw_sha256,
            "compact_adp_h5_sha256": compact_sha256,
            "cadence_reference_table_sha256": cadence_sha256,
            "cadence_reference_manifest_sha256": cadence_manifest_sha256,
            "orbitid_reconciliation_policy": "strict",
        }
        for field, expected in expected_attrs.items():
            if str(h5.attrs.get(field, "")) != expected:
                raise ValueError(f"native HDF5 {field} binding mismatch")
        if int(h5.attrs.get("periodogram_n", -1)) != PERIODOGRAM_N:
            raise ValueError("native HDF5 periodogram length is not 4096")
        observed_target_tics = (
            {int(key) for key in h5["targets"]} if "targets" in h5 else set()
        )
        if observed_target_tics != expected_candidate_tics:
            raise ValueError("native HDF5 target coverage is not exact")
        if "injections" not in h5 or len(h5["injections"]) != 0:
            raise ValueError("prospective S63 native HDF5 must be real-only")
    return {"n_targets": candidate_tics, "periodogram_n": PERIODOGRAM_N, "passed": True}


def _git_identity(repo: Path, expected_git_sha: str) -> dict[str, Any]:
    repo = Path(repo).expanduser().resolve(strict=True)
    expected = str(expected_git_sha).strip().lower()
    if GIT_SHA_RE.fullmatch(expected) is None:
        raise ValueError("expected Git SHA must be a full lowercase 40-character SHA")
    observed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip().lower()
    if observed != expected:
        raise ValueError(f"checkout Git SHA {observed} != expected {expected}")
    dirty = subprocess.run(
        ["git", "status", "--short", "--untracked-files=all"],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if dirty:
        raise ValueError("launch manifest requires a fully clean checkout")
    return {
        "repo": str(repo),
        "sha": observed,
        "checkout_clean": True,
        "untracked_files_checked": True,
    }


def _publish_json_atomic(payload: Mapping[str, Any], output: Path) -> None:
    output = Path(output).expanduser().resolve()
    claim = output.with_suffix(output.suffix + ".claim")
    temporary = output.with_suffix(output.suffix + ".tmp")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() or temporary.exists():
        raise FileExistsError(f"refusing to overwrite launch output/pending file: {output}")
    descriptor: int | None = None
    try:
        descriptor = os.open(claim, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        os.write(descriptor, f"pid={os.getpid()}\n".encode("ascii"))
        os.close(descriptor)
        descriptor = None
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        temporary.unlink(missing_ok=True)
        claim.unlink(missing_ok=True)


def build_launch_manifest(args: argparse.Namespace) -> dict[str, Any]:
    git = _git_identity(args.repo, args.expected_git_sha)
    producer_git_sha = git["sha"]
    artifacts = {
        "preregistered_contract": _artifact(args.preregistered_contract, label="preregistered contract"),
        "selection_policy": _artifact(args.selection_policy, label="selection policy"),
        "accepted_stage1_validation": _artifact(args.accepted_validation, label="accepted Stage-1 validation"),
        "cadence_table": _artifact(args.cadence_table, label="cadence table"),
        "cadence_manifest": _artifact(args.cadence_manifest, label="cadence manifest"),
        "compact_h5": _artifact(args.compact_h5, label="compact ADP HDF5"),
        "compact_manifest": _artifact(args.compact_manifest, label="compact ADP manifest"),
        "reserved_tics": _artifact(args.reserved_tics, label="reserved TIC inventory"),
        "teacher_v3_corpus": _artifact(args.teacher_v3_corpus, label="frozen Teacher-v3 corpus"),
        "model_ready_allowlist": _artifact(args.model_ready_allowlist, label="model-ready allowlist"),
        "model_ready_summary": _artifact(args.model_ready_summary, label="model-ready summary"),
        "primary_cohort": _artifact(args.primary_cohort, label="primary prospective cohort"),
        "repeated_host_cohort": _artifact(args.repeated_host_cohort, label="repeated-host cohort"),
        "reserved_not_model_ready": _artifact(args.reserved_not_model_ready, label="reserved-not-model-ready inventory"),
        "cohort_summary": _artifact(args.cohort_summary, label="prospective cohort summary"),
        "bls_peaks": _artifact(args.bls_peaks, label="merged BLS peaks"),
        "bls_summary": _artifact(args.bls_summary, label="merged BLS summary"),
        "candidates": _artifact(args.candidates, label="S63 candidates"),
        "candidate_summary": _artifact(args.candidate_summary, label="S63 candidate summary"),
        "raw_source_h5": _artifact(args.raw_source_h5, label="raw-source HDF5"),
        "raw_source_summary": _artifact(args.raw_source_summary, label="raw-source summary"),
        "native_h5": _artifact(args.native_h5, label="native HDF5"),
        "native_summary": _artifact(args.native_summary, label="native summary"),
    }
    preregistered, _ = load_prospective_plan(
        args.preregistered_contract,
        expected_sha256=artifacts["preregistered_contract"]["sha256"],
    )
    selection_policy = _validate_selection_policy(args.selection_policy)
    stage1 = validate_accepted_stage1(
        args.accepted_validation,
        expected_producer_git_sha=producer_git_sha,
    )
    cadence = _validate_cadence(
        args.cadence_table,
        args.cadence_manifest,
        expected_producer_git_sha=producer_git_sha,
    )
    compact = _validate_compact(
        args.compact_h5,
        args.compact_manifest,
        accepted_validation_path=args.accepted_validation,
        expected_producer_git_sha=producer_git_sha,
    )
    reserved_tics = _read_tic_inventory(args.reserved_tics, label="reserved TIC inventory")
    reservation = preregistered["s63_identity_reservation"]
    if reservation.get("reserved_tics_sha256") != artifacts["reserved_tics"]["sha256"]:
        raise ValueError("prospective plan does not bind the sealed S63 TIC file")
    if int(reservation.get("n_requested_tics", -1)) != len(reserved_tics):
        raise ValueError("prospective plan S63 reservation count mismatch")
    frozen_training = preregistered["frozen_training_identity"]
    if frozen_training.get("morphology_corpus_sha256") != artifacts[
        "teacher_v3_corpus"
    ]["sha256"]:
        raise ValueError("prospective plan does not bind the frozen Teacher-v3 corpus")
    teacher_v3_corpus_tics = _read_tic_inventory(
        args.teacher_v3_corpus,
        label="frozen Teacher-v3 corpus",
        allow_duplicate_rows=True,
    )
    if int(frozen_training.get("n_corpus_tics", -1)) != len(
        teacher_v3_corpus_tics
    ):
        raise ValueError("prospective plan Teacher-v3 corpus count mismatch")
    model_ready_tics = _read_tic_inventory(
        args.model_ready_allowlist, label="model-ready allowlist"
    )
    if not set(model_ready_tics).issubset(set(reserved_tics)):
        raise ValueError("model-ready allowlist is not a subset of sealed S63 identities")
    if len(model_ready_tics) != compact["n_targets"]:
        raise ValueError("model-ready allowlist does not exactly cover the compact export")
    if set(model_ready_tics) != set(compact.pop("_target_tics")):
        raise ValueError("model-ready allowlist identities differ from compact targets")
    model_ready = _validate_model_ready_summary(
        summary_path=args.model_ready_summary,
        expected_hashes={
            "prospective_plan_sha256": artifacts["preregistered_contract"]["sha256"],
            "accepted_stage1_validation_sha256": artifacts["accepted_stage1_validation"]["sha256"],
            "compact_h5_sha256": artifacts["compact_h5"]["sha256"],
            "compact_manifest_sha256": artifacts["compact_manifest"]["sha256"],
            "reserved_tics_sha256": artifacts["reserved_tics"]["sha256"],
            "model_ready_allowlist_sha256": artifacts["model_ready_allowlist"]["sha256"],
        },
        n_reserved=len(reserved_tics),
        n_model_ready=len(model_ready_tics),
        expected_producer_git_sha=producer_git_sha,
    )
    primary_tics = _read_tic_inventory(
        args.primary_cohort, label="primary prospective cohort"
    )
    repeated_tics = _read_tic_inventory(
        args.repeated_host_cohort, label="repeated-host cohort"
    )
    reserved_not_ready_tics = _read_tic_inventory(
        args.reserved_not_model_ready,
        label="reserved-not-model-ready inventory",
    )
    if set(primary_tics) & set(repeated_tics):
        raise ValueError("primary and repeated-host prospective cohorts overlap")
    if set(primary_tics) | set(repeated_tics) != set(model_ready_tics):
        raise ValueError("prospective cohorts do not exactly partition model-ready TICs")
    expected_repeated = set(model_ready_tics) & set(teacher_v3_corpus_tics)
    expected_primary = set(model_ready_tics) - set(teacher_v3_corpus_tics)
    if set(primary_tics) != expected_primary or set(repeated_tics) != expected_repeated:
        raise ValueError(
            "prospective cohort identities disagree with the frozen Teacher-v3 corpus"
        )
    if set(reserved_not_ready_tics) & set(model_ready_tics):
        raise ValueError("reserved-not-model-ready inventory overlaps model-ready TICs")
    if set(reserved_not_ready_tics) | set(model_ready_tics) != set(reserved_tics):
        raise ValueError("model-ready plus excluded inventories do not recover the seal")
    cohort_summary = _load_json(args.cohort_summary, label="prospective cohort summary")
    require_producer_git_sha(
        cohort_summary, producer_git_sha, label="prospective cohort summary"
    )
    if cohort_summary.get("schema_version") != "twirl_teacher_v3_s63_cohorts_v1":
        raise ValueError("prospective cohort summary schema mismatch")
    partition_checks = cohort_summary.get("partition_checks")
    if not isinstance(partition_checks, Mapping) or not partition_checks or not all(
        value is True for value in partition_checks.values()
    ):
        raise ValueError("prospective cohort partition checks did not all pass")
    cohort_outputs = cohort_summary.get("outputs")
    if not isinstance(cohort_outputs, Mapping):
        raise ValueError("prospective cohort summary lacks output provenance")
    cohort_sources = cohort_summary.get("source_hashes")
    if not isinstance(cohort_sources, Mapping):
        raise ValueError("prospective cohort summary lacks source hashes")
    expected_cohort_sources = {
        "model_ready_allowlist_sha256": artifacts["model_ready_allowlist"]["sha256"],
        "teacher_v3_corpus_sha256": artifacts["teacher_v3_corpus"]["sha256"],
        "reserved_tics_sha256": artifacts["reserved_tics"]["sha256"],
    }
    for field, expected in expected_cohort_sources.items():
        if _valid_sha256(
            cohort_sources.get(field), label=f"cohort source {field}"
        ) != expected:
            raise ValueError(f"prospective cohort source hash mismatch for {field}")
    for key, artifact_name in (
        ("primary_cohort", "primary_cohort"),
        ("repeated_host_cohort", "repeated_host_cohort"),
        ("reserved_not_model_ready", "reserved_not_model_ready"),
    ):
        declaration = cohort_outputs.get(key)
        if not isinstance(declaration, Mapping):
            raise ValueError(f"cohort summary lacks outputs.{key}")
        if _valid_sha256(
            declaration.get("sha256"), label=f"cohort outputs.{key}.sha256"
        ) != artifacts[artifact_name]["sha256"]:
            raise ValueError(f"cohort summary outputs.{key} hash mismatch")
    bls = _validate_bls(
        peaks_path=args.bls_peaks,
        summary_path=args.bls_summary,
        compact_sha256=artifacts["compact_h5"]["sha256"],
        cadence_sha256=artifacts["cadence_table"]["sha256"],
        cadence_manifest_sha256=artifacts["cadence_manifest"]["sha256"],
        allowlist_sha256=artifacts["model_ready_allowlist"]["sha256"],
        allowlist_tics=model_ready_tics,
        expected_producer_git_sha=producer_git_sha,
    )
    candidates = _validate_candidates(
        candidates_path=args.candidates,
        summary_path=args.candidate_summary,
        model_ready_tics=set(model_ready_tics),
        expected_hashes={
            "prospective_plan": artifacts["preregistered_contract"]["sha256"],
            "bls_peaks": artifacts["bls_peaks"]["sha256"],
            "allowlist": artifacts["model_ready_allowlist"]["sha256"],
        },
        expected_producer_git_sha=producer_git_sha,
    )
    raw = _validate_raw(
        raw_h5=args.raw_source_h5,
        raw_summary_path=args.raw_source_summary,
        candidate_tics=set(candidates["_candidate_tics"]),
        candidate_sha256=artifacts["candidates"]["sha256"],
        compact_sha256=artifacts["compact_h5"]["sha256"],
        expected_producer_git_sha=producer_git_sha,
    )
    native = _validate_native(
        native_h5=args.native_h5,
        native_summary_path=args.native_summary,
        candidate_sha256=artifacts["candidates"]["sha256"],
        raw_sha256=artifacts["raw_source_h5"]["sha256"],
        compact_sha256=artifacts["compact_h5"]["sha256"],
        cadence_sha256=artifacts["cadence_table"]["sha256"],
        cadence_manifest_sha256=artifacts["cadence_manifest"]["sha256"],
        candidate_tics=candidates["n_candidate_tics"],
        expected_candidate_tics=set(candidates.pop("_candidate_tics")),
        expected_producer_git_sha=producer_git_sha,
    )
    manifest = {
        "contract_version": CONTRACT_VERSION,
        "authorization_contract": AUTHORIZATION_CONTRACT,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": SECTOR,
        "orbits": list(ORBITS),
        "status": "preprocessing_complete_scoring_not_started",
        "score_or_queue_read": False,
        "git": git,
        "producer_git_sha": producer_git_sha,
        "artifacts": artifacts,
        "checks": {
            "accepted_stage1": stage1,
            "selection_policy": selection_policy,
            "cadence": cadence,
            "compact": compact,
            "inventories": {
                "n_reserved_tics": len(reserved_tics),
                "n_model_ready_tics": len(model_ready_tics),
                "cohort_counts": {
                    "primary_corpus_disjoint": len(primary_tics),
                    "repeated_host_side": len(repeated_tics),
                    "reserved_not_model_ready": len(reserved_not_ready_tics),
                },
                "passed": True,
            },
            "model_ready_derivation": model_ready,
            "bls": bls,
            "candidates": candidates,
            "raw_source": raw,
            "native": native,
        },
        "locked_search": {
            "apertures": list(ADP_ONLY_APERTURES),
            "n_periods": 50_000,
            "n_peaks": 10,
            "orbitid_policy": "strict",
            "native_periodogram_n": PERIODOGRAM_N,
        },
        "passed": True,
    }
    # Every hash published by this manifest is resolved by construction.
    for label, artifact in artifacts.items():
        _valid_sha256(artifact["sha256"], label=f"manifest {label} hash")
    return manifest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted-validation", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--preregistered-contract", type=Path)
    parser.add_argument("--selection-policy", type=Path)
    parser.add_argument("--cadence-table", type=Path)
    parser.add_argument("--cadence-manifest", type=Path)
    parser.add_argument("--compact-h5", type=Path)
    parser.add_argument("--compact-manifest", type=Path)
    parser.add_argument("--reserved-tics", type=Path)
    parser.add_argument("--teacher-v3-corpus", type=Path)
    parser.add_argument("--model-ready-allowlist", type=Path)
    parser.add_argument("--model-ready-summary", type=Path)
    parser.add_argument("--primary-cohort", type=Path)
    parser.add_argument("--repeated-host-cohort", type=Path)
    parser.add_argument("--reserved-not-model-ready", type=Path)
    parser.add_argument("--cohort-summary", type=Path)
    parser.add_argument("--bls-peaks", type=Path)
    parser.add_argument("--bls-summary", type=Path)
    parser.add_argument("--candidates", type=Path)
    parser.add_argument("--candidate-summary", type=Path)
    parser.add_argument("--raw-source-h5", type=Path)
    parser.add_argument("--raw-source-summary", type=Path)
    parser.add_argument("--native-h5", type=Path)
    parser.add_argument("--native-summary", type=Path)
    parser.add_argument("--repo", type=Path, default=ROOT)
    parser.add_argument("--expected-git-sha")
    parser.add_argument("--out", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if args.preflight_only:
        expected_preflight_sha = None
        if args.expected_git_sha is not None:
            expected_preflight_sha = validate_producer_git_sha(
                args.expected_git_sha,
                label="preflight producer Git SHA",
            )
        print(
            json.dumps(
                validate_accepted_stage1(
                    args.accepted_validation,
                    expected_producer_git_sha=expected_preflight_sha,
                ),
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
        )
        return 0
    required = (
        "preregistered_contract",
        "selection_policy",
        "cadence_table",
        "cadence_manifest",
        "compact_h5",
        "compact_manifest",
        "reserved_tics",
        "teacher_v3_corpus",
        "model_ready_allowlist",
        "model_ready_summary",
        "primary_cohort",
        "repeated_host_cohort",
        "reserved_not_model_ready",
        "cohort_summary",
        "bls_peaks",
        "bls_summary",
        "candidates",
        "candidate_summary",
        "raw_source_h5",
        "raw_source_summary",
        "native_h5",
        "native_summary",
        "expected_git_sha",
        "out",
    )
    missing = [name for name in required if getattr(args, name) in (None, "")]
    if missing:
        parser.error("full launch manifest requires: " + ", ".join(missing))
    manifest = build_launch_manifest(args)
    _publish_json_atomic(manifest, args.out)
    print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
