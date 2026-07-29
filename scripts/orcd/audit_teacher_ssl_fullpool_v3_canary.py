#!/usr/bin/env python3
"""Fail-closed model-facing audit for the fresh S60/native_4 v3r1 canary."""
from __future__ import annotations

import argparse
from collections.abc import Mapping
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.harmonic_dataset import (  # noqa: E402
    HarmonicNativeDataset,
    collate_native_samples,
)
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    MODEL_INPUT_CONTRACT_VERSION,
)
from twirl.vetting.ssl_full_pool_eligibility import (  # noqa: E402
    PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_FROZEN_POOL_CSV_SHA256,
    PRODUCTION_FROZEN_POOL_SUMMARY_SHA256,
    PRODUCTION_FULL_IDENTITY_SHA256,
    PRODUCTION_FULL_OBSERVATIONS,
    PRODUCTION_GLOBAL_BLS_SHA256,
    PRODUCTION_GLOBAL_BLS_SUMMARY_SHA256,
    derive_anchor_eligibility,
    observation_identity_sha256,
)
from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION,
    FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D,
    FULL_POOL_NATIVE_CONTRACT_VERSION,
    FULL_POOL_NATIVE_DETREND_CONFIG_SHA256,
    FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION,
    FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION,
    FULL_POOL_NATIVE_DETREND_TIME_DATASET,
    FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON,
    FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256,
    FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY,
    FULL_POOL_NATIVE_RELEASE_BINDING,
    FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION,
    FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY,
    full_pool_native_group_failures,
    shard_for_tic,
    verify_full_pool_native_shard,
)
from twirl.vetting.ssl_full_pool_numeric import (  # noqa: E402
    FULL_POOL_NUMERIC_ENVELOPE_V1,
    MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
    audit_collated_sample,
)
from twirl.vetting.teacher_native_registry import file_sha256  # noqa: E402
from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    FULLPOOL_SSL_PROFILE,
    fullpool_dataset_rows,
    load_frozen_ssl_full_pool,
    load_global_full_pool_bls,
)


CANARY_SCHEMA = (
    "twirl_teacher_ssl_fullpool_v3r1_s60_shard4_model_facing_canary_v2"
)
CANARY_SECTOR = 60
CANARY_SHARD_INDEX = 4
CANARY_N_SHARDS = 16
CANARY_TARGET_TIC = 704_538_814
AUDIT_BATCH_SIZE = 8


def _artifact_binding(path: Path) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve(strict=True)
    before = resolved.stat()
    if not resolved.is_file() or before.st_size <= 0:
        raise ValueError(f"authority is missing or empty: {resolved}")
    digest = file_sha256(resolved)
    after = resolved.stat()
    identity_before = (
        before.st_size,
        before.st_mtime_ns,
        before.st_dev,
        before.st_ino,
    )
    identity_after = (
        after.st_size,
        after.st_mtime_ns,
        after.st_dev,
        after.st_ino,
    )
    if identity_before != identity_after:
        raise RuntimeError(f"authority changed while hashing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": int(after.st_size),
        "sha256": digest,
    }


def _assert_unchanged(binding: Mapping[str, Any]) -> None:
    if _artifact_binding(Path(str(binding["path"]))) != dict(binding):
        raise RuntimeError(f"authority changed during canary audit: {binding['path']}")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid JSON authority {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"JSON authority must be an object: {path}")
    return payload


def _code_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    revision = completed.stdout.strip()
    if (
        completed.returncode != 0
        or len(revision) != 40
        or any(value not in "0123456789abcdef" for value in revision)
    ):
        raise RuntimeError("cannot bind canary audit to a Git revision")
    return revision


def _validate_sha_sidecar(
    sidecar_path: Path,
    *,
    artifact_path: Path,
    artifact_sha256: str,
) -> None:
    lines = Path(sidecar_path).read_text(encoding="ascii").splitlines()
    if len(lines) != 1:
        raise ValueError("native canary SHA sidecar must contain exactly one line")
    parts = lines[0].split(maxsplit=1)
    if len(parts) != 2 or parts[0] != artifact_sha256:
        raise ValueError("native canary SHA sidecar digest differs")
    declared = parts[1].lstrip("*").strip()
    if Path(declared).expanduser().resolve() != Path(artifact_path).resolve():
        raise ValueError("native canary SHA sidecar path differs")


def _identity(frame: pd.DataFrame) -> str:
    return observation_identity_sha256(frame.loc[:, ["sector", "tic"]])


def _validate_summary(
    summary: Mapping[str, Any],
    *,
    native_h5: Path,
    native_sha256: str,
    expected_code_revision: str,
    source_rows: pd.DataFrame,
    eligible_rows: pd.DataFrame,
    excluded_rows: pd.DataFrame,
) -> None:
    expected = {
        "schema_version": FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION,
        "sector": CANARY_SECTOR,
        "shard_index": CANARY_SHARD_INDEX,
        "n_shards": CANARY_N_SHARDS,
        "n_source_shard_observations": len(source_rows),
        "n_shard_observations": len(eligible_rows),
        "n_shard_excluded_observations": len(excluded_rows),
        "source_shard_observation_identity_sha256": _identity(source_rows),
        "shard_observation_identity_sha256": _identity(eligible_rows),
        "shard_excluded_observation_identity_sha256": _identity(excluded_rows),
        "expected_git_sha": expected_code_revision,
        "builder_contract_version": FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION,
        "detrend_contract_version": FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION,
        "detrend_config_sha256": FULL_POOL_NATIVE_DETREND_CONFIG_SHA256,
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
        "native_model_full_identity_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
        "native_model_eligible_identity_sha256": (
            PRODUCTION_ELIGIBLE_IDENTITY_SHA256
        ),
        "out_h5": str(native_h5.resolve()),
        "out_h5_sha256": native_sha256,
        "native_contract_version": FULL_POOL_NATIVE_CONTRACT_VERSION,
        "real_only": True,
    }
    mismatches = {
        name: {"expected": value, "observed": summary.get(name)}
        for name, value in expected.items()
        if summary.get(name) != value
    }
    if mismatches:
        raise ValueError(f"native canary summary differs: {mismatches}")
    if summary.get("verification", {}).get("passed") is not True:
        raise ValueError("native canary summary verification did not pass")


def _validate_h5_inventory(
    native_h5: Path,
    *,
    expected_rows: pd.DataFrame,
) -> list[str]:
    import h5py

    expected_groups = {
        f"targets/{int(tic):016d}"
        for tic in expected_rows["tic"].astype(int)
    }
    failures: list[str] = []
    with h5py.File(native_h5, "r") as h5:
        if "targets" not in h5:
            raise ValueError("native canary lacks the targets group")
        observed_groups = {f"targets/{key}" for key in h5["targets"].keys()}
        if observed_groups != expected_groups:
            raise ValueError(
                "native canary group inventory differs from the primary "
                f"eligibility authority: missing={sorted(expected_groups - observed_groups)[:10]}, "
                f"extra={sorted(observed_groups - expected_groups)[:10]}"
            )
        target_group = f"targets/{CANARY_TARGET_TIC:016d}"
        if target_group not in observed_groups:
            raise ValueError("native canary omitted required TIC 704538814")
        root_policy = str(h5.attrs.get("orbitid_reconciliation_policy", ""))
        root_contract = str(h5.attrs.get("contract_version", ""))
        expected_by_group = {
            f"targets/{int(row.tic):016d}": (CANARY_SECTOR, int(row.tic))
            for row in expected_rows.itertuples(index=False)
        }
        for group_path, expected_key in expected_by_group.items():
            group = h5[group_path]
            actual_key = (
                int(group.attrs.get("sector", -1)),
                int(group.attrs.get("tic", -1)),
            )
            if actual_key != expected_key:
                failures.append(
                    f"/{group_path}: identity {actual_key} != {expected_key}"
                )
            group_failures, _ = full_pool_native_group_failures(
                group,
                context=f"/{group_path}",
                root_policy=root_policy,
                root_contract=root_contract,
            )
            failures.extend(group_failures)
    return failures


def _audit_dataset(
    dataset: HarmonicNativeDataset,
    selected: pd.DataFrame,
) -> list[dict[str, Any]]:
    reports: list[dict[str, Any]] = []
    samples: list[dict[str, Any]] = []
    rows: list[Any] = []

    def flush() -> None:
        if not samples:
            return
        batch = collate_native_samples(samples)
        for sample_index, row in enumerate(rows):
            result = audit_collated_sample(batch, sample_index=sample_index)
            if (
                result.get("review_id") != str(row.ssl_observation_id)
                or int(result.get("tic", -1)) != int(row.tic)
            ):
                raise RuntimeError(
                    "canary collator changed the observation identity"
                )
            reports.append(
                {
                    "ssl_observation_id": str(row.ssl_observation_id),
                    "sector": int(row.sector),
                    "tic": int(row.tic),
                    "passed": bool(result["passed"]),
                    "n_failures": int(result["n_failures"]),
                    "failures": list(result["failures"]),
                    "action": str(result["action"]),
                }
            )
        samples.clear()
        rows.clear()

    for index, row in enumerate(selected.itertuples(index=False)):
        try:
            sample = dataset[index]
        except Exception as exc:
            flush()
            reports.append(
                {
                    "ssl_observation_id": str(row.ssl_observation_id),
                    "sector": int(row.sector),
                    "tic": int(row.tic),
                    "passed": False,
                    "n_failures": 1,
                    "failures": [
                        {
                            "code": "transform_exception",
                            "type": type(exc).__name__,
                            "message": str(exc),
                        }
                    ],
                    "action": "audit_only_no_clip_no_exclusion",
                }
            )
            continue
        samples.append(sample)
        rows.append(row)
        if len(samples) == AUDIT_BATCH_SIZE:
            flush()
    flush()
    return reports


def _json_bytes(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _publish_new(path: Path, payload: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
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
        os.link(temporary, path)
    except FileExistsError as exc:
        raise FileExistsError(
            f"refusing to overwrite or reuse canary audit artifact: {path}"
        ) from exc
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def audit_canary(
    *,
    pool_path: Path,
    pool_summary_path: Path,
    bls_path: Path,
    bls_summary_path: Path,
    native_h5_path: Path,
    native_summary_path: Path,
    native_sha_sidecar_path: Path,
    report_out: Path,
    expected_code_revision: str,
) -> dict[str, Any]:
    if (
        len(expected_code_revision) != 40
        or any(value not in "0123456789abcdef" for value in expected_code_revision)
        or expected_code_revision != _code_revision()
    ):
        raise ValueError("canary audit repository is not the expected revision")

    paths = {
        "frozen_pool": pool_path,
        "frozen_pool_summary": pool_summary_path,
        "global_bls": bls_path,
        "global_bls_summary": bls_summary_path,
        "native_h5": native_h5_path,
        "native_summary": native_summary_path,
        "native_sha_sidecar": native_sha_sidecar_path,
    }
    bindings = {name: _artifact_binding(path) for name, path in paths.items()}
    expected_primary_hashes = {
        "frozen_pool": PRODUCTION_FROZEN_POOL_CSV_SHA256,
        "frozen_pool_summary": PRODUCTION_FROZEN_POOL_SUMMARY_SHA256,
        "global_bls": PRODUCTION_GLOBAL_BLS_SHA256,
        "global_bls_summary": PRODUCTION_GLOBAL_BLS_SUMMARY_SHA256,
    }
    for name, expected in expected_primary_hashes.items():
        if bindings[name]["sha256"] != expected:
            raise ValueError(f"canary primary authority hash differs: {name}")

    pool, _ = load_frozen_ssl_full_pool(
        pool_path=Path(bindings["frozen_pool"]["path"]),
        summary_path=Path(bindings["frozen_pool_summary"]["path"]),
    )
    if (
        len(pool) != PRODUCTION_FULL_OBSERVATIONS
        or _identity(pool) != PRODUCTION_FULL_IDENTITY_SHA256
    ):
        raise ValueError("canary frozen-pool identity differs from production")
    bls, _, _ = load_global_full_pool_bls(
        summary_path=Path(bindings["global_bls_summary"]["path"]),
        frozen_pool=pool,
        frozen_pool_summary_path=Path(
            bindings["frozen_pool_summary"]["path"]
        ),
        output_path_override=Path(bindings["global_bls"]["path"]),
    )
    decisions = derive_anchor_eligibility(bls, frozen_pool=pool)
    eligible_all = decisions.loc[
        decisions["native_model_eligible"].astype(bool)
    ]
    if (
        len(decisions) != PRODUCTION_FULL_OBSERVATIONS
        or _identity(decisions) != PRODUCTION_FULL_IDENTITY_SHA256
        or len(eligible_all) != PRODUCTION_ELIGIBLE_OBSERVATIONS
        or _identity(eligible_all) != PRODUCTION_ELIGIBLE_IDENTITY_SHA256
    ):
        raise ValueError("derived canary eligibility identity differs")

    shard_mask = decisions["sector"].astype(int).eq(CANARY_SECTOR) & (
        decisions["tic"]
        .astype(int)
        .map(
            lambda tic: shard_for_tic(
                sector=CANARY_SECTOR,
                tic=int(tic),
                n_shards=CANARY_N_SHARDS,
            )
        )
        .eq(CANARY_SHARD_INDEX)
    )
    source_rows = decisions.loc[shard_mask].sort_values(
        ["sector", "tic"], kind="stable"
    )
    eligible_rows = source_rows.loc[
        source_rows["native_model_eligible"].astype(bool)
    ].copy()
    excluded_rows = source_rows.loc[
        ~source_rows["native_model_eligible"].astype(bool)
    ].copy()
    if eligible_rows.empty:
        raise ValueError("S60 shard-4 primary eligible inventory is empty")
    target = eligible_rows.loc[
        eligible_rows["tic"].astype(int).eq(CANARY_TARGET_TIC)
    ]
    if len(target) != 1:
        raise ValueError("S60 shard-4 lacks exactly one eligible TIC 704538814")
    target_row = target.iloc[0]
    if (
        str(target_row["aperture"]) != "DET_FLUX_ADP_SML"
        or int(target_row["peak_rank"]) != 1
        or str(target_row["status"]).lower() != "ok"
        or not np.isfinite(float(target_row["period_d"]))
        or float(target_row["period_d"]) <= 0
        or not np.isfinite(float(target_row["t0_bjd"]))
        or not np.isfinite(float(target_row["duration_min"]))
        or float(target_row["duration_min"]) <= 0
    ):
        raise ValueError("TIC 704538814 is not the eligible ADP-small rank-1 anchor")

    native_h5 = Path(bindings["native_h5"]["path"])
    native_sha = str(bindings["native_h5"]["sha256"])
    _validate_sha_sidecar(
        Path(bindings["native_sha_sidecar"]["path"]),
        artifact_path=native_h5,
        artifact_sha256=native_sha,
    )
    native_summary = _read_json(Path(bindings["native_summary"]["path"]))
    _validate_summary(
        native_summary,
        native_h5=native_h5,
        native_sha256=native_sha,
        expected_code_revision=expected_code_revision,
        source_rows=source_rows,
        eligible_rows=eligible_rows,
        excluded_rows=excluded_rows,
    )
    native_verification = verify_full_pool_native_shard(
        native_h5,
        expected_git_sha=expected_code_revision,
        expected_sector=CANARY_SECTOR,
        expected_shard_index=CANARY_SHARD_INDEX,
        expected_n_shards=CANARY_N_SHARDS,
        expected_observations=len(eligible_rows),
        expected_source_observations=len(source_rows),
        expected_excluded_observations=len(excluded_rows),
        expected_source_identity_sha256=_identity(source_rows),
        expected_output_identity_sha256=_identity(eligible_rows),
        expected_excluded_identity_sha256=_identity(excluded_rows),
    )
    if native_verification.get("passed") is not True:
        raise ValueError(
            "S60 shard-4 native canary verification failed: "
            + "; ".join(native_verification.get("failures", [])[:20])
        )
    group_failures = _validate_h5_inventory(
        native_h5,
        expected_rows=eligible_rows,
    )
    if group_failures:
        raise ValueError(
            "S60 shard-4 group provenance failed: "
            + "; ".join(group_failures[:20])
        )

    selected = pd.DataFrame(
        {
            "ssl_observation_id": eligible_rows["observation_id"].astype(str),
            "sector": eligible_rows["sector"].astype(np.int64),
            "tic": eligible_rows["tic"].astype(np.int64),
            "period_d": eligible_rows["period_d"].astype(float),
            "t0_bjd": eligible_rows["t0_bjd"].astype(float),
            "duration_min": eligible_rows["duration_min"].astype(float),
            "native_h5_path": str(native_h5),
            "native_group_path": eligible_rows["tic"].map(
                lambda tic: f"targets/{int(tic):016d}"
            ),
            "native_h5_sha256": native_sha,
            "native_contract_version": FULL_POOL_NATIVE_CONTRACT_VERSION,
            "ssl_held_out_fold": np.int8(0),
        }
    ).sort_values(["sector", "tic"], kind="stable").reset_index(drop=True)
    dataset_rows = fullpool_dataset_rows(selected)
    dataset = HarmonicNativeDataset(
        dataset_rows,
        native_h5=None,
        metadata=np.empty((len(dataset_rows), 0), dtype=np.float32),
        cache_size=0,
        profile=FULLPOOL_SSL_PROFILE,
    )
    try:
        row_reports = _audit_dataset(dataset, selected)
    finally:
        dataset.close()

    if (
        len(row_reports) != len(selected)
        or len({row["ssl_observation_id"] for row in row_reports})
        != len(selected)
        or observation_identity_sha256(pd.DataFrame(row_reports))
        != _identity(selected)
    ):
        raise RuntimeError("canary model-facing audit lost exact coverage")
    failed_rows = [row for row in row_reports if row["passed"] is not True]
    n_failures = sum(int(row["n_failures"]) for row in row_reports)
    target_report = [
        row for row in row_reports if int(row["tic"]) == CANARY_TARGET_TIC
    ]
    if len(target_report) != 1:
        raise RuntimeError("canary model-facing audit lost TIC 704538814")

    for binding in bindings.values():
        _assert_unchanged(binding)
    report = {
        "schema_version": CANARY_SCHEMA,
        "stage": "pre_production_native_v3r1_model_facing_canary",
        "code_revision": expected_code_revision,
        "sector": CANARY_SECTOR,
        "shard_index": CANARY_SHARD_INDEX,
        "n_shards": CANARY_N_SHARDS,
        "required_tic": CANARY_TARGET_TIC,
        "required_tic_passed": target_report[0]["passed"] is True,
        "passed": not failed_rows and n_failures == 0,
        "n_failures": n_failures,
        "counts": {
            "source_observations": int(len(source_rows)),
            "eligible_observations": int(len(eligible_rows)),
            "excluded_observations": int(len(excluded_rows)),
            "audited_observations": int(len(row_reports)),
            "passed_observations": int(len(row_reports) - len(failed_rows)),
            "failed_observations": int(len(failed_rows)),
        },
        "identities": {
            "full_observations_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
            "eligible_observations_sha256": (
                PRODUCTION_ELIGIBLE_IDENTITY_SHA256
            ),
            "source_shard_observations_sha256": _identity(source_rows),
            "eligible_shard_observations_sha256": _identity(eligible_rows),
            "excluded_shard_observations_sha256": _identity(excluded_rows),
        },
        "native_contract_version": FULL_POOL_NATIVE_CONTRACT_VERSION,
        "native_release_binding": FULL_POOL_NATIVE_RELEASE_BINDING,
        "builder_contract_version": FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION,
        "detrend_contract_version": FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION,
        "detrend_config_sha256": FULL_POOL_NATIVE_DETREND_CONFIG_SHA256,
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
        "rank_warning_ledger_sha256": (
            FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256
        ),
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "numeric_envelope_contract": (
            FULL_POOL_NUMERIC_ENVELOPE_V1.contract_version
        ),
        "numeric_envelope_sha256": MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
        "native_h5": bindings["native_h5"],
        "authority_bindings": bindings,
        "native_verification": native_verification,
        "rows": row_reports,
        "action": "audit_only_no_clip_no_exclusion",
    }
    payload = _json_bytes(report)
    _publish_new(report_out, payload)
    digest = hashlib.sha256(payload).hexdigest()
    _publish_new(
        Path(str(report_out) + ".sha256"),
        f"{digest}  {Path(report_out).name}\n".encode("ascii"),
    )
    if report["passed"] is not True:
        raise RuntimeError(
            "S60 shard-4 model-facing canary failed: "
            f"failed_rows={len(failed_rows)} n_failures={n_failures}"
        )
    return report


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-pool", type=Path, required=True)
    parser.add_argument("--frozen-pool-summary", type=Path, required=True)
    parser.add_argument("--bls-peaks", type=Path, required=True)
    parser.add_argument("--bls-summary", type=Path, required=True)
    parser.add_argument("--native-h5", type=Path, required=True)
    parser.add_argument("--native-summary", type=Path, required=True)
    parser.add_argument("--native-sha-sidecar", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    parser.add_argument("--expected-code-revision", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    report = audit_canary(
        pool_path=args.frozen_pool,
        pool_summary_path=args.frozen_pool_summary,
        bls_path=args.bls_peaks,
        bls_summary_path=args.bls_summary,
        native_h5_path=args.native_h5,
        native_summary_path=args.native_summary,
        native_sha_sidecar_path=args.native_sha_sidecar,
        report_out=args.report_out,
        expected_code_revision=args.expected_code_revision,
    )
    print(
        json.dumps(
            {
                "event": "teacher_ssl_fullpool_v3r1_canary_numeric_complete",
                "counts": report["counts"],
                "report": str(args.report_out.expanduser().resolve()),
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
