"""Prospective, hash-bound Teacher-v3 inference for real Sector 63 targets.

This module is deliberately separate from :mod:`harmonic_inference`: the
historical Teacher-v1 scorer retains its original release contract, while the
prospective S63 run is bound to the frozen seven-sector Teacher-v3 release.
The neural-network architecture and low-level inference helpers are reused;
no second model implementation lives here.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    MODEL_VERSION,
    MORPHOLOGY_CLASSES,
    HarmonicModelConfig,
    build_harmonic_cnn,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    MetadataNormalization,
    build_metadata_matrix,
    collate_native_samples,
)
from twirl.vetting.harmonic_inference import (
    HARMONIC_DISPLAY_LABELS,
    HARMONIC_FACTORS,
    _normalization,
    _softmax,
    _stage_table,
    _to_device,
    prepare_inference_rows,
    rank_planet_enrichment,
)
from twirl.vetting.harmonic_inputs import (
    RAW_PAIR_CONTRACT_VERSION,
    native_group_path,
    verify_raw_pair_contract,
)
from twirl.vetting.teacher_v3 import TEACHER_V3_RUN_NAME
from twirl.vetting.teacher_v3_prospective import (
    PROSPECTIVE_PLAN_SCHEMA,
    load_selection_policy,
)
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_CALIBRATION_SCHEMA,
    TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
    TEACHER_V3_CHECKPOINT_NAMESPACE,
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
    verify_teacher_v3_checkpoint_manifest,
)


S63 = 63
S63_NATIVE_PERIODOGRAM_N = 4096
S63_SCORE_SCHEMA = "twirl_teacher_v3_s63_candidate_scores_v1"
S63_RANKING_SCHEMA = "twirl_teacher_v3_s63_planet_ranking_v1"
S63_SUMMARY_SCHEMA = "twirl_teacher_v3_s63_scoring_summary_v1"
S63_SCORE_POLICY = "prospective_human_review_enrichment_not_object_confirmation"
S63_LAUNCH_CONTRACT = "twirl_teacher_v3_s63_prospective_launch_v1"
S63_AUTHORIZATION_CONTRACT = "s63_teacher_v3_prospective_v1"

_RELEASE_SUMMARY_SCHEMA = "twirl_teacher_v3_training_summary_v1"
_PRETEST_FREEZE_SCHEMA = "twirl_teacher_v3_model_freeze_v1"
_SELECTION_POLICY = "fixed_before_test"
_GIT_SHA_PATTERN = re.compile(r"[0-9a-f]{40}")
_PROVENANCE_FIELDS: tuple[str, ...] = (
    "checkpoint_namespace",
    "input_contract_version",
    "native_h5_sha256",
    "native_manifest_sha256",
    "training_table_sha256",
    "native_registry_sha256",
    "native_registry_summary_sha256",
)


def file_sha256(path: Path) -> str:
    """Return the SHA-256 of one file without loading it into memory."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _require_sha256(value: str, *, name: str) -> str:
    digest = str(value)
    if len(digest) != 64 or any(
        character not in "0123456789abcdef" for character in digest
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return digest


def _read_json_with_hash(
    path: Path,
    *,
    expected_sha256: str,
    artifact: str,
) -> tuple[dict[str, Any], str]:
    path = Path(path).resolve()
    expected = _require_sha256(expected_sha256, name=f"{artifact} SHA-256")
    before = file_sha256(path)
    if before != expected:
        raise RuntimeError(
            f"trusted {artifact} SHA-256 mismatch: {path}; "
            f"observed={before}, expected={expected}"
        )
    try:
        payload = json.loads(path.read_text())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid {artifact} JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{artifact} must be a JSON mapping")
    if file_sha256(path) != before:
        raise RuntimeError(f"{artifact} changed while it was read: {path}")
    return payload, before


def _expect_fields(
    payload: Mapping[str, Any],
    expected: Mapping[str, Any],
    *,
    artifact: str,
) -> None:
    failures = [
        f"{name}={payload.get(name)!r}, expected {value!r}"
        for name, value in expected.items()
        if payload.get(name) != value
    ]
    if failures:
        raise ValueError(f"{artifact} contract failed: " + "; ".join(failures))


def _release_provenance(manifest: Mapping[str, Any]) -> dict[str, str]:
    provenance = {
        name: str(manifest.get(name, "")) for name in _PROVENANCE_FIELDS
    }
    if provenance["checkpoint_namespace"] != TEACHER_V3_CHECKPOINT_NAMESPACE:
        raise ValueError("selected manifest has the wrong checkpoint namespace")
    if provenance["input_contract_version"] != RAW_PAIR_CONTRACT_VERSION:
        raise ValueError("selected manifest has the wrong native input contract")
    for name in _PROVENANCE_FIELDS[2:]:
        _require_sha256(provenance[name], name=f"selected manifest {name}")
    return provenance


def verify_teacher_v3_release(
    *,
    release_summary_path: Path,
    expected_release_summary_sha256: str,
    pretest_freeze_path: Path,
    expected_pretest_freeze_sha256: str,
    checkpoint_manifest_path: Path,
    expected_checkpoint_manifest_sha256: str,
    calibration_path: Path,
    expected_calibration_sha256: str,
) -> dict[str, Any]:
    """Verify the exact frozen Teacher-v3 primary ensemble and release chain.

    All four document hashes are explicit trust anchors supplied by the
    caller.  The selected manifest then binds the five checkpoint paths and
    hashes, which are independently rehashed and inspected by the canonical
    Teacher-v3 verifier.
    """

    paths = {
        "release_summary": Path(release_summary_path).resolve(),
        "pretest_freeze": Path(pretest_freeze_path).resolve(),
        "checkpoint_manifest": Path(checkpoint_manifest_path).resolve(),
        "calibration": Path(calibration_path).resolve(),
    }
    expected_hashes = {
        "release_summary": expected_release_summary_sha256,
        "pretest_freeze": expected_pretest_freeze_sha256,
        "checkpoint_manifest": expected_checkpoint_manifest_sha256,
        "calibration": expected_calibration_sha256,
    }
    documents: dict[str, dict[str, Any]] = {}
    observed_hashes: dict[str, str] = {}
    for name, path in paths.items():
        documents[name], observed_hashes[name] = _read_json_with_hash(
            path,
            expected_sha256=expected_hashes[name],
            artifact=name.replace("_", " "),
        )

    summary = documents["release_summary"]
    freeze = documents["pretest_freeze"]
    manifest = documents["checkpoint_manifest"]
    calibration = documents["calibration"]
    _expect_fields(
        summary,
        {
            "schema_version": _RELEASE_SUMMARY_SCHEMA,
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
            "profile_selection_policy": _SELECTION_POLICY,
            "student_training_blocked": True,
            "automatic_production_promotion": False,
        },
        artifact="Teacher-v3 release summary",
    )
    _expect_fields(
        freeze,
        {
            "schema_version": _PRETEST_FREEZE_SCHEMA,
            "run_id": TEACHER_V3_RUN_ID,
            "model_version": MODEL_VERSION,
            "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
            "profile_selection_policy": _SELECTION_POLICY,
            "test_rows_used_for_selection_or_calibration": False,
        },
        artifact="Teacher-v3 pretest freeze",
    )
    _expect_fields(
        manifest,
        {
            "schema_version": TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "selected_profile": TEACHER_V3_PRIMARY_PROFILE,
            "selection_policy": _SELECTION_POLICY,
        },
        artifact="Teacher-v3 selected checkpoint manifest",
    )
    _expect_fields(
        calibration,
        {
            "schema_version": TEACHER_V3_CALIBRATION_SCHEMA,
            "run_id": TEACHER_V3_RUN_ID,
            "profile": TEACHER_V3_PRIMARY_PROFILE,
            "scope": "concatenated_five_fold_development_oof_logits",
        },
        artifact="Teacher-v3 pooled calibration",
    )

    manifest_sha256 = observed_hashes["checkpoint_manifest"]
    freeze_sha256 = observed_hashes["pretest_freeze"]
    calibration_sha256 = observed_hashes["calibration"]
    if str(summary.get("selected_checkpoint_manifest_sha256", "")) != (
        manifest_sha256
    ):
        raise ValueError("release summary does not bind the selected manifest")
    if str(summary.get("pretest_model_freeze_sha256", "")) != freeze_sha256:
        raise ValueError("release summary does not bind the pretest freeze")
    if str(freeze.get("primary_checkpoint_manifest_sha256", "")) != (
        manifest_sha256
    ):
        raise ValueError("pretest freeze does not bind the selected manifest")
    primary_calibrations = [
        value
        for value in summary.get("calibration", [])
        if isinstance(value, Mapping)
        and value.get("profile") == TEACHER_V3_PRIMARY_PROFILE
    ]
    if len(primary_calibrations) != 1 or str(
        primary_calibrations[0].get("calibration_sha256", "")
    ) != calibration_sha256:
        raise ValueError("release summary does not bind the primary calibration")

    provenance = _release_provenance(manifest)
    for artifact_name, payload in (
        ("release summary", summary),
        ("pretest freeze", freeze),
        ("pooled calibration", calibration),
    ):
        failures = [
            f"{name}={payload.get(name)!r}, expected {expected!r}"
            for name, expected in provenance.items()
            if str(payload.get(name, "")) != expected
        ]
        if failures:
            raise ValueError(
                f"Teacher-v3 {artifact_name} provenance differs from the "
                "selected manifest: " + "; ".join(failures)
            )

    try:
        temperature = float(calibration["temperature"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Teacher-v3 calibration temperature is invalid") from exc
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("Teacher-v3 calibration temperature must be positive")

    manifest_audit = verify_teacher_v3_checkpoint_manifest(
        manifest_path=paths["checkpoint_manifest"],
        expected_manifest_sha256=manifest_sha256,
        expected_profile=TEACHER_V3_PRIMARY_PROFILE,
        expected_selection_policy=_SELECTION_POLICY,
        expected_provenance=provenance,
        calibration_path=paths["calibration"],
        expected_calibration_sha256=calibration_sha256,
    )
    if int(manifest_audit.get("n_verified_checkpoints", -1)) != 5:
        raise RuntimeError("Teacher-v3 verifier did not prove five checkpoints")

    records = manifest.get("checkpoints")
    if not isinstance(records, list) or len(records) != 5:
        raise ValueError("selected manifest must contain exactly five checkpoints")
    checkpoint_paths: list[Path] = []
    checkpoint_sha256_by_fold: dict[str, str] = {}
    for expected_fold, record in enumerate(
        sorted(records, key=lambda value: int(value.get("fold", -1)))
    ):
        if not isinstance(record, Mapping):
            raise ValueError("selected manifest contains a malformed checkpoint")
        fold = int(record.get("fold", -1))
        if fold != expected_fold:
            raise ValueError("selected manifest checkpoint folds must be 0..4")
        relative = (
            Path(TEACHER_V3_PRIMARY_PROFILE)
            / f"fold_{fold}"
            / "teacher.pt"
        )
        if str(record.get("path", "")) != relative.as_posix():
            raise ValueError(f"selected manifest fold {fold} path is not canonical")
        digest = _require_sha256(
            str(record.get("sha256", "")),
            name=f"selected manifest fold {fold} SHA-256",
        )
        checkpoint_path = paths["checkpoint_manifest"].parent / relative
        if file_sha256(checkpoint_path) != digest:
            raise RuntimeError(f"Teacher-v3 fold {fold} checkpoint hash changed")
        checkpoint_paths.append(checkpoint_path)
        checkpoint_sha256_by_fold[str(fold)] = digest
    if manifest_audit.get("checkpoint_sha256_by_fold") != (
        checkpoint_sha256_by_fold
    ):
        raise RuntimeError("Teacher-v3 verifier checkpoint hashes disagree")

    for name, path in paths.items():
        if file_sha256(path) != observed_hashes[name]:
            raise RuntimeError(f"Teacher-v3 {name} changed during verification")
    return {
        "schema_version": "twirl_teacher_v3_verified_release_v1",
        "release_documents": {
            name: {"path": str(paths[name]), "sha256": observed_hashes[name]}
            for name in paths
        },
        "input_provenance": provenance,
        "calibration_temperature": temperature,
        "checkpoint_paths": [str(path) for path in checkpoint_paths],
        "checkpoint_sha256_by_fold": checkpoint_sha256_by_fold,
        "manifest_audit": manifest_audit,
    }


def _verify_launch_artifact(
    *,
    artifacts: Mapping[str, Any],
    name: str,
    actual_path: Path,
    expected_sha256: str,
) -> dict[str, Any]:
    record = artifacts.get(name)
    if not isinstance(record, Mapping):
        raise ValueError(f"S63 launch manifest lacks artifacts.{name}")
    actual_path = Path(actual_path).resolve()
    expected = _require_sha256(
        expected_sha256, name=f"S63 launch artifacts.{name} SHA-256"
    )
    if str(record.get("sha256", "")) != expected:
        raise ValueError(f"S63 launch artifacts.{name} hash binding differs")
    try:
        declared_size = int(record["size_bytes"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            f"S63 launch artifacts.{name} has invalid size_bytes"
        ) from exc
    if declared_size != int(actual_path.stat().st_size):
        raise ValueError(f"S63 launch artifacts.{name} size binding differs")
    declared_text = str(record.get("path", "")).strip()
    if not declared_text:
        raise ValueError(f"S63 launch artifacts.{name} has no path")
    declared_path = Path(declared_text).expanduser()
    canonical_match = declared_path.resolve(strict=False) == actual_path
    # Cross-node staging may relocate a hash-identical artifact.  If the
    # canonical path from the manifest exists on this host, however, silently
    # substituting another path is forbidden.
    if declared_path.exists() and not canonical_match:
        raise ValueError(
            f"S63 launch artifacts.{name} canonical path differs from input"
        )
    return {
        "declared_path": declared_text,
        "actual_path": str(actual_path),
        "sha256": expected,
        "size_bytes": declared_size,
        "canonical_path_match": bool(canonical_match),
        "relocated_by_hash": bool(not canonical_match),
    }


def verify_execution_checkout(
    *,
    repo: Path,
    expected_git_sha: str,
    launch_git: Mapping[str, Any],
) -> dict[str, Any]:
    """Bind execution to the exact clean Git checkout sealed at launch."""

    try:
        repo = Path(repo).expanduser().resolve(strict=True)
    except (FileNotFoundError, OSError) as exc:
        raise ValueError(f"execution repository does not exist: {repo}") from exc
    if not repo.is_dir():
        raise ValueError(f"execution repository is not a directory: {repo}")
    expected = str(expected_git_sha).strip()
    if _GIT_SHA_PATTERN.fullmatch(expected) is None:
        raise ValueError(
            "operator expected Git SHA must be a full lowercase 40-character SHA"
        )
    if not isinstance(launch_git, Mapping):
        raise ValueError("S63 launch manifest lacks a git mapping")
    if launch_git.get("checkout_clean") is not True:
        raise ValueError("S63 launch git.checkout_clean must be true")
    if launch_git.get("untracked_files_checked") is not True:
        raise ValueError("S63 launch git.untracked_files_checked must be true")
    launch_sha = str(launch_git.get("sha", "")).strip()
    if _GIT_SHA_PATTERN.fullmatch(launch_sha) is None:
        raise ValueError(
            "S63 launch git.sha must be a full lowercase 40-character SHA"
        )

    def run_git(*arguments: str) -> str:
        try:
            completed = subprocess.run(
                ["git", *arguments],
                cwd=repo,
                check=True,
                capture_output=True,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            detail = getattr(exc, "stderr", "") or str(exc)
            raise ValueError(
                f"could not audit execution repository with git: {detail.strip()}"
            ) from exc
        return completed.stdout

    observed = run_git("rev-parse", "HEAD").strip()
    if _GIT_SHA_PATTERN.fullmatch(observed) is None:
        raise ValueError("git rev-parse HEAD did not return a full Git SHA")
    status = run_git("status", "--porcelain", "--untracked-files=all")
    if status:
        raise ValueError(
            "execution checkout is not fully clean, including untracked files"
        )
    if not observed == expected == launch_sha:
        raise ValueError(
            "execution checkout Git SHA binding differs: "
            f"observed={observed}, operator_expected={expected}, "
            f"launch={launch_sha}"
        )
    return {
        "schema_version": "twirl_execution_checkout_audit_v1",
        "passed": True,
        "repo": str(repo),
        "launch_declared_repo": str(launch_git.get("repo", "")),
        "observed_git_sha": observed,
        "operator_expected_git_sha": expected,
        "launch_git_sha": launch_sha,
        "checkout_clean": True,
        "untracked_files_checked": True,
        "status_command": "git status --porcelain --untracked-files=all",
    }


def verify_s63_launch_authorization(
    *,
    repo: Path,
    expected_git_sha: str,
    launch_manifest_path: Path,
    expected_launch_manifest_sha256: str,
    preregistered_contract_path: Path,
    selection_policy_path: Path,
    release_verification: Mapping[str, Any],
    candidates_path: Path,
    expected_candidates_sha256: str,
    reserved_tics_path: Path,
    expected_reserved_tics_sha256: str,
    native_h5: Path,
    expected_native_h5_sha256: str,
) -> dict[str, Any]:
    """Verify the immutable preprocessing gate and its preregistered plan."""

    launch_manifest_path = Path(launch_manifest_path).resolve()
    preregistered_contract_path = Path(preregistered_contract_path).resolve()
    selection_policy_path = Path(selection_policy_path).resolve()
    launch, launch_sha256 = _read_json_with_hash(
        launch_manifest_path,
        expected_sha256=expected_launch_manifest_sha256,
        artifact="S63 launch manifest",
    )
    _expect_fields(
        launch,
        {
            "contract_version": S63_LAUNCH_CONTRACT,
            "authorization_contract": S63_AUTHORIZATION_CONTRACT,
            "sector": S63,
            "orbits": [133, 134],
            "status": "preprocessing_complete_scoring_not_started",
            "score_or_queue_read": False,
            "passed": True,
        },
        artifact="S63 prospective launch manifest",
    )
    execution_checkout = verify_execution_checkout(
        repo=repo,
        expected_git_sha=expected_git_sha,
        launch_git=launch.get("git"),
    )
    artifacts = launch.get("artifacts")
    if not isinstance(artifacts, Mapping):
        raise ValueError("S63 launch manifest lacks an artifacts mapping")
    artifact_audit = {
        "candidates": _verify_launch_artifact(
            artifacts=artifacts,
            name="candidates",
            actual_path=candidates_path,
            expected_sha256=expected_candidates_sha256,
        ),
        "reserved_tics": _verify_launch_artifact(
            artifacts=artifacts,
            name="reserved_tics",
            actual_path=reserved_tics_path,
            expected_sha256=expected_reserved_tics_sha256,
        ),
        "native_h5": _verify_launch_artifact(
            artifacts=artifacts,
            name="native_h5",
            actual_path=native_h5,
            expected_sha256=expected_native_h5_sha256,
        ),
    }
    prereg_record = artifacts.get("preregistered_contract")
    if not isinstance(prereg_record, Mapping):
        raise ValueError("S63 launch manifest lacks its preregistered contract")
    prereg_sha256 = _require_sha256(
        str(prereg_record.get("sha256", "")),
        name="S63 preregistered contract SHA-256",
    )
    preregistered, observed_prereg_sha256 = _read_json_with_hash(
        preregistered_contract_path,
        expected_sha256=prereg_sha256,
        artifact="S63 preregistered contract",
    )
    artifact_audit["preregistered_contract"] = _verify_launch_artifact(
        artifacts=artifacts,
        name="preregistered_contract",
        actual_path=preregistered_contract_path,
        expected_sha256=observed_prereg_sha256,
    )
    selection_record = artifacts.get("selection_policy")
    if not isinstance(selection_record, Mapping):
        raise ValueError("S63 launch manifest lacks artifacts.selection_policy")
    selection_sha256 = _require_sha256(
        str(selection_record.get("sha256", "")),
        name="S63 selection policy SHA-256",
    )
    selection_policy, observed_selection_sha256 = load_selection_policy(
        selection_policy_path,
        expected_sha256=selection_sha256,
    )
    artifact_audit["selection_policy"] = _verify_launch_artifact(
        artifacts=artifacts,
        name="selection_policy",
        actual_path=selection_policy_path,
        expected_sha256=observed_selection_sha256,
    )
    _expect_fields(
        preregistered,
        {
            "schema_version": PROSPECTIVE_PLAN_SCHEMA,
            "status": "plan_frozen_launch_manifest_pending",
            "sector": S63,
        },
        artifact="S63 preregistered prospective plan",
    )

    release_documents = release_verification.get("release_documents")
    checkpoint_hashes = release_verification.get("checkpoint_sha256_by_fold")
    if not isinstance(release_documents, Mapping) or not isinstance(
        checkpoint_hashes, Mapping
    ):
        raise ValueError("verified Teacher-v3 release lacks immutable hashes")
    expected_documents = {
        "summary_sha256": "release_summary",
        "pretest_model_freeze_sha256": "pretest_freeze",
        "selected_checkpoint_manifest_sha256": "checkpoint_manifest",
        "pooled_oof_calibration_sha256": "calibration",
    }
    teacher = preregistered.get("frozen_teacher_v3")
    if not isinstance(teacher, Mapping):
        raise ValueError("S63 preregistered plan lacks frozen_teacher_v3")
    _expect_fields(
        teacher,
        {
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "selected_profile": TEACHER_V3_PRIMARY_PROFILE,
            "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
            "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
            "selection_policy": _SELECTION_POLICY,
            "temperature_calibration_scope": (
                "concatenated_five_fold_development_oof_logits"
            ),
            "automatic_promotion": False,
            "student_training_authorized": False,
        },
        artifact="S63 preregistered Teacher-v3 identity",
    )
    document_bindings: dict[str, str] = {}
    for plan_field, release_name in expected_documents.items():
        record = release_documents.get(release_name)
        if not isinstance(record, Mapping):
            raise ValueError(f"verified Teacher-v3 release lacks {release_name}")
        expected = _require_sha256(
            str(record.get("sha256", "")),
            name=f"verified Teacher-v3 {release_name} SHA-256",
        )
        observed = str(teacher.get(plan_field, ""))
        if observed != expected:
            raise ValueError(
                f"S63 preregistered {plan_field} differs from verified release"
            )
        document_bindings[plan_field] = expected
    try:
        plan_temperature = float(teacher["temperature"])
        release_temperature = float(
            release_verification["calibration_temperature"]
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("S63 preregistered Teacher-v3 temperature is invalid") from exc
    if not np.isclose(plan_temperature, release_temperature):
        raise ValueError(
            "S63 preregistered temperature differs from verified calibration"
        )

    records = teacher.get("checkpoints")
    if not isinstance(records, list) or len(records) != 5:
        raise ValueError("S63 preregistered plan must bind five checkpoints")
    observed_checkpoint_hashes: dict[str, str] = {}
    for expected_fold, record in enumerate(
        sorted(records, key=lambda value: int(value.get("fold", -1)))
    ):
        if not isinstance(record, Mapping):
            raise ValueError("S63 preregistered checkpoint record is malformed")
        fold = int(record.get("fold", -1))
        if fold != expected_fold:
            raise ValueError("S63 preregistered checkpoint folds must be 0..4")
        digest = _require_sha256(
            str(record.get("sha256", "")),
            name=f"S63 preregistered fold-{fold} SHA-256",
        )
        if digest != str(checkpoint_hashes.get(str(fold), "")):
            raise ValueError(
                f"S63 preregistered fold-{fold} differs from verified release"
            )
        observed_checkpoint_hashes[str(fold)] = digest

    reservation = preregistered.get("s63_identity_reservation")
    if not isinstance(reservation, Mapping):
        raise ValueError("S63 preregistered plan lacks identity reservation")
    _expect_fields(
        reservation,
        {
            "identity_only": True,
            "light_curves_opened_when_reserved": False,
            "reserved_tics_sha256": str(expected_reserved_tics_sha256),
        },
        artifact="S63 preregistered identity reservation",
    )

    checks = launch.get("checks")
    if not isinstance(checks, Mapping):
        raise ValueError("S63 launch manifest lacks preprocessing checks")
    required_checks = (
        "accepted_stage1",
        "selection_policy",
        "cadence",
        "compact",
        "inventories",
        "model_ready_derivation",
        "bls",
        "candidates",
        "raw_source",
        "native",
    )
    for name in required_checks:
        check = checks.get(name)
        if not isinstance(check, Mapping) or check.get("passed") is not True:
            raise ValueError(f"S63 launch check {name} did not pass")

    if file_sha256(launch_manifest_path) != launch_sha256:
        raise RuntimeError("S63 launch manifest changed during authorization")
    if file_sha256(preregistered_contract_path) != observed_prereg_sha256:
        raise RuntimeError("S63 preregistered contract changed during authorization")
    if file_sha256(selection_policy_path) != observed_selection_sha256:
        raise RuntimeError("S63 selection policy changed during authorization")
    return {
        "schema_version": "twirl_teacher_v3_s63_launch_authorization_v1",
        "passed": True,
        "score_or_queue_read_before_authorization": False,
        "launch_manifest": {
            "path": str(launch_manifest_path),
            "sha256": launch_sha256,
        },
        "preregistered_contract": {
            "path": str(preregistered_contract_path),
            "sha256": observed_prereg_sha256,
        },
        "selection_policy": {
            "path": str(selection_policy_path),
            "sha256": observed_selection_sha256,
            "schema_version": str(selection_policy["schema_version"]),
        },
        "execution_checkout": execution_checkout,
        "artifact_bindings": artifact_audit,
        "teacher_release_document_bindings": document_bindings,
        "teacher_checkpoint_sha256_by_fold": observed_checkpoint_hashes,
    }


def _as_bool(values: pd.Series, *, name: str) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    normalized = values.fillna("").astype(str).str.strip().str.lower()
    valid = normalized.isin(
        {"1", "1.0", "true", "t", "yes", "y", "0", "0.0", "false", "f", "no", "n", ""}
    )
    if not valid.all():
        raise ValueError(f"candidate {name} contains invalid Boolean values")
    return normalized.isin({"1", "1.0", "true", "t", "yes", "y"})


def _read_candidate_table(path: Path) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"unsupported candidate table format: {path}")


def verify_s63_scoring_inputs(
    *,
    candidates_path: Path,
    expected_candidates_sha256: str,
    reserved_tics_path: Path,
    expected_reserved_tics_sha256: str,
    native_h5: Path,
    expected_native_h5_sha256: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Prove that one native file exactly represents real S63 candidates."""

    import h5py

    candidates_path = Path(candidates_path).resolve()
    reserved_tics_path = Path(reserved_tics_path).resolve()
    native_h5 = Path(native_h5).resolve()
    expected_hashes = {
        "candidate_table": _require_sha256(
            expected_candidates_sha256, name="candidate table SHA-256"
        ),
        "reserved_tics": _require_sha256(
            expected_reserved_tics_sha256, name="S63 reservation SHA-256"
        ),
        "native_h5": _require_sha256(
            expected_native_h5_sha256, name="native HDF5 SHA-256"
        ),
    }
    paths = {
        "candidate_table": candidates_path,
        "reserved_tics": reserved_tics_path,
        "native_h5": native_h5,
    }
    for name, path in paths.items():
        observed = file_sha256(path)
        if observed != expected_hashes[name]:
            raise RuntimeError(
                f"trusted {name} SHA-256 mismatch: observed={observed}, "
                f"expected={expected_hashes[name]}"
            )

    candidates = _read_candidate_table(candidates_path)
    required = {
        "review_id",
        "sector",
        "tic",
        "period_d",
        "t0_bjd",
        "duration_min",
        "native_input_include",
    }
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise KeyError(f"S63 candidate table lacks columns: {missing}")
    if candidates.empty:
        raise ValueError("S63 candidate table is empty")
    rows = prepare_inference_rows(candidates, allow_injections=False)
    sectors = set(
        pd.to_numeric(rows["sector"], errors="raise").astype(int).tolist()
    )
    if sectors != {S63}:
        raise ValueError(f"prospective Teacher-v3 scoring is S63-only; got {sectors}")
    review_ids = rows["review_id"].astype(str)
    if review_ids.duplicated().any():
        raise ValueError("S63 candidate review_id values must be unique")
    if not _as_bool(
        candidates["native_input_include"], name="native_input_include"
    ).all():
        raise ValueError("every S63 scoring row must be native_input_include=true")
    if "is_injected_row" in candidates and _as_bool(
        candidates["is_injected_row"], name="is_injected_row"
    ).any():
        raise ValueError("S63 production scoring accepts real candidates only")
    if "source_kind" in candidates and candidates["source_kind"].fillna("").astype(
        str
    ).str.lower().str.contains("inject").any():
        raise ValueError("S63 production scoring contains an injected source_kind")

    candidate_tic_values = pd.to_numeric(
        rows["tic"], errors="raise"
    ).astype(np.int64)
    if candidate_tic_values.duplicated().any():
        duplicated = sorted(
            candidate_tic_values.loc[candidate_tic_values.duplicated(False)]
            .astype(int)
            .unique()
            .tolist()
        )
        raise ValueError(
            "S63 prospective scoring requires exactly one candidate row per TIC; "
            f"duplicates={duplicated[:5]}"
        )
    expected_group_by_tic = {
        int(row["tic"]): native_group_path(row)
        for row in rows.to_dict("records")
    }
    if any(
        not value.startswith("targets/")
        for value in expected_group_by_tic.values()
    ):
        raise ValueError("S63 candidates must map only to real target groups")
    if "native_group_path" in candidates:
        observed_groups = rows["native_group_path"].astype(str)
        derived_groups = pd.Series(
            [
                expected_group_by_tic[int(tic)]
                for tic in pd.to_numeric(rows["tic"], errors="raise")
            ],
            index=rows.index,
        )
        if not observed_groups.equals(derived_groups):
            raise ValueError("candidate native_group_path differs from its TIC identity")

    try:
        reserved_values = [
            int(value)
            for value in reserved_tics_path.read_text().splitlines()
            if value.strip()
        ]
    except (OSError, ValueError) as exc:
        raise ValueError(f"invalid S63 reserved-TIC inventory: {exc}") from exc
    if (
        not reserved_values
        or any(value <= 0 for value in reserved_values)
        or reserved_values != sorted(set(reserved_values))
    ):
        raise ValueError("S63 reserved-TIC inventory must be sorted, unique, and positive")
    candidate_tics = set(expected_group_by_tic)
    unreserved = sorted(candidate_tics - set(reserved_values))
    if unreserved:
        raise ValueError(
            f"S63 candidates include {len(unreserved)} unreserved hosts; "
            f"first={unreserved[:5]}"
        )

    native_contract = verify_raw_pair_contract(
        native_h5,
        require_errors=True,
        require_periodograms=True,
    )
    if not native_contract["passed"]:
        raise ValueError(
            "S63 native HDF5 contract failed: "
            + "; ".join(native_contract["failures"][:10])
        )
    expected_groups = set(expected_group_by_tic.values())
    with h5py.File(native_h5, "r") as h5:
        if str(h5.attrs.get("contract_version", "")) != RAW_PAIR_CONTRACT_VERSION:
            raise ValueError("S63 native HDF5 has the wrong contract_version")
        if str(h5.attrs.get("training_table_sha256", "")) != expected_hashes[
            "candidate_table"
        ]:
            raise ValueError("S63 native HDF5 is not bound to the candidate table")
        try:
            periodogram_n = int(h5.attrs.get("periodogram_n", -1))
        except (TypeError, ValueError):
            periodogram_n = -1
        if periodogram_n != S63_NATIVE_PERIODOGRAM_N:
            raise ValueError(
                f"S63 native periodogram_n={periodogram_n}, "
                f"expected {S63_NATIVE_PERIODOGRAM_N}"
            )
        injected_groups = list(h5["injections"].keys()) if "injections" in h5 else []
        if injected_groups:
            raise ValueError("S63 native HDF5 contains injection groups")
        observed_groups = {
            f"targets/{name}" for name in h5.get("targets", {})
        }
        missing_groups = sorted(expected_groups - observed_groups)
        extra_groups = sorted(observed_groups - expected_groups)
        if missing_groups or extra_groups:
            raise ValueError(
                "S63 native/candidate identity mismatch: "
                f"missing={missing_groups[:5]}, extra={extra_groups[:5]}"
            )
        for tic, group_path in expected_group_by_tic.items():
            group = h5[group_path]
            try:
                group_tic = int(group.attrs["tic"])
                group_sector = int(group.attrs["sector"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"{group_path} lacks exact TIC/sector identity attrs"
                ) from exc
            if group_tic != tic or group_sector != S63:
                raise ValueError(
                    f"{group_path} identity is tic={group_tic}, sector={group_sector}"
                )

    for name, path in paths.items():
        if file_sha256(path) != expected_hashes[name]:
            raise RuntimeError(f"{name} changed during S63 input verification")
    return candidates, {
        "schema_version": "twirl_teacher_v3_s63_input_verification_v1",
        "passed": True,
        "sector": S63,
        "real_only": True,
        "candidate_table_sha256": expected_hashes["candidate_table"],
        "reserved_tics_sha256": expected_hashes["reserved_tics"],
        "native_h5_sha256": expected_hashes["native_h5"],
        "n_candidate_rows": int(len(candidates)),
        "n_candidate_tics": int(len(candidate_tics)),
        "n_exact_native_groups": int(len(expected_groups)),
        "n_reserved_tics": int(len(reserved_values)),
        "native_contract": native_contract,
    }


def _load_verified_models(
    *,
    checkpoint_paths: Sequence[Path],
    checkpoint_sha256_by_fold: Mapping[str, str],
    input_provenance: Mapping[str, str],
    calibration_temperature: float,
    device: Any,
) -> tuple[list[Any], list[dict[str, Any]], list[Any]]:
    import torch

    if len(checkpoint_paths) != 5:
        raise ValueError("Teacher-v3 inference requires exactly five checkpoints")
    models: list[Any] = []
    checkpoints: list[dict[str, Any]] = []
    normalizations: list[Any] = []
    reference_model_config: Mapping[str, Any] | None = None
    for fold, path in enumerate(checkpoint_paths):
        path = Path(path).resolve()
        expected_hash = str(checkpoint_sha256_by_fold.get(str(fold), ""))
        if file_sha256(path) != expected_hash:
            raise RuntimeError(f"Teacher-v3 fold {fold} checkpoint hash changed")
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        expected_fields = {
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
            "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
            "profile": TEACHER_V3_PRIMARY_PROFILE,
            "fold": fold,
            "temperature_calibration_scope": "pooled_oof_development",
        }
        _expect_fields(
            checkpoint,
            expected_fields,
            artifact=f"Teacher-v3 fold-{fold} checkpoint",
        )
        for name, expected in input_provenance.items():
            if str(checkpoint.get(name, "")) != str(expected):
                raise ValueError(
                    f"Teacher-v3 fold-{fold} checkpoint {name} differs from release"
                )
        if not np.isclose(
            float(checkpoint.get("temperature", float("nan"))),
            float(calibration_temperature),
        ):
            raise ValueError(
                f"Teacher-v3 fold-{fold} temperature differs from calibration"
            )
        model_config = checkpoint.get("model_config")
        if not isinstance(model_config, Mapping):
            raise ValueError(f"Teacher-v3 fold-{fold} lacks a model_config")
        if reference_model_config is None:
            reference_model_config = dict(model_config)
        elif dict(model_config) != dict(reference_model_config):
            raise ValueError("Teacher-v3 fold checkpoints disagree on model_config")
        model = build_harmonic_cnn(
            HarmonicModelConfig(**dict(model_config)),
            profile=TEACHER_V3_PRIMARY_PROFILE,
        ).to(device)
        model.load_state_dict(checkpoint["model_state_dict"], strict=True)
        model.eval()
        if file_sha256(path) != expected_hash:
            raise RuntimeError(f"Teacher-v3 fold {fold} changed while it was loaded")
        models.append(model)
        checkpoints.append(checkpoint)
        normalizations.append(_normalization(checkpoint))
    return models, checkpoints, normalizations


def _validated_probability(values: np.ndarray, *, name: str) -> np.ndarray:
    probability = np.asarray(values, dtype=np.float64)
    if probability.ndim != 3 or probability.shape[0] != 5:
        raise RuntimeError(f"{name} ensemble probability has shape {probability.shape}")
    if not np.isfinite(probability).all():
        raise FloatingPointError(f"{name} ensemble probability is non-finite")
    if not np.allclose(probability.sum(axis=2), 1.0, atol=1.0e-6):
        raise FloatingPointError(f"{name} ensemble probability is not normalized")
    return probability


def _validate_checkpoint_metadata_inputs(
    rows: pd.DataFrame,
    normalizations: Sequence[MetadataNormalization],
) -> dict[str, Any]:
    """Require every checkpoint metadata feature to have S63 support.

    Row-level non-finite values retain the frozen training behavior and are
    imputed to the fold-local normalization center by ``build_metadata_matrix``.
    An absent or wholly non-finite feature is different: it would silently
    replace that learned feature with its training mean for every S63 target.
    """

    required_columns = tuple(
        dict.fromkeys(
            column
            for normalization in normalizations
            for column in normalization.columns
        )
    )
    missing = sorted(set(required_columns) - set(rows.columns))
    if missing:
        raise KeyError(
            "Teacher-v3 checkpoint metadata columns are absent from the "
            f"S63 candidate table: {missing}"
        )

    finite_rows_by_column: dict[str, int] = {}
    wholly_nonfinite: list[str] = []
    partially_nonfinite: list[str] = []
    for column in required_columns:
        values = pd.to_numeric(rows[column], errors="coerce").to_numpy(dtype=float)
        finite_count = int(np.isfinite(values).sum())
        finite_rows_by_column[column] = finite_count
        if finite_count == 0:
            wholly_nonfinite.append(column)
        elif finite_count < len(rows):
            partially_nonfinite.append(column)
    if wholly_nonfinite:
        raise ValueError(
            "Teacher-v3 checkpoint metadata columns are wholly non-finite "
            f"across S63: {sorted(wholly_nonfinite)}"
        )

    return {
        "schema_version": "twirl_teacher_v3_s63_metadata_input_validation_v1",
        "passed": True,
        "n_rows": int(len(rows)),
        "n_required_columns": int(len(required_columns)),
        "required_columns": list(required_columns),
        "finite_rows_by_column": finite_rows_by_column,
        "partially_nonfinite_columns": sorted(partially_nonfinite),
        "row_level_center_imputation_allowed": True,
        "whole_column_center_imputation_allowed": False,
    }


def score_teacher_v3_ensemble(
    *,
    candidates: pd.DataFrame,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
    checkpoint_sha256_by_fold: Mapping[str, str],
    input_provenance: Mapping[str, str],
    calibration_temperature: float,
    batch_size: int = 32,
    workers: int = 4,
    require_cuda: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run the frozen five-fold Teacher-v3 primary ensemble on real S63 rows."""

    import torch
    from torch.utils.data import DataLoader

    if int(batch_size) < 1:
        raise ValueError("batch_size must be positive")
    if int(workers) < 0:
        raise ValueError("workers cannot be negative")
    rows = prepare_inference_rows(candidates, allow_injections=False)
    sectors = set(pd.to_numeric(rows["sector"], errors="raise").astype(int))
    if sectors != {S63}:
        raise ValueError("Teacher-v3 prospective ensemble accepts S63 only")
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("Teacher-v3 S63 scoring requires CUDA")
    models, checkpoints, normalizations = _load_verified_models(
        checkpoint_paths=checkpoint_paths,
        checkpoint_sha256_by_fold=checkpoint_sha256_by_fold,
        input_provenance=input_provenance,
        calibration_temperature=calibration_temperature,
        device=device,
    )
    metadata_input_validation = _validate_checkpoint_metadata_inputs(
        rows,
        normalizations,
    )

    metadata_matrices: list[np.ndarray] = []
    for normalization in normalizations:
        metadata, _ = build_metadata_matrix(
            rows,
            fit_mask=np.zeros(len(rows), dtype=bool),
            normalization=normalization,
        )
        metadata_matrices.append(metadata)
    dataset = HarmonicNativeDataset(
        rows,
        native_h5=Path(native_h5),
        metadata=metadata_matrices[0],
        cache_size=0,
        profile=TEACHER_V3_PRIMARY_PROFILE,
    )
    loader = DataLoader(
        dataset,
        batch_size=int(batch_size),
        shuffle=False,
        num_workers=int(workers),
        pin_memory=device.type == "cuda",
        persistent_workers=int(workers) > 0,
        collate_fn=collate_native_samples,
    )
    row_lookup = {
        value: index for index, value in enumerate(rows["review_id"].astype(str))
    }
    morphology_members: list[list[np.ndarray]] = [[] for _ in range(5)]
    preserve_members: list[list[np.ndarray]] = [[] for _ in range(5)]
    harmonic_members: list[list[np.ndarray]] = [[] for _ in range(5)]
    output_order: list[str] = []
    with torch.no_grad():
        for batch_index, raw_batch in enumerate(loader, start=1):
            review_ids = [str(value) for value in raw_batch["review_id"]]
            indices = np.asarray([row_lookup[value] for value in review_ids], dtype=int)
            output_order.extend(review_ids)
            batch = _to_device(raw_batch, device)
            for member, (model, metadata) in enumerate(
                zip(models, metadata_matrices)
            ):
                batch["metadata"] = torch.as_tensor(
                    metadata[indices], dtype=torch.float32, device=device
                )
                with torch.autocast(
                    device_type=device.type,
                    dtype=torch.bfloat16,
                    enabled=device.type == "cuda",
                ):
                    output = model(batch)
                morphology_members[member].append(
                    _softmax(
                        output["morphology_logits"].float().cpu().numpy(),
                        temperature=float(calibration_temperature),
                    )
                )
                preserve_members[member].append(
                    _softmax(output["preserve_logits"].float().cpu().numpy())
                )
                harmonic_members[member].append(
                    _softmax(output["harmonic_logits"].float().cpu().numpy())
                )
            if batch_index % 100 == 0:
                print(
                    "[teacher-v3-s63] "
                    f"batches={batch_index:,} rows={len(output_order):,}/{len(rows):,}",
                    flush=True,
                )
    if output_order != rows["review_id"].astype(str).tolist():
        raise RuntimeError("Teacher-v3 S63 inference output ordering changed")

    morphology = _validated_probability(
        np.stack(
            [np.concatenate(values, axis=0) for values in morphology_members],
            axis=0,
        ),
        name="morphology",
    )
    preserve = _validated_probability(
        np.stack(
            [np.concatenate(values, axis=0) for values in preserve_members], axis=0
        ),
        name="preserve",
    )
    harmonic = _validated_probability(
        np.stack(
            [np.concatenate(values, axis=0) for values in harmonic_members], axis=0
        ),
        name="harmonic",
    )
    mean_morphology = morphology.mean(axis=0)
    mean_preserve = preserve.mean(axis=0)
    mean_harmonic = harmonic.mean(axis=0)
    scored = candidates.copy().reset_index(drop=True)
    scored.insert(0, "score_schema_version", S63_SCORE_SCHEMA)
    for class_index, label in enumerate(MORPHOLOGY_CLASSES):
        scored[f"p_{label}"] = mean_morphology[:, class_index]
        scored[f"std_p_{label}"] = morphology[:, :, class_index].std(axis=0)
        for member in range(5):
            scored[f"member_{member}_p_{label}"] = morphology[
                member, :, class_index
            ]
    scored["p_preserve"] = mean_preserve[:, 1]
    scored["std_p_preserve"] = preserve[:, :, 1].std(axis=0)
    for class_index, label in enumerate(HARMONIC_DISPLAY_LABELS):
        scored[f"p_harmonic_{label}"] = mean_harmonic[:, class_index]
        scored[f"std_p_harmonic_{label}"] = harmonic[:, :, class_index].std(axis=0)
    sorted_probability = np.sort(mean_morphology, axis=1)
    scored["morphology_entropy"] = -np.sum(
        mean_morphology * np.log(np.clip(mean_morphology, 1.0e-12, 1.0)), axis=1
    )
    scored["morphology_margin"] = (
        sorted_probability[:, -1] - sorted_probability[:, -2]
    )
    scored["ensemble_disagreement"] = morphology.std(axis=0).mean(axis=1)
    scored["predicted_morphology"] = np.asarray(
        MORPHOLOGY_CLASSES, dtype=object
    )[mean_morphology.argmax(axis=1)]
    factor_index = mean_harmonic.argmax(axis=1)
    scored["predicted_period_factor"] = np.asarray(HARMONIC_FACTORS)[factor_index]
    scored["predicted_period_d"] = (
        pd.to_numeric(scored["period_d"], errors="raise")
        * scored["predicted_period_factor"]
    )
    scored["model_version"] = MODEL_VERSION
    scored["model_profile"] = TEACHER_V3_PRIMARY_PROFILE
    scored["checkpoint_namespace"] = TEACHER_V3_CHECKPOINT_NAMESPACE
    scored["input_contract_version"] = RAW_PAIR_CONTRACT_VERSION
    scored["teacher_run_id"] = TEACHER_V3_RUN_ID
    scored["teacher_training_sectors"] = "56-62"
    scored["scoring_sector"] = S63
    scored["score_policy"] = S63_SCORE_POLICY

    return scored, {
        "schema_version": "twirl_teacher_v3_s63_ensemble_execution_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "profile": TEACHER_V3_PRIMARY_PROFILE,
        "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "calibration_temperature": float(calibration_temperature),
        "sector": S63,
        "n_rows": int(len(scored)),
        "n_unique_tics": int(
            pd.to_numeric(scored["tic"], errors="raise").nunique()
        ),
        "device": str(device),
        "torch_version": str(torch.__version__),
        "torch_cuda_version": str(torch.version.cuda),
        "cuda_device_name": (
            str(torch.cuda.get_device_name(0)) if device.type == "cuda" else ""
        ),
        "batch_size": int(batch_size),
        "workers": int(workers),
        "checkpoint_sha256_by_fold": dict(checkpoint_sha256_by_fold),
        "metadata_columns_by_fold": [
            list(value.columns) for value in normalizations
        ],
        "metadata_input_validation": metadata_input_validation,
        "predicted_class_counts": {
            str(key): int(value)
            for key, value in scored["predicted_morphology"]
            .value_counts()
            .sort_index()
            .items()
        },
    }


def _all_input_paths(
    *,
    release: Mapping[str, Any],
    launch_manifest_path: Path,
    preregistered_contract_path: Path,
    selection_policy_path: Path,
    candidates_path: Path,
    reserved_tics_path: Path,
    native_h5: Path,
) -> dict[str, Path]:
    paths = {
        name: Path(record["path"]).resolve()
        for name, record in release["release_documents"].items()
    }
    paths.update(
        {
            f"checkpoint_fold_{fold}": Path(path).resolve()
            for fold, path in enumerate(release["checkpoint_paths"])
        }
    )
    paths.update(
        {
            "s63_launch_manifest": Path(launch_manifest_path).resolve(),
            "s63_preregistered_contract": Path(
                preregistered_contract_path
            ).resolve(),
            "s63_selection_policy": Path(selection_policy_path).resolve(),
            "candidate_table": Path(candidates_path).resolve(),
            "reserved_tics": Path(reserved_tics_path).resolve(),
            "native_h5": Path(native_h5).resolve(),
        }
    )
    return paths


def _snapshot(
    paths: Mapping[str, Path],
    *,
    trusted_hashes: Mapping[str, str] | None = None,
) -> dict[str, dict[str, Any]]:
    """Snapshot paths, optionally reusing hashes just proven by a verifier."""

    if trusted_hashes is not None and set(trusted_hashes) != set(paths):
        raise ValueError("trusted snapshot hashes do not exactly cover input paths")
    output: dict[str, dict[str, Any]] = {}
    for name, path in paths.items():
        stat = path.stat()
        output[name] = {
            "path": str(path),
            "sha256": (
                _require_sha256(
                    trusted_hashes[name], name=f"trusted {name} SHA-256"
                )
                if trusted_hashes is not None
                else file_sha256(path)
            ),
            "size_bytes": int(stat.st_size),
            "mtime_ns": int(stat.st_mtime_ns),
            "inode": int(stat.st_ino),
        }
    return output


def _trusted_input_hashes(
    *,
    release: Mapping[str, Any],
    launch_authorization: Mapping[str, Any],
    expected_candidates_sha256: str,
    expected_reserved_tics_sha256: str,
    expected_native_h5_sha256: str,
) -> dict[str, str]:
    hashes = {
        name: str(record["sha256"])
        for name, record in release["release_documents"].items()
    }
    hashes.update(
        {
            f"checkpoint_fold_{fold}": str(
                release["checkpoint_sha256_by_fold"][str(fold)]
            )
            for fold in range(5)
        }
    )
    hashes.update(
        {
            "s63_launch_manifest": str(
                launch_authorization["launch_manifest"]["sha256"]
            ),
            "s63_preregistered_contract": str(
                launch_authorization["preregistered_contract"]["sha256"]
            ),
            "s63_selection_policy": str(
                launch_authorization["selection_policy"]["sha256"]
            ),
            "candidate_table": str(expected_candidates_sha256),
            "reserved_tics": str(expected_reserved_tics_sha256),
            "native_h5": str(expected_native_h5_sha256),
        }
    )
    return hashes


def score_teacher_v3_s63_to_disk(
    *,
    repo: Path,
    expected_git_sha: str,
    launch_manifest_path: Path,
    expected_launch_manifest_sha256: str,
    preregistered_contract_path: Path,
    selection_policy_path: Path,
    release_summary_path: Path,
    expected_release_summary_sha256: str,
    pretest_freeze_path: Path,
    expected_pretest_freeze_sha256: str,
    checkpoint_manifest_path: Path,
    expected_checkpoint_manifest_sha256: str,
    calibration_path: Path,
    expected_calibration_sha256: str,
    candidates_path: Path,
    expected_candidates_sha256: str,
    reserved_tics_path: Path,
    expected_reserved_tics_sha256: str,
    native_h5: Path,
    expected_native_h5_sha256: str,
    out_dir: Path,
    batch_size: int = 32,
    workers: int = 4,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Verify, score, and atomically publish the prospective S63 artifacts."""

    out_dir = Path(out_dir).resolve()
    if out_dir.exists():
        raise FileExistsError(
            f"Teacher-v3 S63 output directory must be a fresh path: {out_dir}"
        )
    release = verify_teacher_v3_release(
        release_summary_path=release_summary_path,
        expected_release_summary_sha256=expected_release_summary_sha256,
        pretest_freeze_path=pretest_freeze_path,
        expected_pretest_freeze_sha256=expected_pretest_freeze_sha256,
        checkpoint_manifest_path=checkpoint_manifest_path,
        expected_checkpoint_manifest_sha256=expected_checkpoint_manifest_sha256,
        calibration_path=calibration_path,
        expected_calibration_sha256=expected_calibration_sha256,
    )
    launch_authorization = verify_s63_launch_authorization(
        repo=repo,
        expected_git_sha=expected_git_sha,
        launch_manifest_path=launch_manifest_path,
        expected_launch_manifest_sha256=expected_launch_manifest_sha256,
        preregistered_contract_path=preregistered_contract_path,
        selection_policy_path=selection_policy_path,
        release_verification=release,
        candidates_path=candidates_path,
        expected_candidates_sha256=expected_candidates_sha256,
        reserved_tics_path=reserved_tics_path,
        expected_reserved_tics_sha256=expected_reserved_tics_sha256,
        native_h5=native_h5,
        expected_native_h5_sha256=expected_native_h5_sha256,
    )
    candidates, input_verification = verify_s63_scoring_inputs(
        candidates_path=candidates_path,
        expected_candidates_sha256=expected_candidates_sha256,
        reserved_tics_path=reserved_tics_path,
        expected_reserved_tics_sha256=expected_reserved_tics_sha256,
        native_h5=native_h5,
        expected_native_h5_sha256=expected_native_h5_sha256,
    )
    all_paths = _all_input_paths(
        release=release,
        launch_manifest_path=launch_manifest_path,
        preregistered_contract_path=preregistered_contract_path,
        selection_policy_path=selection_policy_path,
        candidates_path=candidates_path,
        reserved_tics_path=reserved_tics_path,
        native_h5=native_h5,
    )
    trusted_hashes = _trusted_input_hashes(
        release=release,
        launch_authorization=launch_authorization,
        expected_candidates_sha256=expected_candidates_sha256,
        expected_reserved_tics_sha256=expected_reserved_tics_sha256,
        expected_native_h5_sha256=expected_native_h5_sha256,
    )
    before = _snapshot(all_paths, trusted_hashes=trusted_hashes)
    scored, execution = score_teacher_v3_ensemble(
        candidates=candidates,
        native_h5=Path(native_h5).resolve(),
        checkpoint_paths=[Path(value) for value in release["checkpoint_paths"]],
        checkpoint_sha256_by_fold=release["checkpoint_sha256_by_fold"],
        input_provenance=release["input_provenance"],
        calibration_temperature=float(release["calibration_temperature"]),
        batch_size=batch_size,
        workers=workers,
        require_cuda=require_cuda,
    )
    if execution.get("checkpoint_sha256_by_fold") != release[
        "checkpoint_sha256_by_fold"
    ]:
        raise RuntimeError("Teacher-v3 execution checkpoint hashes disagree")
    ranked = rank_planet_enrichment(scored)
    ranked.insert(0, "ranking_schema_version", S63_RANKING_SCHEMA)
    after = _snapshot(all_paths)
    if after != before:
        raise RuntimeError("Teacher-v3 S63 inputs changed during scoring")

    out_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = out_dir.parent / f".{out_dir.name}.tmp-{os.getpid()}"
    if staging.exists():
        raise FileExistsError(f"staging directory already exists: {staging}")
    staging.mkdir(mode=0o700)
    staging.chmod(0o700)
    try:
        score_tmp, staged_score = _stage_table(
            scored,
            out_dir=staging,
            stem="teacher_v3_s63_real_candidate_scores",
        )
        ranking_tmp, staged_ranking = _stage_table(
            ranked,
            out_dir=staging,
            stem="teacher_v3_s63_planet_enrichment_ranked",
        )
        score_tmp.replace(staged_score)
        ranking_tmp.replace(staged_ranking)
        staged_score.chmod(0o600)
        staged_ranking.chmod(0o600)
        initial_checkout_audit = launch_authorization["execution_checkout"]
        try:
            publication_checkout_audit = verify_execution_checkout(
                repo=repo,
                expected_git_sha=expected_git_sha,
                launch_git={
                    "repo": initial_checkout_audit["launch_declared_repo"],
                    "sha": initial_checkout_audit["launch_git_sha"],
                    "checkout_clean": initial_checkout_audit[
                        "checkout_clean"
                    ],
                    "untracked_files_checked": initial_checkout_audit[
                        "untracked_files_checked"
                    ],
                },
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(
                "execution checkout changed during Teacher-v3 S63 scoring"
            ) from exc
        if publication_checkout_audit != initial_checkout_audit:
            raise RuntimeError(
                "execution checkout audit changed during Teacher-v3 S63 scoring"
            )
        final_score = out_dir / staged_score.name
        final_ranking = out_dir / staged_ranking.name
        final_summary = out_dir / "teacher_v3_s63_scoring_summary.json"
        artifacts = {
            name: {
                **record,
                "sha256_before": record["sha256"],
                "sha256_after": after[name]["sha256"],
            }
            for name, record in before.items()
        }
        summary: dict[str, Any] = {
            "schema_version": S63_SUMMARY_SCHEMA,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "strict_provenance_passed": True,
            "atomic_publication": "fresh_staging_directory_renamed_as_one_commit",
            "sector": S63,
            "score_policy": S63_SCORE_POLICY,
            "release_verification": release,
            "launch_authorization": launch_authorization,
            "publication_checkout_audit": publication_checkout_audit,
            "input_verification": input_verification,
            "input_artifacts": artifacts,
            "execution": execution,
            "outputs": {
                "scores": str(final_score),
                "planet_enrichment_ranking": str(final_ranking),
                "summary": str(final_summary),
            },
            "output_sha256": {
                "scores": file_sha256(staged_score),
                "planet_enrichment_ranking": file_sha256(staged_ranking),
            },
            "n_scored_rows": int(len(scored)),
            "n_ranked_unique_tics": int(len(ranked)),
            "output_permissions": {
                "directory_mode": "0700",
                "score_mode": "0600",
                "ranking_mode": "0600",
                "summary_mode": "0600",
            },
        }
        staged_summary = staging / final_summary.name
        staged_summary.write_text(
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        staged_summary.chmod(0o600)
        after_stats = _snapshot(
            all_paths,
            trusted_hashes={
                name: str(record["sha256"]) for name, record in after.items()
            },
        )
        if after_stats != after:
            raise RuntimeError(
                "Teacher-v3 S63 inputs changed before atomic publication"
            )
        if out_dir.exists():
            raise FileExistsError(
                f"Teacher-v3 S63 output path appeared during scoring: {out_dir}"
            )
        staging.replace(out_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return summary


__all__ = [
    "S63_AUTHORIZATION_CONTRACT",
    "S63_LAUNCH_CONTRACT",
    "S63_RANKING_SCHEMA",
    "S63_SCORE_POLICY",
    "S63_SCORE_SCHEMA",
    "S63_SUMMARY_SCHEMA",
    "file_sha256",
    "score_teacher_v3_ensemble",
    "score_teacher_v3_s63_to_disk",
    "verify_execution_checkout",
    "verify_s63_scoring_inputs",
    "verify_s63_launch_authorization",
    "verify_teacher_v3_release",
]
