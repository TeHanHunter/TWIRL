from __future__ import annotations

import json
from pathlib import Path
import stat
import subprocess
from typing import Any

import h5py
import numpy as np
import pandas as pd
import pytest

from twirl.vetting import teacher_v3_inference as inference
from twirl.vetting.harmonic_cnn import MODEL_VERSION
from twirl.vetting.harmonic_dataset import MetadataNormalization
from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION
from twirl.vetting.teacher_v3 import TEACHER_V3_RUN_NAME
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_CALIBRATION_SCHEMA,
    TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
    TEACHER_V3_CHECKPOINT_NAMESPACE,
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
)


ROOT = Path(__file__).resolve().parents[1]
SELECTION_POLICY_FIXTURE = (
    ROOT
    / "reports"
    / "stage5_validation"
    / "teacher_v3_s63_prospective_v1"
    / "preregistered"
    / "selection_policy_v1.json"
)

AT_RISK_S63_METADATA_COLUMNS = (
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


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _at_risk_metadata_normalizations() -> list[MetadataNormalization]:
    normalization = MetadataNormalization(
        columns=AT_RISK_S63_METADATA_COLUMNS,
        center=(0.0,) * len(AT_RISK_S63_METADATA_COLUMNS),
        scale=(1.0,) * len(AT_RISK_S63_METADATA_COLUMNS),
    )
    return [normalization] * 5


def _at_risk_metadata_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            column: np.asarray([1.0, 2.0], dtype=float)
            for column in AT_RISK_S63_METADATA_COLUMNS
        }
    )


@pytest.mark.parametrize("column", AT_RISK_S63_METADATA_COLUMNS)
def test_checkpoint_metadata_rejects_absent_s63_feature(column: str) -> None:
    rows = _at_risk_metadata_rows().drop(columns=column)

    with pytest.raises(KeyError, match=column):
        inference._validate_checkpoint_metadata_inputs(
            rows,
            _at_risk_metadata_normalizations(),
        )


@pytest.mark.parametrize("column", AT_RISK_S63_METADATA_COLUMNS)
def test_checkpoint_metadata_rejects_wholly_nonfinite_s63_feature(
    column: str,
) -> None:
    rows = _at_risk_metadata_rows()
    rows[column] = [np.nan, np.inf]

    with pytest.raises(ValueError, match=column):
        inference._validate_checkpoint_metadata_inputs(
            rows,
            _at_risk_metadata_normalizations(),
        )


def test_checkpoint_metadata_allows_row_level_center_imputation() -> None:
    rows = _at_risk_metadata_rows()
    rows.loc[0, list(AT_RISK_S63_METADATA_COLUMNS)] = np.nan
    normalizations = _at_risk_metadata_normalizations()

    audit = inference._validate_checkpoint_metadata_inputs(rows, normalizations)
    metadata, _ = inference.build_metadata_matrix(
        rows,
        fit_mask=np.zeros(len(rows), dtype=bool),
        normalization=normalizations[0],
    )

    assert audit["passed"] is True
    assert audit["finite_rows_by_column"] == {
        column: 1 for column in AT_RISK_S63_METADATA_COLUMNS
    }
    assert audit["partially_nonfinite_columns"] == sorted(
        AT_RISK_S63_METADATA_COLUMNS
    )
    assert np.array_equal(metadata[0], np.zeros(len(AT_RISK_S63_METADATA_COLUMNS)))
    assert np.array_equal(metadata[1], np.full(len(AT_RISK_S63_METADATA_COLUMNS), 2.0))


def _release_fixture(tmp_path: Path) -> dict[str, Any]:
    root = tmp_path / "release"
    provenance = {
        "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "native_h5_sha256": "1" * 64,
        "native_manifest_sha256": "1" * 64,
        "training_table_sha256": "2" * 64,
        "native_registry_sha256": "3" * 64,
        "native_registry_summary_sha256": "4" * 64,
    }
    calibration = root / TEACHER_V3_PRIMARY_PROFILE / "pooled_oof_calibration.json"
    _write_json(
        calibration,
        {
            "schema_version": TEACHER_V3_CALIBRATION_SCHEMA,
            "run_id": TEACHER_V3_RUN_ID,
            "profile": TEACHER_V3_PRIMARY_PROFILE,
            "scope": "concatenated_five_fold_development_oof_logits",
            "temperature": 1.25,
            **provenance,
        },
    )
    calibration_sha256 = inference.file_sha256(calibration)
    checkpoint_records = []
    checkpoint_paths = []
    for fold in range(5):
        checkpoint = (
            root
            / TEACHER_V3_PRIMARY_PROFILE
            / f"fold_{fold}"
            / "teacher.pt"
        )
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.write_bytes(f"synthetic-fold-{fold}\n".encode())
        checkpoint_paths.append(checkpoint)
        checkpoint_records.append(
            {
                "fold": fold,
                "path": (
                    f"{TEACHER_V3_PRIMARY_PROFILE}/fold_{fold}/teacher.pt"
                ),
                "sha256": inference.file_sha256(checkpoint),
                "pooled_oof_calibration_sha256": calibration_sha256,
            }
        )
    manifest = root / "selected_checkpoint_manifest.json"
    _write_json(
        manifest,
        {
            "schema_version": TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "selected_profile": TEACHER_V3_PRIMARY_PROFILE,
            "selection_policy": "fixed_before_test",
            **provenance,
            "checkpoints": checkpoint_records,
        },
    )
    manifest_sha256 = inference.file_sha256(manifest)
    freeze = root / "pretest_model_freeze.json"
    _write_json(
        freeze,
        {
            "schema_version": "twirl_teacher_v3_model_freeze_v1",
            "run_id": TEACHER_V3_RUN_ID,
            "model_version": MODEL_VERSION,
            "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
            "profile_selection_policy": "fixed_before_test",
            "primary_checkpoint_manifest_sha256": manifest_sha256,
            "test_rows_used_for_selection_or_calibration": False,
            **provenance,
        },
    )
    freeze_sha256 = inference.file_sha256(freeze)
    summary = root / "summary.json"
    _write_json(
        summary,
        {
            "schema_version": "twirl_teacher_v3_training_summary_v1",
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
            "profile_selection_policy": "fixed_before_test",
            "student_training_blocked": True,
            "automatic_production_promotion": False,
            "selected_checkpoint_manifest_sha256": manifest_sha256,
            "pretest_model_freeze_sha256": freeze_sha256,
            "calibration": [
                {
                    "profile": TEACHER_V3_PRIMARY_PROFILE,
                    "calibration_sha256": calibration_sha256,
                }
            ],
            **provenance,
        },
    )
    return {
        "summary": summary,
        "summary_sha256": inference.file_sha256(summary),
        "freeze": freeze,
        "freeze_sha256": freeze_sha256,
        "manifest": manifest,
        "manifest_sha256": manifest_sha256,
        "calibration": calibration,
        "calibration_sha256": calibration_sha256,
        "checkpoints": checkpoint_paths,
        "checkpoint_records": checkpoint_records,
        "provenance": provenance,
    }


def _verify_release(
    case: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> dict[str, Any]:
    def fake_manifest_verifier(**kwargs: Any) -> dict[str, Any]:
        assert kwargs["expected_manifest_sha256"] == case["manifest_sha256"]
        assert kwargs["expected_calibration_sha256"] == case["calibration_sha256"]
        assert kwargs["expected_profile"] == TEACHER_V3_PRIMARY_PROFILE
        assert kwargs["expected_provenance"] == case["provenance"]
        observed = {
            str(record["fold"]): inference.file_sha256(
                case["manifest"].parent / record["path"]
            )
            for record in case["checkpoint_records"]
        }
        assert observed == {
            str(record["fold"]): record["sha256"]
            for record in case["checkpoint_records"]
        }
        return {
            "n_verified_checkpoints": 5,
            "checkpoint_sha256_by_fold": observed,
        }

    monkeypatch.setattr(
        inference,
        "verify_teacher_v3_checkpoint_manifest",
        fake_manifest_verifier,
    )
    return inference.verify_teacher_v3_release(
        release_summary_path=case["summary"],
        expected_release_summary_sha256=case["summary_sha256"],
        pretest_freeze_path=case["freeze"],
        expected_pretest_freeze_sha256=case["freeze_sha256"],
        checkpoint_manifest_path=case["manifest"],
        expected_checkpoint_manifest_sha256=case["manifest_sha256"],
        calibration_path=case["calibration"],
        expected_calibration_sha256=case["calibration_sha256"],
    )


def test_release_verification_binds_four_documents_and_five_checkpoints(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _release_fixture(tmp_path)
    result = _verify_release(case, monkeypatch)

    assert result["input_provenance"] == case["provenance"]
    assert result["calibration_temperature"] == 1.25
    assert result["checkpoint_paths"] == [
        str(path) for path in case["checkpoints"]
    ]
    assert result["checkpoint_sha256_by_fold"] == {
        str(record["fold"]): record["sha256"]
        for record in case["checkpoint_records"]
    }


def test_release_verification_rejects_trust_anchor_or_checkpoint_change(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _release_fixture(tmp_path)
    with pytest.raises(RuntimeError, match="release summary SHA-256 mismatch"):
        inference.verify_teacher_v3_release(
            release_summary_path=case["summary"],
            expected_release_summary_sha256="f" * 64,
            pretest_freeze_path=case["freeze"],
            expected_pretest_freeze_sha256=case["freeze_sha256"],
            checkpoint_manifest_path=case["manifest"],
            expected_checkpoint_manifest_sha256=case["manifest_sha256"],
            calibration_path=case["calibration"],
            expected_calibration_sha256=case["calibration_sha256"],
        )

    case["checkpoints"][2].write_bytes(b"changed\n")
    monkeypatch.setattr(
        inference,
        "verify_teacher_v3_checkpoint_manifest",
        lambda **kwargs: {
            "n_verified_checkpoints": 5,
            "checkpoint_sha256_by_fold": {
                str(record["fold"]): record["sha256"]
                for record in case["checkpoint_records"]
            },
        },
    )
    with pytest.raises(RuntimeError, match="fold 2 checkpoint hash changed"):
        inference.verify_teacher_v3_release(
            release_summary_path=case["summary"],
            expected_release_summary_sha256=case["summary_sha256"],
            pretest_freeze_path=case["freeze"],
            expected_pretest_freeze_sha256=case["freeze_sha256"],
            checkpoint_manifest_path=case["manifest"],
            expected_checkpoint_manifest_sha256=case["manifest_sha256"],
            calibration_path=case["calibration"],
            expected_calibration_sha256=case["calibration_sha256"],
        )


def _s63_inputs(tmp_path: Path) -> dict[str, Any]:
    candidates = tmp_path / "s63_candidates.csv"
    pd.DataFrame(
        {
            "review_id": ["s63:11:0", "s63:12:0"],
            "sector": [63, 63],
            "tic": [11, 12],
            "period_d": [1.1, 2.2],
            "t0_bjd": [2460155.1, 2460155.2],
            "duration_min": [10.0, 12.0],
            "native_input_include": [True, True],
            "source_kind": ["real_candidate", "real_candidate"],
            "native_group_path": [
                "targets/0000000000000011",
                "targets/0000000000000012",
            ],
        }
    ).to_csv(candidates, index=False)
    reserved = tmp_path / "s63_reserved_tics.txt"
    reserved.write_text("11\n12\n99\n")
    native = tmp_path / "s63_native.h5"
    with h5py.File(native, "w") as h5:
        h5.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        h5.attrs["training_table_sha256"] = inference.file_sha256(candidates)
        h5.attrs["periodogram_n"] = inference.S63_NATIVE_PERIODOGRAM_N
        h5.create_group("injections")
        for tic in (11, 12):
            group = h5.create_group(f"targets/{tic:016d}")
            group.attrs["tic"] = tic
            group.attrs["sector"] = 63
    return {"candidates": candidates, "reserved": reserved, "native": native}


def test_s63_inputs_are_real_reserved_and_exactly_group_bound(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _s63_inputs(tmp_path)
    monkeypatch.setattr(
        inference,
        "verify_raw_pair_contract",
        lambda *args, **kwargs: {"passed": True, "failures": []},
    )
    candidates, audit = inference.verify_s63_scoring_inputs(
        candidates_path=case["candidates"],
        expected_candidates_sha256=inference.file_sha256(case["candidates"]),
        reserved_tics_path=case["reserved"],
        expected_reserved_tics_sha256=inference.file_sha256(case["reserved"]),
        native_h5=case["native"],
        expected_native_h5_sha256=inference.file_sha256(case["native"]),
    )
    assert len(candidates) == 2
    assert audit["real_only"] is True
    assert audit["n_exact_native_groups"] == 2

    with h5py.File(case["native"], "a") as h5:
        h5.create_group("injections/synthetic")
    with pytest.raises(ValueError, match="contains injection groups"):
        inference.verify_s63_scoring_inputs(
            candidates_path=case["candidates"],
            expected_candidates_sha256=inference.file_sha256(case["candidates"]),
            reserved_tics_path=case["reserved"],
            expected_reserved_tics_sha256=inference.file_sha256(case["reserved"]),
            native_h5=case["native"],
            expected_native_h5_sha256=inference.file_sha256(case["native"]),
        )


def test_s63_inputs_reject_duplicate_tic_rows_with_unique_review_ids(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _s63_inputs(tmp_path)
    candidates = pd.read_csv(case["candidates"])
    duplicate = candidates.iloc[[0]].copy()
    duplicate["review_id"] = "s63:11:alternate"
    pd.concat([candidates, duplicate], ignore_index=True).to_csv(
        case["candidates"], index=False
    )
    monkeypatch.setattr(
        inference,
        "verify_raw_pair_contract",
        lambda *args, **kwargs: {"passed": True, "failures": []},
    )

    with pytest.raises(ValueError, match="exactly one candidate row per TIC"):
        inference.verify_s63_scoring_inputs(
            candidates_path=case["candidates"],
            expected_candidates_sha256=inference.file_sha256(case["candidates"]),
            reserved_tics_path=case["reserved"],
            expected_reserved_tics_sha256=inference.file_sha256(case["reserved"]),
            native_h5=case["native"],
            expected_native_h5_sha256=inference.file_sha256(case["native"]),
        )


def _fake_release_paths(tmp_path: Path) -> dict[str, Any]:
    paths = {}
    for name in (
        "release_summary",
        "pretest_freeze",
        "checkpoint_manifest",
        "calibration",
    ):
        path = tmp_path / f"{name}.json"
        path.write_text(f"{name}\n")
        paths[name] = path
    checkpoints = []
    hashes = {}
    for fold in range(5):
        path = tmp_path / f"fold_{fold}.pt"
        path.write_text(f"fold {fold}\n")
        checkpoints.append(path)
        hashes[str(fold)] = inference.file_sha256(path)
    return {
        "schema_version": "twirl_teacher_v3_verified_release_v1",
        "release_documents": {
            name: {"path": str(path), "sha256": inference.file_sha256(path)}
            for name, path in paths.items()
        },
        "input_provenance": {
            "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
            "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        },
        "calibration_temperature": 1.25,
        "checkpoint_paths": [str(path) for path in checkpoints],
        "checkpoint_sha256_by_fold": hashes,
        "manifest_audit": {"n_verified_checkpoints": 5},
    }


def _artifact(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": inference.file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def _git_repo(path: Path) -> tuple[Path, str]:
    path.mkdir(parents=True)
    (path / "tracked.txt").write_text("clean\n")
    subprocess.run(["git", "init", "-q"], cwd=path, check=True)
    subprocess.run(["git", "add", "tracked.txt"], cwd=path, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=TWIRL Test",
            "-c",
            "user.email=twirl@example.invalid",
            "-c",
            "commit.gpgsign=false",
            "commit",
            "-qm",
            "test checkout",
        ],
        cwd=path,
        check=True,
    )
    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=path,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return path, git_sha


def _launch_fixture(
    tmp_path: Path,
    *,
    release: dict[str, Any],
    candidates: Path,
    reserved: Path,
    native: Path,
) -> dict[str, Any]:
    repo, git_sha = _git_repo(tmp_path / "execution_repo")
    plan = tmp_path / "prospective_plan.json"
    _write_json(
        plan,
        {
            "schema_version": inference.PROSPECTIVE_PLAN_SCHEMA,
            "status": "plan_frozen_launch_manifest_pending",
            "sector": 63,
            "frozen_teacher_v3": {
                "run_id": TEACHER_V3_RUN_ID,
                "release_name": TEACHER_V3_RUN_NAME,
                "model_version": MODEL_VERSION,
                "selected_profile": TEACHER_V3_PRIMARY_PROFILE,
                "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
                "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
                "selection_policy": "fixed_before_test",
                "temperature_calibration_scope": (
                    "concatenated_five_fold_development_oof_logits"
                ),
                "temperature": release["calibration_temperature"],
                "summary_sha256": release["release_documents"][
                    "release_summary"
                ]["sha256"],
                "pretest_model_freeze_sha256": release["release_documents"][
                    "pretest_freeze"
                ]["sha256"],
                "selected_checkpoint_manifest_sha256": release[
                    "release_documents"
                ]["checkpoint_manifest"]["sha256"],
                "pooled_oof_calibration_sha256": release["release_documents"][
                    "calibration"
                ]["sha256"],
                "checkpoints": [
                    {"fold": fold, "sha256": release["checkpoint_sha256_by_fold"][str(fold)]}
                    for fold in range(5)
                ],
                "automatic_promotion": False,
                "student_training_authorized": False,
            },
            "s63_identity_reservation": {
                "identity_only": True,
                "light_curves_opened_when_reserved": False,
                "reserved_tics_sha256": inference.file_sha256(reserved),
            },
        },
    )
    selection_policy = tmp_path / "selection_policy_v1.json"
    selection_policy.parent.mkdir(parents=True, exist_ok=True)
    selection_policy.write_text(SELECTION_POLICY_FIXTURE.read_text())
    launch = tmp_path / "launch_manifest.json"
    passed_checks = {
        name: {"passed": True}
        for name in (
            "accepted_stage1",
            "cadence",
            "compact",
            "bls",
            "candidates",
            "raw_source",
            "native",
            "inventories",
            "selection_policy",
            "model_ready_derivation",
        )
    }
    _write_json(
        launch,
        {
            "contract_version": inference.S63_LAUNCH_CONTRACT,
            "authorization_contract": inference.S63_AUTHORIZATION_CONTRACT,
            "sector": 63,
            "orbits": [133, 134],
            "status": "preprocessing_complete_scoring_not_started",
            "score_or_queue_read": False,
            "passed": True,
            "git": {
                "repo": str(repo),
                "sha": git_sha,
                "checkout_clean": True,
                "untracked_files_checked": True,
            },
            "artifacts": {
                "preregistered_contract": _artifact(plan),
                "selection_policy": _artifact(selection_policy),
                "candidates": _artifact(candidates),
                "reserved_tics": _artifact(reserved),
                "native_h5": _artifact(native),
            },
            "checks": passed_checks,
        },
    )
    return {
        "plan": plan,
        "selection_policy": selection_policy,
        "launch": launch,
        "repo": repo,
        "git_sha": git_sha,
    }


def _verify_launch(
    *,
    case: dict[str, Any],
    release: dict[str, Any],
    candidates: Path,
    reserved: Path,
    native: Path,
) -> dict[str, Any]:
    return inference.verify_s63_launch_authorization(
        repo=case["repo"],
        expected_git_sha=case["git_sha"],
        launch_manifest_path=case["launch"],
        expected_launch_manifest_sha256=inference.file_sha256(case["launch"]),
        preregistered_contract_path=case["plan"],
        selection_policy_path=case["selection_policy"],
        release_verification=release,
        candidates_path=candidates,
        expected_candidates_sha256=inference.file_sha256(candidates),
        reserved_tics_path=reserved,
        expected_reserved_tics_sha256=inference.file_sha256(reserved),
        native_h5=native,
        expected_native_h5_sha256=inference.file_sha256(native),
    )


def test_launch_authorization_binds_plan_inputs_and_exact_teacher_release(
    tmp_path: Path,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates = tmp_path / "candidate.csv"
    candidates.write_text("tic\n11\n")
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native\n")
    case = _launch_fixture(
        tmp_path,
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )

    audit = _verify_launch(
        case=case,
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    assert audit["passed"] is True
    assert audit["score_or_queue_read_before_authorization"] is False
    assert audit["teacher_checkpoint_sha256_by_fold"] == release[
        "checkpoint_sha256_by_fold"
    ]
    checkout = audit["execution_checkout"]
    assert checkout["passed"] is True
    assert checkout["observed_git_sha"] == case["git_sha"]
    assert checkout["operator_expected_git_sha"] == case["git_sha"]
    assert checkout["launch_git_sha"] == case["git_sha"]
    assert checkout["untracked_files_checked"] is True


def test_launch_authorization_rejects_clean_different_commit(
    tmp_path: Path,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates = tmp_path / "candidate.csv"
    candidates.write_text("tic\n11\n")
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native\n")
    case = _launch_fixture(
        tmp_path,
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    repo = case["repo"]
    (repo / "tracked.txt").write_text("different clean commit\n")
    subprocess.run(["git", "add", "tracked.txt"], cwd=repo, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=TWIRL Test",
            "-c",
            "user.email=twirl@example.invalid",
            "-c",
            "commit.gpgsign=false",
            "commit",
            "-qm",
            "different checkout",
        ],
        cwd=repo,
        check=True,
    )
    case["git_sha"] = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()

    with pytest.raises(ValueError, match="Git SHA binding differs"):
        _verify_launch(
            case=case,
            release=release,
            candidates=candidates,
            reserved=reserved,
            native=native,
        )


@pytest.mark.parametrize("dirty_kind", ("tracked", "untracked"))
def test_launch_authorization_rejects_dirty_or_untracked_checkout(
    tmp_path: Path,
    dirty_kind: str,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates = tmp_path / "candidate.csv"
    candidates.write_text("tic\n11\n")
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native\n")
    case = _launch_fixture(
        tmp_path,
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    repo = case["repo"]
    if dirty_kind == "tracked":
        (repo / "tracked.txt").write_text("dirty\n")
    else:
        (repo / "untracked.txt").write_text("untracked\n")

    with pytest.raises(ValueError, match="not fully clean"):
        _verify_launch(
            case=case,
            release=release,
            candidates=candidates,
            reserved=reserved,
            native=native,
        )


@pytest.mark.parametrize(
    "field", ("checkout_clean", "untracked_files_checked")
)
def test_launch_authorization_requires_declared_full_git_audit(
    tmp_path: Path,
    field: str,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates = tmp_path / "candidate.csv"
    candidates.write_text("tic\n11\n")
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native\n")
    case = _launch_fixture(
        tmp_path,
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    launch = json.loads(case["launch"].read_text())
    launch["git"][field] = False
    _write_json(case["launch"], launch)

    with pytest.raises(ValueError, match=rf"git\.{field} must be true"):
        _verify_launch(
            case=case,
            release=release,
            candidates=candidates,
            reserved=reserved,
            native=native,
        )


def test_launch_authorization_rejects_hash_status_and_teacher_drift(
    tmp_path: Path,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates = tmp_path / "candidate.csv"
    candidates.write_text("tic\n11\n")
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native\n")
    case = _launch_fixture(
        tmp_path,
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    with pytest.raises(RuntimeError, match="launch manifest SHA-256 mismatch"):
        inference.verify_s63_launch_authorization(
            repo=case["repo"],
            expected_git_sha=case["git_sha"],
            launch_manifest_path=case["launch"],
            expected_launch_manifest_sha256="f" * 64,
            preregistered_contract_path=case["plan"],
            selection_policy_path=case["selection_policy"],
            release_verification=release,
            candidates_path=candidates,
            expected_candidates_sha256=inference.file_sha256(candidates),
            reserved_tics_path=reserved,
            expected_reserved_tics_sha256=inference.file_sha256(reserved),
            native_h5=native,
            expected_native_h5_sha256=inference.file_sha256(native),
        )

    launch = json.loads(case["launch"].read_text())
    launch["score_or_queue_read"] = True
    _write_json(case["launch"], launch)
    with pytest.raises(ValueError, match="score_or_queue_read"):
        _verify_launch(
            case=case,
            release=release,
            candidates=candidates,
            reserved=reserved,
            native=native,
        )

    case = _launch_fixture(
        tmp_path / "teacher_drift",
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    plan = json.loads(case["plan"].read_text())
    plan["frozen_teacher_v3"]["checkpoints"][3]["sha256"] = "e" * 64
    _write_json(case["plan"], plan)
    launch = json.loads(case["launch"].read_text())
    launch["artifacts"]["preregistered_contract"] = _artifact(case["plan"])
    _write_json(case["launch"], launch)
    with pytest.raises(ValueError, match="fold-3 differs"):
        _verify_launch(
            case=case,
            release=release,
            candidates=candidates,
            reserved=reserved,
            native=native,
        )

    case = _launch_fixture(
        tmp_path / "missing_policy",
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    launch = json.loads(case["launch"].read_text())
    launch["artifacts"].pop("selection_policy")
    _write_json(case["launch"], launch)
    with pytest.raises(ValueError, match="artifacts.selection_policy"):
        _verify_launch(
            case=case,
            release=release,
            candidates=candidates,
            reserved=reserved,
            native=native,
        )


@pytest.mark.parametrize(
    ("check_name", "failure_mode"),
    (
        ("selection_policy", "missing"),
        ("selection_policy", "failed"),
        ("model_ready_derivation", "missing"),
        ("model_ready_derivation", "failed"),
    ),
)
def test_launch_authorization_requires_new_pre_score_checks(
    tmp_path: Path,
    check_name: str,
    failure_mode: str,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates = tmp_path / "candidate.csv"
    candidates.write_text("tic\n11\n")
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native\n")
    case = _launch_fixture(
        tmp_path,
        release=release,
        candidates=candidates,
        reserved=reserved,
        native=native,
    )
    launch = json.loads(case["launch"].read_text())
    if failure_mode == "missing":
        launch["checks"].pop(check_name)
    else:
        launch["checks"][check_name]["passed"] = False
    _write_json(case["launch"], launch)

    with pytest.raises(
        ValueError, match=rf"launch check {check_name} did not pass"
    ):
        _verify_launch(
            case=case,
            release=release,
            candidates=candidates,
            reserved=reserved,
            native=native,
        )


def test_s63_outputs_publish_as_one_fresh_hash_bound_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates_path = tmp_path / "candidates.csv"
    candidates = pd.DataFrame(
        {
            "review_id": ["s63:11", "s63:12"],
            "sector": [63, 63],
            "tic": [11, 12],
            "period_d": [1.0, 2.0],
            "t0_bjd": [2460155.0, 2460156.0],
            "duration_min": [10.0, 12.0],
            "native_input_include": [True, True],
        }
    )
    candidates.to_csv(candidates_path, index=False)
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n12\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"synthetic-native\n")
    launch = tmp_path / "launch.json"
    launch.write_text("launch\n")
    plan = tmp_path / "plan.json"
    plan.write_text("plan\n")
    policy = tmp_path / "selection_policy.json"
    policy.write_text("policy\n")
    repo, expected_git_sha = _git_repo(tmp_path / "execution_repo")
    execution_checkout = inference.verify_execution_checkout(
        repo=repo,
        expected_git_sha=expected_git_sha,
        launch_git={
            "repo": str(repo),
            "sha": expected_git_sha,
            "checkout_clean": True,
            "untracked_files_checked": True,
        },
    )
    launch_authorization = {
        "passed": True,
        "execution_checkout": execution_checkout,
        "launch_manifest": {
            "path": str(launch),
            "sha256": inference.file_sha256(launch),
        },
        "preregistered_contract": {
            "path": str(plan),
            "sha256": inference.file_sha256(plan),
        },
        "selection_policy": {
            "path": str(policy),
            "sha256": inference.file_sha256(policy),
        },
    }
    monkeypatch.setattr(inference, "verify_teacher_v3_release", lambda **kwargs: release)
    monkeypatch.setattr(
        inference,
        "verify_s63_launch_authorization",
        lambda **kwargs: launch_authorization,
    )
    monkeypatch.setattr(
        inference,
        "verify_s63_scoring_inputs",
        lambda **kwargs: (
            candidates,
            {
                "schema_version": "twirl_teacher_v3_s63_input_verification_v1",
                "passed": True,
            },
        ),
    )

    def fake_score(**kwargs: Any) -> tuple[pd.DataFrame, dict[str, Any]]:
        scored = candidates.copy()
        scored.insert(0, "score_schema_version", inference.S63_SCORE_SCHEMA)
        scored["p_planet_like"] = [0.9, 0.4]
        scored["p_preserve"] = [0.8, 0.7]
        scored["std_p_planet_like"] = [0.01, 0.02]
        scored["sde_max"] = [12.0, 9.0]
        return scored, {
            "schema_version": "twirl_teacher_v3_s63_ensemble_execution_v1",
            "checkpoint_sha256_by_fold": release[
                "checkpoint_sha256_by_fold"
            ],
        }

    monkeypatch.setattr(inference, "score_teacher_v3_ensemble", fake_score)
    out_dir = tmp_path / "published"
    summary = inference.score_teacher_v3_s63_to_disk(
        repo=repo,
        expected_git_sha=expected_git_sha,
        launch_manifest_path=launch,
        expected_launch_manifest_sha256=inference.file_sha256(launch),
        preregistered_contract_path=plan,
        selection_policy_path=policy,
        release_summary_path=tmp_path / "unused-summary",
        expected_release_summary_sha256="0" * 64,
        pretest_freeze_path=tmp_path / "unused-freeze",
        expected_pretest_freeze_sha256="0" * 64,
        checkpoint_manifest_path=tmp_path / "unused-manifest",
        expected_checkpoint_manifest_sha256="0" * 64,
        calibration_path=tmp_path / "unused-calibration",
        expected_calibration_sha256="0" * 64,
        candidates_path=candidates_path,
        expected_candidates_sha256=inference.file_sha256(candidates_path),
        reserved_tics_path=reserved,
        expected_reserved_tics_sha256=inference.file_sha256(reserved),
        native_h5=native,
        expected_native_h5_sha256=inference.file_sha256(native),
        out_dir=out_dir,
        require_cuda=False,
    )

    assert summary["schema_version"] == inference.S63_SUMMARY_SCHEMA
    assert summary["strict_provenance_passed"] is True
    assert out_dir.is_dir()
    assert not list(tmp_path.glob(".published.tmp-*"))
    summary_path = out_dir / "teacher_v3_s63_scoring_summary.json"
    on_disk = json.loads(summary_path.read_text())
    assert on_disk["output_sha256"] == summary["output_sha256"]
    assert on_disk["launch_authorization"]["execution_checkout"] == (
        launch_authorization["execution_checkout"]
    )
    assert on_disk["publication_checkout_audit"] == execution_checkout
    score_path = Path(on_disk["outputs"]["scores"])
    ranking_path = Path(on_disk["outputs"]["planet_enrichment_ranking"])
    assert score_path.is_file() and ranking_path.is_file()
    assert inference.file_sha256(score_path) == on_disk["output_sha256"]["scores"]
    assert inference.file_sha256(ranking_path) == on_disk["output_sha256"][
        "planet_enrichment_ranking"
    ]
    assert stat.S_IMODE(out_dir.stat().st_mode) == 0o700
    for private_path in (score_path, ranking_path, summary_path):
        assert stat.S_IMODE(private_path.stat().st_mode) == 0o600
    with pytest.raises(FileExistsError, match="fresh path"):
        inference.score_teacher_v3_s63_to_disk(
            repo=repo,
            expected_git_sha=expected_git_sha,
            launch_manifest_path=launch,
            expected_launch_manifest_sha256=inference.file_sha256(launch),
            preregistered_contract_path=plan,
            selection_policy_path=policy,
            release_summary_path=tmp_path / "unused-summary",
            expected_release_summary_sha256="0" * 64,
            pretest_freeze_path=tmp_path / "unused-freeze",
            expected_pretest_freeze_sha256="0" * 64,
            checkpoint_manifest_path=tmp_path / "unused-manifest",
            expected_checkpoint_manifest_sha256="0" * 64,
            calibration_path=tmp_path / "unused-calibration",
            expected_calibration_sha256="0" * 64,
            candidates_path=candidates_path,
            expected_candidates_sha256=inference.file_sha256(candidates_path),
            reserved_tics_path=reserved,
            expected_reserved_tics_sha256=inference.file_sha256(reserved),
            native_h5=native,
            expected_native_h5_sha256=inference.file_sha256(native),
            out_dir=out_dir,
            require_cuda=False,
        )


def test_checkout_change_during_scoring_prevents_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates_path = tmp_path / "candidates.csv"
    candidates = pd.DataFrame(
        {
            "review_id": ["s63:11"],
            "sector": [63],
            "tic": [11],
            "period_d": [1.0],
            "t0_bjd": [2460155.0],
            "duration_min": [10.0],
            "native_input_include": [True],
        }
    )
    candidates.to_csv(candidates_path, index=False)
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native\n")
    launch = tmp_path / "launch.json"
    launch.write_text("launch\n")
    plan = tmp_path / "plan.json"
    plan.write_text("plan\n")
    policy = tmp_path / "selection_policy.json"
    policy.write_text("policy\n")
    repo, expected_git_sha = _git_repo(tmp_path / "execution_repo")
    launch_git = {
        "repo": str(repo),
        "sha": expected_git_sha,
        "checkout_clean": True,
        "untracked_files_checked": True,
    }
    launch_authorization = {
        "passed": True,
        "execution_checkout": inference.verify_execution_checkout(
            repo=repo,
            expected_git_sha=expected_git_sha,
            launch_git=launch_git,
        ),
        "launch_manifest": {
            "path": str(launch),
            "sha256": inference.file_sha256(launch),
        },
        "preregistered_contract": {
            "path": str(plan),
            "sha256": inference.file_sha256(plan),
        },
        "selection_policy": {
            "path": str(policy),
            "sha256": inference.file_sha256(policy),
        },
    }
    monkeypatch.setattr(inference, "verify_teacher_v3_release", lambda **kwargs: release)
    monkeypatch.setattr(
        inference,
        "verify_s63_launch_authorization",
        lambda **kwargs: launch_authorization,
    )
    monkeypatch.setattr(
        inference,
        "verify_s63_scoring_inputs",
        lambda **kwargs: (candidates, {"passed": True}),
    )

    def mutating_score(**kwargs: Any) -> tuple[pd.DataFrame, dict[str, Any]]:
        (repo / "tracked.txt").write_text("changed during scoring\n")
        scored = candidates.copy()
        scored["p_planet_like"] = 0.9
        scored["p_preserve"] = 0.8
        return scored, {
            "checkpoint_sha256_by_fold": release[
                "checkpoint_sha256_by_fold"
            ]
        }

    monkeypatch.setattr(inference, "score_teacher_v3_ensemble", mutating_score)
    out_dir = tmp_path / "must_not_publish_checkout_drift"
    with pytest.raises(
        RuntimeError, match="execution checkout changed during"
    ):
        inference.score_teacher_v3_s63_to_disk(
            repo=repo,
            expected_git_sha=expected_git_sha,
            launch_manifest_path=launch,
            expected_launch_manifest_sha256=inference.file_sha256(launch),
            preregistered_contract_path=plan,
            selection_policy_path=policy,
            release_summary_path=tmp_path / "unused-summary",
            expected_release_summary_sha256="0" * 64,
            pretest_freeze_path=tmp_path / "unused-freeze",
            expected_pretest_freeze_sha256="0" * 64,
            checkpoint_manifest_path=tmp_path / "unused-manifest",
            expected_checkpoint_manifest_sha256="0" * 64,
            calibration_path=tmp_path / "unused-calibration",
            expected_calibration_sha256="0" * 64,
            candidates_path=candidates_path,
            expected_candidates_sha256=inference.file_sha256(candidates_path),
            reserved_tics_path=reserved,
            expected_reserved_tics_sha256=inference.file_sha256(reserved),
            native_h5=native,
            expected_native_h5_sha256=inference.file_sha256(native),
            out_dir=out_dir,
            require_cuda=False,
        )
    assert not out_dir.exists()
    assert not list(tmp_path.glob(".must_not_publish_checkout_drift.tmp-*"))


def test_s63_input_change_prevents_any_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    release = _fake_release_paths(tmp_path)
    candidates_path = tmp_path / "candidates.csv"
    candidates = pd.DataFrame(
        {
            "review_id": ["s63:11"],
            "sector": [63],
            "tic": [11],
            "period_d": [1.0],
            "t0_bjd": [2460155.0],
            "duration_min": [10.0],
            "native_input_include": [True],
        }
    )
    candidates.to_csv(candidates_path, index=False)
    reserved = tmp_path / "reserved.txt"
    reserved.write_text("11\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"before\n")
    launch = tmp_path / "launch.json"
    launch.write_text("launch\n")
    plan = tmp_path / "plan.json"
    plan.write_text("plan\n")
    policy = tmp_path / "selection_policy.json"
    policy.write_text("policy\n")
    expected_git_sha = "b" * 40
    launch_authorization = {
        "passed": True,
        "execution_checkout": {
            "passed": True,
            "observed_git_sha": expected_git_sha,
            "operator_expected_git_sha": expected_git_sha,
            "launch_git_sha": expected_git_sha,
            "checkout_clean": True,
            "untracked_files_checked": True,
        },
        "launch_manifest": {
            "path": str(launch),
            "sha256": inference.file_sha256(launch),
        },
        "preregistered_contract": {
            "path": str(plan),
            "sha256": inference.file_sha256(plan),
        },
        "selection_policy": {
            "path": str(policy),
            "sha256": inference.file_sha256(policy),
        },
    }
    monkeypatch.setattr(inference, "verify_teacher_v3_release", lambda **kwargs: release)
    monkeypatch.setattr(
        inference,
        "verify_s63_launch_authorization",
        lambda **kwargs: launch_authorization,
    )
    monkeypatch.setattr(
        inference,
        "verify_s63_scoring_inputs",
        lambda **kwargs: (candidates, {"passed": True}),
    )

    def mutating_score(**kwargs: Any) -> tuple[pd.DataFrame, dict[str, Any]]:
        native.write_bytes(b"after\n")
        scored = candidates.copy()
        scored["p_planet_like"] = 0.9
        scored["p_preserve"] = 0.8
        return scored, {
            "checkpoint_sha256_by_fold": release[
                "checkpoint_sha256_by_fold"
            ]
        }

    monkeypatch.setattr(inference, "score_teacher_v3_ensemble", mutating_score)
    out_dir = tmp_path / "must_not_publish"
    with pytest.raises(RuntimeError, match="inputs changed during scoring"):
        inference.score_teacher_v3_s63_to_disk(
            repo=tmp_path,
            expected_git_sha=expected_git_sha,
            launch_manifest_path=launch,
            expected_launch_manifest_sha256=inference.file_sha256(launch),
            preregistered_contract_path=plan,
            selection_policy_path=policy,
            release_summary_path=tmp_path / "unused-summary",
            expected_release_summary_sha256="0" * 64,
            pretest_freeze_path=tmp_path / "unused-freeze",
            expected_pretest_freeze_sha256="0" * 64,
            checkpoint_manifest_path=tmp_path / "unused-manifest",
            expected_checkpoint_manifest_sha256="0" * 64,
            calibration_path=tmp_path / "unused-calibration",
            expected_calibration_sha256="0" * 64,
            candidates_path=candidates_path,
            expected_candidates_sha256=inference.file_sha256(candidates_path),
            reserved_tics_path=reserved,
            expected_reserved_tics_sha256=inference.file_sha256(reserved),
            native_h5=native,
            expected_native_h5_sha256=inference.file_sha256(native),
            out_dir=out_dir,
            require_cuda=False,
        )
    assert not out_dir.exists()
