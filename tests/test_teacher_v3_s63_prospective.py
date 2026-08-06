from __future__ import annotations

import hashlib
import json
from pathlib import Path
import stat

import numpy as np
import pandas as pd
import pytest

from twirl.vetting import teacher_v3_prospective as prospective
from twirl.vetting.adp_only import ADP_ONLY_APERTURES
from twirl.vetting.harmonic_cnn import MODEL_VERSION
from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION
from twirl.vetting.teacher_v3 import TEACHER_V3_RUN_NAME
from twirl.vetting.teacher_v3_prospective import (
    CANDIDATE_CONTRACT_VERSION,
    PRIMARY_COHORT,
    PRIMARY_QUOTAS,
    PROSPECTIVE_PLAN_SCHEMA,
    PROSPECTIVE_USE_SCOPE,
    PUBLIC_REVIEW_QUEUE_COLUMNS,
    QUEUE_POLICY_VERSION,
    REPEATED_HOST_COHORT,
    REPEATED_HOST_QUOTAS,
    S63_SCORE_POLICY,
    S63_TARGET_METADATA_CONTRACT,
    TEACHER_V3_CANDIDATE_METADATA_COLUMNS,
    build_s63_prospective_review_queue,
    build_s63_rank_one_candidates,
    derive_s63_prospective_cohorts,
    file_sha256,
    load_prospective_plan,
    load_selection_policy,
    read_tic_inventory,
    tic_inventory_sha256,
    validate_s63_bls_summary,
    validate_queue_launch_bindings,
    verify_s63_prospective_review_queue,
    verify_s63_prospective_review_bundle,
    write_s63_prospective_cohorts,
    write_s63_prospective_review_queue,
    write_s63_rank_one_candidates,
)
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_CHECKPOINT_NAMESPACE,
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
)


ROOT = Path(__file__).resolve().parents[1]
SELECTION_POLICY_PATH = (
    ROOT
    / "reports"
    / "stage5_validation"
    / "teacher_v3_s63_prospective_v1"
    / "preregistered"
    / "selection_policy_v1.json"
)


def _hash(label: str) -> str:
    return hashlib.sha256(label.encode("utf-8")).hexdigest()


def _git_sha(label: str = "git") -> str:
    return hashlib.sha1(label.encode("utf-8")).hexdigest()


def _plan() -> dict[str, object]:
    checkpoints = [
        {"fold": fold, "sha256": _hash(f"checkpoint-{fold}")}
        for fold in range(5)
    ]
    return {
        "schema_version": PROSPECTIVE_PLAN_SCHEMA,
        "sector": 63,
        "search_contract": {
            "apertures": list(ADP_ONLY_APERTURES),
            "anchor_aperture": ADP_ONLY_APERTURES[0],
            "context_aperture": ADP_ONLY_APERTURES[1],
            "n_periods": 50_000,
            "n_retained_peaks_per_aperture": 10,
            "candidate_peak_rank": 1,
            "orbitid_policy": "strict",
            "absolute_probability_threshold": None,
        },
        "blinded_queue": {
            "selection_seed": 630056,
            "primary": {
                "n_tics": PRIMARY_QUOTAS.total,
                "quotas": PRIMARY_QUOTAS.__dict__,
            },
            "repeated_host_side_cohort": {
                "n_tics": REPEATED_HOST_QUOTAS.total,
                "quotas": REPEATED_HOST_QUOTAS.__dict__,
            },
        },
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
            "temperature": 1.5,
            "summary_sha256": _hash("release-summary"),
            "pretest_model_freeze_sha256": _hash("pretest-freeze"),
            "selected_checkpoint_manifest_sha256": _hash(
                "checkpoint-manifest"
            ),
            "pooled_oof_calibration_sha256": _hash("calibration"),
            "checkpoints": checkpoints,
            "automatic_promotion": False,
            "student_training_authorized": False,
        },
        "frozen_training_identity": {
            "morphology_corpus_sha256": _hash("corpus"),
            "n_corpus_tics": 3,
        },
        "s63_identity_reservation": {
            "reserved_tics_sha256": _hash("reservation"),
            "n_requested_tics": 5,
        },
    }


def _candidate_peaks(tics: list[int]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for tic in tics:
        for aperture, scale in zip(ADP_ONLY_APERTURES, (1.0, 2.0), strict=True):
            row: dict[str, object] = {
                "sector": 63,
                "tic": tic,
                "aperture": aperture,
                "peak_rank": 1,
                "status": "ok",
                "source_product_tag": "A2v1",
                "period_d": 1.0 + tic / 10_000.0,
                "t0_bjd": 2459000.0 + tic / 100_000.0,
                "duration_min": 10.0,
                "depth": 0.01 * scale,
                "depth_snr": 8.0 / scale,
                "sde": 20.0 / scale,
                "log_power": 4.0 / scale,
                "target_metadata_contract_version": S63_TARGET_METADATA_CONTRACT,
            }
            for metadata_index, field in enumerate(
                TEACHER_V3_CANDIDATE_METADATA_COLUMNS, start=1
            ):
                row[field] = metadata_index / 100.0
            rows.append(row)
    return pd.DataFrame(rows)


def _candidate_hashes() -> dict[str, str]:
    return {
        "prospective_plan_sha256": _hash("plan"),
        "model_ready_allowlist_sha256": _hash("allowlist"),
        "bls_merged_table_sha256": _hash("bls"),
    }


def _queue_hashes() -> dict[str, str]:
    return {
        "prospective_plan_sha256": _hash("plan"),
        "selection_policy_sha256": file_sha256(SELECTION_POLICY_PATH),
        "launch_manifest_sha256": _hash("launch"),
        "reserved_tics_sha256": _hash("reservation"),
        "teacher_v3_corpus_sha256": _hash("corpus"),
        "model_ready_allowlist_sha256": _hash("allowlist"),
        "primary_cohort_sha256": _hash("primary"),
        "repeated_host_cohort_sha256": _hash("repeated"),
        "cohort_summary_sha256": _hash("cohort-summary"),
        "candidate_table_sha256": _hash("candidate"),
        "teacher_score_table_sha256": _hash("scores"),
        "teacher_score_summary_sha256": _hash("score-summary"),
    }


def _selection_policy() -> dict[str, object]:
    policy, _ = load_selection_policy(
        SELECTION_POLICY_PATH,
        expected_sha256=file_sha256(SELECTION_POLICY_PATH),
    )
    return policy


def _launch_and_scoring_summary(
    hashes: dict[str, str],
) -> tuple[dict[str, object], dict[str, object]]:
    plan = _plan()
    teacher = plan["frozen_teacher_v3"]
    assert isinstance(teacher, dict)
    checkpoint_hashes = {
        str(record["fold"]): str(record["sha256"])
        for record in teacher["checkpoints"]
    }
    release_documents = {
        "release_summary": {
            "path": "/synthetic/release-summary.json",
            "sha256": teacher["summary_sha256"],
        },
        "pretest_freeze": {
            "path": "/synthetic/pretest-freeze.json",
            "sha256": teacher["pretest_model_freeze_sha256"],
        },
        "checkpoint_manifest": {
            "path": "/synthetic/checkpoint-manifest.json",
            "sha256": teacher["selected_checkpoint_manifest_sha256"],
        },
        "calibration": {
            "path": "/synthetic/calibration.json",
            "sha256": teacher["pooled_oof_calibration_sha256"],
        },
    }
    artifact_map = {
        "preregistered_contract": "prospective_plan_sha256",
        "selection_policy": "selection_policy_sha256",
        "reserved_tics": "reserved_tics_sha256",
        "teacher_v3_corpus": "teacher_v3_corpus_sha256",
        "model_ready_allowlist": "model_ready_allowlist_sha256",
        "primary_cohort": "primary_cohort_sha256",
        "repeated_host_cohort": "repeated_host_cohort_sha256",
        "cohort_summary": "cohort_summary_sha256",
        "candidates": "candidate_table_sha256",
    }
    artifacts = {
        name: {
            "path": f"/synthetic/{name}",
            "sha256": hashes[hash_name],
            "size_bytes": 1,
        }
        for name, hash_name in artifact_map.items()
    }
    checks = {
        name: {"passed": True}
        for name in (
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
    }
    launch = {
        "contract_version": "twirl_teacher_v3_s63_prospective_launch_v1",
        "authorization_contract": "s63_teacher_v3_prospective_v1",
        "sector": 63,
        "orbits": [133, 134],
        "status": "preprocessing_complete_scoring_not_started",
        "score_or_queue_read": False,
        "passed": True,
        "producer_git_sha": _git_sha(),
        "git": {
            "repo": "/synthetic/repo",
            "sha": _git_sha(),
            "checkout_clean": True,
            "untracked_files_checked": True,
        },
        "artifacts": artifacts,
        "checks": checks,
    }
    score_inputs = {
        name: {
            "sha256": hashes[hash_name],
            "sha256_before": hashes[hash_name],
            "sha256_after": hashes[hash_name],
        }
        for name, hash_name in (
            ("s63_launch_manifest", "launch_manifest_sha256"),
            ("s63_preregistered_contract", "prospective_plan_sha256"),
            ("candidate_table", "candidate_table_sha256"),
        )
    }
    scoring = {
        "schema_version": "twirl_teacher_v3_s63_scoring_summary_v1",
        "sector": 63,
        "strict_provenance_passed": True,
        "score_policy": S63_SCORE_POLICY,
        "release_verification": {
            "schema_version": "twirl_teacher_v3_verified_release_v1",
            "release_documents": release_documents,
            "input_provenance": {
                "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
                "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
            },
            "calibration_temperature": teacher["temperature"],
            "checkpoint_sha256_by_fold": checkpoint_hashes,
            "manifest_audit": {
                "n_verified_checkpoints": 5,
                "checkpoint_sha256_by_fold": checkpoint_hashes,
            },
        },
        "execution": {
            "schema_version": "twirl_teacher_v3_s63_ensemble_execution_v1",
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "profile": TEACHER_V3_PRIMARY_PROFILE,
            "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
            "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
            "sector": 63,
            "checkpoint_sha256_by_fold": checkpoint_hashes,
        },
        "output_sha256": {"scores": hashes["teacher_score_table_sha256"]},
        "launch_authorization": {
            "passed": True,
            "execution_checkout": _execution_git_audit(launch),
            "teacher_checkpoint_sha256_by_fold": checkpoint_hashes,
            "teacher_release_document_bindings": {
                "summary_sha256": teacher["summary_sha256"],
                "pretest_model_freeze_sha256": teacher[
                    "pretest_model_freeze_sha256"
                ],
                "selected_checkpoint_manifest_sha256": teacher[
                    "selected_checkpoint_manifest_sha256"
                ],
                "pooled_oof_calibration_sha256": teacher[
                    "pooled_oof_calibration_sha256"
                ],
            },
            "launch_manifest": {"sha256": hashes["launch_manifest_sha256"]},
            "preregistered_contract": {
                "sha256": hashes["prospective_plan_sha256"]
            },
        },
        "publication_checkout_audit": _execution_git_audit(launch),
        "input_artifacts": score_inputs,
    }
    return launch, scoring


def _execution_git_audit(
    launch: dict[str, object],
) -> dict[str, object]:
    git = launch["git"]
    assert isinstance(git, dict)
    sha = str(git["sha"])
    return {
        "schema_version": "twirl_execution_checkout_audit_v1",
        "passed": True,
        "repo": "/synthetic/repo",
        "launch_declared_repo": str(git["repo"]),
        "observed_git_sha": sha,
        "operator_expected_git_sha": sha,
        "launch_git_sha": sha,
        "checkout_clean": True,
        "untracked_files_checked": True,
        "status_command": "git status --porcelain --untracked-files=all",
    }


def _queue_launch_binding() -> dict[str, object]:
    hashes = _queue_hashes()
    launch, scoring = _launch_and_scoring_summary(hashes)
    return validate_queue_launch_bindings(
        launch_manifest=launch,
        prospective_plan=_plan(),
        scoring_summary=scoring,
        artifact_hashes=hashes,
        git_audit=_execution_git_audit(launch),
    )


def _score_rows(primary: list[int], repeated: list[int]) -> pd.DataFrame:
    tics = np.asarray(primary + repeated, dtype=np.int64)
    rng = np.random.default_rng(630056)
    probabilities = rng.dirichlet(np.ones(4), len(tics))
    payload: dict[str, object] = {
            "score_schema_version": "twirl_teacher_v3_s63_candidate_scores_v1",
            "tic": tics,
            "sector": 63,
            "review_id": [f"candidate-{tic}" for tic in tics],
            "period_d": rng.uniform(0.2, 5.0, len(tics)),
            "t0_bjd": 2459000.0 + rng.uniform(0.0, 1.0, len(tics)),
            "duration_min": rng.uniform(3.0, 80.0, len(tics)),
            "depth": rng.uniform(0.001, 0.5, len(tics)),
            "depth_snr": rng.uniform(2.0, 50.0, len(tics)),
            "sde_max": rng.uniform(5.0, 100.0, len(tics)),
            "rep_aperture": ADP_ONLY_APERTURES[0],
            "rep_peak_rank": 1,
            "anchor_aperture": ADP_ONLY_APERTURES[0],
            "candidate_provenance_contract_version": CANDIDATE_CONTRACT_VERSION,
            "prospective_use_scope": PROSPECTIVE_USE_SCOPE,
            "science_ready": False,
            "promotion_enabled": False,
            "source_kind": "real_candidate",
            "p_planet_like": probabilities[:, 0],
            "p_eclipse_contact": probabilities[:, 1],
            "p_smooth_variable": probabilities[:, 2],
            "p_other": probabilities[:, 3],
            "p_preserve": rng.uniform(0.0, 1.0, len(tics)),
            "p_harmonic_P_over_2": rng.uniform(0.0, 1.0, len(tics)),
            "morphology_entropy": rng.uniform(0.0, np.log(4.0), len(tics)),
            "ensemble_disagreement": rng.uniform(0.0, 0.3, len(tics)),
            "model_version": MODEL_VERSION,
            "model_profile": TEACHER_V3_PRIMARY_PROFILE,
            "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
            "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
            "teacher_run_id": TEACHER_V3_RUN_ID,
            "teacher_training_sectors": "56-62",
            "scoring_sector": 63,
            "score_policy": S63_SCORE_POLICY,
        }
    for index, column in enumerate(TEACHER_V3_CANDIDATE_METADATA_COLUMNS, start=1):
        payload[column] = rng.normal(index, 0.01, len(tics))
    return pd.DataFrame(payload).assign(
        p_harmonic_P=lambda frame: 1.0 - frame["p_harmonic_P_over_2"]
    )


def test_plan_and_tic_inventories_are_byte_hash_bound(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(_plan(), sort_keys=True) + "\n", encoding="utf-8")
    plan_sha = file_sha256(plan_path)
    plan, observed = load_prospective_plan(plan_path, expected_sha256=plan_sha)
    assert plan["sector"] == 63
    assert observed == plan_sha

    tics_path = tmp_path / "tics.txt"
    tics_path.write_text("7\n9\n", encoding="ascii")
    assert read_tic_inventory(tics_path, expected_sha256=file_sha256(tics_path)) == (7, 9)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        read_tic_inventory(tics_path, expected_sha256="0" * 64)

    policy, policy_sha = load_selection_policy(
        SELECTION_POLICY_PATH,
        expected_sha256=file_sha256(SELECTION_POLICY_PATH),
    )
    assert policy["status"] == "frozen_before_s63_scores"
    assert policy["deterministic_bucket_selection"]["inclusion_probability"] is None
    assert policy["erratum_to_prospective_plan"]["prospective_plan_bytes_unchanged"] is True
    assert policy["blinding"]["cohort_annotation_withheld"] is True
    assert policy["blinding"]["cohort_joinable_from_visible_tic"] is True
    assert policy_sha == file_sha256(SELECTION_POLICY_PATH)
    tampered = dict(policy)
    tampered["score_data_opened_when_frozen"] = True
    tampered_path = tmp_path / "tampered_policy.json"
    tampered_path.write_text(json.dumps(tampered) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="score_data_opened_when_frozen"):
        load_selection_policy(
            tampered_path,
            expected_sha256=file_sha256(tampered_path),
        )
    missing_erratum = json.loads(json.dumps(policy))
    missing_erratum.pop("erratum_to_prospective_plan")
    missing_erratum_path = tmp_path / "missing_erratum_policy.json"
    missing_erratum_path.write_text(
        json.dumps(missing_erratum) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="prospective-plan erratum"):
        load_selection_policy(
            missing_erratum_path,
            expected_sha256=file_sha256(missing_erratum_path),
        )


def test_cohorts_exactly_partition_model_ready_and_reservation_complement() -> None:
    source_hashes = {
        "model_ready_allowlist_sha256": _hash("allowlist"),
        "teacher_v3_corpus_sha256": _hash("corpus"),
        "reserved_tics_sha256": _hash("reserved"),
    }
    primary, repeated, not_ready, summary = derive_s63_prospective_cohorts(
        model_ready_tics=[1, 2, 3, 4],
        frozen_corpus_tics=[2, 4, 9, 9],
        reserved_tics=[1, 2, 3, 4, 5],
        source_hashes=source_hashes,
    )
    assert primary["tic"].tolist() == [1, 3]
    assert repeated["tic"].tolist() == [2, 4]
    assert not_ready["tic"].tolist() == [5]
    assert set(primary["prospective_cohort"]) == {PRIMARY_COHORT}
    assert set(repeated["prospective_cohort"]) == {REPEATED_HOST_COHORT}
    assert summary["partition_checks"]["reservation_complement_exact"] is True
    assert summary["tic_inventory_sha256"]["model_ready"] == tic_inventory_sha256({1, 2, 3, 4})

    with pytest.raises(ValueError, match="outside the S63 reservation"):
        derive_s63_prospective_cohorts(
            model_ready_tics=[1, 6],
            frozen_corpus_tics=[1],
            reserved_tics=[1, 2],
            source_hashes=source_hashes,
        )


def test_cohort_writer_publishes_one_fresh_atomic_directory(tmp_path: Path) -> None:
    model = tmp_path / "model.txt"
    corpus = tmp_path / "corpus.csv"
    reserved = tmp_path / "reserved.txt"
    model.write_text("1\n2\n", encoding="ascii")
    pd.DataFrame({"tic": [2, 9]}).to_csv(corpus, index=False)
    reserved.write_text("1\n2\n3\n", encoding="ascii")
    source_hashes = {
        "model_ready_allowlist_sha256": file_sha256(model),
        "teacher_v3_corpus_sha256": file_sha256(corpus),
        "reserved_tics_sha256": file_sha256(reserved),
    }
    output = tmp_path / "cohorts"
    summary = write_s63_prospective_cohorts(
        model_ready_path=model,
        teacher_v3_corpus_path=corpus,
        reserved_tics_path=reserved,
        out_dir=output,
        source_hashes=source_hashes,
        producer_git_sha=_git_sha(),
        expected_reserved_count=3,
        expected_corpus_count=2,
    )
    assert (output / "cohort_summary.json").is_file()
    assert summary["outputs"]["primary_cohort"]["sha256"] == file_sha256(
        output / "primary_teacher_v3_disjoint_tics.txt"
    )
    assert summary["producer_git_sha"] == _git_sha()
    with pytest.raises(FileExistsError, match="fresh directory"):
        write_s63_prospective_cohorts(
            model_ready_path=model,
            teacher_v3_corpus_path=corpus,
            reserved_tics_path=reserved,
            out_dir=output,
            source_hashes=source_hashes,
            producer_git_sha=_git_sha(),
        )


def test_candidate_contract_uses_rank_one_small_with_primary_context(tmp_path: Path) -> None:
    candidates, summary = build_s63_rank_one_candidates(
        _candidate_peaks([101, 102, 103]),
        model_ready_tics=[101, 102, 103],
        artifact_hashes=_candidate_hashes(),
    )
    assert candidates["tic"].tolist() == [101, 102, 103]
    assert candidates["rep_peak_rank"].eq(1).all()
    assert candidates["rep_aperture"].eq(ADP_ONLY_APERTURES[0]).all()
    assert candidates["adp_period_d"].notna().all()
    assert candidates["science_ready"].eq(False).all()  # noqa: E712
    assert candidates["promotion_enabled"].eq(False).all()  # noqa: E712
    assert set(candidates["candidate_provenance_contract_version"]) == {
        CANDIDATE_CONTRACT_VERSION
    }
    assert set(TEACHER_V3_CANDIDATE_METADATA_COLUMNS).issubset(candidates.columns)
    for column in TEACHER_V3_CANDIDATE_METADATA_COLUMNS:
        assert pd.to_numeric(candidates[column], errors="coerce").notna().any()
    assert summary["prospective_use_scope"] == PROSPECTIVE_USE_SCOPE
    written = write_s63_rank_one_candidates(
        candidates,
        summary,
        out_path=tmp_path / "candidates.parquet",
        producer_git_sha=_git_sha(),
    )
    assert written["candidate_table_sha256"] == file_sha256(
        tmp_path / "candidates.parquet"
    )
    assert written["producer_git_sha"] == _git_sha()
    with pytest.raises(FileExistsError, match="candidate output must be fresh"):
        write_s63_rank_one_candidates(
            candidates,
            summary,
            out_path=tmp_path / "candidates.parquet",
            producer_git_sha=_git_sha(),
        )

    missing_primary = _candidate_peaks([101, 102, 103]).drop(index=3)
    with pytest.raises(ValueError, match="lack successful rank-one ADP-primary"):
        build_s63_rank_one_candidates(
            missing_primary,
            model_ready_tics=[101, 102, 103],
            artifact_hashes=_candidate_hashes(),
        )

    absent_metadata = _candidate_peaks([101, 102, 103])
    absent_metadata["adp_trend_ptp"] = np.nan
    with pytest.raises(ValueError, match="adp_trend_ptp is wholly absent"):
        build_s63_rank_one_candidates(
            absent_metadata,
            model_ready_tics=[101, 102, 103],
            artifact_hashes=_candidate_hashes(),
        )

    nonnumeric_metadata = _candidate_peaks([101, 102, 103])
    nonnumeric_metadata["adp_sml_own_even_depth"] = nonnumeric_metadata[
        "adp_sml_own_even_depth"
    ].astype(object)
    nonnumeric_metadata.loc[0, "adp_sml_own_even_depth"] = "not-a-number"
    with pytest.raises(ValueError, match="contains nonnumeric values"):
        build_s63_rank_one_candidates(
            nonnumeric_metadata,
            model_ready_tics=[101, 102, 103],
            artifact_hashes=_candidate_hashes(),
        )

    inconsistent_metadata = _candidate_peaks([101, 102, 103])
    inconsistent_metadata.loc[
        inconsistent_metadata["aperture"].eq(ADP_ONLY_APERTURES[1])
        & inconsistent_metadata["tic"].eq(101),
        "adp_sml_own_even_depth",
    ] = 999.0
    with pytest.raises(ValueError, match="aperture copies.*differ"):
        build_s63_rank_one_candidates(
            inconsistent_metadata,
            model_ready_tics=[101, 102, 103],
            artifact_hashes=_candidate_hashes(),
        )


def test_bls_summary_is_bound_without_tier1_or_science_claim() -> None:
    model_tics = [101, 102]
    model_inventory_sha = tic_inventory_sha256(set(model_tics))
    payload = {
        "passed": True,
        "sector": 63,
        "n_shards": 1,
        "shard_index": 0,
        "n_targets": 2,
        "n_targets_total": 2,
        "peak_table_sha256": _hash("peaks"),
        "target_allowlist_sha256": _hash("allowlist"),
        "target_allowlist_tics_sha256": model_inventory_sha,
        "orbitid_policy": "strict",
        "config": {
            "apertures": list(ADP_ONLY_APERTURES),
            "n_periods": 50_000,
            "n_peaks": 10,
            "source_product_tag": "A2v1",
        },
    }
    validated = validate_s63_bls_summary(
        payload,
        peak_table_sha256=_hash("peaks"),
        model_ready_allowlist_sha256=_hash("allowlist"),
        model_ready_tics_sha256=model_inventory_sha,
        n_model_ready_tics=2,
    )
    assert validated["science_ready"] is False
    assert validated["promotion_enabled"] is False
    assert "tier1" not in validated


def test_review_queue_is_deterministic_exact_and_blinded(tmp_path: Path) -> None:
    primary = list(range(10_000, 11_200))
    repeated = list(range(20_000, 20_140))
    scores = _score_rows(primary, repeated)
    queue, hidden, public_summary, hidden_summary = build_s63_prospective_review_queue(
        scores,
        primary_tics=primary,
        repeated_host_tics=repeated,
        artifact_hashes=_queue_hashes(),
        selection_policy=_selection_policy(),
        launch_binding=_queue_launch_binding(),
        producer_git_sha=_git_sha(),
    )
    queue_again, hidden_again, public_again, hidden_again_summary = (
        build_s63_prospective_review_queue(
        scores,
        primary_tics=primary,
        repeated_host_tics=repeated,
        artifact_hashes=_queue_hashes(),
        selection_policy=_selection_policy(),
        launch_binding=_queue_launch_binding(),
        producer_git_sha=_git_sha(),
        )
    )
    pd.testing.assert_frame_equal(queue, queue_again)
    pd.testing.assert_frame_equal(hidden, hidden_again)
    assert public_summary.keys() == public_again.keys()
    assert hidden_summary.keys() == hidden_again_summary.keys()
    assert len(queue) == 1100
    assert queue["tic"].nunique() == 1100
    assert tuple(queue.columns) == PUBLIC_REVIEW_QUEUE_COLUMNS
    assert public_summary["schema_version"] == (
        "twirl_teacher_v3_s63_public_review_queue_v1"
    )
    assert hidden_summary["queue_policy_version"] == QUEUE_POLICY_VERSION
    assert public_summary["producer_git_sha"] == _git_sha()
    assert hidden_summary["producer_git_sha"] == _git_sha()
    assert hidden_summary["cohort_annotation_withheld"] is True
    assert hidden_summary["cohort_joinable_from_visible_tic"] is True
    assert hidden_summary["cohort_analysis_deferred_until_review_complete"] is True
    assert "cohort_hidden_from_reviewer" not in hidden_summary
    assert "technically joinable" in public_summary["identity_joinability_notice"]
    assert hidden_summary["cohort_counts"] == {
        PRIMARY_COHORT: 1000,
        REPEATED_HOST_COHORT: 100,
    }
    assert hidden_summary["bucket_counts_by_cohort"][PRIMARY_COHORT] == {
        key: value for key, value in PRIMARY_QUOTAS.__dict__.items()
    }
    assert hidden_summary["bucket_counts_by_cohort"][REPEATED_HOST_COHORT] == {
        key: value for key, value in REPEATED_HOST_QUOTAS.__dict__.items()
    }
    for forbidden in (
        "prospective_cohort",
        "selection_bucket",
        "selection_score",
        "selection_rank",
        "p_planet_like",
        "ensemble_disagreement",
        "control_sde_stratum",
    ):
        assert forbidden not in queue
        assert forbidden in hidden
    assert hidden["teacher_scores_are_labels"].eq(False).all()  # noqa: E712
    assert hidden["pseudo_labels_authorized"].eq(False).all()  # noqa: E712
    assert hidden["execution_git_checkout_clean"].eq(True).all()  # noqa: E712
    assert hidden["execution_git_untracked_files_checked"].eq(True).all()  # noqa: E712
    assert hidden["execution_git_sha"].eq(_git_sha()).all()
    assert hidden["queue_order_hash"].is_monotonic_increasing
    assert "queue_order_hash" not in queue
    fractional_queue = queue.copy()
    fractional_hidden = hidden.copy()
    fractional_queue["row_id"] = fractional_queue["row_id"].astype(float)
    fractional_hidden["row_id"] = fractional_hidden["row_id"].astype(float)
    fractional_queue.loc[0, "row_id"] = 0.5
    fractional_hidden.loc[0, "row_id"] = 0.5
    with pytest.raises(ValueError, match="ordered range"):
        verify_s63_prospective_review_queue(
            fractional_queue,
            fractional_hidden,
            artifact_hashes=_queue_hashes(),
        )
    for cohort in (PRIMARY_COHORT, REPEATED_HOST_COHORT):
        control = hidden.loc[
            hidden["prospective_cohort"].eq(cohort)
            & hidden["selection_bucket"].eq("stratified_control")
        ]
        assert set(control["control_sde_stratum"].astype(int)) == {1, 2, 3, 4}
        assert np.allclose(
            control["selection_inclusion_probability"].astype(float),
            control["control_stratum_draw_count"].astype(float)
            / control["control_stratum_population"].astype(float),
        )
        deterministic = hidden.loc[
            hidden["prospective_cohort"].eq(cohort)
            & ~hidden["selection_bucket"].eq("stratified_control")
        ]
        assert deterministic["selection_inclusion_probability"].isna().all()
        assert deterministic["selection_weight"].isna().all()
        assert deterministic["deterministic_cut_fraction"].notna().all()
        source = scores.loc[
            scores["tic"].isin(primary if cohort == PRIMARY_COHORT else repeated)
        ].copy()
        expected_duty = (
            source["duration_min"] / (source["period_d"] * 1440.0)
        ).rank(method="average", pct=True)
        expected_by_tic = dict(zip(source["tic"], expected_duty, strict=True))
        observed = hidden.loc[
            hidden["prospective_cohort"].eq(cohort),
            ["tic", "broad_duty_cycle_percentile"],
        ]
        assert np.allclose(
            observed["broad_duty_cycle_percentile"],
            observed["tic"].map(expected_by_tic),
            rtol=0.0,
            atol=1e-15,
        )
    written = write_s63_prospective_review_queue(
        queue,
        hidden,
        public_summary,
        hidden_summary,
        public_dir=tmp_path / "public_queue",
        hidden_dir=tmp_path / "hidden_provenance",
    )
    assert written["passed"] is True
    assert written["completion"]["passed"] is True
    completion = verify_s63_prospective_review_bundle(
        public_dir=tmp_path / "public_queue",
        hidden_dir=tmp_path / "hidden_provenance",
    )
    assert completion["public_completion_marker_sha256"] == written["completion"][
        "public_completion_marker_sha256"
    ]
    public_marker = json.loads(
        (tmp_path / "public_queue" / "bundle_complete.json").read_text()
    )
    hidden_marker = json.loads(
        (tmp_path / "hidden_provenance" / "bundle_complete.json").read_text()
    )
    assert public_marker["bundle_id"] == hidden_marker["bundle_id"]
    assert set(public_marker["artifacts"]) == {
        "public_review_queue",
        "public_summary",
    }
    serialized_public_marker = json.dumps(public_marker, sort_keys=True)
    assert "hidden_selection_provenance.parquet" not in serialized_public_marker
    assert "hidden_summary.json" not in serialized_public_marker
    assert file_sha256(
        tmp_path / "hidden_provenance" / "hidden_selection_provenance.parquet"
    ) not in serialized_public_marker
    assert file_sha256(
        tmp_path / "hidden_provenance" / "hidden_summary.json"
    ) not in serialized_public_marker
    persisted_hidden = pd.read_parquet(
        tmp_path / "hidden_provenance" / "hidden_selection_provenance.parquet"
    )
    assert len(persisted_hidden) == 1100
    public_on_disk = json.loads(
        (tmp_path / "public_queue" / "public_summary.json").read_text()
    )
    public_keys: list[str] = []

    def collect_keys(value: object) -> None:
        if isinstance(value, dict):
            for key, nested in value.items():
                public_keys.append(str(key).lower())
                collect_keys(nested)
        elif isinstance(value, list):
            for nested in value:
                collect_keys(nested)

    collect_keys(public_on_disk)
    for forbidden in ("score", "bucket", "cohort", "hidden"):
        assert not any(forbidden in key for key in public_keys)
    assert "technically joinable" in public_on_disk["identity_joinability_notice"]
    assert not (tmp_path / "public_queue" / "hidden_selection_provenance.parquet").exists()
    hidden_on_disk = json.loads(
        (tmp_path / "hidden_provenance" / "hidden_summary.json").read_text()
    )
    assert hidden_on_disk["artifact_hashes"]["teacher_score_table_sha256"] == _hash(
        "scores"
    )
    assert stat.S_IMODE((tmp_path / "public_queue").stat().st_mode) == 0o700
    assert stat.S_IMODE(
        (tmp_path / "public_queue" / "review_queue_1100.csv").stat().st_mode
    ) == 0o600
    assert stat.S_IMODE(
        (tmp_path / "public_queue" / "public_summary.json").stat().st_mode
    ) == 0o600
    assert stat.S_IMODE((tmp_path / "hidden_provenance").stat().st_mode) == 0o700
    assert stat.S_IMODE(
        (
            tmp_path
            / "hidden_provenance"
            / "hidden_selection_provenance.parquet"
        ).stat().st_mode
    ) == 0o600
    assert stat.S_IMODE(
        (tmp_path / "hidden_provenance" / "hidden_summary.json").stat().st_mode
    ) == 0o600
    for directory in (tmp_path / "public_queue", tmp_path / "hidden_provenance"):
        assert stat.S_IMODE((directory / "bundle_complete.json").stat().st_mode) == 0o600

    leaked_queue = queue.copy()
    leaked_queue["member_0_morphology_logit"] = 0.5
    leaked_queue["member_0_p_planet_like"] = 0.5
    leaked_queue["predicted_morphology"] = "planet_like"
    leaked_queue["morphology_margin"] = 0.1
    with pytest.raises(ValueError, match="exact writer allowlist"):
        write_s63_prospective_review_queue(
            leaked_queue,
            hidden,
            public_summary,
            hidden_summary,
            public_dir=tmp_path / "leaky_public",
            hidden_dir=tmp_path / "leaky_hidden",
        )
    (tmp_path / "public_queue" / "bundle_complete.json").unlink()
    with pytest.raises(FileNotFoundError, match="completion marker"):
        verify_s63_prospective_review_bundle(
            public_dir=tmp_path / "public_queue",
            hidden_dir=tmp_path / "hidden_provenance",
        )


def test_review_queue_fails_closed_on_shortage_and_hash_drift() -> None:
    primary = list(range(10_000, 11_200))
    repeated = list(range(20_000, 20_099))
    scores = _score_rows(primary, repeated)
    with pytest.raises(ValueError, match="requires 100"):
        build_s63_prospective_review_queue(
            scores,
            primary_tics=primary,
            repeated_host_tics=repeated,
            artifact_hashes=_queue_hashes(),
            selection_policy=_selection_policy(),
            launch_binding=_queue_launch_binding(),
            producer_git_sha=_git_sha(),
        )

    repeated = list(range(20_000, 20_140))
    scores = _score_rows(primary, repeated)
    with pytest.raises(ValueError, match="producer_git_sha"):
        build_s63_prospective_review_queue(
            scores,
            primary_tics=primary,
            repeated_host_tics=repeated,
            artifact_hashes=_queue_hashes(),
            selection_policy=_selection_policy(),
            launch_binding=_queue_launch_binding(),
            producer_git_sha=_git_sha("wrong producer"),
        )
    queue, hidden, _, _ = build_s63_prospective_review_queue(
        scores,
        primary_tics=primary,
        repeated_host_tics=repeated,
        artifact_hashes=_queue_hashes(),
        selection_policy=_selection_policy(),
        launch_binding=_queue_launch_binding(),
        producer_git_sha=_git_sha(),
    )
    tampered = hidden.copy()
    tampered.loc[0, "candidate_table_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="inconsistent candidate_table_sha256"):
        verify_s63_prospective_review_queue(
            queue,
            tampered,
            artifact_hashes=_queue_hashes(),
        )
    reordered = queue.iloc[[1, 0, *range(2, len(queue))]].reset_index(drop=True)
    with pytest.raises(ValueError, match="ordered range 0..1099"):
        verify_s63_prospective_review_queue(
            reordered,
            hidden,
            artifact_hashes=_queue_hashes(),
        )


def test_queue_launch_binding_is_transitive_and_hash_exact() -> None:
    hashes = _queue_hashes()
    launch, scoring = _launch_and_scoring_summary(hashes)
    audit = validate_queue_launch_bindings(
        launch_manifest=launch,
        prospective_plan=_plan(),
        scoring_summary=scoring,
        artifact_hashes=hashes,
        git_audit=_execution_git_audit(launch),
    )
    assert audit["passed"] is True
    assert audit["score_bound_transitively_through_launch_authorization"] is True
    assert audit["launch_manifest_sha256"] == hashes["launch_manifest_sha256"]

    tampered = json.loads(json.dumps(launch))
    tampered["artifacts"]["selection_policy"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="selection_policy_sha256"):
        validate_queue_launch_bindings(
            launch_manifest=tampered,
            prospective_plan=_plan(),
            scoring_summary=scoring,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(tampered),
        )
    tampered_score = json.loads(json.dumps(scoring))
    tampered_score["output_sha256"]["scores"] = "0" * 64
    with pytest.raises(ValueError, match="score table"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=_plan(),
            scoring_summary=tampered_score,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(launch),
        )
    tampered_audit = _execution_git_audit(launch)
    tampered_audit["operator_expected_git_sha"] = _git_sha("drift")
    with pytest.raises(ValueError, match="operator_expected_git_sha"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=_plan(),
            scoring_summary=scoring,
            artifact_hashes=hashes,
            git_audit=tampered_audit,
        )
    failed_policy_check = json.loads(json.dumps(launch))
    failed_policy_check["checks"]["selection_policy"]["passed"] = False
    with pytest.raises(ValueError, match="selection_policy"):
        validate_queue_launch_bindings(
            launch_manifest=failed_policy_check,
            prospective_plan=_plan(),
            scoring_summary=scoring,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(failed_policy_check),
        )

    producer_drift = json.loads(json.dumps(launch))
    producer_drift["producer_git_sha"] = _git_sha("wrong launch producer")
    with pytest.raises(ValueError, match="producer_git_sha"):
        validate_queue_launch_bindings(
            launch_manifest=producer_drift,
            prospective_plan=_plan(),
            scoring_summary=scoring,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(producer_drift),
        )

    model_drift = json.loads(json.dumps(scoring))
    model_drift["execution"]["model_version"] = "wrong-model"
    with pytest.raises(ValueError, match="model_version"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=_plan(),
            scoring_summary=model_drift,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(launch),
        )

    checkpoint_drift = json.loads(json.dumps(scoring))
    checkpoint_drift["execution"]["checkpoint_sha256_by_fold"]["2"] = "0" * 64
    with pytest.raises(ValueError, match="checkpoint hashes"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=_plan(),
            scoring_summary=checkpoint_drift,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(launch),
        )

    publication_drift = json.loads(json.dumps(scoring))
    publication_drift["publication_checkout_audit"]["repo"] = "/other/checkout"
    with pytest.raises(ValueError, match="audits are not identical"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=_plan(),
            scoring_summary=publication_drift,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(launch),
        )

    policy_drift = json.loads(json.dumps(scoring))
    policy_drift["score_policy"] = "automatic-object-confirmation"
    with pytest.raises(ValueError, match="scoring policy"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=_plan(),
            scoring_summary=policy_drift,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(launch),
        )

    document_drift = json.loads(json.dumps(scoring))
    document_drift["release_verification"]["release_documents"][
        "release_summary"
    ]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="release_summary differs"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=_plan(),
            scoring_summary=document_drift,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(launch),
        )

    plan_drift = _plan()
    plan_drift["frozen_teacher_v3"]["run_id"] = "wrong-run"
    with pytest.raises(ValueError, match="run_id"):
        validate_queue_launch_bindings(
            launch_manifest=launch,
            prospective_plan=plan_drift,
            scoring_summary=scoring,
            artifact_hashes=hashes,
            git_audit=_execution_git_audit(launch),
        )


def test_split_publication_rejects_nesting_and_rolls_back(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    queue = pd.DataFrame(
        [{column: "" for column in PUBLIC_REVIEW_QUEUE_COLUMNS}],
        columns=PUBLIC_REVIEW_QUEUE_COLUMNS,
    )
    queue.loc[0, ["tic", "sector", "row_id"]] = [1, 63, 0]
    hidden = pd.DataFrame({"tic": [1], "selection_bucket": ["hidden"]})
    public_summary = {"schema_version": "reviewer_safe_test", "n_queue_rows": 1}
    hidden_summary = {"schema_version": "hidden_test", "artifact_hashes": _queue_hashes()}
    with pytest.raises(ValueError, match="must not be nested"):
        write_s63_prospective_review_queue(
            queue,
            hidden,
            public_summary,
            hidden_summary,
            public_dir=tmp_path / "nested",
            hidden_dir=tmp_path / "nested" / "private",
        )

    original_replace = prospective._replace_directory
    calls = 0

    def fail_second_replace(source: Path, target: Path) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("synthetic second-directory commit failure")
        original_replace(source, target)

    monkeypatch.setattr(prospective, "_replace_directory", fail_second_replace)
    public_dir = tmp_path / "public"
    hidden_dir = tmp_path / "private"
    with pytest.raises(OSError, match="synthetic second-directory"):
        write_s63_prospective_review_queue(
            queue,
            hidden,
            public_summary,
            hidden_summary,
            public_dir=public_dir,
            hidden_dir=hidden_dir,
        )
    assert not public_dir.exists()
    assert not hidden_dir.exists()
    assert not list(tmp_path.glob(".*.tmp-*"))
    assert not list(tmp_path.glob(".*.publish.claim"))
