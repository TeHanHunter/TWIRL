"""Prospective S63 Teacher-v3 enrichment and blinded-review contracts.

This module is deliberately separate from the S56 Tier-1 and active-learning
contracts.  It treats Teacher v3 only as an enrichment ranker for a frozen
human morphology review queue; it does not assert science readiness, create
pseudo-labels, or turn model probabilities into labels.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
from typing import Any, Mapping, Sequence
import uuid

import numpy as np
import pandas as pd

from twirl.vetting.adp_only import ADP_ONLY_APERTURES, classify_period_relation
from twirl.vetting.harmonic_cnn import MODEL_VERSION
from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION, native_group_path
from twirl.vetting.label_io import candidate_key
from twirl.vetting.teacher_v3 import TEACHER_V3_RUN_NAME
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_CHECKPOINT_NAMESPACE,
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
)


S63_SECTOR = 63
PROSPECTIVE_PLAN_SCHEMA = "twirl_teacher_v3_s63_prospective_plan_v1"
SELECTION_POLICY_SCHEMA = "twirl_teacher_v3_s63_selection_policy_v1"
COHORT_CONTRACT_VERSION = "twirl_teacher_v3_s63_cohorts_v1"
CANDIDATE_CONTRACT_VERSION = "twirl_teacher_v3_s63_rank1_candidates_v1"
QUEUE_POLICY_VERSION = "twirl_teacher_v3_s63_blinded_enrichment_v1"
PROSPECTIVE_USE_SCOPE = "human_morphology_enrichment_only"
S63_LAUNCH_CONTRACT = "twirl_teacher_v3_s63_prospective_launch_v1"
S63_AUTHORIZATION_CONTRACT = "s63_teacher_v3_prospective_v1"
S63_SCORE_SUMMARY_SCHEMA = "twirl_teacher_v3_s63_scoring_summary_v1"
S63_SCORE_TABLE_SCHEMA = "twirl_teacher_v3_s63_candidate_scores_v1"
S63_SCORE_POLICY = "prospective_human_review_enrichment_not_object_confirmation"
EXECUTION_CHECKOUT_AUDIT_SCHEMA = "twirl_execution_checkout_audit_v1"
QUEUE_LAUNCH_BINDING_SCHEMA = "twirl_teacher_v3_s63_queue_launch_binding_v1"
QUEUE_BUNDLE_COMPLETION_SCHEMA = (
    "twirl_teacher_v3_s63_queue_bundle_complete_v1"
)
QUEUE_BUNDLE_COMPLETION_FILENAME = "bundle_complete.json"

PRIMARY_COHORT = "teacher_v3_disjoint"
REPEATED_HOST_COHORT = "teacher_v3_repeated_host"
RESERVATION_EXCLUSION = "reserved_not_model_ready"


@dataclass(frozen=True)
class ProspectiveQuotas:
    """Locked six-bucket quota for one prospective cohort."""

    planet_preserve: int
    eclipse_contact: int
    smooth_variable: int
    broad_dip: int
    disagreement_harmonic: int
    stratified_control: int

    @property
    def total(self) -> int:
        return int(sum(asdict(self).values()))


PRIMARY_QUOTAS = ProspectiveQuotas(300, 200, 150, 100, 150, 100)
REPEATED_HOST_QUOTAS = ProspectiveQuotas(30, 20, 15, 10, 15, 10)

REQUIRED_ARTIFACT_HASHES: tuple[str, ...] = (
    "prospective_plan_sha256",
    "selection_policy_sha256",
    "launch_manifest_sha256",
    "reserved_tics_sha256",
    "teacher_v3_corpus_sha256",
    "model_ready_allowlist_sha256",
    "primary_cohort_sha256",
    "repeated_host_cohort_sha256",
    "cohort_summary_sha256",
    "candidate_table_sha256",
    "teacher_score_table_sha256",
    "teacher_score_summary_sha256",
)

REVIEWER_COLUMNS: tuple[str, ...] = (
    "row_id",
    "review_id",
    "candidate_key",
    "tic",
    "sector",
    "cam",
    "ccd",
    "tmag",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "rep_aperture",
    "rep_peak_rank",
    "anchor_aperture",
    "n_apertures_agree",
    "apertures_agree",
    "aperture_period_relation",
    "aperture_period_rel_delta",
    "aperture_depth_ratio_primary_over_small",
    "aperture_disagreement_flag",
    "twirl_vet_sheet_name",
    "twirl_vet_sheet_pdf_name",
)
PUBLIC_REVIEW_QUEUE_COLUMNS: tuple[str, ...] = (
    *REVIEWER_COLUMNS,
    "label",
    "label_source",
    "labeler",
    "notes",
    "updated_utc",
)

TEACHER_V3_APERTURE_METADATA_FIELDS: tuple[str, ...] = (
    "own_even_depth",
    "own_odd_depth",
    "own_even_odd_depth_delta",
    "own_even_odd_sigma_delta",
    "trend_ptp",
)
TEACHER_V3_CANDIDATE_METADATA_COLUMNS: tuple[str, ...] = tuple(
    f"{prefix}_{field}"
    for prefix in ("adp_sml", "adp")
    for field in TEACHER_V3_APERTURE_METADATA_FIELDS
)
S63_TARGET_METADATA_CONTRACT = (
    "twirl_teacher_v3_s63_rank1_candidate_metadata_v1"
)

_TIC_PATTERN = re.compile(r"[1-9][0-9]*")
_GIT_SHA_PATTERN = re.compile(r"[0-9a-f]{40}")
_MAX_INT64 = int(np.iinfo(np.int64).max)


def file_sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Return the byte-exact SHA-256 digest of one artifact."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_sha256(value: Any, *, context: str) -> str:
    """Normalize and validate a SHA-256 binding."""

    normalized = str(value).strip().lower()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise ValueError(f"{context} must be a valid SHA-256 digest")
    return normalized


def validate_git_sha(value: Any, *, context: str) -> str:
    """Normalize and validate one full lowercase Git object identity."""

    normalized = str(value).strip()
    if _GIT_SHA_PATTERN.fullmatch(normalized) is None:
        raise ValueError(f"{context} must be a full lowercase 40-character Git SHA")
    return normalized


def tic_inventory_sha256(tics: Sequence[int] | set[int]) -> str:
    """Hash a sorted, newline-delimited unique TIC inventory."""

    normalized = _normalize_tic_values(tics, context="TIC inventory")
    if len(normalized) != len(tics):
        raise ValueError("TIC inventory contains duplicate TICs")
    payload = "".join(f"{tic}\n" for tic in normalized).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _normalize_tic_values(values: Sequence[Any], *, context: str) -> list[int]:
    output: list[int] = []
    for value in values:
        if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
            tic = int(value)
        else:
            text = str(value).strip()
            if not _TIC_PATTERN.fullmatch(text):
                raise ValueError(f"{context} contains an invalid TIC value: {value!r}")
            tic = int(text)
        if tic <= 0 or tic > _MAX_INT64:
            raise ValueError(f"{context} contains an out-of-range TIC: {tic}")
        output.append(tic)
    return sorted(output)


def _read_table(path: Path) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"unsupported table format: {path}")


def _hash_bound_read(path: Path, expected_sha256: str, reader: Any) -> Any:
    path = Path(path)
    expected = validate_sha256(expected_sha256, context=f"expected hash for {path}")
    before = file_sha256(path)
    if before != expected:
        raise ValueError(f"SHA-256 mismatch for {path}: expected {expected}, observed {before}")
    value = reader(path)
    after = file_sha256(path)
    if after != before:
        raise RuntimeError(f"artifact changed while it was read: {path}")
    return value


def read_tic_inventory(
    path: Path,
    *,
    expected_sha256: str,
    allow_duplicate_rows: bool = False,
) -> tuple[int, ...]:
    """Read a hash-bound TXT/CSV/Parquet TIC inventory."""

    def reader(source: Path) -> tuple[int, ...]:
        if source.suffix.lower() == ".txt":
            raw = source.read_text(encoding="utf-8").splitlines()
            values = [line.strip() for line in raw if line.strip()]
        else:
            frame = _read_table(source)
            tic_column = "tic" if "tic" in frame else "tic_id" if "tic_id" in frame else None
            if tic_column is None:
                raise KeyError(f"TIC table is missing tic/tic_id: {source}")
            values = frame[tic_column].tolist()
        normalized = _normalize_tic_values(values, context=str(source))
        if not allow_duplicate_rows and len(normalized) != len(set(normalized)):
            raise ValueError(f"TIC inventory contains duplicates: {source}")
        return tuple(sorted(set(normalized)))

    return _hash_bound_read(path, expected_sha256, reader)


def load_prospective_plan(
    path: Path,
    *,
    expected_sha256: str,
) -> tuple[dict[str, Any], str]:
    """Read and validate the immutable S63 prospective plan."""

    payload = _hash_bound_read(
        path,
        expected_sha256,
        lambda source: json.loads(source.read_text(encoding="utf-8")),
    )
    if not isinstance(payload, dict):
        raise ValueError("prospective plan must be a JSON object")
    if payload.get("schema_version") != PROSPECTIVE_PLAN_SCHEMA:
        raise ValueError("prospective plan schema is incompatible")
    if payload.get("sector") != S63_SECTOR:
        raise ValueError("prospective plan must be bounded to Sector 63")
    search = payload.get("search_contract")
    if not isinstance(search, Mapping):
        raise ValueError("prospective plan is missing search_contract")
    expected_search = {
        "apertures": list(ADP_ONLY_APERTURES),
        "anchor_aperture": ADP_ONLY_APERTURES[0],
        "context_aperture": ADP_ONLY_APERTURES[1],
        "n_periods": 50_000,
        "n_retained_peaks_per_aperture": 10,
        "candidate_peak_rank": 1,
        "orbitid_policy": "strict",
        "absolute_probability_threshold": None,
    }
    for key, expected in expected_search.items():
        if search.get(key) != expected:
            raise ValueError(f"prospective plan search_contract has incompatible {key}")
    queue = payload.get("blinded_queue")
    if not isinstance(queue, Mapping) or queue.get("selection_seed") != 630056:
        raise ValueError("prospective plan has an incompatible blinded queue seed")
    expected_quotas = (
        ("primary", PRIMARY_QUOTAS),
        ("repeated_host_side_cohort", REPEATED_HOST_QUOTAS),
    )
    for name, quotas in expected_quotas:
        entry = queue.get(name)
        if not isinstance(entry, Mapping):
            raise ValueError(f"prospective plan is missing {name} queue")
        if entry.get("n_tics") != quotas.total or entry.get("quotas") != asdict(quotas):
            raise ValueError(f"prospective plan has incompatible {name} quotas")
    teacher = payload.get("frozen_teacher_v3")
    if not isinstance(teacher, Mapping):
        raise ValueError("prospective plan is missing frozen Teacher-v3 identity")
    if teacher.get("automatic_promotion") is not False:
        raise ValueError("prospective plan must disable automatic promotion")
    if teacher.get("student_training_authorized") is not False:
        raise ValueError("prospective plan must not authorize student training")
    training = payload.get("frozen_training_identity")
    if not isinstance(training, Mapping):
        raise ValueError("prospective plan is missing frozen training identity")
    validate_sha256(
        training.get("morphology_corpus_sha256"),
        context="prospective plan morphology corpus",
    )
    if type(training.get("n_corpus_tics")) is not int or training["n_corpus_tics"] <= 0:
        raise ValueError("prospective plan has an invalid frozen corpus TIC count")
    reservation = payload.get("s63_identity_reservation")
    if not isinstance(reservation, Mapping):
        raise ValueError("prospective plan is missing the S63 identity reservation")
    validate_sha256(
        reservation.get("reserved_tics_sha256"),
        context="prospective plan S63 reservation",
    )
    if (
        type(reservation.get("n_requested_tics")) is not int
        or reservation["n_requested_tics"] <= 0
    ):
        raise ValueError("prospective plan has an invalid reservation TIC count")
    return payload, validate_sha256(expected_sha256, context="prospective plan")


def _mapping(value: Any, *, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{context} must be a JSON object")
    return value


def validate_selection_policy(policy: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the score-independent queue policy against executable choices."""

    policy = _mapping(policy, context="S63 selection policy")
    exact_top = {
        "schema_version": SELECTION_POLICY_SCHEMA,
        "status": "frozen_before_s63_scores",
        "sector": S63_SECTOR,
        "score_data_opened_when_frozen": False,
    }
    for key, expected in exact_top.items():
        if policy.get(key) != expected:
            raise ValueError(f"S63 selection policy has incompatible {key}")

    cohort = _mapping(
        policy.get("cohort_evaluation"), context="selection policy cohort_evaluation"
    )
    if cohort.get("independence_rule") != (
        "Compute every percentile and selection score independently within each "
        "complete cohort before any bucket depletion"
    ):
        raise ValueError("selection policy does not freeze independent cohort scoring")

    formulas = _mapping(policy.get("score_formulas"), context="selection policy formulas")
    expected_formulas = {
        "planet_preserve": (
            "0.55*p_planet_like + 0.35*p_preserve + "
            "0.10*I(period_d >= 1.0 and sde_max <= 40.0)"
        ),
        "eclipse_contact": "p_eclipse_contact",
        "smooth_variable": "p_smooth_variable",
        "broad_dip": (
            "p_preserve * broad_duty_cycle_percentile * "
            "(1 - 0.5*p_eclipse_contact) * (1 - 0.5*p_smooth_variable)"
        ),
        "disagreement_harmonic": (
            "0.45*morphology_entropy + 0.35*ensemble_disagreement + "
            "0.20*max(non-P p_harmonic_* probabilities)"
        ),
    }
    for name, expected in expected_formulas.items():
        entry = _mapping(formulas.get(name), context=f"selection formula {name}")
        if entry.get("formula") != expected:
            raise ValueError(f"selection policy formula drifted for {name}")
    expected_weights = {
        "planet_preserve": {
            "p_planet_like": 0.55,
            "p_preserve": 0.35,
            "long_period_low_sde_indicator": 0.1,
        },
        "eclipse_contact": {"p_eclipse_contact": 1.0},
        "smooth_variable": {"p_smooth_variable": 1.0},
        "disagreement_harmonic": {
            "morphology_entropy": 0.45,
            "ensemble_disagreement": 0.35,
            "maximum_non_P_harmonic_probability": 0.2,
        },
    }
    for name, expected in expected_weights.items():
        if formulas[name].get("weights") != expected:
            raise ValueError(f"selection policy weights drifted for {name}")
    if formulas["broad_dip"].get("duty_cycle") != (
        "duration_min / (period_d * 1440.0)"
    ):
        raise ValueError("selection policy broad-duty definition drifted")

    duty = _mapping(
        policy.get("broad_duty_cycle_percentile"),
        context="selection policy duty percentile",
    )
    if duty.get("cohorts_processed_separately") is not True or duty.get(
        "algorithm"
    ) != "pandas Series.rank(method='average', pct=True) on duty_cycle":
        raise ValueError("selection policy broad-duty percentile algorithm drifted")
    if duty.get("population") != (
        "all validated unique-TIC candidates in the current cohort before any "
        "bucket depletion"
    ):
        raise ValueError("selection policy broad-duty population drifted")

    deterministic = _mapping(
        policy.get("deterministic_bucket_selection"),
        context="selection policy deterministic buckets",
    )
    expected_order = [
        "planet_preserve",
        "eclipse_contact",
        "smooth_variable",
        "broad_dip",
        "disagreement_harmonic",
    ]
    if deterministic.get("depletion_order") != expected_order:
        raise ValueError("selection policy bucket-depletion order drifted")
    if deterministic.get("ordering_within_each_remaining_pool") != [
        "bucket selection score descending",
        "sde_max descending",
        "tic ascending",
    ]:
        raise ValueError("selection policy deterministic tie-breaks drifted")
    for key in ("probability_threshold", "inclusion_probability", "selection_weight"):
        if key not in deterministic or deterministic[key] is not None:
            raise ValueError(f"deterministic selection {key} must be explicitly null")

    control = _mapping(
        policy.get("control_selection"), context="selection policy control"
    )
    exact_control = {
        "selection_seed": 630056,
        "pre_qcut_order": "tic ascending",
        "tie_rule": (
            "rank(method='first') so equal SDE values are ordered by ascending TIC"
        ),
        "qcut_algorithm": (
            "pandas qcut of the first-ranked SDE values into exactly four integer "
            "quartiles labeled 1,2,3,4; duplicate edges or fewer than four strata "
            "fail closed"
        ),
        "draw_allocation": (
            "floor(control_quota/4) per quartile, with one additional draw assigned "
            "in quartile order 1,2,... for the control_quota modulo 4 remainder"
        ),
        "within_stratum_order": (
            "ascending lowercase SHA-256 hex digest of the ASCII string "
            "'{seed}|{cohort}|control|q{stratum}|{tic}'"
        ),
        "sampling_without_replacement": True,
        "inclusion_probability": (
            "realized stratum draw count divided by stratum population"
        ),
        "selection_weight": (
            "stratum population divided by realized stratum draw count"
        ),
    }
    for key, expected in exact_control.items():
        if control.get(key) != expected:
            raise ValueError(f"selection policy control {key} drifted")

    combined = _mapping(
        policy.get("combined_queue_order"), context="selection policy queue order"
    )
    if combined.get("algorithm") != (
        "ascending lowercase SHA-256 hex digest of the ASCII string "
        "'{seed}|combined-review-order|{tic}'"
    ):
        raise ValueError("selection policy combined queue order drifted")
    if combined.get("explicitly_annotates_cohort_or_bucket") is not False:
        raise ValueError("selection policy queue order must not annotate selection")
    if combined.get("scope") != (
        "The seeded order itself does not encode an explicit cohort or bucket "
        "annotation; visible TIC identity remains externally joinable"
    ):
        raise ValueError("selection policy queue-order joinability scope drifted")

    quotas = _mapping(policy.get("quotas"), context="selection policy quotas")
    expected_quota_payloads = {
        "primary": {**asdict(PRIMARY_QUOTAS), "total": PRIMARY_QUOTAS.total},
        "repeated_host": {
            **asdict(REPEATED_HOST_QUOTAS),
            "total": REPEATED_HOST_QUOTAS.total,
        },
    }
    if dict(quotas) != expected_quota_payloads:
        raise ValueError("selection policy quotas drifted")
    if not str(policy.get("shortage_policy", "")).startswith("Fail closed"):
        raise ValueError("selection policy does not fail closed on shortages")
    erratum = _mapping(
        policy.get("erratum_to_prospective_plan"),
        context="selection policy prospective-plan erratum",
    )
    expected_erratum = {
        "affected_field": (
            "blinded_queue.reviewer_visible_fields_exclude entry 'cohort membership'"
        ),
        "interpretation": (
            "The original wording means explicit cohort annotation is withheld; "
            "TIC and candidate identity remain visible and technically joinable to "
            "the frozen Teacher-v3 corpus"
        ),
        "prospective_plan_bytes_unchanged": True,
        "s63_teacher_scores_opened_before_erratum": False,
    }
    if dict(erratum) != expected_erratum:
        raise ValueError("selection policy prospective-plan erratum drifted")
    blinding = _mapping(
        policy.get("blinding"), context="selection policy blinding"
    )
    expected_exclusions = [
        "teacher probabilities and logits",
        "ensemble disagreement and entropy",
        "all selection scores and ranks",
        "selection bucket",
        "control stratum",
        "explicit cohort annotation",
        "all hidden artifact hashes",
    ]
    if blinding.get("reviewer_visible_excludes") != expected_exclusions:
        raise ValueError("selection policy reviewer exclusions drifted")
    if blinding.get("reviewer_visible_identity") != [
        "tic",
        "candidate_key",
        "vet-sheet identity",
    ]:
        raise ValueError("selection policy visible identity drifted")
    for key in (
        "cohort_annotation_withheld",
        "cohort_joinable_from_visible_tic",
        "cohort_analysis_deferred_until_review_complete",
    ):
        if blinding.get(key) is not True:
            raise ValueError(f"selection policy blinding {key} must be true")
    if blinding.get("review_completion") != (
        "All 1100 frozen rows require an accepted human morphology decision before "
        "one-time unblinding and cohort-wise analysis"
    ):
        raise ValueError("selection policy review-completion boundary drifted")
    return dict(policy)


def load_selection_policy(
    path: Path,
    *,
    expected_sha256: str,
) -> tuple[dict[str, Any], str]:
    """Read the exact pre-score selection policy and validate its semantics."""

    payload = _hash_bound_read(
        path,
        expected_sha256,
        lambda source: json.loads(source.read_text(encoding="utf-8")),
    )
    validated = validate_selection_policy(payload)
    return validated, validate_sha256(expected_sha256, context="selection policy")


def _validate_execution_checkout_audit(
    payload: Mapping[str, Any],
    *,
    launch_git: Mapping[str, Any],
    context: str,
) -> dict[str, Any]:
    audit = _mapping(payload, context=context)
    expected_fields = {
        "schema_version": EXECUTION_CHECKOUT_AUDIT_SCHEMA,
        "passed": True,
        "checkout_clean": True,
        "untracked_files_checked": True,
        "status_command": "git status --porcelain --untracked-files=all",
    }
    for key, expected in expected_fields.items():
        if audit.get(key) != expected:
            raise ValueError(f"{context} has incompatible {key}")
    if not str(audit.get("repo", "")).strip():
        raise ValueError(f"{context} has no repository path")
    if audit.get("launch_declared_repo") != str(launch_git.get("repo", "")):
        raise ValueError(f"{context} launch repository drifted")
    launch_git_sha = str(launch_git.get("sha", "")).strip()
    for key in (
        "observed_git_sha",
        "operator_expected_git_sha",
        "launch_git_sha",
    ):
        observed = str(audit.get(key, "")).strip()
        if observed != launch_git_sha:
            raise ValueError(f"{context} {key} differs from launch git.sha")
    return dict(audit)


def _validate_frozen_teacher_v3_identity(
    prospective_plan: Mapping[str, Any],
) -> dict[str, Any]:
    plan = _mapping(prospective_plan, context="S63 prospective plan")
    if plan.get("schema_version") != PROSPECTIVE_PLAN_SCHEMA:
        raise ValueError("queue prospective plan schema drifted")
    if plan.get("sector") != S63_SECTOR:
        raise ValueError("queue prospective plan has the wrong sector")
    teacher = _mapping(
        plan.get("frozen_teacher_v3"), context="frozen Teacher-v3 identity"
    )
    expected_identity = {
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
        "automatic_promotion": False,
        "student_training_authorized": False,
    }
    for key, expected in expected_identity.items():
        if teacher.get(key) != expected:
            raise ValueError(f"frozen Teacher-v3 identity has incompatible {key}")
    document_fields = (
        "summary_sha256",
        "pretest_model_freeze_sha256",
        "selected_checkpoint_manifest_sha256",
        "pooled_oof_calibration_sha256",
    )
    documents = {
        key: validate_sha256(
            teacher.get(key), context=f"frozen Teacher-v3 {key}"
        )
        for key in document_fields
    }
    try:
        temperature = float(teacher["temperature"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("frozen Teacher-v3 temperature is invalid") from exc
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("frozen Teacher-v3 temperature must be positive")
    records = teacher.get("checkpoints")
    if not isinstance(records, list) or len(records) != 5:
        raise ValueError("frozen Teacher-v3 identity must bind five checkpoints")
    checkpoints: dict[str, str] = {}
    for expected_fold, record in enumerate(
        sorted(records, key=lambda value: int(value.get("fold", -1)))
    ):
        record = _mapping(
            record, context=f"frozen Teacher-v3 checkpoint {expected_fold}"
        )
        if record.get("fold") != expected_fold:
            raise ValueError("frozen Teacher-v3 checkpoint folds must be 0..4")
        checkpoints[str(expected_fold)] = validate_sha256(
            record.get("sha256"),
            context=f"frozen Teacher-v3 checkpoint {expected_fold}",
        )
    return {
        **expected_identity,
        "temperature": temperature,
        "documents": documents,
        "checkpoint_sha256_by_fold": checkpoints,
    }


def _validate_teacher_v3_scoring_release(
    scoring: Mapping[str, Any],
    *,
    frozen_teacher: Mapping[str, Any],
) -> dict[str, Any]:
    if scoring.get("score_policy") != S63_SCORE_POLICY:
        raise ValueError("Teacher-v3 S63 scoring policy drifted")
    execution = _mapping(scoring.get("execution"), context="scoring execution")
    expected_execution = {
        "schema_version": "twirl_teacher_v3_s63_ensemble_execution_v1",
        "run_id": frozen_teacher["run_id"],
        "release_name": frozen_teacher["release_name"],
        "model_version": frozen_teacher["model_version"],
        "profile": frozen_teacher["selected_profile"],
        "checkpoint_namespace": frozen_teacher["checkpoint_namespace"],
        "input_contract_version": frozen_teacher["input_contract_version"],
        "sector": S63_SECTOR,
    }
    for key, expected in expected_execution.items():
        if execution.get(key) != expected:
            raise ValueError(f"Teacher-v3 scoring execution has incompatible {key}")
    checkpoint_hashes = {
        str(key): validate_sha256(
            value, context=f"scoring execution checkpoint {key}"
        )
        for key, value in _mapping(
            execution.get("checkpoint_sha256_by_fold"),
            context="scoring execution checkpoints",
        ).items()
    }
    if checkpoint_hashes != frozen_teacher["checkpoint_sha256_by_fold"]:
        raise ValueError("scoring execution checkpoint hashes differ from the freeze")
    release = _mapping(
        scoring.get("release_verification"),
        context="scoring verified Teacher-v3 release",
    )
    if release.get("schema_version") != "twirl_teacher_v3_verified_release_v1":
        raise ValueError("verified Teacher-v3 release schema drifted")
    release_documents = _mapping(
        release.get("release_documents"), context="verified release documents"
    )
    expected_documents = {
        "release_summary": "summary_sha256",
        "pretest_freeze": "pretest_model_freeze_sha256",
        "checkpoint_manifest": "selected_checkpoint_manifest_sha256",
        "calibration": "pooled_oof_calibration_sha256",
    }
    for release_name, plan_field in expected_documents.items():
        record = _mapping(
            release_documents.get(release_name),
            context=f"verified release document {release_name}",
        )
        if validate_sha256(
            record.get("sha256"), context=f"verified {release_name}"
        ) != frozen_teacher["documents"][plan_field]:
            raise ValueError(
                f"verified Teacher-v3 {release_name} differs from the freeze"
            )
    provenance = _mapping(
        release.get("input_provenance"), context="verified release input provenance"
    )
    for key, expected in {
        "checkpoint_namespace": frozen_teacher["checkpoint_namespace"],
        "input_contract_version": frozen_teacher["input_contract_version"],
    }.items():
        if provenance.get(key) != expected:
            raise ValueError(f"verified Teacher-v3 release has incompatible {key}")
    release_checkpoints = {
        str(key): validate_sha256(
            value, context=f"verified release checkpoint {key}"
        )
        for key, value in _mapping(
            release.get("checkpoint_sha256_by_fold"),
            context="verified release checkpoints",
        ).items()
    }
    if release_checkpoints != frozen_teacher["checkpoint_sha256_by_fold"]:
        raise ValueError("verified release checkpoint hashes differ from the freeze")
    try:
        release_temperature = float(release["calibration_temperature"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("verified release calibration temperature is invalid") from exc
    if not np.isclose(
        release_temperature,
        float(frozen_teacher["temperature"]),
        rtol=0.0,
        atol=1.0e-15,
    ):
        raise ValueError("verified release temperature differs from the freeze")
    manifest_audit = _mapping(
        release.get("manifest_audit"), context="verified release manifest audit"
    )
    if manifest_audit.get("n_verified_checkpoints") != 5:
        raise ValueError("verified release manifest did not prove five checkpoints")
    manifest_checkpoints = {
        str(key): validate_sha256(
            value, context=f"manifest-audit checkpoint {key}"
        )
        for key, value in _mapping(
            manifest_audit.get("checkpoint_sha256_by_fold"),
            context="verified manifest checkpoints",
        ).items()
    }
    if manifest_checkpoints != frozen_teacher["checkpoint_sha256_by_fold"]:
        raise ValueError("manifest-audit checkpoint hashes differ from the freeze")
    return {
        "score_policy": S63_SCORE_POLICY,
        "identity": dict(expected_execution),
        "release_document_sha256": {
            name: frozen_teacher["documents"][field]
            for name, field in expected_documents.items()
        },
        "checkpoint_sha256_by_fold": dict(checkpoint_hashes),
    }


def validate_queue_launch_bindings(
    *,
    launch_manifest: Mapping[str, Any],
    prospective_plan: Mapping[str, Any],
    scoring_summary: Mapping[str, Any],
    artifact_hashes: Mapping[str, str],
    git_audit: Mapping[str, Any],
) -> dict[str, Any]:
    """Prove queue inputs descend from one authorized preprocessing launch."""

    launch = _mapping(launch_manifest, context="S63 launch manifest")
    expected_launch_fields = {
        "contract_version": S63_LAUNCH_CONTRACT,
        "authorization_contract": S63_AUTHORIZATION_CONTRACT,
        "sector": S63_SECTOR,
        "orbits": [133, 134],
        "status": "preprocessing_complete_scoring_not_started",
        "score_or_queue_read": False,
        "passed": True,
    }
    for key, expected in expected_launch_fields.items():
        if launch.get(key) != expected:
            raise ValueError(f"S63 queue launch manifest has incompatible {key}")
    frozen_teacher = _validate_frozen_teacher_v3_identity(prospective_plan)
    launch_git = _mapping(launch.get("git"), context="launch git identity")
    if launch_git.get("checkout_clean") is not True:
        raise ValueError("S63 queue launch git.checkout_clean must be true")
    if launch_git.get("untracked_files_checked") is not True:
        raise ValueError("S63 queue launch git.untracked_files_checked must be true")
    launch_git_sha = str(launch_git.get("sha", "")).strip()
    if _GIT_SHA_PATTERN.fullmatch(launch_git_sha) is None:
        raise ValueError("S63 queue launch git.sha must be a full lowercase Git SHA")
    launch_producer_git_sha = validate_git_sha(
        launch.get("producer_git_sha"), context="launch producer_git_sha"
    )
    if launch_producer_git_sha != launch_git_sha:
        raise ValueError("launch producer_git_sha differs from launch git.sha")
    execution_checkout = _validate_execution_checkout_audit(
        git_audit,
        launch_git=launch_git,
        context="queue execution checkout audit",
    )
    hashes = {
        key: validate_sha256(value, context=key)
        for key, value in artifact_hashes.items()
    }
    missing = sorted(set(REQUIRED_ARTIFACT_HASHES) - set(hashes))
    if missing:
        raise KeyError(f"queue launch bindings are missing hashes: {missing}")
    artifacts = _mapping(launch.get("artifacts"), context="launch artifacts")
    launch_bindings = {
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
    verified_launch_artifacts: dict[str, str] = {}
    for artifact_name, hash_name in launch_bindings.items():
        record = _mapping(
            artifacts.get(artifact_name),
            context=f"launch artifacts.{artifact_name}",
        )
        observed = validate_sha256(
            record.get("sha256"), context=f"launch artifacts.{artifact_name}"
        )
        if observed != hashes[hash_name]:
            raise ValueError(
                f"launch artifacts.{artifact_name} does not match {hash_name}"
            )
        if type(record.get("size_bytes")) is not int or record["size_bytes"] <= 0:
            raise ValueError(f"launch artifacts.{artifact_name} has invalid size")
        if not str(record.get("path", "")).strip():
            raise ValueError(f"launch artifacts.{artifact_name} has no path")
        verified_launch_artifacts[artifact_name] = observed
    checks = _mapping(launch.get("checks"), context="launch checks")
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
        check = _mapping(checks.get(name), context=f"launch checks.{name}")
        if check.get("passed") is not True:
            raise ValueError(f"launch preprocessing check {name} did not pass")

    scoring = _mapping(scoring_summary, context="Teacher-v3 S63 scoring summary")
    if scoring.get("schema_version") != S63_SCORE_SUMMARY_SCHEMA:
        raise ValueError("Teacher-v3 S63 scoring summary schema drifted")
    if scoring.get("sector") != S63_SECTOR:
        raise ValueError("Teacher-v3 S63 scoring summary has the wrong sector")
    if scoring.get("strict_provenance_passed") is not True:
        raise ValueError("Teacher-v3 S63 scoring provenance did not pass")
    teacher_release_binding = _validate_teacher_v3_scoring_release(
        scoring,
        frozen_teacher=frozen_teacher,
    )
    score_outputs = _mapping(
        scoring.get("output_sha256"), context="scoring output hashes"
    )
    observed_score_sha = validate_sha256(
        score_outputs.get("scores"), context="scoring output score table"
    )
    if observed_score_sha != hashes["teacher_score_table_sha256"]:
        raise ValueError("scoring summary does not bind the queue score table")
    authorization = _mapping(
        scoring.get("launch_authorization"), context="scoring launch authorization"
    )
    if authorization.get("passed") is not True:
        raise ValueError("scoring launch authorization did not pass")
    scoring_execution_checkout = _validate_execution_checkout_audit(
        authorization.get("execution_checkout"),
        launch_git=launch_git,
        context="scoring execution checkout audit",
    )
    publication_checkout = _validate_execution_checkout_audit(
        scoring.get("publication_checkout_audit"),
        launch_git=launch_git,
        context="scoring publication checkout audit",
    )
    if not (
        scoring_execution_checkout == publication_checkout == execution_checkout
    ):
        raise ValueError(
            "scoring launch, publication, and queue checkout audits are not identical"
        )
    authorized_launch = _mapping(
        authorization.get("launch_manifest"),
        context="scoring authorized launch manifest",
    )
    if validate_sha256(
        authorized_launch.get("sha256"), context="scoring authorized launch"
    ) != hashes["launch_manifest_sha256"]:
        raise ValueError("score table was not produced from the queue launch manifest")
    authorized_plan = _mapping(
        authorization.get("preregistered_contract"),
        context="scoring authorized prospective plan",
    )
    if validate_sha256(
        authorized_plan.get("sha256"), context="scoring authorized plan"
    ) != hashes["prospective_plan_sha256"]:
        raise ValueError("score table was not produced from the frozen prospective plan")
    authorized_checkpoints = {
        str(key): validate_sha256(
            value, context=f"scoring authorized checkpoint {key}"
        )
        for key, value in _mapping(
            authorization.get("teacher_checkpoint_sha256_by_fold"),
            context="scoring authorized checkpoints",
        ).items()
    }
    if authorized_checkpoints != frozen_teacher["checkpoint_sha256_by_fold"]:
        raise ValueError("scoring authorization checkpoints differ from the freeze")
    document_bindings = _mapping(
        authorization.get("teacher_release_document_bindings"),
        context="scoring authorized release documents",
    )
    for key, expected in frozen_teacher["documents"].items():
        if validate_sha256(
            document_bindings.get(key), context=f"authorized release {key}"
        ) != expected:
            raise ValueError(f"scoring authorization release {key} drifted")
    scoring_inputs = _mapping(
        scoring.get("input_artifacts"), context="scoring input artifacts"
    )
    for input_name, expected_hash_name in (
        ("s63_launch_manifest", "launch_manifest_sha256"),
        ("s63_preregistered_contract", "prospective_plan_sha256"),
        ("candidate_table", "candidate_table_sha256"),
    ):
        record = _mapping(
            scoring_inputs.get(input_name),
            context=f"scoring input_artifacts.{input_name}",
        )
        before = validate_sha256(
            record.get("sha256_before", record.get("sha256")),
            context=f"scoring input {input_name} before",
        )
        after = validate_sha256(
            record.get("sha256_after", record.get("sha256")),
            context=f"scoring input {input_name} after",
        )
        if before != hashes[expected_hash_name] or after != before:
            raise ValueError(f"scoring input {input_name} hash binding drifted")
    return {
        "schema_version": QUEUE_LAUNCH_BINDING_SCHEMA,
        "passed": True,
        "launch_manifest_sha256": hashes["launch_manifest_sha256"],
        "teacher_score_table_sha256": observed_score_sha,
        "teacher_score_summary_sha256": hashes["teacher_score_summary_sha256"],
        "execution_checkout": dict(execution_checkout),
        "scoring_execution_checkout": dict(scoring_execution_checkout),
        "publication_checkout_audit": dict(publication_checkout),
        "teacher_v3_release": teacher_release_binding,
        "producer_git_sha": launch_producer_git_sha,
        "verified_launch_artifacts": verified_launch_artifacts,
        "score_bound_transitively_through_launch_authorization": True,
    }


def _cohort_frame(tics: set[int], cohort: str) -> pd.DataFrame:
    frame = pd.DataFrame({"tic": np.asarray(sorted(tics), dtype=np.int64)})
    frame["sector"] = np.int16(S63_SECTOR)
    frame["prospective_cohort"] = cohort
    frame["cohort_contract_version"] = COHORT_CONTRACT_VERSION
    return frame.loc[:, ["sector", "tic", "prospective_cohort", "cohort_contract_version"]]


def derive_s63_prospective_cohorts(
    *,
    model_ready_tics: Sequence[int],
    frozen_corpus_tics: Sequence[int],
    reserved_tics: Sequence[int],
    source_hashes: Mapping[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Partition every model-ready S63 TIC by frozen corpus membership.

    The third table is the reservation complement: reserved TICs that were not
    model-ready.  Together, it and the two model cohorts must exactly partition
    the original reservation.
    """

    required = {
        "model_ready_allowlist_sha256",
        "teacher_v3_corpus_sha256",
        "reserved_tics_sha256",
    }
    missing = sorted(required - set(source_hashes))
    if missing:
        raise KeyError(f"cohort source hashes are missing: {missing}")
    hashes = {
        key: validate_sha256(source_hashes[key], context=key) for key in sorted(required)
    }
    model_values = _normalize_tic_values(model_ready_tics, context="model-ready allowlist")
    corpus_values = _normalize_tic_values(frozen_corpus_tics, context="Teacher-v3 corpus")
    reserved_values = _normalize_tic_values(reserved_tics, context="S63 reservation")
    if len(model_values) != len(set(model_values)):
        raise ValueError("model-ready allowlist contains duplicate TICs")
    if len(reserved_values) != len(set(reserved_values)):
        raise ValueError("S63 reservation contains duplicate TICs")
    model = set(model_values)
    corpus = set(corpus_values)
    reserved = set(reserved_values)
    if not model:
        raise ValueError("model-ready allowlist is empty")
    if not model.issubset(reserved):
        examples = sorted(model - reserved)[:10]
        raise ValueError(f"model-ready TICs fall outside the S63 reservation: {examples}")

    primary_tics = model - corpus
    repeated_tics = model & corpus
    not_ready_tics = reserved - model
    if (primary_tics & repeated_tics) or (primary_tics | repeated_tics) != model:
        raise RuntimeError("prospective cohorts do not exactly partition model-ready TICs")
    partitions = (primary_tics, repeated_tics, not_ready_tics)
    if any(left & right for index, left in enumerate(partitions) for right in partitions[index + 1 :]):
        raise RuntimeError("S63 reservation partitions overlap")
    if set().union(*partitions) != reserved:
        raise RuntimeError("S63 reservation complement is incomplete")

    primary = _cohort_frame(primary_tics, PRIMARY_COHORT)
    repeated = _cohort_frame(repeated_tics, REPEATED_HOST_COHORT)
    not_ready = _cohort_frame(not_ready_tics, RESERVATION_EXCLUSION)
    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": COHORT_CONTRACT_VERSION,
        "sector": S63_SECTOR,
        "source_hashes": hashes,
        "counts": {
            "reserved": len(reserved),
            "model_ready": len(model),
            "primary_teacher_v3_disjoint": len(primary_tics),
            "repeated_host": len(repeated_tics),
            "reserved_not_model_ready": len(not_ready_tics),
        },
        "tic_inventory_sha256": {
            "reserved": tic_inventory_sha256(reserved),
            "model_ready": tic_inventory_sha256(model),
            "primary_teacher_v3_disjoint": tic_inventory_sha256(primary_tics),
            "repeated_host": tic_inventory_sha256(repeated_tics),
            "reserved_not_model_ready": tic_inventory_sha256(not_ready_tics),
        },
        "partition_checks": {
            "model_ready_subset_of_reservation": True,
            "primary_repeated_disjoint": True,
            "primary_union_repeated_equals_model_ready": True,
            "reservation_complement_exact": True,
        },
        "human_label_interpretation": "human-confirmed Planet-like transit morphology",
        "confirmed_exoplanet_claimed": False,
    }
    return primary, repeated, not_ready, summary


def write_s63_prospective_cohorts(
    *,
    model_ready_path: Path,
    teacher_v3_corpus_path: Path,
    reserved_tics_path: Path,
    out_dir: Path,
    source_hashes: Mapping[str, str],
    producer_git_sha: str,
    expected_reserved_count: int | None = None,
    expected_corpus_count: int | None = None,
) -> dict[str, Any]:
    """Build and write canonical cohort inventories with exact source binding."""

    model = read_tic_inventory(
        model_ready_path,
        expected_sha256=source_hashes["model_ready_allowlist_sha256"],
    )
    corpus = read_tic_inventory(
        teacher_v3_corpus_path,
        expected_sha256=source_hashes["teacher_v3_corpus_sha256"],
        allow_duplicate_rows=True,
    )
    reserved = read_tic_inventory(
        reserved_tics_path,
        expected_sha256=source_hashes["reserved_tics_sha256"],
    )
    if expected_reserved_count is not None and len(reserved) != int(expected_reserved_count):
        raise ValueError(
            "S63 reservation count disagrees with the prospective plan: "
            f"{len(reserved)} != {int(expected_reserved_count)}"
        )
    if expected_corpus_count is not None and len(corpus) != int(expected_corpus_count):
        raise ValueError(
            "Teacher-v3 corpus TIC count disagrees with the prospective plan: "
            f"{len(corpus)} != {int(expected_corpus_count)}"
        )
    primary, repeated, not_ready, summary = derive_s63_prospective_cohorts(
        model_ready_tics=model,
        frozen_corpus_tics=corpus,
        reserved_tics=reserved,
        source_hashes=source_hashes,
    )
    summary["producer_git_sha"] = validate_git_sha(
        producer_git_sha, context="cohort producer_git_sha"
    )
    out_dir = Path(out_dir).expanduser().resolve(strict=False)
    if out_dir.exists():
        raise FileExistsError(f"cohort output must be a fresh directory: {out_dir}")
    claims = _claim_publication_targets((out_dir,))
    staging = out_dir.parent / f".{out_dir.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    try:
        if out_dir.exists():
            raise FileExistsError(f"cohort output appeared during build: {out_dir}")
        staging.mkdir()
        final_outputs = {
            "primary_cohort": out_dir / "primary_teacher_v3_disjoint_tics.txt",
            "repeated_host_cohort": out_dir / "repeated_host_tics.txt",
            "reserved_not_model_ready": out_dir / "reserved_not_model_ready_tics.txt",
        }
        staged_outputs = {
            key: staging / final_path.name for key, final_path in final_outputs.items()
        }
        for key, frame in (
            ("primary_cohort", primary),
            ("repeated_host_cohort", repeated),
            ("reserved_not_model_ready", not_ready),
        ):
            staged_outputs[key].write_text(
                "".join(f"{tic}\n" for tic in frame["tic"].astype(np.int64)),
                encoding="ascii",
            )
        summary["outputs"] = {
            key: {
                "path": str(final_outputs[key]),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for key, path in staged_outputs.items()
        }
        final_summary = out_dir / "cohort_summary.json"
        summary["outputs"]["summary"] = {"path": str(final_summary)}
        staged_summary = staging / final_summary.name
        staged_summary.write_text(
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        staged_summary_sha256 = file_sha256(staged_summary)
        staged_summary_size = staged_summary.stat().st_size
        if out_dir.exists():
            raise FileExistsError(f"cohort output appeared before publication: {out_dir}")
        _replace_directory(staging, out_dir)
        summary["outputs"]["summary"].update(
            {
                "sha256": staged_summary_sha256,
                "size_bytes": staged_summary_size,
            }
        )
        return summary
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    finally:
        _release_publication_claims(claims)


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], *, context: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise KeyError(f"{context} is missing columns: {missing}")


def validate_s63_bls_summary(
    summary: Mapping[str, Any],
    *,
    peak_table_sha256: str,
    model_ready_allowlist_sha256: str,
    model_ready_tics_sha256: str,
    n_model_ready_tics: int,
) -> dict[str, Any]:
    """Validate the prospective BLS product without asserting S56 Tier-1 QA."""

    if not isinstance(summary, Mapping):
        raise TypeError("S63 BLS summary must be a JSON object")
    if summary.get("passed") is not True:
        raise ValueError("S63 BLS summary did not pass its own validation")
    if summary.get("sector") != S63_SECTOR:
        raise ValueError("S63 BLS summary has the wrong sector")
    if summary.get("n_shards") != 1 or summary.get("shard_index") != 0:
        raise ValueError("S63 candidate construction requires a complete merged BLS table")
    if summary.get("n_targets") != n_model_ready_tics or summary.get("n_targets_total") != n_model_ready_tics:
        raise ValueError("S63 BLS target count does not match the model-ready allowlist")
    exact_hashes = {
        "peak_table_sha256": peak_table_sha256,
        "target_allowlist_sha256": model_ready_allowlist_sha256,
        "target_allowlist_tics_sha256": model_ready_tics_sha256,
    }
    for field, expected in exact_hashes.items():
        observed = validate_sha256(summary.get(field), context=f"S63 BLS {field}")
        if observed != validate_sha256(expected, context=f"expected {field}"):
            raise ValueError(f"S63 BLS {field} does not match its frozen input")
    config = summary.get("config")
    if not isinstance(config, Mapping):
        raise ValueError("S63 BLS summary is missing config")
    expected_config = {
        "apertures": list(ADP_ONLY_APERTURES),
        "n_periods": 50_000,
        "n_peaks": 10,
        "source_product_tag": "A2v1",
    }
    for field, expected in expected_config.items():
        if config.get(field) != expected:
            raise ValueError(f"S63 BLS config has incompatible {field}")
    if summary.get("orbitid_policy") != "strict":
        raise ValueError("S63 prospective BLS must use strict orbit IDs")
    return {
        "status": "pass",
        "scope": PROSPECTIVE_USE_SCOPE,
        "science_ready": False,
        "promotion_enabled": False,
        "sector": S63_SECTOR,
        **{key: validate_sha256(value, context=key) for key, value in exact_hashes.items()},
    }


def _strict_numeric(frame: pd.DataFrame, column: str, *, positive: bool = False) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="coerce")
    valid = np.isfinite(values.to_numpy(dtype=float))
    if positive:
        valid &= values.to_numpy(dtype=float) > 0
    if not bool(np.all(valid)):
        raise ValueError(f"S63 BLS rows contain invalid {column}")
    return values


def build_s63_rank_one_candidates(
    peaks: pd.DataFrame,
    *,
    model_ready_tics: Sequence[int],
    artifact_hashes: Mapping[str, str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build one rank-one ADP-small candidate per model-ready S63 TIC."""

    required_hashes = {
        "prospective_plan_sha256",
        "model_ready_allowlist_sha256",
        "bls_merged_table_sha256",
    }
    missing_hashes = sorted(required_hashes - set(artifact_hashes))
    if missing_hashes:
        raise KeyError(f"candidate artifact hashes are missing: {missing_hashes}")
    hashes = {
        key: validate_sha256(artifact_hashes[key], context=key)
        for key in sorted(required_hashes)
    }
    required_columns = (
        "sector",
        "tic",
        "aperture",
        "peak_rank",
        "status",
        "source_product_tag",
        "period_d",
        "t0_bjd",
        "duration_min",
        "depth",
        "depth_snr",
        "sde",
        "log_power",
        "target_metadata_contract_version",
        *TEACHER_V3_CANDIDATE_METADATA_COLUMNS,
    )
    _require_columns(peaks, required_columns, context="S63 BLS table")
    if peaks.empty:
        raise ValueError("S63 BLS table is empty")
    frame = peaks.copy()
    frame["sector"] = pd.to_numeric(frame["sector"], errors="coerce").astype("Int64")
    frame["tic"] = pd.to_numeric(frame["tic"], errors="coerce").astype("Int64")
    if frame[["sector", "tic"]].isna().any().any():
        raise ValueError("S63 BLS table has invalid sector/TIC values")
    if set(frame["sector"].astype(int)) != {S63_SECTOR}:
        raise ValueError("candidate table must contain only Sector 63 BLS rows")
    if set(frame["aperture"].astype(str)) != set(ADP_ONLY_APERTURES):
        raise ValueError("S63 BLS table must contain exactly the ADP-small/primary pair")
    if set(frame["source_product_tag"].astype(str)) != {"A2v1"}:
        raise ValueError("S63 BLS rows must come from A2v1")
    model_values = _normalize_tic_values(model_ready_tics, context="model-ready allowlist")
    if len(model_values) != len(set(model_values)):
        raise ValueError("model-ready allowlist contains duplicate TICs")
    model = set(model_values)
    observed = set(frame["tic"].astype(np.int64))
    if observed != model:
        raise ValueError(
            "S63 BLS TIC coverage differs from model-ready allowlist; "
            f"missing={sorted(model - observed)[:10]}, extra={sorted(observed - model)[:10]}"
        )
    rank = pd.to_numeric(frame["peak_rank"], errors="coerce")
    rank_one = rank.eq(1)
    good = frame["status"].astype(str).eq("ok")
    small = frame.loc[
        rank_one & good & frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[0])
    ].copy()
    primary = frame.loc[
        rank_one & good & frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[1])
    ].copy()
    for name, selected in (("ADP-small", small), ("ADP-primary", primary)):
        if selected["tic"].duplicated().any():
            raise ValueError(f"S63 BLS table has duplicate rank-one {name} rows")
        selected_tics = set(selected["tic"].astype(np.int64))
        if selected_tics != model:
            raise ValueError(
                f"model-ready TICs lack successful rank-one {name} results: "
                f"{sorted(model - selected_tics)[:10]}"
            )
        for column in ("period_d", "t0_bjd", "duration_min"):
            _strict_numeric(selected, column, positive=column != "t0_bjd")
        if set(selected["target_metadata_contract_version"].astype(str)) != {
            S63_TARGET_METADATA_CONTRACT
        }:
            raise ValueError(f"rank-one {name} target metadata contract drifted")
        for column in TEACHER_V3_CANDIDATE_METADATA_COLUMNS:
            raw = selected[column]
            numeric = pd.to_numeric(raw, errors="coerce")
            nonempty = raw.notna() & raw.astype(str).str.strip().ne("")
            if (nonempty & numeric.isna()).any():
                raise ValueError(
                    f"rank-one {name} metadata {column} contains nonnumeric values"
                )
            if not numeric.notna().any():
                raise ValueError(
                    f"rank-one {name} metadata {column} is wholly absent"
                )
            selected[column] = numeric
    small = small.sort_values("tic", kind="stable").reset_index(drop=True)
    primary = primary.sort_values("tic", kind="stable").reset_index(drop=True)
    for column in TEACHER_V3_CANDIDATE_METADATA_COLUMNS:
        if not np.array_equal(
            small[column].to_numpy(dtype=float),
            primary[column].to_numpy(dtype=float),
            equal_nan=True,
        ):
            raise ValueError(
                f"rank-one aperture copies of target metadata {column} differ"
            )
    peak_fields = (
        "peak_rank",
        "period_d",
        "t0_bjd",
        "duration_min",
        "depth",
        "depth_snr",
        "sde",
        "log_power",
    )
    context = primary.loc[:, ["tic", *peak_fields]].rename(
        columns={column: f"adp_{column}" for column in peak_fields}
    )
    out = small.merge(context, on="tic", how="left", validate="one_to_one")
    for column in peak_fields:
        out[f"adp_sml_{column}"] = pd.to_numeric(out[column], errors="coerce")
    out["tic"] = out["tic"].astype(np.int64)
    out["sector"] = np.int16(S63_SECTOR)
    out["rep_aperture"] = ADP_ONLY_APERTURES[0]
    out["rep_peak_rank"] = np.int8(1)
    out["sde_max"] = out["adp_sml_sde"]
    out["anchor_aperture"] = ADP_ONLY_APERTURES[0]
    out["anchor_period_d"] = out["adp_sml_period_d"]
    out["anchor_t0_bjd"] = out["adp_sml_t0_bjd"]
    out["anchor_duration_min"] = out["adp_sml_duration_min"]
    out["anchor_sde"] = out["adp_sml_sde"]
    relation = classify_period_relation(out["adp_sml_period_d"], out["adp_period_d"])
    out["aperture_period_relation"] = relation
    out["aperture_period_rel_delta"] = (
        np.abs(out["adp_period_d"] - out["adp_sml_period_d"])
        / out["adp_sml_period_d"]
    )
    out["aperture_depth_ratio_primary_over_small"] = (
        pd.to_numeric(out["adp_depth"], errors="coerce")
        / pd.to_numeric(out["adp_sml_depth"], errors="coerce").replace(0, np.nan)
    )
    agreement = relation.isin(("exact", "harmonic"))
    out["aperture_disagreement_flag"] = ~agreement
    out["n_apertures_agree"] = np.where(agreement, 2, 1)
    out["apertures_agree"] = np.where(
        agreement,
        ",".join(ADP_ONLY_APERTURES),
        ADP_ONLY_APERTURES[0],
    )
    out["review_id"] = [
        f"s0063-A2v1-adp-small-{tic:016d}-r1" for tic in out["tic"]
    ]
    out["source_kind"] = "real_candidate"
    out["source_bucket"] = ""
    out["source_product_tag"] = "A2v1"
    out["candidate_provenance_contract_version"] = CANDIDATE_CONTRACT_VERSION
    out["prospective_use_scope"] = PROSPECTIVE_USE_SCOPE
    out["science_ready"] = False
    out["promotion_enabled"] = False
    out["native_input_include"] = True
    out["native_group_path"] = [native_group_path(row) for row in out.to_dict("records")]
    out["candidate_key"] = out.apply(candidate_key, axis=1)
    for key, value in hashes.items():
        out[key] = value
    if out["tic"].nunique() != len(model) or len(out) != len(model):
        raise RuntimeError("S63 candidate table is not one unique row per model-ready TIC")
    out = out.sort_values("tic", kind="stable").reset_index(drop=True)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": CANDIDATE_CONTRACT_VERSION,
        "sector": S63_SECTOR,
        "prospective_use_scope": PROSPECTIVE_USE_SCOPE,
        "science_ready": False,
        "promotion_enabled": False,
        "candidate_policy": "one successful ADP-small rank-one row per model-ready TIC",
        "context_policy": "successful ADP-primary rank-one context required",
        "teacher_v3_metadata_columns": list(TEACHER_V3_CANDIDATE_METADATA_COLUMNS),
        "n_rows": len(out),
        "n_unique_tics": int(out["tic"].nunique()),
        "artifact_hashes": hashes,
        "candidate_tics_sha256": tic_inventory_sha256(set(out["tic"].astype(int))),
    }
    return out, summary


def _numeric(frame: pd.DataFrame, column: str, default: float = np.nan) -> pd.Series:
    if column not in frame:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _probability(frame: pd.DataFrame, column: str) -> pd.Series:
    values = _numeric(frame, column)
    if values.isna().any() or (~values.between(0.0, 1.0)).any():
        raise ValueError(f"Teacher-v3 score table has invalid {column}")
    return values


def _selection_scores(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    p_planet = _probability(out, "p_planet_like")
    p_eb = _probability(out, "p_eclipse_contact")
    p_variable = _probability(out, "p_smooth_variable")
    p_other = _probability(out, "p_other")
    morphology_sum = p_planet + p_eb + p_variable + p_other
    if not np.allclose(morphology_sum, 1.0, rtol=0.0, atol=1.0e-5):
        raise ValueError("Teacher-v3 morphology probabilities do not sum to one")
    p_preserve = _probability(out, "p_preserve")
    period = _numeric(out, "period_d")
    duration = _numeric(out, "duration_min")
    sde = _numeric(out, "sde_max")
    if not (
        np.isfinite(period).all()
        and np.isfinite(duration).all()
        and np.isfinite(sde).all()
        and period.gt(0).all()
        and duration.gt(0).all()
    ):
        raise ValueError("Teacher-v3 score table has invalid rank-one BLS evidence")
    duty_cycle = duration / (period * 1440.0)
    broad_rank = duty_cycle.rank(pct=True, method="average")
    out["broad_duty_cycle"] = duty_cycle
    out["broad_duty_cycle_percentile"] = broad_rank
    long_low_sde = ((period >= 1.0) & (sde <= 40.0)).astype(float)
    entropy = _numeric(out, "morphology_entropy")
    disagreement = _numeric(out, "ensemble_disagreement")
    if (
        entropy.isna().any()
        or disagreement.isna().any()
        or entropy.lt(0.0).any()
        or entropy.gt(np.log(4.0) + 1.0e-5).any()
        or disagreement.lt(0.0).any()
    ):
        raise ValueError("Teacher-v3 score table lacks finite disagreement quantities")
    harmonic_columns = [
        column for column in out if column.startswith("p_harmonic_") and column != "p_harmonic_P"
    ]
    if not harmonic_columns:
        raise KeyError("Teacher-v3 score table has no non-P harmonic probabilities")
    all_harmonic_columns = [column for column in out if column.startswith("p_harmonic_")]
    harmonic_all = out.loc[:, all_harmonic_columns].apply(pd.to_numeric, errors="coerce")
    if (
        harmonic_all.isna().any().any()
        or ((harmonic_all < 0.0) | (harmonic_all > 1.0)).any().any()
        or not np.allclose(harmonic_all.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-5)
    ):
        raise ValueError("Teacher-v3 score table has invalid harmonic probabilities")
    harmonic = harmonic_all.loc[:, harmonic_columns]
    rare_harmonic = harmonic.max(axis=1)
    out["_score_planet_preserve"] = 0.55 * p_planet + 0.35 * p_preserve + 0.10 * long_low_sde
    out["_score_eclipse_contact"] = p_eb
    out["_score_smooth_variable"] = p_variable
    out["_score_broad_dip"] = p_preserve * broad_rank * (1.0 - 0.5 * p_eb) * (1.0 - 0.5 * p_variable)
    out["_score_disagreement_harmonic"] = 0.45 * entropy + 0.35 * disagreement + 0.20 * rare_harmonic
    return out


def _stable_hash(seed: int, namespace: str, tic: int) -> str:
    return hashlib.sha256(f"{seed}|{namespace}|{tic}".encode("ascii")).hexdigest()


def _top_unique(
    candidates: pd.DataFrame,
    *,
    score_column: str,
    count: int,
    used: set[int],
) -> pd.DataFrame:
    available = candidates.loc[~candidates["tic"].isin(used)].copy()
    available = available.loc[np.isfinite(pd.to_numeric(available[score_column], errors="coerce"))]
    available = available.sort_values(
        [score_column, "sde_max", "tic"],
        ascending=[False, False, True],
        kind="stable",
    )
    if len(available) < count:
        raise ValueError(
            f"selection bucket {score_column} requested {count} unique TICs but found {len(available)}"
        )
    selected = available.head(int(count)).copy()
    selected["selection_rank"] = np.arange(1, count + 1, dtype=int)
    return selected


def _stratified_control(
    candidates: pd.DataFrame,
    *,
    count: int,
    used: set[int],
    seed: int,
    cohort: str,
) -> pd.DataFrame:
    available = candidates.loc[~candidates["tic"].isin(used)].copy()
    available = available.sort_values("tic", kind="stable")
    if len(available) < count:
        raise ValueError(f"stratified control requested {count} TICs but found {len(available)}")
    ranks = pd.to_numeric(available["sde_max"], errors="coerce").rank(method="first")
    try:
        strata = pd.qcut(ranks, q=4, labels=False, duplicates="raise")
    except ValueError as exc:
        raise ValueError("stratified control cannot form four ADP-small SDE quartiles") from exc
    available["control_sde_stratum"] = strata.astype(int) + 1
    draw_by_stratum = {
        stratum: count // 4 + (1 if stratum <= count % 4 else 0)
        for stratum in range(1, 5)
    }
    parts: list[pd.DataFrame] = []
    for stratum in range(1, 5):
        pool = available.loc[available["control_sde_stratum"].eq(stratum)].copy()
        draw_count = draw_by_stratum[stratum]
        if len(pool) < draw_count:
            raise ValueError(
                f"control SDE stratum {stratum} requested {draw_count} TICs but found {len(pool)}"
            )
        pool["_control_draw_hash"] = [
            _stable_hash(seed, f"{cohort}|control|q{stratum}", int(tic))
            for tic in pool["tic"]
        ]
        selected = pool.sort_values("_control_draw_hash", kind="stable").head(draw_count).copy()
        selected["control_stratum_population"] = len(pool)
        selected["control_stratum_draw_count"] = draw_count
        selected["selection_inclusion_probability"] = draw_count / len(pool)
        selected["selection_weight"] = len(pool) / draw_count
        parts.append(selected)
    output = pd.concat(parts, ignore_index=True, sort=False)
    output = output.sort_values(
        ["control_sde_stratum", "_control_draw_hash"], kind="stable"
    ).reset_index(drop=True)
    output["selection_rank"] = np.arange(1, len(output) + 1, dtype=int)
    return output


def _select_cohort(
    candidates: pd.DataFrame,
    *,
    cohort: str,
    quotas: ProspectiveQuotas,
    seed: int,
) -> pd.DataFrame:
    if len(candidates) < quotas.total:
        raise ValueError(
            f"{cohort} has {len(candidates)} candidate TICs but requires {quotas.total}"
        )
    used: set[int] = set()
    parts: list[pd.DataFrame] = []
    specs = (
        ("planet_preserve", quotas.planet_preserve, "_score_planet_preserve"),
        ("eclipse_contact", quotas.eclipse_contact, "_score_eclipse_contact"),
        ("smooth_variable", quotas.smooth_variable, "_score_smooth_variable"),
        ("broad_dip", quotas.broad_dip, "_score_broad_dip"),
        ("disagreement_harmonic", quotas.disagreement_harmonic, "_score_disagreement_harmonic"),
    )
    for bucket, count, score_column in specs:
        if count <= 0:
            continue
        pool_size = int((~candidates["tic"].isin(used)).sum())
        selected = _top_unique(
            candidates,
            score_column=score_column,
            count=count,
            used=used,
        )
        selected["selection_bucket"] = bucket
        selected["selection_score"] = selected[score_column]
        selected["deterministic_pool_size"] = pool_size
        selected["deterministic_cut_fraction"] = count / pool_size
        selected["selection_inclusion_probability"] = np.nan
        selected["selection_weight"] = np.nan
        selected["control_sde_stratum"] = pd.NA
        selected["control_stratum_population"] = pd.NA
        selected["control_stratum_draw_count"] = pd.NA
        used.update(int(value) for value in selected["tic"])
        parts.append(selected)
    pool_size = int((~candidates["tic"].isin(used)).sum())
    controls = _stratified_control(
        candidates,
        count=quotas.stratified_control,
        used=used,
        seed=seed,
        cohort=cohort,
    )
    controls["selection_bucket"] = "stratified_control"
    controls["selection_score"] = np.nan
    controls["control_candidate_pool_size"] = pool_size
    controls["deterministic_pool_size"] = pd.NA
    controls["deterministic_cut_fraction"] = np.nan
    used.update(int(value) for value in controls["tic"])
    parts.append(controls)
    selected = pd.concat(parts, ignore_index=True, sort=False)
    if len(selected) != quotas.total or selected["tic"].nunique() != quotas.total:
        raise RuntimeError(f"{cohort} selection did not produce exact unique-TIC quotas")
    selected["prospective_cohort"] = cohort
    selected["selection_score_population_size"] = len(candidates)
    selected["selection_score_population_scope"] = (
        f"{cohort}_complete_cohort_before_depletion"
    )
    return selected


def _validate_score_candidates(scores: pd.DataFrame, expected_tics: set[int]) -> pd.DataFrame:
    required = (
        "score_schema_version",
        "tic",
        "sector",
        "review_id",
        "period_d",
        "t0_bjd",
        "duration_min",
        "sde_max",
        "rep_aperture",
        "rep_peak_rank",
        "candidate_provenance_contract_version",
        "prospective_use_scope",
        "science_ready",
        "promotion_enabled",
        "p_planet_like",
        "p_eclipse_contact",
        "p_smooth_variable",
        "p_other",
        "p_preserve",
        "morphology_entropy",
        "ensemble_disagreement",
        "model_version",
        "model_profile",
        "checkpoint_namespace",
        "input_contract_version",
        "teacher_run_id",
        "teacher_training_sectors",
        "scoring_sector",
        "score_policy",
        *TEACHER_V3_CANDIDATE_METADATA_COLUMNS,
    )
    _require_columns(scores, required, context="Teacher-v3 S63 score table")
    frame = scores.copy()
    frame["tic"] = pd.to_numeric(frame["tic"], errors="coerce").astype("Int64")
    frame["sector"] = pd.to_numeric(frame["sector"], errors="coerce").astype("Int64")
    if frame[["tic", "sector"]].isna().any().any():
        raise ValueError("Teacher-v3 score table has invalid sector/TIC values")
    if frame["tic"].duplicated().any():
        raise ValueError("Teacher-v3 score table must have one rank-one row per TIC")
    frame["tic"] = frame["tic"].astype(np.int64)
    if set(frame["tic"]) != expected_tics:
        raise ValueError("Teacher-v3 score TICs do not exactly equal the cohort union")
    if set(frame["sector"].astype(int)) != {S63_SECTOR}:
        raise ValueError("Teacher-v3 prospective score table must contain only S63")
    if set(frame["rep_aperture"].astype(str)) != {ADP_ONLY_APERTURES[0]}:
        raise ValueError("Teacher-v3 prospective scores must use ADP-small anchors")
    if not pd.to_numeric(frame["rep_peak_rank"], errors="coerce").eq(1).all():
        raise ValueError("Teacher-v3 prospective scores must use rank-one candidates")
    if set(frame["candidate_provenance_contract_version"].astype(str)) != {CANDIDATE_CONTRACT_VERSION}:
        raise ValueError("Teacher-v3 scores have an incompatible S63 candidate contract")
    if set(frame["prospective_use_scope"].astype(str)) != {PROSPECTIVE_USE_SCOPE}:
        raise ValueError("Teacher-v3 scores have an incompatible use scope")
    if frame["science_ready"].astype(bool).any() or frame["promotion_enabled"].astype(bool).any():
        raise ValueError("S63 prospective scores must not claim readiness or promotion")
    expected_score_identity = {
        "score_schema_version": S63_SCORE_TABLE_SCHEMA,
        "model_version": MODEL_VERSION,
        "model_profile": TEACHER_V3_PRIMARY_PROFILE,
        "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "teacher_run_id": TEACHER_V3_RUN_ID,
        "teacher_training_sectors": "56-62",
        "scoring_sector": S63_SECTOR,
        "score_policy": S63_SCORE_POLICY,
    }
    for column, expected in expected_score_identity.items():
        if set(frame[column]) != {expected}:
            raise ValueError(f"Teacher-v3 score table has incompatible {column}")
    if "source_kind" in frame and frame["source_kind"].astype(str).str.contains("inject", case=False).any():
        raise ValueError("S63 prospective score table contains injected rows")
    return frame


def build_s63_prospective_review_queue(
    scores: pd.DataFrame,
    *,
    primary_tics: Sequence[int],
    repeated_host_tics: Sequence[int],
    artifact_hashes: Mapping[str, str],
    selection_policy: Mapping[str, Any],
    launch_binding: Mapping[str, Any],
    producer_git_sha: str,
    selection_seed: int = 630056,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, Any]]:
    """Select and blind the frozen 1000+100 S63 morphology review queue."""

    missing_hashes = sorted(set(REQUIRED_ARTIFACT_HASHES) - set(artifact_hashes))
    if missing_hashes:
        raise KeyError(f"queue artifact hashes are missing: {missing_hashes}")
    hashes = {
        key: validate_sha256(value, context=key)
        for key, value in sorted(artifact_hashes.items())
    }
    validate_selection_policy(selection_policy)
    if hashes["selection_policy_sha256"] != validate_sha256(
        artifact_hashes["selection_policy_sha256"], context="selection policy"
    ):
        raise RuntimeError("selection-policy hash normalization changed")
    launch_binding = _mapping(
        launch_binding, context="queue launch-binding verification"
    )
    if (
        launch_binding.get("schema_version") != QUEUE_LAUNCH_BINDING_SCHEMA
        or launch_binding.get("passed") is not True
    ):
        raise ValueError("queue construction requires a passed launch binding")
    for key in (
        "launch_manifest_sha256",
        "teacher_score_table_sha256",
        "teacher_score_summary_sha256",
    ):
        if validate_sha256(
            launch_binding.get(key), context=f"queue launch binding {key}"
        ) != hashes[key]:
            raise ValueError(f"queue launch binding disagrees on {key}")
    execution_checkout = _mapping(
        launch_binding.get("execution_checkout"),
        context="queue launch binding execution checkout",
    )
    scoring_execution_checkout = _mapping(
        launch_binding.get("scoring_execution_checkout"),
        context="queue launch binding scoring checkout",
    )
    publication_checkout = _mapping(
        launch_binding.get("publication_checkout_audit"),
        context="queue launch binding publication checkout",
    )
    if not (
        dict(execution_checkout)
        == dict(scoring_execution_checkout)
        == dict(publication_checkout)
    ):
        raise ValueError("queue launch binding checkout audits are not identical")
    release_binding = _mapping(
        launch_binding.get("teacher_v3_release"),
        context="queue launch binding Teacher-v3 release",
    )
    if release_binding.get("score_policy") != S63_SCORE_POLICY:
        raise ValueError("queue launch binding has an incompatible score policy")
    release_identity = _mapping(
        release_binding.get("identity"),
        context="queue launch binding Teacher-v3 identity",
    )
    expected_release_identity = {
        "schema_version": "twirl_teacher_v3_s63_ensemble_execution_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "profile": TEACHER_V3_PRIMARY_PROFILE,
        "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "sector": S63_SECTOR,
    }
    if dict(release_identity) != expected_release_identity:
        raise ValueError("queue launch binding Teacher-v3 identity drifted")
    release_checkpoints = _mapping(
        release_binding.get("checkpoint_sha256_by_fold"),
        context="queue launch binding Teacher-v3 checkpoints",
    )
    if set(release_checkpoints) != {"0", "1", "2", "3", "4"}:
        raise ValueError("queue launch binding does not contain five checkpoints")
    for fold, digest in release_checkpoints.items():
        validate_sha256(digest, context=f"queue launch binding checkpoint {fold}")
    for key, expected in {
        "schema_version": EXECUTION_CHECKOUT_AUDIT_SCHEMA,
        "passed": True,
        "checkout_clean": True,
        "untracked_files_checked": True,
        "status_command": "git status --porcelain --untracked-files=all",
    }.items():
        if execution_checkout.get(key) != expected:
            raise ValueError(
                f"queue launch binding execution checkout has incompatible {key}"
            )
    execution_repo = str(execution_checkout.get("repo", "")).strip()
    execution_declared_repo = str(
        execution_checkout.get("launch_declared_repo", "")
    )
    if not execution_repo:
        raise ValueError("queue launch binding execution checkout has no repo")
    execution_git_sha = str(execution_checkout.get("observed_git_sha", "")).strip()
    if _GIT_SHA_PATTERN.fullmatch(execution_git_sha) is None:
        raise ValueError("queue launch binding execution checkout has invalid Git SHA")
    for key in ("operator_expected_git_sha", "launch_git_sha"):
        if execution_checkout.get(key) != execution_git_sha:
            raise ValueError(
                f"queue launch binding execution checkout disagrees on {key}"
            )
    producer_git_sha = validate_git_sha(
        producer_git_sha, context="queue producer_git_sha"
    )
    if producer_git_sha != execution_git_sha:
        raise ValueError(
            "queue producer_git_sha differs from the audited execution checkout"
        )
    if launch_binding.get("producer_git_sha") != producer_git_sha:
        raise ValueError("queue launch binding producer_git_sha drifted")
    if selection_seed != 630056:
        raise ValueError("S63 prospective selection seed is frozen at 630056")
    primary_values = _normalize_tic_values(primary_tics, context="primary cohort")
    repeated_values = _normalize_tic_values(repeated_host_tics, context="repeated-host cohort")
    if len(primary_values) != len(set(primary_values)) or len(repeated_values) != len(set(repeated_values)):
        raise ValueError("prospective cohort inventories must contain unique TICs")
    primary = set(primary_values)
    repeated = set(repeated_values)
    if primary & repeated:
        raise ValueError("primary and repeated-host cohorts overlap")
    candidates = _validate_score_candidates(scores, primary | repeated)
    # Percentile-derived quantities are deliberately computed on each complete
    # cohort independently and before the first deterministic bucket is removed.
    primary_candidates = _selection_scores(
        candidates.loc[candidates["tic"].isin(primary)].copy()
    )
    repeated_candidates = _selection_scores(
        candidates.loc[candidates["tic"].isin(repeated)].copy()
    )
    primary_selected = _select_cohort(
        primary_candidates,
        cohort=PRIMARY_COHORT,
        quotas=PRIMARY_QUOTAS,
        seed=selection_seed,
    )
    repeated_selected = _select_cohort(
        repeated_candidates,
        cohort=REPEATED_HOST_COHORT,
        quotas=REPEATED_HOST_QUOTAS,
        seed=selection_seed,
    )
    hidden = pd.concat([primary_selected, repeated_selected], ignore_index=True, sort=False)
    hidden["source_candidate_review_id"] = hidden["review_id"].astype(str)
    hidden["queue_order_hash"] = [
        _stable_hash(selection_seed, "combined-review-order", int(tic)) for tic in hidden["tic"]
    ]
    hidden = hidden.sort_values("queue_order_hash", kind="stable").reset_index(drop=True)
    hidden["row_id"] = np.arange(len(hidden), dtype=int)
    hidden["review_id"] = [f"s0063-teacher-v3-review-{index:04d}" for index in hidden["row_id"]]
    hidden["source_bucket"] = ""
    hidden["candidate_key"] = hidden.apply(candidate_key, axis=1)
    hidden["queue_policy_version"] = QUEUE_POLICY_VERSION
    hidden["human_label_interpretation"] = "human-confirmed Planet-like transit morphology"
    hidden["teacher_scores_are_labels"] = False
    hidden["pseudo_labels_authorized"] = False
    hidden["execution_git_audit_schema"] = EXECUTION_CHECKOUT_AUDIT_SCHEMA
    hidden["execution_git_repo"] = execution_repo
    hidden["execution_git_launch_declared_repo"] = execution_declared_repo
    hidden["execution_git_sha"] = execution_git_sha
    hidden["execution_git_operator_expected_sha"] = execution_git_sha
    hidden["execution_git_launch_sha"] = execution_git_sha
    hidden["execution_git_checkout_clean"] = True
    hidden["execution_git_untracked_files_checked"] = True
    hidden["execution_git_status_command"] = execution_checkout["status_command"]
    hidden["producer_git_sha"] = producer_git_sha
    for key, value in hashes.items():
        hidden[key] = value
    hidden["twirl_vet_sheet_name"] = hidden["review_id"] + "_twirl_twoap_current_adp.png"
    hidden["twirl_vet_sheet_pdf_name"] = ""
    queue = hidden.reindex(columns=REVIEWER_COLUMNS).copy()
    for column in ("label", "label_source", "labeler", "notes", "updated_utc"):
        queue[column] = ""
    queue = queue.loc[:, PUBLIC_REVIEW_QUEUE_COLUMNS]
    public_summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": "twirl_teacher_v3_s63_public_review_queue_v1",
        "sector": S63_SECTOR,
        "queue_policy_version": QUEUE_POLICY_VERSION,
        "producer_git_sha": producer_git_sha,
        "n_queue_rows": len(queue),
        "n_unique_tics": int(queue["tic"].nunique()),
        "blinding_contract": (
            "teacher outputs and explicit selection annotations withheld; "
            "scientific identity remains visible"
        ),
        "ordering_contract": "seeded_identity_hash_after_selection",
        "identity_joinability_notice": (
            "Cohort is not explicitly annotated, but it is technically joinable "
            "from visible TIC and candidate identity; cohort-wise analysis is "
            "deferred until all 1100 labels are complete"
        ),
        "human_label_interpretation": (
            "human-confirmed Planet-like transit morphology"
        ),
        "confirmed_exoplanet_claimed": False,
        "review_completion_requirement": (
            "all 1100 frozen rows require an accepted human morphology decision "
            "before one-time unblinding and cohort-wise analysis"
        ),
    }
    hidden_summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": "twirl_teacher_v3_s63_hidden_selection_summary_v1",
        "queue_policy_version": QUEUE_POLICY_VERSION,
        "producer_git_sha": producer_git_sha,
        "selection_policy_schema_version": SELECTION_POLICY_SCHEMA,
        "sector": S63_SECTOR,
        "selection_seed": selection_seed,
        "prospective_use_scope": PROSPECTIVE_USE_SCOPE,
        "science_ready": False,
        "promotion_enabled": False,
        "scores_hidden_from_reviewer": True,
        "buckets_hidden_from_reviewer": True,
        "cohort_annotation_withheld": True,
        "cohort_joinable_from_visible_tic": True,
        "cohort_analysis_deferred_until_review_complete": True,
        "n_queue_rows": len(queue),
        "n_unique_tics": int(queue["tic"].nunique()),
        "cohort_counts": {
            str(key): int(value)
            for key, value in hidden["prospective_cohort"].value_counts().sort_index().items()
        },
        "bucket_counts_by_cohort": {
            cohort: {
                str(key): int(value)
                for key, value in group["selection_bucket"].value_counts().sort_index().items()
            }
            for cohort, group in hidden.groupby("prospective_cohort", sort=True)
        },
        "artifact_hashes": hashes,
        "launch_binding": dict(launch_binding),
        "review_completion_requirement": "all 1100 frozen rows require an accepted human morphology decision before one-time unblinding and cohort-wise analysis",
        "confirmed_exoplanet_claimed": False,
    }
    verify_s63_prospective_review_queue(
        queue,
        hidden,
        artifact_hashes=hashes,
    )
    return queue, hidden, public_summary, hidden_summary


def verify_s63_prospective_review_queue(
    queue: pd.DataFrame,
    hidden: pd.DataFrame,
    *,
    artifact_hashes: Mapping[str, str],
) -> None:
    """Fail closed on queue blinding, quota, identity, and hash drift."""

    if len(queue) != 1100 or len(hidden) != 1100:
        raise ValueError("S63 prospective review queue must contain exactly 1100 rows")
    if tuple(queue.columns) != PUBLIC_REVIEW_QUEUE_COLUMNS:
        raise ValueError("reviewer queue columns differ from the exact public allowlist")
    for name, frame in (("reviewer queue", queue), ("hidden provenance", hidden)):
        _require_columns(frame, ("row_id", "review_id", "candidate_key", "tic", "sector"), context=name)
        if frame["row_id"].duplicated().any() or frame["review_id"].duplicated().any():
            raise ValueError(f"{name} has duplicate row/review identities")
        if frame["tic"].nunique() != 1100:
            raise ValueError(f"{name} does not contain one unique TIC per row")
        row_ids = pd.to_numeric(frame["row_id"], errors="coerce").to_numpy(
            dtype=float
        )
        if not np.array_equal(row_ids, np.arange(1100, dtype=float)):
            raise ValueError(f"{name} row_id must be the ordered range 0..1099")
    visible_forbidden = {
        "prospective_cohort",
        "selection_bucket",
        "selection_score",
        "selection_rank",
        "control_sde_stratum",
        "control_stratum_population",
        "control_stratum_draw_count",
        "selection_inclusion_probability",
        "selection_weight",
        "deterministic_pool_size",
        "deterministic_cut_fraction",
        "selection_score_population_size",
        "selection_score_population_scope",
        "queue_order_hash",
        "ensemble_disagreement",
        "morphology_entropy",
    }
    visible_forbidden.update(column for column in queue if column.startswith(("p_", "std_p_", "_score_")))
    leaked = sorted(visible_forbidden & set(queue.columns))
    if leaked:
        raise ValueError(f"reviewer queue leaks model/selection metadata: {leaked}")
    identity_columns = ["row_id", "review_id", "candidate_key", "tic", "sector"]
    public_identity = queue.loc[:, identity_columns].reset_index(drop=True).astype(str)
    hidden_identity = hidden.loc[:, identity_columns].reset_index(drop=True).astype(str)
    if not public_identity.equals(hidden_identity):
        raise ValueError(
            "reviewer queue and hidden provenance are not row-for-row identical"
        )
    expected_order_hash = [
        _stable_hash(630056, "combined-review-order", int(tic))
        for tic in hidden["tic"]
    ]
    if hidden["queue_order_hash"].astype(str).tolist() != expected_order_hash:
        raise ValueError("hidden queue order is not the frozen seeded TIC-identity hash")
    if hidden["queue_order_hash"].astype(str).tolist() != sorted(expected_order_hash):
        raise ValueError("reviewer queue order is not ascending identity hash")
    expected_counts = {
        PRIMARY_COHORT: (PRIMARY_QUOTAS, 1000),
        REPEATED_HOST_COHORT: (REPEATED_HOST_QUOTAS, 100),
    }
    for cohort, (quotas, total) in expected_counts.items():
        subset = hidden.loc[hidden["prospective_cohort"].eq(cohort)]
        if len(subset) != total or subset["tic"].nunique() != total:
            raise ValueError(f"{cohort} does not contain its exact unique-TIC quota")
        counts = subset["selection_bucket"].value_counts().to_dict()
        if counts != {key: value for key, value in asdict(quotas).items()}:
            raise ValueError(f"{cohort} selection buckets do not match frozen quotas")
        population_size = pd.to_numeric(
            subset["selection_score_population_size"], errors="coerce"
        )
        if population_size.isna().any() or population_size.nunique() != 1:
            raise ValueError(f"{cohort} score-population size is inconsistent")
        if int(population_size.iloc[0]) < total:
            raise ValueError(f"{cohort} score population is smaller than its queue")
        if set(subset["selection_score_population_scope"].astype(str)) != {
            f"{cohort}_complete_cohort_before_depletion"
        }:
            raise ValueError(f"{cohort} score-population scope drifted")
        duty_percentile = pd.to_numeric(
            subset["broad_duty_cycle_percentile"], errors="coerce"
        )
        if duty_percentile.isna().any() or not duty_percentile.between(
            0.0, 1.0, inclusive="right"
        ).all():
            raise ValueError(f"{cohort} broad-duty percentiles are invalid")
        controls = subset.loc[subset["selection_bucket"].eq("stratified_control")]
        deterministic = subset.loc[
            ~subset["selection_bucket"].eq("stratified_control")
        ]
        if deterministic["selection_inclusion_probability"].notna().any():
            raise ValueError(
                f"{cohort} deterministic rows must have null inclusion probability"
            )
        if deterministic["selection_weight"].notna().any():
            raise ValueError(
                f"{cohort} deterministic rows must have null selection weight"
            )
        if (
            deterministic["deterministic_pool_size"].isna().any()
            or deterministic["deterministic_cut_fraction"].isna().any()
        ):
            raise ValueError(f"{cohort} deterministic cut accounting is incomplete")
        if not pd.to_numeric(
            deterministic["deterministic_cut_fraction"], errors="coerce"
        ).between(0.0, 1.0, inclusive="right").all():
            raise ValueError(f"{cohort} deterministic cut fractions are invalid")
        if controls["deterministic_cut_fraction"].notna().any():
            raise ValueError(f"{cohort} controls must not have deterministic cut fractions")
        if set(pd.to_numeric(controls["control_sde_stratum"]).astype(int)) != {1, 2, 3, 4}:
            raise ValueError(f"{cohort} control does not cover four SDE quartiles")
        for _, stratum in controls.groupby("control_sde_stratum", sort=True):
            population = pd.to_numeric(stratum["control_stratum_population"], errors="raise")
            draw = pd.to_numeric(stratum["control_stratum_draw_count"], errors="raise")
            fraction = pd.to_numeric(
                stratum["selection_inclusion_probability"], errors="raise"
            )
            if population.nunique() != 1 or draw.nunique() != 1 or draw.iloc[0] != len(stratum):
                raise ValueError(f"{cohort} control stratum accounting is inconsistent")
            if not np.allclose(fraction, draw.iloc[0] / population.iloc[0], rtol=0.0, atol=1e-15):
                raise ValueError(f"{cohort} control inclusion fractions are inconsistent")
    required_hidden = {
        "selection_score",
        "selection_rank",
        "selection_bucket",
        "prospective_cohort",
        "control_sde_stratum",
        "selection_inclusion_probability",
        "deterministic_cut_fraction",
        "selection_score_population_scope",
        "execution_git_audit_schema",
        "execution_git_repo",
        "execution_git_launch_declared_repo",
        "execution_git_sha",
        "execution_git_operator_expected_sha",
        "execution_git_launch_sha",
        "execution_git_checkout_clean",
        "execution_git_untracked_files_checked",
        "execution_git_status_command",
        "producer_git_sha",
        *REQUIRED_ARTIFACT_HASHES,
    }
    missing_hidden = sorted(required_hidden - set(hidden))
    if missing_hidden:
        raise KeyError(f"hidden provenance is missing required fields: {missing_hidden}")
    for key in REQUIRED_ARTIFACT_HASHES:
        expected = validate_sha256(artifact_hashes[key], context=key)
        observed = set(hidden[key].astype(str))
        if observed != {expected}:
            raise ValueError(f"hidden provenance has inconsistent {key}")
    expected_execution_values = {
        "execution_git_audit_schema": EXECUTION_CHECKOUT_AUDIT_SCHEMA,
        "execution_git_checkout_clean": True,
        "execution_git_untracked_files_checked": True,
        "execution_git_status_command": "git status --porcelain --untracked-files=all",
    }
    for key, expected in expected_execution_values.items():
        if set(hidden[key]) != {expected}:
            raise ValueError(f"hidden provenance has inconsistent {key}")
    if not all(str(value).strip() for value in hidden["execution_git_repo"]):
        raise ValueError("hidden provenance has an empty execution Git repository")
    observed_git_shas = set(hidden["execution_git_sha"].astype(str))
    if len(observed_git_shas) != 1 or _GIT_SHA_PATTERN.fullmatch(
        next(iter(observed_git_shas), "")
    ) is None:
        raise ValueError("hidden provenance has an invalid execution Git SHA")
    for key in ("execution_git_operator_expected_sha", "execution_git_launch_sha"):
        if set(hidden[key].astype(str)) != observed_git_shas:
            raise ValueError(f"hidden provenance has inconsistent {key}")
    if set(hidden["producer_git_sha"].astype(str)) != observed_git_shas:
        raise ValueError("hidden provenance has inconsistent producer_git_sha")


def _write_table(frame: pd.DataFrame, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".csv":
        frame.to_csv(path, index=False, float_format="%.17g")
    elif path.suffix.lower() == ".parquet":
        frame.to_parquet(path, compression="zstd", index=False)
    else:
        raise ValueError(f"unsupported output table format: {path}")
    return path


def write_s63_rank_one_candidates(
    candidates: pd.DataFrame,
    summary: Mapping[str, Any],
    *,
    out_path: Path,
    producer_git_sha: str,
) -> dict[str, Any]:
    """Failure-atomically publish a fresh candidate table/summary pair."""

    output = Path(out_path).expanduser().resolve(strict=False)
    summary_path = output.with_suffix(".summary.json")
    for path in (output, summary_path):
        if path.exists():
            raise FileExistsError(f"candidate output must be fresh: {path}")
    claims = _claim_publication_targets((output, summary_path))
    token = f"{os.getpid()}-{uuid.uuid4().hex}"
    staged_output = output.parent / f".{output.stem}.tmp-{token}{output.suffix}"
    staged_summary = summary_path.parent / f".{summary_path.name}.tmp-{token}"
    committed: list[Path] = []
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_table(candidates, staged_output)
        payload = dict(summary)
        payload["producer_git_sha"] = validate_git_sha(
            producer_git_sha, context="candidate producer_git_sha"
        )
        payload["candidate_table_sha256"] = file_sha256(staged_output)
        payload["outputs"] = {
            "candidate_table": str(output),
            "candidate_summary": str(summary_path),
        }
        staged_summary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        staged_summary_sha256 = file_sha256(staged_summary)
        for path in (output, summary_path):
            if path.exists():
                raise FileExistsError(f"candidate output appeared before publication: {path}")
        os.replace(staged_output, output)
        committed.append(output)
        os.replace(staged_summary, summary_path)
        committed.append(summary_path)
        payload["candidate_summary_sha256"] = staged_summary_sha256
        return payload
    except Exception:
        for path in reversed(committed):
            path.unlink(missing_ok=True)
        staged_output.unlink(missing_ok=True)
        staged_summary.unlink(missing_ok=True)
        raise
    finally:
        _release_publication_claims(claims)


def _replace_directory(source: Path, target: Path) -> None:
    """Commit one fully staged directory; split out for rollback testing."""

    source.replace(target)


def _assert_non_nested_fresh_directories(public_dir: Path, hidden_dir: Path) -> None:
    if public_dir == hidden_dir:
        raise ValueError("public and hidden queue directories must be distinct")
    if public_dir in hidden_dir.parents or hidden_dir in public_dir.parents:
        raise ValueError("public and hidden queue directories must not be nested")
    for name, path in (("public", public_dir), ("hidden", hidden_dir)):
        if path.exists():
            raise FileExistsError(f"{name} queue directory must be a fresh path: {path}")


def _claim_publication_targets(targets: Sequence[Path]) -> list[tuple[Path, int]]:
    claims: list[tuple[Path, int]] = []
    try:
        for target in sorted(targets, key=str):
            target.parent.mkdir(parents=True, exist_ok=True)
            claim = target.parent / f".{target.name}.publish.claim"
            descriptor = os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
            claims.append((claim, descriptor))
        return claims
    except Exception:
        for claim, descriptor in claims:
            os.close(descriptor)
            claim.unlink(missing_ok=True)
        raise


def _release_publication_claims(claims: Sequence[tuple[Path, int]]) -> None:
    for claim, descriptor in claims:
        os.close(descriptor)
        claim.unlink(missing_ok=True)


def _assert_public_summary_safe(payload: Mapping[str, Any]) -> None:
    forbidden_key_fragments = ("score", "bucket", "cohort", "hidden")

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            for key, nested in value.items():
                normalized = str(key).lower()
                if any(fragment in normalized for fragment in forbidden_key_fragments):
                    raise ValueError(
                        f"public queue summary exposes forbidden field: {key}"
                    )
                visit(nested)
        elif isinstance(value, list):
            for nested in value:
                visit(nested)

    visit(payload)


def _write_private_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite completion marker: {path}")
    temporary = path.parent / f".{path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(
                json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, 0o600)
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        temporary.unlink(missing_ok=True)


def verify_s63_prospective_review_bundle(
    *,
    public_dir: Path,
    hidden_dir: Path,
) -> dict[str, Any]:
    """Require the final cross-directory completion marker and exact hashes."""

    public_dir = Path(public_dir).expanduser().resolve(strict=True)
    hidden_dir = Path(hidden_dir).expanduser().resolve(strict=True)
    for name, directory in (("public", public_dir), ("hidden", hidden_dir)):
        if not directory.is_dir():
            raise ValueError(f"{name} queue bundle path is not a directory")
        if directory.stat().st_mode & 0o777 != 0o700:
            raise PermissionError(f"{name} queue bundle directory must be mode 0700")
    public_marker = public_dir / QUEUE_BUNDLE_COMPLETION_FILENAME
    hidden_marker = hidden_dir / QUEUE_BUNDLE_COMPLETION_FILENAME
    if not public_marker.is_file() or not hidden_marker.is_file():
        raise FileNotFoundError("queue bundle completion marker is missing")
    try:
        public_payload = json.loads(public_marker.read_bytes())
        hidden_payload = json.loads(hidden_marker.read_bytes())
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("queue bundle completion marker is invalid JSON") from exc
    if not isinstance(public_payload, dict) or not isinstance(hidden_payload, dict):
        raise ValueError("queue bundle completion markers must be objects")
    expected_marker_fields = {
        "schema_version": QUEUE_BUNDLE_COMPLETION_SCHEMA,
        "passed": True,
        "sector": S63_SECTOR,
        "publication_complete": True,
    }
    for name, payload, visibility in (
        ("public", public_payload, "reviewer_safe_unpublished"),
        ("hidden", hidden_payload, "private_provenance"),
    ):
        for key, expected in expected_marker_fields.items():
            if payload.get(key) != expected:
                raise ValueError(
                    f"{name} queue completion marker has incompatible {key}"
                )
        if payload.get("visibility") != visibility:
            raise ValueError(
                f"{name} queue completion marker has incompatible visibility"
            )
    bundle_id = str(public_payload.get("bundle_id", ""))
    if not re.fullmatch(r"[0-9a-f]{32}", bundle_id):
        raise ValueError("public queue completion marker has invalid bundle_id")
    if hidden_payload.get("bundle_id") != bundle_id:
        raise ValueError("queue completion marker bundle IDs differ")
    artifact_paths = {
        "public_review_queue": public_dir / "review_queue_1100.csv",
        "public_summary": public_dir / "public_summary.json",
        "hidden_selection_provenance": (
            hidden_dir / "hidden_selection_provenance.parquet"
        ),
        "hidden_summary": hidden_dir / "hidden_summary.json",
    }
    public_artifact_names = {"public_review_queue", "public_summary"}
    hidden_artifact_names = set(artifact_paths)
    public_records = _mapping(
        public_payload.get("artifacts"), context="public queue completion artifacts"
    )
    hidden_records = _mapping(
        hidden_payload.get("artifacts"), context="hidden queue completion artifacts"
    )
    if set(public_records) != public_artifact_names:
        raise ValueError("public completion marker contains non-public artifacts")
    if set(hidden_records) != hidden_artifact_names:
        raise ValueError("hidden completion marker does not bind the full bundle")

    def verify_records(
        records: Mapping[str, Any], names: set[str], *, context: str
    ) -> None:
        for artifact_name in sorted(names):
            path = artifact_paths[artifact_name]
            if not path.is_file():
                raise FileNotFoundError(
                    f"queue bundle artifact is missing: {artifact_name}"
                )
            if path.stat().st_mode & 0o777 != 0o600:
                raise PermissionError(
                    f"queue bundle artifact {artifact_name} must be mode 0600"
                )
            record = _mapping(
                records.get(artifact_name),
                context=f"{context} completion {artifact_name}",
            )
            if record.get("filename") != path.name:
                raise ValueError(
                    f"{context} completion filename drifted for {artifact_name}"
                )
            if validate_sha256(
                record.get("sha256"),
                context=f"{context} completion {artifact_name}",
            ) != file_sha256(path):
                raise ValueError(
                    f"{context} completion hash drifted for {artifact_name}"
                )
            if record.get("size_bytes") != path.stat().st_size:
                raise ValueError(
                    f"{context} completion size drifted for {artifact_name}"
                )

    verify_records(public_records, public_artifact_names, context="public")
    verify_records(hidden_records, hidden_artifact_names, context="hidden")
    for marker_path in (public_marker, hidden_marker):
        if marker_path.stat().st_mode & 0o777 != 0o600:
            raise PermissionError("queue completion markers must be mode 0600")
    return {
        "schema_version": QUEUE_BUNDLE_COMPLETION_SCHEMA,
        "passed": True,
        "bundle_id": bundle_id,
        "public_completion_marker_sha256": file_sha256(public_marker),
        "hidden_completion_marker_sha256": file_sha256(hidden_marker),
        "artifacts": dict(hidden_records),
    }


def write_s63_prospective_review_queue(
    queue: pd.DataFrame,
    hidden: pd.DataFrame,
    public_summary: Mapping[str, Any],
    hidden_summary: Mapping[str, Any],
    *,
    public_dir: Path,
    hidden_dir: Path,
) -> dict[str, Any]:
    """Failure-atomically publish isolated public and hidden directories."""

    public_dir = Path(public_dir).expanduser().resolve(strict=False)
    hidden_dir = Path(hidden_dir).expanduser().resolve(strict=False)
    if tuple(queue.columns) != PUBLIC_REVIEW_QUEUE_COLUMNS:
        unexpected = sorted(set(queue.columns) - set(PUBLIC_REVIEW_QUEUE_COLUMNS))
        missing = sorted(set(PUBLIC_REVIEW_QUEUE_COLUMNS) - set(queue.columns))
        raise ValueError(
            "public queue violates the exact writer allowlist; "
            f"unexpected={unexpected}, missing={missing}"
        )
    _assert_non_nested_fresh_directories(public_dir, hidden_dir)
    _assert_public_summary_safe(public_summary)
    claims = _claim_publication_targets((public_dir, hidden_dir))
    token = f"{os.getpid()}-{uuid.uuid4().hex}"
    public_staging = public_dir.parent / f".{public_dir.name}.tmp-{token}"
    hidden_staging = hidden_dir.parent / f".{hidden_dir.name}.tmp-{token}"
    committed: list[Path] = []
    try:
        _assert_non_nested_fresh_directories(public_dir, hidden_dir)
        public_staging.mkdir(mode=0o700)
        hidden_staging.mkdir(mode=0o700)
        os.chmod(public_staging, 0o700)
        os.chmod(hidden_staging, 0o700)

        staged_queue = _write_table(queue, public_staging / "review_queue_1100.csv")
        os.chmod(staged_queue, 0o600)
        queue_record = {
            "path": str(public_dir / staged_queue.name),
            "sha256": file_sha256(staged_queue),
            "size_bytes": staged_queue.stat().st_size,
        }
        public_payload = dict(public_summary)
        public_payload["publication_contract"] = (
            "private_split_bundle_requires_matching_completion_bundle_ids"
        )
        public_payload["outputs"] = {"review_queue": queue_record}
        _assert_public_summary_safe(public_payload)
        if str(hidden_dir) in json.dumps(public_payload, sort_keys=True):
            raise ValueError("public queue summary exposes the hidden directory path")
        staged_public_summary = public_staging / "public_summary.json"
        staged_public_summary.write_text(
            json.dumps(public_payload, indent=2, sort_keys=True, allow_nan=False)
            + "\n",
            encoding="utf-8",
        )
        os.chmod(staged_public_summary, 0o600)

        staged_hidden = _write_table(
            hidden,
            hidden_staging / "hidden_selection_provenance.parquet",
        )
        os.chmod(staged_hidden, 0o600)
        hidden_payload = dict(hidden_summary)
        hidden_payload["publication_contract"] = (
            "private_split_bundle_requires_matching_completion_bundle_ids"
        )
        hidden_payload["outputs"] = {
            "public_review_queue": queue_record,
            "public_summary": {
                "path": str(public_dir / staged_public_summary.name),
                "sha256": file_sha256(staged_public_summary),
                "size_bytes": staged_public_summary.stat().st_size,
            },
            "hidden_selection_provenance": {
                "path": str(hidden_dir / staged_hidden.name),
                "sha256": file_sha256(staged_hidden),
                "size_bytes": staged_hidden.stat().st_size,
            },
        }
        staged_hidden_summary = hidden_staging / "hidden_summary.json"
        staged_hidden_summary.write_text(
            json.dumps(hidden_payload, indent=2, sort_keys=True, allow_nan=False)
            + "\n",
            encoding="utf-8",
        )
        os.chmod(staged_hidden_summary, 0o600)

        _assert_non_nested_fresh_directories(public_dir, hidden_dir)
        _replace_directory(public_staging, public_dir)
        committed.append(public_dir)
        _replace_directory(hidden_staging, hidden_dir)
        committed.append(hidden_dir)
        os.chmod(public_dir, 0o700)
        os.chmod(hidden_dir, 0o700)
        final_artifacts = {
            "public_review_queue": public_dir / "review_queue_1100.csv",
            "public_summary": public_dir / "public_summary.json",
            "hidden_selection_provenance": (
                hidden_dir / "hidden_selection_provenance.parquet"
            ),
            "hidden_summary": hidden_dir / "hidden_summary.json",
        }
        bundle_id = uuid.uuid4().hex
        common_completion = {
            "schema_version": QUEUE_BUNDLE_COMPLETION_SCHEMA,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "passed": True,
            "sector": S63_SECTOR,
            "publication_complete": True,
            "bundle_id": bundle_id,
        }
        public_completion_payload = {
            **common_completion,
            "visibility": "reviewer_safe_unpublished",
            "private_completion_marker_required": True,
            "artifacts": {
                name: {
                    "filename": path.name,
                    "sha256": file_sha256(path),
                    "size_bytes": path.stat().st_size,
                }
                for name, path in final_artifacts.items()
                if name in {"public_review_queue", "public_summary"}
            },
        }
        hidden_completion_payload = {
            **common_completion,
            "visibility": "private_provenance",
            "artifacts": {
                name: {
                    "filename": path.name,
                    "sha256": file_sha256(path),
                    "size_bytes": path.stat().st_size,
                }
                for name, path in final_artifacts.items()
            },
        }
        # The public marker is the final commit signal. A process or node crash
        # before it appears leaves a bundle that every validator rejects.
        _write_private_json_atomic(
            hidden_dir / QUEUE_BUNDLE_COMPLETION_FILENAME,
            hidden_completion_payload,
        )
        _write_private_json_atomic(
            public_dir / QUEUE_BUNDLE_COMPLETION_FILENAME,
            public_completion_payload,
        )
        completion = verify_s63_prospective_review_bundle(
            public_dir=public_dir,
            hidden_dir=hidden_dir,
        )
        return {
            "schema_version": "twirl_teacher_v3_s63_split_queue_publication_v1",
            "passed": True,
            "public": {
                "directory": str(public_dir),
                "summary": str(public_dir / "public_summary.json"),
                "summary_sha256": file_sha256(public_dir / "public_summary.json"),
            },
            "hidden": {
                "directory": str(hidden_dir),
                "summary": str(hidden_dir / "hidden_summary.json"),
                "summary_sha256": file_sha256(hidden_dir / "hidden_summary.json"),
            },
            "completion": completion,
            "rollback_required": False,
        }
    except Exception:
        # Targets were proven absent and claimed before publication, so only
        # directories committed by this transaction can be removed here.
        for committed_dir in reversed(committed):
            shutil.rmtree(committed_dir, ignore_errors=True)
        shutil.rmtree(public_staging, ignore_errors=True)
        shutil.rmtree(hidden_staging, ignore_errors=True)
        raise
    finally:
        _release_publication_claims(claims)


__all__ = [
    "CANDIDATE_CONTRACT_VERSION",
    "COHORT_CONTRACT_VERSION",
    "PRIMARY_COHORT",
    "PRIMARY_QUOTAS",
    "PROSPECTIVE_PLAN_SCHEMA",
    "PROSPECTIVE_USE_SCOPE",
    "PUBLIC_REVIEW_QUEUE_COLUMNS",
    "ProspectiveQuotas",
    "QUEUE_POLICY_VERSION",
    "REPEATED_HOST_COHORT",
    "REPEATED_HOST_QUOTAS",
    "S63_SECTOR",
    "S63_SCORE_POLICY",
    "S63_TARGET_METADATA_CONTRACT",
    "SELECTION_POLICY_SCHEMA",
    "TEACHER_V3_APERTURE_METADATA_FIELDS",
    "TEACHER_V3_CANDIDATE_METADATA_COLUMNS",
    "build_s63_prospective_review_queue",
    "build_s63_rank_one_candidates",
    "derive_s63_prospective_cohorts",
    "file_sha256",
    "load_prospective_plan",
    "load_selection_policy",
    "read_tic_inventory",
    "tic_inventory_sha256",
    "validate_s63_bls_summary",
    "validate_git_sha",
    "validate_queue_launch_bindings",
    "validate_selection_policy",
    "verify_s63_prospective_review_queue",
    "verify_s63_prospective_review_bundle",
    "write_s63_prospective_cohorts",
    "write_s63_prospective_review_queue",
    "write_s63_rank_one_candidates",
]
