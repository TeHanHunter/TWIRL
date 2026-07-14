"""Data roles and locked-evaluation policy for the S56 Teacher-v2 cycle.

Teacher v2 deliberately separates three kinds of supervision:

* fresh, BLS-recovered injections supervise compact-transit recognition and
  the harmonic head;
* clean real human labels supervise morphology and preserve/reject;
* injection truth and Franklin labels remain audit/provenance columns unless a
  later, explicit policy promotes them.

None of the helpers in this module select model inputs from truth or label
columns.  They build target/role manifests consumed by the training driver.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.adjudication_audit import (
    MORPHOLOGY_CLASSES_V1,
    TRANSIT_HARMONIC_FACTORS_V1,
    TRANSIT_HARMONIC_TARGETS_V1,
)
from twirl.vetting.harmonic_cnn import PRESERVE_CLASSES, build_grouped_test_and_cv_folds
from twirl.vetting.harmonic_inputs import native_group_path
from twirl.vetting.label_io import candidate_key as current_candidate_key
from twirl.vetting.recovery50_teacher import leakage_columns


TEACHER_V2_MODEL_VERSION = "s56_harmonic_cnn_teacher_v2"
TEACHER_V2_ROLE_POLICY = "s56_a2v1_teacher_v2_hybrid_roles_v1"
COMPACT_CLASSES: tuple[str, str] = ("background", "compact_transit")
FRANKLIN_ACTIVE_LABELS: frozenset[str] = frozenset(
    {
        "planet_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "wide_transit_like",
        "instrumental_or_systematic",
        "uncertain",
        "no_visible_signal",
        "centroid_contaminant",
        "skip",
    }
)

TEACHER_V2_FORBIDDEN_FEATURE_TOKENS: tuple[str, ...] = (
    "truth",
    "recovery",
    "label",
    "note",
    "adjudicat",
    "corrected_period",
    "effective_period",
    "refold",
    "period_factor",
    "source_kind",
    "source_bucket",
    "cohort",
    "selection_weight",
    "pseudo",
    "path",
    "h5",
    "role",
    "split",
    "fold",
    "injection_id",
    "candidate_key",
)


def _as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.fillna("").astype(str).str.strip().str.lower().isin(
        {"1", "1.0", "true", "t", "yes", "y"}
    )


def _stable_digest(value: str, *, seed: int) -> str:
    return hashlib.sha256(f"{seed}|{value}".encode("utf-8")).hexdigest()


def _stable_fold(value: str, *, seed: int, n_folds: int = 5) -> int:
    digest = hashlib.sha256(f"{seed}|{value}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % int(n_folds)


def _require_unique(frame: pd.DataFrame, column: str, *, name: str) -> None:
    if column not in frame:
        raise KeyError(f"{name} is missing {column}")
    text = frame[column].fillna("").astype(str)
    if text.eq("").any():
        raise ValueError(f"{name} contains empty {column} values")
    if text.duplicated().any():
        duplicate = text[text.duplicated(keep=False)].head(5).tolist()
        raise ValueError(f"{name} contains duplicate {column} values: {duplicate}")


def normalize_franklin_labels(labels: pd.DataFrame) -> pd.DataFrame:
    """Validate Franklin's handoff using ``candidate_key`` as identity.

    The legacy ``row_id`` is retained as provenance but is never used for a
    join.  This makes the normalization safe against the known queue-row drift.
    """

    work = labels.copy()
    _require_unique(work, "candidate_key", name="Franklin labels")
    required = {"tic", "sector", "label"}
    missing = sorted(required - set(work.columns))
    if missing:
        raise KeyError(f"Franklin labels are missing columns: {missing}")
    work["tic"] = pd.to_numeric(work["tic"], errors="coerce")
    work["sector"] = pd.to_numeric(work["sector"], errors="coerce")
    if work[["tic", "sector"]].isna().any().any():
        raise ValueError("Franklin labels contain invalid TIC or sector values")
    work["tic"] = work["tic"].astype(np.int64)
    work["sector"] = work["sector"].astype(np.int16)
    work["label"] = work["label"].fillna("").astype(str).str.strip()
    unknown = sorted(set(work["label"]) - set(FRANKLIN_ACTIVE_LABELS))
    if unknown:
        raise ValueError(f"Franklin labels contain unknown labels: {unknown}")
    work = work.rename(columns={"row_id": "legacy_row_id"})
    work.insert(0, "franklin_identity", "candidate_key:" + work["candidate_key"])
    work["teacher_v2_role"] = "franklin_audit_only"
    work["teacher_v2_training_include"] = False
    return work.reset_index(drop=True)


def join_franklin_labels_to_queue(
    labels: pd.DataFrame,
    queue: pd.DataFrame,
) -> pd.DataFrame:
    """Recover Franklin queue metadata by exact app-written candidate key."""

    normalized = normalize_franklin_labels(labels)
    required = {"review_id", "tic", "sector", "period_d", "t0_bjd", "duration_min"}
    missing = sorted(required - set(queue.columns))
    if missing:
        raise KeyError(f"Franklin source queue is missing columns: {missing}")
    source = queue.copy()
    for column in required:
        source[column] = source[column].fillna("").astype(str)
    source["_franklin_candidate_key"] = source.apply(
        lambda row: "|".join(
            str(row[column])
            for column in ("review_id", "tic", "sector", "period_d", "t0_bjd")
        ),
        axis=1,
    )
    _require_unique(
        source,
        "_franklin_candidate_key",
        name="Franklin source queue",
    )
    metadata_columns = [
        column
        for column in source.columns
        if column
        not in {
            "row_id",
            "candidate_key",
            "label",
            "label_source",
            "labeler",
            "notes",
            "updated_utc",
        }
    ]
    metadata = source[metadata_columns].rename(
        columns={
            "_franklin_candidate_key": "candidate_key",
            "review_id": "franklin_queue_review_id",
        }
    )
    joined = normalized.merge(
        metadata,
        on="candidate_key",
        how="left",
        validate="one_to_one",
        suffixes=("", "_queue"),
        indicator=True,
    )
    missing_identity = joined["_merge"].ne("both")
    if missing_identity.any():
        keys = joined.loc[missing_identity, "candidate_key"].head(5).tolist()
        raise ValueError(
            f"{int(missing_identity.sum())} Franklin labels do not match candidate_key; "
            f"first={keys}"
        )
    joined = joined.drop(columns=["_merge"])
    for column in ("period_d", "t0_bjd", "duration_min"):
        joined[column] = pd.to_numeric(joined[column], errors="raise")
    if "tic_queue" in joined:
        queue_tic = pd.to_numeric(joined["tic_queue"], errors="raise").astype(np.int64)
        if not queue_tic.equals(joined["tic"].astype(np.int64)):
            raise ValueError("Franklin candidate_key join changed TIC identity")
        joined = joined.drop(columns=["tic_queue"])
    if "sector_queue" in joined:
        queue_sector = pd.to_numeric(joined["sector_queue"], errors="raise").astype(np.int16)
        if not queue_sector.equals(joined["sector"].astype(np.int16)):
            raise ValueError("Franklin candidate_key join changed sector identity")
        joined = joined.drop(columns=["sector_queue"])
    joined["teacher_v2_training_include"] = False
    joined["teacher_v2_role"] = "franklin_audit_only"
    return joined.reset_index(drop=True)


def build_franklin_a2v1_rereview_queue(
    franklin: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    seed: int = 56,
    n_controls: int = 100,
    expected_signal_rows: int = 34,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Build a blinded current-A2v1 queue for Franklin-label transfer QA."""

    current = candidates.sort_values(
        ["tic", "rep_peak_rank"], kind="stable"
    ).drop_duplicates("tic", keep="first")
    current = current.copy()
    current["tic"] = pd.to_numeric(current["tic"], errors="raise").astype(np.int64)
    prior = franklin.rename(
        columns={"candidate_key": "franklin_candidate_key"}
    ).copy()
    prior["tic"] = pd.to_numeric(prior["tic"], errors="raise").astype(np.int64)
    joined = prior.merge(
        current,
        on="tic",
        how="left",
        validate="many_to_one",
        suffixes=("_franklin", ""),
        indicator=True,
    )
    signal_labels = {
        "planet_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "wide_transit_like",
    }
    control_labels = {"instrumental_or_systematic", "uncertain", "no_visible_signal"}
    signals = joined.loc[joined["label"].isin(signal_labels)].copy()
    if len(signals) != int(expected_signal_rows):
        raise ValueError(
            f"Franklin re-review expected {expected_signal_rows} signal rows, got {len(signals)}"
        )
    missing_signals = signals["_merge"].ne("both")
    if missing_signals.any():
        raise ValueError(
            f"{int(missing_signals.sum())} Franklin signal rows have no current A2v1 candidate"
        )
    controls = joined.loc[
        joined["label"].isin(control_labels) & joined["_merge"].eq("both")
    ].copy()
    if len(controls) < int(n_controls):
        raise ValueError(
            f"Franklin re-review has only {len(controls)} eligible controls; "
            f"requested={n_controls}"
        )
    controls = controls.sample(n=int(n_controls), random_state=int(seed))
    selected = pd.concat([signals, controls], ignore_index=True, sort=False)
    selected = selected.sample(frac=1.0, random_state=int(seed)).reset_index(drop=True)
    selected["source_a2v1_review_id"] = selected["review_id"].astype(str)
    selected["review_id"] = selected["franklin_identity"].map(
        lambda value: "s0056-a2v1-franklin-rereview-"
        + _stable_digest(str(value), seed=seed)[:20]
    )
    selected.insert(0, "row_id", np.arange(len(selected), dtype=np.int64))
    selected["source_kind"] = "real_candidate"
    selected["source_bucket"] = ""
    selected["candidate_key"] = selected.apply(current_candidate_key, axis=1)
    selected["vet_sheet_version"] = "S56-ADP-HV1"
    selected["twirl_vet_sheet_name"] = (
        selected["review_id"].str.replace(":", "_", regex=False)
        + "_twirl_twoap_current_a2v1_adp.png"
    )
    selected["twirl_vet_sheet_pdf_name"] = ""
    private_columns = [
        "row_id",
        "review_id",
        "candidate_key",
        "franklin_identity",
        "franklin_candidate_key",
        "franklin_queue_review_id",
        "label",
        "label_source",
        "labeler",
        "notes",
        "updated_utc",
        "source_a2v1_review_id",
        "tic",
        "period_d",
        "t0_bjd",
        "duration_min",
        "period_d_franklin",
        "t0_bjd_franklin",
        "duration_min_franklin",
    ]
    private_columns = [column for column in private_columns if column in selected]
    private = selected[private_columns].copy()
    hidden_tokens = ("label", "franklin", "legacy", "prior", "truth", "recovery")
    public_columns = [
        column
        for column in (
            "row_id",
            "review_id",
            "tic",
            "sector",
            "cam",
            "ccd",
            "tmag",
            "source_kind",
            "source_bucket",
            "period_d",
            "t0_bjd",
            "duration_min",
            "depth",
            "depth_snr",
            "sde_max",
            "rep_aperture",
            "rep_peak_rank",
            "n_apertures_agree",
            "apertures_agree",
            "aperture_period_rel_delta",
            "aperture_disagreement_flag",
            "candidate_key",
            "vet_sheet_version",
            "twirl_vet_sheet_name",
            "twirl_vet_sheet_pdf_name",
        )
        if column in selected
    ]
    public = selected[public_columns].copy()
    if any(any(token in column.lower() for token in hidden_tokens) for column in public):
        raise RuntimeError("Franklin public re-review queue exposes hidden provenance")
    summary = {
        "n_rows": int(len(public)),
        "n_signal_rows": int(len(signals)),
        "n_control_rows": int(len(controls)),
        "n_unique_tics": int(public["tic"].nunique()),
        "seed": int(seed),
        "prior_labels_hidden": True,
        "vet_sheet_version": "S56-ADP-HV1",
    }
    return public, private, summary


def build_franklin_current_a2v1_audit_candidates(
    franklin: pd.DataFrame,
    candidates: pd.DataFrame,
    compatibility: pd.DataFrame,
    *,
    seed: int = 56,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Attach every Franklin row to the current A2v1 top peak for score audit."""

    current = candidates.sort_values(
        ["tic", "rep_peak_rank"], kind="stable"
    ).drop_duplicates("tic", keep="first")
    current = current.copy()
    current["tic"] = pd.to_numeric(current["tic"], errors="raise").astype(np.int64)
    prior_columns = [
        column
        for column in (
            "franklin_identity",
            "tic",
            "label",
            "label_source",
            "labeler",
            "notes",
            "updated_utc",
            "candidate_key",
            "franklin_queue_review_id",
        )
        if column in franklin
    ]
    prior = franklin[prior_columns].rename(
        columns={
            "label": "legacy_human_label",
            "label_source": "legacy_label_source",
            "labeler": "legacy_labeler",
            "notes": "legacy_notes",
            "updated_utc": "legacy_updated_utc",
            "candidate_key": "franklin_candidate_key",
        }
    )
    prior["tic"] = pd.to_numeric(prior["tic"], errors="raise").astype(np.int64)
    joined = prior.merge(current, on="tic", how="left", validate="many_to_one", indicator=True)
    missing = joined["_merge"].ne("both")
    audit = joined.loc[~missing].drop(columns=["_merge"]).copy().reset_index(drop=True)
    audit["source_a2v1_review_id"] = audit["review_id"].astype(str)
    audit["review_id"] = audit["franklin_identity"].map(
        lambda value: "s0056-a2v1-franklin-audit-"
        + _stable_digest(str(value), seed=seed)[:20]
    )
    compatibility_key = (
        "source_uid" if "source_uid" in compatibility else "franklin_identity"
    )
    transfer = compatibility[
        [compatibility_key, "a2v1_transfer_ok", "a2v1_match_review_id"]
    ].rename(columns={compatibility_key: "franklin_identity"})
    audit = audit.merge(transfer, on="franklin_identity", how="left", validate="one_to_one")
    audit["a2v1_transfer_ok"] = _as_bool(audit["a2v1_transfer_ok"])
    audit["teacher_v2_role"] = "franklin_current_a2v1_score_audit_only"
    audit["teacher_v2_training_include"] = False
    audit["is_injected_row"] = False
    summary = {
        "n_franklin_rows": int(len(franklin)),
        "n_current_a2v1_rows": int(len(audit)),
        "n_missing_current_a2v1": int(missing.sum()),
        "n_ephemeris_compatible": int(audit["a2v1_transfer_ok"].sum()),
        "training_include": False,
        "legacy_labels_are_descriptive_not_active_targets": True,
    }
    return audit, summary


def assign_s56_injection_roles(
    schedule: pd.DataFrame,
    *,
    seed: int = 560202,
) -> pd.DataFrame:
    """Assign exactly six development and two locked injections per cell."""

    work = schedule.copy()
    required = {"injection_id", "tic", "grid_cell_id"}
    missing = sorted(required - set(work.columns))
    if missing:
        raise KeyError(f"S56 injection schedule is missing columns: {missing}")
    _require_unique(work, "injection_id", name="S56 injection schedule")
    work["tic"] = pd.to_numeric(work["tic"], errors="raise").astype(np.int64)
    if work["tic"].duplicated().any():
        raise ValueError("Teacher-v2 injections require one unique host TIC each")
    support = work.groupby("grid_cell_id", dropna=False).size()
    if len(support) != 2500 or not support.eq(8).all():
        raise ValueError(
            "S56 Teacher-v2 split requires exactly 2,500 cells with eight injections each"
        )
    work["_split_order"] = work["injection_id"].astype(str).map(
        lambda value: _stable_digest(value, seed=seed)
    )
    work["_cell_rank"] = (
        work.sort_values(["grid_cell_id", "_split_order"], kind="stable")
        .groupby("grid_cell_id", sort=False)
        .cumcount()
        .reindex(work.index)
    )
    work["teacher_v2_partition"] = np.where(
        work["_cell_rank"].lt(2), "locked_holdout", "development"
    )
    work["fixed_split"] = np.where(
        work["teacher_v2_partition"].eq("locked_holdout"), "test", "development"
    )
    work["cv_fold"] = [
        -1
        if partition == "locked_holdout"
        else _stable_fold(str(injection_id), seed=seed + 1)
        for injection_id, partition in zip(
            work["injection_id"], work["teacher_v2_partition"]
        )
    ]
    cell_partition = pd.crosstab(
        work["grid_cell_id"], work["teacher_v2_partition"]
    ).reindex(columns=["development", "locked_holdout"], fill_value=0)
    if not (
        cell_partition["development"].eq(6).all()
        and cell_partition["locked_holdout"].eq(2).all()
    ):
        raise RuntimeError("S56 injection role assignment violated the 6/2 cell contract")
    return work.drop(columns=["_split_order", "_cell_rank"]).reset_index(drop=True)


def _nearest_harmonic_index(value: Any, *, tolerance: float = 0.02) -> int:
    try:
        factor = float(value)
    except (TypeError, ValueError):
        return -1
    factors = np.asarray(TRANSIT_HARMONIC_FACTORS_V1, dtype=float)
    error = np.abs(factor / factors - 1.0)
    index = int(np.argmin(error))
    return index if np.isfinite(error[index]) and error[index] <= tolerance else -1


def build_s56_injection_training_rows(
    schedule: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    native_h5: Path,
    seed: int = 560202,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build compact/harmonic rows from genuinely recovered S56 injections.

    Every truth-matched top-five candidate is positive.  Every unmatched peak
    from the same recovered injected light curve is a hard negative.  A paired
    pre-injection view at each positive ephemeris is a shape-only negative.
    Injections with no truth-matched top-five candidate produce no training row.
    """

    roles = assign_s56_injection_roles(schedule, seed=seed)
    frame = candidates.copy()
    required = {
        "injection_id",
        "review_id",
        "tic",
        "period_d",
        "t0_bjd",
        "duration_min",
        "rep_peak_rank",
        "is_injected_signal_peak",
        "nearest_harmonic_factor",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError(f"S56 candidate table is missing columns: {missing}")
    _require_unique(frame, "review_id", name="S56 top-five candidates")
    match = _as_bool(frame["is_injected_signal_peak"])
    recovered_ids = frozenset(frame.loc[match, "injection_id"].astype(str))
    frame = frame.loc[frame["injection_id"].astype(str).isin(recovered_ids)].copy()
    frame["_truth_match"] = _as_bool(frame["is_injected_signal_peak"])
    frame = frame.merge(
        roles[
            [
                "injection_id",
                "grid_cell_id",
                "teacher_v2_partition",
                "fixed_split",
                "cv_fold",
            ]
        ],
        on="injection_id",
        how="left",
        validate="many_to_one",
    )
    if frame["teacher_v2_partition"].isna().any():
        raise ValueError("candidate table contains injections absent from the S56 schedule")
    frame["input_variant"] = "injected"
    frame["teacher_v2_role"] = np.where(
        frame["_truth_match"], "compact_positive", "same_lc_unmatched_peak"
    )
    frame["compact_target_v2"] = np.where(
        frame["_truth_match"], COMPACT_CLASSES[1], COMPACT_CLASSES[0]
    )
    frame["compact_target_index"] = frame["_truth_match"].astype(np.int8)
    frame["harmonic_target_index"] = -1
    positive_index = frame.index[frame["_truth_match"]]
    harmonic_index = frame.loc[positive_index, "nearest_harmonic_factor"].map(
        _nearest_harmonic_index
    )
    if harmonic_index.lt(0).any():
        bad = frame.loc[harmonic_index.index[harmonic_index.lt(0)], "review_id"].head(5)
        raise ValueError(f"truth-matched candidates have unsupported harmonic factors: {bad.tolist()}")
    frame.loc[positive_index, "harmonic_target_index"] = harmonic_index.astype(np.int8)
    frame["harmonic_target_v2"] = ""
    frame.loc[positive_index, "harmonic_target_v2"] = [
        TRANSIT_HARMONIC_TARGETS_V1[int(index)] for index in harmonic_index
    ]

    positives = frame.loc[frame["_truth_match"]].copy()
    paired = positives.copy()
    paired["review_id"] = paired["review_id"].astype(str) + ":paired_original"
    paired["input_variant"] = "paired_original"
    paired["teacher_v2_role"] = "paired_pre_injection"
    paired["compact_target_v2"] = COMPACT_CLASSES[0]
    paired["compact_target_index"] = np.int8(0)
    paired["harmonic_target_index"] = np.int8(-1)
    paired["harmonic_target_v2"] = ""
    training = pd.concat([frame, paired], ignore_index=True, sort=False)
    if training["review_id"].duplicated().any():
        raise RuntimeError("Teacher-v2 injection training review IDs are not unique")
    training["morphology_target_v1"] = ""
    training["morphology_target_index"] = np.int8(-1)
    training["preserve_target_v1"] = ""
    training["preserve_target_index"] = np.int8(-1)
    training["morphology_weight"] = np.float32(0.0)
    training["preserve_weight"] = np.float32(0.0)
    training["harmonic_weight"] = np.float32(0.0)
    training["compact_weight"] = np.float32(0.0)
    training["pretrain_target"] = -1
    training["native_h5_path"] = str(Path(native_h5))
    training["native_group_path"] = [
        native_group_path(row) for row in training.to_dict("records")
    ]
    training["teacher_v2_role_policy"] = TEACHER_V2_ROLE_POLICY
    training["teacher_v2_training_include"] = True
    training["is_injected_row"] = True
    training["source_kind"] = "fresh_s56_a2v1_injection"
    training = training.drop(columns=["_truth_match"])

    role_manifest = roles.copy()
    role_manifest["bls_top5_recovered"] = role_manifest["injection_id"].astype(str).isin(
        recovered_ids
    )
    role_manifest["teacher_v2_classifier_include"] = role_manifest[
        "bls_top5_recovered"
    ]
    role_manifest["teacher_v2_role_policy"] = TEACHER_V2_ROLE_POLICY
    return training.reset_index(drop=True), role_manifest.reset_index(drop=True)


def build_real_human_training_rows(
    human_rows: pd.DataFrame,
    *,
    native_h5: Path,
    compatibility: pd.DataFrame | None = None,
    seed: int = 56,
) -> pd.DataFrame:
    """Build real-only morphology/preserve rows after A2v1 compatibility QA."""

    work = human_rows.copy()
    injected = (
        _as_bool(work["is_injected_row"])
        if "is_injected_row" in work
        else work.get("source_kind", pd.Series("", index=work.index))
        .fillna("")
        .astype(str)
        .str.contains("inject", case=False)
    )
    work = work.loc[~injected].copy()
    if compatibility is not None:
        key = "source_uid" if "source_uid" in work and "source_uid" in compatibility else "review_id"
        if key not in work or key not in compatibility:
            raise KeyError("A2v1 compatibility table has no usable human-row identity")
        transfer = compatibility[[key, "a2v1_transfer_ok"]].drop_duplicates(key)
        work = work.drop(columns=["a2v1_transfer_ok"], errors="ignore")
        work = work.merge(transfer, on=key, how="left", validate="one_to_one")
        work = work.loc[_as_bool(work["a2v1_transfer_ok"])].copy()
    label = work.get("human_label", pd.Series("", index=work.index)).fillna("").astype(str)
    work = work.loc[~label.isin({"", "skip"})].copy()
    unresolved = work.get("ephemeris_status", pd.Series("", index=work.index)).fillna("").astype(str).eq(
        "unresolved"
    )
    signal = work["human_label"].isin(
        {"planet_like", "eclipsing_binary_or_pceb", "wide_transit_like"}
    )
    work = work.loc[~(unresolved & signal)].copy()

    morphology = work.get("morphology_target_v1", pd.Series("", index=work.index)).fillna("").astype(str)
    preserve = work.get("preserve_target_v1", pd.Series("", index=work.index)).fillna("").astype(str)
    harmonic = work.get("harmonic_target_v1", pd.Series("", index=work.index)).fillna("").astype(str)
    work["morphology_target_index"] = [
        MORPHOLOGY_CLASSES_V1.index(value) if value in MORPHOLOGY_CLASSES_V1 else -1
        for value in morphology
    ]
    work["preserve_target_index"] = [
        PRESERVE_CLASSES.index(value) if value in PRESERVE_CLASSES else -1 for value in preserve
    ]
    harmonic_include = (
        _as_bool(work["harmonic_include_v1"])
        if "harmonic_include_v1" in work
        else harmonic.ne("")
    )
    work["harmonic_target_index"] = [
        TRANSIT_HARMONIC_TARGETS_V1.index(value)
        if include and value in TRANSIT_HARMONIC_TARGETS_V1
        else -1
        for value, include in zip(harmonic, harmonic_include)
    ]
    active = (
        work["morphology_target_index"].ge(0)
        | work["preserve_target_index"].ge(0)
        | work["harmonic_target_index"].ge(0)
    )
    work = work.loc[active].copy().reset_index(drop=True)
    if "fixed_split" not in work or not set(work["fixed_split"].astype(str)).issubset(
        {"development", "test"}
    ):
        split = build_grouped_test_and_cv_folds(work, seed=seed)
        work["fixed_split"] = split["fixed_split"].to_numpy()
        work["cv_fold"] = split["cv_fold"].to_numpy()
    elif "cv_fold" not in work:
        raise KeyError("human rows with fixed_split must also provide cv_fold")
    work["teacher_v2_partition"] = np.where(
        work["fixed_split"].eq("test"), "locked_holdout", "development"
    )
    work["compact_target_v2"] = ""
    work["compact_target_index"] = np.int8(-1)
    work["input_variant"] = "observed"
    work["morphology_weight"] = np.float32(0.0)
    work["preserve_weight"] = np.float32(0.0)
    work["harmonic_weight"] = np.float32(0.0)
    work["compact_weight"] = np.float32(0.0)
    work["pretrain_target"] = -1
    work["native_h5_path"] = str(Path(native_h5))
    work["native_group_path"] = [native_group_path(row) for row in work.to_dict("records")]
    work["teacher_v2_role"] = "real_human_morphology"
    work["teacher_v2_role_policy"] = TEACHER_V2_ROLE_POLICY
    work["teacher_v2_training_include"] = True
    work["is_injected_row"] = False
    return work.reset_index(drop=True)


def normalize_real_adp_candidates(
    peaks: pd.DataFrame,
    *,
    small_peaks_per_tic: int = 5,
) -> pd.DataFrame:
    """Normalize ADP-small BLS peaks with ADP-primary rank-one context."""

    frame = peaks.copy()
    required = {
        "tic",
        "sector",
        "aperture",
        "status",
        "peak_rank",
        "period_d",
        "t0_bjd",
        "duration_min",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError(f"real ADP BLS table is missing columns: {missing}")
    aperture = frame["aperture"].fillna("").astype(str)
    forbidden = sorted(
        set(aperture.dropna()) - {"DET_FLUX_ADP_SML", "DET_FLUX_ADP"}
    )
    if forbidden:
        raise ValueError(f"real BLS table contains non-ADP apertures: {forbidden}")
    rank = pd.to_numeric(frame["peak_rank"], errors="coerce")
    ok = frame["status"].fillna("").astype(str).eq("ok")
    small = frame.loc[
        ok
        & aperture.eq("DET_FLUX_ADP_SML")
        & rank.between(1, int(small_peaks_per_tic))
    ].copy()
    primary = frame.loc[
        ok & aperture.eq("DET_FLUX_ADP") & rank.eq(1)
    ].sort_values("tic", kind="stable").drop_duplicates("tic", keep="first")
    metric_columns = tuple(
        column
        for column in (
            "peak_rank",
            "period_d",
            "t0_bjd",
            "duration_min",
            "depth",
            "depth_snr",
            "sde",
            "log_power",
        )
        if column in frame
    )
    primary_context = primary[["tic", *metric_columns]].rename(
        columns={column: f"adp_{column}" for column in metric_columns}
    )
    out = small.merge(primary_context, on="tic", how="left", validate="many_to_one")
    for column in metric_columns:
        out[f"adp_sml_{column}"] = pd.to_numeric(out[column], errors="coerce")
    out["tic"] = pd.to_numeric(out["tic"], errors="raise").astype(np.int64)
    out["sector"] = pd.to_numeric(out["sector"], errors="raise").astype(np.int16)
    out["rep_peak_rank"] = pd.to_numeric(out["peak_rank"], errors="raise").astype(np.int16)
    out["rep_aperture"] = "DET_FLUX_ADP_SML"
    out["sde_max"] = pd.to_numeric(out.get("sde"), errors="coerce")
    out["anchor_aperture"] = "DET_FLUX_ADP_SML"
    out["anchor_period_d"] = pd.to_numeric(out["period_d"], errors="coerce")
    out["anchor_t0_bjd"] = pd.to_numeric(out["t0_bjd"], errors="coerce")
    out["anchor_duration_min"] = pd.to_numeric(out["duration_min"], errors="coerce")
    out["anchor_sde"] = out["sde_max"]
    small_period = pd.to_numeric(out["adp_sml_period_d"], errors="coerce")
    primary_period = pd.to_numeric(out.get("adp_period_d"), errors="coerce")
    out["aperture_period_rel_delta"] = np.abs(primary_period - small_period) / small_period
    small_depth = pd.to_numeric(out.get("adp_sml_depth"), errors="coerce")
    primary_depth = pd.to_numeric(out.get("adp_depth"), errors="coerce")
    out["aperture_depth_ratio_primary_over_small"] = primary_depth / small_depth.replace(0, np.nan)
    ratios = primary_period / small_period
    factors = np.asarray(TRANSIT_HARMONIC_FACTORS_V1, dtype=float)
    relation_error = np.min(
        np.abs(ratios.to_numpy(dtype=float)[:, None] / factors[None, :] - 1.0),
        axis=1,
    )
    out["aperture_disagreement_flag"] = ~(np.isfinite(relation_error) & (relation_error <= 0.02))
    out["n_apertures_agree"] = np.where(out["aperture_disagreement_flag"], 1, 2)
    out["apertures_agree"] = np.where(
        out["aperture_disagreement_flag"],
        "DET_FLUX_ADP_SML",
        "DET_FLUX_ADP_SML,DET_FLUX_ADP",
    )
    out["review_id"] = [
        f"s{sector:04d}-a2v1-real-{tic:016d}-r{candidate_rank}"
        for sector, tic, candidate_rank in zip(
            out["sector"], out["tic"], out["rep_peak_rank"]
        )
    ]
    _require_unique(out, "review_id", name="normalized real ADP candidates")
    out["source_kind"] = "real_candidate"
    out["is_injected_row"] = False
    out["native_input_include"] = True
    out["native_group_path"] = [native_group_path(row) for row in out.to_dict("records")]
    return out.sort_values(["tic", "rep_peak_rank"], kind="stable").reset_index(drop=True)


def mark_native_input_availability(
    candidates: pd.DataFrame,
    *,
    available_tics: Iterable[int],
) -> pd.DataFrame:
    """Mark candidates whose raw/error host exists in the native source."""

    if "tic" not in candidates:
        raise KeyError("candidate table is missing tic")
    work = candidates.copy()
    tic = pd.to_numeric(work["tic"], errors="raise").astype(np.int64)
    available = frozenset(int(value) for value in available_tics)
    raw_available = tic.isin(available)
    upstream_include = (
        _as_bool(work["native_input_include"])
        if "native_input_include" in work
        else pd.Series(True, index=work.index)
    )
    work["native_raw_source_available"] = raw_available
    work["native_input_include"] = upstream_include & raw_available
    work["native_input_status"] = np.select(
        [~upstream_include.to_numpy(), ~raw_available.to_numpy()],
        ["excluded_upstream", "missing_raw_source"],
        default="available",
    )
    return work


def _phase_delta_min(first_t0: float, second_t0: float, period_d: float) -> float:
    if not all(np.isfinite([first_t0, second_t0, period_d])) or period_d <= 0:
        return np.nan
    delta_d = ((first_t0 - second_t0 + 0.5 * period_d) % period_d) - 0.5 * period_d
    return float(delta_d * 1440.0)


def transfer_human_labels_to_a2v1_candidates(
    human_rows: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    period_tolerance: float = 0.02,
    minimum_window_overlap: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Match old reviewed ephemerides to current A2v1 ADP-small top peaks."""

    human = human_rows.copy()
    injected = (
        _as_bool(human["is_injected_row"])
        if "is_injected_row" in human
        else human.get("source_kind", pd.Series("", index=human.index))
        .fillna("")
        .astype(str)
        .str.contains("inject", case=False)
    )
    human = human.loc[~injected].copy().reset_index(drop=True)
    identity = "source_uid" if "source_uid" in human else "review_id"
    _require_unique(human, identity, name="real human labels")
    candidate_groups = {
        int(tic): group.copy()
        for tic, group in candidates.groupby("tic", sort=False)
    }
    supervision_columns = [
        column
        for column in (
            identity,
            "human_label",
            "human_label_raw",
            "human_label_source",
            "human_labeler",
            "human_notes",
            "human_notes_raw",
            "human_updated_utc",
            "human_updated_utc_raw",
            "morphology_target_v1",
            "morphology_include_v1",
            "preserve_target_v1",
            "preserve_include_v1",
            "harmonic_target_v1",
            "harmonic_include_v1",
            "broad_preserve_only",
            "model_target_policy_version",
            "ephemeris_status",
            "effective_period_factor",
            "period_factor_source",
            "period_task",
        )
        if column in human
    ]
    transfer_rows: list[dict[str, Any]] = []
    compatibility_rows: list[dict[str, Any]] = []
    factors = np.asarray(TRANSIT_HARMONIC_FACTORS_V1, dtype=float)
    for record in human.to_dict("records"):
        tic = int(float(record["tic"]))
        old_period = float(record.get("effective_period_d", np.nan))
        if not np.isfinite(old_period) or old_period <= 0:
            old_period = float(record.get("period_d", np.nan))
        old_t0 = float(record.get("t0_bjd", np.nan))
        old_duration = float(record.get("duration_min", np.nan))
        label = str(record.get("human_label", ""))
        options: list[dict[str, Any]] = []
        for candidate in candidate_groups.get(tic, pd.DataFrame()).to_dict("records"):
            new_period = float(candidate["period_d"])
            ratio = old_period / new_period if new_period > 0 else np.nan
            errors = np.abs(ratio / factors - 1.0) if np.isfinite(ratio) else np.full(len(factors), np.inf)
            factor_index = int(np.argmin(errors))
            period_error = float(errors[factor_index])
            phase_period = min(old_period, new_period)
            delta_min = _phase_delta_min(float(candidate["t0_bjd"]), old_t0, phase_period)
            new_duration = float(candidate["duration_min"])
            half_sum = 0.5 * (old_duration + new_duration)
            overlap_min = max(0.0, half_sum - abs(delta_min)) if np.isfinite(delta_min) else 0.0
            overlap_fraction = (
                min(1.0, overlap_min / min(old_duration, new_duration))
                if min(old_duration, new_duration) > 0
                else 0.0
            )
            period_ok = np.isfinite(period_error) and period_error <= float(period_tolerance)
            window_ok = overlap_fraction >= float(minimum_window_overlap)
            # Smooth variability has no discrete transit window. All other
            # labels transfer only when the reviewed event phase also agrees.
            transfer_ok = period_ok and (window_ok or label == "stellar_variability")
            options.append(
                {
                    **candidate,
                    "a2v1_period_factor": float(factors[factor_index]),
                    "a2v1_period_relative_error": period_error,
                    "a2v1_window_delta_min": delta_min,
                    "a2v1_window_overlap_fraction": overlap_fraction,
                    "a2v1_transfer_ok": bool(transfer_ok),
                }
            )
        valid = [option for option in options if option["a2v1_transfer_ok"]]
        valid.sort(
            key=lambda option: (
                option["a2v1_period_relative_error"],
                -option["a2v1_window_overlap_fraction"],
                int(option["rep_peak_rank"]),
            )
        )
        selected = valid[0] if valid else None
        compatibility_rows.append(
            {
                identity: str(record[identity]),
                "tic": tic,
                "human_label": label,
                "a2v1_transfer_ok": selected is not None,
                "a2v1_match_review_id": selected["review_id"] if selected else "",
                "a2v1_period_factor": selected["a2v1_period_factor"] if selected else np.nan,
                "a2v1_period_relative_error": selected["a2v1_period_relative_error"] if selected else np.nan,
                "a2v1_window_overlap_fraction": selected["a2v1_window_overlap_fraction"] if selected else np.nan,
                "a2v1_match_rank": int(selected["rep_peak_rank"]) if selected else -1,
                "a2v1_transfer_policy": "period_harmonic_2pct_and_phase_overlap_50pct; variable_period_only",
            }
        )
        if selected is None:
            continue
        output = dict(selected)
        output["source_review_id"] = str(record.get("review_id", record[identity]))
        output["review_id"] = "s0056-a2v1-human-" + _stable_digest(
            str(record[identity]), seed=56
        )[:20]
        for column in supervision_columns:
            output[column] = record.get(column, "")
        output["a2v1_transfer_ok"] = True
        output["source_kind"] = "real_candidate"
        output["is_injected_row"] = False
        transfer_rows.append(output)
    transferred = pd.DataFrame(transfer_rows)
    compatibility = pd.DataFrame(compatibility_rows)
    if not transferred.empty:
        _require_unique(transferred, "review_id", name="A2v1 human transfers")
        if transferred[identity].duplicated().any():
            raise RuntimeError("A2v1 human transfer produced duplicate source identities")
    return transferred.reset_index(drop=True), compatibility.reset_index(drop=True)


def build_global_tic_split_registry(
    *,
    s56_injection_roles: pd.DataFrame,
    s56_human_rows: pd.DataFrame,
    franklin_rows: pd.DataFrame,
    s57_schedule: pd.DataFrame,
) -> pd.DataFrame:
    """Create one auditable TIC registry spanning training, audit, and test."""

    injection = s56_injection_roles[["tic", "teacher_v2_partition"]].copy()
    injection["tic"] = pd.to_numeric(injection["tic"], errors="raise").astype(np.int64)
    human_tics = set(pd.to_numeric(s56_human_rows["tic"], errors="raise").astype(np.int64))
    injection_tics = set(injection["tic"])
    overlap = human_tics & injection_tics
    if overlap:
        raise ValueError(f"S56 human and fresh-injection TICs overlap: {sorted(overlap)[:5]}")
    franklin_tics = set(pd.to_numeric(franklin_rows["tic"], errors="raise").astype(np.int64))
    s57_tics = set(pd.to_numeric(s57_schedule["tic"], errors="raise").astype(np.int64))
    all_tics = sorted(human_tics | injection_tics | franklin_tics | s57_tics)
    partition_lookup = injection.set_index("tic")["teacher_v2_partition"].to_dict()
    rows: list[dict[str, Any]] = []
    for tic in all_tics:
        in_human = tic in human_tics
        in_injection = tic in injection_tics
        in_franklin = tic in franklin_tics
        in_s57 = tic in s57_tics
        if in_human:
            assigned = "s56_human"
        elif in_injection:
            assigned = f"s56_injection_{partition_lookup[tic]}"
        elif in_franklin:
            assigned = "franklin_audit_only"
        else:
            assigned = "s57_external_test"
        rows.append(
            {
                "tic": tic,
                "assigned_role": assigned,
                "in_s56_human": in_human,
                "in_s56_fresh_injection": in_injection,
                "in_franklin_audit": in_franklin,
                "in_s57_schedule": in_s57,
                "s57_external_eligible": in_s57
                and not (in_human or in_injection or in_franklin),
            }
        )
    return pd.DataFrame(rows)


def build_s57_external_role_manifest(
    schedule: pd.DataFrame,
    registry: pd.DataFrame,
    *,
    minimum_injections: int = 10_000,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Filter the immutable S57 schedule through the frozen global registry."""

    eligible = set(
        pd.to_numeric(
            registry.loc[_as_bool(registry["s57_external_eligible"]), "tic"],
            errors="raise",
        ).astype(np.int64)
    )
    work = schedule.copy()
    work["tic"] = pd.to_numeric(work["tic"], errors="raise").astype(np.int64)
    work["teacher_v2_external_include"] = work["tic"].isin(eligible)
    selected = work.loc[work["teacher_v2_external_include"]].copy().reset_index(drop=True)
    support = selected.groupby("grid_cell_id", dropna=False).size()
    expected_cells = schedule["grid_cell_id"].nunique(dropna=False)
    all_cells = len(support) == expected_cells and support.gt(0).all()
    passed = len(selected) >= int(minimum_injections) and all_cells
    summary = {
        "n_schedule": int(len(schedule)),
        "n_external": int(len(selected)),
        "n_external_unique_tics": int(selected["tic"].nunique()),
        "n_cells": int(len(support)),
        "expected_cells": int(expected_cells),
        "cell_support_min": int(support.min()) if len(support) else 0,
        "cell_support_max": int(support.max()) if len(support) else 0,
        "minimum_injections": int(minimum_injections),
        "passed": bool(passed),
        "regenerate_required": not bool(passed),
    }
    return selected, summary


def build_s56_real_workload_candidates(
    candidates: pd.DataFrame,
    registry: pd.DataFrame,
    *,
    minimum_tics: int = 1_000,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build a real workload pool disjoint from every registered model/audit host."""

    required_candidates = {"tic", "review_id"}
    missing = sorted(required_candidates - set(candidates.columns))
    if missing:
        raise KeyError(f"S56 real candidates are missing columns: {missing}")
    if "tic" not in registry:
        raise KeyError("global TIC registry is missing tic")
    work = candidates.copy()
    work["tic"] = pd.to_numeric(work["tic"], errors="raise").astype(np.int64)
    n_input_rows = int(len(work))
    n_input_tics = int(work["tic"].nunique())
    if work["review_id"].fillna("").astype(str).duplicated().any():
        raise ValueError("S56 real candidates contain duplicate review IDs")
    native_include = (
        _as_bool(work["native_input_include"])
        if "native_input_include" in work
        else pd.Series(True, index=work.index)
    )
    n_native_unavailable_tics = int(work.loc[~native_include, "tic"].nunique())
    work = work.loc[native_include].copy()
    excluded = frozenset(
        pd.to_numeric(registry["tic"], errors="raise").astype(np.int64)
    )
    selected = work.loc[~work["tic"].isin(excluded)].copy().reset_index(drop=True)
    n_tics = int(selected["tic"].nunique())
    passed = n_tics >= int(minimum_tics)
    summary = {
        "n_input_rows": n_input_rows,
        "n_input_tics": n_input_tics,
        "n_native_unavailable_tics": n_native_unavailable_tics,
        "n_registry_tics_excluded": int(work.loc[work["tic"].isin(excluded), "tic"].nunique()),
        "n_workload_rows": int(len(selected)),
        "n_workload_tics": n_tics,
        "minimum_tics": int(minimum_tics),
        "passed": bool(passed),
        "host_disjoint_from_s56_injections_human_franklin_and_s57": True,
    }
    if not passed:
        raise ValueError(
            f"untouched S56 real workload pool has only {n_tics} TICs; "
            f"minimum={minimum_tics}"
        )
    return selected, summary


def assert_teacher_v2_feature_columns(columns: Iterable[str]) -> tuple[str, ...]:
    """Reject truth, target, provenance, and adjudication columns as inputs."""

    names = tuple(str(column) for column in columns)
    leaks = set(leakage_columns(names))
    for column in names:
        normalized = column.lower()
        if any(token in normalized for token in TEACHER_V2_FORBIDDEN_FEATURE_TOKENS):
            leaks.add(column)
    if leaks:
        raise ValueError(f"Teacher-v2 feature columns contain leakage: {sorted(leaks)}")
    return names


@dataclass(frozen=True)
class WorkloadThreshold:
    threshold: float
    n_tics: int
    max_pass_tics: int
    n_pass_tics: int
    pass_fraction: float
    score_column: str


def freeze_real_tic_workload_threshold(
    scores: pd.DataFrame,
    *,
    score_column: str = "p_compact_transit",
    max_fraction: float = 0.05,
) -> WorkloadThreshold:
    """Freeze a threshold passing no more than ``max_fraction`` of real TICs."""

    if not 0.0 < float(max_fraction) < 1.0:
        raise ValueError("max_fraction must be between zero and one")
    required = {"tic", score_column}
    missing = sorted(required - set(scores.columns))
    if missing:
        raise KeyError(f"real score table is missing columns: {missing}")
    work = scores.loc[:, ["tic", score_column]].copy()
    work["tic"] = pd.to_numeric(work["tic"], errors="coerce")
    work[score_column] = pd.to_numeric(work[score_column], errors="coerce")
    if work.isna().any().any():
        raise ValueError("real score table contains invalid TIC or compact scores")
    per_tic = work.groupby("tic", sort=False)[score_column].max().to_numpy(dtype=float)
    n_tics = int(len(per_tic))
    allowed = int(np.floor(float(max_fraction) * n_tics))
    if n_tics == 0:
        raise ValueError("cannot freeze a workload threshold from zero TICs")
    if allowed == 0:
        threshold = float("inf")
    else:
        ordered = np.sort(per_tic)[::-1]
        boundary = float(ordered[allowed - 1])
        threshold = boundary
        if int(np.count_nonzero(per_tic >= boundary)) > allowed:
            threshold = float(np.nextafter(boundary, np.inf))
    passed = int(np.count_nonzero(per_tic >= threshold))
    if passed > allowed:
        raise RuntimeError("frozen workload threshold exceeds the requested TIC budget")
    return WorkloadThreshold(
        threshold=threshold,
        n_tics=n_tics,
        max_pass_tics=allowed,
        n_pass_tics=passed,
        pass_fraction=float(passed / n_tics),
        score_column=score_column,
    )


def injection_recall_at_threshold(
    scores: pd.DataFrame,
    *,
    threshold: float,
    partition: str,
) -> dict[str, Any]:
    """Measure recovery per injected host using only truth-matched candidates."""

    required = {
        "injection_id",
        "teacher_v2_partition",
        "is_injected_signal_peak",
        "p_compact_transit",
    }
    missing = sorted(required - set(scores.columns))
    if missing:
        raise KeyError(f"injection score table is missing columns: {missing}")
    selected = scores.loc[scores["teacher_v2_partition"].astype(str).eq(partition)].copy()
    if "teacher_v2_role" in selected:
        selected = selected.loc[
            ~selected["teacher_v2_role"].astype(str).eq("paired_pre_injection")
        ].copy()
    selected["_truth_match"] = _as_bool(selected["is_injected_signal_peak"])
    selected["_pass"] = pd.to_numeric(
        selected["p_compact_transit"], errors="coerce"
    ).ge(float(threshold))
    selected["_success"] = selected["_truth_match"] & selected["_pass"]
    per_injection = selected.groupby("injection_id", sort=False)["_success"].any()
    return {
        "partition": partition,
        "threshold": float(threshold),
        "n_injections": int(len(per_injection)),
        "n_recovered": int(per_injection.sum()),
        "recovery_fraction": float(per_injection.mean()) if len(per_injection) else np.nan,
    }


def injection_recall_at_fold_thresholds(
    scores: pd.DataFrame,
    *,
    thresholds: Mapping[int, float],
    partition: str,
) -> dict[str, Any]:
    """Measure OOF injection recall at fold-matched real-workload thresholds."""

    required = {
        "injection_id",
        "teacher_v2_partition",
        "is_injected_signal_peak",
        "p_compact_transit",
        "cv_fold",
    }
    missing = sorted(required - set(scores.columns))
    if missing:
        raise KeyError(f"injection OOF score table is missing columns: {missing}")
    expected = set(range(5))
    supplied = {int(key) for key in thresholds}
    if supplied != expected:
        raise ValueError(f"fold thresholds must cover 0..4; got {sorted(supplied)}")
    selected = scores.loc[
        scores["teacher_v2_partition"].astype(str).eq(partition)
    ].copy()
    if "teacher_v2_role" in selected:
        selected = selected.loc[
            ~selected["teacher_v2_role"].astype(str).eq("paired_pre_injection")
        ].copy()
    fold = pd.to_numeric(selected["cv_fold"], errors="coerce")
    if fold.isna().any() or not fold.between(0, 4).all():
        raise ValueError("development OOF injection scores require cv_fold in 0..4")
    selected["_truth_match"] = _as_bool(selected["is_injected_signal_peak"])
    selected["_threshold"] = fold.astype(int).map(
        {int(key): float(value) for key, value in thresholds.items()}
    )
    selected["_pass"] = pd.to_numeric(
        selected["p_compact_transit"], errors="coerce"
    ).ge(selected["_threshold"])
    selected["_success"] = selected["_truth_match"] & selected["_pass"]
    per_injection = selected.groupby("injection_id", sort=False)["_success"].any()
    return {
        "partition": partition,
        "fold_thresholds": {
            str(key): float(value) for key, value in sorted(thresholds.items())
        },
        "n_injections": int(len(per_injection)),
        "n_recovered": int(per_injection.sum()),
        "recovery_fraction": (
            float(per_injection.mean()) if len(per_injection) else np.nan
        ),
    }


__all__ = [
    "COMPACT_CLASSES",
    "FRANKLIN_ACTIVE_LABELS",
    "TEACHER_V2_MODEL_VERSION",
    "TEACHER_V2_ROLE_POLICY",
    "WorkloadThreshold",
    "assert_teacher_v2_feature_columns",
    "assign_s56_injection_roles",
    "build_global_tic_split_registry",
    "build_franklin_a2v1_rereview_queue",
    "build_franklin_current_a2v1_audit_candidates",
    "build_real_human_training_rows",
    "build_s56_injection_training_rows",
    "build_s57_external_role_manifest",
    "build_s56_real_workload_candidates",
    "freeze_real_tic_workload_threshold",
    "injection_recall_at_threshold",
    "injection_recall_at_fold_thresholds",
    "normalize_real_adp_candidates",
    "normalize_franklin_labels",
    "join_franklin_labels_to_queue",
    "mark_native_input_availability",
    "transfer_human_labels_to_a2v1_candidates",
]
