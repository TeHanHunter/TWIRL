"""Frozen seven-sector corpus contract for the operational Teacher v3 run."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.vetting.adjudication_audit import (
    HARMONIC_CNN_TARGET_POLICY,
    MORPHOLOGY_TARGET_MAP_V1,
    PRESERVE_LABELS_V1,
    REJECT_LABELS_V1,
    TRANSIT_HARMONIC_LABELS_V1,
)
from twirl.vetting.harmonic_cnn import MODEL_VERSION
from twirl.vetting.label_io import candidate_key
from twirl.vetting.multisector_signal_review import FINAL_LABELS


TEACHER_V3_RUN_NAME = "teacher_v3"
TEACHER_V3_CORPUS_VERSION = "teacher_v3_s56_s62_morphology_corpus_v1"
TEACHER_V3_SPLIT_SEED = 560062
TEACHER_V3_SOURCE_TAGS: tuple[str, ...] = (
    "s56_adjudicated_training_v1",
    "s57_s59_tehan_accepted_v1",
    "s60_s62_tehan_accepted_v1",
)

_SIGNAL_UPDATE_COLUMNS: Mapping[str, str] = {
    "human_label": "final_human_label",
    "human_label_source": "human_label_source",
    "human_labeler": "final_labeler",
    "human_notes": "final_notes",
    "human_updated_utc": "final_updated_utc",
    "morphology_target_v1": "morphology_target_v1",
    "morphology_include_v1": "morphology_include_v1",
    "preserve_target_v1": "preserve_target_v1",
    "preserve_include_v1": "preserve_include_v1",
    "harmonic_target_v1": "harmonic_target_v1",
    "harmonic_include_v1": "harmonic_include_v1",
    "harmonic_supervision_verified": "harmonic_supervision_verified",
    "broad_preserve_only": "broad_preserve_only",
    "model_target_policy_version": "model_target_policy_version",
}


def file_sha256(path: Path) -> str:
    """Hash one immutable input or output."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value).strip()


def _truth(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "1.0", "true", "t", "yes", "y"})
    )


def _require(
    frame: pd.DataFrame,
    columns: set[str],
    *,
    name: str,
) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise KeyError(f"{name} is missing columns: {missing}")


def _normalize_source(
    frame: pd.DataFrame,
    *,
    source_tag: str,
) -> pd.DataFrame:
    required = {
        "review_id",
        "source_uid",
        "tic",
        "sector",
        "period_d",
        "t0_bjd",
        "duration_min",
        "human_label",
        "morphology_target_v1",
        "morphology_include_v1",
        "preserve_target_v1",
        "preserve_include_v1",
        "harmonic_target_v1",
        "harmonic_include_v1",
        "broad_preserve_only",
        "model_target_policy_version",
        "is_injected_row",
    }
    _require(frame, required, name=source_tag)
    work = frame.copy()
    sector_text = work["sector"].fillna("").astype(str).str.strip()
    allowed_sectors = {
        TEACHER_V3_SOURCE_TAGS[0]: {56},
        TEACHER_V3_SOURCE_TAGS[1]: {57, 58, 59},
        TEACHER_V3_SOURCE_TAGS[2]: {60, 61, 62},
    }[source_tag]
    if allowed_sectors == {56}:
        sector_text = sector_text.mask(sector_text.eq(""), "56")
    if sector_text.eq("").any():
        raise ValueError(f"{source_tag} contains blank sector values")
    sectors = pd.to_numeric(sector_text, errors="raise").astype(int)
    if not set(sectors).issubset(allowed_sectors):
        raise ValueError(
            f"{source_tag} contains sectors outside {sorted(allowed_sectors)}"
        )
    work["sector"] = sectors
    work["teacher_v3_source_tag"] = source_tag
    return work


def _validate_targets(frame: pd.DataFrame) -> None:
    labels = frame["human_label"].fillna("").astype(str)
    invalid_labels = sorted(set(labels) - FINAL_LABELS)
    if invalid_labels:
        raise ValueError(f"Teacher v3 corpus has invalid labels: {invalid_labels}")
    policy = (
        frame["model_target_policy_version"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    if not policy.eq(HARMONIC_CNN_TARGET_POLICY).all():
        observed = sorted(policy[policy.ne(HARMONIC_CNN_TARGET_POLICY)].unique())
        raise ValueError(
            "Teacher v3 corpus has unexpected target policies: "
            f"{observed}"
        )

    morphology = frame["morphology_target_v1"].fillna("").astype(str)
    expected_morphology = labels.map(MORPHOLOGY_TARGET_MAP_V1).fillna("")
    morphology_include = _truth(frame["morphology_include_v1"])
    if morphology.ne(expected_morphology).any() or morphology_include.ne(
        expected_morphology.ne("")
    ).any():
        raise ValueError("Teacher v3 morphology targets disagree with labels")

    preserve = frame["preserve_target_v1"].fillna("").astype(str)
    expected_preserve = pd.Series("", index=frame.index, dtype=str)
    expected_preserve.loc[labels.isin(PRESERVE_LABELS_V1)] = "preserve"
    expected_preserve.loc[labels.isin(REJECT_LABELS_V1)] = "reject"
    preserve_include = _truth(frame["preserve_include_v1"])
    if preserve.ne(expected_preserve).any() or preserve_include.ne(
        expected_preserve.ne("")
    ).any():
        raise ValueError("Teacher v3 preserve targets disagree with labels")

    broad = _truth(frame["broad_preserve_only"])
    if broad.ne(labels.eq("wide_transit_like")).any():
        raise ValueError("Teacher v3 broad-preserve flags disagree with labels")
    harmonic = _truth(frame["harmonic_include_v1"])
    if (
        harmonic
        & ~labels.isin(TRANSIT_HARMONIC_LABELS_V1)
    ).any():
        raise ValueError(
            "Teacher v3 harmonic supervision survived a non-transit relabel"
        )


def build_teacher_v3_corpus(
    *,
    s56_training: pd.DataFrame,
    s57_s59_labels: pd.DataFrame,
    s60_s62_labels: pd.DataFrame,
    signal_freeze: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Combine the accepted sources and apply the exact signal re-review."""

    sources = (
        _normalize_source(
            s56_training,
            source_tag=TEACHER_V3_SOURCE_TAGS[0],
        ),
        _normalize_source(
            s57_s59_labels,
            source_tag=TEACHER_V3_SOURCE_TAGS[1],
        ),
        _normalize_source(
            s60_s62_labels,
            source_tag=TEACHER_V3_SOURCE_TAGS[2],
        ),
    )
    base = pd.concat(sources, ignore_index=True, sort=False)
    base["source_uid"] = (
        base["source_uid"].fillna("").astype(str).str.strip()
    )
    if base["source_uid"].eq("").any() or base["source_uid"].duplicated().any():
        raise ValueError(
            "Teacher v3 source rows have blank or duplicate source_uid values"
        )

    _require(
        signal_freeze,
        {
            "row_id",
            "review_id",
            "selected_source_uid",
            "sector",
            "tic",
            "observation_candidate_key",
            "morphology_review_status",
            "rereview_period_factor",
            "rereview_period_status",
            "rereview_period_is_audit_only",
            *_SIGNAL_UPDATE_COLUMNS.values(),
        },
        name="accepted signal re-review",
    )
    decisions = signal_freeze.copy()
    decisions["selected_source_uid"] = (
        decisions["selected_source_uid"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    if (
        decisions["selected_source_uid"].eq("").any()
        or decisions["selected_source_uid"].duplicated().any()
    ):
        raise ValueError(
            "accepted signal review has blank or duplicate selected_source_uid"
        )
    base_uids = set(base["source_uid"])
    decisions["matched_existing_source"] = decisions[
        "selected_source_uid"
    ].isin(base_uids)
    unmatched = decisions.loc[~decisions["matched_existing_source"]]
    if not unmatched["sector"].astype(int).eq(56).all():
        raise ValueError(
            "only additional S56 observations may be absent from the accepted "
            "source tables"
        )

    for column in (
        "signal_rereview_applied",
        "signal_rereview_row_id",
        "signal_rereview_status",
        "signal_rereview_period_factor",
        "signal_rereview_period_status",
        "signal_rereview_period_is_audit_only",
    ):
        base[column] = False if column.endswith(("applied", "audit_only")) else ""

    base_by_uid = pd.Series(base.index, index=base["source_uid"])
    for decision in decisions.loc[
        decisions["matched_existing_source"]
    ].to_dict("records"):
        index = int(base_by_uid[_text(decision["selected_source_uid"])])
        for destination, source in _SIGNAL_UPDATE_COLUMNS.items():
            base.at[index, destination] = decision.get(source, "")
        base.at[index, "signal_rereview_applied"] = True
        base.at[index, "signal_rereview_row_id"] = decision["row_id"]
        base.at[index, "signal_rereview_status"] = decision[
            "morphology_review_status"
        ]
        base.at[index, "signal_rereview_period_factor"] = decision[
            "rereview_period_factor"
        ]
        base.at[index, "signal_rereview_period_status"] = decision[
            "rereview_period_status"
        ]
        base.at[index, "signal_rereview_period_is_audit_only"] = True

    additional = unmatched.copy()
    additional["source_uid"] = additional["selected_source_uid"]
    additional["source_kind"] = "real_candidate"
    additional["source_bucket"] = "review_candidate"
    additional["is_injected_row"] = False
    additional["teacher_v3_source_tag"] = "s56_signal_rereview_addition_v1"
    additional["signal_rereview_applied"] = True
    additional["signal_rereview_row_id"] = additional["row_id"]
    additional["signal_rereview_status"] = additional[
        "morphology_review_status"
    ]
    additional["signal_rereview_period_factor"] = additional[
        "rereview_period_factor"
    ]
    additional["signal_rereview_period_status"] = additional[
        "rereview_period_status"
    ]
    additional["signal_rereview_period_is_audit_only"] = True
    for destination, source in _SIGNAL_UPDATE_COLUMNS.items():
        additional[destination] = additional[source]
    additional["human_labeler"] = additional["final_labeler"]
    additional["human_notes"] = additional["final_notes"]
    additional["human_updated_utc"] = additional["final_updated_utc"]
    additional["candidate_key"] = additional["observation_candidate_key"]
    additional["pipeline_candidate_key"] = additional[
        "observation_candidate_key"
    ]

    corpus = pd.concat([base, additional], ignore_index=True, sort=False)
    corpus["teacher_v3_observation_id"] = corpus["source_uid"]
    if (
        corpus["teacher_v3_observation_id"].eq("").any()
        or corpus["teacher_v3_observation_id"].duplicated().any()
    ):
        raise ValueError("Teacher v3 observation identities are not unique")
    sectors = pd.to_numeric(corpus["sector"], errors="raise").astype(int)
    if set(sectors) != set(range(56, 63)):
        raise ValueError(
            f"Teacher v3 corpus sectors must be S56-S62; got {sorted(set(sectors))}"
        )
    corpus["sector"] = sectors
    corpus["tic"] = pd.to_numeric(corpus["tic"], errors="raise").astype(np.int64)
    computed_key = corpus.apply(candidate_key, axis=1)
    existing_key = (
        corpus.get(
            "observation_candidate_key",
            pd.Series("", index=corpus.index),
        )
        .fillna("")
        .astype(str)
        .str.strip()
    )
    corpus["observation_candidate_key"] = existing_key.where(
        existing_key.ne(""),
        computed_key,
    )
    morphology = (
        corpus["morphology_target_v1"].fillna("").astype(str).str.strip()
    )
    broad = _truth(corpus["broad_preserve_only"])
    corpus["split_stratum"] = morphology.where(
        morphology.ne(""),
        np.where(broad, "broad_preserve", "auxiliary"),
    )
    corpus["teacher_v3_training_include"] = (
        _truth(corpus["morphology_include_v1"])
        | _truth(corpus["preserve_include_v1"])
        | _truth(corpus["harmonic_include_v1"])
    )
    corpus["teacher_run_name"] = TEACHER_V3_RUN_NAME
    corpus["teacher_v3_corpus_version"] = TEACHER_V3_CORPUS_VERSION
    corpus["teacher_architecture_version"] = MODEL_VERSION
    _validate_targets(corpus)

    applied = _truth(corpus["signal_rereview_applied"])
    if int(applied.sum()) != len(decisions):
        raise ValueError(
            "accepted signal decisions were not applied exactly once: "
            f"{int(applied.sum())} != {len(decisions)}"
        )
    real = ~_truth(corpus["is_injected_row"])
    active = _truth(corpus["teacher_v3_training_include"])
    summary: dict[str, Any] = {
        "corpus_contract_version": TEACHER_V3_CORPUS_VERSION,
        "teacher_run_name": TEACHER_V3_RUN_NAME,
        "teacher_architecture_version": MODEL_VERSION,
        "split_seed": TEACHER_V3_SPLIT_SEED,
        "split_group": "tic",
        "split_identity_columns": [
            "sector",
            "tic",
            "teacher_v3_observation_id",
        ],
        "n_rows": int(len(corpus)),
        "n_active_training_rows": int(active.sum()),
        "n_real_rows": int(real.sum()),
        "n_injected_rows": int((~real).sum()),
        "n_unique_tics": int(corpus["tic"].nunique()),
        "n_real_unique_tics": int(corpus.loc[real, "tic"].nunique()),
        "n_signal_review_rows": int(applied.sum()),
        "n_signal_review_matched_existing": int(
            decisions["matched_existing_source"].sum()
        ),
        "n_signal_review_additional_s56": int(
            (~decisions["matched_existing_source"]).sum()
        ),
        "sector_counts": {
            str(key): int(value)
            for key, value in corpus["sector"].value_counts().sort_index().items()
        },
        "human_label_counts": {
            str(key): int(value)
            for key, value in corpus["human_label"]
            .value_counts()
            .sort_index()
            .items()
        },
        "real_human_label_counts": {
            str(key): int(value)
            for key, value in corpus.loc[real, "human_label"]
            .value_counts()
            .sort_index()
            .items()
        },
        "real_unique_tics_by_label": {
            str(key): int(value)
            for key, value in corpus.loc[real]
            .groupby("human_label")["tic"]
            .nunique()
            .sort_index()
            .items()
        },
        "split_stratum_counts": {
            str(key): int(value)
            for key, value in corpus["split_stratum"]
            .value_counts()
            .sort_index()
            .items()
        },
    }
    return corpus, summary


def write_teacher_v3_corpus(
    *,
    s56_training_path: Path,
    s57_s59_labels_path: Path,
    s60_s62_labels_path: Path,
    signal_freeze_path: Path,
    output_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    """Build and publish the immutable corpus plus a hash-bound summary."""

    inputs = {
        "s56_training": Path(s56_training_path),
        "s57_s59_labels": Path(s57_s59_labels_path),
        "s60_s62_labels": Path(s60_s62_labels_path),
        "signal_freeze": Path(signal_freeze_path),
    }
    input_hashes = {
        name: file_sha256(path) for name, path in inputs.items()
    }
    frames = {
        name: pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)
        for name, path in inputs.items()
    }
    corpus, summary = build_teacher_v3_corpus(
        s56_training=frames["s56_training"],
        s57_s59_labels=frames["s57_s59_labels"],
        s60_s62_labels=frames["s60_s62_labels"],
        signal_freeze=frames["signal_freeze"],
    )
    if {
        name: file_sha256(path) for name, path in inputs.items()
    } != input_hashes:
        raise RuntimeError("Teacher v3 source changed while corpus was built")

    output_path = Path(output_path)
    summary_path = Path(summary_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    corpus_payload = corpus.to_csv(
        index=False,
        lineterminator="\n",
        float_format="%.15g",
    ).encode("utf-8")
    corpus_hash = hashlib.sha256(corpus_payload).hexdigest()
    summary = {
        **summary,
        "inputs": {
            name: {
                "path": str(path.resolve()),
                "sha256": input_hashes[name],
            }
            for name, path in inputs.items()
        },
        "output_corpus": {
            "path": str(output_path.resolve()),
            "sha256": corpus_hash,
        },
    }
    summary_payload = (
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    for path, payload in (
        (output_path, corpus_payload),
        (summary_path, summary_payload),
    ):
        if path.exists() and path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace immutable Teacher v3 output: {path}"
            )
    if not output_path.exists():
        output_path.write_bytes(corpus_payload)
    if not summary_path.exists():
        summary_path.write_bytes(summary_payload)
    return summary


__all__ = [
    "TEACHER_V3_CORPUS_VERSION",
    "TEACHER_V3_RUN_NAME",
    "TEACHER_V3_SPLIT_SEED",
    "build_teacher_v3_corpus",
    "file_sha256",
    "write_teacher_v3_corpus",
]
