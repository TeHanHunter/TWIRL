"""Development-only architecture selection and workload-threshold freezing."""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES
from twirl.vetting.harmonic_training import classification_metrics, expected_calibration_error
from twirl.vetting.teacher_v2 import (
    TEACHER_V2_MODEL_VERSION,
    freeze_real_tic_workload_threshold,
    injection_recall_at_fold_thresholds,
)
from twirl.vetting.teacher_v2_training import compact_metrics


SELECTION_FORMULA = (
    "0.60*development_injection_recall_at_frozen_real_5pct_workload + "
    "0.25*real_human_morphology_macro_f1 + 0.10*compact_candidate_AP + "
    "0.05*(1-compact_ECE)"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def _profile_development_metrics(
    predictions: pd.DataFrame,
    *,
    fold_thresholds: Mapping[int, float],
) -> dict[str, Any]:
    if predictions["review_id"].duplicated().any():
        raise ValueError("development OOF predictions contain duplicate review IDs")
    real = predictions["teacher_v2_role"].astype(str).eq("real_human_morphology")
    morphology_probability = predictions[
        [f"p_{label}" for label in MORPHOLOGY_CLASSES]
    ].to_numpy(dtype=float)
    morphology = classification_metrics(
        predictions.loc[real, "morphology_target"].to_numpy(dtype=int),
        morphology_probability[real.to_numpy()],
        classes=MORPHOLOGY_CLASSES,
    )
    candidate = ~predictions["teacher_v2_role"].astype(str).eq("paired_pre_injection")
    compact_truth = predictions.loc[candidate, "compact_target"].to_numpy(dtype=int)
    compact_score = predictions.loc[candidate, "p_compact_transit"].to_numpy(dtype=float)
    compact = compact_metrics(compact_truth, compact_score)
    active = compact_truth >= 0
    compact_probability = np.column_stack(
        [1.0 - compact_score[active], compact_score[active]]
    )
    compact_calibration = expected_calibration_error(
        compact_truth[active], compact_probability
    )
    injection = injection_recall_at_fold_thresholds(
        predictions,
        thresholds=fold_thresholds,
        partition="development",
    )
    macro_f1 = float(morphology["macro_f1"])
    ap = float(compact["average_precision"])
    ece = float(compact_calibration["ece"])
    values = [macro_f1, ap, ece, float(injection["recovery_fraction"])]
    if not all(np.isfinite(values)):
        selection_score = -np.inf
    else:
        selection_score = (
            0.60 * float(injection["recovery_fraction"])
            + 0.25 * macro_f1
            + 0.10 * ap
            + 0.05 * (1.0 - min(max(ece, 0.0), 1.0))
        )
    return {
        "selection_score": float(selection_score),
        "injection_recall_at_workload": injection,
        "real_human_morphology": morphology,
        "compact_candidates_only": compact,
        "compact_calibration": compact_calibration,
    }


def freeze_teacher_v2_selection(
    *,
    profile_prediction_paths: Mapping[str, Sequence[Path]],
    profile_real_score_paths: Mapping[str, Path],
    profile_checkpoint_paths: Mapping[str, Sequence[Path]],
    max_real_tic_fraction: float = 0.05,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Select once using development OOF labels and unlabeled real workload."""

    profiles = sorted(profile_prediction_paths)
    if set(profiles) != set(profile_real_score_paths) or set(profiles) != set(
        profile_checkpoint_paths
    ):
        raise ValueError("profile predictions, real scores, and checkpoints must have identical keys")
    rows: list[dict[str, Any]] = []
    details: dict[str, Any] = {}
    for profile in profiles:
        prediction_parts = [_read(Path(path)) for path in profile_prediction_paths[profile]]
        predictions = pd.concat(prediction_parts, ignore_index=True)
        real_scores = _read(Path(profile_real_score_paths[profile]))
        threshold = freeze_real_tic_workload_threshold(
            real_scores, max_fraction=max_real_tic_fraction
        )
        fold_workloads = {
            fold: freeze_real_tic_workload_threshold(
                real_scores,
                score_column=f"member_{fold}_p_compact_transit",
                max_fraction=max_real_tic_fraction,
            )
            for fold in range(5)
        }
        metrics = _profile_development_metrics(
            predictions,
            fold_thresholds={
                fold: workload.threshold for fold, workload in fold_workloads.items()
            },
        )
        morphology = metrics["real_human_morphology"]
        row = {
            "profile": profile,
            "selection_score": metrics["selection_score"],
            "frozen_threshold": threshold.threshold,
            "real_tic_pass_fraction": threshold.pass_fraction,
            "development_injection_recall": metrics[
                "injection_recall_at_workload"
            ]["recovery_fraction"],
            "real_morphology_macro_f1": morphology["macro_f1"],
            "compact_candidate_average_precision": metrics[
                "compact_candidates_only"
            ]["average_precision"],
            "compact_ece": metrics["compact_calibration"]["ece"],
        }
        rows.append(row)
        details[profile] = {
            "workload_threshold": asdict(threshold),
            "oof_fold_workload_thresholds": {
                str(fold): asdict(workload)
                for fold, workload in fold_workloads.items()
            },
            "development_metrics": metrics,
            "development_prediction_paths": [
                str(Path(path)) for path in profile_prediction_paths[profile]
            ],
            "real_score_path": str(Path(profile_real_score_paths[profile])),
        }
    ranking = pd.DataFrame(rows).sort_values(
        ["selection_score", "development_injection_recall"],
        ascending=False,
        kind="stable",
    )
    if ranking.empty or not np.isfinite(float(ranking.iloc[0]["selection_score"])):
        raise RuntimeError("no Teacher-v2 profile produced a finite development selection score")
    selected = str(ranking.iloc[0]["profile"])
    checkpoints = [Path(path) for path in profile_checkpoint_paths[selected]]
    if len(checkpoints) != 5 or not all(path.exists() for path in checkpoints):
        raise FileNotFoundError("selected Teacher-v2 ensemble is missing fold checkpoints")
    frozen = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_version": TEACHER_V2_MODEL_VERSION,
        "selection_formula": SELECTION_FORMULA,
        "selected_profile": selected,
        "frozen_real_tic_workload_fraction": float(max_real_tic_fraction),
        "frozen_compact_threshold": float(ranking.iloc[0]["frozen_threshold"]),
        "selected_checkpoints": [str(path) for path in checkpoints],
        "selected_checkpoint_sha256": {str(path): _sha256(path) for path in checkpoints},
        "profile_details": details,
        "profile_ranking": ranking.to_dict("records"),
        "architecture_frozen": True,
        "threshold_frozen": True,
        "s56_holdout_opened": False,
        "s57_opened": False,
    }
    return ranking.reset_index(drop=True), frozen


__all__ = ["SELECTION_FORMULA", "freeze_teacher_v2_selection"]
