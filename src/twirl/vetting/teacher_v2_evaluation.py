"""Locked human-holdout evaluation for the S56 Teacher-v2 cycle."""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES
from twirl.vetting.harmonic_training import classification_metrics


def _probability_columns() -> list[str]:
    return [f"p_{label}" for label in MORPHOLOGY_CLASSES]


def _validate_scores(scores: pd.DataFrame, *, name: str) -> pd.DataFrame:
    required = {"review_id", *_probability_columns()}
    missing = sorted(required - set(scores.columns))
    if missing:
        raise KeyError(f"{name} scores are missing columns: {missing}")
    out = scores[["review_id", *_probability_columns()]].copy()
    out["review_id"] = out["review_id"].fillna("").astype(str)
    if out["review_id"].eq("").any() or out["review_id"].duplicated().any():
        raise ValueError(f"{name} score review IDs must be nonempty and unique")
    probability = out[_probability_columns()].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(probability.to_numpy(dtype=float)).all():
        raise ValueError(f"{name} scores contain non-finite morphology probabilities")
    if len(probability) and not np.allclose(
        probability.to_numpy(dtype=float).sum(axis=1), 1.0, atol=1.0e-4
    ):
        raise ValueError(f"{name} morphology probabilities do not sum to one")
    out[_probability_columns()] = probability
    return out


def evaluate_locked_human_holdout(
    human_rows: pd.DataFrame,
    teacher_v2_scores: pd.DataFrame,
    teacher_v1_scores: pd.DataFrame,
    *,
    maximum_macro_f1_decline: float = 0.03,
    minimum_nonplanet_recall: float = 0.60,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Compare v1 and v2 on exactly the same locked real human examples."""

    required = {
        "review_id",
        "fixed_split",
        "teacher_v2_role",
        "morphology_target_index",
        "human_label",
    }
    missing = sorted(required - set(human_rows.columns))
    if missing:
        raise KeyError(f"human role table is missing columns: {missing}")
    truth = human_rows.loc[
        human_rows["fixed_split"].astype(str).eq("test")
        & human_rows["teacher_v2_role"].astype(str).eq("real_human_morphology")
        & pd.to_numeric(
            human_rows["morphology_target_index"], errors="coerce"
        ).ge(0)
    ].copy()
    truth["review_id"] = truth["review_id"].fillna("").astype(str)
    if truth.empty:
        raise ValueError("locked human holdout contains no active morphology rows")
    if truth["review_id"].eq("").any() or truth["review_id"].duplicated().any():
        raise ValueError("locked human holdout review IDs must be nonempty and unique")

    v2 = _validate_scores(teacher_v2_scores, name="Teacher-v2").rename(
        columns={column: f"v2_{column}" for column in _probability_columns()}
    )
    v1 = _validate_scores(teacher_v1_scores, name="Teacher-v1").rename(
        columns={column: f"v1_{column}" for column in _probability_columns()}
    )
    expected = set(truth["review_id"])
    for name, scores in (("Teacher-v2", v2), ("Teacher-v1", v1)):
        observed = set(scores["review_id"])
        missing_ids = expected - observed
        if missing_ids:
            raise ValueError(
                f"{name} is missing {len(missing_ids)} locked holdout rows; "
                f"first={sorted(missing_ids)[:5]}"
            )
    joined = truth.merge(v2, on="review_id", how="left", validate="one_to_one")
    joined = joined.merge(v1, on="review_id", how="left", validate="one_to_one")
    target = pd.to_numeric(
        joined["morphology_target_index"], errors="raise"
    ).to_numpy(dtype=int)
    metrics: dict[str, dict[str, Any]] = {}
    for version in ("v1", "v2"):
        probability = joined[
            [f"{version}_p_{label}" for label in MORPHOLOGY_CLASSES]
        ].to_numpy(dtype=float)
        metrics[version] = classification_metrics(
            target, probability, classes=MORPHOLOGY_CLASSES
        )
        joined[f"{version}_predicted_morphology"] = np.asarray(
            MORPHOLOGY_CLASSES, dtype=object
        )[probability.argmax(axis=1)]

    macro_delta = float(metrics["v2"]["macro_f1"] - metrics["v1"]["macro_f1"])
    required_recall: dict[str, float] = {}
    for label in ("eclipse_contact", "smooth_variable", "other"):
        value = float(metrics["v2"]["per_class"][label]["recall"])
        required_recall[label] = value
    finite_recall = all(np.isfinite(value) for value in required_recall.values())
    acceptance = {
        "macro_f1_non_regression": bool(
            np.isfinite(macro_delta) and macro_delta >= -float(maximum_macro_f1_decline)
        ),
        "eb_variable_other_recall": bool(
            finite_recall
            and all(
                value >= float(minimum_nonplanet_recall)
                for value in required_recall.values()
            )
        ),
    }
    acceptance["passed"] = bool(all(acceptance.values()))
    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_same_row_locked_human_examples": int(len(joined)),
        "teacher_v1": metrics["v1"],
        "teacher_v2": metrics["v2"],
        "teacher_v2_minus_v1_macro_f1": macro_delta,
        "required_teacher_v2_recall": required_recall,
        "maximum_allowed_macro_f1_decline": float(maximum_macro_f1_decline),
        "minimum_required_eb_variable_other_recall": float(
            minimum_nonplanet_recall
        ),
        "real_planet_performance_is_descriptive_only": True,
        "acceptance": acceptance,
    }
    return joined.reset_index(drop=True), summary


__all__ = ["evaluate_locked_human_holdout"]
