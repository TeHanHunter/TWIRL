from __future__ import annotations

import numpy as np
import pandas as pd

from twirl.vetting.self_training import (
    SelfTrainingConfig,
    load_model,
    save_model,
    train_teacher_student,
)


def _synthetic_candidates(n: int = 120) -> pd.DataFrame:
    rng = np.random.default_rng(123)
    rows = []
    for i in range(n):
        is_good = i < n // 2
        sde = rng.normal(42.0 if is_good else 9.0, 3.0)
        period = rng.normal(1.3 if is_good else 0.16, 0.04)
        duration = rng.normal(10.0 if is_good else 29.0, 1.5)
        rows.append(
            {
                "tic": 100000 + i,
                "sector": 56,
                "tmag": rng.uniform(14.0, 19.5),
                "period_d": period,
                "duration_min": duration,
                "depth": rng.normal(0.25 if is_good else 0.03, 0.01),
                "sde_max": sde,
                "n_apertures_agree": 2 if is_good else 1,
                "vet_class": "planet_candidate" if is_good else "alias_artifact",
                "class_rank": i + 1,
            }
        )
    return pd.DataFrame(rows)


def test_teacher_student_selects_pseudo_labels_and_review_queue(tmp_path) -> None:
    candidates = _synthetic_candidates()
    labels = pd.DataFrame(
        {
            "tic": [100000, 100001, 100002, 100060, 100061, 100062],
            "label": [
                "transit_like",
                "transit_like",
                "transit_like",
                "false_positive",
                "false_positive",
                "false_positive",
            ],
            "label_source": ["human"] * 6,
        }
    )

    cfg = SelfTrainingConfig(
        pseudo_min_confidence=0.80,
        pseudo_min_margin=0.30,
        review_queue_size=10,
    )
    result = train_teacher_student(candidates, labels, cfg)

    assert result.summary["n_human_labeled"] == 6
    assert result.summary["n_pseudo"] > 0
    assert set(result.student.classes) == {"false_positive", "transit_like"}
    assert "student_label" in result.scored.columns
    assert len(result.review_queue) <= 10
    assert not set(labels["tic"]).intersection(set(result.pseudo_labels["tic"]))

    model_path = tmp_path / "student_model.npz"
    save_model(result.student, model_path)
    loaded = load_model(model_path)

    p0 = result.student.predict_proba(candidates.head(5))
    p1 = loaded.predict_proba(candidates.head(5))
    assert np.allclose(p0, p1)


def test_teacher_requires_two_classes() -> None:
    candidates = _synthetic_candidates()
    labels = pd.DataFrame(
        {
            "tic": [100000, 100001, 100002],
            "label": ["transit_like", "transit_like", "transit_like"],
        }
    )

    try:
        train_teacher_student(candidates, labels, SelfTrainingConfig())
    except ValueError as exc:
        assert "at least two classes" in str(exc)
    else:
        raise AssertionError("single-class labels should not train")
