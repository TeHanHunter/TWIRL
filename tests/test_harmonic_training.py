from __future__ import annotations

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_training import (
    classification_metrics,
    expected_calibration_error,
    injection_truth_human_audit,
)


def test_multiclass_metrics_report_balanced_accuracy_and_confusion() -> None:
    truth = np.asarray([0, 0, 1, 1, 2, 2, 3, 3])
    predicted = np.asarray([0, 0, 1, 0, 2, 2, 3, 1])
    probability = np.full((len(truth), 4), 0.01)
    probability[np.arange(len(truth)), predicted] = 0.97

    metrics = classification_metrics(
        truth,
        probability,
        classes=("planet", "eb", "variable", "other"),
    )

    assert np.isclose(metrics["accuracy"], 0.75)
    assert np.isclose(metrics["balanced_accuracy"], 0.75)
    assert metrics["confusion_matrix"][1][0] == 1
    assert metrics["per_class"]["planet"]["recall"] == 1.0


def test_calibration_error_is_finite_and_small_for_correct_probabilities() -> None:
    truth = np.asarray([0, 1, 0, 1])
    probability = np.asarray(
        [[0.9, 0.1], [0.1, 0.9], [0.8, 0.2], [0.2, 0.8]],
        dtype=float,
    )

    calibration = expected_calibration_error(truth, probability, n_bins=5)

    assert np.isfinite(calibration["ece"])
    assert 0.0 <= calibration["ece"] <= 0.2


def test_injection_truth_audit_keeps_human_and_truth_separate(tmp_path) -> None:
    rows = pd.DataFrame(
        {
            "review_id": ["inj:1", "inj:2", "real:3"],
            "source_kind": ["injection_recovery", "injection_recovery", "real_candidate"],
            "injection_id": ["one", "two", ""],
            "tic": [1, 2, 3],
            "human_label": ["planet_like", "uncertain", "planet_like"],
            "truth_period_d": [2.0, 1.37, np.nan],
            "period_d": [1.0, 1.0, 1.0],
            "truth_radius_rearth": [2.0, 4.0, np.nan],
            "truth_model_depth": [0.1, 0.2, np.nan],
            "tmag": [18.0, 19.0, 17.0],
            "adp_sml_sde": [20.0, 3.0, 10.0],
            "boundary_period_bin_key": ["p1", "p2", ""],
            "boundary_radius_bin_key": ["r1", "r2", ""],
            "boundary_tmag_bin_key": ["m1", "m2", ""],
        }
    )

    summary = injection_truth_human_audit(rows, out_dir=tmp_path)
    audit = pd.read_csv(tmp_path / "injection_truth_human_rows.csv")

    assert summary["n_injections"] == 2
    assert summary["n_human_visible_planet"] == 1
    assert summary["n_bls_truth_period_or_harmonic_match"] == 1
    assert audit.loc[audit["injection_id"].eq("one"), "truth_to_bls_period_factor"].iloc[0] == 2.0
    assert audit.loc[audit["injection_id"].eq("two"), "human_label"].iloc[0] == "uncertain"
