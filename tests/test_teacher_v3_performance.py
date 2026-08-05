from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from twirl.plotting.teacher_v3_performance import (
    build_teacher_v3_performance_figure,
)
from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES
from twirl.vetting.harmonic_training import (
    _calibration_by_source,
    _subset_metrics,
)
import twirl.vetting.teacher_v3_training as teacher_v3_training
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_BASELINE_PROFILE,
    TEACHER_V3_PRIMARY_LABEL_POLICY,
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
    TEACHER_V3_UNCERTAIN_MASKED_POLICY,
    write_teacher_v3_fixed_test_bootstrap,
)


def _probability(truth: np.ndarray, *, strength: float) -> np.ndarray:
    probability = np.full(
        (len(truth), len(MORPHOLOGY_CLASSES)),
        (1.0 - strength) / (len(MORPHOLOGY_CLASSES) - 1),
        dtype=float,
    )
    probability[np.arange(len(truth)), truth] = strength
    probability[1, :] = [0.10, 0.65, 0.15, 0.10]
    return probability


def _write_profile(
    *,
    profile_dir: Path,
    rows: pd.DataFrame,
    probability: np.ndarray,
    profile: str,
    label_policy: str,
    common_support: bool,
) -> dict[str, object]:
    profile_dir.mkdir(parents=True, exist_ok=True)
    truth = rows["morphology_target_index"].to_numpy(dtype=int)
    predictions = rows.copy()
    predictions["morphology_prediction_index"] = probability.argmax(axis=1)
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        predictions[f"p_{label}"] = probability[:, index]
    predictions_path = profile_dir / "fixed_test_predictions.csv"
    predictions.to_csv(predictions_path, index=False)
    metrics = {
        "morphology_by_source": _subset_metrics(
            rows,
            truth,
            probability,
        ),
        "calibration": _calibration_by_source(
            rows,
            truth,
            probability,
        ),
        "label_policy": label_policy,
    }
    (profile_dir / "fixed_test_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    calibration_path = profile_dir / "pooled_oof_calibration.json"
    calibration_path.write_text(
        json.dumps(
            {
                "run_id": TEACHER_V3_RUN_ID,
                "profile": profile,
                "scope": (
                    "concatenated_five_fold_development_oof_logits"
                ),
            },
            sort_keys=True,
        )
        + "\n"
    )
    if common_support:
        selected = ~rows["human_label"].eq("uncertain").to_numpy()
        write_teacher_v3_fixed_test_bootstrap(
            rows=rows.loc[selected].reset_index(drop=True),
            truth=truth[selected],
            probability=probability[selected],
            profile=profile,
            label_policy=label_policy,
            out_dir=profile_dir,
            predictions_path=predictions_path,
            draws=50,
            seed=123,
            filename_prefix="fixed_test_non_uncertain",
            evaluation_scope="fixed_test_real_non_uncertain",
        )
    else:
        write_teacher_v3_fixed_test_bootstrap(
            rows=rows,
            truth=truth,
            probability=probability,
            profile=profile,
            label_policy=label_policy,
            out_dir=profile_dir,
            predictions_path=predictions_path,
            draws=50,
            seed=123,
            evaluation_scope="fixed_test_real_non_uncertain",
        )
    return {
        "profile": profile,
        "calibration_path": str(calibration_path),
        "calibration_sha256": teacher_v3_training._file_sha256(
            calibration_path
        ),
    }


def _synthetic_run(tmp_path: Path) -> Path:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    labels = [
        "planet_like",
        "planet_like",
        "eclipsing_binary_or_pceb",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "stellar_variability",
        "instrumental_or_systematic",
        "uncertain",
        "planet_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "instrumental_or_systematic",
    ]
    truth = np.asarray([0, 0, 1, 1, 2, 2, 3, 3, 0, 1, 2, 3])
    rows = pd.DataFrame(
        {
            "review_id": [f"row_{index}" for index in range(len(labels))],
            "sector": [56] * len(labels),
            "tic": [1, 1, 2, 3, 4, 5, 6, 7, 101, 102, 103, 104],
            "human_label": labels,
            "is_injected_row": [False] * 8 + [True] * 4,
            "fixed_split": ["test"] * len(labels),
            "morphology_target_index": truth,
        }
    )
    baseline_calibration = _write_profile(
        profile_dir=run_dir / TEACHER_V3_BASELINE_PROFILE,
        rows=rows,
        probability=_probability(truth, strength=0.55),
        profile=TEACHER_V3_BASELINE_PROFILE,
        label_policy=TEACHER_V3_PRIMARY_LABEL_POLICY,
        common_support=True,
    )
    primary_calibration = _write_profile(
        profile_dir=run_dir / TEACHER_V3_PRIMARY_PROFILE,
        rows=rows,
        probability=_probability(truth, strength=0.82),
        profile=TEACHER_V3_PRIMARY_PROFILE,
        label_policy=TEACHER_V3_PRIMARY_LABEL_POLICY,
        common_support=True,
    )
    sensitivity_rows = rows.loc[
        ~rows["human_label"].eq("uncertain")
    ].reset_index(drop=True)
    sensitivity_truth = sensitivity_rows[
        "morphology_target_index"
    ].to_numpy(dtype=int)
    sensitivity_calibration = _write_profile(
        profile_dir=(
            run_dir
            / "sensitivities"
            / "uncertain_masked"
            / TEACHER_V3_PRIMARY_PROFILE
        ),
        rows=sensitivity_rows,
        probability=_probability(sensitivity_truth, strength=0.78),
        profile=TEACHER_V3_PRIMARY_PROFILE,
        label_policy=TEACHER_V3_UNCERTAIN_MASKED_POLICY,
        common_support=False,
    )
    freeze_path = run_dir / "pretest_model_freeze.json"
    freeze_path.write_text(
        json.dumps(
            {
                "schema_version": "twirl_teacher_v3_model_freeze_v1",
                "run_id": TEACHER_V3_RUN_ID,
                "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
                "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
                "profile_selection_policy": "fixed_before_test",
                "test_rows_used_for_selection_or_calibration": False,
                "bootstrap": {
                    "unit": "tic",
                    "draws": 50,
                    "seed": 123,
                    "confidence_level": 0.95,
                },
            },
            sort_keys=True,
        )
        + "\n"
    )
    selected_manifest = run_dir / "selected_checkpoint_manifest.json"
    selected_manifest.write_text('{"selected": "primary"}\n')
    metadata_manifest = (
        run_dir / "metadata_baseline_checkpoint_manifest.json"
    )
    metadata_manifest.write_text('{"selected": "metadata"}\n')
    sensitivity_manifest = (
        run_dir
        / "sensitivities"
        / "uncertain_masked"
        / "uncertain_masked_checkpoint_manifest.json"
    )
    sensitivity_manifest.write_text('{"selected": "uncertain_masked"}\n')
    native_manifest = run_dir / "native_input_manifest.json"
    native_manifest.write_text('{"native": "synthetic"}\n')
    freeze_payload = json.loads(freeze_path.read_text())
    freeze_payload.update(
        {
            "primary_checkpoint_manifest": str(selected_manifest),
            "primary_checkpoint_manifest_sha256": (
                teacher_v3_training._file_sha256(selected_manifest)
            ),
            "metadata_baseline": {
                "checkpoint_manifest": str(metadata_manifest),
                "checkpoint_manifest_sha256": (
                    teacher_v3_training._file_sha256(metadata_manifest)
                ),
            },
            "label_policy_sensitivity": {
                "checkpoint_manifest": str(sensitivity_manifest),
                "checkpoint_manifest_sha256": (
                    teacher_v3_training._file_sha256(
                        sensitivity_manifest
                    )
                ),
            },
        }
    )
    freeze_path.write_text(
        json.dumps(freeze_payload, sort_keys=True) + "\n"
    )
    baseline_metrics = json.loads(
        (
            run_dir
            / TEACHER_V3_BASELINE_PROFILE
            / "fixed_test_metrics.json"
        ).read_text()
    )
    primary_metrics = json.loads(
        (
            run_dir
            / TEACHER_V3_PRIMARY_PROFILE
            / "fixed_test_metrics.json"
        ).read_text()
    )
    sensitivity_metrics = json.loads(
        (
            run_dir
            / "sensitivities"
            / "uncertain_masked"
            / TEACHER_V3_PRIMARY_PROFILE
            / "fixed_test_metrics.json"
        ).read_text()
    )
    summary = {
        "schema_version": "twirl_teacher_v3_training_summary_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
        "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
        "fixed_test_status": "opened_once_after_pretest_model_freeze",
        "pretest_model_freeze": str(freeze_path),
        "pretest_model_freeze_sha256": teacher_v3_training._file_sha256(
            freeze_path
        ),
        "bootstrap": {
            "unit": "tic",
            "draws": 50,
            "seed": 123,
            "confidence_level": 0.95,
        },
        "calibration": [
            baseline_calibration,
            primary_calibration,
        ],
        "test_metrics": {
            TEACHER_V3_BASELINE_PROFILE: baseline_metrics,
            TEACHER_V3_PRIMARY_PROFILE: primary_metrics,
        },
        "selected_checkpoint_manifest": str(selected_manifest),
        "selected_checkpoint_manifest_sha256": (
            teacher_v3_training._file_sha256(selected_manifest)
        ),
        "metadata_baseline_checkpoint_manifest": str(metadata_manifest),
        "metadata_baseline_checkpoint_manifest_sha256": (
            teacher_v3_training._file_sha256(metadata_manifest)
        ),
        "native_input_manifest": str(native_manifest),
        "native_manifest_sha256": teacher_v3_training._file_sha256(
            native_manifest
        ),
        "sensitivity_analyses": {
            "uncertain_masked": {
                "label_policy": TEACHER_V3_UNCERTAIN_MASKED_POLICY,
                "training_execution": "independent_five_fold_retraining",
                "calibration": sensitivity_calibration,
                "test_metrics": sensitivity_metrics,
                "checkpoint_manifest": str(sensitivity_manifest),
                "checkpoint_manifest_sha256": (
                    teacher_v3_training._file_sha256(
                        sensitivity_manifest
                    )
                ),
            }
        },
        "automatic_production_promotion": False,
    }
    (run_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return run_dir


def test_teacher_v3_performance_figure_writes_verified_products(
    tmp_path: Path,
) -> None:
    run_dir = _synthetic_run(tmp_path)

    first = build_teacher_v3_performance_figure(
        run_dir=run_dir,
        out_dir=tmp_path / "figure_one",
        expected_draws=50,
        expected_seed=123,
    )
    second = build_teacher_v3_performance_figure(
        run_dir=run_dir,
        out_dir=tmp_path / "figure_two",
        expected_draws=50,
        expected_seed=123,
    )

    first_dir = tmp_path / "figure_one"
    assert (first_dir / "teacher_v3_performance.png").stat().st_size > 0
    assert (first_dir / "teacher_v3_performance.pdf").stat().st_size > 0
    plotted = pd.read_csv(
        first_dir / "teacher_v3_performance_metrics.csv"
    )
    assert set(plotted["model_key"]) == {
        "metadata_baseline",
        "primary",
        "uncertain_masked",
    }
    assert "expected_calibration_error" in set(plotted["metric"])
    confusion = pd.read_csv(
        first_dir / "teacher_v3_primary_confusion_matrix.csv"
    )
    assert confusion["count"].sum() == 8
    assert first["automatic_production_promotion"] is False
    assert teacher_v3_training._file_sha256(
        first_dir / "teacher_v3_performance.png"
    ) == teacher_v3_training._file_sha256(
        tmp_path / "figure_two" / "teacher_v3_performance.png"
    )
    assert first["bootstrap"]["unit"] == second["bootstrap"]["unit"] == "tic"
    assert first["bootstrap"]["seed"] == second["bootstrap"]["seed"] == 123
    assert first["common_support"]["n_rows"] == 7
    assert len(first["common_support"]["review_id_tic_pairs_sha256"]) == 64


def test_teacher_v3_figure_rejects_bootstrap_hash_and_seed_drift(
    tmp_path: Path,
) -> None:
    seed_case = tmp_path / "seed_case"
    seed_case.mkdir()
    run_dir = _synthetic_run(seed_case)
    bootstrap_summary = (
        run_dir
        / TEACHER_V3_PRIMARY_PROFILE
        / "fixed_test_non_uncertain_bootstrap.summary.json"
    )
    payload = json.loads(bootstrap_summary.read_text())
    payload["seed"] = 999
    bootstrap_summary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    with pytest.raises(RuntimeError, match="bootstrap summary is invalid"):
        build_teacher_v3_performance_figure(
            run_dir=run_dir,
            out_dir=seed_case / "figure",
            expected_draws=50,
            expected_seed=123,
        )

    hash_case = tmp_path / "hash_case"
    hash_case.mkdir()
    run_dir = _synthetic_run(hash_case)
    metrics_path = (
        run_dir
        / TEACHER_V3_BASELINE_PROFILE
        / "fixed_test_non_uncertain_bootstrap_metrics.csv"
    )
    metrics_path.write_text(metrics_path.read_text() + "\n")
    with pytest.raises(RuntimeError, match="metrics_sha256"):
        build_teacher_v3_performance_figure(
            run_dir=run_dir,
            out_dir=hash_case / "figure",
            expected_draws=50,
            expected_seed=123,
        )


def test_teacher_v3_figure_rejects_nonidentical_common_support(
    tmp_path: Path,
) -> None:
    run_dir = _synthetic_run(tmp_path)
    sensitivity_dir = (
        run_dir
        / "sensitivities"
        / "uncertain_masked"
        / TEACHER_V3_PRIMARY_PROFILE
    )
    predictions_path = sensitivity_dir / "fixed_test_predictions.csv"
    predictions = pd.read_csv(predictions_path)
    predictions.loc[0, "review_id"] = "different_review_id"
    predictions.to_csv(predictions_path, index=False)
    bootstrap_summary_path = (
        sensitivity_dir / "fixed_test_bootstrap.summary.json"
    )
    bootstrap_summary = json.loads(bootstrap_summary_path.read_text())
    bootstrap_summary["fixed_test_predictions_sha256"] = (
        teacher_v3_training._file_sha256(predictions_path)
    )
    bootstrap_summary_path.write_text(
        json.dumps(bootstrap_summary, indent=2, sort_keys=True) + "\n"
    )

    with pytest.raises(RuntimeError, match="exact real non-uncertain"):
        build_teacher_v3_performance_figure(
            run_dir=run_dir,
            out_dir=tmp_path / "figure",
            expected_draws=50,
            expected_seed=123,
        )


def test_teacher_v3_figure_rejects_freeze_summary_manifest_drift(
    tmp_path: Path,
) -> None:
    run_dir = _synthetic_run(tmp_path)
    freeze_path = run_dir / "pretest_model_freeze.json"
    freeze = json.loads(freeze_path.read_text())
    freeze["primary_checkpoint_manifest_sha256"] = "0" * 64
    freeze_path.write_text(
        json.dumps(freeze, indent=2, sort_keys=True) + "\n"
    )
    summary_path = run_dir / "summary.json"
    summary = json.loads(summary_path.read_text())
    summary["pretest_model_freeze_sha256"] = (
        teacher_v3_training._file_sha256(freeze_path)
    )
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )

    with pytest.raises(RuntimeError, match="not identical in freeze and summary"):
        build_teacher_v3_performance_figure(
            run_dir=run_dir,
            out_dir=tmp_path / "figure",
            expected_draws=50,
            expected_seed=123,
        )
