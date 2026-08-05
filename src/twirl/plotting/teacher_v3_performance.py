"""Deterministic publication products for the frozen Teacher v3 evaluation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.plotting.style import apply_twirl_style, get_ordered_palette
from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_BASELINE_PROFILE,
    TEACHER_V3_BOOTSTRAP_SCHEMA,
    TEACHER_V3_PRIMARY_LABEL_POLICY,
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
    TEACHER_V3_UNCERTAIN_MASKED_POLICY,
)


TEACHER_V3_PERFORMANCE_FIGURE_SCHEMA = (
    "twirl_teacher_v3_performance_figure_v1"
)
DISPLAY_MODELS: tuple[tuple[str, str], ...] = (
    ("metadata_baseline", "Metadata baseline"),
    ("primary", "Primary"),
    ("uncertain_masked", "Primary, uncertain masked"),
)
DISPLAY_METRICS: tuple[tuple[str, str], ...] = (
    ("balanced_accuracy", "Balanced accuracy"),
    ("macro_f1", "Macro F1"),
    ("planet_recall", "Planet recall"),
    ("eb_recall", "EB recall"),
    ("expected_calibration_error", "ECE ↓"),
)
CLASS_LABELS: Mapping[str, str] = {
    "planet_like": "Planet",
    "eclipse_contact": "EB",
    "smooth_variable": "Variable",
    "other": "Other",
}


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _truth(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "1.0", "true", "t", "yes", "y"})
    )


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"expected a JSON mapping: {path}")
    return payload


def _evaluation_artifacts(
    run_dir: Path,
) -> dict[str, dict[str, Path]]:
    primary_dir = run_dir / TEACHER_V3_PRIMARY_PROFILE
    baseline_dir = run_dir / TEACHER_V3_BASELINE_PROFILE
    sensitivity_dir = (
        run_dir
        / "sensitivities"
        / "uncertain_masked"
        / TEACHER_V3_PRIMARY_PROFILE
    )
    return {
        "metadata_baseline": {
            "directory": baseline_dir,
            "metrics": (
                baseline_dir
                / "fixed_test_non_uncertain_bootstrap_metrics.csv"
            ),
            "reliability": (
                baseline_dir / "fixed_test_non_uncertain_reliability.csv"
            ),
            "bootstrap_summary": (
                baseline_dir
                / "fixed_test_non_uncertain_bootstrap.summary.json"
            ),
            "predictions": baseline_dir / "fixed_test_predictions.csv",
            "fixed_metrics": baseline_dir / "fixed_test_metrics.json",
        },
        "primary": {
            "directory": primary_dir,
            "metrics": (
                primary_dir
                / "fixed_test_non_uncertain_bootstrap_metrics.csv"
            ),
            "reliability": (
                primary_dir / "fixed_test_non_uncertain_reliability.csv"
            ),
            "bootstrap_summary": (
                primary_dir
                / "fixed_test_non_uncertain_bootstrap.summary.json"
            ),
            "predictions": primary_dir / "fixed_test_predictions.csv",
            "fixed_metrics": primary_dir / "fixed_test_metrics.json",
        },
        "uncertain_masked": {
            "directory": sensitivity_dir,
            "metrics": sensitivity_dir / "fixed_test_bootstrap_metrics.csv",
            "reliability": sensitivity_dir / "fixed_test_reliability.csv",
            "bootstrap_summary": (
                sensitivity_dir / "fixed_test_bootstrap.summary.json"
            ),
            "predictions": sensitivity_dir / "fixed_test_predictions.csv",
            "fixed_metrics": sensitivity_dir / "fixed_test_metrics.json",
        },
    }


def _validate_run_summary(
    run_dir: Path,
    *,
    expected_draws: int,
    expected_seed: int,
) -> tuple[dict[str, Any], Path]:
    summary_path = run_dir / "summary.json"
    summary = _load_json(summary_path)
    expected = {
        "schema_version": "twirl_teacher_v3_training_summary_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
        "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
        "fixed_test_status": "opened_once_after_pretest_model_freeze",
    }
    failures = [
        f"{name}={summary.get(name)!r}, expected {value!r}"
        for name, value in expected.items()
        if summary.get(name) != value
    ]
    sensitivity_analyses = summary.get("sensitivity_analyses")
    if not isinstance(sensitivity_analyses, dict):
        failures.append("summary sensitivity analyses are missing")
        sensitivity_analyses = {}
    sensitivity = sensitivity_analyses.get("uncertain_masked", {})
    if not isinstance(sensitivity, dict):
        failures.append("summary uncertain-masked sensitivity is invalid")
        sensitivity = {}
    if sensitivity.get("label_policy") != TEACHER_V3_UNCERTAIN_MASKED_POLICY:
        failures.append("summary lacks the uncertain-masked label policy")
    if sensitivity.get("training_execution") != (
        "independent_five_fold_retraining"
    ):
        failures.append(
            "uncertain-masked result is not independent five-fold retraining"
        )
    expected_bootstrap = {
        "unit": "tic",
        "draws": int(expected_draws),
        "seed": int(expected_seed),
        "confidence_level": 0.95,
    }
    bootstrap = summary.get("bootstrap")
    if not isinstance(bootstrap, dict):
        failures.append("summary bootstrap contract is missing")
    else:
        for name, value in expected_bootstrap.items():
            if bootstrap.get(name) != value:
                failures.append(
                    f"summary bootstrap {name}={bootstrap.get(name)!r}, "
                    f"expected {value!r}"
                )

    freeze_path = run_dir / "pretest_model_freeze.json"
    freeze: dict[str, Any] = {}
    if (
        not freeze_path.is_file()
        or _file_sha256(freeze_path)
        != str(summary.get("pretest_model_freeze_sha256", ""))
    ):
        failures.append("pretest model-freeze path/hash is invalid")
    else:
        freeze = _load_json(freeze_path)
        expected_freeze = {
            "schema_version": "twirl_teacher_v3_model_freeze_v1",
            "run_id": TEACHER_V3_RUN_ID,
            "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
            "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
            "profile_selection_policy": "fixed_before_test",
            "test_rows_used_for_selection_or_calibration": False,
        }
        for name, value in expected_freeze.items():
            if freeze.get(name) != value:
                failures.append(
                    f"pretest freeze {name}={freeze.get(name)!r}, "
                    f"expected {value!r}"
                )
        freeze_bootstrap = freeze.get("bootstrap")
        if not isinstance(freeze_bootstrap, dict):
            failures.append("pretest freeze lacks its bootstrap contract")
        else:
            for name, value in expected_bootstrap.items():
                if freeze_bootstrap.get(name) != value:
                    failures.append(
                        f"pretest freeze bootstrap {name} does not match"
                    )

    bound_artifacts = (
        (
            run_dir / "selected_checkpoint_manifest.json",
            summary.get("selected_checkpoint_manifest_sha256", ""),
            "primary checkpoint manifest",
        ),
        (
            run_dir / "metadata_baseline_checkpoint_manifest.json",
            summary.get(
                "metadata_baseline_checkpoint_manifest_sha256",
                "",
            ),
            "metadata checkpoint manifest",
        ),
        (
            run_dir
            / "sensitivities"
            / "uncertain_masked"
            / "uncertain_masked_checkpoint_manifest.json",
            sensitivity.get("checkpoint_manifest_sha256", ""),
            "uncertain-masked checkpoint manifest",
        ),
        (
            run_dir / "native_input_manifest.json",
            summary.get("native_manifest_sha256", ""),
            "native-input manifest",
        ),
    )
    for path, expected_sha256, label in bound_artifacts:
        if (
            not path.is_file()
            or _file_sha256(path) != str(expected_sha256)
        ):
            failures.append(f"{label} path/hash is invalid")

    freeze_sensitivity = freeze.get("label_policy_sensitivity", {})
    if not isinstance(freeze_sensitivity, dict):
        failures.append("pretest freeze sensitivity binding is invalid")
        freeze_sensitivity = {}
    freeze_metadata = freeze.get("metadata_baseline", {})
    if not isinstance(freeze_metadata, dict):
        failures.append("pretest freeze metadata binding is invalid")
        freeze_metadata = {}
    manifest_bindings = (
        (
            run_dir / "selected_checkpoint_manifest.json",
            summary.get("selected_checkpoint_manifest", ""),
            summary.get("selected_checkpoint_manifest_sha256", ""),
            freeze.get("primary_checkpoint_manifest", ""),
            freeze.get("primary_checkpoint_manifest_sha256", ""),
            "primary checkpoint manifest",
        ),
        (
            run_dir / "metadata_baseline_checkpoint_manifest.json",
            summary.get("metadata_baseline_checkpoint_manifest", ""),
            summary.get(
                "metadata_baseline_checkpoint_manifest_sha256",
                "",
            ),
            freeze_metadata.get("checkpoint_manifest", ""),
            freeze_metadata.get("checkpoint_manifest_sha256", ""),
            "metadata checkpoint manifest",
        ),
        (
            run_dir
            / "sensitivities"
            / "uncertain_masked"
            / "uncertain_masked_checkpoint_manifest.json",
            sensitivity.get("checkpoint_manifest", ""),
            sensitivity.get("checkpoint_manifest_sha256", ""),
            freeze_sensitivity.get("checkpoint_manifest", ""),
            freeze_sensitivity.get("checkpoint_manifest_sha256", ""),
            "uncertain-masked checkpoint manifest",
        ),
    )
    for (
        expected_path,
        summary_path_value,
        summary_sha256,
        freeze_path_value,
        freeze_sha256,
        label,
    ) in manifest_bindings:
        try:
            summary_path_matches = (
                Path(str(summary_path_value)).expanduser().resolve()
                == expected_path.resolve()
            )
            freeze_path_matches = (
                Path(str(freeze_path_value)).expanduser().resolve()
                == expected_path.resolve()
            )
        except (OSError, RuntimeError, ValueError):
            summary_path_matches = False
            freeze_path_matches = False
        if not summary_path_matches or not freeze_path_matches:
            failures.append(
                f"{label} path is not identical in freeze and summary"
            )
        observed_sha256 = (
            _file_sha256(expected_path) if expected_path.is_file() else ""
        )
        if (
            str(summary_sha256) != str(freeze_sha256)
            or str(summary_sha256) != observed_sha256
        ):
            failures.append(
                f"{label} hash is not identical in freeze and summary"
            )

    calibration_records = summary.get("calibration")
    if not isinstance(calibration_records, list):
        failures.append("summary calibration records are missing")
        calibration_records = []
    calibration_by_profile = {
        str(record.get("profile", "")): record
        for record in calibration_records
        if isinstance(record, dict)
    }
    expected_calibration_paths = {
        TEACHER_V3_BASELINE_PROFILE: (
            run_dir
            / TEACHER_V3_BASELINE_PROFILE
            / "pooled_oof_calibration.json"
        ),
        TEACHER_V3_PRIMARY_PROFILE: (
            run_dir
            / TEACHER_V3_PRIMARY_PROFILE
            / "pooled_oof_calibration.json"
        ),
    }
    for profile, path in expected_calibration_paths.items():
        record = calibration_by_profile.get(profile, {})
        if (
            not path.is_file()
            or _file_sha256(path)
            != str(record.get("calibration_sha256", ""))
        ):
            failures.append(f"{profile} pooled calibration hash is invalid")
    sensitivity_calibration = sensitivity.get("calibration")
    sensitivity_calibration_path = (
        run_dir
        / "sensitivities"
        / "uncertain_masked"
        / TEACHER_V3_PRIMARY_PROFILE
        / "pooled_oof_calibration.json"
    )
    if not isinstance(sensitivity_calibration, dict) or (
        not sensitivity_calibration_path.is_file()
        or _file_sha256(sensitivity_calibration_path)
        != str(sensitivity_calibration.get("calibration_sha256", ""))
    ):
        failures.append("uncertain-masked pooled calibration hash is invalid")
    if failures:
        raise RuntimeError(
            "Teacher v3 run is not figure-ready: " + "; ".join(failures)
        )
    return summary, summary_path


def _load_interval_products(
    *,
    artifacts: Mapping[str, Mapping[str, Path]],
    expected_draws: int,
    expected_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, Any]]]:
    display_lookup = dict(DISPLAY_MODELS)
    metrics_parts: list[pd.DataFrame] = []
    reliability_parts: list[pd.DataFrame] = []
    bootstrap_summaries: dict[str, dict[str, Any]] = {}
    expected_profiles = {
        "metadata_baseline": TEACHER_V3_BASELINE_PROFILE,
        "primary": TEACHER_V3_PRIMARY_PROFILE,
        "uncertain_masked": TEACHER_V3_PRIMARY_PROFILE,
    }
    expected_policies = {
        "metadata_baseline": TEACHER_V3_PRIMARY_LABEL_POLICY,
        "primary": TEACHER_V3_PRIMARY_LABEL_POLICY,
        "uncertain_masked": TEACHER_V3_UNCERTAIN_MASKED_POLICY,
    }
    for model_key, paths in artifacts.items():
        bootstrap_summary = _load_json(paths["bootstrap_summary"])
        expected_summary = {
            "schema_version": TEACHER_V3_BOOTSTRAP_SCHEMA,
            "run_id": TEACHER_V3_RUN_ID,
            "profile": expected_profiles[model_key],
            "label_policy": expected_policies[model_key],
            "evaluation_scope": "fixed_test_real_non_uncertain",
            "bootstrap_unit": "tic",
            "draws": int(expected_draws),
            "seed": int(expected_seed),
            "confidence_level": 0.95,
        }
        summary_failures = [
            f"{name}={bootstrap_summary.get(name)!r}, expected {value!r}"
            for name, value in expected_summary.items()
            if bootstrap_summary.get(name) != value
        ]
        hash_bindings = (
            ("metrics", "metrics_sha256"),
            ("reliability", "reliability_sha256"),
            ("predictions", "fixed_test_predictions_sha256"),
        )
        for path_name, digest_name in hash_bindings:
            if _file_sha256(paths[path_name]) != str(
                bootstrap_summary.get(digest_name, "")
            ):
                summary_failures.append(
                    f"{digest_name} does not bind {path_name}"
                )
        if summary_failures:
            raise RuntimeError(
                f"{model_key} bootstrap summary is invalid: "
                + "; ".join(summary_failures)
            )
        metrics = pd.read_csv(paths["metrics"])
        reliability = pd.read_csv(paths["reliability"])
        for name, frame in (
            ("metrics", metrics),
            ("reliability", reliability),
        ):
            if frame.empty:
                raise ValueError(f"{model_key} {name} table is empty")
            if set(pd.to_numeric(frame["draws"], errors="raise")) != {
                int(expected_draws)
            }:
                raise ValueError(
                    f"{model_key} {name} does not use {expected_draws} draws"
                )
            if set(frame["bootstrap_unit"].astype(str)) != {"tic"}:
                raise ValueError(
                    f"{model_key} {name} is not TIC-clustered"
                )
            if set(
                pd.to_numeric(
                    frame["confidence_level"],
                    errors="raise",
                )
            ) != {0.95}:
                raise ValueError(
                    f"{model_key} {name} does not report 95% intervals"
                )
            if set(frame["evaluation_scope"].astype(str)) != {
                "fixed_test_real_non_uncertain"
            }:
                raise ValueError(
                    f"{model_key} {name} is not on common real test support"
                )
            if set(pd.to_numeric(frame["seed"], errors="raise")) != {
                int(expected_seed)
            }:
                raise ValueError(
                    f"{model_key} {name} does not use seed {expected_seed}"
                )
            if set(frame["profile"].astype(str)) != {
                expected_profiles[model_key]
            }:
                raise ValueError(
                    f"{model_key} {name} has the wrong profile"
                )
            if set(frame["label_policy"].astype(str)) != {
                expected_policies[model_key]
            }:
                raise ValueError(
                    f"{model_key} {name} has the wrong label policy"
                )
        expected_rows = int(bootstrap_summary["n_rows"])
        expected_tics = int(bootstrap_summary["n_tics"])
        if set(pd.to_numeric(metrics["n_rows"], errors="raise")) != {
            expected_rows
        } or set(pd.to_numeric(metrics["n_tics"], errors="raise")) != {
            expected_tics
        }:
            raise ValueError(
                f"{model_key} metric support disagrees with its bootstrap "
                "summary"
            )
        metrics["model_key"] = model_key
        metrics["display_model"] = display_lookup[model_key]
        reliability["model_key"] = model_key
        reliability["display_model"] = display_lookup[model_key]
        metrics_parts.append(metrics)
        reliability_parts.append(reliability)
        bootstrap_summaries[model_key] = bootstrap_summary
    metrics = pd.concat(metrics_parts, ignore_index=True)
    reliability = pd.concat(reliability_parts, ignore_index=True)
    required_metrics = {name for name, _ in DISPLAY_METRICS}
    for model_key in dict(DISPLAY_MODELS):
        observed = set(
            metrics.loc[metrics["model_key"].eq(model_key), "metric"]
            .astype(str)
        )
        missing = sorted(required_metrics - observed)
        if missing:
            raise ValueError(
                f"{model_key} lacks plotted metrics: {missing}"
            )
    return metrics, reliability, bootstrap_summaries


def _support_digest(values: list[str]) -> str:
    payload = "\n".join(sorted(values)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _validate_common_non_uncertain_support(
    *,
    artifacts: Mapping[str, Mapping[str, Path]],
    bootstrap_summaries: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Prove that all three interval products use the identical real support."""

    support_by_model: dict[str, dict[str, Any]] = {}
    reference_pairs: set[tuple[str, int]] | None = None
    for model_key, paths in artifacts.items():
        predictions = pd.read_csv(paths["predictions"], low_memory=False)
        required = {
            "review_id",
            "tic",
            "human_label",
            "is_injected_row",
            "fixed_split",
            "morphology_target_index",
        }
        missing = sorted(required - set(predictions.columns))
        if missing:
            raise KeyError(
                f"{model_key} fixed-test predictions lack support columns: "
                f"{missing}"
            )
        if predictions["review_id"].fillna("").astype(str).duplicated().any():
            raise ValueError(
                f"{model_key} fixed-test review IDs are not unique"
            )
        split = (
            predictions["fixed_split"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
        )
        truth = pd.to_numeric(
            predictions["morphology_target_index"],
            errors="raise",
        ).to_numpy(dtype=int)
        active = (
            split.eq("test").to_numpy()
            & (~_truth(predictions["is_injected_row"]).to_numpy())
            & (
                ~predictions["human_label"]
                .fillna("")
                .astype(str)
                .str.strip()
                .eq("uncertain")
                .to_numpy()
            )
            & (truth >= 0)
        )
        selected = predictions.loc[active, ["review_id", "tic"]].copy()
        selected["review_id"] = (
            selected["review_id"].fillna("").astype(str).str.strip()
        )
        selected["tic"] = pd.to_numeric(
            selected["tic"],
            errors="raise",
        ).astype(np.int64)
        if selected.empty or selected["review_id"].eq("").any():
            raise ValueError(
                f"{model_key} has no valid common-support review rows"
            )
        pairs = {
            (str(row.review_id), int(row.tic))
            for row in selected.itertuples(index=False)
        }
        if len(pairs) != len(selected):
            raise ValueError(
                f"{model_key} common support contains duplicate review/TIC pairs"
            )
        summary = bootstrap_summaries[model_key]
        if (
            int(summary.get("n_rows", -1)) != len(pairs)
            or int(summary.get("n_tics", -1))
            != selected["tic"].nunique()
        ):
            raise RuntimeError(
                f"{model_key} bootstrap support counts do not match "
                "fixed-test predictions"
            )
        if reference_pairs is None:
            reference_pairs = pairs
        elif pairs != reference_pairs:
            missing_pairs = sorted(reference_pairs - pairs)[:5]
            extra_pairs = sorted(pairs - reference_pairs)[:5]
            raise RuntimeError(
                "Teacher v3 figure inputs do not share the exact real "
                "non-uncertain review-ID/TIC support; "
                f"{model_key} missing={missing_pairs}, extra={extra_pairs}"
            )
        review_ids = sorted(review_id for review_id, _ in pairs)
        unique_tics = sorted({tic for _, tic in pairs})
        pair_values = [
            f"{review_id}\t{tic}"
            for review_id, tic in sorted(pairs)
        ]
        support_by_model[model_key] = {
            "n_rows": int(len(pairs)),
            "n_tics": int(len(unique_tics)),
            "review_ids_sha256": _support_digest(review_ids),
            "tics_sha256": _support_digest(
                [str(tic) for tic in unique_tics]
            ),
            "review_id_tic_pairs_sha256": _support_digest(pair_values),
        }
    assert reference_pairs is not None
    first = support_by_model[next(iter(dict(DISPLAY_MODELS)))]
    for model_key, support in support_by_model.items():
        if support != first:
            raise RuntimeError(
                f"{model_key} common-support digests disagree"
            )
    return {
        "evaluation_scope": "fixed_test_real_non_uncertain",
        **first,
        "verified_models": list(dict(DISPLAY_MODELS)),
    }


def _primary_confusion(
    predictions_path: Path,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    predictions = pd.read_csv(predictions_path, low_memory=False)
    required = {
        "is_injected_row",
        "morphology_target_index",
        "morphology_prediction_index",
    }
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise KeyError(f"primary fixed-test predictions lack: {missing}")
    real = ~_truth(predictions["is_injected_row"]).to_numpy()
    truth = pd.to_numeric(
        predictions["morphology_target_index"],
        errors="raise",
    ).to_numpy(dtype=int)
    estimated = pd.to_numeric(
        predictions["morphology_prediction_index"],
        errors="raise",
    ).to_numpy(dtype=int)
    active = real & (truth >= 0)
    truth = truth[active]
    estimated = estimated[active]
    if not len(truth):
        raise ValueError("primary result has no real fixed-test targets")
    n_classes = len(MORPHOLOGY_CLASSES)
    matrix = np.zeros((n_classes, n_classes), dtype=int)
    for actual, prediction in zip(truth, estimated):
        if not 0 <= int(prediction) < n_classes:
            raise ValueError("primary predictions contain an invalid class")
        matrix[int(actual), int(prediction)] += 1
    support = matrix.sum(axis=1)
    normalized = np.divide(
        matrix,
        support[:, None],
        out=np.full_like(matrix, np.nan, dtype=float),
        where=support[:, None] > 0,
    )
    records = [
        {
            "actual_class": actual_label,
            "predicted_class": predicted_label,
            "count": int(matrix[actual_index, predicted_index]),
            "actual_class_support": int(support[actual_index]),
            "row_fraction": float(
                normalized[actual_index, predicted_index]
            ),
            "evaluation_scope": "fixed_test_real_labels",
            "model_key": "primary",
        }
        for actual_index, actual_label in enumerate(MORPHOLOGY_CLASSES)
        for predicted_index, predicted_label in enumerate(MORPHOLOGY_CLASSES)
    ]
    return pd.DataFrame(records), matrix, normalized


def _flatten_source_metrics(
    artifacts: Mapping[str, Mapping[str, Path]],
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for model_key, paths in artifacts.items():
        payload = _load_json(paths["fixed_metrics"])
        morphology = payload["morphology_by_source"]
        calibration = payload["calibration"]
        for source_scope in ("all", "real", "injected"):
            values = morphology[source_scope]
            base = {
                "model_key": model_key,
                "source_scope": source_scope,
                "n": int(values["n"]),
                "accuracy": float(values["accuracy"]),
                "balanced_accuracy": float(values["balanced_accuracy"]),
                "macro_f1": float(values["macro_f1"]),
                "expected_calibration_error": float(
                    calibration[source_scope]["ece"]
                ),
            }
            for class_name in MORPHOLOGY_CLASSES:
                class_metrics = values["per_class"].get(class_name, {})
                for metric in ("n", "precision", "recall", "f1"):
                    base[f"{class_name}_{metric}"] = class_metrics.get(
                        metric,
                        float("nan"),
                    )
            records.append(base)
    return pd.DataFrame(records)


def _metric_support_label(
    metric_row: pd.Series,
    label: str,
) -> str:
    support = (
        f"n={int(metric_row['support_rows'])}; "
        f"{int(metric_row['support_tics'])} TICs"
    )
    warning_text = str(metric_row["support_warning"]).strip().lower()
    warning = (
        " †" if warning_text not in {"", "nan", "none"} else ""
    )
    return f"{label}{warning}\n({support})"


def _load_verified_figure_inputs(
    *,
    run_dir: Path,
    expected_draws: int,
    expected_seed: int,
) -> tuple[
    dict[str, Any],
    Path,
    dict[str, dict[str, Path]],
    pd.DataFrame,
    pd.DataFrame,
    dict[str, Any],
    dict[str, str],
]:
    summary, summary_path = _validate_run_summary(
        run_dir,
        expected_draws=expected_draws,
        expected_seed=expected_seed,
    )
    artifacts = _evaluation_artifacts(run_dir)
    for paths in artifacts.values():
        for name, path in paths.items():
            if name != "directory" and not path.is_file():
                raise FileNotFoundError(path)
    summary_test_metrics = summary.get("test_metrics")
    if not isinstance(summary_test_metrics, dict):
        raise RuntimeError(
            "completed Teacher v3 summary lacks fixed-test metrics"
        )
    sensitivity_payload = summary["sensitivity_analyses"]["uncertain_masked"]
    sensitivity_test_metrics = sensitivity_payload.get("test_metrics")
    if not isinstance(sensitivity_test_metrics, dict):
        raise RuntimeError(
            "completed Teacher v3 summary lacks uncertain-masked test metrics"
        )
    expected_fixed_metrics = {
        "metadata_baseline": summary_test_metrics.get(
            TEACHER_V3_BASELINE_PROFILE,
            {},
        ),
        "primary": summary_test_metrics.get(
            TEACHER_V3_PRIMARY_PROFILE,
            {},
        ),
        "uncertain_masked": sensitivity_test_metrics,
    }
    for model_key, expected_payload in expected_fixed_metrics.items():
        observed_payload = _load_json(
            artifacts[model_key]["fixed_metrics"]
        )
        expected_digest = hashlib.sha256(
            json.dumps(
                expected_payload,
                sort_keys=True,
                allow_nan=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()
        observed_digest = hashlib.sha256(
            json.dumps(
                observed_payload,
                sort_keys=True,
                allow_nan=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()
        if expected_digest != observed_digest:
            raise RuntimeError(
                f"{model_key} fixed-test metrics do not match the completed "
                "training summary"
            )
    metrics, reliability, bootstrap_summaries = _load_interval_products(
        artifacts=artifacts,
        expected_draws=expected_draws,
        expected_seed=expected_seed,
    )
    support = _validate_common_non_uncertain_support(
        artifacts=artifacts,
        bootstrap_summaries=bootstrap_summaries,
    )
    input_paths = {
        "training_summary": summary_path,
        **{
            f"{model_key}_{name}": path
            for model_key, paths in artifacts.items()
            for name, path in paths.items()
            if name != "directory"
        },
    }
    input_hashes = {
        name: _file_sha256(path)
        for name, path in sorted(input_paths.items())
    }
    return (
        summary,
        summary_path,
        artifacts,
        metrics,
        reliability,
        support,
        input_hashes,
    )


def validate_teacher_v3_performance_inputs(
    *,
    run_dir: Path,
    expected_draws: int = 2000,
    expected_seed: int = 560063,
) -> dict[str, Any]:
    """Validate the completed full run without creating figure products."""

    run_dir = Path(run_dir).resolve()
    (
        _summary,
        summary_path,
        _artifacts,
        _metrics,
        _reliability,
        support,
        input_hashes,
    ) = _load_verified_figure_inputs(
        run_dir=run_dir,
        expected_draws=int(expected_draws),
        expected_seed=int(expected_seed),
    )
    return {
        "schema_version": "twirl_teacher_v3_performance_input_audit_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "run_dir": str(run_dir),
        "summary_path": str(summary_path),
        "summary_sha256": input_hashes["training_summary"],
        "bootstrap": {
            "unit": "tic",
            "draws": int(expected_draws),
            "seed": int(expected_seed),
            "confidence_level": 0.95,
        },
        "common_support": support,
        "n_verified_input_artifacts": int(len(input_hashes)),
        "verified": True,
    }


def build_teacher_v3_performance_figure(
    *,
    run_dir: Path,
    out_dir: Path,
    expected_draws: int = 2000,
    expected_seed: int = 560063,
) -> dict[str, Any]:
    """Render the fixed-test comparison and write every plotted value."""

    import matplotlib.pyplot as plt

    run_dir = Path(run_dir).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    (
        summary,
        summary_path,
        artifacts,
        metrics,
        reliability,
        common_support,
        verified_input_hashes,
    ) = _load_verified_figure_inputs(
        run_dir=run_dir,
        expected_draws=int(expected_draws),
        expected_seed=int(expected_seed),
    )
    confusion, confusion_counts, confusion_fraction = _primary_confusion(
        artifacts["primary"]["predictions"]
    )
    source_metrics = _flatten_source_metrics(artifacts)

    selected_metrics = metrics.loc[
        metrics["metric"].isin({name for name, _ in DISPLAY_METRICS})
    ].copy()
    metric_order = {name: index for index, (name, _) in enumerate(DISPLAY_METRICS)}
    model_order = {
        name: index for index, (name, _) in enumerate(DISPLAY_MODELS)
    }
    selected_metrics["_metric_order"] = selected_metrics["metric"].map(
        metric_order
    )
    selected_metrics["_model_order"] = selected_metrics["model_key"].map(
        model_order
    )
    selected_metrics = selected_metrics.sort_values(
        ["_metric_order", "_model_order"],
        kind="stable",
    ).drop(columns=["_metric_order", "_model_order"])

    metrics_path = out_dir / "teacher_v3_performance_metrics.csv"
    reliability_path = out_dir / "teacher_v3_performance_reliability.csv"
    confusion_path = out_dir / "teacher_v3_primary_confusion_matrix.csv"
    source_metrics_path = out_dir / "teacher_v3_source_metrics.csv"
    selected_metrics.to_csv(
        metrics_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    reliability.to_csv(
        reliability_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    confusion.to_csv(
        confusion_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    source_metrics.to_csv(
        source_metrics_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )

    template = apply_twirl_style("full_page")
    colors = get_ordered_palette(3, "viridis")
    fig, (metric_axis, confusion_axis) = plt.subplots(
        1,
        2,
        figsize=template["figsize"],
        gridspec_kw={"width_ratios": [1.45, 1.0]},
    )
    y_base = np.arange(len(DISPLAY_METRICS), dtype=float)[::-1]
    offsets = (-0.18, 0.0, 0.18)
    markers = ("s", "o", "^")
    for model_index, (model_key, display_name) in enumerate(DISPLAY_MODELS):
        part = (
            selected_metrics.loc[
                selected_metrics["model_key"].eq(model_key)
            ]
            .set_index("metric")
            .loc[[name for name, _ in DISPLAY_METRICS]]
        )
        point = part["point_estimate"].to_numpy(dtype=float)
        low = part["ci_low"].to_numpy(dtype=float)
        high = part["ci_high"].to_numpy(dtype=float)
        xerr = np.vstack(
            [
                np.maximum(point - low, 0.0),
                np.maximum(high - point, 0.0),
            ]
        )
        metric_axis.errorbar(
            point,
            y_base + offsets[model_index],
            xerr=xerr,
            fmt=markers[model_index],
            markersize=4.2,
            markerfacecolor=(
                "white" if model_key == "metadata_baseline" else colors[model_index]
            ),
            markeredgecolor=colors[model_index],
            color=colors[model_index],
            elinewidth=1.0,
            capsize=2.0,
            label=display_name,
            zorder=3 + model_index,
        )
    support_source = (
        selected_metrics.loc[selected_metrics["model_key"].eq("primary")]
        .set_index("metric")
    )
    metric_axis.set_yticks(y_base)
    metric_axis.set_yticklabels(
        [
            _metric_support_label(support_source.loc[name], label)
            for name, label in DISPLAY_METRICS
        ]
    )
    metric_axis.axhline(
        y_base[-1] + 0.5,
        color="0.75",
        linewidth=0.7,
        linestyle="--",
        zorder=1,
    )
    metric_axis.set_xlim(-0.02, 1.02)
    metric_axis.set_xlabel("Metric value")
    metric_axis.set_title(
        "(a) Common non-uncertain test support",
        loc="left",
    )
    metric_axis.legend(
        loc="lower left",
        bbox_to_anchor=(0.0, 1.08),
        ncol=3,
        frameon=True,
        borderaxespad=0.0,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    metric_axis.text(
        0.0,
        -0.16,
        "95% TIC-cluster percentile intervals "
        f"({int(expected_draws):,} draws). "
        "† n<20: descriptive only.",
        transform=metric_axis.transAxes,
        ha="left",
        va="top",
        fontsize=template["annotation_size"],
    )

    image = confusion_axis.imshow(
        confusion_fraction,
        vmin=0.0,
        vmax=1.0,
        cmap="viridis",
        aspect="equal",
    )
    class_ticks = [
        CLASS_LABELS[class_name] for class_name in MORPHOLOGY_CLASSES
    ]
    support = confusion_counts.sum(axis=1)
    confusion_axis.set_xticks(np.arange(len(MORPHOLOGY_CLASSES)))
    confusion_axis.set_xticklabels(class_ticks, rotation=35, ha="right")
    confusion_axis.set_yticks(np.arange(len(MORPHOLOGY_CLASSES)))
    confusion_axis.set_yticklabels(
        [
            f"{label} (n={int(n)})"
            for label, n in zip(class_ticks, support)
        ]
    )
    confusion_axis.set_xlabel("Predicted")
    confusion_axis.set_ylabel("Human label")
    confusion_axis.set_title(
        "(b) Primary: all real test labels",
        loc="left",
    )
    for row in range(len(MORPHOLOGY_CLASSES)):
        for column in range(len(MORPHOLOGY_CLASSES)):
            fraction = confusion_fraction[row, column]
            count = confusion_counts[row, column]
            text_color = "white" if np.isfinite(fraction) and fraction > 0.55 else "black"
            label = (
                "—"
                if not np.isfinite(fraction)
                else f"{fraction:.2f}\n({int(count)})"
            )
            confusion_axis.text(
                column,
                row,
                label,
                ha="center",
                va="center",
                color=text_color,
                fontsize=template["annotation_size"],
            )
    colorbar = fig.colorbar(
        image,
        ax=confusion_axis,
        fraction=0.046,
        pad=0.04,
    )
    colorbar.set_label("Fraction within human-label class")
    fig.subplots_adjust(
        left=0.15,
        right=0.97,
        bottom=0.22,
        top=0.78,
        wspace=0.34,
    )

    png_path = out_dir / "teacher_v3_performance.png"
    pdf_path = out_dir / "teacher_v3_performance.pdf"
    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
        metadata={"Software": "TWIRL"},
    )
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        metadata={
            "Creator": "TWIRL",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)

    input_paths = {
        "training_summary": summary_path,
        **{
            f"{model_key}_{name}": path
            for model_key, paths in artifacts.items()
            for name, path in paths.items()
            if name != "directory"
        },
    }
    current_input_hashes = {
        name: _file_sha256(path)
        for name, path in sorted(input_paths.items())
    }
    if current_input_hashes != verified_input_hashes:
        raise RuntimeError(
            "Teacher v3 figure input changed during rendering"
        )
    output_paths = {
        "png": png_path,
        "pdf": pdf_path,
        "metrics": metrics_path,
        "reliability": reliability_path,
        "confusion": confusion_path,
        "source_metrics": source_metrics_path,
    }
    provenance = {
        "schema_version": TEACHER_V3_PERFORMANCE_FIGURE_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "figure_question": (
            "How does the fixed Teacher v3 primary compare with metadata-only "
            "and genuinely retrained uncertain-masked controls?"
        ),
        "panel_a_scope": (
            "real non-uncertain fixed-test observations shared by all three "
            "models"
        ),
        "panel_b_scope": (
            "all active real fixed-test labels under the accepted "
            "uncertain-as-other primary policy"
        ),
        "bootstrap": {
            "unit": "tic",
            "draws": int(expected_draws),
            "seed": int(expected_seed),
            "confidence_level": 0.95,
            "interval": "percentile",
        },
        "common_support": common_support,
        "small_support_policy": (
            "Class estimates with fewer than 20 observations are marked with "
            "a dagger and described as descriptive only."
        ),
        "style": {
            "authority": "twirl.plotting.style",
            "template": "full_page",
            "ordered_palette": "viridis",
        },
        "input_artifacts": {
            name: {
                "path": str(path),
                "sha256": verified_input_hashes[name],
            }
            for name, path in sorted(input_paths.items())
        },
        "output_artifacts": {
            name: {
                "path": str(path),
                "sha256": _file_sha256(path),
            }
            for name, path in sorted(output_paths.items())
        },
        "automatic_production_promotion": bool(
            summary.get("automatic_production_promotion", False)
        ),
    }
    provenance_path = out_dir / "teacher_v3_performance.provenance.json"
    provenance_path.write_text(
        json.dumps(
            provenance,
            indent=2,
            sort_keys=True,
            allow_nan=True,
        )
        + "\n"
    )
    return {
        **provenance,
        "provenance_path": str(provenance_path),
        "provenance_sha256": _file_sha256(provenance_path),
    }


__all__ = [
    "TEACHER_V3_PERFORMANCE_FIGURE_SCHEMA",
    "build_teacher_v3_performance_figure",
    "validate_teacher_v3_performance_inputs",
]
