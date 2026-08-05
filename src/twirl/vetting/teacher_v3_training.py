"""Seven-sector Teacher v3 training on observation-addressed native inputs.

Teacher v3 is a new data release, not a new neural-network architecture.  It
retains :data:`s56_harmonic_cnn_v1`, fixes the primary and metadata-baseline
profiles before opening the test partition, and fits one temperature per
profile from the concatenated five-fold development OOF logits.
"""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    HARMONIC_CLASSES,
    MODEL_VERSION,
    MORPHOLOGY_CLASSES,
    PRESERVE_CLASSES,
    HarmonicModelConfig,
    HarmonicTrainConfig,
    build_harmonic_cnn,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    MetadataNormalization,
    build_metadata_matrix,
    prepare_harmonic_training_rows,
)
from twirl.vetting.harmonic_inputs import (
    CHANNEL_CONTRACT,
    RAW_PAIR_CONTRACT_VERSION,
    verify_raw_pair_contract,
)
from twirl.vetting.harmonic_training import (
    _calibration_by_source,
    _evaluate_loader,
    _file_sha256,
    _loader,
    _softmax,
    _subset_metrics,
    _train_one_fold,
    classification_metrics,
    expected_calibration_error,
    fit_temperature,
    injection_truth_human_audit,
    validate_teacher_input_provenance,
)
from twirl.vetting.teacher_native_registry import (
    file_sha256,
    read_table,
    validate_native_input_registry,
    validate_native_input_registry_path,
)
from twirl.vetting.teacher_v3 import (
    TEACHER_V3_CORPUS_VERSION,
    TEACHER_V3_RUN_NAME,
)


TEACHER_V3_RUN_ID = "teacher_v3_s56_s62_a2v1_current_adp"
TEACHER_V3_CHECKPOINT_NAMESPACE = (
    "twirl_teacher_v3_s56_s62_native_v1"
)
TEACHER_V3_METADATA_CHECKPOINT_NAMESPACE = (
    "twirl_teacher_v3_s56_s62_metadata_baseline_v1"
)
TEACHER_V3_NATIVE_MANIFEST_SCHEMA = (
    "twirl_teacher_v3_native_input_manifest_v1"
)
TEACHER_V3_CALIBRATION_SCHEMA = (
    "twirl_teacher_v3_pooled_oof_calibration_v1"
)
TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA = (
    "twirl_teacher_v3_selected_checkpoint_manifest_v1"
)
TEACHER_V3_BOOTSTRAP_SCHEMA = (
    "twirl_teacher_v3_tic_cluster_bootstrap_v1"
)
TEACHER_V3_UNCERTAIN_MASKED_CHECKPOINT_NAMESPACE = (
    "twirl_teacher_v3_s56_s62_native_uncertain_masked_v1"
)
TEACHER_V3_PRIMARY_PROFILE = "shape_plus_periodogram_bls"
TEACHER_V3_BASELINE_PROFILE = "metadata_only"
TEACHER_V3_PRIMARY_LABEL_POLICY = "uncertain_as_other"
TEACHER_V3_UNCERTAIN_MASKED_POLICY = "uncertain_masked_retrained"
TEACHER_V3_PROFILES: tuple[str, str] = (
    TEACHER_V3_BASELINE_PROFILE,
    TEACHER_V3_PRIMARY_PROFILE,
)
TEACHER_V3_SECTORS: tuple[int, ...] = tuple(range(56, 63))
TEACHER_V3_BOOTSTRAP_METRICS: tuple[str, ...] = (
    "macro_f1",
    "balanced_accuracy",
    "planet_recall",
    "eb_recall",
    "variable_recall",
    "other_recall",
    "expected_calibration_error",
)


def _json_payload(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
            allow_nan=True,
        )
        + "\n"
    ).encode("utf-8")


def _write_json(path: Path, payload: Mapping[str, Any]) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    serialized = _json_payload(payload)
    path.write_bytes(serialized)
    return hashlib.sha256(serialized).hexdigest()


def require_fresh_teacher_v3_output_dir(
    path: Path,
    *,
    artifact: str,
) -> Path:
    """Create or accept only an empty, non-symlink output directory."""

    requested = Path(path).expanduser()
    if requested.is_symlink():
        raise RuntimeError(
            f"{artifact} output may not be a symlink: {requested}"
        )
    path = requested.resolve()
    if path.exists():
        if not path.is_dir():
            raise NotADirectoryError(path)
        first = next(path.iterdir(), None)
        if first is not None:
            raise FileExistsError(
                f"{artifact} requires a fresh empty output directory; "
                f"found prior content: {first}"
            )
    else:
        path.mkdir(parents=True, exist_ok=False)
    return path


def _read_training_table_with_stable_hash(
    path: Path,
    *,
    artifact: str,
) -> tuple[pd.DataFrame, str, dict[str, Any]]:
    """Read one frozen table while proving that its bytes stayed unchanged."""

    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    stat_before = path.stat()
    sha256_before = _file_sha256(path)
    frame = pd.read_csv(path, low_memory=False)
    sha256_after = _file_sha256(path)
    stat_after = path.stat()
    if sha256_before != sha256_after:
        raise RuntimeError(
            f"{artifact} changed while it was initially read: {path}"
        )
    audit = {
        "path": str(path),
        "sha256_before_read": sha256_before,
        "sha256_after_read": sha256_after,
        "size_bytes_before_read": int(stat_before.st_size),
        "size_bytes_after_read": int(stat_after.st_size),
        "mtime_ns_before_read": int(stat_before.st_mtime_ns),
        "mtime_ns_after_read": int(stat_after.st_mtime_ns),
        "stable_during_initial_read": True,
    }
    return frame, sha256_before, audit


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


def _string_list_sha256(values: Sequence[str]) -> str:
    payload = "\n".join(sorted(str(value) for value in values)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def prepare_teacher_v3_uncertain_masked_rows(
    rows: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Remove every ``uncertain`` row before sensitivity-model fitting.

    This is deliberately a row-level training-policy change, not a post-hoc
    metric mask.  Removed rows cannot affect metadata normalization, encoder
    pretraining, morphology/preserve/harmonic losses, OOF calibration, or the
    fixed-test sensitivity estimate.  The immutable TIC split assignments of
    retained rows are never changed.
    """

    required = {"review_id", "human_label", "fixed_split", "cv_fold", "tic"}
    missing = sorted(required - set(rows.columns))
    if missing:
        raise KeyError(
            "Teacher v3 uncertain-masked rows are missing columns: "
            f"{missing}"
        )
    uncertain = (
        rows["human_label"].fillna("").astype(str).str.strip().eq("uncertain")
    )
    if not uncertain.any():
        raise ValueError(
            "Teacher v3 uncertain-masked sensitivity found no uncertain rows"
        )
    masked = rows.loc[uncertain].copy()
    retained = rows.loc[~uncertain].copy().reset_index(drop=True)
    split = retained["fixed_split"].fillna("").astype(str).str.strip()
    folds = pd.to_numeric(retained["cv_fold"], errors="raise").astype(int)
    if set(split) != {"development", "test"}:
        raise ValueError(
            "uncertain-masked sensitivity must retain development and test rows"
        )
    if set(folds[split.eq("development")]) != set(range(5)):
        raise ValueError(
            "uncertain-masked sensitivity must retain all five development folds"
        )
    if not folds[split.eq("test")].eq(-1).all():
        raise ValueError(
            "uncertain-masked fixed-test rows must retain cv_fold=-1"
        )
    masked_injected = (
        _truth(masked["is_injected_row"])
        if "is_injected_row" in masked
        else pd.Series(False, index=masked.index)
    )
    audit = {
        "label_policy": TEACHER_V3_UNCERTAIN_MASKED_POLICY,
        "operation": "rows_removed_before_all_model_fitting_and_evaluation",
        "n_input_rows": int(len(rows)),
        "n_retained_rows": int(len(retained)),
        "n_masked_rows": int(uncertain.sum()),
        "n_masked_development_rows": int(
            masked["fixed_split"].astype(str).eq("development").sum()
        ),
        "n_masked_fixed_test_rows": int(
            masked["fixed_split"].astype(str).eq("test").sum()
        ),
        "n_masked_tics": int(
            pd.to_numeric(masked["tic"], errors="raise").nunique()
        ),
        "n_masked_real_rows": int((~masked_injected).sum()),
        "n_masked_injected_rows": int(masked_injected.sum()),
        "masked_rows_by_sector": {
            str(int(key)): int(value)
            for key, value in pd.to_numeric(
                masked.get("sector", pd.Series(dtype=float)),
                errors="coerce",
            )
            .dropna()
            .astype(int)
            .value_counts()
            .sort_index()
            .items()
        },
        "masked_review_ids_sha256": _string_list_sha256(
            masked["review_id"].astype(str).tolist()
        ),
    }
    return retained, audit


def build_teacher_v3_uncertain_masked_provenance(
    *,
    input_provenance: Mapping[str, str],
    sensitivity_audit: Mapping[str, Any],
) -> dict[str, str]:
    """Create a distinct checkpoint namespace for the retrained sensitivity."""

    required = {
        "native_h5_sha256",
        "native_manifest_sha256",
        "training_table_sha256",
    }
    missing = sorted(required - set(input_provenance))
    if missing:
        raise KeyError(
            "Teacher v3 sensitivity provenance is missing fields: "
            f"{missing}"
        )
    return {
        **{str(key): str(value) for key, value in input_provenance.items()},
        "checkpoint_namespace": (
            TEACHER_V3_UNCERTAIN_MASKED_CHECKPOINT_NAMESPACE
        ),
        "label_policy": TEACHER_V3_UNCERTAIN_MASKED_POLICY,
        "masked_review_ids_sha256": str(
            sensitivity_audit["masked_review_ids_sha256"]
        ),
    }


def _morphology_metric_values(
    truth: np.ndarray,
    probability: np.ndarray,
    *,
    require_supported_classes: Sequence[int] = (),
) -> dict[str, float]:
    classification = classification_metrics(
        truth,
        probability,
        classes=MORPHOLOGY_CLASSES,
    )
    per_class = classification["per_class"]
    supported = all(
        int(per_class.get(MORPHOLOGY_CLASSES[index], {}).get("n", 0)) > 0
        for index in require_supported_classes
    )
    return {
        "macro_f1": (
            float(classification["macro_f1"])
            if supported
            else float("nan")
        ),
        "balanced_accuracy": (
            float(classification["balanced_accuracy"])
            if supported
            else float("nan")
        ),
        "planet_recall": float(
            per_class.get("planet_like", {}).get("recall", float("nan"))
        ),
        "eb_recall": float(
            per_class.get("eclipse_contact", {}).get(
                "recall",
                float("nan"),
            )
        ),
        "variable_recall": float(
            per_class.get("smooth_variable", {}).get(
                "recall",
                float("nan"),
            )
        ),
        "other_recall": float(
            per_class.get("other", {}).get("recall", float("nan"))
        ),
        "expected_calibration_error": float(
            expected_calibration_error(truth, probability)["ece"]
        ),
    }


def teacher_v3_tic_cluster_bootstrap(
    *,
    rows: pd.DataFrame,
    truth: np.ndarray,
    probability: np.ndarray,
    profile: str,
    label_policy: str,
    draws: int = 2000,
    seed: int = 560063,
    confidence_level: float = 0.95,
    n_calibration_bins: int = 10,
    evaluation_scope: str = "fixed_test_real_labels",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Estimate fixed-test real-label uncertainty by resampling whole TICs."""

    required = {"review_id", "tic", "is_injected_row", "fixed_split"}
    missing = sorted(required - set(rows.columns))
    if missing:
        raise KeyError(
            f"Teacher v3 bootstrap rows are missing columns: {missing}"
        )
    truth = np.asarray(truth, dtype=int)
    probability = np.asarray(probability, dtype=float)
    if truth.shape != (len(rows),):
        raise ValueError("bootstrap truth must have one value per row")
    if probability.shape != (len(rows), len(MORPHOLOGY_CLASSES)):
        raise ValueError(
            "bootstrap probabilities must have shape "
            f"(n,{len(MORPHOLOGY_CLASSES)})"
        )
    if int(draws) <= 0:
        raise ValueError("bootstrap draws must be positive")
    if not 0.0 < float(confidence_level) < 1.0:
        raise ValueError("bootstrap confidence_level must lie in (0,1)")
    if int(n_calibration_bins) <= 1:
        raise ValueError("n_calibration_bins must exceed one")
    split = rows["fixed_split"].fillna("").astype(str).str.strip().str.lower()
    if set(split) != {"test"}:
        raise ValueError(
            "Teacher v3 bootstrap accepts fixed-test rows only"
        )
    injected = _truth(rows["is_injected_row"]).to_numpy()
    active = (~injected) & (truth >= 0)
    if not np.any(active):
        raise ValueError(
            "Teacher v3 bootstrap has no active real fixed-test rows"
        )
    active_probability = probability[active]
    if (
        not np.isfinite(active_probability).all()
        or (active_probability < 0.0).any()
        or not np.allclose(
            active_probability.sum(axis=1),
            1.0,
            rtol=1.0e-5,
            atol=1.0e-7,
        )
    ):
        raise ValueError(
            "Teacher v3 bootstrap probabilities must be finite, nonnegative, "
            "and row-normalized"
        )

    selected_positions = np.flatnonzero(active)
    evaluation_rows = rows.iloc[selected_positions].copy()
    evaluation_rows["_source_position"] = selected_positions
    evaluation_rows["tic"] = pd.to_numeric(
        evaluation_rows["tic"],
        errors="raise",
    ).astype(np.int64)
    evaluation_rows = evaluation_rows.sort_values(
        ["tic", "review_id"],
        kind="stable",
    ).reset_index(drop=True)
    ordered_positions = evaluation_rows["_source_position"].to_numpy(dtype=int)
    evaluation_truth = truth[ordered_positions]
    evaluation_probability = probability[ordered_positions]
    evaluation_tics = evaluation_rows["tic"].to_numpy(dtype=np.int64)
    unique_tics = np.sort(np.unique(evaluation_tics))
    cluster_indices = [
        np.flatnonzero(evaluation_tics == tic)
        for tic in unique_tics
    ]
    if not len(cluster_indices):
        raise RuntimeError("Teacher v3 bootstrap has no TIC clusters")

    full_class_support = np.bincount(
        evaluation_truth,
        minlength=len(MORPHOLOGY_CLASSES),
    )
    required_classes = tuple(
        int(index)
        for index, count in enumerate(full_class_support)
        if int(count) > 0
    )
    point = _morphology_metric_values(
        evaluation_truth,
        evaluation_probability,
        require_supported_classes=required_classes,
    )

    confidence = evaluation_probability.max(axis=1)
    correct = evaluation_probability.argmax(axis=1) == evaluation_truth
    bin_index = np.minimum(
        np.floor(confidence * int(n_calibration_bins)).astype(int),
        int(n_calibration_bins) - 1,
    )
    edges = np.linspace(0.0, 1.0, int(n_calibration_bins) + 1)
    reliability_accuracy = np.full(
        (int(draws), int(n_calibration_bins)),
        np.nan,
        dtype=float,
    )
    reliability_confidence = np.full_like(
        reliability_accuracy,
        np.nan,
    )
    bootstrap_values = {
        metric: np.full(int(draws), np.nan, dtype=float)
        for metric in TEACHER_V3_BOOTSTRAP_METRICS
    }
    rng = np.random.default_rng(int(seed))
    for draw in range(int(draws)):
        sampled_clusters = rng.integers(
            0,
            len(cluster_indices),
            size=len(cluster_indices),
        )
        sampled_indices = np.concatenate(
            [cluster_indices[int(index)] for index in sampled_clusters]
        )
        draw_metrics = _morphology_metric_values(
            evaluation_truth[sampled_indices],
            evaluation_probability[sampled_indices],
            require_supported_classes=required_classes,
        )
        for metric in TEACHER_V3_BOOTSTRAP_METRICS:
            bootstrap_values[metric][draw] = draw_metrics[metric]
        sampled_bins = bin_index[sampled_indices]
        for calibration_bin in range(int(n_calibration_bins)):
            in_bin = sampled_bins == calibration_bin
            if not np.any(in_bin):
                continue
            selected = sampled_indices[in_bin]
            reliability_accuracy[draw, calibration_bin] = float(
                np.mean(correct[selected])
            )
            reliability_confidence[draw, calibration_bin] = float(
                np.mean(confidence[selected])
            )

    alpha = 0.5 * (1.0 - float(confidence_level))

    def interval(values: np.ndarray) -> tuple[float, float, int]:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if not len(finite):
            return float("nan"), float("nan"), 0
        low, high = np.quantile(
            finite,
            [alpha, 1.0 - alpha],
            method="linear",
        )
        return float(low), float(high), int(len(finite))

    class_metric_index = {
        "planet_recall": 0,
        "eb_recall": 1,
        "variable_recall": 2,
        "other_recall": 3,
    }
    metric_records: list[dict[str, Any]] = []
    for metric in TEACHER_V3_BOOTSTRAP_METRICS:
        class_index = class_metric_index.get(metric)
        if class_index is None:
            support = np.ones(len(evaluation_truth), dtype=bool)
        else:
            support = evaluation_truth == class_index
        low, high, finite_draws = interval(bootstrap_values[metric])
        n_support = int(np.count_nonzero(support))
        n_support_tics = int(np.unique(evaluation_tics[support]).size)
        metric_records.append(
            {
                "profile": str(profile),
                "label_policy": str(label_policy),
                "evaluation_scope": str(evaluation_scope),
                "metric": metric,
                "point_estimate": float(point[metric]),
                "ci_low": low,
                "ci_high": high,
                "confidence_level": float(confidence_level),
                "bootstrap_unit": "tic",
                "draws": int(draws),
                "finite_draws": finite_draws,
                "seed": int(seed),
                "n_rows": int(len(evaluation_truth)),
                "n_tics": int(len(unique_tics)),
                "support_rows": n_support,
                "support_tics": n_support_tics,
                "support_warning": (
                    "small_real_test_support"
                    if n_support < 20
                    else ""
                ),
            }
        )

    reliability_records: list[dict[str, Any]] = []
    for calibration_bin in range(int(n_calibration_bins)):
        in_bin = bin_index == calibration_bin
        accuracy_low, accuracy_high, accuracy_draws = interval(
            reliability_accuracy[:, calibration_bin]
        )
        confidence_low, confidence_high, confidence_draws = interval(
            reliability_confidence[:, calibration_bin]
        )
        reliability_records.append(
            {
                "profile": str(profile),
                "label_policy": str(label_policy),
                "evaluation_scope": str(evaluation_scope),
                "bin_index": calibration_bin,
                "bin_low": float(edges[calibration_bin]),
                "bin_high": float(edges[calibration_bin + 1]),
                "n_rows": int(np.count_nonzero(in_bin)),
                "n_tics": int(np.unique(evaluation_tics[in_bin]).size),
                "mean_confidence": (
                    float(np.mean(confidence[in_bin]))
                    if np.any(in_bin)
                    else float("nan")
                ),
                "accuracy": (
                    float(np.mean(correct[in_bin]))
                    if np.any(in_bin)
                    else float("nan")
                ),
                "accuracy_ci_low": accuracy_low,
                "accuracy_ci_high": accuracy_high,
                "confidence_ci_low": confidence_low,
                "confidence_ci_high": confidence_high,
                "confidence_level": float(confidence_level),
                "bootstrap_unit": "tic",
                "draws": int(draws),
                "accuracy_finite_draws": accuracy_draws,
                "confidence_finite_draws": confidence_draws,
                "seed": int(seed),
                "support_warning": (
                    "empty_bin"
                    if not np.any(in_bin)
                    else (
                        "small_real_test_support"
                        if int(np.count_nonzero(in_bin)) < 20
                        else ""
                    )
                ),
            }
        )

    metrics = pd.DataFrame(metric_records)
    reliability = pd.DataFrame(reliability_records)
    summary = {
        "schema_version": TEACHER_V3_BOOTSTRAP_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "profile": str(profile),
        "label_policy": str(label_policy),
        "evaluation_scope": str(evaluation_scope),
        "bootstrap_unit": "tic",
        "draws": int(draws),
        "seed": int(seed),
        "confidence_level": float(confidence_level),
        "n_calibration_bins": int(n_calibration_bins),
        "n_rows": int(len(evaluation_truth)),
        "n_tics": int(len(unique_tics)),
        "class_support_rows": {
            label: int(full_class_support[index])
            for index, label in enumerate(MORPHOLOGY_CLASSES)
        },
        "class_support_tics": {
            label: int(
                np.unique(
                    evaluation_tics[evaluation_truth == index]
                ).size
            )
            for index, label in enumerate(MORPHOLOGY_CLASSES)
        },
        "interpretation": (
            "Percentile intervals resample complete TIC clusters from the "
            "real fixed-test partition. Class metrics with fewer than 20 "
            "support rows are explicitly marked as small-support estimates."
        ),
    }
    return metrics, reliability, summary


def write_teacher_v3_fixed_test_bootstrap(
    *,
    rows: pd.DataFrame,
    truth: np.ndarray,
    probability: np.ndarray,
    profile: str,
    label_policy: str,
    out_dir: Path,
    predictions_path: Path,
    draws: int = 2000,
    seed: int = 560063,
    confidence_level: float = 0.95,
    filename_prefix: str = "fixed_test",
    evaluation_scope: str = "fixed_test_real_labels",
) -> dict[str, Any]:
    """Write deterministic interval and reliability products for one profile."""

    out_dir = Path(out_dir)
    predictions_path = Path(predictions_path)
    metrics, reliability, summary = teacher_v3_tic_cluster_bootstrap(
        rows=rows,
        truth=truth,
        probability=probability,
        profile=profile,
        label_policy=label_policy,
        draws=draws,
        seed=seed,
        confidence_level=confidence_level,
        evaluation_scope=evaluation_scope,
    )
    prefix = str(filename_prefix).strip()
    if not prefix or "/" in prefix or "\\" in prefix:
        raise ValueError("bootstrap filename_prefix must be one safe basename")
    metrics_path = out_dir / f"{prefix}_bootstrap_metrics.csv"
    reliability_path = out_dir / f"{prefix}_reliability.csv"
    metrics.to_csv(
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
    summary.update(
        {
            "fixed_test_predictions": str(predictions_path),
            "fixed_test_predictions_sha256": _file_sha256(predictions_path),
            "metrics_path": str(metrics_path),
            "metrics_sha256": _file_sha256(metrics_path),
            "reliability_path": str(reliability_path),
            "reliability_sha256": _file_sha256(reliability_path),
        }
    )
    summary_path = out_dir / f"{prefix}_bootstrap.summary.json"
    summary_sha256 = _write_json(summary_path, summary)
    return {
        **summary,
        "summary_path": str(summary_path),
        "summary_sha256": summary_sha256,
    }


def validate_teacher_v3_training_table(source: pd.DataFrame) -> None:
    """Validate the frozen release identity before touching native inputs."""

    required = {
        "review_id",
        "teacher_v3_observation_id",
        "teacher_run_name",
        "teacher_v3_corpus_version",
        "teacher_architecture_version",
        "teacher_v3_training_include",
        "sector",
        "tic",
        "fixed_split",
        "cv_fold",
        "native_h5_path",
        "native_group_path",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise KeyError(f"Teacher v3 training table is missing columns: {missing}")
    if source.empty:
        raise ValueError("Teacher v3 training table is empty")
    for name in ("review_id", "teacher_v3_observation_id"):
        values = source[name].fillna("").astype(str).str.strip()
        if values.eq("").any() or values.duplicated().any():
            raise ValueError(f"{name} must be nonblank and globally unique")
    identities = {
        "teacher_run_name": TEACHER_V3_RUN_NAME,
        "teacher_v3_corpus_version": TEACHER_V3_CORPUS_VERSION,
        "teacher_architecture_version": MODEL_VERSION,
    }
    for name, expected in identities.items():
        observed = set(source[name].fillna("").astype(str).str.strip())
        if observed != {expected}:
            raise ValueError(
                f"{name} must be exactly {expected!r}; observed={sorted(observed)}"
            )
    sectors = {
        int(value)
        for value in pd.to_numeric(source["sector"], errors="raise")
    }
    if sectors != set(TEACHER_V3_SECTORS):
        raise ValueError(
            f"Teacher v3 sectors must be {list(TEACHER_V3_SECTORS)}; "
            f"observed={sorted(sectors)}"
        )
    active = _truth(source["teacher_v3_training_include"])
    if not active.any():
        raise ValueError("Teacher v3 has no active training rows")
    for name in ("native_h5_path", "native_group_path"):
        blank = source.loc[active, name].fillna("").astype(str).str.strip().eq("")
        if blank.any():
            raise ValueError(
                f"{name} is blank for {int(blank.sum())} active Teacher v3 rows"
            )
    split = source["fixed_split"].fillna("").astype(str).str.strip().str.lower()
    if set(split) != {"development", "test"}:
        raise ValueError(
            "fixed_split must contain exactly development and test"
        )
    folds = pd.to_numeric(source["cv_fold"], errors="raise").astype(int)
    if set(folds[split.eq("development")]) != set(range(5)):
        raise ValueError("development rows must cover all five folds 0..4")
    if not folds[split.eq("test")].eq(-1).all():
        raise ValueError("fixed test rows must have cv_fold=-1")


def _native_mapping_rows(rows: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "sector",
        "tic",
        "native_h5_path",
        "native_group_path",
        "is_injected_row",
    ]
    missing = sorted(set(columns) - set(rows.columns))
    if missing:
        raise KeyError(f"prepared Teacher v3 rows are missing columns: {missing}")
    mappings = rows.loc[:, columns].copy()
    mappings["sector"] = pd.to_numeric(
        mappings["sector"], errors="raise"
    ).astype(np.int64)
    mappings["tic"] = pd.to_numeric(
        mappings["tic"], errors="raise"
    ).astype(np.int64)
    mappings["native_h5_path"] = (
        mappings["native_h5_path"]
        .fillna("")
        .astype(str)
        .str.strip()
        .map(lambda value: str(Path(value).expanduser().resolve()))
    )
    mappings["native_group_path"] = (
        mappings["native_group_path"].fillna("").astype(str).str.strip()
    )
    mappings["is_injected_row"] = _truth(mappings["is_injected_row"])
    if mappings["native_h5_path"].eq("").any():
        raise ValueError("prepared Teacher v3 rows contain blank native paths")
    if mappings["native_group_path"].eq("").any():
        raise ValueError("prepared Teacher v3 rows contain blank native groups")
    identity = mappings.loc[
        :,
        ["native_h5_path", "native_group_path", "sector", "tic"],
    ].drop_duplicates()
    collisions = identity.duplicated(
        ["native_h5_path", "native_group_path"],
        keep=False,
    )
    if collisions.any():
        examples = identity.loc[collisions].head(5).to_dict("records")
        raise ValueError(
            "one native HDF5 group maps to more than one (sector,TIC); "
            f"first={examples}"
        )
    return mappings.drop_duplicates().reset_index(drop=True)


def _validate_real_registry(
    *,
    mappings: pd.DataFrame,
    registry_path: Path,
    registry_summary_path: Path,
) -> dict[str, Any]:
    audit = validate_native_input_registry_path(
        registry_path=registry_path,
        summary_path=registry_summary_path,
        expected_contract_version=RAW_PAIR_CONTRACT_VERSION,
    )
    registry = validate_native_input_registry(
        read_table(registry_path),
        path_base=Path(registry_path).parent,
        expected_contract_version=RAW_PAIR_CONTRACT_VERSION,
    )
    expected = (
        mappings.loc[
            ~mappings["is_injected_row"],
            ["sector", "tic", "native_h5_path", "native_group_path"],
        ]
        .drop_duplicates()
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    observed = (
        registry.loc[
            :,
            ["sector", "tic", "native_h5_path", "native_group_path"],
        ]
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    if not expected.equals(observed):
        comparison = expected.merge(
            observed,
            on=["sector", "tic", "native_h5_path", "native_group_path"],
            how="outer",
            indicator=True,
        )
        examples = comparison.loc[comparison["_merge"].ne("both")].head(5)
        raise ValueError(
            "real Teacher v3 rows do not exactly match the frozen native "
            f"registry; first={examples.to_dict('records')}"
        )
    return audit


def build_teacher_v3_native_manifest(
    *,
    rows: pd.DataFrame,
    registry_path: Path,
    registry_summary_path: Path,
) -> dict[str, Any]:
    """Fully validate and hash all row-addressed Teacher v3 HDF5 inputs."""

    import h5py

    mappings = _native_mapping_rows(rows)
    registry_path = Path(registry_path).resolve()
    registry_summary_path = Path(registry_summary_path).resolve()
    registry_audit = _validate_real_registry(
        mappings=mappings,
        registry_path=registry_path,
        registry_summary_path=registry_summary_path,
    )
    records: list[dict[str, Any]] = []
    for path_text, expected_rows in mappings.groupby(
        "native_h5_path",
        sort=True,
    ):
        path = Path(path_text)
        if not path.is_file():
            raise FileNotFoundError(f"missing Teacher v3 native HDF5: {path}")
        stat_before = path.stat()
        verification = verify_raw_pair_contract(
            path,
            require_errors=True,
            require_periodograms=True,
        )
        if not verification["passed"]:
            raise RuntimeError(
                f"Teacher v3 native contract failed for {path}: "
                f"{verification['failures'][:10]}"
            )
        expected_groups = set(expected_rows["native_group_path"])
        with h5py.File(path, "r") as h5:
            contract = str(h5.attrs.get("contract_version", ""))
            observed_groups = {
                f"{root}/{key}"
                for root in ("targets", "injections")
                if root in h5
                for key in h5[root]
            }
            if observed_groups != expected_groups:
                missing = sorted(expected_groups - observed_groups)[:10]
                extra = sorted(observed_groups - expected_groups)[:10]
                raise ValueError(
                    f"native group set mismatch in {path}: "
                    f"missing={missing}, extra={extra}"
                )
            for row in (
                expected_rows.loc[
                    :,
                    ["native_group_path", "sector", "tic"],
                ]
                .drop_duplicates()
                .to_dict("records")
            ):
                group = h5[str(row["native_group_path"])]
                try:
                    observed_sector = int(group.attrs["sector"])
                    observed_tic = int(group.attrs["tic"])
                except (KeyError, TypeError, ValueError) as exc:
                    raise ValueError(
                        f"{path}:{row['native_group_path']} lacks scalar "
                        "sector/TIC identity attributes"
                    ) from exc
                if (
                    observed_sector,
                    observed_tic,
                ) != (int(row["sector"]), int(row["tic"])):
                    raise ValueError(
                        f"{path}:{row['native_group_path']} declares "
                        f"(sector,TIC)=({observed_sector},{observed_tic}), "
                        f"expected=({row['sector']},{row['tic']})"
                    )
            source_table_sha256 = str(
                h5.attrs.get("training_table_sha256", "")
            )
        digest = file_sha256(path)
        stat_after = path.stat()
        if (
            stat_before.st_size != stat_after.st_size
            or stat_before.st_mtime_ns != stat_after.st_mtime_ns
        ):
            raise RuntimeError(
                f"native HDF5 changed while it was validated: {path}"
            )
        records.append(
            {
                "path": str(path),
                "sha256": digest,
                "size_bytes": int(stat_after.st_size),
                "mtime_ns": int(stat_after.st_mtime_ns),
                "contract_version": contract,
                "source_training_table_sha256": source_table_sha256,
                "n_training_rows": int(len(expected_rows)),
                "n_unique_groups": int(len(expected_groups)),
                "group_counts": verification["counts"],
                "external_quality_counts": verification[
                    "external_quality_counts"
                ],
            }
        )
    observed_sectors = set(
        pd.to_numeric(mappings["sector"], errors="raise").astype(int)
    )
    if observed_sectors != set(TEACHER_V3_SECTORS):
        raise ValueError(
            "native manifest does not cover exactly sectors 56--62"
        )
    if len(records) != len(TEACHER_V3_SECTORS):
        raise ValueError(
            "Teacher v3 requires exactly one native HDF5 per sector; "
            f"found {len(records)}"
        )
    manifest: dict[str, Any] = {
        "schema_version": TEACHER_V3_NATIVE_MANIFEST_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "model_version": MODEL_VERSION,
        "native_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "sectors": list(TEACHER_V3_SECTORS),
        "n_training_rows": int(len(rows)),
        "n_unique_native_groups": int(
            mappings[
                ["native_h5_path", "native_group_path"]
            ].drop_duplicates().shape[0]
        ),
        "registry": {
            "path": str(registry_path),
            "sha256": _file_sha256(registry_path),
            "summary_path": str(registry_summary_path),
            "summary_sha256": _file_sha256(registry_summary_path),
            "audit": registry_audit,
        },
        "native_files": records,
    }
    return manifest


def _manifest_digest(manifest: Mapping[str, Any]) -> str:
    return hashlib.sha256(_json_payload(manifest)).hexdigest()


def build_teacher_v3_input_provenance(
    *,
    training_table: Path,
    native_manifest: Mapping[str, Any],
) -> dict[str, str]:
    """Bind checkpoints to the table and canonical seven-file manifest."""

    manifest_sha256 = _manifest_digest(native_manifest)
    return {
        "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        # Kept for compatibility with the v1 checkpoint validator.  For v3
        # this is the digest of the canonical multi-file manifest, not one H5.
        "native_h5_sha256": manifest_sha256,
        "native_manifest_sha256": manifest_sha256,
        "training_table_sha256": _file_sha256(training_table),
        "native_registry_sha256": str(
            native_manifest["registry"]["sha256"]
        ),
        "native_registry_summary_sha256": str(
            native_manifest["registry"]["summary_sha256"]
        ),
    }


def validate_teacher_v3_input_provenance(
    payload: Mapping[str, Any],
    *,
    expected: Mapping[str, str],
    artifact: str,
) -> None:
    """Validate the legacy core plus every Teacher-v3-specific binding."""

    validate_teacher_input_provenance(
        payload,
        expected=expected,
        artifact=artifact,
    )
    extended_fields = (
        "native_manifest_sha256",
        "native_registry_sha256",
        "native_registry_summary_sha256",
        "native_input_mode",
        "label_policy",
        "masked_review_ids_sha256",
    )
    failures = [
        f"{name}={payload.get(name)!r}, expected {expected.get(name)!r}"
        for name in extended_fields
        if name in expected
        and str(payload.get(name, "")) != str(expected.get(name, ""))
    ]
    if failures:
        raise RuntimeError(
            f"stale or unprovenanced {artifact}: " + "; ".join(failures)
        )


def _oof_frame_hash(frame: pd.DataFrame) -> str:
    columns = [
        "review_id",
        "morphology_target",
        *[f"logit_{label}" for label in MORPHOLOGY_CLASSES],
    ]
    ordered = frame.loc[:, columns].sort_values(
        "review_id",
        kind="stable",
    )
    payload = ordered.to_csv(
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def calibrate_teacher_v3_profile_oof(
    *,
    rows: pd.DataFrame,
    out_dir: Path,
    profile: str,
    input_provenance: Mapping[str, str],
) -> dict[str, Any]:
    """Fit and apply one temperature to concatenated five-fold OOF logits."""

    import torch

    out_dir = Path(out_dir)
    parts: list[pd.DataFrame] = []
    for fold in range(5):
        path = out_dir / profile / f"fold_{fold}" / "validation_predictions.csv"
        part = pd.read_csv(path)
        part["cv_fold"] = fold
        parts.append(part)
    predictions = pd.concat(parts, ignore_index=True)
    if predictions["review_id"].duplicated().any():
        raise RuntimeError(
            f"duplicate development OOF predictions for profile {profile}"
        )
    expected_ids = set(
        rows.loc[rows["fixed_split"].eq("development"), "review_id"].astype(str)
    )
    observed_ids = set(predictions["review_id"].astype(str))
    if observed_ids != expected_ids:
        raise RuntimeError(
            f"{profile} OOF predictions do not exactly cover development rows"
        )
    logit_columns = [f"logit_{label}" for label in MORPHOLOGY_CLASSES]
    missing = sorted(set(logit_columns) - set(predictions.columns))
    if missing:
        raise RuntimeError(
            f"{profile} OOF predictions lack raw logits: {missing}"
        )
    logits = predictions.loc[:, logit_columns].to_numpy(dtype=float)
    truth = predictions["morphology_target"].to_numpy(dtype=int)
    temperature = fit_temperature(logits, truth)
    probability = _softmax(logits, temperature=temperature)
    predictions["morphology_prediction"] = probability.argmax(axis=1)
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        predictions[f"p_{label}"] = probability[:, index]
    oof_hash = _oof_frame_hash(predictions)
    active = truth >= 0
    if not np.any(active):
        raise RuntimeError(f"{profile} OOF predictions have no morphology targets")
    uncalibrated = _softmax(logits)
    calibration = {
        "schema_version": TEACHER_V3_CALIBRATION_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "profile": profile,
        "scope": "concatenated_five_fold_development_oof_logits",
        "n_oof_rows": int(len(predictions)),
        "oof_logits_sha256": oof_hash,
        "temperature": float(temperature),
        "uncalibrated_cross_entropy": float(
            -np.mean(
                np.log(
                    np.clip(
                        uncalibrated[
                            np.flatnonzero(active),
                            truth[active],
                        ],
                        1.0e-12,
                        1.0,
                    )
                )
            )
        ),
        "calibrated_cross_entropy": float(
            -np.mean(
                np.log(
                    np.clip(
                        probability[
                            np.flatnonzero(active),
                            truth[active],
                        ],
                        1.0e-12,
                        1.0,
                    )
                )
            )
        ),
        **dict(input_provenance),
    }
    calibration_path = out_dir / profile / "pooled_oof_calibration.json"
    calibration_sha256 = _write_json(calibration_path, calibration)

    row_lookup = rows.set_index("review_id", drop=False)
    for fold in range(5):
        fold_mask = predictions["cv_fold"].eq(fold)
        fold_predictions = predictions.loc[fold_mask].drop(
            columns=["cv_fold"]
        )
        fold_path = (
            out_dir
            / profile
            / f"fold_{fold}"
            / "validation_predictions.csv"
        )
        fold_predictions.to_csv(fold_path, index=False)
        checkpoint_path = (
            out_dir / profile / f"fold_{fold}" / "teacher.pt"
        )
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        validate_teacher_v3_input_provenance(
            checkpoint,
            expected=input_provenance,
            artifact=f"Teacher v3 {profile} fold-{fold} checkpoint",
        )
        if checkpoint.get("temperature_calibration_scope") != (
            "pending_pooled_oof_development"
        ):
            raise RuntimeError(
                f"{checkpoint_path} was not held for pooled OOF calibration"
            )
        checkpoint["temperature"] = float(temperature)
        checkpoint["temperature_calibration_scope"] = (
            "pooled_oof_development"
        )
        checkpoint["pooled_oof_calibration_sha256"] = calibration_sha256
        checkpoint["run_id"] = TEACHER_V3_RUN_ID
        checkpoint["release_name"] = TEACHER_V3_RUN_NAME
        temporary = checkpoint_path.with_suffix(".pt.tmp")
        temporary.unlink(missing_ok=True)
        torch.save(checkpoint, temporary)
        temporary.replace(checkpoint_path)

        fold_truth = fold_predictions[
            "morphology_target"
        ].to_numpy(dtype=int)
        fold_probability = fold_predictions.loc[
            :,
            [f"p_{label}" for label in MORPHOLOGY_CLASSES],
        ].to_numpy(dtype=float)
        fold_rows = row_lookup.loc[
            fold_predictions["review_id"].astype(str)
        ].reset_index(drop=True)
        metrics_path = (
            out_dir / profile / f"fold_{fold}" / "metrics.json"
        )
        metrics = json.loads(metrics_path.read_text())
        metrics["morphology"] = classification_metrics(
            fold_truth,
            fold_probability,
            classes=MORPHOLOGY_CLASSES,
        )
        metrics["morphology_by_source"] = _subset_metrics(
            fold_rows,
            fold_truth,
            fold_probability,
        )
        metrics["calibration"] = expected_calibration_error(
            fold_truth,
            fold_probability,
        )
        metrics["temperature"] = float(temperature)
        metrics["temperature_calibration_scope"] = (
            "pooled_oof_development"
        )
        metrics["pooled_oof_calibration_sha256"] = calibration_sha256
        _write_json(metrics_path, metrics)
    development_rows = row_lookup.loc[
        predictions["review_id"].astype(str)
    ].reset_index(drop=True)
    source_metrics = _subset_metrics(
        development_rows,
        truth,
        probability,
    )
    return {
        "profile": profile,
        "temperature": float(temperature),
        "calibration_path": str(calibration_path),
        "calibration_sha256": calibration_sha256,
        "oof_logits_sha256": oof_hash,
        "metrics": {
            "morphology_by_source": source_metrics,
            "calibration": _calibration_by_source(
                development_rows,
                truth,
                probability,
            ),
        },
    }


def _evaluate_fixed_profile(
    *,
    rows: pd.DataFrame,
    out_dir: Path,
    profile: str,
    train_config: HarmonicTrainConfig,
    workers: int,
    input_provenance: Mapping[str, str],
    require_cuda: bool,
    label_policy: str = TEACHER_V3_PRIMARY_LABEL_POLICY,
    bootstrap_draws: int = 2000,
    bootstrap_seed: int = 560063,
    bootstrap_confidence_level: float = 0.95,
) -> dict[str, Any]:
    import torch

    test_mask = rows["fixed_split"].eq("test").to_numpy()
    test_rows = rows.loc[test_mask].reset_index(drop=True)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required but is unavailable")
    morphology_members: list[np.ndarray] = []
    preserve_members: list[np.ndarray] = []
    harmonic_members: list[np.ndarray] = []
    truth: np.ndarray | None = None
    preserve_truth: np.ndarray | None = None
    harmonic_truth: np.ndarray | None = None
    review_ids: list[str] | None = None
    temperatures: set[float] = set()
    for fold in range(5):
        checkpoint_path = (
            Path(out_dir) / profile / f"fold_{fold}" / "teacher.pt"
        )
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        validate_teacher_v3_input_provenance(
            checkpoint,
            expected=input_provenance,
            artifact=f"Teacher v3 {profile} fold-{fold} checkpoint",
        )
        if checkpoint.get("temperature_calibration_scope") != (
            "pooled_oof_development"
        ):
            raise RuntimeError(
                f"{checkpoint_path} lacks pooled OOF calibration"
            )
        temperatures.add(float(checkpoint["temperature"]))
        norm = checkpoint["metadata_normalization"]
        normalization = MetadataNormalization(
            columns=tuple(norm["columns"]),
            center=tuple(norm["center"]),
            scale=tuple(norm["scale"]),
        )
        metadata, _ = build_metadata_matrix(
            rows,
            fit_mask=np.zeros(len(rows), dtype=bool),
            normalization=normalization,
        )
        dataset = HarmonicNativeDataset(
            rows,
            native_h5=None,
            metadata=metadata,
            profile=profile,
        )
        loader = _loader(
            dataset,
            np.flatnonzero(test_mask),
            batch_size=train_config.batch_size,
            shuffle=False,
            workers=workers,
            seed=train_config.seed,
        )
        model = build_harmonic_cnn(
            HarmonicModelConfig(**checkpoint["model_config"]),
            profile=profile,
        ).to(device)
        model.load_state_dict(checkpoint["model_state_dict"])
        scored = _evaluate_loader(model, loader, device=device)
        morphology_members.append(
            _softmax(
                scored["morphology_logits"],
                temperature=float(checkpoint["temperature"]),
            )
        )
        preserve_members.append(_softmax(scored["preserve_logits"]))
        harmonic_members.append(_softmax(scored["harmonic_logits"]))
        if truth is None:
            truth = scored["morphology_target"]
            preserve_truth = scored["preserve_target"]
            harmonic_truth = scored["harmonic_target"]
            review_ids = scored["review_id"]
        elif review_ids != scored["review_id"]:
            raise RuntimeError(
                f"{profile} test ordering changed between ensemble members"
            )
    if len(temperatures) != 1:
        raise RuntimeError(
            f"{profile} fold checkpoints disagree on pooled temperature"
        )
    assert truth is not None
    assert preserve_truth is not None
    assert harmonic_truth is not None
    assert review_ids is not None
    morphology_probability = np.mean(morphology_members, axis=0)
    preserve_probability = np.mean(preserve_members, axis=0)
    harmonic_probability = np.mean(harmonic_members, axis=0)
    metrics: dict[str, Any] = {
        "morphology_by_source": _subset_metrics(
            test_rows,
            truth,
            morphology_probability,
        ),
        "calibration": _calibration_by_source(
            test_rows,
            truth,
            morphology_probability,
        ),
        "preserve": classification_metrics(
            preserve_truth,
            preserve_probability,
            classes=PRESERVE_CLASSES,
        ),
        "harmonic": classification_metrics(
            harmonic_truth,
            harmonic_probability,
            classes=HARMONIC_CLASSES,
        ),
        "temperature": float(next(iter(temperatures))),
        "temperature_calibration_scope": "pooled_oof_development",
        "label_policy": str(label_policy),
        "fixed_test_role": "opened_after_pretest_model_freeze",
    }
    predictions = test_rows.loc[
        :,
        [
            "review_id",
            "sector",
            "tic",
            "human_label",
            "is_injected_row",
            "fixed_split",
        ],
    ].copy()
    predictions["morphology_target_index"] = truth
    predictions["morphology_prediction_index"] = (
        morphology_probability.argmax(axis=1)
    )
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        predictions[f"p_{label}"] = morphology_probability[:, index]
    predictions["p_preserve"] = preserve_probability[:, 1]
    predictions["predicted_period_factor"] = np.asarray(
        [0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0]
    )[harmonic_probability.argmax(axis=1)]
    predictions["label_policy"] = str(label_policy)
    predictions_path = (
        Path(out_dir) / profile / "fixed_test_predictions.csv"
    )
    predictions.to_csv(
        predictions_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    metrics["tic_cluster_bootstrap"] = (
        write_teacher_v3_fixed_test_bootstrap(
            rows=test_rows,
            truth=truth,
            probability=morphology_probability,
            profile=profile,
            label_policy=label_policy,
            out_dir=Path(out_dir) / profile,
            predictions_path=predictions_path,
            draws=bootstrap_draws,
            seed=bootstrap_seed,
            confidence_level=bootstrap_confidence_level,
            evaluation_scope=(
                "fixed_test_real_non_uncertain"
                if label_policy == TEACHER_V3_UNCERTAIN_MASKED_POLICY
                else "fixed_test_real_labels"
            ),
        )
    )
    _write_json(
        Path(out_dir) / profile / "fixed_test_metrics.json",
        metrics,
    )
    return metrics


def write_teacher_v3_non_uncertain_comparison_bootstrap(
    *,
    predictions_path: Path,
    profile: str,
    out_dir: Path,
    draws: int = 2000,
    seed: int = 560063,
    confidence_level: float = 0.95,
) -> dict[str, Any]:
    """Evaluate a canonical model on the sensitivity model's test support."""

    predictions_path = Path(predictions_path)
    predictions = pd.read_csv(predictions_path, low_memory=False)
    required = {
        "review_id",
        "tic",
        "human_label",
        "is_injected_row",
        "fixed_split",
        "morphology_target_index",
        *{f"p_{label}" for label in MORPHOLOGY_CLASSES},
    }
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise KeyError(
            "Teacher v3 fixed-test predictions lack common-support columns: "
            f"{missing}"
        )
    selected = predictions.loc[
        ~predictions["human_label"]
        .fillna("")
        .astype(str)
        .str.strip()
        .eq("uncertain")
    ].copy()
    if selected.empty or len(selected) == len(predictions):
        raise ValueError(
            "common-support evaluation must remove at least one uncertain row"
        )
    truth = pd.to_numeric(
        selected["morphology_target_index"],
        errors="raise",
    ).to_numpy(dtype=int)
    probability = selected.loc[
        :,
        [f"p_{label}" for label in MORPHOLOGY_CLASSES],
    ].to_numpy(dtype=float)
    return write_teacher_v3_fixed_test_bootstrap(
        rows=selected.reset_index(drop=True),
        truth=truth,
        probability=probability,
        profile=profile,
        label_policy=TEACHER_V3_PRIMARY_LABEL_POLICY,
        out_dir=Path(out_dir),
        predictions_path=predictions_path,
        draws=draws,
        seed=seed,
        confidence_level=confidence_level,
        filename_prefix="fixed_test_non_uncertain",
        evaluation_scope="fixed_test_real_non_uncertain",
    )


def _profile_development_record(
    result: Mapping[str, Any],
) -> dict[str, Any]:
    metrics = result["metrics"]["morphology_by_source"]
    all_metrics = metrics["all"]
    real_metrics = metrics["real"]
    per_class = real_metrics["per_class"]
    return {
        "profile": str(result["profile"]),
        "fixed_role": (
            "primary"
            if result["profile"] == TEACHER_V3_PRIMARY_PROFILE
            else "metadata_baseline"
        ),
        "temperature": float(result["temperature"]),
        "validation_macro_f1": float(all_metrics["macro_f1"]),
        "validation_balanced_accuracy": float(
            all_metrics["balanced_accuracy"]
        ),
        "validation_real_planet_recall": float(
            per_class.get("planet_like", {}).get("recall", float("nan"))
        ),
        "validation_real_eb_recall": float(
            per_class.get("eclipse_contact", {}).get(
                "recall",
                float("nan"),
            )
        ),
        "validation_real_variable_recall": float(
            per_class.get("smooth_variable", {}).get(
                "recall",
                float("nan"),
            )
        ),
        "validation_real_other_recall": float(
            per_class.get("other", {}).get("recall", float("nan"))
        ),
        "validation_ece": float(
            result["metrics"]["calibration"]["all"]["ece"]
        ),
    }


def _write_selected_checkpoint_manifest(
    *,
    out_dir: Path,
    input_provenance: Mapping[str, str],
    selected_profile: str = TEACHER_V3_PRIMARY_PROFILE,
    filename: str = "selected_checkpoint_manifest.json",
    selection_policy: str = "fixed_before_test",
) -> Path:
    import torch

    records: list[dict[str, Any]] = []
    for fold in range(5):
        path = (
            Path(out_dir)
            / selected_profile
            / f"fold_{fold}"
            / "teacher.pt"
        )
        checkpoint = torch.load(
            path,
            map_location="cpu",
            weights_only=False,
        )
        validate_teacher_v3_input_provenance(
            checkpoint,
            expected=input_provenance,
            artifact=f"Teacher v3 selected fold-{fold} checkpoint",
        )
        if checkpoint.get("temperature_calibration_scope") != (
            "pooled_oof_development"
        ):
            raise RuntimeError(
                f"selected Teacher v3 checkpoint lacks pooled calibration: {path}"
            )
        records.append(
            {
                "fold": fold,
                "path": str(path.relative_to(out_dir)),
                "sha256": _file_sha256(path),
                "pooled_oof_calibration_sha256": checkpoint[
                    "pooled_oof_calibration_sha256"
                ],
            }
        )
    payload = {
        "schema_version": TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "selected_profile": selected_profile,
        "selection_policy": selection_policy,
        **dict(input_provenance),
        "checkpoints": records,
    }
    path = Path(out_dir) / str(filename)
    _write_json(path, payload)
    return path


def verify_teacher_v3_checkpoint_manifest(
    *,
    manifest_path: Path,
    expected_manifest_sha256: str,
    expected_profile: str,
    expected_selection_policy: str,
    expected_provenance: Mapping[str, str],
    calibration_path: Path,
    expected_calibration_sha256: str,
) -> dict[str, Any]:
    """Verify one frozen five-fold ensemble without touching fixed-test rows."""

    import torch

    manifest_path = Path(manifest_path).resolve()
    calibration_path = Path(calibration_path).resolve()
    manifest_sha256 = _file_sha256(manifest_path)
    if manifest_sha256 != str(expected_manifest_sha256):
        raise RuntimeError(
            "Teacher v3 checkpoint-manifest hash changed after freeze: "
            f"{manifest_path}"
        )
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict):
        raise RuntimeError(
            f"Teacher v3 checkpoint manifest must be a JSON mapping: "
            f"{manifest_path}"
        )
    expected_manifest = {
        "schema_version": TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "selected_profile": str(expected_profile),
        "selection_policy": str(expected_selection_policy),
    }
    failures = [
        f"{name}={manifest.get(name)!r}, expected {expected!r}"
        for name, expected in expected_manifest.items()
        if manifest.get(name) != expected
    ]
    try:
        validate_teacher_v3_input_provenance(
            manifest,
            expected=expected_provenance,
            artifact=f"Teacher v3 {expected_profile} checkpoint manifest",
        )
    except RuntimeError as exc:
        failures.append(str(exc))

    calibration_sha256 = _file_sha256(calibration_path)
    if calibration_sha256 != str(expected_calibration_sha256):
        failures.append(
            "pooled calibration hash changed after freeze: "
            f"{calibration_path}"
        )
    calibration = json.loads(calibration_path.read_text())
    if not isinstance(calibration, dict):
        failures.append("pooled calibration must be a JSON mapping")
        calibration = {}
    expected_calibration = {
        "schema_version": TEACHER_V3_CALIBRATION_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "profile": str(expected_profile),
        "scope": "concatenated_five_fold_development_oof_logits",
    }
    for name, expected in expected_calibration.items():
        if calibration.get(name) != expected:
            failures.append(
                f"pooled calibration {name}={calibration.get(name)!r}, "
                f"expected {expected!r}"
            )
    for name, expected in expected_provenance.items():
        if str(calibration.get(name, "")) != str(expected):
            failures.append(
                f"pooled calibration {name}={calibration.get(name)!r}, "
                f"expected {expected!r}"
            )
    try:
        calibration_temperature = float(calibration["temperature"])
    except (KeyError, TypeError, ValueError):
        calibration_temperature = float("nan")
        failures.append("pooled calibration temperature is invalid")
    if (
        not np.isfinite(calibration_temperature)
        or calibration_temperature <= 0
    ):
        failures.append(
            "pooled calibration temperature must be finite and positive"
        )

    records = manifest.get("checkpoints")
    if not isinstance(records, list) or len(records) != 5:
        failures.append("checkpoint manifest must contain exactly five folds")
        records = []
    observed_folds: list[int] = []
    checkpoint_hashes: dict[str, str] = {}
    for record in records:
        if not isinstance(record, Mapping):
            failures.append("checkpoint manifest contains a malformed record")
            continue
        try:
            fold = int(record.get("fold", -1))
        except (TypeError, ValueError):
            failures.append("checkpoint manifest contains a non-integer fold")
            continue
        observed_folds.append(fold)
        expected_relative = (
            Path(str(expected_profile)) / f"fold_{fold}" / "teacher.pt"
        ).as_posix()
        if str(record.get("path", "")) != expected_relative:
            failures.append(
                f"fold {fold} checkpoint path is not exactly {expected_relative}"
            )
            continue
        checkpoint_path = manifest_path.parent / expected_relative
        if not checkpoint_path.is_file():
            failures.append(f"fold {fold} checkpoint is missing")
            continue
        checkpoint_sha256 = _file_sha256(checkpoint_path)
        if checkpoint_sha256 != str(record.get("sha256", "")):
            failures.append(f"fold {fold} checkpoint SHA-256 mismatch")
            continue
        if str(
            record.get("pooled_oof_calibration_sha256", "")
        ) != calibration_sha256:
            failures.append(f"fold {fold} manifest calibration hash mismatch")
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        if _file_sha256(checkpoint_path) != checkpoint_sha256:
            failures.append(f"fold {fold} checkpoint changed while it was read")
        try:
            validate_teacher_v3_input_provenance(
                checkpoint,
                expected=expected_provenance,
                artifact=f"Teacher v3 {expected_profile} fold-{fold}",
            )
        except RuntimeError as exc:
            failures.append(str(exc))
        expected_checkpoint = {
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
            "profile": str(expected_profile),
            "fold": fold,
            "temperature_calibration_scope": "pooled_oof_development",
            "pooled_oof_calibration_sha256": calibration_sha256,
        }
        for name, expected in expected_checkpoint.items():
            if checkpoint.get(name) != expected:
                failures.append(
                    f"fold {fold} checkpoint {name}={checkpoint.get(name)!r}, "
                    f"expected {expected!r}"
                )
        try:
            temperature_matches = bool(
                np.isclose(
                    float(checkpoint.get("temperature", float("nan"))),
                    calibration_temperature,
                )
            )
        except (TypeError, ValueError):
            temperature_matches = False
        if not temperature_matches:
            failures.append(
                f"fold {fold} checkpoint temperature disagrees with calibration"
            )
        checkpoint_hashes[str(fold)] = checkpoint_sha256
    if sorted(observed_folds) != list(range(5)):
        failures.append(
            f"checkpoint folds must be 0..4; observed={observed_folds}"
        )
    if _file_sha256(manifest_path) != manifest_sha256:
        failures.append("checkpoint manifest changed while it was verified")
    if _file_sha256(calibration_path) != calibration_sha256:
        failures.append("pooled calibration changed while it was verified")
    if failures:
        raise RuntimeError(
            f"Teacher v3 {expected_profile} frozen-ensemble verification "
            "failed: "
            + "; ".join(failures)
        )
    return {
        "manifest_path": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "profile": str(expected_profile),
        "selection_policy": str(expected_selection_policy),
        "calibration_path": str(calibration_path),
        "calibration_sha256": calibration_sha256,
        "checkpoint_sha256_by_fold": checkpoint_hashes,
        "n_verified_checkpoints": 5,
    }


def build_teacher_v3_metadata_provenance(
    *,
    training_table: Path,
) -> dict[str, str]:
    """Bind the native-independent metadata baseline to the frozen table."""

    not_applicable = hashlib.sha256(
        b"teacher_v3_metadata_only_no_native_input_v1"
    ).hexdigest()
    return {
        "checkpoint_namespace": TEACHER_V3_METADATA_CHECKPOINT_NAMESPACE,
        "input_contract_version": "teacher_v3_metadata_only_v1",
        "native_h5_sha256": not_applicable,
        "native_manifest_sha256": not_applicable,
        "native_input_mode": "metadata_only_no_hdf5",
        "training_table_sha256": _file_sha256(training_table),
    }


def import_verified_teacher_v3_metadata_baseline(
    *,
    source_dir: Path,
    out_dir: Path,
    training_table: Path,
    expected_train_config: HarmonicTrainConfig,
) -> tuple[dict[str, Any], dict[str, str], dict[str, Any]]:
    """Verify and copy the sealed CPU baseline without opening its test set."""

    import torch

    source_dir = Path(source_dir).resolve()
    out_dir = Path(out_dir).resolve()
    training_table = Path(training_table).resolve()
    if source_dir == out_dir:
        raise ValueError("metadata baseline source and combined output must differ")
    summary_path = source_dir / "summary.json"
    manifest_path = source_dir / "metadata_baseline_checkpoint_manifest.json"
    freeze_path = source_dir / "pretest_model_freeze.json"
    calibration_path = (
        source_dir
        / TEACHER_V3_BASELINE_PROFILE
        / "pooled_oof_calibration.json"
    )
    core_paths = {
        "summary": summary_path,
        "checkpoint_manifest": manifest_path,
        "pretest_model_freeze": freeze_path,
        "pooled_oof_calibration": calibration_path,
    }
    if any(not path.is_file() for path in core_paths.values()):
        raise FileNotFoundError(
            "sealed Teacher v3 metadata baseline lacks its summary, "
            "checkpoint manifest, pretest freeze, or pooled calibration"
        )
    core_hashes_before = {
        name: _file_sha256(path)
        for name, path in core_paths.items()
    }
    summary = json.loads(summary_path.read_text())
    if not isinstance(summary, dict):
        raise RuntimeError("metadata baseline summary must be a JSON mapping")
    provenance = build_teacher_v3_metadata_provenance(
        training_table=training_table,
    )
    expected_summary = {
        "schema_version": (
            "twirl_teacher_v3_metadata_baseline_training_summary_v1"
        ),
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "fixed_role": "metadata_baseline",
        "fixed_test_status": "sealed_pending_primary_profile_freeze",
        "native_inputs_used": False,
        "encoder_pretraining_used": False,
    }
    failures = [
        f"{name}={summary.get(name)!r}, expected {expected!r}"
        for name, expected in expected_summary.items()
        if summary.get(name) != expected
    ]
    if summary.get("profiles") != [TEACHER_V3_BASELINE_PROFILE]:
        failures.append("profiles must contain only metadata_only")
    if summary.get("train_config") != asdict(expected_train_config):
        failures.append("metadata baseline train_config does not match")
    if summary.get("test_metrics") not in ({}, None):
        failures.append("sealed metadata baseline contains fixed-test metrics")
    for name, expected in provenance.items():
        if str(summary.get(name, "")) != str(expected):
            failures.append(
                f"summary {name}={summary.get(name)!r}, expected {expected!r}"
            )

    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict):
        raise RuntimeError(
            "metadata checkpoint manifest must be a JSON mapping"
        )
    manifest_sha256 = core_hashes_before["checkpoint_manifest"]
    if str(
        summary.get("metadata_baseline_checkpoint_manifest_sha256", "")
    ) != manifest_sha256:
        failures.append(
            "summary does not bind the metadata checkpoint-manifest hash"
        )
    if manifest.get("schema_version") != (
        TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA
    ):
        failures.append("metadata checkpoint manifest has wrong schema")
    for name, expected in {
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
    }.items():
        if manifest.get(name) != expected:
            failures.append(
                f"metadata checkpoint manifest {name} is invalid"
            )
    if manifest.get("selected_profile") != TEACHER_V3_BASELINE_PROFILE:
        failures.append("metadata checkpoint manifest has wrong profile")
    if "fixed_test_sealed" not in str(manifest.get("selection_policy", "")):
        failures.append("metadata checkpoint manifest does not seal the test")
    for name, expected in provenance.items():
        if str(manifest.get(name, "")) != str(expected):
            failures.append(
                f"manifest {name}={manifest.get(name)!r}, expected {expected!r}"
            )

    freeze = json.loads(freeze_path.read_text())
    if not isinstance(freeze, dict):
        raise RuntimeError("metadata pretest freeze must be a JSON mapping")
    freeze_sha256 = core_hashes_before["pretest_model_freeze"]
    if str(summary.get("pretest_model_freeze_sha256", "")) != freeze_sha256:
        failures.append("summary does not bind the pretest model-freeze hash")
    expected_freeze = {
        "schema_version": "twirl_teacher_v3_model_freeze_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "model_version": MODEL_VERSION,
        "profile": TEACHER_V3_BASELINE_PROFILE,
        "fixed_role": "metadata_baseline",
        "native_inputs_used": False,
        "encoder_pretraining_used": False,
        "temperature_calibration": (
            "one_temperature_from_concatenated_five_fold_"
            "development_oof_logits"
        ),
        "test_rows_used_for_selection_or_calibration": False,
    }
    for name, expected in expected_freeze.items():
        if freeze.get(name) != expected:
            failures.append(
                f"metadata freeze {name}={freeze.get(name)!r}, "
                f"expected {expected!r}"
            )
    for name, expected in provenance.items():
        if str(freeze.get(name, "")) != str(expected):
            failures.append(
                f"metadata freeze {name}={freeze.get(name)!r}, "
                f"expected {expected!r}"
            )

    calibration = json.loads(calibration_path.read_text())
    if not isinstance(calibration, dict):
        raise RuntimeError(
            "metadata pooled calibration must be a JSON mapping"
        )
    calibration_sha256 = core_hashes_before["pooled_oof_calibration"]
    expected_calibration = {
        "schema_version": TEACHER_V3_CALIBRATION_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "profile": TEACHER_V3_BASELINE_PROFILE,
        "scope": "concatenated_five_fold_development_oof_logits",
    }
    for name, expected in expected_calibration.items():
        if calibration.get(name) != expected:
            failures.append(
                f"metadata calibration {name}={calibration.get(name)!r}, "
                f"expected {expected!r}"
            )
    for name, expected in provenance.items():
        if str(calibration.get(name, "")) != str(expected):
            failures.append(
                f"metadata calibration {name}={calibration.get(name)!r}, "
                f"expected {expected!r}"
            )
    summary_calibration = summary.get("calibration")
    if not isinstance(summary_calibration, dict):
        failures.append("summary calibration must be a mapping")
        summary_calibration = {}
    if str(
        summary_calibration.get("calibration_sha256", "")
    ) != calibration_sha256:
        failures.append(
            "summary does not bind the pooled calibration hash"
        )
    for name in ("profile", "temperature", "oof_logits_sha256"):
        if summary_calibration.get(name) != calibration.get(name):
            failures.append(
                f"summary calibration {name} does not match its artifact"
            )

    records = manifest.get("checkpoints")
    if not isinstance(records, list) or len(records) != 5:
        failures.append("metadata checkpoint manifest must contain five folds")
        records = []
    observed_folds: list[int] = []
    for record in records:
        fold = int(record.get("fold", -1))
        observed_folds.append(fold)
        expected_relative = (
            Path(TEACHER_V3_BASELINE_PROFILE)
            / f"fold_{fold}"
            / "teacher.pt"
        )
        if str(record.get("path", "")) != expected_relative.as_posix():
            failures.append(f"metadata fold {fold} has an unsafe path")
            continue
        checkpoint_path = source_dir / expected_relative
        if (
            not checkpoint_path.is_file()
            or _file_sha256(checkpoint_path) != str(record.get("sha256", ""))
        ):
            failures.append(f"metadata fold {fold} checkpoint hash mismatch")
            continue
        if str(
            record.get("pooled_oof_calibration_sha256", "")
        ) != calibration_sha256:
            failures.append(
                f"metadata fold {fold} calibration hash mismatch"
            )
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        try:
            validate_teacher_v3_input_provenance(
                checkpoint,
                expected=provenance,
                artifact=f"sealed metadata baseline fold-{fold}",
            )
        except RuntimeError as exc:
            failures.append(str(exc))
        if checkpoint.get("profile") != TEACHER_V3_BASELINE_PROFILE:
            failures.append(f"metadata fold {fold} has wrong profile")
        for name, expected in {
            "run_id": TEACHER_V3_RUN_ID,
            "release_name": TEACHER_V3_RUN_NAME,
            "model_version": MODEL_VERSION,
        }.items():
            if checkpoint.get(name) != expected:
                failures.append(
                    f"metadata fold {fold} has wrong {name}"
                )
        if int(checkpoint.get("fold", -1)) != fold:
            failures.append(f"metadata fold {fold} has wrong fold identity")
        if checkpoint.get("temperature_calibration_scope") != (
            "pooled_oof_development"
        ):
            failures.append(f"metadata fold {fold} lacks pooled calibration")
        if str(
            checkpoint.get("pooled_oof_calibration_sha256", "")
        ) != calibration_sha256:
            failures.append(
                f"metadata fold {fold} internal calibration hash mismatch"
            )
        try:
            temperature_matches = bool(
                np.isclose(
                    float(checkpoint.get("temperature", float("nan"))),
                    float(calibration.get("temperature", float("nan"))),
                )
            )
        except (TypeError, ValueError):
            temperature_matches = False
        if not temperature_matches:
            failures.append(
                f"metadata fold {fold} temperature does not match calibration"
            )
        if not bool(checkpoint.get("encoder_pretraining_skipped", False)):
            failures.append(f"metadata fold {fold} used encoder pretraining")
        if checkpoint.get("train_config") != asdict(expected_train_config):
            failures.append(f"metadata fold {fold} has wrong train_config")
    if sorted(observed_folds) != list(range(5)):
        failures.append(
            f"metadata folds must be 0..4; observed={observed_folds}"
        )
    if failures:
        raise RuntimeError(
            "metadata baseline verification failed: " + "; ".join(failures)
        )

    source_profile = source_dir / TEACHER_V3_BASELINE_PROFILE
    destination_profile = out_dir / TEACHER_V3_BASELINE_PROFILE
    source_hashes = {
        path.relative_to(source_profile).as_posix(): _file_sha256(path)
        for path in sorted(source_profile.rglob("*"))
        if path.is_file()
    }
    if any(
        any(part.startswith("fixed_test_") for part in Path(name).parts)
        for name in source_hashes
    ):
        raise RuntimeError(
            "sealed metadata baseline directory already contains fixed-test "
            "artifacts"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    if destination_profile.exists():
        raise FileExistsError(
            "full Teacher v3 import will not resume an existing metadata "
            f"destination: {destination_profile}"
        )
    shutil.copytree(source_profile, destination_profile)
    copied_hashes = {
        path.relative_to(destination_profile).as_posix(): _file_sha256(path)
        for path in sorted(destination_profile.rglob("*"))
        if path.is_file()
    }
    if copied_hashes != source_hashes:
        raise RuntimeError("metadata baseline changed while it was copied")
    source_hashes_after = {
        path.relative_to(source_profile).as_posix(): _file_sha256(path)
        for path in sorted(source_profile.rglob("*"))
        if path.is_file()
    }
    core_hashes_after = {
        name: _file_sha256(path)
        for name, path in core_paths.items()
    }
    if source_hashes_after != source_hashes:
        raise RuntimeError(
            "metadata baseline source profile mutated while it was copied"
        )
    if core_hashes_after != core_hashes_before:
        raise RuntimeError(
            "metadata baseline control artifact mutated during import"
        )
    if _file_sha256(training_table) != provenance[
        "training_table_sha256"
    ]:
        raise RuntimeError(
            "Teacher v3 training table mutated during metadata-baseline import"
        )

    imported_calibration = dict(summary_calibration)
    imported_calibration["calibration_path"] = str(
        destination_profile / "pooled_oof_calibration.json"
    )
    imported_calibration["calibration_sha256"] = calibration_sha256
    development_records = summary.get(
        "development_fixed_profile_comparison",
        [],
    )
    if not isinstance(development_records, list) or len(
        development_records
    ) != 1:
        raise RuntimeError(
            "sealed metadata baseline lacks one development comparison record"
        )
    reuse_audit = {
        "schema_version": "twirl_teacher_v3_metadata_baseline_reuse_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "source_dir": str(source_dir),
        "source_summary": str(summary_path),
        "source_summary_sha256": core_hashes_before["summary"],
        "source_checkpoint_manifest": str(manifest_path),
        "source_checkpoint_manifest_sha256": manifest_sha256,
        "source_pretest_model_freeze": str(freeze_path),
        "source_pretest_model_freeze_sha256": freeze_sha256,
        "source_calibration_sha256": calibration_sha256,
        "source_profile_tree_sha256": hashlib.sha256(
            _json_payload(source_hashes)
        ).hexdigest(),
        "n_verified_checkpoints": 5,
        "n_copied_files": int(len(copied_hashes)),
        "fixed_test_state_at_import": "sealed",
        "test_rows_used_for_import_decision": False,
        "destination_profile": str(destination_profile),
    }
    _write_json(out_dir / "metadata_baseline_reuse.json", reuse_audit)
    return (
        imported_calibration,
        provenance,
        {
            **reuse_audit,
            "development_record": dict(development_records[0]),
        },
    )


def run_teacher_v3_metadata_baseline(
    *,
    training_table: Path,
    out_dir: Path,
    train_config: HarmonicTrainConfig = HarmonicTrainConfig(seed=560062),
    workers: int = 8,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Train a genuine seven-sector metadata baseline without native HDF5s."""

    training_table = Path(training_table).resolve()
    out_dir = require_fresh_teacher_v3_output_dir(
        Path(out_dir),
        artifact="Teacher v3 metadata baseline",
    )
    source, initial_training_table_sha256, initial_read_audit = (
        _read_training_table_with_stable_hash(
            training_table,
            artifact="Teacher v3 metadata-baseline training table",
        )
    )
    validate_teacher_v3_training_table(source)
    rows = prepare_harmonic_training_rows(
        source,
        seed=train_config.seed,
    )
    expected_active = int(
        _truth(source["teacher_v3_training_include"]).sum()
    )
    if len(rows) != expected_active:
        raise RuntimeError(
            "Teacher v3 active-row contract disagrees with model targets: "
            f"{len(rows)} != {expected_active}"
        )
    provenance = build_teacher_v3_metadata_provenance(
        training_table=training_table,
    )
    if provenance["training_table_sha256"] != (
        initial_training_table_sha256
    ):
        raise RuntimeError(
            "Teacher v3 metadata-baseline table changed after its initial read"
        )
    rows.to_csv(out_dir / "training_rows_with_fixed_splits.csv", index=False)
    injection_audit = injection_truth_human_audit(
        source,
        out_dir=out_dir / "injection_truth_human_audit",
    )
    for fold in range(5):
        _train_one_fold(
            rows=rows,
            native_h5=None,
            out_dir=out_dir,
            profile=TEACHER_V3_BASELINE_PROFILE,
            fold=fold,
            train_config=train_config,
            workers=workers,
            pretrain_epochs=0,
            require_cuda=require_cuda,
            input_provenance=provenance,
            defer_temperature_calibration=True,
            pretraining_cache_namespace="not_applicable_metadata_only",
            skip_encoder_pretraining=True,
        )
    calibration = calibrate_teacher_v3_profile_oof(
        rows=rows,
        out_dir=out_dir,
        profile=TEACHER_V3_BASELINE_PROFILE,
        input_provenance=provenance,
    )
    development = _profile_development_record(calibration)
    pd.DataFrame([development]).to_csv(
        out_dir / "development_fixed_profile_comparison.csv",
        index=False,
    )
    freeze = {
        "schema_version": "twirl_teacher_v3_model_freeze_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "model_version": MODEL_VERSION,
        "profile": TEACHER_V3_BASELINE_PROFILE,
        "fixed_role": "metadata_baseline",
        "native_inputs_used": False,
        "encoder_pretraining_used": False,
        "temperature_calibration": (
            "one_temperature_from_concatenated_five_fold_"
            "development_oof_logits"
        ),
        "test_rows_used_for_selection_or_calibration": False,
        "training_table_initial_read": initial_read_audit,
        **provenance,
    }
    freeze_path = out_dir / "pretest_model_freeze.json"
    freeze_sha256 = _write_json(freeze_path, freeze)
    if _file_sha256(training_table) != initial_training_table_sha256:
        raise RuntimeError(
            "Teacher v3 training table changed during metadata training"
        )
    manifest = _write_selected_checkpoint_manifest(
        out_dir=out_dir,
        input_provenance=provenance,
        selected_profile=TEACHER_V3_BASELINE_PROFILE,
        filename="metadata_baseline_checkpoint_manifest.json",
        selection_policy=(
            "fixed_metadata_baseline; fixed_test_sealed_until_primary_freeze"
        ),
    )
    if _file_sha256(training_table) != initial_training_table_sha256:
        raise RuntimeError(
            "Teacher v3 training table changed before metadata summary"
        )
    summary: dict[str, Any] = {
        "schema_version": (
            "twirl_teacher_v3_metadata_baseline_training_summary_v1"
        ),
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profiles": [TEACHER_V3_BASELINE_PROFILE],
        "fixed_role": "metadata_baseline",
        "native_inputs_used": False,
        "encoder_pretraining_used": False,
        "train_config": asdict(train_config),
        "n_training_rows": int(len(rows)),
        "n_development_rows": int(
            rows["fixed_split"].eq("development").sum()
        ),
        "n_fixed_test_rows": int(rows["fixed_split"].eq("test").sum()),
        "training_table_initial_read": initial_read_audit,
        **provenance,
        "calibration": calibration,
        "development_fixed_profile_comparison": [development],
        "fixed_test_status": "sealed_pending_primary_profile_freeze",
        "test_metrics": {},
        "pretest_model_freeze": str(freeze_path),
        "pretest_model_freeze_sha256": freeze_sha256,
        "metadata_baseline_checkpoint_manifest": str(manifest),
        "metadata_baseline_checkpoint_manifest_sha256": _file_sha256(
            manifest
        ),
        "injection_truth_human_audit": injection_audit,
        "student_training_blocked": True,
        "automatic_production_promotion": False,
    }
    _write_json(out_dir / "summary.json", summary)
    return summary


def _assert_native_inputs_unchanged(
    manifest: Mapping[str, Any],
) -> None:
    for record in manifest["native_files"]:
        path = Path(record["path"])
        if (
            path.stat().st_size != int(record["size_bytes"])
            or path.stat().st_mtime_ns != int(record["mtime_ns"])
            or file_sha256(path) != str(record["sha256"])
        ):
            raise RuntimeError(
                f"Teacher v3 native input changed during training: {path}"
            )
    registry = manifest["registry"]
    for path_name, digest_name in (
        ("path", "sha256"),
        ("summary_path", "summary_sha256"),
    ):
        path = Path(registry[path_name])
        if _file_sha256(path) != str(registry[digest_name]):
            raise RuntimeError(
                f"Teacher v3 native registry input changed: {path}"
            )


def _assert_teacher_v3_inputs_unchanged(
    *,
    native_manifest: Mapping[str, Any],
    training_table: Path,
    training_table_sha256: str,
) -> None:
    """Recheck every frozen table/registry/HDF5 byte binding."""

    _assert_native_inputs_unchanged(native_manifest)
    training_table = Path(training_table).resolve()
    if _file_sha256(training_table) != str(training_table_sha256):
        raise RuntimeError(
            "Teacher v3 training table changed after its initial read: "
            f"{training_table}"
        )


def run_teacher_v3_training(
    *,
    training_table: Path,
    native_registry: Path,
    native_registry_summary: Path,
    out_dir: Path,
    train_config: HarmonicTrainConfig = HarmonicTrainConfig(seed=560062),
    workers: int = 8,
    pretrain_epochs: int = 20,
    require_cuda: bool = True,
    metadata_baseline_dir: Path | None = None,
    bootstrap_draws: int = 2000,
    bootstrap_seed: int = 560063,
    bootstrap_confidence_level: float = 0.95,
) -> dict[str, Any]:
    """Train both predeclared primary policies, then open the fixed test once."""

    training_table = Path(training_table).resolve()
    out_dir = require_fresh_teacher_v3_output_dir(
        Path(out_dir),
        artifact="full Teacher v3 training",
    )
    if (
        int(bootstrap_draws) != 2000
        or int(bootstrap_seed) != 560063
        or not np.isclose(float(bootstrap_confidence_level), 0.95)
    ):
        raise ValueError(
            "frozen Teacher v3 evaluation requires 2,000 TIC-cluster draws, "
            "seed 560063, and 95% intervals"
        )
    source, initial_training_table_sha256, initial_read_audit = (
        _read_training_table_with_stable_hash(
            training_table,
            artifact="Teacher v3 full-run training table",
        )
    )
    validate_teacher_v3_training_table(source)
    rows = prepare_harmonic_training_rows(
        source,
        seed=train_config.seed,
    )
    expected_active = int(
        _truth(source["teacher_v3_training_include"]).sum()
    )
    if len(rows) != expected_active:
        raise RuntimeError(
            "Teacher v3 active-row contract disagrees with model targets: "
            f"{len(rows)} != {expected_active}"
        )
    native_manifest = build_teacher_v3_native_manifest(
        rows=rows,
        registry_path=native_registry,
        registry_summary_path=native_registry_summary,
    )
    native_manifest_path = out_dir / "native_input_manifest.json"
    native_manifest_sha256 = _write_json(
        native_manifest_path,
        native_manifest,
    )
    input_provenance = build_teacher_v3_input_provenance(
        training_table=training_table,
        native_manifest=native_manifest,
    )
    if input_provenance["training_table_sha256"] != (
        initial_training_table_sha256
    ):
        raise RuntimeError(
            "Teacher v3 training table changed after its initial read"
        )
    if input_provenance["native_manifest_sha256"] != (
        native_manifest_sha256
    ):
        raise RuntimeError("Teacher v3 native manifest serialization changed")
    rows.to_csv(out_dir / "training_rows_with_fixed_splits.csv", index=False)
    sensitivity_rows, sensitivity_audit = (
        prepare_teacher_v3_uncertain_masked_rows(rows)
    )
    sensitivity_root = (
        out_dir / "sensitivities" / "uncertain_masked"
    )
    sensitivity_root.mkdir(parents=True, exist_ok=True)
    sensitivity_rows.to_csv(
        sensitivity_root / "training_rows_with_fixed_splits.csv",
        index=False,
        lineterminator="\n",
    )
    sensitivity_provenance = (
        build_teacher_v3_uncertain_masked_provenance(
            input_provenance=input_provenance,
            sensitivity_audit=sensitivity_audit,
        )
    )
    metadata_provenance = build_teacher_v3_metadata_provenance(
        training_table=training_table,
    )
    injection_audit = injection_truth_human_audit(
        source,
        out_dir=out_dir / "injection_truth_human_audit",
    )

    fold_results: list[dict[str, Any]] = []
    metadata_reuse: dict[str, Any] | None = None
    if metadata_baseline_dir is not None:
        (
            metadata_calibration,
            imported_metadata_provenance,
            metadata_reuse,
        ) = import_verified_teacher_v3_metadata_baseline(
            source_dir=metadata_baseline_dir,
            out_dir=out_dir,
            training_table=training_table,
            expected_train_config=train_config,
        )
        if imported_metadata_provenance != metadata_provenance:
            raise RuntimeError(
                "imported metadata baseline provenance changed after validation"
            )
    else:
        for fold in range(5):
            fold_results.append(
                _train_one_fold(
                    rows=rows,
                    native_h5=None,
                    out_dir=out_dir,
                    profile=TEACHER_V3_BASELINE_PROFILE,
                    fold=fold,
                    train_config=train_config,
                    workers=workers,
                    pretrain_epochs=0,
                    require_cuda=require_cuda,
                    input_provenance=metadata_provenance,
                    defer_temperature_calibration=True,
                    pretraining_cache_namespace="not_applicable_metadata_only",
                    skip_encoder_pretraining=True,
                )
            )
        metadata_calibration = calibrate_teacher_v3_profile_oof(
            rows=rows,
            out_dir=out_dir,
            profile=TEACHER_V3_BASELINE_PROFILE,
            input_provenance=metadata_provenance,
        )

    for fold in range(5):
        fold_results.append(
            _train_one_fold(
                rows=rows,
                native_h5=None,
                out_dir=out_dir,
                profile=TEACHER_V3_PRIMARY_PROFILE,
                fold=fold,
                train_config=train_config,
                workers=workers,
                pretrain_epochs=pretrain_epochs,
                require_cuda=require_cuda,
                input_provenance=input_provenance,
                defer_temperature_calibration=True,
                pretraining_cache_namespace=(
                    "teacher_v3_s56_s62_native_v1"
                ),
                skip_encoder_pretraining=False,
            )
        )
    primary_calibration = calibrate_teacher_v3_profile_oof(
        rows=rows,
        out_dir=out_dir,
        profile=TEACHER_V3_PRIMARY_PROFILE,
        input_provenance=input_provenance,
    )

    for fold in range(5):
        fold_results.append(
            _train_one_fold(
                rows=sensitivity_rows,
                native_h5=None,
                out_dir=sensitivity_root,
                profile=TEACHER_V3_PRIMARY_PROFILE,
                fold=fold,
                train_config=train_config,
                workers=workers,
                pretrain_epochs=pretrain_epochs,
                require_cuda=require_cuda,
                input_provenance=sensitivity_provenance,
                defer_temperature_calibration=True,
                pretraining_cache_namespace=(
                    "teacher_v3_s56_s62_native_uncertain_masked_v1"
                ),
                skip_encoder_pretraining=False,
            )
        )
    sensitivity_calibration = calibrate_teacher_v3_profile_oof(
        rows=sensitivity_rows,
        out_dir=sensitivity_root,
        profile=TEACHER_V3_PRIMARY_PROFILE,
        input_provenance=sensitivity_provenance,
    )

    calibration_results = [
        metadata_calibration,
        primary_calibration,
    ]
    development = pd.DataFrame(
        [
            _profile_development_record(result)
            for result in calibration_results
        ]
    )
    development.to_csv(
        out_dir / "development_fixed_profile_comparison.csv",
        index=False,
    )
    sensitivity_development = _profile_development_record(
        sensitivity_calibration
    )
    sensitivity_development.update(
        {
            "fixed_role": "label_policy_sensitivity",
            "label_policy": TEACHER_V3_UNCERTAIN_MASKED_POLICY,
            "training_execution": "independent_five_fold_retraining",
        }
    )
    pd.DataFrame([sensitivity_development]).to_csv(
        sensitivity_root / "development_metrics.csv",
        index=False,
        lineterminator="\n",
    )

    metadata_selection_policy = (
        "fixed_metadata_baseline; fixed_test_sealed_until_primary_freeze"
    )
    primary_selection_policy = "fixed_before_test"
    sensitivity_selection_policy = (
        "predeclared_uncertain_masked_retraining; fixed_before_test"
    )
    metadata_manifest = _write_selected_checkpoint_manifest(
        out_dir=out_dir,
        input_provenance=metadata_provenance,
        selected_profile=TEACHER_V3_BASELINE_PROFILE,
        filename="metadata_baseline_checkpoint_manifest.json",
        selection_policy=metadata_selection_policy,
    )
    selected_manifest = _write_selected_checkpoint_manifest(
        out_dir=out_dir,
        input_provenance=input_provenance,
        selection_policy=primary_selection_policy,
    )
    sensitivity_manifest = _write_selected_checkpoint_manifest(
        out_dir=sensitivity_root,
        input_provenance=sensitivity_provenance,
        selected_profile=TEACHER_V3_PRIMARY_PROFILE,
        filename="uncertain_masked_checkpoint_manifest.json",
        selection_policy=sensitivity_selection_policy,
    )
    metadata_manifest_sha256 = _file_sha256(metadata_manifest)
    selected_manifest_sha256 = _file_sha256(selected_manifest)
    sensitivity_manifest_sha256 = _file_sha256(sensitivity_manifest)
    freeze = {
        "schema_version": "twirl_teacher_v3_model_freeze_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "model_version": MODEL_VERSION,
        "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
        "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
        "profile_selection_policy": "fixed_before_test",
        "temperature_calibration": (
            "one_temperature_per_profile_from_concatenated_"
            "five_fold_development_oof_logits"
        ),
        "label_policy_sensitivity": {
            **sensitivity_audit,
            "training_execution": "independent_five_fold_retraining",
            "calibration_scope": (
                "concatenated_five_fold_development_oof_logits"
            ),
            "checkpoint_manifest": str(sensitivity_manifest),
            "checkpoint_manifest_sha256": sensitivity_manifest_sha256,
        },
        "metadata_baseline": {
            "provenance": metadata_provenance,
            "reused_sealed_cpu_run": bool(metadata_reuse is not None),
            "reuse_audit": metadata_reuse,
            "checkpoint_manifest": str(metadata_manifest),
            "checkpoint_manifest_sha256": metadata_manifest_sha256,
        },
        "primary_checkpoint_manifest": str(selected_manifest),
        "primary_checkpoint_manifest_sha256": selected_manifest_sha256,
        "bootstrap": {
            "unit": "tic",
            "draws": int(bootstrap_draws),
            "seed": int(bootstrap_seed),
            "confidence_level": float(bootstrap_confidence_level),
        },
        "test_rows_used_for_selection_or_calibration": False,
        "training_table_initial_read": initial_read_audit,
        **input_provenance,
    }
    _assert_teacher_v3_inputs_unchanged(
        native_manifest=native_manifest,
        training_table=training_table,
        training_table_sha256=initial_training_table_sha256,
    )
    freeze_path = out_dir / "pretest_model_freeze.json"
    freeze_sha256 = _write_json(freeze_path, freeze)

    _assert_teacher_v3_inputs_unchanged(
        native_manifest=native_manifest,
        training_table=training_table,
        training_table_sha256=initial_training_table_sha256,
    )

    def verify_frozen_ensembles() -> dict[str, dict[str, Any]]:
        return {
            TEACHER_V3_BASELINE_PROFILE: (
                verify_teacher_v3_checkpoint_manifest(
                    manifest_path=metadata_manifest,
                    expected_manifest_sha256=metadata_manifest_sha256,
                    expected_profile=TEACHER_V3_BASELINE_PROFILE,
                    expected_selection_policy=metadata_selection_policy,
                    expected_provenance=metadata_provenance,
                    calibration_path=Path(
                        metadata_calibration["calibration_path"]
                    ),
                    expected_calibration_sha256=str(
                        metadata_calibration["calibration_sha256"]
                    ),
                )
            ),
            TEACHER_V3_PRIMARY_PROFILE: (
                verify_teacher_v3_checkpoint_manifest(
                    manifest_path=selected_manifest,
                    expected_manifest_sha256=selected_manifest_sha256,
                    expected_profile=TEACHER_V3_PRIMARY_PROFILE,
                    expected_selection_policy=primary_selection_policy,
                    expected_provenance=input_provenance,
                    calibration_path=Path(
                        primary_calibration["calibration_path"]
                    ),
                    expected_calibration_sha256=str(
                        primary_calibration["calibration_sha256"]
                    ),
                )
            ),
            "uncertain_masked": verify_teacher_v3_checkpoint_manifest(
                manifest_path=sensitivity_manifest,
                expected_manifest_sha256=sensitivity_manifest_sha256,
                expected_profile=TEACHER_V3_PRIMARY_PROFILE,
                expected_selection_policy=sensitivity_selection_policy,
                expected_provenance=sensitivity_provenance,
                calibration_path=Path(
                    sensitivity_calibration["calibration_path"]
                ),
                expected_calibration_sha256=str(
                    sensitivity_calibration["calibration_sha256"]
                ),
            ),
        }

    pretest_checkpoint_verification = verify_frozen_ensembles()
    test_metrics = {
        TEACHER_V3_BASELINE_PROFILE: _evaluate_fixed_profile(
            rows=rows,
            out_dir=out_dir,
            profile=TEACHER_V3_BASELINE_PROFILE,
            train_config=train_config,
            workers=workers,
            input_provenance=metadata_provenance,
            require_cuda=require_cuda,
            label_policy=TEACHER_V3_PRIMARY_LABEL_POLICY,
            bootstrap_draws=bootstrap_draws,
            bootstrap_seed=bootstrap_seed,
            bootstrap_confidence_level=bootstrap_confidence_level,
        ),
        TEACHER_V3_PRIMARY_PROFILE: _evaluate_fixed_profile(
            rows=rows,
            out_dir=out_dir,
            profile=TEACHER_V3_PRIMARY_PROFILE,
            train_config=train_config,
            workers=workers,
            input_provenance=input_provenance,
            require_cuda=require_cuda,
            label_policy=TEACHER_V3_PRIMARY_LABEL_POLICY,
            bootstrap_draws=bootstrap_draws,
            bootstrap_seed=bootstrap_seed,
            bootstrap_confidence_level=bootstrap_confidence_level,
        ),
    }
    sensitivity_test_metrics = _evaluate_fixed_profile(
        rows=sensitivity_rows,
        out_dir=sensitivity_root,
        profile=TEACHER_V3_PRIMARY_PROFILE,
        train_config=train_config,
        workers=workers,
        input_provenance=sensitivity_provenance,
        require_cuda=require_cuda,
        label_policy=TEACHER_V3_UNCERTAIN_MASKED_POLICY,
        bootstrap_draws=bootstrap_draws,
        bootstrap_seed=bootstrap_seed,
        bootstrap_confidence_level=bootstrap_confidence_level,
    )
    for profile in TEACHER_V3_PROFILES:
        profile_dir = out_dir / profile
        test_metrics[profile][
            "non_uncertain_common_support_bootstrap"
        ] = write_teacher_v3_non_uncertain_comparison_bootstrap(
            predictions_path=profile_dir / "fixed_test_predictions.csv",
            profile=profile,
            out_dir=profile_dir,
            draws=bootstrap_draws,
            seed=bootstrap_seed,
            confidence_level=bootstrap_confidence_level,
        )
        _write_json(
            profile_dir / "fixed_test_metrics.json",
            test_metrics[profile],
        )
    _assert_teacher_v3_inputs_unchanged(
        native_manifest=native_manifest,
        training_table=training_table,
        training_table_sha256=initial_training_table_sha256,
    )
    posttest_checkpoint_verification = verify_frozen_ensembles()
    if posttest_checkpoint_verification != pretest_checkpoint_verification:
        raise RuntimeError(
            "Teacher v3 frozen-ensemble audit changed across fixed-test "
            "evaluation"
        )
    summary: dict[str, Any] = {
        "schema_version": "twirl_teacher_v3_training_summary_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profiles": list(TEACHER_V3_PROFILES),
        "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
        "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
        "profile_selection_policy": "fixed_before_test",
        "train_config": asdict(train_config),
        "pretraining_epochs": int(pretrain_epochs),
        "h200_trained_fold_count": (
            10 if metadata_reuse is not None else 15
        ),
        "metadata_baseline_reused": bool(metadata_reuse is not None),
        "n_training_rows": int(len(rows)),
        "n_development_rows": int(
            rows["fixed_split"].eq("development").sum()
        ),
        "n_fixed_test_rows": int(rows["fixed_split"].eq("test").sum()),
        "training_table_initial_read": initial_read_audit,
        "input_channel_contract": {
            name: list(channels)
            for name, channels in CHANNEL_CONTRACT.items()
        },
        **input_provenance,
        "native_input_manifest": str(native_manifest_path),
        "pretest_model_freeze": str(freeze_path),
        "pretest_model_freeze_sha256": freeze_sha256,
        "calibration": calibration_results,
        "development_fixed_profile_comparison": (
            development.to_dict("records")
        ),
        "test_metrics": test_metrics,
        "fixed_test_status": "opened_once_after_pretest_model_freeze",
        "bootstrap": {
            "unit": "tic",
            "draws": int(bootstrap_draws),
            "seed": int(bootstrap_seed),
            "confidence_level": float(bootstrap_confidence_level),
        },
        "sensitivity_analyses": {
            "uncertain_masked": {
                **sensitivity_audit,
                "training_execution": "independent_five_fold_retraining",
                "calibration": sensitivity_calibration,
                "development_metrics": sensitivity_development,
                "test_metrics": sensitivity_test_metrics,
                "checkpoint_manifest": str(sensitivity_manifest),
                "checkpoint_manifest_sha256": sensitivity_manifest_sha256,
            }
        },
        "selected_checkpoint_manifest": str(selected_manifest),
        "selected_checkpoint_manifest_sha256": selected_manifest_sha256,
        "metadata_baseline_checkpoint_manifest": str(metadata_manifest),
        "metadata_baseline_checkpoint_manifest_sha256": (
            metadata_manifest_sha256
        ),
        "checkpoint_freeze_verification": {
            "before_fixed_test": pretest_checkpoint_verification,
            "after_fixed_test": posttest_checkpoint_verification,
        },
        "metadata_baseline_reuse": metadata_reuse,
        "injection_truth_human_audit": injection_audit,
        "student_training_blocked": True,
        "automatic_production_promotion": False,
    }
    _write_json(out_dir / "summary.json", summary)
    return summary


__all__ = [
    "TEACHER_V3_BASELINE_PROFILE",
    "TEACHER_V3_BOOTSTRAP_SCHEMA",
    "TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA",
    "TEACHER_V3_CHECKPOINT_NAMESPACE",
    "TEACHER_V3_METADATA_CHECKPOINT_NAMESPACE",
    "TEACHER_V3_PRIMARY_LABEL_POLICY",
    "TEACHER_V3_PRIMARY_PROFILE",
    "TEACHER_V3_PROFILES",
    "TEACHER_V3_RUN_ID",
    "TEACHER_V3_UNCERTAIN_MASKED_POLICY",
    "build_teacher_v3_input_provenance",
    "build_teacher_v3_metadata_provenance",
    "build_teacher_v3_native_manifest",
    "build_teacher_v3_uncertain_masked_provenance",
    "calibrate_teacher_v3_profile_oof",
    "import_verified_teacher_v3_metadata_baseline",
    "prepare_teacher_v3_uncertain_masked_rows",
    "run_teacher_v3_training",
    "run_teacher_v3_metadata_baseline",
    "teacher_v3_tic_cluster_bootstrap",
    "verify_teacher_v3_checkpoint_manifest",
    "validate_teacher_v3_input_provenance",
    "validate_teacher_v3_training_table",
    "write_teacher_v3_fixed_test_bootstrap",
    "write_teacher_v3_non_uncertain_comparison_bootstrap",
]
