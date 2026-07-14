"""Grouped-fold training for the hybrid S56 Teacher-v2 CNN."""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    HARMONIC_CLASSES,
    MODEL_PROFILES,
    MORPHOLOGY_CLASSES,
    PRESERVE_CLASSES,
    HarmonicModelConfig,
    profile_branches,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    MetadataNormalization,
    collate_native_samples,
)
from twirl.vetting.harmonic_inputs import (
    CHANNEL_CONTRACT,
    RAW_PAIR_CONTRACT_VERSION,
    verify_raw_pair_contract,
)
from twirl.vetting.harmonic_training import (
    classification_metrics,
    expected_calibration_error,
    fit_temperature,
)
from twirl.vetting.teacher_v2 import (
    COMPACT_CLASSES,
    TEACHER_V2_MODEL_VERSION,
    TEACHER_V2_ROLE_POLICY,
    assert_teacher_v2_feature_columns,
)
from twirl.vetting.teacher_v2_cnn import (
    TeacherV2TrainConfig,
    build_teacher_v2_cnn,
    teacher_v2_multitask_loss,
)


TEACHER_V2_METADATA_COLUMNS: tuple[str, ...] = (
    "period_d",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "tmag",
    "adp_sml_peak_rank",
    "adp_sml_period_d",
    "adp_sml_duration_min",
    "adp_sml_depth",
    "adp_sml_depth_snr",
    "adp_sml_sde",
    "adp_sml_log_power",
    "adp_peak_rank",
    "adp_period_d",
    "adp_duration_min",
    "adp_depth",
    "adp_depth_snr",
    "adp_sde",
    "adp_log_power",
    "aperture_period_rel_delta",
    "aperture_depth_ratio_primary_over_small",
    "aperture_disagreement_flag",
)

DEFAULT_TEACHER_V2_PROFILES: tuple[str, ...] = (
    "metadata_only",
    "single_period_native_fold",
    "seven_harmonic_shape",
    "shape_plus_raw_chronology",
    "shape_plus_periodogram_bls",
    "full_combined",
)


def _softmax(values: np.ndarray, *, temperature: float = 1.0) -> np.ndarray:
    scaled = np.asarray(values, dtype=np.float64) / max(float(temperature), 1.0e-6)
    scaled -= np.max(scaled, axis=1, keepdims=True)
    exponential = np.exp(scaled)
    return exponential / np.maximum(exponential.sum(axis=1, keepdims=True), 1.0e-12)


def build_teacher_v2_metadata_matrix(
    rows: pd.DataFrame,
    *,
    fit_mask: np.ndarray,
    normalization: MetadataNormalization | None = None,
) -> tuple[np.ndarray, MetadataNormalization]:
    """Build the explicitly shared injection/real scalar feature branch."""

    columns = (
        normalization.columns
        if normalization is not None
        else tuple(column for column in TEACHER_V2_METADATA_COLUMNS if column in rows)
    )
    assert_teacher_v2_feature_columns(columns)
    values = (
        np.column_stack(
            [pd.to_numeric(rows[column], errors="coerce").to_numpy(dtype=float) for column in columns]
        )
        if columns
        else np.empty((len(rows), 0), dtype=float)
    )
    fit_mask = np.asarray(fit_mask, dtype=bool)
    if fit_mask.shape != (len(rows),):
        raise ValueError("fit_mask must have one value per Teacher-v2 row")
    if normalization is None:
        if not fit_mask.any():
            raise ValueError("metadata normalization requires at least one fit row")
        center = np.nanmedian(values[fit_mask], axis=0) if values.shape[1] else np.empty(0)
        scale = np.nanstd(values[fit_mask], axis=0) if values.shape[1] else np.empty(0)
        center = np.where(np.isfinite(center), center, 0.0)
        scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
        normalization = MetadataNormalization(
            columns=tuple(columns),
            center=tuple(float(value) for value in center),
            scale=tuple(float(value) for value in scale),
        )
    center = np.asarray(normalization.center, dtype=float)
    scale = np.asarray(normalization.scale, dtype=float)
    values = np.where(np.isfinite(values), values, center)
    values = ((values - center) / scale).astype(np.float32)
    # Paired originals have no independently searched BLS metadata. Keep their
    # scalar branch neutral; profiles with BLS/periodograms exclude these rows
    # from compact-head optimization below.
    paired = rows.get("teacher_v2_role", pd.Series("", index=rows.index)).astype(str).eq(
        "paired_pre_injection"
    )
    values[paired.to_numpy()] = 0.0
    return values, normalization


def _inverse_sqrt_weights(target: np.ndarray, *, cap: float) -> np.ndarray:
    target = np.asarray(target, dtype=int)
    weights = np.zeros(len(target), dtype=np.float64)
    active = target >= 0
    if not active.any():
        return weights.astype(np.float32)
    counts = pd.Series(target[active]).value_counts()
    largest = float(counts.max())
    for value, count in counts.items():
        weights[target == int(value)] = min(float(cap), np.sqrt(largest / float(count)))
    weights[active] /= np.mean(weights[active])
    return weights.astype(np.float32)


def _equalize_aggregate(weights: np.ndarray, masks: Sequence[np.ndarray]) -> None:
    active_masks = [np.asarray(mask, dtype=bool) for mask in masks if np.any(mask)]
    totals = [float(weights[mask].sum()) for mask in active_masks]
    totals = [value for value in totals if value > 0]
    if len(totals) < 2:
        return
    target = float(np.mean(totals))
    for mask in active_masks:
        total = float(weights[mask].sum())
        if total > 0:
            weights[mask] *= target / total


def attach_teacher_v2_fold_weights(
    rows: pd.DataFrame,
    *,
    fit_mask: np.ndarray,
    profile: str,
    class_weight_cap: float = 4.0,
) -> pd.DataFrame:
    """Fit task weights from one grouped training partition only."""

    fit_mask = np.asarray(fit_mask, dtype=bool)
    if fit_mask.shape != (len(rows),):
        raise ValueError("fit_mask must have one value per Teacher-v2 row")
    out = rows.copy()
    for name in ("morphology", "preserve", "harmonic", "compact"):
        out[f"{name}_weight"] = np.float32(0.0)
    fit = out.loc[fit_mask].copy()
    for name in ("morphology", "preserve", "harmonic", "compact"):
        target = fit[f"{name}_target_index"].to_numpy(dtype=int)
        weights = _inverse_sqrt_weights(target, cap=class_weight_cap)
        if name == "morphology":
            human = fit.get("human_label", pd.Series("", index=fit.index)).fillna("").astype(str)
            other = target == MORPHOLOGY_CLASSES.index("other")
            _equalize_aggregate(
                weights,
                [
                    other & human.eq("instrumental_or_systematic").to_numpy(),
                    other
                    & human.isin({"uncertain", "no_visible_signal"}).to_numpy(),
                ],
            )
        elif name == "compact":
            role = fit["teacher_v2_role"].fillna("").astype(str)
            _equalize_aggregate(weights, [target == 0, target == 1])
            negative = target == 0
            _equalize_aggregate(
                weights,
                [
                    negative & role.eq("same_lc_unmatched_peak").to_numpy(),
                    negative & role.eq("paired_pre_injection").to_numpy(),
                ],
            )
            if {"metadata", "periodogram"} & profile_branches(profile):
                weights[role.eq("paired_pre_injection").to_numpy()] = 0.0
        positive = weights > 0
        if positive.any():
            weights[positive] /= np.mean(weights[positive])
        out.loc[fit_mask, f"{name}_weight"] = weights.astype(np.float32)
    return out


def _balanced_source_sampler_weights(rows: pd.DataFrame) -> np.ndarray:
    human = rows["teacher_v2_role"].astype(str).eq("real_human_morphology").to_numpy()
    injection = ~human
    weights = np.zeros(len(rows), dtype=np.float64)
    for mask in (human, injection):
        if np.any(mask):
            weights[mask] = 0.5 / np.count_nonzero(mask)
    weights /= weights.sum()
    return weights


def _loader(
    dataset: Any,
    indices: np.ndarray,
    *,
    batch_size: int,
    workers: int,
    seed: int,
    balanced_sources: bool,
) -> Any:
    import torch
    from torch.utils.data import DataLoader, Subset, WeightedRandomSampler

    subset_indices = np.asarray(indices, dtype=int)
    subset = Subset(dataset, subset_indices.tolist())
    sampler = None
    if balanced_sources:
        source_weights = _balanced_source_sampler_weights(
            dataset.rows.iloc[subset_indices].reset_index(drop=True)
        )
        sampler = WeightedRandomSampler(
            torch.as_tensor(source_weights, dtype=torch.double),
            num_samples=len(subset_indices),
            replacement=True,
            generator=torch.Generator().manual_seed(seed),
        )
    kwargs: dict[str, Any] = {}
    if workers > 0:
        kwargs.update({"persistent_workers": True, "prefetch_factor": 2})
    return DataLoader(
        subset,
        batch_size=int(batch_size),
        sampler=sampler,
        shuffle=False,
        num_workers=max(0, int(workers)),
        pin_memory=True,
        collate_fn=collate_native_samples,
        **kwargs,
    )


def _to_device(batch: Mapping[str, Any], device: Any) -> dict[str, Any]:
    import torch

    return {
        key: value.to(device, non_blocking=True) if isinstance(value, torch.Tensor) else value
        for key, value in batch.items()
    }


def _evaluate(model: Any, loader: Any, *, device: Any) -> dict[str, Any]:
    import torch

    names = ("morphology", "preserve", "harmonic", "compact")
    collected: dict[str, list[Any]] = {"review_id": []}
    for name in names:
        collected[f"{name}_logits"] = []
        collected[f"{name}_target"] = []
    model.eval()
    with torch.no_grad():
        for raw_batch in loader:
            batch = _to_device(raw_batch, device)
            with torch.autocast(
                device_type=device.type,
                dtype=torch.bfloat16,
                enabled=device.type == "cuda",
            ):
                output = model(batch)
            collected["review_id"].extend(raw_batch["review_id"])
            for name in names:
                collected[f"{name}_logits"].append(
                    output[f"{name}_logits"].float().cpu().numpy()
                )
                collected[f"{name}_target"].append(
                    raw_batch[f"{name}_target"].cpu().numpy()
                )
    return {
        "review_id": collected["review_id"],
        **{
            name: np.concatenate(values, axis=0) if values else np.empty((0, 0))
            for name, values in collected.items()
            if name != "review_id"
        },
    }


def compact_metrics(truth: np.ndarray, probability: np.ndarray) -> dict[str, Any]:
    truth = np.asarray(truth, dtype=int)
    probability = np.asarray(probability, dtype=float)
    active = truth >= 0
    truth = truth[active]
    probability = probability[active]
    if not len(truth):
        return {"n": 0, "average_precision": np.nan, "roc_auc": np.nan}
    predicted = probability >= 0.5
    positive = truth == 1
    order = np.argsort(-probability, kind="stable")
    ordered_positive = positive[order]
    cumulative_positive = np.cumsum(ordered_positive)
    precision_curve = cumulative_positive / np.arange(1, len(truth) + 1)
    average_precision = (
        float(np.sum(precision_curve * ordered_positive) / np.count_nonzero(positive))
        if np.any(positive)
        else np.nan
    )
    n_positive = int(np.count_nonzero(positive))
    n_negative = int(len(truth) - n_positive)
    ranks = pd.Series(probability).rank(method="average").to_numpy(dtype=float)
    roc_auc = (
        float(
            (ranks[positive].sum() - n_positive * (n_positive + 1) / 2)
            / (n_positive * n_negative)
        )
        if n_positive and n_negative
        else np.nan
    )
    true_positive = int(np.count_nonzero(predicted & positive))
    precision = true_positive / max(int(np.count_nonzero(predicted)), 1)
    recall = true_positive / max(int(np.count_nonzero(positive)), 1)
    return {
        "n": int(len(truth)),
        "n_positive": n_positive,
        "average_precision": average_precision,
        "roc_auc": roc_auc,
        "precision_at_0p5": float(precision),
        "recall_at_0p5": float(recall),
    }


def _development_score(
    rows: pd.DataFrame,
    evaluated: Mapping[str, Any],
) -> tuple[float, dict[str, Any]]:
    morphology_probability = _softmax(evaluated["morphology_logits"])
    morphology = classification_metrics(
        evaluated["morphology_target"], morphology_probability, classes=MORPHOLOGY_CLASSES
    )
    candidate_mask = ~rows["teacher_v2_role"].astype(str).eq("paired_pre_injection").to_numpy()
    compact_probability = _softmax(evaluated["compact_logits"])[:, 1]
    compact = compact_metrics(
        evaluated["compact_target"][candidate_mask], compact_probability[candidate_mask]
    )
    macro_f1 = float(morphology["macro_f1"])
    average_precision = float(compact["average_precision"])
    if not np.isfinite(macro_f1):
        macro_f1 = 0.0
    if not np.isfinite(average_precision):
        average_precision = 0.0
    score = 0.70 * average_precision + 0.30 * macro_f1
    return score, {
        "selection_score": score,
        "morphology": morphology,
        "compact_candidates_only": compact,
    }


def train_teacher_v2_fold(
    *,
    rows: pd.DataFrame,
    out_dir: Path,
    profile: str,
    fold: int,
    train_config: TeacherV2TrainConfig,
    workers: int,
    require_cuda: bool,
) -> dict[str, Any]:
    import torch

    seed = int(train_config.seed + 100 * fold)
    torch.manual_seed(seed)
    np.random.seed(seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("Teacher-v2 training requires CUDA")
    train_mask = rows["fixed_split"].eq("development") & rows["cv_fold"].ne(fold)
    validation_mask = rows["fixed_split"].eq("development") & rows["cv_fold"].eq(fold)
    weighted = attach_teacher_v2_fold_weights(
        rows,
        fit_mask=train_mask.to_numpy(),
        profile=profile,
        class_weight_cap=train_config.class_weight_cap,
    )
    metadata, normalization = build_teacher_v2_metadata_matrix(
        weighted, fit_mask=train_mask.to_numpy()
    )
    model_config = HarmonicModelConfig(metadata_dim=metadata.shape[1])
    dataset = HarmonicNativeDataset(
        weighted,
        native_h5=None,
        metadata=metadata,
        cache_size=2048,
        profile=profile,
    )
    train_loader = _loader(
        dataset,
        np.flatnonzero(train_mask.to_numpy()),
        batch_size=train_config.batch_size,
        workers=workers,
        seed=seed,
        balanced_sources=True,
    )
    validation_loader = _loader(
        dataset,
        np.flatnonzero(validation_mask.to_numpy()),
        batch_size=train_config.batch_size,
        workers=workers,
        seed=seed,
        balanced_sources=False,
    )
    model = build_teacher_v2_cnn(model_config, profile=profile).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=train_config.learning_rate,
        weight_decay=train_config.weight_decay,
    )
    history: list[dict[str, Any]] = []
    best_score = -np.inf
    best_state: dict[str, Any] | None = None
    best_epoch = -1
    stale = 0
    for epoch in range(1, int(train_config.epochs) + 1):
        model.train()
        total_loss = 0.0
        seen = 0
        for raw_batch in train_loader:
            batch = _to_device(raw_batch, device)
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(
                device_type=device.type,
                dtype=torch.bfloat16,
                enabled=device.type == "cuda",
            ):
                output = model(batch)
                loss, _ = teacher_v2_multitask_loss(output, batch, config=train_config)
            loss.backward()
            optimizer.step()
            size = len(raw_batch["review_id"])
            total_loss += float(loss.detach()) * size
            seen += size
        evaluated = _evaluate(model, validation_loader, device=device)
        validation_rows = weighted.loc[validation_mask].reset_index(drop=True)
        score, metrics = _development_score(validation_rows, evaluated)
        history.append(
            {
                "epoch": epoch,
                "train_loss": total_loss / max(seen, 1),
                "validation_selection_score": score,
                "validation_macro_f1": metrics["morphology"]["macro_f1"],
                "validation_compact_ap": metrics["compact_candidates_only"]["average_precision"],
            }
        )
        if np.isfinite(score) and score > best_score + 1.0e-6:
            best_score = score
            best_epoch = epoch
            best_state = {
                name: value.detach().cpu().clone() for name, value in model.state_dict().items()
            }
            stale = 0
        else:
            stale += 1
        print(
            f"[teacher-v2 {profile} fold={fold}] epoch={epoch} "
            f"loss={history[-1]['train_loss']:.5f} score={score:.4f} stale={stale}",
            flush=True,
        )
        if stale >= int(train_config.patience):
            break
    if best_state is None:
        raise RuntimeError(f"Teacher-v2 {profile} fold {fold} produced no finite checkpoint")
    model.load_state_dict(best_state)
    evaluated = _evaluate(model, validation_loader, device=device)
    validation_rows = weighted.loc[validation_mask].reset_index(drop=True)
    morphology_temperature = fit_temperature(
        evaluated["morphology_logits"], evaluated["morphology_target"]
    )
    compact_temperature = fit_temperature(
        evaluated["compact_logits"], evaluated["compact_target"]
    )
    morphology_probability = _softmax(
        evaluated["morphology_logits"], temperature=morphology_temperature
    )
    compact_probability = _softmax(
        evaluated["compact_logits"], temperature=compact_temperature
    )
    candidate_mask = ~validation_rows["teacher_v2_role"].astype(str).eq(
        "paired_pre_injection"
    ).to_numpy()
    metrics = {
        "morphology": classification_metrics(
            evaluated["morphology_target"], morphology_probability, classes=MORPHOLOGY_CLASSES
        ),
        "morphology_calibration": expected_calibration_error(
            evaluated["morphology_target"], morphology_probability
        ),
        "compact_candidates_only": compact_metrics(
            evaluated["compact_target"][candidate_mask],
            compact_probability[candidate_mask, 1],
        ),
        "compact_all_roles": compact_metrics(
            evaluated["compact_target"], compact_probability[:, 1]
        ),
        "paired_original_rejection_at_0p5": float(
            np.mean(
                compact_probability[
                    validation_rows["teacher_v2_role"].astype(str).eq(
                        "paired_pre_injection"
                    ),
                    1,
                ]
                < 0.5
            )
        )
        if validation_rows["teacher_v2_role"].astype(str).eq("paired_pre_injection").any()
        else np.nan,
    }
    fold_dir = out_dir / profile / f"fold_{fold}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = {
        "model_version": TEACHER_V2_MODEL_VERSION,
        "role_policy": TEACHER_V2_ROLE_POLICY,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "input_channel_contract": {
            name: list(channels) for name, channels in CHANNEL_CONTRACT.items()
        },
        "profile": profile,
        "fold": int(fold),
        "model_config": asdict(model_config),
        "train_config": asdict(train_config),
        "metadata_normalization": asdict(normalization),
        "morphology_temperature": float(morphology_temperature),
        "compact_temperature": float(compact_temperature),
        "best_epoch": int(best_epoch),
        "model_state_dict": best_state,
    }
    torch.save(checkpoint, fold_dir / "teacher_v2.pt")
    pd.DataFrame(history).to_csv(fold_dir / "history.csv", index=False)
    (fold_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    predictions = validation_rows.loc[
        :,
        [
            "review_id",
            "tic",
            "teacher_v2_role",
            "teacher_v2_partition",
            "fixed_split",
            "cv_fold",
            "is_injected_row",
            "injection_id",
            "is_injected_signal_peak",
            "human_label",
        ],
    ].copy()
    predictions["morphology_target"] = evaluated["morphology_target"]
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        predictions[f"p_{label}"] = morphology_probability[:, index]
    predictions["compact_target"] = evaluated["compact_target"]
    predictions["p_compact_transit"] = compact_probability[:, 1]
    predictions.to_parquet(
        fold_dir / "validation_predictions.parquet", compression="zstd", index=False
    )
    return {
        "profile": profile,
        "fold": int(fold),
        "best_epoch": int(best_epoch),
        "best_selection_score": float(best_score),
        "checkpoint": str(fold_dir / "teacher_v2.pt"),
        "metrics": metrics,
    }


def _validate_training_rows(rows: pd.DataFrame) -> dict[str, Any]:
    required = {
        "review_id",
        "tic",
        "teacher_v2_role",
        "teacher_v2_role_policy",
        "fixed_split",
        "cv_fold",
        "native_h5_path",
        "native_group_path",
        "morphology_target_index",
        "preserve_target_index",
        "harmonic_target_index",
        "compact_target_index",
    }
    missing = sorted(required - set(rows.columns))
    if missing:
        raise KeyError(f"Teacher-v2 training table is missing columns: {missing}")
    if rows["review_id"].fillna("").astype(str).duplicated().any():
        raise ValueError("Teacher-v2 training review IDs are not unique")
    if not rows["teacher_v2_role_policy"].astype(str).eq(TEACHER_V2_ROLE_POLICY).all():
        raise ValueError("Teacher-v2 training table has an unexpected role policy")
    development = set(rows.loc[rows["fixed_split"].eq("development"), "tic"])
    test = set(rows.loc[rows["fixed_split"].eq("test"), "tic"])
    if development & test:
        raise ValueError("Teacher-v2 development and test TICs overlap")
    dev_folds = pd.to_numeric(
        rows.loc[rows["fixed_split"].eq("development"), "cv_fold"], errors="coerce"
    )
    if dev_folds.isna().any() or not dev_folds.between(0, 4).all():
        raise ValueError("Teacher-v2 development rows require cv_fold in 0..4")
    verifications: dict[str, Any] = {}
    for value in sorted(rows["native_h5_path"].astype(str).unique()):
        verification = verify_raw_pair_contract(
            Path(value), require_errors=True, require_periodograms=True
        )
        if not verification["passed"]:
            raise RuntimeError(f"native input failed verification: {value}")
        verifications[value] = verification
    return {
        "n_rows": int(len(rows)),
        "n_unique_tics": int(pd.to_numeric(rows["tic"], errors="coerce").nunique()),
        "native_inputs": verifications,
    }


def run_teacher_v2_training(
    *,
    training_table: Path,
    out_dir: Path,
    profiles: Sequence[str] = DEFAULT_TEACHER_V2_PROFILES,
    train_config: TeacherV2TrainConfig = TeacherV2TrainConfig(),
    workers: int = 8,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Train development folds only; holdout and S57 remain unopened."""

    unknown = sorted(set(profiles) - set(MODEL_PROFILES))
    if unknown:
        raise ValueError(f"unknown Teacher-v2 profiles: {unknown}")
    rows = (
        pd.read_parquet(training_table)
        if Path(training_table).suffix.lower() == ".parquet"
        else pd.read_csv(training_table, low_memory=False)
    )
    verification = _validate_training_rows(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for profile in profiles:
        for fold in range(5):
            results.append(
                train_teacher_v2_fold(
                    rows=rows,
                    out_dir=out_dir,
                    profile=profile,
                    fold=fold,
                    train_config=train_config,
                    workers=workers,
                    require_cuda=require_cuda,
                )
            )
    profile_rows: list[dict[str, Any]] = []
    for profile in profiles:
        selected = [result for result in results if result["profile"] == profile]
        profile_rows.append(
            {
                "profile": profile,
                "fold_selection_score_mean": float(
                    np.mean([result["best_selection_score"] for result in selected])
                ),
                "fold_compact_ap_mean": float(
                    np.mean(
                        [
                            result["metrics"]["compact_candidates_only"]["average_precision"]
                            for result in selected
                        ]
                    )
                ),
                "fold_morphology_macro_f1_mean": float(
                    np.mean([result["metrics"]["morphology"]["macro_f1"] for result in selected])
                ),
            }
        )
    ranking = pd.DataFrame(profile_rows).sort_values(
        "fold_selection_score_mean", ascending=False, kind="stable"
    )
    ranking.to_csv(out_dir / "development_profile_diagnostics.csv", index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_version": TEACHER_V2_MODEL_VERSION,
        "role_policy": TEACHER_V2_ROLE_POLICY,
        "profiles": list(profiles),
        "train_config": asdict(train_config),
        "training_verification": verification,
        "development_diagnostics": ranking.to_dict("records"),
        "architecture_selection_frozen": False,
        "holdout_opened": False,
        "s57_opened": False,
    }
    (out_dir / "training_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    return summary


__all__ = [
    "DEFAULT_TEACHER_V2_PROFILES",
    "TEACHER_V2_METADATA_COLUMNS",
    "attach_teacher_v2_fold_weights",
    "build_teacher_v2_metadata_matrix",
    "compact_metrics",
    "run_teacher_v2_training",
    "train_teacher_v2_fold",
]
