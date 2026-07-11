"""Training and evaluation loop for the S56 post-adjudication harmonic CNN."""
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
    MODEL_VERSION,
    MORPHOLOGY_CLASSES,
    PRESERVE_CLASSES,
    HarmonicModelConfig,
    HarmonicTrainConfig,
    build_harmonic_cnn,
    multitask_loss,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    MetadataNormalization,
    attach_fold_training_weights,
    build_injection_pretraining_rows,
    build_metadata_matrix,
    collate_native_samples,
    prepare_harmonic_training_rows,
)
from twirl.vetting.harmonic_inputs import (
    CHANNEL_CONTRACT,
    RAW_PAIR_CONTRACT_VERSION,
    verify_raw_pair_contract,
)


DEFAULT_PROFILES: tuple[str, ...] = (
    "metadata_only",
    "single_period_native_fold",
    "seven_harmonic_shape",
    "shape_plus_raw_chronology",
    "shape_plus_periodogram_bls",
    "full_combined",
)


def classification_metrics(
    truth: np.ndarray,
    probability: np.ndarray,
    *,
    classes: Sequence[str],
) -> dict[str, Any]:
    truth = np.asarray(truth, dtype=int)
    probability = np.asarray(probability, dtype=float)
    active = truth >= 0
    truth = truth[active]
    probability = probability[active]
    if not len(truth):
        return {
            "n": 0,
            "accuracy": float("nan"),
            "balanced_accuracy": float("nan"),
            "macro_f1": float("nan"),
            "per_class": {},
            "confusion_matrix": np.zeros((len(classes), len(classes)), dtype=int).tolist(),
        }
    predicted = probability.argmax(axis=1)
    matrix = np.zeros((len(classes), len(classes)), dtype=int)
    for actual, estimate in zip(truth, predicted):
        matrix[int(actual), int(estimate)] += 1
    per_class: dict[str, Any] = {}
    recalls: list[float] = []
    f1_values: list[float] = []
    for index, label in enumerate(classes):
        true_positive = int(matrix[index, index])
        actual = int(matrix[index].sum())
        estimated = int(matrix[:, index].sum())
        recall = true_positive / actual if actual else float("nan")
        precision = true_positive / estimated if estimated else 0.0
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if actual and precision + recall > 0
            else (0.0 if actual else float("nan"))
        )
        if np.isfinite(recall):
            recalls.append(float(recall))
        if np.isfinite(f1):
            f1_values.append(float(f1))
        per_class[str(label)] = {
            "n": actual,
            "precision": float(precision),
            "recall": float(recall),
            "f1": float(f1),
        }
    return {
        "n": int(len(truth)),
        "accuracy": float(np.mean(predicted == truth)),
        "balanced_accuracy": float(np.mean(recalls)) if recalls else float("nan"),
        "macro_f1": float(np.mean(f1_values)) if f1_values else float("nan"),
        "per_class": per_class,
        "confusion_matrix": matrix.tolist(),
        "predicted_class_counts": {
            str(classes[index]): int(np.count_nonzero(predicted == index))
            for index in range(len(classes))
        },
    }


def expected_calibration_error(
    truth: np.ndarray,
    probability: np.ndarray,
    *,
    n_bins: int = 10,
) -> dict[str, Any]:
    truth = np.asarray(truth, dtype=int)
    probability = np.asarray(probability, dtype=float)
    active = truth >= 0
    truth = truth[active]
    probability = probability[active]
    if not len(truth):
        return {"ece": float("nan"), "bins": []}
    confidence = probability.max(axis=1)
    correct = probability.argmax(axis=1) == truth
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    bins: list[dict[str, Any]] = []
    for low, high in zip(edges[:-1], edges[1:]):
        mask = (confidence >= low) & (confidence < high if high < 1 else confidence <= high)
        n = int(mask.sum())
        if not n:
            bins.append({"low": float(low), "high": float(high), "n": 0})
            continue
        accuracy = float(np.mean(correct[mask]))
        mean_confidence = float(np.mean(confidence[mask]))
        ece += n / len(truth) * abs(accuracy - mean_confidence)
        bins.append(
            {
                "low": float(low),
                "high": float(high),
                "n": n,
                "accuracy": accuracy,
                "confidence": mean_confidence,
            }
        )
    return {"ece": float(ece), "bins": bins}


def injection_truth_human_audit(
    rows: pd.DataFrame,
    *,
    out_dir: Path,
) -> dict[str, Any]:
    """Write truth-versus-human recovery tables without exposing truth to inputs."""

    source_kind = rows.get("source_kind", pd.Series("", index=rows.index)).fillna("").astype(str)
    injected = rows.loc[source_kind.str.contains("inject", case=False)].copy()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if injected.empty:
        summary = {"n_injections": 0, "trend_tables": {}}
        (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return summary
    truth_period = pd.to_numeric(injected.get("truth_period_d"), errors="coerce")
    bls_period = pd.to_numeric(injected.get("period_d"), errors="coerce")
    ratio = truth_period / bls_period
    factors = np.asarray([0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0])
    error = np.abs(ratio.to_numpy(dtype=float)[:, None] / factors[None, :] - 1.0)
    nearest_index = np.argmin(np.where(np.isfinite(error), error, np.inf), axis=1)
    nearest_error = error[np.arange(len(injected)), nearest_index]
    injected["truth_to_bls_period_factor"] = factors[nearest_index]
    injected["truth_to_bls_period_rel_error"] = nearest_error
    injected["bls_truth_period_or_harmonic_match"] = np.isfinite(nearest_error) & (
        nearest_error <= 0.02
    )
    injected["human_visible_planet"] = injected["human_label"].fillna("").astype(str).eq(
        "planet_like"
    )
    audit_columns = [
        column
        for column in (
            "review_id",
            "injection_id",
            "tic",
            "human_label",
            "truth_period_d",
            "truth_radius_rearth",
            "truth_model_depth",
            "truth_depth",
            "tmag",
            "period_d",
            "sde_max",
            "adp_sml_sde",
            "recovery_status_DET_FLUX_ADP_SML",
            "topn_exact_recovered_DET_FLUX_ADP_SML",
            "topn_harmonic_match_DET_FLUX_ADP_SML",
            "truth_grid_cell_id",
            "boundary_period_bin_key",
            "boundary_radius_bin_key",
            "boundary_tmag_bin_key",
            "truth_to_bls_period_factor",
            "truth_to_bls_period_rel_error",
            "bls_truth_period_or_harmonic_match",
            "human_visible_planet",
        )
        if column in injected
    ]
    injected.loc[:, audit_columns].to_csv(out_dir / "injection_truth_human_rows.csv", index=False)
    pd.crosstab(
        injected["human_label"].fillna("").astype(str),
        injected["bls_truth_period_or_harmonic_match"],
    ).to_csv(out_dir / "human_label_by_bls_truth_period_match.csv")
    recovery_tables: dict[str, str] = {}
    for column in (
        "recovery_status_DET_FLUX_ADP_SML",
        "topn_exact_recovered_DET_FLUX_ADP_SML",
        "topn_harmonic_match_DET_FLUX_ADP_SML",
    ):
        if column not in injected:
            continue
        path = out_dir / f"human_label_by_{column}.csv"
        pd.crosstab(
            injected[column].fillna("").astype(str),
            injected["human_label"].fillna("").astype(str),
        ).to_csv(path)
        recovery_tables[column] = str(path)

    specifications: dict[str, tuple[Sequence[float], Sequence[str]]] = {
        "truth_period_d": (
            (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf),
            ("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10"),
        ),
        "truth_radius_rearth": (
            (0.0, 1.0, 2.0, 4.0, 8.0, 12.0, np.inf),
            ("<1", "1-2", "2-4", "4-8", "8-12", ">12"),
        ),
        "tmag": (
            (-np.inf, 16.0, 17.0, 18.0, 19.0, 20.0, np.inf),
            ("<16", "16-17", "17-18", "18-19", "19-20", ">20"),
        ),
        "truth_model_depth": (
            (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, np.inf),
            ("<1%", "1-3%", "3-10%", "10-30%", "30-100%", ">100%"),
        ),
        "adp_sml_sde": (
            (-np.inf, 5.0, 10.0, 20.0, 40.0, 80.0, np.inf),
            ("<5", "5-10", "10-20", "20-40", "40-80", ">80"),
        ),
    }
    trend_tables: dict[str, str] = {}
    for column, (bins, labels) in specifications.items():
        if column not in injected:
            continue
        work = injected.copy()
        work["parameter_bin"] = pd.cut(
            pd.to_numeric(work[column], errors="coerce"),
            bins=bins,
            labels=labels,
        )
        table = (
            work.groupby("parameter_bin", observed=False)
            .agg(
                n=("human_visible_planet", "size"),
                n_human_visible=("human_visible_planet", "sum"),
                human_visible_fraction=("human_visible_planet", "mean"),
                bls_truth_match_fraction=("bls_truth_period_or_harmonic_match", "mean"),
            )
            .reset_index()
        )
        path = out_dir / f"injection_recovery_by_{column}.csv"
        table.to_csv(path, index=False)
        trend_tables[column] = str(path)
    boundary_columns = (
        "boundary_period_bin_key",
        "boundary_radius_bin_key",
        "boundary_tmag_bin_key",
    )
    if all(column in injected for column in boundary_columns):
        injected["recovery50_cell"] = injected.loc[:, boundary_columns].fillna("").astype(str).agg(
            "|".join, axis=1
        )
        cell_column = "recovery50_cell"
    else:
        cell_column = (
            "truth_grid_cell_id"
            if "truth_grid_cell_id" in injected and injected["truth_grid_cell_id"].notna().any()
            else ""
        )
    if cell_column:
        cell_table = (
            injected.groupby(cell_column, dropna=False)
            .agg(
                n=("human_visible_planet", "size"),
                n_human_visible=("human_visible_planet", "sum"),
                human_visible_fraction=("human_visible_planet", "mean"),
                bls_truth_match_fraction=("bls_truth_period_or_harmonic_match", "mean"),
            )
            .reset_index()
        )
        cell_path = out_dir / "injection_recovery_by_recovery50_cell.csv"
        cell_table.to_csv(cell_path, index=False)
        trend_tables["recovery50_cell"] = str(cell_path)
    summary = {
        "n_injections": int(len(injected)),
        "human_label_counts": {
            str(key): int(value)
            for key, value in injected["human_label"].fillna("").astype(str).value_counts().items()
        },
        "n_human_visible_planet": int(injected["human_visible_planet"].sum()),
        "n_bls_truth_period_or_harmonic_match": int(
            injected["bls_truth_period_or_harmonic_match"].sum()
        ),
        "recovery_tables": recovery_tables,
        "trend_tables": trend_tables,
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    return summary


def _softmax(values: np.ndarray, *, temperature: float = 1.0) -> np.ndarray:
    scaled = np.asarray(values, dtype=np.float64) / max(float(temperature), 1.0e-6)
    scaled -= np.nanmax(scaled, axis=1, keepdims=True)
    exponential = np.exp(scaled)
    return exponential / np.maximum(exponential.sum(axis=1, keepdims=True), 1.0e-12)


def _cross_entropy_from_probability(truth: np.ndarray, probability: np.ndarray) -> float:
    truth = np.asarray(truth, dtype=int)
    probability = np.asarray(probability, dtype=float)
    active = truth >= 0
    if not np.any(active):
        return float("nan")
    selected = probability[active, truth[active]]
    return float(-np.mean(np.log(np.clip(selected, 1.0e-12, 1.0))))


def fit_temperature(logits: np.ndarray, truth: np.ndarray) -> float:
    """Fit one positive morphology temperature on validation logits."""

    import torch
    from torch.nn import functional as functional

    truth = np.asarray(truth, dtype=int)
    active = truth >= 0
    if np.count_nonzero(active) < 2:
        return 1.0
    values = torch.as_tensor(np.asarray(logits)[active], dtype=torch.float64)
    targets = torch.as_tensor(truth[active], dtype=torch.long)
    log_temperature = torch.nn.Parameter(torch.zeros((), dtype=torch.float64))
    optimizer = torch.optim.LBFGS([log_temperature], lr=0.05, max_iter=100)

    def closure() -> Any:
        optimizer.zero_grad()
        temperature = log_temperature.exp().clamp(0.05, 20.0)
        loss = functional.cross_entropy(values / temperature, targets)
        loss.backward()
        return loss

    try:
        optimizer.step(closure)
    except Exception:
        return 1.0
    return float(log_temperature.detach().exp().clamp(0.05, 20.0))


def _loader(
    dataset: Any,
    indices: np.ndarray,
    *,
    batch_size: int,
    shuffle: bool,
    workers: int,
    seed: int,
) -> Any:
    import torch
    from torch.utils.data import DataLoader, Subset

    kwargs: dict[str, Any] = {}
    if workers > 0:
        kwargs.update({"persistent_workers": True, "prefetch_factor": 2})
    return DataLoader(
        Subset(dataset, [int(value) for value in indices]),
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=workers,
        pin_memory=True,
        collate_fn=collate_native_samples,
        generator=torch.Generator().manual_seed(seed),
        **kwargs,
    )


def _to_device(batch: Mapping[str, Any], device: Any) -> dict[str, Any]:
    import torch

    return {
        key: value.to(device, non_blocking=True) if isinstance(value, torch.Tensor) else value
        for key, value in batch.items()
    }


def _evaluate_loader(model: Any, loader: Any, *, device: Any) -> dict[str, Any]:
    import torch

    model.eval()
    collected: dict[str, list[Any]] = {
        "review_id": [],
        "morphology_logits": [],
        "preserve_logits": [],
        "harmonic_logits": [],
        "morphology_target": [],
        "preserve_target": [],
        "harmonic_target": [],
    }
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
            for name in ("morphology_logits", "preserve_logits", "harmonic_logits"):
                collected[name].append(output[name].float().cpu().numpy())
            for name in ("morphology_target", "preserve_target", "harmonic_target"):
                collected[name].append(raw_batch[name].cpu().numpy())
    return {
        "review_id": collected["review_id"],
        **{
            name: np.concatenate(values, axis=0) if values else np.empty((0, 0))
            for name, values in collected.items()
            if name != "review_id"
        },
    }


def _train_pretraining_encoder(
    *,
    rows: pd.DataFrame,
    native_h5: Path,
    metadata: np.ndarray,
    model_config: HarmonicModelConfig,
    train_config: HarmonicTrainConfig,
    device: Any,
    workers: int,
    epochs: int,
    seed: int,
) -> dict[str, Any]:
    import torch
    from torch import nn

    pretrain_rows = build_injection_pretraining_rows(rows, seed=seed)
    pretrain_metadata, _ = build_metadata_matrix(
        pretrain_rows,
        fit_mask=np.ones(len(pretrain_rows), dtype=bool),
        normalization=MetadataNormalization(
            columns=tuple(), center=tuple(), scale=tuple()
        ) if metadata.shape[1] == 0 else None,
    )
    if metadata.shape[1] and pretrain_metadata.shape[1] != metadata.shape[1]:
        # The pretraining profile does not consume metadata, but model tensor
        # dimensions must still match the fine-tuning checkpoint.
        pretrain_metadata = np.zeros((len(pretrain_rows), metadata.shape[1]), dtype=np.float32)
    dataset = HarmonicNativeDataset(
        pretrain_rows,
        native_h5=native_h5,
        metadata=pretrain_metadata,
        profile="shape_plus_raw_chronology",
    )
    loader = _loader(
        dataset,
        np.arange(len(dataset)),
        batch_size=train_config.batch_size,
        shuffle=True,
        workers=workers,
        seed=seed,
    )
    model = build_harmonic_cnn(model_config, profile="shape_plus_raw_chronology").to(device)
    head = nn.Linear(2 * model_config.embedding_dim, 2).to(device)
    optimizer = torch.optim.AdamW(
        list(model.parameters()) + list(head.parameters()),
        lr=train_config.learning_rate,
        weight_decay=train_config.weight_decay,
    )
    counts = pretrain_rows["pretrain_target"].value_counts()
    class_weight = torch.tensor(
        [len(pretrain_rows) / max(int(counts.get(index, 0)), 1) for index in range(2)],
        dtype=torch.float32,
        device=device,
    )
    history: list[dict[str, float]] = []
    for epoch in range(1, int(epochs) + 1):
        model.train()
        head.train()
        total = 0.0
        seen = 0
        for raw_batch in loader:
            batch = _to_device(raw_batch, device)
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(
                device_type=device.type,
                dtype=torch.bfloat16,
                enabled=device.type == "cuda",
            ):
                embedding = model(batch)["embedding"]
                logits = head(embedding)
                loss = nn.functional.cross_entropy(
                    logits,
                    batch["pretrain_target"],
                    weight=class_weight,
                )
            loss.backward()
            optimizer.step()
            total += float(loss.detach()) * len(raw_batch["review_id"])
            seen += len(raw_batch["review_id"])
        history.append({"epoch": float(epoch), "loss": total / max(seen, 1)})
    return {
        "state_dict": {name: value.detach().cpu() for name, value in model.state_dict().items()},
        "history": history,
        "n_rows": len(pretrain_rows),
        "n_visible_injected": int(pretrain_rows["pretrain_target"].eq(1).sum()),
        "model_config": asdict(model_config),
        "seed": seed,
    }


def _train_one_fold(
    *,
    rows: pd.DataFrame,
    native_h5: Path,
    out_dir: Path,
    profile: str,
    fold: int,
    train_config: HarmonicTrainConfig,
    workers: int,
    pretrain_epochs: int,
    require_cuda: bool,
) -> dict[str, Any]:
    import torch

    seed = train_config.seed + 100 * fold
    torch.manual_seed(seed)
    np.random.seed(seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required but is unavailable")
    train_mask = rows["fixed_split"].eq("development") & rows["cv_fold"].ne(fold)
    validation_mask = rows["fixed_split"].eq("development") & rows["cv_fold"].eq(fold)
    fold_rows = attach_fold_training_weights(rows, fit_mask=train_mask.to_numpy())
    metadata, normalization = build_metadata_matrix(rows, fit_mask=train_mask.to_numpy())
    model_config = HarmonicModelConfig(metadata_dim=metadata.shape[1])
    pretrain_path = Path(out_dir) / "encoder_pretraining" / f"fold_{fold}.pt"
    if pretrain_path.exists():
        pretrain = torch.load(pretrain_path, map_location="cpu", weights_only=False)
        if pretrain.get("model_config") != asdict(model_config) or int(pretrain.get("seed", -1)) != seed:
            raise RuntimeError(f"stale encoder pretraining checkpoint: {pretrain_path}")
    else:
        pretrain = _train_pretraining_encoder(
            rows=fold_rows.loc[train_mask].copy(),
            native_h5=native_h5,
            metadata=metadata[train_mask],
            model_config=model_config,
            train_config=train_config,
            device=device,
            workers=workers,
            epochs=pretrain_epochs,
            seed=seed,
        )
        pretrain_path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(pretrain, pretrain_path)
    dataset = HarmonicNativeDataset(
        fold_rows,
        native_h5=native_h5,
        metadata=metadata,
        profile=profile,
    )
    train_loader = _loader(
        dataset,
        np.flatnonzero(train_mask.to_numpy()),
        batch_size=train_config.batch_size,
        shuffle=True,
        workers=workers,
        seed=seed,
    )
    validation_loader = _loader(
        dataset,
        np.flatnonzero(validation_mask.to_numpy()),
        batch_size=train_config.batch_size,
        shuffle=False,
        workers=workers,
        seed=seed,
    )
    model = build_harmonic_cnn(model_config, profile=profile).to(device)
    model.load_state_dict(pretrain["state_dict"], strict=True)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=train_config.learning_rate,
        weight_decay=train_config.weight_decay,
    )
    best_state: dict[str, Any] | None = None
    best_macro_f1 = -np.inf
    best_epoch = 0
    stale = 0
    history: list[dict[str, Any]] = []
    for epoch in range(1, train_config.epochs + 1):
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
                outputs = model(batch)
                loss, _ = multitask_loss(outputs, batch, config=train_config)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.detach()) * len(raw_batch["review_id"])
            seen += len(raw_batch["review_id"])
        validation = _evaluate_loader(model, validation_loader, device=device)
        probability = _softmax(validation["morphology_logits"])
        validation_loss = _cross_entropy_from_probability(
            validation["morphology_target"], probability
        )
        metrics = classification_metrics(
            validation["morphology_target"], probability, classes=MORPHOLOGY_CLASSES
        )
        macro_f1 = float(metrics["macro_f1"])
        history.append(
            {
                "epoch": epoch,
                "train_loss": total_loss / max(seen, 1),
                "validation_morphology_loss": validation_loss,
                "validation_macro_f1": macro_f1,
                "validation_balanced_accuracy": metrics["balanced_accuracy"],
            }
        )
        if np.isfinite(macro_f1) and macro_f1 > best_macro_f1 + 1.0e-6:
            best_macro_f1 = macro_f1
            best_epoch = epoch
            best_state = {name: value.detach().cpu().clone() for name, value in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
        print(
            f"[{profile} fold={fold}] epoch={epoch} loss={history[-1]['train_loss']:.5f} "
            f"val_macro_f1={macro_f1:.4f} stale={stale}",
            flush=True,
        )
        if stale >= train_config.patience:
            break
    if best_state is None:
        raise RuntimeError(f"{profile} fold {fold} never produced finite validation metrics")
    model.load_state_dict(best_state)
    validation = _evaluate_loader(model, validation_loader, device=device)
    temperature = fit_temperature(
        validation["morphology_logits"], validation["morphology_target"]
    )
    morphology_probability = _softmax(
        validation["morphology_logits"], temperature=temperature
    )
    preserve_probability = _softmax(validation["preserve_logits"])
    harmonic_probability = _softmax(validation["harmonic_logits"])
    validation_metrics = {
        "morphology": classification_metrics(
            validation["morphology_target"], morphology_probability, classes=MORPHOLOGY_CLASSES
        ),
        "morphology_by_source": _subset_metrics(
            fold_rows.loc[validation_mask].reset_index(drop=True),
            validation["morphology_target"],
            morphology_probability,
        ),
        "preserve": classification_metrics(
            validation["preserve_target"], preserve_probability, classes=PRESERVE_CLASSES
        ),
        "harmonic": classification_metrics(
            validation["harmonic_target"], harmonic_probability, classes=HARMONIC_CLASSES
        ),
        "calibration": expected_calibration_error(
            validation["morphology_target"], morphology_probability
        ),
    }
    fold_dir = Path(out_dir) / profile / f"fold_{fold}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = {
        "model_version": MODEL_VERSION,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "input_channel_contract": {
            name: list(channels) for name, channels in CHANNEL_CONTRACT.items()
        },
        "profile": profile,
        "fold": fold,
        "model_config": asdict(model_config),
        "train_config": asdict(train_config),
        "metadata_normalization": asdict(normalization),
        "temperature": temperature,
        "best_epoch": best_epoch,
        "model_state_dict": best_state,
    }
    torch.save(checkpoint, fold_dir / "teacher.pt")
    (fold_dir / "metrics.json").write_text(
        json.dumps(validation_metrics, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    pd.DataFrame(history).to_csv(fold_dir / "history.csv", index=False)
    pd.DataFrame(pretrain["history"]).to_csv(fold_dir / "pretraining_history.csv", index=False)
    predictions = pd.DataFrame(
        {
            "review_id": validation["review_id"],
            "morphology_target": validation["morphology_target"],
            "morphology_prediction": morphology_probability.argmax(axis=1),
            **{
                f"p_{label}": morphology_probability[:, index]
                for index, label in enumerate(MORPHOLOGY_CLASSES)
            },
        }
    )
    predictions.to_csv(fold_dir / "validation_predictions.csv", index=False)
    return {
        "profile": profile,
        "fold": fold,
        "best_epoch": best_epoch,
        "temperature": temperature,
        "metrics": validation_metrics,
        "checkpoint": str(fold_dir / "teacher.pt"),
        "pretraining_n": pretrain["n_rows"],
    }


def _subset_metrics(
    rows: pd.DataFrame,
    truth: np.ndarray,
    probability: np.ndarray,
) -> dict[str, Any]:
    injected = rows["is_injected_row"]
    if injected.dtype != bool:
        injected = injected.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})
    output: dict[str, Any] = {
        "all": classification_metrics(truth, probability, classes=MORPHOLOGY_CLASSES),
        "real": classification_metrics(
            truth[~injected.to_numpy()], probability[~injected.to_numpy()], classes=MORPHOLOGY_CLASSES
        ),
        "injected": classification_metrics(
            truth[injected.to_numpy()], probability[injected.to_numpy()], classes=MORPHOLOGY_CLASSES
        ),
    }
    human = rows["human_label"].fillna("").astype(str)
    for name, label in (
        ("other_flat", "uncertain"),
        ("other_systematic", "instrumental_or_systematic"),
    ):
        mask = human.eq(label).to_numpy()
        output[name] = classification_metrics(
            truth[mask], probability[mask], classes=MORPHOLOGY_CLASSES
        )
    return output


def _calibration_by_source(
    rows: pd.DataFrame,
    truth: np.ndarray,
    probability: np.ndarray,
) -> dict[str, Any]:
    injected = rows["is_injected_row"]
    if injected.dtype != bool:
        injected = injected.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})
    mask = injected.to_numpy()
    return {
        "all": expected_calibration_error(truth, probability),
        "real": expected_calibration_error(truth[~mask], probability[~mask]),
        "injected": expected_calibration_error(truth[mask], probability[mask]),
    }


def run_harmonic_teacher_training(
    *,
    training_table: Path,
    native_h5: Path,
    out_dir: Path,
    profiles: Sequence[str] = DEFAULT_PROFILES,
    train_config: HarmonicTrainConfig = HarmonicTrainConfig(),
    workers: int = 8,
    pretrain_epochs: int = 20,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Train all development folds, select once, then open the fixed test set."""

    import torch

    unknown = sorted(set(profiles) - set(MODEL_PROFILES))
    if unknown:
        raise ValueError(f"unknown profiles: {unknown}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(training_table, low_memory=False)
    injection_audit = injection_truth_human_audit(
        source,
        out_dir=out_dir / "injection_truth_human_audit",
    )
    verification = verify_raw_pair_contract(
        native_h5,
        require_errors=True,
        require_periodograms=True,
    )
    if not verification["passed"]:
        raise RuntimeError(f"native input contract failed: {verification['failures'][:10]}")
    rows = prepare_harmonic_training_rows(source, seed=train_config.seed)
    rows.to_csv(out_dir / "training_rows_with_fixed_splits.csv", index=False)
    fold_results: list[dict[str, Any]] = []
    for profile in profiles:
        for fold in range(5):
            fold_results.append(
                _train_one_fold(
                    rows=rows,
                    native_h5=native_h5,
                    out_dir=out_dir,
                    profile=profile,
                    fold=fold,
                    train_config=train_config,
                    workers=workers,
                    pretrain_epochs=pretrain_epochs,
                    require_cuda=require_cuda,
                )
            )
    ranking_rows: list[dict[str, Any]] = []
    for profile in profiles:
        selected = [result for result in fold_results if result["profile"] == profile]
        prediction_parts = [
            pd.read_csv(out_dir / profile / f"fold_{fold}" / "validation_predictions.csv")
            for fold in range(5)
        ]
        predictions = pd.concat(prediction_parts, ignore_index=True)
        if predictions["review_id"].duplicated().any():
            raise RuntimeError(f"duplicate development predictions for profile {profile}")
        development_rows = predictions.loc[:, ["review_id"]].merge(
            rows,
            on="review_id",
            how="left",
            validate="one_to_one",
        )
        probability = predictions.loc[
            :, [f"p_{label}" for label in MORPHOLOGY_CLASSES]
        ].to_numpy(dtype=float)
        truth = predictions["morphology_target"].to_numpy(dtype=int)
        source_metrics = _subset_metrics(development_rows, truth, probability)
        calibration = expected_calibration_error(truth, probability)
        real = source_metrics["real"]["per_class"]
        planet_recall = float(real.get("planet_like", {}).get("recall", 0.0))
        eb_recall = float(real.get("eclipse_contact", {}).get("recall", 0.0))
        variable_recall = float(real.get("smooth_variable", {}).get("recall", 0.0))
        other_recall = float(real.get("other", {}).get("recall", 0.0))
        macro_f1 = float(source_metrics["all"]["macro_f1"])
        balanced_accuracy = float(source_metrics["all"]["balanced_accuracy"])
        ece = float(calibration["ece"])
        selection_score = (
            0.30 * macro_f1
            + 0.20 * balanced_accuracy
            + 0.15 * planet_recall
            + 0.10 * eb_recall
            + 0.10 * variable_recall
            + 0.05 * other_recall
            + 0.10 * (1.0 - min(max(ece, 0.0), 1.0))
        )
        ranking_rows.append(
            {
                "profile": profile,
                "development_selection_score": selection_score,
                "validation_macro_f1": macro_f1,
                "validation_balanced_accuracy": balanced_accuracy,
                "validation_real_planet_recall": planet_recall,
                "validation_real_eb_recall": eb_recall,
                "validation_real_variable_recall": variable_recall,
                "validation_real_other_recall": other_recall,
                "validation_ece": ece,
                "fold_macro_f1_mean": float(
                    np.mean(
                        [result["metrics"]["morphology"]["macro_f1"] for result in selected]
                    )
                ),
            }
        )
    ranking = pd.DataFrame(ranking_rows).sort_values(
        ["development_selection_score", "validation_macro_f1"],
        ascending=False,
        kind="stable",
    )
    ranking.to_csv(out_dir / "development_profile_ranking.csv", index=False)
    selected_profile = str(ranking.iloc[0]["profile"])

    test_mask = rows["fixed_split"].eq("test").to_numpy()
    test_rows = rows.loc[test_mask].reset_index(drop=True)
    ensemble_probabilities: list[np.ndarray] = []
    ensemble_preserve: list[np.ndarray] = []
    ensemble_harmonic: list[np.ndarray] = []
    truth: np.ndarray | None = None
    preserve_truth: np.ndarray | None = None
    harmonic_truth: np.ndarray | None = None
    review_ids: list[str] | None = None
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    for fold in range(5):
        checkpoint_path = out_dir / selected_profile / f"fold_{fold}" / "teacher.pt"
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
        norm_dict = checkpoint["metadata_normalization"]
        normalization = MetadataNormalization(
            columns=tuple(norm_dict["columns"]),
            center=tuple(norm_dict["center"]),
            scale=tuple(norm_dict["scale"]),
        )
        metadata, _ = build_metadata_matrix(
            rows,
            fit_mask=np.zeros(len(rows), dtype=bool),
            normalization=normalization,
        )
        dataset = HarmonicNativeDataset(
            rows,
            native_h5=native_h5,
            metadata=metadata,
            profile=selected_profile,
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
            HarmonicModelConfig(**checkpoint["model_config"]), profile=selected_profile
        ).to(device)
        model.load_state_dict(checkpoint["model_state_dict"])
        scored = _evaluate_loader(model, loader, device=device)
        ensemble_probabilities.append(
            _softmax(scored["morphology_logits"], temperature=float(checkpoint["temperature"]))
        )
        ensemble_preserve.append(_softmax(scored["preserve_logits"]))
        ensemble_harmonic.append(_softmax(scored["harmonic_logits"]))
        if truth is None:
            truth = scored["morphology_target"]
            preserve_truth = scored["preserve_target"]
            harmonic_truth = scored["harmonic_target"]
            review_ids = scored["review_id"]
        elif review_ids != scored["review_id"]:
            raise RuntimeError("test score ordering changed between ensemble members")
    morphology_probability = np.mean(ensemble_probabilities, axis=0)
    preserve_probability = np.mean(ensemble_preserve, axis=0)
    harmonic_probability = np.mean(ensemble_harmonic, axis=0)
    assert truth is not None and preserve_truth is not None and harmonic_truth is not None
    test_metrics = {
        "morphology_by_source": _subset_metrics(test_rows, truth, morphology_probability),
        "preserve": classification_metrics(
            preserve_truth, preserve_probability, classes=PRESERVE_CLASSES
        ),
        "harmonic": classification_metrics(
            harmonic_truth, harmonic_probability, classes=HARMONIC_CLASSES
        ),
        "calibration": _calibration_by_source(test_rows, truth, morphology_probability),
    }
    broad_mask = test_rows["human_label"].fillna("").astype(str).eq("wide_transit_like").to_numpy()
    test_metrics["broad_preserve_audit"] = {
        "n": int(broad_mask.sum()),
        "mean_p_preserve": (
            float(np.mean(preserve_probability[broad_mask, 1]))
            if np.any(broad_mask)
            else float("nan")
        ),
        "preserve_recall_at_0p5": (
            float(np.mean(preserve_probability[broad_mask, 1] >= 0.5))
            if np.any(broad_mask)
            else float("nan")
        ),
    }
    active_harmonic = harmonic_truth >= 0
    if np.any(active_harmonic):
        counts = np.bincount(harmonic_truth[active_harmonic], minlength=len(HARMONIC_CLASSES))
        test_metrics["harmonic_majority_baseline_accuracy"] = float(counts.max() / counts.sum())
    else:
        test_metrics["harmonic_majority_baseline_accuracy"] = float("nan")
    prediction_table = test_rows.loc[:, ["review_id", "tic", "human_label", "is_injected_row"]].copy()
    prediction_table["morphology_target_index"] = truth
    prediction_table["morphology_prediction_index"] = morphology_probability.argmax(axis=1)
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        prediction_table[f"p_{label}"] = morphology_probability[:, index]
    prediction_table["p_preserve"] = preserve_probability[:, 1]
    prediction_table["predicted_period_factor"] = np.asarray(
        [0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0]
    )[harmonic_probability.argmax(axis=1)]
    prediction_table["predicted_period_d"] = (
        pd.to_numeric(test_rows["period_d"], errors="coerce").to_numpy()
        * prediction_table["predicted_period_factor"].to_numpy()
    )
    prediction_table.to_csv(out_dir / "fixed_test_predictions.csv", index=False)

    ranking_by_profile = ranking.set_index("profile")
    combined = ranking_by_profile.loc["full_combined"] if "full_combined" in ranking_by_profile.index else None
    metadata = ranking_by_profile.loc["metadata_only"] if "metadata_only" in ranking_by_profile.index else None
    real_metrics = test_metrics["morphology_by_source"]["real"]
    harmonic_accuracy = test_metrics["harmonic"]["accuracy"]
    majority = test_metrics["harmonic_majority_baseline_accuracy"]
    promotion = {
        "test_balanced_accuracy": bool(real_metrics["balanced_accuracy"] >= 0.75),
        "real_planet_recall": bool(
            real_metrics["per_class"].get("planet_like", {}).get("recall", 0.0) >= 0.80
        ),
        "eb_recall": bool(
            real_metrics["per_class"].get("eclipse_contact", {}).get("recall", 0.0) >= 0.60
        ),
        "variable_recall": bool(
            real_metrics["per_class"].get("smooth_variable", {}).get("recall", 0.0) >= 0.60
        ),
        "other_recall": bool(
            real_metrics["per_class"].get("other", {}).get("recall", 0.0) >= 0.60
        ),
        "ece": bool(test_metrics["calibration"]["all"]["ece"] <= 0.10),
        "harmonic_above_majority": bool(
            np.isfinite(harmonic_accuracy)
            and np.isfinite(majority)
            and harmonic_accuracy >= majority + 0.10
        ),
    }
    acceptance = {
        "nondegenerate_all_classes": all(
            test_metrics["morphology_by_source"]["all"]["predicted_class_counts"].get(label, 0) > 0
            for label in MORPHOLOGY_CLASSES
        ),
        "finite_calibration": bool(np.isfinite(test_metrics["calibration"]["all"]["ece"])),
        "combined_exceeds_metadata_macro_f1": bool(
            combined is not None
            and metadata is not None
            and combined["validation_macro_f1"] > metadata["validation_macro_f1"]
        ),
        "promotion_gates": promotion,
        "promotion_passed": all(promotion.values()),
    }
    summary = {
        "model_version": MODEL_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_profile": selected_profile,
        "profiles": list(profiles),
        "native_contract_verification": verification,
        "input_channel_contract": {
            name: list(channels) for name, channels in CHANNEL_CONTRACT.items()
        },
        "injection_truth_human_audit": injection_audit,
        "n_training_rows": len(rows),
        "n_fixed_test_rows": len(test_rows),
        "development_ranking": ranking.to_dict("records"),
        "development_selection_formula": (
            "0.30*macro_f1 + 0.20*balanced_accuracy + 0.15*real_planet_recall + "
            "0.10*real_eb_recall + 0.10*real_variable_recall + "
            "0.05*real_other_recall + 0.10*(1-ECE)"
        ),
        "test_metrics": test_metrics,
        "acceptance": acceptance,
        "student_training_blocked": True,
        "student_block_reason": "requires at least 50 unique real Planet-like examples",
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    return summary


__all__ = [
    "DEFAULT_PROFILES",
    "classification_metrics",
    "expected_calibration_error",
    "fit_temperature",
    "injection_truth_human_audit",
    "run_harmonic_teacher_training",
]
