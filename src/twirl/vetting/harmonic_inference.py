"""Ensemble inference for the S56-trained harmonic-CNN teacher."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    MODEL_VERSION,
    MORPHOLOGY_CLASSES,
    HarmonicModelConfig,
    build_harmonic_cnn,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    MetadataNormalization,
    build_metadata_matrix,
    collate_native_samples,
)
from twirl.vetting.harmonic_inputs import (
    A2V1_TEACHER_INPUT_CONTRACT,
    RAW_PAIR_CONTRACT_VERSION,
    native_group_path,
)
from twirl.vetting.recovery50_teacher import leakage_columns


SELECTED_TEACHER_PROFILE = "shape_plus_periodogram_bls"
TEACHER_TRAINING_SECTOR = 56
TEACHER_SCORE_POLICY = "active_learning_priority_not_scientific_classification"
HARMONIC_DISPLAY_LABELS: tuple[str, ...] = (
    "P_over_4",
    "P_over_3",
    "P_over_2",
    "P",
    "2P",
    "3P",
    "4P",
)
HARMONIC_FACTORS: tuple[float, ...] = (0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def rank_planet_enrichment(scored: pd.DataFrame) -> pd.DataFrame:
    """Return the highest-scoring Planet candidate per TIC with stable ranks."""

    required = {"tic", "p_planet_like", "p_preserve"}
    missing = sorted(required - set(scored.columns))
    if missing:
        raise KeyError(
            f"teacher score table is missing Planet-ranking columns: {missing}"
        )
    work = scored.copy()
    work["tic"] = pd.to_numeric(work["tic"], errors="coerce")
    work["p_planet_like"] = pd.to_numeric(work["p_planet_like"], errors="coerce")
    work["p_preserve"] = pd.to_numeric(work["p_preserve"], errors="coerce")
    if work[["tic", "p_planet_like", "p_preserve"]].isna().any().any():
        raise ValueError("teacher score table has invalid Planet-ranking values")
    order = ["p_planet_like", "p_preserve"]
    ascending = [False, False]
    if "std_p_planet_like" in work:
        order.append("std_p_planet_like")
        ascending.append(True)
    if "sde_max" in work:
        order.append("sde_max")
        ascending.append(False)
    order.extend(["tic"])
    ascending.extend([True])
    ranked = work.sort_values(order, ascending=ascending, kind="stable")
    ranked = ranked.drop_duplicates("tic", keep="first").reset_index(drop=True)
    ranked.insert(0, "planet_rank", np.arange(1, len(ranked) + 1, dtype=np.int64))
    ranked["ranking_policy"] = (
        "p_planet_like desc, p_preserve desc, ensemble uncertainty asc, SDE desc; "
        "one candidate per TIC"
    )
    return ranked


def _softmax(values: np.ndarray, *, temperature: float = 1.0) -> np.ndarray:
    scaled = np.asarray(values, dtype=np.float64) / max(float(temperature), 1.0e-6)
    scaled -= np.max(scaled, axis=1, keepdims=True)
    probability = np.exp(scaled)
    probability /= np.sum(probability, axis=1, keepdims=True)
    return probability


def _to_device(batch: dict[str, Any], device: Any) -> dict[str, Any]:
    import torch

    return {
        key: value.to(device, non_blocking=True)
        if isinstance(value, torch.Tensor)
        else value
        for key, value in batch.items()
    }


def _normalization(checkpoint: dict[str, Any]) -> MetadataNormalization:
    payload = checkpoint["metadata_normalization"]
    columns = tuple(str(value) for value in payload["columns"])
    leaks = leakage_columns(columns)
    if leaks:
        raise ValueError(f"teacher checkpoint metadata contains leakage: {leaks}")
    return MetadataNormalization(
        columns=columns,
        center=tuple(float(value) for value in payload["center"]),
        scale=tuple(float(value) for value in payload["scale"]),
    )


def prepare_inference_rows(
    candidates: pd.DataFrame,
    *,
    allow_injections: bool = False,
) -> pd.DataFrame:
    """Attach only structural fields required by ``HarmonicNativeDataset``.

    Injection groups are rejected by default so production candidate scoring
    cannot accidentally consume simulation products.  Validation workflows
    must opt in explicitly and still use the same leakage-audited model path.
    """

    rows = candidates.copy().reset_index(drop=True)
    required = {"tic", "period_d", "t0_bjd", "duration_min"}
    missing = sorted(required - set(rows.columns))
    if missing:
        raise KeyError(f"candidate table is missing columns: {missing}")
    for column in required:
        rows[column] = pd.to_numeric(rows[column], errors="coerce")
    invalid = (
        rows["tic"].isna()
        | rows["period_d"].isna()
        | rows["t0_bjd"].isna()
        | rows["duration_min"].isna()
        | rows["period_d"].le(0)
        | rows["duration_min"].le(0)
    )
    if invalid.any():
        raise ValueError(
            f"candidate table has {int(invalid.sum())} invalid ephemerides"
        )
    rows["tic"] = rows["tic"].astype(np.int64)
    if "sector" in rows:
        sectors = (
            pd.to_numeric(rows["sector"], errors="coerce").dropna().astype(int).unique()
        )
        if len(sectors) != 1:
            raise ValueError(
                "teacher inference requires exactly one sector per score artifact; "
                f"got {sorted(sectors.tolist())}"
            )
        rows["sector"] = int(sectors[0])
    if "review_id" not in rows:
        rows["review_id"] = [f"teacher-score-{index}" for index in range(len(rows))]
    rows["review_id"] = rows["review_id"].fillna("").astype(str)
    if rows["review_id"].eq("").any() or rows["review_id"].duplicated().any():
        raise ValueError("candidate review_id values must be nonempty and unique")
    if "native_group_path" not in rows:
        rows["native_group_path"] = [
            native_group_path(row) for row in rows.to_dict("records")
        ]
    has_injections = rows["native_group_path"].astype(str).str.startswith("injections/")
    if has_injections.any() and not allow_injections:
        raise ValueError("production teacher scoring accepts real target groups only")
    rows["morphology_target_index"] = -1
    rows["preserve_target_index"] = -1
    rows["harmonic_target_index"] = -1
    rows["morphology_weight"] = np.float32(0.0)
    rows["preserve_weight"] = np.float32(0.0)
    rows["harmonic_weight"] = np.float32(0.0)
    rows["pretrain_target"] = -1
    return rows


def _load_ensemble(
    checkpoint_paths: Sequence[Path],
    *,
    profile: str,
    device: Any,
) -> tuple[list[Any], list[dict[str, Any]], list[MetadataNormalization]]:
    import torch

    if len(checkpoint_paths) != 5:
        raise ValueError(f"expected five fold checkpoints, got {len(checkpoint_paths)}")
    models: list[Any] = []
    checkpoints: list[dict[str, Any]] = []
    normalizations: list[MetadataNormalization] = []
    seen_folds: set[int] = set()
    for path in checkpoint_paths:
        checkpoint = torch.load(Path(path), map_location="cpu", weights_only=False)
        if checkpoint.get("model_version") != MODEL_VERSION:
            raise ValueError(
                f"{path} has model_version={checkpoint.get('model_version')!r}"
            )
        if checkpoint.get("input_contract_version") != RAW_PAIR_CONTRACT_VERSION:
            raise ValueError(f"{path} has the wrong native input contract")
        if checkpoint.get("profile") != profile:
            raise ValueError(
                f"{path} has profile={checkpoint.get('profile')!r}, expected {profile!r}"
            )
        fold = int(checkpoint.get("fold", -1))
        if fold in seen_folds or fold not in range(5):
            raise ValueError(f"invalid or duplicate fold {fold} in {path}")
        seen_folds.add(fold)
        model = build_harmonic_cnn(
            HarmonicModelConfig(**checkpoint["model_config"]), profile=profile
        ).to(device)
        model.load_state_dict(checkpoint["model_state_dict"], strict=True)
        model.eval()
        models.append(model)
        checkpoints.append(checkpoint)
        normalizations.append(_normalization(checkpoint))
    order = np.argsort([int(value["fold"]) for value in checkpoints])
    return (
        [models[int(index)] for index in order],
        [checkpoints[int(index)] for index in order],
        [normalizations[int(index)] for index in order],
    )


def score_harmonic_teacher_ensemble(
    *,
    candidates: pd.DataFrame,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
    profile: str = SELECTED_TEACHER_PROFILE,
    batch_size: int = 32,
    workers: int = 4,
    require_cuda: bool = True,
    allow_injections: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Score candidates once through the data path and five times through the CNN."""

    import h5py
    import torch
    from torch.utils.data import DataLoader

    rows = prepare_inference_rows(candidates, allow_injections=allow_injections)
    with h5py.File(Path(native_h5), "r") as h5:
        if str(h5.attrs.get("contract_version", "")) != RAW_PAIR_CONTRACT_VERSION:
            raise ValueError("native HDF5 has the wrong contract_version")
        required_groups = rows["native_group_path"].astype(str).drop_duplicates()
        missing_groups = [value for value in required_groups if value not in h5]
        if missing_groups:
            raise ValueError(
                f"native HDF5 is missing {len(missing_groups)} candidate groups; "
                f"first={missing_groups[:5]}"
            )
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("teacher ensemble scoring requires CUDA")
    models, checkpoints, normalizations = _load_ensemble(
        checkpoint_paths, profile=profile, device=device
    )
    metadata_matrices: list[np.ndarray] = []
    for normalization in normalizations:
        work = rows.copy()
        for column in normalization.columns:
            if column not in work:
                work[column] = np.nan
        metadata, _ = build_metadata_matrix(
            work,
            fit_mask=np.zeros(len(work), dtype=bool),
            normalization=normalization,
        )
        metadata_matrices.append(metadata)

    # Shape tensors are identical across ensemble members. Build each batch
    # once, then swap only the fold-specific normalized scalar metadata.
    dataset = HarmonicNativeDataset(
        rows,
        native_h5=Path(native_h5),
        metadata=metadata_matrices[0],
        cache_size=0,
        profile=profile,
    )
    loader = DataLoader(
        dataset,
        batch_size=int(batch_size),
        shuffle=False,
        num_workers=max(0, int(workers)),
        pin_memory=device.type == "cuda",
        persistent_workers=int(workers) > 0,
        collate_fn=collate_native_samples,
    )
    row_lookup = {
        value: index for index, value in enumerate(rows["review_id"].astype(str))
    }
    morphology_members: list[list[np.ndarray]] = [[] for _ in models]
    preserve_members: list[list[np.ndarray]] = [[] for _ in models]
    harmonic_members: list[list[np.ndarray]] = [[] for _ in models]
    output_order: list[str] = []
    with torch.no_grad():
        for batch_index, raw_batch in enumerate(loader, start=1):
            ids = [str(value) for value in raw_batch["review_id"]]
            indices = np.asarray([row_lookup[value] for value in ids], dtype=int)
            output_order.extend(ids)
            batch = _to_device(raw_batch, device)
            for member, (model, checkpoint, metadata) in enumerate(
                zip(models, checkpoints, metadata_matrices)
            ):
                batch["metadata"] = torch.as_tensor(
                    metadata[indices], dtype=torch.float32, device=device
                )
                with torch.autocast(
                    device_type=device.type,
                    dtype=torch.bfloat16,
                    enabled=device.type == "cuda",
                ):
                    output = model(batch)
                morphology_members[member].append(
                    _softmax(
                        output["morphology_logits"].float().cpu().numpy(),
                        temperature=float(checkpoint["temperature"]),
                    )
                )
                preserve_members[member].append(
                    _softmax(output["preserve_logits"].float().cpu().numpy())
                )
                harmonic_members[member].append(
                    _softmax(output["harmonic_logits"].float().cpu().numpy())
                )
            if batch_index % 100 == 0:
                print(
                    f"[teacher-inference] batches={batch_index:,} rows={len(output_order):,}/{len(rows):,}",
                    flush=True,
                )
    if output_order != rows["review_id"].astype(str).tolist():
        raise RuntimeError("teacher inference output ordering changed")

    morphology = np.stack(
        [np.concatenate(values, axis=0) for values in morphology_members], axis=0
    )
    preserve = np.stack(
        [np.concatenate(values, axis=0) for values in preserve_members], axis=0
    )
    harmonic = np.stack(
        [np.concatenate(values, axis=0) for values in harmonic_members], axis=0
    )
    mean_morphology = morphology.mean(axis=0)
    mean_preserve = preserve.mean(axis=0)
    mean_harmonic = harmonic.mean(axis=0)
    scored = candidates.copy().reset_index(drop=True)
    for class_index, label in enumerate(MORPHOLOGY_CLASSES):
        scored[f"p_{label}"] = mean_morphology[:, class_index]
        scored[f"std_p_{label}"] = morphology[:, :, class_index].std(axis=0)
        for member in range(5):
            scored[f"member_{member}_p_{label}"] = morphology[member, :, class_index]
    scored["p_preserve"] = mean_preserve[:, 1]
    scored["std_p_preserve"] = preserve[:, :, 1].std(axis=0)
    for class_index, label in enumerate(HARMONIC_DISPLAY_LABELS):
        scored[f"p_harmonic_{label}"] = mean_harmonic[:, class_index]
        scored[f"std_p_harmonic_{label}"] = harmonic[:, :, class_index].std(axis=0)
    sorted_probability = np.sort(mean_morphology, axis=1)
    scored["morphology_entropy"] = -np.sum(
        mean_morphology * np.log(np.clip(mean_morphology, 1.0e-12, 1.0)), axis=1
    )
    scored["morphology_margin"] = sorted_probability[:, -1] - sorted_probability[:, -2]
    scored["ensemble_disagreement"] = morphology.std(axis=0).mean(axis=1)
    scored["predicted_morphology"] = np.asarray(MORPHOLOGY_CLASSES, dtype=object)[
        mean_morphology.argmax(axis=1)
    ]
    factor_index = mean_harmonic.argmax(axis=1)
    scored["predicted_period_factor"] = np.asarray(HARMONIC_FACTORS)[factor_index]
    scored["predicted_period_d"] = (
        pd.to_numeric(scored["period_d"], errors="coerce")
        * scored["predicted_period_factor"]
    )
    scored["model_version"] = MODEL_VERSION
    scored["model_profile"] = profile
    scored["input_contract_version"] = A2V1_TEACHER_INPUT_CONTRACT
    scored["teacher_training_sector"] = TEACHER_TRAINING_SECTOR
    scored["score_policy"] = TEACHER_SCORE_POLICY
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_version": MODEL_VERSION,
        "profile": profile,
        "teacher_training_sector": TEACHER_TRAINING_SECTOR,
        "score_policy": TEACHER_SCORE_POLICY,
        "input_contract_version": A2V1_TEACHER_INPUT_CONTRACT,
        "native_h5_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "n_rows": int(len(scored)),
        "n_unique_tics": int(pd.to_numeric(scored["tic"], errors="coerce").nunique()),
        "sectors": sorted(
            pd.to_numeric(scored.get("sector", pd.Series(dtype=float)), errors="coerce")
            .dropna()
            .astype(int)
            .unique()
            .tolist()
        ),
        "device": str(device),
        "torch_version": str(torch.__version__),
        "torch_cuda_version": str(torch.version.cuda),
        "cuda_device_name": (
            str(torch.cuda.get_device_name(0)) if device.type == "cuda" else ""
        ),
        "batch_size": int(batch_size),
        "workers": int(workers),
        "allow_injections": bool(allow_injections),
        "checkpoints": [str(Path(path)) for path in checkpoint_paths],
        "checkpoint_sha256": {
            str(Path(path)): _sha256(Path(path)) for path in checkpoint_paths
        },
        "checkpoint_folds": [int(value["fold"]) for value in checkpoints],
        "metadata_columns_by_fold": [list(value.columns) for value in normalizations],
        "predicted_class_counts": {
            str(key): int(value)
            for key, value in scored["predicted_morphology"]
            .value_counts()
            .sort_index()
            .items()
        },
    }
    return scored, summary


def score_harmonic_teacher_to_disk(
    *,
    candidates_path: Path,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
    out_dir: Path,
    profile: str = SELECTED_TEACHER_PROFILE,
    batch_size: int = 32,
    workers: int = 4,
    require_cuda: bool = True,
) -> dict[str, Any]:
    candidates_path = Path(candidates_path)
    candidates = (
        pd.read_parquet(candidates_path)
        if candidates_path.suffix.lower() == ".parquet"
        else pd.read_csv(candidates_path, low_memory=False)
    )
    scored, summary = score_harmonic_teacher_ensemble(
        candidates=candidates,
        native_h5=native_h5,
        checkpoint_paths=checkpoint_paths,
        profile=profile,
        batch_size=batch_size,
        workers=workers,
        require_cuda=require_cuda,
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    score_path = out_dir / "teacher_v1_real_candidate_scores.parquet"
    try:
        scored.to_parquet(score_path, compression="zstd", index=False)
    except (ImportError, ModuleNotFoundError, ValueError):
        score_path = score_path.with_suffix(".csv")
        scored.to_csv(score_path, index=False)
    ranked = rank_planet_enrichment(scored)
    ranking_path = out_dir / "teacher_v1_planet_enrichment_ranked.parquet"
    try:
        ranked.to_parquet(ranking_path, compression="zstd", index=False)
    except (ImportError, ModuleNotFoundError, ValueError):
        ranking_path = ranking_path.with_suffix(".csv")
        ranked.to_csv(ranking_path, index=False)
    summary["outputs"] = {
        "scores": str(score_path),
        "planet_enrichment_ranking": str(ranking_path),
        "summary": str(out_dir / "summary.json"),
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


__all__ = [
    "HARMONIC_DISPLAY_LABELS",
    "SELECTED_TEACHER_PROFILE",
    "prepare_inference_rows",
    "score_harmonic_teacher_ensemble",
    "score_harmonic_teacher_to_disk",
]
