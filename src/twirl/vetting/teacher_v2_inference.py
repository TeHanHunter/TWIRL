"""Frozen five-member ensemble inference for S56 Teacher v2."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    HARMONIC_CLASSES,
    MORPHOLOGY_CLASSES,
    PRESERVE_CLASSES,
    HarmonicModelConfig,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    MetadataNormalization,
    collate_native_samples,
)
from twirl.vetting.harmonic_inputs import HARMONIC_FACTORS, RAW_PAIR_CONTRACT_VERSION, native_group_path
from twirl.vetting.teacher_v2 import COMPACT_CLASSES, TEACHER_V2_MODEL_VERSION
from twirl.vetting.teacher_v2_cnn import build_teacher_v2_cnn
from twirl.vetting.teacher_v2_training import build_teacher_v2_metadata_matrix


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _softmax(values: np.ndarray, *, temperature: float = 1.0) -> np.ndarray:
    scaled = np.asarray(values, dtype=np.float64) / max(float(temperature), 1.0e-6)
    scaled -= np.max(scaled, axis=1, keepdims=True)
    exponential = np.exp(scaled)
    return exponential / np.maximum(exponential.sum(axis=1, keepdims=True), 1.0e-12)


def prepare_teacher_v2_inference_rows(
    candidates: pd.DataFrame,
    *,
    native_h5: Path,
) -> pd.DataFrame:
    rows = candidates.copy().reset_index(drop=True)
    required = {"tic", "period_d", "t0_bjd", "duration_min"}
    missing = sorted(required - set(rows.columns))
    if missing:
        raise KeyError(f"Teacher-v2 candidate table is missing columns: {missing}")
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
        raise ValueError(f"Teacher-v2 candidates contain {int(invalid.sum())} invalid ephemerides")
    rows["tic"] = rows["tic"].astype(np.int64)
    if "review_id" not in rows:
        rows["review_id"] = [f"teacher-v2-score-{index}" for index in range(len(rows))]
    rows["review_id"] = rows["review_id"].fillna("").astype(str)
    if rows["review_id"].eq("").any() or rows["review_id"].duplicated().any():
        raise ValueError("Teacher-v2 inference review IDs must be nonempty and unique")
    if "native_group_path" not in rows:
        rows["native_group_path"] = [native_group_path(row) for row in rows.to_dict("records")]
    rows["native_h5_path"] = str(Path(native_h5))
    for name in ("morphology", "preserve", "harmonic", "compact"):
        rows[f"{name}_target_index"] = -1
        rows[f"{name}_weight"] = np.float32(0.0)
    rows["pretrain_target"] = -1
    rows["teacher_v2_role"] = "inference"
    return rows


def _normalization(checkpoint: Mapping[str, Any]) -> MetadataNormalization:
    payload = checkpoint["metadata_normalization"]
    return MetadataNormalization(
        columns=tuple(str(value) for value in payload["columns"]),
        center=tuple(float(value) for value in payload["center"]),
        scale=tuple(float(value) for value in payload["scale"]),
    )


def _load_ensemble(
    checkpoint_paths: Sequence[Path],
    *,
    profile: str,
    device: Any,
) -> tuple[list[Any], list[dict[str, Any]], list[MetadataNormalization]]:
    import torch

    if len(checkpoint_paths) != 5:
        raise ValueError(f"Teacher-v2 ensemble requires five checkpoints, got {len(checkpoint_paths)}")
    payload: list[tuple[int, Any, dict[str, Any], MetadataNormalization]] = []
    for path in checkpoint_paths:
        checkpoint = torch.load(Path(path), map_location="cpu", weights_only=False)
        if checkpoint.get("model_version") != TEACHER_V2_MODEL_VERSION:
            raise ValueError(f"{path} is not a Teacher-v2 checkpoint")
        if checkpoint.get("input_contract_version") != RAW_PAIR_CONTRACT_VERSION:
            raise ValueError(f"{path} has the wrong native input contract")
        if checkpoint.get("profile") != profile:
            raise ValueError(f"{path} has profile={checkpoint.get('profile')!r}, expected {profile!r}")
        fold = int(checkpoint.get("fold", -1))
        model = build_teacher_v2_cnn(
            HarmonicModelConfig(**checkpoint["model_config"]), profile=profile
        ).to(device)
        model.load_state_dict(checkpoint["model_state_dict"], strict=True)
        model.eval()
        payload.append((fold, model, checkpoint, _normalization(checkpoint)))
    folds = [value[0] for value in payload]
    if sorted(folds) != list(range(5)):
        raise ValueError(f"Teacher-v2 checkpoint folds must be 0..4; got {folds}")
    payload.sort(key=lambda value: value[0])
    return (
        [value[1] for value in payload],
        [value[2] for value in payload],
        [value[3] for value in payload],
    )


def _to_device(batch: Mapping[str, Any], device: Any) -> dict[str, Any]:
    import torch

    return {
        key: value.to(device, non_blocking=True) if isinstance(value, torch.Tensor) else value
        for key, value in batch.items()
    }


def score_teacher_v2_ensemble(
    *,
    candidates: pd.DataFrame,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
    profile: str,
    batch_size: int = 32,
    workers: int = 8,
    require_cuda: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Score one sector without consulting truth or human labels as inputs."""

    import h5py
    import torch
    from torch.utils.data import DataLoader

    rows = prepare_teacher_v2_inference_rows(candidates, native_h5=native_h5)
    with h5py.File(native_h5, "r") as h5:
        if str(h5.attrs.get("contract_version", "")) != RAW_PAIR_CONTRACT_VERSION:
            raise ValueError("Teacher-v2 native HDF5 has the wrong contract")
        missing_groups = [
            value
            for value in rows["native_group_path"].astype(str).drop_duplicates()
            if value not in h5
        ]
    if missing_groups:
        raise ValueError(f"native HDF5 is missing {len(missing_groups)} groups; first={missing_groups[:5]}")
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("Teacher-v2 inference requires CUDA")
    models, checkpoints, normalizations = _load_ensemble(
        checkpoint_paths, profile=profile, device=device
    )
    metadata_matrices = [
        build_teacher_v2_metadata_matrix(
            rows,
            fit_mask=np.zeros(len(rows), dtype=bool),
            normalization=normalization,
        )[0]
        for normalization in normalizations
    ]
    dataset = HarmonicNativeDataset(
        rows,
        native_h5=None,
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
    row_lookup = {value: index for index, value in enumerate(rows["review_id"].astype(str))}
    head_names = ("morphology", "preserve", "harmonic", "compact")
    members: dict[str, list[list[np.ndarray]]] = {
        name: [[] for _ in models] for name in head_names
    }
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
                temperatures = {
                    "morphology": float(checkpoint["morphology_temperature"]),
                    "compact": float(checkpoint["compact_temperature"]),
                    "preserve": 1.0,
                    "harmonic": 1.0,
                }
                for name in head_names:
                    members[name][member].append(
                        _softmax(
                            output[f"{name}_logits"].float().cpu().numpy(),
                            temperature=temperatures[name],
                        )
                    )
            if batch_index % 100 == 0:
                print(
                    f"[teacher-v2-inference] batches={batch_index:,} "
                    f"rows={len(output_order):,}/{len(rows):,}",
                    flush=True,
                )
    if output_order != rows["review_id"].astype(str).tolist():
        raise RuntimeError("Teacher-v2 inference output order changed")
    arrays = {
        name: np.stack(
            [np.concatenate(values, axis=0) for values in member_values], axis=0
        )
        for name, member_values in members.items()
    }
    means = {name: values.mean(axis=0) for name, values in arrays.items()}
    scored = candidates.copy().reset_index(drop=True)
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        scored[f"p_{label}"] = means["morphology"][:, index]
        scored[f"std_p_{label}"] = arrays["morphology"][:, :, index].std(axis=0)
    scored["p_preserve"] = means["preserve"][:, PRESERVE_CLASSES.index("preserve")]
    scored["std_p_preserve"] = arrays["preserve"][:, :, PRESERVE_CLASSES.index("preserve")].std(axis=0)
    scored["p_compact_transit"] = means["compact"][:, COMPACT_CLASSES.index("compact_transit")]
    scored["std_p_compact_transit"] = arrays["compact"][:, :, COMPACT_CLASSES.index("compact_transit")].std(axis=0)
    for member in range(5):
        scored[f"member_{member}_p_compact_transit"] = arrays["compact"][member, :, 1]
    for index, label in enumerate(HARMONIC_CLASSES):
        scored[f"p_harmonic_{label}"] = means["harmonic"][:, index]
    scored["predicted_morphology"] = np.asarray(MORPHOLOGY_CLASSES, dtype=object)[
        means["morphology"].argmax(axis=1)
    ]
    factor_index = means["harmonic"].argmax(axis=1)
    scored["predicted_period_factor"] = np.asarray(HARMONIC_FACTORS)[factor_index]
    scored["predicted_period_d"] = pd.to_numeric(scored["period_d"], errors="coerce") * scored[
        "predicted_period_factor"
    ]
    compact = means["compact"]
    scored["compact_entropy"] = -np.sum(
        compact * np.log(np.clip(compact, 1.0e-12, 1.0)), axis=1
    )
    scored["model_version"] = TEACHER_V2_MODEL_VERSION
    scored["model_profile"] = profile
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_version": TEACHER_V2_MODEL_VERSION,
        "profile": profile,
        "n_rows": int(len(scored)),
        "n_unique_tics": int(pd.to_numeric(scored["tic"], errors="coerce").nunique()),
        "device": str(device),
        "torch_version": str(torch.__version__),
        "cuda_version": str(torch.version.cuda),
        "checkpoints": [str(Path(path)) for path in checkpoint_paths],
        "checkpoint_sha256": {str(Path(path)): _sha256(Path(path)) for path in checkpoint_paths},
    }
    return scored, summary


__all__ = [
    "prepare_teacher_v2_inference_rows",
    "score_teacher_v2_ensemble",
]
