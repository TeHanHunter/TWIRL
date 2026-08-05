"""Development-only Teacher v4-SSL training on frozen Teacher v3 inputs.

This module is intentionally a sibling of :mod:`teacher_v3_training`.  It
validates and binds the frozen Teacher v3 release, but it never calls the
Teacher v3 full runner and it never constructs model tensors or evaluates
predictions for the already-opened fixed test.  Each SSL encoder is fitted on
real observations from four development folds.  Supervised fine-tuning uses
the same fold-specific early-stopping protocol as the Teacher v3 baseline, so
its out-of-fold metrics are matched development estimates rather than an
untouched prospective result.
"""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
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
    multitask_loss,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    attach_fold_training_weights,
    build_metadata_matrix,
    prepare_harmonic_training_rows,
)
from twirl.vetting.harmonic_training import (
    _calibration_by_source,
    _cross_entropy_from_probability,
    _evaluate_loader,
    _file_sha256,
    _loader,
    _softmax,
    _subset_metrics,
    _to_device,
    classification_metrics,
    expected_calibration_error,
    fit_temperature,
)
from twirl.vetting.harmonic_ssl import (
    HARMONIC_SSL_CONTRACT_VERSION,
    EventPreservingAugmentationConfig,
    VICRegConfig,
    augment_ssl_batch,
    build_ssl_cache_identity,
    select_ssl_fold_rows,
    validate_ssl_cache_identity,
    vicreg_loss,
)
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
    _assert_teacher_v3_inputs_unchanged,
    _read_training_table_with_stable_hash,
    _truth,
    build_teacher_v3_input_provenance,
    build_teacher_v3_native_manifest,
    prepare_teacher_v3_uncertain_masked_rows,
    validate_teacher_v3_training_table,
)


TEACHER_SSL_RUN_ID = "teacher_ssl_v1_s56_s62_a2v1_current_adp"
TEACHER_SSL_ENCODER_NAME = "teacher_ssl_v1"
TEACHER_SSL_MODEL_FACING_NAME = "Teacher v4-SSL"
TEACHER_SSL_PROFILE = TEACHER_V3_PRIMARY_PROFILE
TEACHER_SSL_CHECKPOINT_NAMESPACE = (
    "twirl_teacher_ssl_v1_s56_s62_a2v1_current_adp"
)
TEACHER_SSL_SUMMARY_SCHEMA = "twirl_teacher_ssl_training_summary_v1"
TEACHER_SSL_CHECKPOINT_SCHEMA = "twirl_teacher_ssl_checkpoint_v1"
TEACHER_SSL_OOF_SCHEMA = "twirl_teacher_ssl_development_oof_v1"
TEACHER_SSL_SECTORS: tuple[int, ...] = tuple(range(56, 63))
TEACHER_V3_EXPECTED_SUMMARY_SHA256 = (
    "0bdcb064a7e67f2304ba58b1e79c462d"
    "aa7cc8aad5ad64d2f57a86c4dae46e99"
)

_TRANSFER_PREFIXES: tuple[str, ...] = (
    "harmonic_encoder.",
    "local_encoder.",
    "periodogram_encoder.",
    "branch_gate.",
    "fusion.",
)


def _json_safe(value: Any) -> Any:
    """Convert nested metrics into strict, portable JSON values."""

    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def _json_bytes(payload: Mapping[str, Any]) -> bytes:
    safe = _json_safe(payload)
    return (
        json.dumps(safe, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _write_json(path: Path, payload: Mapping[str, Any]) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    serialized = _json_bytes(payload)
    path.write_bytes(serialized)
    return hashlib.sha256(serialized).hexdigest()


def _string_set_sha256(values: Sequence[Any]) -> str:
    payload = "\n".join(sorted({str(value) for value in values})) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def require_fresh_teacher_ssl_output_dir(path: Path) -> Path:
    """Create or accept only an empty, non-symlink pilot output directory."""

    requested = Path(path).expanduser()
    if requested.is_symlink():
        raise RuntimeError(
            f"Teacher v4-SSL output may not be a symlink: {requested}"
        )
    resolved = requested.resolve()
    if resolved.exists():
        if not resolved.is_dir():
            raise NotADirectoryError(resolved)
        prior = next(resolved.iterdir(), None)
        if prior is not None:
            raise FileExistsError(
                "Teacher v4-SSL requires a fresh empty output directory; "
                f"found prior content: {prior}"
            )
    else:
        resolved.mkdir(parents=True, exist_ok=False)
    return resolved


def _is_injected(rows: pd.DataFrame) -> pd.Series:
    if "is_injected_row" not in rows:
        return pd.Series(False, index=rows.index, dtype=bool)
    return _truth(rows["is_injected_row"])


def _assert_development_only(
    rows: pd.DataFrame,
    *,
    artifact: str,
    held_fold: int | None = None,
    real_only: bool = False,
) -> None:
    """Fail closed if any forbidden row reaches a training dataset."""

    required = {"review_id", "sector", "tic", "fixed_split", "cv_fold"}
    missing = sorted(required - set(rows.columns))
    if missing:
        raise KeyError(f"{artifact} rows lack columns: {missing}")
    if rows.empty:
        raise ValueError(f"{artifact} rows are empty")
    split = rows["fixed_split"].fillna("").astype(str).str.strip().str.lower()
    if not split.eq("development").all():
        raise RuntimeError(f"{artifact} contains fixed-test rows")
    sectors = pd.to_numeric(rows["sector"], errors="raise").astype(int)
    if not sectors.isin(TEACHER_SSL_SECTORS).all():
        raise RuntimeError(f"{artifact} contains Sector 63 or another sector")
    folds = pd.to_numeric(rows["cv_fold"], errors="raise").astype(int)
    if not folds.isin(range(5)).all():
        raise RuntimeError(f"{artifact} contains invalid development folds")
    if held_fold is not None and folds.eq(int(held_fold)).any():
        raise RuntimeError(
            f"{artifact} contains held development fold {held_fold}"
        )
    if real_only and _is_injected(rows).any():
        raise RuntimeError(f"{artifact} contains injected rows")


def _active_release_rows(source: pd.DataFrame, *, seed: int) -> pd.DataFrame:
    rows = prepare_harmonic_training_rows(source, seed=seed)
    expected = int(_truth(source["teacher_v3_training_include"]).sum())
    if len(rows) != expected:
        raise RuntimeError(
            "Teacher v4-SSL active-row contract disagrees with Teacher v3: "
            f"{len(rows)} != {expected}"
        )
    return rows


def _split_registry_sha256(
    *,
    training_table: Path,
    rows: pd.DataFrame,
) -> tuple[str, str]:
    """Hash the frozen split registry, with a canonical assignment fallback."""

    training_table = Path(training_table).resolve()
    candidates = (
        training_table.parent.parent
        / "frozen"
        / "tic_split_registry.csv",
        training_table.parent / "tic_split_registry.csv",
    )
    for candidate in candidates:
        if candidate.is_file():
            return _file_sha256(candidate), str(candidate.resolve())
    required = ["tic", "fixed_split", "cv_fold"]
    assignment = (
        rows.loc[:, required]
        .drop_duplicates()
        .sort_values(required, kind="stable")
    )
    if assignment["tic"].duplicated().any():
        raise RuntimeError(
            "TIC assignments disagree while deriving split-registry identity"
        )
    serialized = assignment.to_csv(
        index=False,
        lineterminator="\n",
    ).encode("utf-8")
    return hashlib.sha256(serialized).hexdigest(), (
        "canonical_unique_tic_assignment_from_training_table"
    )


def _code_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(__file__).resolve().parents[3],
        check=False,
        capture_output=True,
        text=True,
    )
    revision = completed.stdout.strip()
    if completed.returncode != 0 or len(revision) != 40:
        raise RuntimeError("cannot bind Teacher v4-SSL to a git revision")
    return revision


def _transfer_ssl_encoder_state(
    *,
    ssl_state: Mapping[str, Any],
    model: Any,
) -> list[str]:
    """Transfer only shared shape/fusion weights into a fresh supervised model."""

    destination = model.state_dict()
    selected: dict[str, Any] = {}
    for name, value in ssl_state.items():
        if not name.startswith(_TRANSFER_PREFIXES):
            continue
        if name not in destination:
            raise RuntimeError(f"SSL state has unknown shared parameter {name}")
        if tuple(value.shape) != tuple(destination[name].shape):
            raise RuntimeError(
                f"SSL shared parameter shape changed for {name}: "
                f"{tuple(value.shape)} != {tuple(destination[name].shape)}"
            )
        selected[name] = value
    expected = {
        name
        for name in destination
        if name.startswith(_TRANSFER_PREFIXES)
    }
    if set(selected) != expected:
        missing = sorted(expected - set(selected))
        raise RuntimeError(
            "SSL checkpoint lacks shared shape/fusion parameters: "
            f"{missing[:10]}"
        )
    destination.update(selected)
    model.load_state_dict(destination, strict=True)
    return sorted(selected)


def _atomic_torch_save(payload: Mapping[str, Any], path: Path) -> str:
    import torch

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        torch.save(dict(payload), temporary)
        digest = _file_sha256(temporary)
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    if _file_sha256(path) != digest:
        raise RuntimeError(f"checkpoint changed while being installed: {path}")
    return digest


def _projection_head(input_dim: int = 128) -> Any:
    """Return the frozen 128 -> 128 -> 64 SSL projection architecture."""

    from torch import nn

    return nn.Sequential(
        nn.Linear(int(input_dim), 128),
        nn.GELU(),
        nn.Linear(128, 64),
    )


def _loss_and_components(result: Any) -> tuple[Any, dict[str, float]]:
    """Normalize the core VICReg return shape for artifact logging."""

    if isinstance(result, tuple) and len(result) == 2:
        loss, raw_components = result
    elif isinstance(result, Mapping) and "loss" in result:
        loss = result["loss"]
        raw_components = {
            key: value for key, value in result.items() if key != "loss"
        }
    else:
        loss = result
        raw_components = {}
    components: dict[str, float] = {}
    for name, value in raw_components.items():
        if str(name) == "loss":
            continue
        try:
            components[str(name)] = float(value.detach())
        except AttributeError:
            components[str(name)] = float(value)
    return loss, components


def _evaluate_loader_with_embeddings(
    model: Any,
    loader: Any,
    *,
    device: Any,
) -> dict[str, Any]:
    import torch

    model.eval()
    collected: dict[str, list[Any]] = {
        "review_id": [],
        "tic": [],
        "embedding": [],
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
            collected["tic"].append(raw_batch["tic"].cpu().numpy())
            for name in (
                "embedding",
                "morphology_logits",
                "preserve_logits",
                "harmonic_logits",
            ):
                collected[name].append(output[name].float().cpu().numpy())
            for name in (
                "morphology_target",
                "preserve_target",
                "harmonic_target",
            ):
                collected[name].append(raw_batch[name].cpu().numpy())
    return {
        "review_id": collected["review_id"],
        **{
            name: (
                np.concatenate(values, axis=0)
                if values
                else np.empty((0, 0), dtype=float)
            )
            for name, values in collected.items()
            if name != "review_id"
        },
    }


def _write_ssl_neighbor_probe(
    *,
    rows: pd.DataFrame,
    embeddings: np.ndarray,
    held_fold: int,
    path: Path,
    neighbors: int = 15,
) -> dict[str, Any]:
    """Write held-query neighbors from the pre-finetune SSL representation."""

    values = np.asarray(embeddings, dtype=np.float32)
    if values.shape != (len(rows), 128):
        raise ValueError(
            "pre-finetune SSL embeddings must have shape "
            f"({len(rows)}, 128), observed {values.shape}"
        )
    norm = np.linalg.norm(values, axis=1, keepdims=True)
    values = values / np.maximum(norm, 1.0e-8)
    folds = pd.to_numeric(rows["cv_fold"], errors="raise").astype(int)
    labels = rows["human_label"].fillna("").astype(str)
    query_index = np.flatnonzero(folds.eq(int(held_fold)).to_numpy())
    reference_index = np.flatnonzero(
        folds.ne(int(held_fold)).to_numpy()
        & labels.ne("uncertain").to_numpy()
    )
    if not len(query_index) or not len(reference_index):
        raise RuntimeError("SSL neighbor probe has no query or reference rows")
    tics = pd.to_numeric(rows["tic"], errors="raise").astype(np.int64)
    query_tics = set(tics.iloc[query_index].astype(int))
    reference_tics = set(tics.iloc[reference_index].astype(int))
    if query_tics & reference_tics:
        raise RuntimeError("SSL neighbor probe leaks held-query TICs into references")
    k = min(int(neighbors), len(reference_index))
    records: list[dict[str, Any]] = []
    for start in range(0, len(query_index), 256):
        query_chunk = query_index[start : start + 256]
        similarity = values[query_chunk] @ values[reference_index].T
        partition = np.argpartition(
            -similarity,
            kth=k - 1,
            axis=1,
        )[:, :k]
        for row_number, source_index in enumerate(query_chunk):
            local = partition[row_number]
            order = local[np.argsort(-similarity[row_number, local])]
            selected = reference_index[order]
            neighbor_labels = labels.iloc[selected].tolist()
            query_label = str(labels.iloc[source_index])
            records.append(
                {
                    "review_id": str(rows.iloc[source_index]["review_id"]),
                    "tic": int(rows.iloc[source_index]["tic"]),
                    "sector": int(rows.iloc[source_index]["sector"]),
                    "cv_fold": int(held_fold),
                    "human_label": query_label,
                    "neighbor_review_ids": json.dumps(
                        rows.iloc[selected]["review_id"].astype(str).tolist()
                    ),
                    "neighbor_human_labels": json.dumps(neighbor_labels),
                    "neighbor_cosine_similarity": json.dumps(
                        [
                            float(value)
                            for value in similarity[row_number, order]
                        ]
                    ),
                    "same_label_fraction": (
                        float(
                            np.mean(
                                np.asarray(neighbor_labels, dtype=str)
                                == query_label
                            )
                        )
                        if query_label != "uncertain"
                        else np.nan
                    ),
                }
            )
    table = pd.DataFrame(records)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False, lineterminator="\n")
    return {
        "path": str(path),
        "sha256": _file_sha256(path),
        "n_query_rows": int(len(table)),
        "n_reference_rows": int(len(reference_index)),
        "neighbors": int(k),
        "mean_same_label_fraction_non_uncertain": float(
            table["same_label_fraction"].mean()
        ),
    }


def _run_ssl_fold(
    *,
    rows: pd.DataFrame,
    fold: int,
    out_dir: Path,
    input_provenance: Mapping[str, str],
    split_registry_sha256: str,
    code_revision: str,
    ssl_epochs: int,
    batch_size: int,
    workers: int,
    seed: int,
    require_cuda: bool,
    augmentation_config: EventPreservingAugmentationConfig,
    vicreg_config: VICRegConfig,
) -> dict[str, Any]:
    """Fit one fold-local SSL encoder without constructing forbidden rows."""

    import torch

    fold_seed = int(seed) + 1000 * int(fold)
    torch.manual_seed(fold_seed)
    np.random.seed(fold_seed)
    selected = select_ssl_fold_rows(rows, held_out_fold=int(fold))
    if isinstance(selected, tuple):
        ssl_rows, selection_audit = selected
    else:
        ssl_rows = selected
        selection_audit = {}
    ssl_rows = pd.DataFrame(ssl_rows).reset_index(drop=True)
    _assert_development_only(
        ssl_rows,
        artifact=f"SSL fold-{fold}",
        held_fold=fold,
        real_only=True,
    )
    metadata = np.empty((len(ssl_rows), 0), dtype=np.float32)
    dataset = HarmonicNativeDataset(
        ssl_rows,
        native_h5=None,
        metadata=metadata,
        profile=TEACHER_SSL_PROFILE,
    )
    loader = _loader(
        dataset,
        np.arange(len(dataset), dtype=int),
        batch_size=int(batch_size),
        shuffle=True,
        workers=int(workers),
        seed=fold_seed,
    )
    device = torch.device(
        "cuda:0" if torch.cuda.is_available() else "cpu"
    )
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required for Teacher v4-SSL but unavailable")
    model_config = HarmonicModelConfig(metadata_dim=0)
    model = build_harmonic_cnn(
        model_config,
        profile=TEACHER_SSL_PROFILE,
    ).to(device)
    projector = _projection_head(2 * model_config.embedding_dim).to(device)
    optimizer = torch.optim.AdamW(
        list(model.parameters()) + list(projector.parameters()),
        lr=3.0e-4,
        weight_decay=1.0e-4,
    )
    identity = build_ssl_cache_identity(
        training_table_sha256=str(
            input_provenance["training_table_sha256"]
        ),
        native_registry_sha256=str(
            input_provenance["native_registry_sha256"]
        ),
        split_registry_sha256=str(split_registry_sha256),
        selected_rows=ssl_rows,
        profile=TEACHER_SSL_PROFILE,
        held_out_fold=int(fold),
        seed=fold_seed,
        model_config=asdict(model_config),
        augmentation_config=augmentation_config,
        vicreg_config=vicreg_config,
        code_revision=str(code_revision),
    )
    validate_ssl_cache_identity(
        identity,
        identity.to_manifest(),
    )
    history: list[dict[str, Any]] = []
    global_step = 0
    for epoch in range(1, int(ssl_epochs) + 1):
        model.train()
        projector.train()
        epoch_totals: dict[str, float] = {"loss": 0.0}
        seen = 0
        for raw_batch in loader:
            if len(raw_batch["review_id"]) < 2:
                continue
            batch = _to_device(raw_batch, device)
            view_seed = fold_seed + epoch * 1_000_003 + global_step * 17
            view_a = augment_ssl_batch(
                batch,
                duration_min=batch["duration_min"],
                config=augmentation_config,
                seed=view_seed,
                view_index=0,
            )
            view_b = augment_ssl_batch(
                batch,
                duration_min=batch["duration_min"],
                config=augmentation_config,
                seed=view_seed,
                view_index=1,
            )
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(
                device_type=device.type,
                dtype=torch.bfloat16,
                enabled=device.type == "cuda",
            ):
                first = projector(model(view_a)["embedding"])
                second = projector(model(view_b)["embedding"])
            loss, components = _loss_and_components(
                vicreg_loss(
                    first.float(),
                    second.float(),
                    config=vicreg_config,
                )
            )
            loss.backward()
            optimizer.step()
            count = len(raw_batch["review_id"])
            epoch_totals["loss"] += float(loss.detach()) * count
            for name, value in components.items():
                epoch_totals[name] = (
                    epoch_totals.get(name, 0.0) + float(value) * count
                )
            seen += count
            global_step += 1
        if seen < 2:
            raise RuntimeError(
                f"SSL fold {fold} epoch {epoch} had no valid VICReg batch"
            )
        record: dict[str, Any] = {
            "epoch": epoch,
            **{
                name: value / max(seen, 1)
                for name, value in epoch_totals.items()
            },
        }
        history.append(record)
        print(
            f"[teacher_ssl fold={fold}] epoch={epoch} "
            f"loss={record['loss']:.6f}",
            flush=True,
        )
    fold_dir = Path(out_dir) / "encoder_pretraining" / f"fold_{fold}"
    representation_rows = rows.loc[
        rows["fixed_split"].astype(str).eq("development")
        & ~_is_injected(rows)
    ].reset_index(drop=True)
    _assert_development_only(
        representation_rows,
        artifact=f"SSL fold-{fold} representation export",
        real_only=True,
    )
    all_metadata = np.empty(
        (len(representation_rows), 0),
        dtype=np.float32,
    )
    embedding_dataset = HarmonicNativeDataset(
        representation_rows,
        native_h5=None,
        metadata=all_metadata,
        profile=TEACHER_SSL_PROFILE,
    )
    embedding_loader = _loader(
        embedding_dataset,
        np.arange(len(representation_rows), dtype=int),
        batch_size=int(batch_size),
        shuffle=False,
        workers=int(workers),
        seed=fold_seed,
    )
    representation = _evaluate_loader_with_embeddings(
        model,
        embedding_loader,
        device=device,
    )
    expected_ids = representation_rows["review_id"].astype(str).tolist()
    if representation["review_id"] != expected_ids:
        raise RuntimeError(
            f"SSL fold {fold} representation row order changed"
        )
    fold_dir.mkdir(parents=True, exist_ok=True)
    embedding_path = fold_dir / "pre_finetune_embeddings.npz"
    with embedding_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=np.asarray(expected_ids, dtype=str),
            tic=representation["tic"].astype(np.int64),
            sector=pd.to_numeric(
                representation_rows["sector"],
                errors="raise",
            ).to_numpy(dtype=np.int16),
            cv_fold=pd.to_numeric(
                representation_rows["cv_fold"],
                errors="raise",
            ).to_numpy(dtype=np.int8),
            role=np.where(
                representation_rows["cv_fold"].eq(int(fold)),
                "held_query",
                "ssl_reference",
            ),
            embedding=representation["embedding"].astype(np.float32),
        )
    neighbor_probe = _write_ssl_neighbor_probe(
        rows=representation_rows,
        embeddings=representation["embedding"],
        held_fold=fold,
        path=fold_dir / "held_query_neighbors.csv",
    )
    checkpoint_path = fold_dir / "encoder.pt"
    checkpoint = {
        "schema_version": TEACHER_SSL_CHECKPOINT_SCHEMA,
        "run_id": TEACHER_SSL_RUN_ID,
        "model_facing_name": TEACHER_SSL_MODEL_FACING_NAME,
        "model_version": MODEL_VERSION,
        "ssl_contract_version": HARMONIC_SSL_CONTRACT_VERSION,
        "profile": TEACHER_SSL_PROFILE,
        "fold": int(fold),
        "model_config": asdict(model_config),
        "projection_architecture": [128, 128, 64],
        "augmentation_config": asdict(augmentation_config),
        "vicreg_config": asdict(vicreg_config),
        "cache_identity": identity.to_manifest(),
        "selection_audit": dict(selection_audit),
        "n_ssl_rows": int(len(ssl_rows)),
        "n_ssl_tics": int(
            pd.to_numeric(ssl_rows["tic"], errors="raise").nunique()
        ),
        "ssl_review_ids_sha256": _string_set_sha256(
            ssl_rows["review_id"].astype(str)
        ),
        "ssl_epochs": int(ssl_epochs),
        "batch_size": int(batch_size),
        "seed": fold_seed,
        "pre_finetune_embeddings": str(embedding_path),
        "pre_finetune_embeddings_sha256": _file_sha256(embedding_path),
        "held_query_neighbor_probe": neighbor_probe,
        "encoder_state_dict": {
            name: value.detach().cpu()
            for name, value in model.state_dict().items()
        },
        "projection_state_dict": {
            name: value.detach().cpu()
            for name, value in projector.state_dict().items()
        },
        **dict(input_provenance),
    }
    validate_ssl_cache_identity(
        identity,
        checkpoint["cache_identity"],
    )
    checkpoint_sha256 = _atomic_torch_save(checkpoint, checkpoint_path)
    fold_dir.mkdir(parents=True, exist_ok=True)
    history_path = fold_dir / "history.csv"
    pd.DataFrame(history).to_csv(
        history_path,
        index=False,
        lineterminator="\n",
    )
    return {
        "fold": int(fold),
        "checkpoint": str(checkpoint_path),
        "checkpoint_sha256": checkpoint_sha256,
        "history": str(history_path),
        "history_sha256": _file_sha256(history_path),
        "cache_identity": identity.to_manifest(),
        "n_ssl_rows": int(len(ssl_rows)),
        "n_ssl_tics": int(ssl_rows["tic"].nunique()),
        "pre_finetune_embeddings": str(embedding_path),
        "pre_finetune_embeddings_sha256": _file_sha256(embedding_path),
        "held_query_neighbor_probe": neighbor_probe,
        "encoder_state_dict": checkpoint["encoder_state_dict"],
    }


def _run_finetune_fold(
    *,
    rows: pd.DataFrame,
    fold: int,
    ssl_result: Mapping[str, Any],
    out_dir: Path,
    input_provenance: Mapping[str, str],
    fine_tune_epochs: int,
    batch_size: int,
    workers: int,
    seed: int,
    require_cuda: bool,
    artifact_context: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Fine-tune one fresh metadata/head model on uncertain-masked labels."""

    import torch

    context = {
        "checkpoint_schema": TEACHER_SSL_CHECKPOINT_SCHEMA,
        "run_id": TEACHER_SSL_RUN_ID,
        "model_facing_name": TEACHER_SSL_MODEL_FACING_NAME,
        "checkpoint_namespace": TEACHER_SSL_CHECKPOINT_NAMESPACE,
        "label_policy": "uncertain_masked",
        "log_name": "teacher_v4_ssl",
    }
    if artifact_context is not None:
        context.update(dict(artifact_context))
    required_context = {
        "checkpoint_schema",
        "run_id",
        "model_facing_name",
        "checkpoint_namespace",
        "label_policy",
        "log_name",
    }
    missing_context = sorted(required_context - set(context))
    if missing_context:
        raise KeyError(
            "fine-tuning artifact context lacks keys: "
            f"{missing_context}"
        )
    log_name = str(context["log_name"])

    _assert_development_only(
        rows,
        artifact="uncertain-masked fine-tuning",
    )
    fold_seed = int(seed) + 10_000 + 1000 * int(fold)
    torch.manual_seed(fold_seed)
    np.random.seed(fold_seed)
    train_mask = rows["cv_fold"].ne(int(fold)).to_numpy()
    validation_mask = rows["cv_fold"].eq(int(fold)).to_numpy()
    if not np.any(train_mask) or not np.any(validation_mask):
        raise RuntimeError(f"fine-tuning fold {fold} has an empty partition")
    train_tics = set(
        pd.to_numeric(rows.loc[train_mask, "tic"], errors="raise").astype(int)
    )
    validation_tics = set(
        pd.to_numeric(
            rows.loc[validation_mask, "tic"],
            errors="raise",
        ).astype(int)
    )
    if train_tics & validation_tics:
        raise RuntimeError(f"fine-tuning fold {fold} leaks TICs")
    weighted_rows = attach_fold_training_weights(
        rows,
        fit_mask=train_mask,
    )
    metadata, normalization = build_metadata_matrix(
        weighted_rows,
        fit_mask=train_mask,
    )
    dataset = HarmonicNativeDataset(
        weighted_rows,
        native_h5=None,
        metadata=metadata,
        profile=TEACHER_SSL_PROFILE,
    )
    train_loader = _loader(
        dataset,
        np.flatnonzero(train_mask),
        batch_size=int(batch_size),
        shuffle=True,
        workers=int(workers),
        seed=fold_seed,
    )
    validation_loader = _loader(
        dataset,
        np.flatnonzero(validation_mask),
        batch_size=int(batch_size),
        shuffle=False,
        workers=int(workers),
        seed=fold_seed,
    )
    device = torch.device(
        "cuda:0" if torch.cuda.is_available() else "cpu"
    )
    if require_cuda and device.type != "cuda":
        raise RuntimeError(
            "CUDA was required for Teacher v4-SSL fine-tuning but unavailable"
        )
    model_config = HarmonicModelConfig(metadata_dim=metadata.shape[1])
    model = build_harmonic_cnn(
        model_config,
        profile=TEACHER_SSL_PROFILE,
    ).to(device)
    transferred = _transfer_ssl_encoder_state(
        ssl_state=ssl_result["encoder_state_dict"],
        model=model,
    )
    train_config = HarmonicTrainConfig(
        epochs=int(fine_tune_epochs),
        batch_size=int(batch_size),
        seed=fold_seed,
    )
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
    for epoch in range(1, int(fine_tune_epochs) + 1):
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
                loss, _ = multitask_loss(
                    output,
                    batch,
                    config=train_config,
                )
            loss.backward()
            optimizer.step()
            count = len(raw_batch["review_id"])
            total_loss += float(loss.detach()) * count
            seen += count
        validation = _evaluate_loader(
            model,
            validation_loader,
            device=device,
        )
        probability = _softmax(validation["morphology_logits"])
        metrics = classification_metrics(
            validation["morphology_target"],
            probability,
            classes=MORPHOLOGY_CLASSES,
        )
        macro_f1 = float(metrics["macro_f1"])
        record = {
            "epoch": epoch,
            "train_loss": total_loss / max(seen, 1),
            "validation_morphology_loss": (
                _cross_entropy_from_probability(
                    validation["morphology_target"],
                    probability,
                )
            ),
            "validation_macro_f1": macro_f1,
            "validation_balanced_accuracy": float(
                metrics["balanced_accuracy"]
            ),
        }
        history.append(record)
        if np.isfinite(macro_f1) and macro_f1 > best_macro_f1 + 1.0e-6:
            best_macro_f1 = macro_f1
            best_epoch = epoch
            best_state = {
                name: value.detach().cpu().clone()
                for name, value in model.state_dict().items()
            }
            stale = 0
        else:
            stale += 1
        print(
            f"[{log_name} fold={fold}] epoch={epoch} "
            f"loss={record['train_loss']:.6f} "
            f"val_macro_f1={macro_f1:.4f} stale={stale}",
            flush=True,
        )
        if stale >= train_config.patience:
            break
    if best_state is None:
        raise RuntimeError(
            f"Teacher v4-SSL fold {fold} produced no finite validation metric"
        )
    model.load_state_dict(best_state, strict=True)
    validation = _evaluate_loader_with_embeddings(
        model,
        validation_loader,
        device=device,
    )
    raw_probability = _softmax(validation["morphology_logits"])
    validation_rows = weighted_rows.loc[validation_mask].reset_index(drop=True)
    if validation["review_id"] != validation_rows["review_id"].astype(str).tolist():
        raise RuntimeError(
            f"Teacher v4-SSL fold {fold} validation order changed"
        )
    fold_metrics = {
        "morphology_by_source": _subset_metrics(
            validation_rows,
            validation["morphology_target"],
            raw_probability,
        ),
        "morphology_calibration": expected_calibration_error(
            validation["morphology_target"],
            raw_probability,
        ),
        "preserve": classification_metrics(
            validation["preserve_target"],
            _softmax(validation["preserve_logits"]),
            classes=PRESERVE_CLASSES,
        ),
        "harmonic": classification_metrics(
            validation["harmonic_target"],
            _softmax(validation["harmonic_logits"]),
            classes=HARMONIC_CLASSES,
        ),
        "temperature_calibration_scope": (
            "pending_pooled_development_oof"
        ),
    }
    fold_dir = Path(out_dir) / TEACHER_SSL_PROFILE / f"fold_{fold}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = fold_dir / "teacher.pt"
    checkpoint = {
        "schema_version": str(context["checkpoint_schema"]),
        "run_id": str(context["run_id"]),
        "model_facing_name": str(context["model_facing_name"]),
        "model_version": MODEL_VERSION,
        "ssl_contract_version": HARMONIC_SSL_CONTRACT_VERSION,
        "checkpoint_namespace": str(context["checkpoint_namespace"]),
        "label_policy": str(context["label_policy"]),
        "evaluation_context": _json_safe(dict(context)),
        "profile": TEACHER_SSL_PROFILE,
        "fold": int(fold),
        "model_config": asdict(model_config),
        "train_config": asdict(train_config),
        "metadata_normalization": asdict(normalization),
        "ssl_encoder_checkpoint": str(ssl_result["checkpoint"]),
        "ssl_encoder_checkpoint_sha256": str(
            ssl_result["checkpoint_sha256"]
        ),
        "ssl_cache_identity": dict(ssl_result["cache_identity"]),
        "transferred_parameter_names": transferred,
        "fresh_metadata_encoder_and_heads": True,
        "temperature": 1.0,
        "temperature_calibration_scope": (
            "pending_pooled_development_oof"
        ),
        "best_epoch": int(best_epoch),
        "model_state_dict": best_state,
        **dict(input_provenance),
    }
    checkpoint_sha256 = _atomic_torch_save(checkpoint, checkpoint_path)
    history_path = fold_dir / "history.csv"
    pd.DataFrame(history).to_csv(
        history_path,
        index=False,
        lineterminator="\n",
    )
    metrics_path = fold_dir / "metrics.json"
    _write_json(metrics_path, fold_metrics)
    prediction_payload: dict[str, Any] = {
        "review_id": validation["review_id"],
        "tic": pd.to_numeric(
            validation_rows["tic"],
            errors="raise",
        ).astype(np.int64),
        "sector": pd.to_numeric(
            validation_rows["sector"],
            errors="raise",
        ).astype(np.int16),
        "is_injected_row": _is_injected(validation_rows).to_numpy(),
        "human_label": validation_rows["human_label"].astype(str),
        "cv_fold": np.full(len(validation_rows), int(fold), dtype=np.int8),
        "morphology_target": validation["morphology_target"],
        "preserve_target": validation["preserve_target"],
        "harmonic_target": validation["harmonic_target"],
        **{
            f"logit_{label}": validation["morphology_logits"][:, index]
            for index, label in enumerate(MORPHOLOGY_CLASSES)
        },
        **{
            f"preserve_logit_{label}": validation["preserve_logits"][:, index]
            for index, label in enumerate(PRESERVE_CLASSES)
        },
        **{
            f"harmonic_logit_{label}": validation["harmonic_logits"][:, index]
            for index, label in enumerate(HARMONIC_CLASSES)
        },
    }
    predictions = pd.DataFrame(prediction_payload)
    predictions_path = fold_dir / "validation_predictions.csv"
    predictions.to_csv(
        predictions_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    embedding_path = fold_dir / "validation_embeddings.npz"
    with embedding_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=np.asarray(validation["review_id"], dtype=str),
            tic=validation["tic"].astype(np.int64),
            embedding=validation["embedding"].astype(np.float32),
        )
    return {
        "fold": int(fold),
        "checkpoint": str(checkpoint_path),
        "checkpoint_sha256": checkpoint_sha256,
        "history": str(history_path),
        "history_sha256": _file_sha256(history_path),
        "metrics": str(metrics_path),
        "metrics_sha256": _file_sha256(metrics_path),
        "predictions": str(predictions_path),
        "predictions_sha256": _file_sha256(predictions_path),
        "embeddings": str(embedding_path),
        "embeddings_sha256": _file_sha256(embedding_path),
        "best_epoch": int(best_epoch),
    }


def _pool_development_oof(
    *,
    rows: pd.DataFrame,
    fold_results: list[dict[str, Any]],
    out_dir: Path,
    input_provenance: Mapping[str, str],
    artifact_context: Mapping[str, Any] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Pool five inductive folds, fit one temperature, and update checkpoints."""

    import torch

    context = {
        "oof_schema": TEACHER_SSL_OOF_SCHEMA,
        "scope": "five_fold_uncertain_masked_development_oof",
        "label_policy": "uncertain_masked",
    }
    if artifact_context is not None:
        context.update(dict(artifact_context))

    parts = [pd.read_csv(result["predictions"]) for result in fold_results]
    oof = pd.concat(parts, ignore_index=True)
    if oof["review_id"].astype(str).duplicated().any():
        raise RuntimeError("Teacher v4-SSL OOF predictions contain duplicates")
    expected_ids = set(rows["review_id"].astype(str))
    observed_ids = set(oof["review_id"].astype(str))
    if observed_ids != expected_ids:
        raise RuntimeError(
            "Teacher v4-SSL OOF predictions do not exactly cover the "
            "uncertain-masked development rows"
        )
    logit_columns = [f"logit_{label}" for label in MORPHOLOGY_CLASSES]
    logits = oof.loc[:, logit_columns].to_numpy(dtype=float)
    truth = oof["morphology_target"].to_numpy(dtype=int)
    temperature = fit_temperature(logits, truth)
    probability = _softmax(logits, temperature=temperature)
    oof["morphology_prediction"] = probability.argmax(axis=1)
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        oof[f"p_{label}"] = probability[:, index]
    preserve_probability = _softmax(
        oof.loc[
            :,
            [f"preserve_logit_{label}" for label in PRESERVE_CLASSES],
        ].to_numpy(dtype=float)
    )
    harmonic_probability = _softmax(
        oof.loc[
            :,
            [f"harmonic_logit_{label}" for label in HARMONIC_CLASSES],
        ].to_numpy(dtype=float)
    )
    row_lookup = rows.set_index("review_id", drop=False)
    aligned_rows = row_lookup.loc[
        oof["review_id"].astype(str)
    ].reset_index(drop=True)
    metrics = {
        "schema_version": str(context["oof_schema"]),
        "scope": str(context["scope"]),
        "label_policy": str(context["label_policy"]),
        "temperature": float(temperature),
        "morphology_by_source": _subset_metrics(
            aligned_rows,
            truth,
            probability,
        ),
        "morphology_calibration_by_source": _calibration_by_source(
            aligned_rows,
            truth,
            probability,
        ),
        "preserve": classification_metrics(
            oof["preserve_target"].to_numpy(dtype=int),
            preserve_probability,
            classes=PRESERVE_CLASSES,
        ),
        "harmonic": classification_metrics(
            oof["harmonic_target"].to_numpy(dtype=int),
            harmonic_probability,
            classes=HARMONIC_CLASSES,
        ),
        "n_rows": int(len(oof)),
        "n_real_rows": int((~_is_injected(aligned_rows)).sum()),
        "n_tics": int(oof["tic"].nunique()),
        "review_ids_sha256": _string_set_sha256(
            oof["review_id"].astype(str)
        ),
        **dict(input_provenance),
    }
    profile_dir = Path(out_dir) / TEACHER_SSL_PROFILE
    oof_path = profile_dir / "development_oof_predictions.csv"
    oof.to_csv(
        oof_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    metrics["oof_predictions"] = str(oof_path)
    metrics["oof_predictions_sha256"] = _file_sha256(oof_path)
    calibration_path = profile_dir / "pooled_oof_calibration.json"
    calibration_sha256 = _write_json(calibration_path, metrics)

    embedding_parts: list[np.ndarray] = []
    embedding_ids: list[np.ndarray] = []
    embedding_tics: list[np.ndarray] = []
    for result in fold_results:
        with np.load(result["embeddings"]) as payload:
            embedding_ids.append(payload["review_id"].astype(str))
            embedding_tics.append(payload["tic"].astype(np.int64))
            embedding_parts.append(payload["embedding"].astype(np.float32))
    embedding_path = profile_dir / "development_oof_embeddings.npz"
    with embedding_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=np.concatenate(embedding_ids),
            tic=np.concatenate(embedding_tics),
            embedding=np.concatenate(embedding_parts, axis=0),
        )
    metrics["oof_embeddings"] = str(embedding_path)
    metrics["oof_embeddings_sha256"] = _file_sha256(embedding_path)
    metrics["calibration_sha256"] = calibration_sha256

    for result in fold_results:
        checkpoint_path = Path(result["checkpoint"])
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        if checkpoint.get("temperature_calibration_scope") != (
            "pending_pooled_development_oof"
        ):
            raise RuntimeError(
                f"{checkpoint_path} was not pending pooled calibration"
            )
        checkpoint["temperature"] = float(temperature)
        checkpoint["temperature_calibration_scope"] = (
            "pooled_development_oof"
        )
        checkpoint["pooled_oof_calibration_sha256"] = calibration_sha256
        result["checkpoint_sha256"] = _atomic_torch_save(
            checkpoint,
            checkpoint_path,
        )
    return oof, metrics


def _common_support_metrics(
    predictions: pd.DataFrame,
    *,
    rows: pd.DataFrame,
) -> dict[str, Any]:
    probability = predictions.loc[
        :,
        [f"p_{label}" for label in MORPHOLOGY_CLASSES],
    ].to_numpy(dtype=float)
    truth = predictions["morphology_target"].to_numpy(dtype=int)
    lookup = rows.set_index("review_id", drop=False)
    aligned = lookup.loc[
        predictions["review_id"].astype(str)
    ].reset_index(drop=True)
    return _subset_metrics(aligned, truth, probability)


def _teacher_v3_baseline_comparison(
    *,
    current_oof: pd.DataFrame,
    rows: pd.DataFrame,
    baseline_root: Path,
    expected_summary_sha256: str | None,
    input_provenance: Mapping[str, str],
    out_dir: Path,
) -> dict[str, Any]:
    """Compare only uncertain-masked development OOF predictions."""

    root = Path(baseline_root).resolve()
    summary_path = root / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    observed_summary_sha256 = _file_sha256(summary_path)
    if (
        expected_summary_sha256 is not None
        and observed_summary_sha256 != str(expected_summary_sha256)
    ):
        raise RuntimeError(
            "Teacher v3 baseline summary hash mismatch: "
            f"{observed_summary_sha256} != {expected_summary_sha256}"
        )
    if expected_summary_sha256 is None:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("run_id") != TEACHER_V3_RUN_ID:
            raise RuntimeError("baseline root is not the frozen Teacher v3 run")
        if str(summary.get("training_table_sha256", "")) != str(
            input_provenance["training_table_sha256"]
        ):
            raise RuntimeError(
                "Teacher v3 baseline used a different training table"
            )
    profile_root = (
        root
        / "sensitivities"
        / "uncertain_masked"
        / TEACHER_SSL_PROFILE
    )
    parts: list[pd.DataFrame] = []
    for fold in range(5):
        path = profile_root / f"fold_{fold}" / "validation_predictions.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        part = pd.read_csv(path)
        part["cv_fold"] = fold
        parts.append(part)
    baseline = pd.concat(parts, ignore_index=True)
    if baseline["review_id"].astype(str).duplicated().any():
        raise RuntimeError("Teacher v3 uncertain-masked OOF has duplicates")
    current_ids = set(current_oof["review_id"].astype(str))
    baseline_ids = set(baseline["review_id"].astype(str))
    if current_ids != baseline_ids:
        raise RuntimeError(
            "Teacher v3 and Teacher v4-SSL OOF supports are not exactly equal"
        )
    current = current_oof.sort_values("review_id", kind="stable").reset_index(
        drop=True
    )
    baseline = baseline.sort_values("review_id", kind="stable").reset_index(
        drop=True
    )
    if not np.array_equal(
        current["morphology_target"].to_numpy(dtype=int),
        baseline["morphology_target"].to_numpy(dtype=int),
    ):
        raise RuntimeError("common-support morphology truth changed")
    required_probability = [f"p_{label}" for label in MORPHOLOGY_CLASSES]
    if not set(required_probability).issubset(baseline.columns):
        raise RuntimeError(
            "Teacher v3 OOF predictions lack calibrated probabilities"
        )
    records: list[dict[str, Any]] = []
    model_metrics: dict[str, Any] = {}
    for name, predictions in (
        ("teacher_v3_uncertain_masked", baseline),
        ("teacher_v4_ssl", current),
    ):
        metrics = _common_support_metrics(predictions, rows=rows)
        model_metrics[name] = metrics
        for source in ("all", "real", "injected"):
            values = metrics.get(source, {})
            records.append(
                {
                    "model": name,
                    "source": source,
                    "n": int(values.get("n", 0)),
                    "balanced_accuracy": values.get(
                        "balanced_accuracy",
                        np.nan,
                    ),
                    "macro_f1": values.get("macro_f1", np.nan),
                }
            )
    table_path = Path(out_dir) / "teacher_v3_common_oof_comparison.csv"
    pd.DataFrame(records).to_csv(
        table_path,
        index=False,
        lineterminator="\n",
    )
    return {
        "scope": "exact_common_uncertain_masked_development_oof",
        "n_rows": int(len(current)),
        "n_tics": int(current["tic"].nunique()),
        "review_ids_sha256": _string_set_sha256(
            current["review_id"].astype(str)
        ),
        "baseline_summary": str(summary_path),
        "baseline_summary_sha256": observed_summary_sha256,
        "baseline_summary_sha256_verified": bool(
            expected_summary_sha256 is not None
        ),
        "metrics": model_metrics,
        "comparison_table": str(table_path),
        "comparison_table_sha256": _file_sha256(table_path),
    }


def run_teacher_ssl_pilot(
    *,
    training_table: Path,
    native_registry: Path,
    native_registry_summary: Path,
    out_dir: Path,
    ssl_epochs: int = 20,
    fine_tune_epochs: int = 100,
    batch_size: int = 32,
    workers: int = 8,
    seed: int = 560064,
    require_cuda: bool = True,
    baseline_teacher_v3_root: Path | None = None,
    expected_baseline_summary_sha256: str | None = (
        TEACHER_V3_EXPECTED_SUMMARY_SHA256
    ),
) -> dict[str, Any]:
    """Run the checksum-bound, development-only Teacher v4-SSL pilot."""

    for name, value in (
        ("ssl_epochs", ssl_epochs),
        ("fine_tune_epochs", fine_tune_epochs),
        ("batch_size", batch_size),
    ):
        if int(value) < 1:
            raise ValueError(f"{name} must be positive")
    if int(workers) < 0:
        raise ValueError("workers must be nonnegative")
    training_table = Path(training_table).resolve()
    native_registry = Path(native_registry).resolve()
    native_registry_summary = Path(native_registry_summary).resolve()
    out_dir = require_fresh_teacher_ssl_output_dir(out_dir)
    source, initial_table_sha256, initial_read_audit = (
        _read_training_table_with_stable_hash(
            training_table,
            artifact="Teacher v4-SSL training table",
        )
    )
    validate_teacher_v3_training_table(source)
    rows = _active_release_rows(source, seed=int(seed))
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
    if input_provenance["training_table_sha256"] != initial_table_sha256:
        raise RuntimeError("Teacher v4-SSL training table changed on startup")
    if (
        input_provenance["native_manifest_sha256"]
        != native_manifest_sha256
    ):
        raise RuntimeError("Teacher v4-SSL native manifest hash changed")
    split_registry_sha256, split_registry_source = _split_registry_sha256(
        training_table=training_table,
        rows=rows,
    )
    code_revision = _code_revision()
    ssl_development = rows.loc[
        rows["fixed_split"].eq("development") & ~_is_injected(rows)
    ].reset_index(drop=True)
    _assert_development_only(
        ssl_development,
        artifact="global SSL pool",
        real_only=True,
    )
    uncertain_masked, uncertain_audit = (
        prepare_teacher_v3_uncertain_masked_rows(rows)
    )
    supervised_development = uncertain_masked.loc[
        uncertain_masked["fixed_split"].eq("development")
    ].reset_index(drop=True)
    _assert_development_only(
        supervised_development,
        artifact="global supervised development pool",
    )
    if supervised_development["human_label"].astype(str).eq("uncertain").any():
        raise RuntimeError("uncertain labels reached supervised fine-tuning")
    augmentation_config = EventPreservingAugmentationConfig()
    vicreg_config = VICRegConfig()
    contract = {
        "schema_version": "twirl_teacher_ssl_run_contract_v1",
        "run_id": TEACHER_SSL_RUN_ID,
        "encoder_name": TEACHER_SSL_ENCODER_NAME,
        "model_facing_name": TEACHER_SSL_MODEL_FACING_NAME,
        "model_version": MODEL_VERSION,
        "ssl_contract_version": HARMONIC_SSL_CONTRACT_VERSION,
        "profile": TEACHER_SSL_PROFILE,
        "sectors": list(TEACHER_SSL_SECTORS),
        "ssl_pool": {
            "real_only": True,
            "development_only": True,
            "held_fold_excluded_per_encoder": True,
            "labels_ignored": True,
            "n_rows": int(len(ssl_development)),
            "n_tics": int(ssl_development["tic"].nunique()),
        },
        "fine_tuning": {
            "uncertain_masked": True,
            "development_only": True,
            "five_fold_oof": True,
            "held_fold_used_for_early_stopping": True,
            "estimate_status": "matched_development_not_untouched",
            "n_rows": int(len(supervised_development)),
            "n_tics": int(supervised_development["tic"].nunique()),
        },
        "fixed_test_container_integrity_validated": True,
        "fixed_test_tensors_constructed": False,
        "fixed_test_labels_used": False,
        "fixed_test_evaluated": False,
        "sector_63_rows_present": False,
        "injections_loaded_by_ssl": False,
        "scalar_metadata_loaded_by_ssl": False,
        "augmentation_config": asdict(augmentation_config),
        "vicreg_config": asdict(vicreg_config),
        "projection_architecture": [128, 128, 64],
        "split_registry_sha256": split_registry_sha256,
        "split_registry_source": split_registry_source,
        "code_revision": code_revision,
        **input_provenance,
    }
    contract_path = out_dir / "run_contract.json"
    contract_sha256 = _write_json(contract_path, contract)

    ssl_results: list[dict[str, Any]] = []
    fine_tune_results: list[dict[str, Any]] = []
    for fold in range(5):
        ssl_result = _run_ssl_fold(
            rows=rows,
            fold=fold,
            out_dir=out_dir,
            input_provenance=input_provenance,
            split_registry_sha256=split_registry_sha256,
            code_revision=code_revision,
            ssl_epochs=int(ssl_epochs),
            batch_size=int(batch_size),
            workers=int(workers),
            seed=int(seed),
            require_cuda=bool(require_cuda),
            augmentation_config=augmentation_config,
            vicreg_config=vicreg_config,
        )
        ssl_results.append(ssl_result)
        fine_tune_results.append(
            _run_finetune_fold(
                rows=supervised_development,
                fold=fold,
                ssl_result=ssl_result,
                out_dir=out_dir,
                input_provenance=input_provenance,
                fine_tune_epochs=int(fine_tune_epochs),
                batch_size=int(batch_size),
                workers=int(workers),
                seed=int(seed),
                require_cuda=bool(require_cuda),
            )
        )
        # The state tensor is a transient handoff, never part of JSON output.
        ssl_result.pop("encoder_state_dict", None)
    current_oof, oof_metrics = _pool_development_oof(
        rows=supervised_development,
        fold_results=fine_tune_results,
        out_dir=out_dir,
        input_provenance=input_provenance,
    )
    baseline_comparison: dict[str, Any] | None = None
    if baseline_teacher_v3_root is not None:
        baseline_comparison = _teacher_v3_baseline_comparison(
            current_oof=current_oof,
            rows=supervised_development,
            baseline_root=baseline_teacher_v3_root,
            expected_summary_sha256=expected_baseline_summary_sha256,
            input_provenance=input_provenance,
            out_dir=out_dir,
        )
    _assert_teacher_v3_inputs_unchanged(
        native_manifest=native_manifest,
        training_table=training_table,
        training_table_sha256=initial_table_sha256,
    )
    summary: dict[str, Any] = {
        "schema_version": TEACHER_SSL_SUMMARY_SCHEMA,
        "run_id": TEACHER_SSL_RUN_ID,
        "encoder_name": TEACHER_SSL_ENCODER_NAME,
        "model_facing_name": TEACHER_SSL_MODEL_FACING_NAME,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_version": MODEL_VERSION,
        "ssl_contract_version": HARMONIC_SSL_CONTRACT_VERSION,
        "profile": TEACHER_SSL_PROFILE,
        "code_revision": code_revision,
        "ssl_epochs": int(ssl_epochs),
        "fine_tune_epochs": int(fine_tune_epochs),
        "batch_size": int(batch_size),
        "workers": int(workers),
        "seed": int(seed),
        "training_table_initial_read": initial_read_audit,
        "native_input_manifest": str(native_manifest_path),
        "native_input_manifest_sha256": native_manifest_sha256,
        "split_registry_sha256": split_registry_sha256,
        "split_registry_source": split_registry_source,
        "run_contract": str(contract_path),
        "run_contract_sha256": contract_sha256,
        "uncertain_masked_audit": uncertain_audit,
        "ssl_folds": ssl_results,
        "fine_tuning_folds": fine_tune_results,
        "development_oof": oof_metrics,
        "teacher_v3_common_oof_comparison": baseline_comparison,
        "fixed_test_status": (
            "containing_files_integrity_validated_"
            "tensors_not_constructed_not_evaluated"
        ),
        "fixed_test_metrics": {},
        "sector_63_status": "reserved_not_loaded",
        "automatic_production_promotion": False,
        "student_training_blocked": True,
        **input_provenance,
    }
    safe_summary = _json_safe(summary)
    if not isinstance(safe_summary, dict):
        raise TypeError("Teacher v4-SSL summary must remain a JSON mapping")
    _write_json(out_dir / "summary.json", safe_summary)
    return safe_summary


__all__ = [
    "TEACHER_SSL_ENCODER_NAME",
    "TEACHER_SSL_MODEL_FACING_NAME",
    "TEACHER_SSL_PROFILE",
    "TEACHER_SSL_RUN_ID",
    "TEACHER_V3_EXPECTED_SUMMARY_SHA256",
    "require_fresh_teacher_ssl_output_dir",
    "run_teacher_ssl_pilot",
]
