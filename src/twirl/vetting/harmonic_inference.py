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
    verify_native_candidate_binding,
    verify_raw_pair_contract,
)
from twirl.vetting.harmonic_training import (
    ENCODER_PRETRAINING_CACHE_SCHEMA,
    TEACHER_NATIVE_V2_CHECKPOINT_NAMESPACE,
    validate_native_v2_provenance_shape,
)
from twirl.vetting.recovery50_teacher import leakage_columns


SELECTED_TEACHER_PROFILE = "shape_plus_periodogram_bls"
TEACHER_TRAINING_SECTOR = 56
TEACHER_SCORE_POLICY = "active_learning_priority_not_scientific_classification"
TEACHER_SCORE_ARTIFACT_CONTRACT = "s56_harmonic_teacher_scores_v2"
RECOVERY_SCORE_ARTIFACT_CONTRACT = "s56_a2v1_teacher_recovery_scores_v2"
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
    ensemble_provenance: tuple[str, str, str, str] | None = None
    for path in checkpoint_paths:
        checkpoint = torch.load(Path(path), map_location="cpu", weights_only=False)
        if checkpoint.get("model_version") != MODEL_VERSION:
            raise ValueError(
                f"{path} has model_version={checkpoint.get('model_version')!r}"
            )
        if checkpoint.get("input_contract_version") != RAW_PAIR_CONTRACT_VERSION:
            raise ValueError(f"{path} has the wrong native input contract")
        try:
            validate_native_v2_provenance_shape(
                checkpoint,
                artifact=f"teacher checkpoint {path}",
            )
        except RuntimeError as exc:
            raise ValueError(str(exc)) from exc
        if checkpoint.get("encoder_pretraining_cache_schema") != (
            ENCODER_PRETRAINING_CACHE_SCHEMA
        ):
            raise ValueError(f"{path} has stale encoder-pretraining provenance")
        pretraining_sha256 = str(checkpoint.get("encoder_pretraining_sha256", ""))
        if len(pretraining_sha256) != 64 or any(
            value not in "0123456789abcdef" for value in pretraining_sha256
        ):
            raise ValueError(f"{path} has an invalid encoder-pretraining SHA256")
        checkpoint_provenance = (
            str(checkpoint["checkpoint_namespace"]),
            str(checkpoint["input_contract_version"]),
            str(checkpoint["native_h5_sha256"]),
            str(checkpoint["training_table_sha256"]),
        )
        if ensemble_provenance is None:
            ensemble_provenance = checkpoint_provenance
        elif checkpoint_provenance != ensemble_provenance:
            raise ValueError(f"{path} has inconsistent ensemble training provenance")
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
    verify_native_contract: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Score candidates once through the data path and five times through the CNN."""

    import h5py
    import torch
    from torch.utils.data import DataLoader

    rows = prepare_inference_rows(candidates, allow_injections=allow_injections)
    native_verification: dict[str, Any] | None = None
    if verify_native_contract:
        native_verification = verify_raw_pair_contract(
            Path(native_h5),
            require_errors=True,
            require_periodograms=True,
        )
        if not native_verification["passed"]:
            raise ValueError(
                "native HDF5 failed the v2 contract: "
                + "; ".join(native_verification["failures"][:10])
            )
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
    scored["checkpoint_namespace"] = TEACHER_NATIVE_V2_CHECKPOINT_NAMESPACE
    scored["input_contract_version"] = A2V1_TEACHER_INPUT_CONTRACT
    scored["teacher_training_sector"] = TEACHER_TRAINING_SECTOR
    scored["score_policy"] = TEACHER_SCORE_POLICY
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_version": MODEL_VERSION,
        "checkpoint_namespace": TEACHER_NATIVE_V2_CHECKPOINT_NAMESPACE,
        "checkpoint_native_training_h5_sha256": checkpoints[0][
            "native_h5_sha256"
        ],
        "checkpoint_training_table_sha256": checkpoints[0][
            "training_table_sha256"
        ],
        "profile": profile,
        "teacher_training_sector": TEACHER_TRAINING_SECTOR,
        "score_policy": TEACHER_SCORE_POLICY,
        "input_contract_version": A2V1_TEACHER_INPUT_CONTRACT,
        "native_h5_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "native_contract_verification": native_verification,
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


def _input_hashes(
    *,
    candidates_path: Path,
    candidate_summary_path: Path | None,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
) -> dict[str, Any]:
    return {
        "candidate_table_sha256": _sha256(candidates_path),
        "candidate_summary_sha256": (
            _sha256(candidate_summary_path)
            if candidate_summary_path is not None
            else ""
        ),
        "native_h5_sha256": _sha256(native_h5),
        "checkpoint_sha256": {
            str(path): _sha256(path) for path in checkpoint_paths
        },
    }


def _input_artifact_summary(
    *,
    candidates_path: Path,
    candidate_summary_path: Path | None,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
    before: dict[str, Any],
    after: dict[str, Any],
) -> dict[str, Any]:
    def entry(path: Path, name: str) -> dict[str, str]:
        return {
            "path": str(path),
            "sha256_before": str(before[name]),
            "sha256_after": str(after[name]),
        }

    artifacts: dict[str, Any] = {
        "candidate_table": entry(candidates_path, "candidate_table_sha256"),
        "native_h5": entry(native_h5, "native_h5_sha256"),
        "checkpoints": [
            {
                "path": str(path),
                "sha256_before": before["checkpoint_sha256"][str(path)],
                "sha256_after": after["checkpoint_sha256"][str(path)],
            }
            for path in checkpoint_paths
        ],
    }
    if candidate_summary_path is not None:
        artifacts["candidate_summary"] = entry(
            candidate_summary_path, "candidate_summary_sha256"
        )
    return artifacts


def _stage_table(
    frame: pd.DataFrame,
    *,
    out_dir: Path,
    stem: str,
) -> tuple[Path, Path]:
    parquet_path = out_dir / f"{stem}.parquet"
    parquet_tmp = out_dir / f".{stem}.tmp.parquet"
    parquet_tmp.unlink(missing_ok=True)
    try:
        frame.to_parquet(parquet_tmp, compression="zstd", index=False)
        return parquet_tmp, parquet_path
    except (ImportError, ModuleNotFoundError, ValueError):
        parquet_tmp.unlink(missing_ok=True)
    csv_path = out_dir / f"{stem}.csv"
    csv_tmp = out_dir / f".{stem}.tmp.csv"
    csv_tmp.unlink(missing_ok=True)
    frame.to_csv(csv_tmp, index=False)
    return csv_tmp, csv_path


def score_harmonic_teacher_to_disk(
    *,
    candidates_path: Path,
    candidate_summary_path: Path | None = None,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
    out_dir: Path,
    profile: str = SELECTED_TEACHER_PROFILE,
    batch_size: int = 32,
    workers: int = 4,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Score and atomically publish a hash-bound teacher artifact set."""

    candidates_path = Path(candidates_path).resolve()
    candidate_summary_path = (
        Path(candidate_summary_path).resolve()
        if candidate_summary_path is not None
        else None
    )
    native_h5 = Path(native_h5).resolve()
    checkpoint_paths = tuple(Path(path).resolve() for path in checkpoint_paths)
    if len(checkpoint_paths) != 5:
        raise ValueError(f"expected five fold checkpoints, got {len(checkpoint_paths)}")
    native_verification = (
        verify_native_candidate_binding(
            native_h5,
            candidate_table=candidates_path,
            candidate_summary=candidate_summary_path,
            require_periodograms=True,
            expected_periodogram_n=4096,
        )
        if candidate_summary_path is not None
        else verify_raw_pair_contract(
            native_h5,
            require_errors=True,
            require_periodograms=True,
        )
    )
    if not native_verification["passed"]:
        raise ValueError(
            "native candidate binding failed: "
            + "; ".join(native_verification["failures"][:10])
        )
    before = _input_hashes(
        candidates_path=candidates_path,
        candidate_summary_path=candidate_summary_path,
        native_h5=native_h5,
        checkpoint_paths=checkpoint_paths,
    )
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
        verify_native_contract=False,
    )
    ranked = rank_planet_enrichment(scored)
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    score_tmp, score_path = _stage_table(
        scored,
        out_dir=out_dir,
        stem="teacher_v1_real_candidate_scores",
    )
    ranking_tmp, ranking_path = _stage_table(
        ranked,
        out_dir=out_dir,
        stem="teacher_v1_planet_enrichment_ranked",
    )
    after = _input_hashes(
        candidates_path=candidates_path,
        candidate_summary_path=candidate_summary_path,
        native_h5=native_h5,
        checkpoint_paths=checkpoint_paths,
    )
    if after != before:
        score_tmp.unlink(missing_ok=True)
        ranking_tmp.unlink(missing_ok=True)
        raise RuntimeError("teacher scoring inputs changed before outputs were published")
    expected_checkpoint_hashes = before["checkpoint_sha256"]
    if summary.get("checkpoint_sha256") != expected_checkpoint_hashes:
        score_tmp.unlink(missing_ok=True)
        ranking_tmp.unlink(missing_ok=True)
        raise RuntimeError("ensemble checkpoint hashes disagree with pre-scoring hashes")

    score_tmp.replace(score_path)
    ranking_tmp.replace(ranking_path)
    for selected in (score_path, ranking_path):
        alternate = selected.with_suffix(
            ".csv" if selected.suffix.lower() == ".parquet" else ".parquet"
        )
        alternate.unlink(missing_ok=True)
    summary_path = out_dir / "summary.json"
    summary.update(
        {
            "artifact_contract_version": TEACHER_SCORE_ARTIFACT_CONTRACT,
            "strict_provenance_passed": True,
            "native_contract_verification": native_verification,
            "input_artifacts": _input_artifact_summary(
                candidates_path=candidates_path,
                candidate_summary_path=candidate_summary_path,
                native_h5=native_h5,
                checkpoint_paths=checkpoint_paths,
                before=before,
                after=after,
            ),
            "outputs": {
                "scores": str(score_path),
                "planet_enrichment_ranking": str(ranking_path),
                "summary": str(summary_path),
            },
            "output_sha256": {
                "scores": _sha256(score_path),
                "planet_enrichment_ranking": _sha256(ranking_path),
            },
            "n_ranked_unique_tics": int(len(ranked)),
        }
    )
    summary_tmp = out_dir / ".summary.json.tmp"
    summary_tmp.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    summary_tmp.replace(summary_path)
    return summary


def verify_teacher_score_cache(
    *,
    summary_path: Path,
    candidates_path: Path,
    candidate_summary_path: Path,
    native_h5: Path,
    checkpoint_paths: Sequence[Path],
) -> dict[str, Any]:
    """Return whether score/ranking outputs exactly match all current inputs."""

    failures: list[str] = []
    summary_path = Path(summary_path).resolve()
    candidates_path = Path(candidates_path).resolve()
    candidate_summary_path = Path(candidate_summary_path).resolve()
    native_h5 = Path(native_h5).resolve()
    checkpoint_paths = tuple(Path(path).resolve() for path in checkpoint_paths)
    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return {"passed": False, "failures": [f"score summary: {exc}"]}
    if summary.get("artifact_contract_version") != TEACHER_SCORE_ARTIFACT_CONTRACT:
        failures.append("wrong artifact_contract_version")
    if summary.get("strict_provenance_passed") is not True:
        failures.append("strict_provenance_passed is not true")
    try:
        current = _input_hashes(
            candidates_path=candidates_path,
            candidate_summary_path=candidate_summary_path,
            native_h5=native_h5,
            checkpoint_paths=checkpoint_paths,
        )
    except OSError as exc:
        failures.append(f"input hash: {exc}")
        current = {}
    artifacts = summary.get("input_artifacts", {})
    expected_inputs: tuple[tuple[str, Path, str], ...] = (
        ("candidate_table", candidates_path, "candidate_table_sha256"),
        ("candidate_summary", candidate_summary_path, "candidate_summary_sha256"),
        ("native_h5", native_h5, "native_h5_sha256"),
    )
    for artifact_name, path, hash_name in expected_inputs:
        record = artifacts.get(artifact_name, {})
        if record.get("path") != str(path):
            failures.append(f"{artifact_name} path mismatch")
        observed_hash = current.get(hash_name, "")
        if (
            record.get("sha256_before") != observed_hash
            or record.get("sha256_after") != observed_hash
        ):
            failures.append(f"{artifact_name} SHA256 mismatch")
    checkpoint_records = artifacts.get("checkpoints", [])
    if len(checkpoint_records) != len(checkpoint_paths):
        failures.append("checkpoint count mismatch")
    else:
        current_checkpoint_hashes = current.get("checkpoint_sha256", {})
        for path, record in zip(checkpoint_paths, checkpoint_records):
            digest = current_checkpoint_hashes.get(str(path), "")
            if record.get("path") != str(path):
                failures.append(f"checkpoint path mismatch: {path}")
            if (
                record.get("sha256_before") != digest
                or record.get("sha256_after") != digest
            ):
                failures.append(f"checkpoint SHA256 mismatch: {path}")
    outputs = summary.get("outputs", {})
    output_hashes = summary.get("output_sha256", {})
    for name in ("scores", "planet_enrichment_ranking"):
        output_path = Path(str(outputs.get(name, "")))
        if not output_path.is_absolute():
            failures.append(f"{name} output path is not absolute")
            continue
        try:
            digest = _sha256(output_path)
        except OSError as exc:
            failures.append(f"{name} output: {exc}")
            continue
        if output_hashes.get(name) != digest:
            failures.append(f"{name} output SHA256 mismatch")
    binding = verify_native_candidate_binding(
        native_h5,
        candidate_table=candidates_path,
        candidate_summary=candidate_summary_path,
        require_periodograms=True,
        expected_periodogram_n=4096,
    )
    if not binding["passed"]:
        failures.extend(
            f"native binding: {failure}" for failure in binding["failures"]
        )
    return {
        "passed": not failures,
        "failures": failures,
        "native_binding": binding,
        "current_input_hashes": current,
    }


__all__ = [
    "HARMONIC_DISPLAY_LABELS",
    "RECOVERY_SCORE_ARTIFACT_CONTRACT",
    "SELECTED_TEACHER_PROFILE",
    "TEACHER_SCORE_ARTIFACT_CONTRACT",
    "prepare_inference_rows",
    "score_harmonic_teacher_ensemble",
    "score_harmonic_teacher_to_disk",
    "verify_teacher_score_cache",
]
