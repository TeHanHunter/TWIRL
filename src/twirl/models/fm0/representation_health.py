"""Development-only representation-health checks for TWIRL-FM0.1.

These checks intentionally use only the immutable FM input release and the
``poc_development`` source partition.  They measure whether an encoder has a
non-collapsed embedding and whether repeated observations of a held-out Gaia
component are closer than unrelated observations.  They do *not* use labels,
BLS products, candidate tables, injection truth, or the sealed source test.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
import csv
import hashlib
import math
from pathlib import Path
from typing import Any

import numpy as np

from .dataset import (
    WINDOW_LENGTH,
    collate_fm0_samples,
    move_batch_to_device,
    prepare_model_window,
)
from .input_release import (
    MANIFEST_COLUMNS,
    deterministic_training_window,
    extract_window,
    load_input_release,
)
from .registry import sha256_file
from .validation import (
    read_json,
    validate_real_run_release,
)


REPRESENTATION_HEALTH_SCHEMA_VERSION = "twirl_fm0_1_representation_health_v1"


def _stable_seed(*parts: object) -> int:
    digest = hashlib.sha256(
        "\x1f".join(str(part) for part in parts).encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def _finite_float(value: float | np.floating[Any]) -> float | None:
    candidate = float(value)
    return candidate if math.isfinite(candidate) else None


def cosine_similarity_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Return one cosine similarity per pair of row embeddings."""

    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    if left_array.ndim != 2 or right_array.shape != left_array.shape:
        raise ValueError("paired embeddings must have the same two-dimensional shape")
    denominator = np.linalg.norm(left_array, axis=1) * np.linalg.norm(
        right_array, axis=1
    )
    result = np.full(left_array.shape[0], np.nan, dtype=np.float64)
    valid = denominator > 0.0
    result[valid] = (
        np.sum(left_array[valid] * right_array[valid], axis=1) / denominator[valid]
    )
    return result


def source_clustered_mean_interval(
    values: np.ndarray,
    cluster_ids: Sequence[str],
    *,
    seed: int,
    n_bootstrap: int = 1_000,
) -> dict[str, Any]:
    """Return a deterministic percentile interval after source clustering.

    Every physical source can produce many correlated windows or sector visits.
    We therefore average within a leakage component before resampling components
    with replacement.  The result is a descriptive 95% interval for the frozen
    development selection, not a sealed-test significance claim.
    """

    array = np.asarray(values, dtype=np.float64)
    identifiers = tuple(str(value) for value in cluster_ids)
    if array.ndim != 1 or array.shape[0] != len(identifiers):
        raise ValueError("clustered values and identifiers must align")
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be positive")
    grouped: dict[str, list[float]] = defaultdict(list)
    for value, identifier in zip(array, identifiers):
        if math.isfinite(float(value)):
            grouped[identifier].append(float(value))
    means = np.asarray(
        [np.mean(grouped[key]) for key in sorted(grouped) if grouped[key]],
        dtype=np.float64,
    )
    if means.size == 0:
        raise ValueError("clustered interval has no finite values")
    rng = np.random.default_rng(int(seed))
    samples = np.mean(
        means[rng.integers(0, means.size, size=(n_bootstrap, means.size))], axis=1
    )
    return {
        "n_source_clusters": int(means.size),
        "mean": _finite_float(np.mean(means)),
        "confidence_level": 0.95,
        "bootstrap_replicates": int(n_bootstrap),
        "lower": _finite_float(np.quantile(samples, 0.025)),
        "upper": _finite_float(np.quantile(samples, 0.975)),
    }


def summarize_embedding_matrix(embeddings: np.ndarray) -> dict[str, Any]:
    """Summarize rank, variance, and near-duplicate dimensions.

    The effective rank is computed from the centered covariance eigenvalue
    distribution.  It is descriptive evidence for representation collapse,
    not a transfer or scientific-performance metric by itself.
    """

    array = np.asarray(embeddings, dtype=np.float64)
    if array.ndim != 2 or array.shape[0] < 2 or array.shape[1] < 1:
        raise ValueError("embeddings must contain at least two rows and one dimension")
    if not np.all(np.isfinite(array)):
        raise ValueError("embeddings contain non-finite values")
    variances = np.var(array, axis=0, ddof=1)
    zero_dimensions = int(np.count_nonzero(variances <= 1.0e-12))
    centered = array - np.mean(array, axis=0, keepdims=True)
    covariance = centered.T @ centered / float(array.shape[0] - 1)
    eigenvalues = np.clip(np.linalg.eigvalsh(covariance), 0.0, None)
    total = float(np.sum(eigenvalues))
    if total <= 0.0:
        effective_rank = 0.0
        leading_share = 1.0
    else:
        weights = eigenvalues / total
        positive = weights[weights > 0.0]
        effective_rank = float(np.exp(-np.sum(positive * np.log(positive))))
        leading_share = float(np.max(weights))

    active = variances > 1.0e-12
    duplicate_pairs = 0
    if int(np.count_nonzero(active)) >= 2:
        correlation = np.corrcoef(array[:, active], rowvar=False)
        upper = np.triu(np.ones(correlation.shape, dtype=bool), k=1)
        duplicate_pairs = int(
            np.count_nonzero(np.abs(correlation[upper]) >= 1.0 - 1.0e-7)
        )
    return {
        "n_embeddings": int(array.shape[0]),
        "embedding_dim": int(array.shape[1]),
        "per_dimension_variance_min": _finite_float(np.min(variances)),
        "per_dimension_variance_median": _finite_float(np.median(variances)),
        "per_dimension_variance_max": _finite_float(np.max(variances)),
        "zero_or_constant_dimensions": zero_dimensions,
        "near_duplicate_dimension_pairs": duplicate_pairs,
        "effective_rank": _finite_float(effective_rank),
        "leading_principal_component_share": _finite_float(leading_share),
    }


def paired_similarity_summary(
    paired_left: np.ndarray,
    paired_right: np.ndarray,
    *,
    cluster_ids: Sequence[str],
) -> dict[str, Any]:
    """Compare same-window masked pairs with deterministically unrelated rows."""

    left = np.asarray(paired_left, dtype=np.float64)
    right = np.asarray(paired_right, dtype=np.float64)
    if left.ndim != 2 or right.shape != left.shape or left.shape[0] < 2:
        raise ValueError("paired similarity requires at least two equal-shaped rows")
    paired = cosine_similarity_rows(left, right)
    unrelated = cosine_similarity_rows(left, np.roll(right, shift=1, axis=0))
    valid = np.isfinite(paired) & np.isfinite(unrelated)
    if not np.any(valid):
        raise ValueError("paired similarity has no finite comparisons")
    difference = paired[valid] - unrelated[valid]
    identifiers = tuple(str(value) for value in cluster_ids)
    if len(identifiers) != left.shape[0]:
        raise ValueError("paired similarity cluster IDs must align with rows")
    return {
        "n_pairs": int(np.count_nonzero(valid)),
        "paired_cosine_mean": _finite_float(np.mean(paired[valid])),
        "unrelated_cosine_mean": _finite_float(np.mean(unrelated[valid])),
        "paired_minus_unrelated_mean": _finite_float(np.mean(difference)),
        "paired_minus_unrelated_median": _finite_float(np.median(difference)),
        "paired_minus_unrelated_source_clustered_95_interval": (
            source_clustered_mean_interval(
                difference,
                tuple(
                    identifier
                    for identifier, include in zip(identifiers, valid)
                    if include
                ),
                seed=_stable_seed("fm0-paired-separation-bootstrap"),
            )
        ),
    }


def same_component_retrieval_summary(
    embeddings: np.ndarray,
    component_ids: Sequence[str],
) -> dict[str, Any]:
    """Measure top-one retrieval of another visit from the same component.

    Queries without a second visit are excluded.  The identity is evaluator
    metadata only; it is never supplied to the encoder.
    """

    array = np.asarray(embeddings, dtype=np.float64)
    groups = tuple(str(value) for value in component_ids)
    if array.ndim != 2 or array.shape[0] != len(groups):
        raise ValueError("retrieval embeddings and component IDs must align")
    if array.shape[0] < 2:
        return {"status": "not_available", "reason": "fewer_than_two_visits"}
    counts: dict[str, int] = defaultdict(int)
    for group in groups:
        counts[group] += 1
    eligible = np.asarray([counts[group] >= 2 for group in groups], dtype=bool)
    if not np.any(eligible):
        return {"status": "not_available", "reason": "no_repeated_components"}
    norm = np.linalg.norm(array, axis=1, keepdims=True)
    safe = np.divide(array, norm, out=np.zeros_like(array), where=norm > 0.0)
    similarity = safe @ safe.T
    np.fill_diagonal(similarity, -np.inf)
    nearest = np.argmax(similarity, axis=1)
    recovered = np.asarray(
        [groups[index] == groups[int(nearest[index])] for index in range(len(groups))],
        dtype=bool,
    )
    eligible_groups = tuple(
        group for group, include in zip(groups, eligible) if include
    )
    return {
        "status": "available",
        "n_visit_embeddings": int(array.shape[0]),
        "n_repeated_component_queries": int(np.count_nonzero(eligible)),
        "top1_same_component_retrieval": _finite_float(np.mean(recovered[eligible])),
        "top1_source_clustered_95_interval": source_clustered_mean_interval(
            recovered[eligible].astype(np.float64),
            eligible_groups,
            seed=_stable_seed("fm0-same-component-retrieval-bootstrap"),
        ),
    }


def _read_development_rows(
    release_root: Path,
    *,
    manifest_sha256: str,
    max_components: int,
) -> list[dict[str, str]]:
    manifest = release_root / "manifest.csv"
    if sha256_file(manifest) != manifest_sha256:
        raise ValueError("FM0 release manifest differs from the run binding")
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise ValueError("FM0 release manifest schema differs from the contract")
        rows = [dict(row) for row in reader]
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["source_partition"] != "poc_development":
            continue
        if row["scientific_training_eligible"] != "True":
            raise ValueError("development selection contains a training-ineligible row")
        grouped[row["leakage_component_id"]].append(row)
    if not grouped:
        raise ValueError("FM0 release has no poc_development rows")
    if max_components <= 1:
        raise ValueError("max_components must be at least two")
    ordered_components = sorted(
        grouped,
        key=lambda component: hashlib.sha256(
            f"twirl_fm0_1_representation_health_v1:{component}".encode("utf-8")
        ).digest(),
    )
    selected = ordered_components[:max_components]
    result: list[dict[str, str]] = []
    for component in selected:
        result.extend(sorted(grouped[component], key=lambda row: row["observation_key"]))
    return result


def _model_window(
    *,
    release_root: Path,
    row: Mapping[str, str],
    variant: str,
    mask_label: str,
    unmasked: bool,
) -> dict[str, np.ndarray]:
    relative = Path(row["relative_path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("release row has an unsafe shard path")
    shard = (release_root / relative).resolve(strict=True)
    if release_root not in shard.parents or sha256_file(shard) != row["sha256"]:
        raise ValueError("release shard is unavailable or differs from its manifest")
    release = load_input_release(shard)
    spec = deterministic_training_window(
        release,
        observation_key=row["observation_key"],
        epoch=0,
        draw_index=0,
    )
    window = extract_window(
        release, segment_id=spec.segment_id, start_offset=spec.start_offset
    )
    sample = prepare_model_window(
        window,
        variant=variant,
        mask_seed=_stable_seed("fm0-representation-health", mask_label, row["observation_key"]),
        window_length=WINDOW_LENGTH,
    )
    if unmasked:
        sample["temporal_mask"] = np.zeros(WINDOW_LENGTH, dtype=bool)
        sample["reconstruction_mask"] = np.zeros_like(
            sample["reconstruction_mask"], dtype=bool
        )
    return sample


def _encode(
    *,
    model: Any,
    samples: Sequence[Mapping[str, np.ndarray]],
    device: Any,
    batch_size: int,
    with_reconstruction: bool,
    huber_delta: float,
) -> tuple[np.ndarray, dict[str, Any] | None]:
    try:
        import torch
        import torch.nn.functional as functional
    except ImportError as exc:  # pragma: no cover - base package keeps Torch optional
        raise RuntimeError("PyTorch is required for FM0 representation evaluation") from exc
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    embeddings: list[np.ndarray] = []
    losses: list[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(samples), batch_size):
            batch = collate_fm0_samples(samples[start : start + batch_size])
            batch = move_batch_to_device(batch, device)
            output = model(batch)
            value = output["z_window"].detach().cpu().numpy()
            embeddings.append(value)
            if with_reconstruction:
                mask = (
                    batch["reconstruction_mask"]
                    & batch["flux_valid"]
                    & batch["view_present"][:, :, None]
                )
                loss = functional.huber_loss(
                    output["reconstruction"],
                    batch["flux"],
                    delta=float(huber_delta),
                    reduction="none",
                )
                for item, item_mask in zip(loss, mask):
                    valid = item[item_mask]
                    if valid.numel():
                        losses.append(valid.detach().cpu().numpy())
    merged = np.concatenate(embeddings, axis=0)
    if not with_reconstruction:
        return merged, None
    all_losses = np.concatenate(losses) if losses else np.empty(0, dtype=np.float32)
    return merged, {
        "masked_valid_target_count": int(all_losses.size),
        "masked_huber_mean": (
            _finite_float(np.mean(all_losses)) if all_losses.size else None
        ),
    }


def _load_trusted_model(run_dir: Path, *, device: Any) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - base package keeps Torch optional
        raise RuntimeError("PyTorch is required for FM0 representation evaluation") from exc
    from .model import FM0ModelConfig, TWIRLFM0

    validation = validate_real_run_release(run_dir, inspect_checkpoint=True)
    contract = read_json(run_dir / "run_contract.json")
    checkpoint = torch.load(run_dir / "checkpoint.pt", map_location=device, weights_only=False)
    model = TWIRLFM0(FM0ModelConfig(**dict(checkpoint["model_config"]))).to(device)
    model.load_state_dict(checkpoint["model_state"], strict=True)
    model.eval()
    return model, contract, validation


def evaluate_representation_health(
    *,
    run_dir: str | Path,
    device: str = "cpu",
    max_components: int = 256,
    batch_size: int = 8,
    random_control_seed: int = 0,
) -> dict[str, Any]:
    """Evaluate one trained encoder only on its held-out development partition."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - base package keeps Torch optional
        raise RuntimeError("PyTorch is required for FM0 representation evaluation") from exc
    root = Path(run_dir).resolve(strict=True)
    target_device = torch.device(device)
    model, contract, run_validation = _load_trusted_model(root, device=target_device)
    release = contract.get("input_release")
    if not isinstance(release, Mapping):
        raise ValueError("real run contract lacks its input-release binding")
    release_root = Path(str(release.get("release_root", ""))).resolve(strict=True)
    manifest_sha256 = str(release.get("manifest_sha256", ""))
    variant = str(contract.get("variant", ""))
    if not variant:
        raise ValueError("real run contract lacks its variant")
    rows = _read_development_rows(
        release_root,
        manifest_sha256=manifest_sha256,
        max_components=max_components,
    )
    by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_component[row["leakage_component_id"]].append(row)
    representatives = [
        sorted(component_rows, key=lambda row: row["observation_key"])[0]
        for _, component_rows in sorted(by_component.items())
    ]
    representative_components = [
        row["leakage_component_id"] for row in representatives
    ]
    paired_left = [
        _model_window(
            release_root=release_root,
            row=row,
            variant=variant,
            mask_label="paired-left",
            unmasked=False,
        )
        for row in representatives
    ]
    paired_right = [
        _model_window(
            release_root=release_root,
            row=row,
            variant=variant,
            mask_label="paired-right",
            unmasked=False,
        )
        for row in representatives
    ]
    unmasked_visits = [
        _model_window(
            release_root=release_root,
            row=row,
            variant=variant,
            mask_label="unmasked",
            unmasked=True,
        )
        for row in rows
    ]
    huber_delta = float(contract["optimization"]["huber_delta"])
    left_embeddings, reconstruction = _encode(
        model=model,
        samples=paired_left,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=True,
        huber_delta=huber_delta,
    )
    right_embeddings, _ = _encode(
        model=model,
        samples=paired_right,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    visit_embeddings, _ = _encode(
        model=model,
        samples=unmasked_visits,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )

    from .model import TWIRLFM0

    torch.manual_seed(int(random_control_seed))
    random_model = TWIRLFM0(model.config).to(target_device)
    random_model.eval()
    random_left, _ = _encode(
        model=random_model,
        samples=paired_left,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    random_right, _ = _encode(
        model=random_model,
        samples=paired_right,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    random_visits, _ = _encode(
        model=random_model,
        samples=unmasked_visits,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    component_ids = [row["leakage_component_id"] for row in rows]
    return {
        "schema_version": REPRESENTATION_HEALTH_SCHEMA_VERSION,
        "passed": True,
        "claim_limit": (
            "development-only representation-health evidence; no morphology, "
            "injection-retention, sealed-test, or foundation-model claim"
        ),
        "run": {
            "run_dir": str(root),
            "variant": variant,
            "architecture": contract.get("architecture"),
            "global_step": int(run_validation["global_step"]),
            "checkpoint_sha256": run_validation["artifact_sha256"]["checkpoint.pt"],
        },
        "evaluation_population": {
            "source_partition": "poc_development",
            "selected_leakage_components": int(len(by_component)),
            "selected_observation_visits": int(len(rows)),
            "selection": "stable SHA-256 ordering of leakage components",
            "one_deterministic_window_per_selected_visit": True,
            "model_visible_identifiers": False,
        },
        "masked_reconstruction": reconstruction,
        "trained_encoder": {
            "embedding_health": summarize_embedding_matrix(left_embeddings),
            "safe_mask_pair_separation": paired_similarity_summary(
                left_embeddings,
                right_embeddings,
                cluster_ids=representative_components,
            ),
            "same_component_cross_visit_retrieval": same_component_retrieval_summary(
                visit_embeddings, component_ids
            ),
        },
        "random_encoder_control": {
            "seed": int(random_control_seed),
            "embedding_health": summarize_embedding_matrix(random_left),
            "safe_mask_pair_separation": paired_similarity_summary(
                random_left,
                random_right,
                cluster_ids=representative_components,
            ),
            "same_component_cross_visit_retrieval": same_component_retrieval_summary(
                random_visits, component_ids
            ),
        },
        "not_evaluated": {
            "frozen_linear_morphology_probe": "requires separately frozen label registry",
            "injection_event_retention": "requires separately frozen development injection canary",
            "shortcut_probes": "requires separately frozen audit metadata registry",
            "sealed_test": "not opened by this development-only evaluator",
        },
    }
