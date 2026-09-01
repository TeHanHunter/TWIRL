"""Fixed-exposure, development-only FM0 context-window diagnostic.

Each selected S66--S77 visit contributes one deterministic, unpadded
2,048-cadence interval.  That same interval and the same two evaluation masks
are tiled without overlap at 256, 512, 1,024, and 2,048 cadences.  Crop
embeddings are averaged back to one row per visit before any metric is
computed, so shorter contexts do not receive extra statistical weight.
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from .dataset import (
    collate_fm0_samples,
    move_batch_to_device,
    prepare_model_window,
    synchronized_temporal_mask,
    variant_view_indices,
)
from .input_release import (
    WindowSpec,
    deterministic_training_window,
    extract_window,
    load_input_release_bytes,
)
from .representation_health import (
    _identity_sha256,
    _load_trusted_model,
    _model_projection_spectrum,
    paired_similarity_summary,
    query_sector_excluded_retrieval_summary,
    summarize_embedding_matrix,
)
from .temporal_panel import DEVELOPMENT_PARTITION
from .temporal_zero_shot import (
    EXPECTED_CHECKPOINT_STEPS,
    TEMPORAL_COHORTS,
    _paired_retrieval_delta,
    _paired_separation_delta,
    _sample_set_sha256,
    load_temporal_panel,
    select_temporal_rows,
)

CONTEXT_WINDOW_DIAGNOSTIC_SCHEMA_VERSION = "twirl_fm0_context_window_diagnostic_v1"
CONTEXT_WINDOW_SELECTION_SALT = "twirl_fm0_context_window_diagnostic_v1"
CONTEXT_LENGTHS = (256, 512, 1_024, 2_048)
FIXED_EXPOSURE_CADENCES = 2_048
MINIMUM_CONTEXT_CADENCES = 256

_CADENCE_KEYS = (
    "flux",
    "flux_valid",
    "flux_error",
    "error_valid",
    "local_time_cadences",
    "delta_time_cadences",
    "time_valid",
    "segment_boundary",
    "temporal_mask",
    "reconstruction_mask",
)
_REPRESENTATIONS = ("h_window", "z_window")


def context_crop_bounds(context_length: int) -> tuple[tuple[int, int], ...]:
    """Return the exact nonoverlapping tiling of the frozen exposure."""

    length = int(context_length)
    if length not in CONTEXT_LENGTHS or FIXED_EXPOSURE_CADENCES % length:
        raise ValueError("context length is outside the frozen diagnostic grid")
    return tuple(
        (start, start + length) for start in range(0, FIXED_EXPOSURE_CADENCES, length)
    )


def _seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def _full_coverage_window_spec(
    release: Any,
    *,
    observation_key: str,
    variant: str,
) -> WindowSpec:
    """Choose an unpadded base with support in every minimum-size tile."""

    initial = deterministic_training_window(
        release,
        observation_key=observation_key,
        epoch=0,
        draw_index=0,
    )
    required = np.asarray(variant_view_indices(variant), dtype=int)
    segments = tuple(int(value) for value in np.unique(release.segment_id))
    ordered_segments = (initial.segment_id,) + tuple(
        value for value in segments if value != initial.segment_id
    )
    candidates: list[tuple[int, int]] = []
    initial_is_eligible = False
    for segment_id in ordered_segments:
        indices = np.flatnonzero(release.segment_id == segment_id)
        if indices.size == 0 or np.any(np.diff(indices) != 1):
            raise ValueError("FM0 release segment is empty or non-contiguous")
        if indices.size < FIXED_EXPOSURE_CADENCES:
            continue
        valid = (
            release.time_valid[indices, None] & release.flux_valid[indices][:, required]
        )
        starts = np.arange(indices.size - FIXED_EXPOSURE_CADENCES + 1)
        cumulative = np.vstack(
            [np.zeros((1, required.size), dtype=np.int64), np.cumsum(valid, axis=0)]
        )
        eligible = np.ones(starts.size, dtype=bool)
        for offset in range(0, FIXED_EXPOSURE_CADENCES, MINIMUM_CONTEXT_CADENCES):
            counts = (
                cumulative[starts + offset + MINIMUM_CONTEXT_CADENCES]
                - cumulative[starts + offset]
            )
            eligible &= np.all(counts > 0, axis=1)
        eligible_starts = starts[eligible]
        if segment_id == initial.segment_id and bool(
            np.any(eligible_starts == initial.start_offset)
        ):
            initial_is_eligible = True
        candidates.extend((segment_id, int(start)) for start in eligible_starts)

    if initial_is_eligible:
        selected = (initial.segment_id, initial.start_offset)
    elif candidates:
        selected = candidates[
            _seed(
                CONTEXT_WINDOW_SELECTION_SALT,
                "base-interval",
                observation_key,
                tuple(int(value) for value in required),
            )
            % len(candidates)
        ]
    else:
        raise ValueError(
            "development visit has no unpadded 2048-cadence interval with "
            "required-view support in every 256-cadence tile: "
            f"{observation_key}"
        )
    return WindowSpec(
        segment_id=int(selected[0]),
        start_offset=int(selected[1]),
        n_observed=FIXED_EXPOSURE_CADENCES,
        n_padded=0,
    )


def _masked_copy(
    base: Mapping[str, np.ndarray],
    *,
    observation_key: str,
    label: str,
) -> dict[str, np.ndarray]:
    sample = {name: np.asarray(value).copy() for name, value in base.items()}
    temporal, reconstruction = synchronized_temporal_mask(
        sample["flux_valid"],
        sample["time_valid"],
        seed=_seed(
            CONTEXT_WINDOW_SELECTION_SALT,
            "paired-mask",
            label,
            observation_key,
        ),
    )
    sample["temporal_mask"] = temporal
    sample["reconstruction_mask"] = reconstruction
    return sample


def _load_visit_base(
    row: Mapping[str, str],
    *,
    variant: str,
) -> dict[str, Any]:
    if row.get("source_partition") != DEVELOPMENT_PARTITION:
        raise ValueError("context diagnostic may open only development rows")
    root = Path(str(row["sector_release_root"])).expanduser().resolve(strict=True)
    relative = Path(str(row["relative_path"]))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("context diagnostic received an unsafe shard path")
    shard = (root / relative).resolve(strict=True)
    if root not in shard.parents:
        raise ValueError("context diagnostic shard escaped its release root")
    payload = shard.read_bytes()
    observed_hash = hashlib.sha256(payload).hexdigest()
    if observed_hash != str(row["sha256"]):
        raise ValueError("context diagnostic shard differs from its panel binding")
    release = load_input_release_bytes(payload)
    spec = _full_coverage_window_spec(
        release,
        observation_key=str(row["observation_key"]),
        variant=variant,
    )
    if spec.n_observed != FIXED_EXPOSURE_CADENCES or spec.n_padded != 0:
        raise ValueError("context diagnostic base interval is padded")
    window = extract_window(
        release,
        segment_id=spec.segment_id,
        start_offset=spec.start_offset,
    )
    unmasked = prepare_model_window(
        window,
        variant=variant,
        mask_seed=_seed(
            CONTEXT_WINDOW_SELECTION_SALT,
            "base-window",
            row["observation_key"],
        ),
        window_length=FIXED_EXPOSURE_CADENCES,
    )
    unmasked["temporal_mask"] = np.zeros(FIXED_EXPOSURE_CADENCES, dtype=bool)
    unmasked["reconstruction_mask"] = np.zeros_like(
        unmasked["reconstruction_mask"], dtype=bool
    )
    left = _masked_copy(
        unmasked,
        observation_key=str(row["observation_key"]),
        label="left",
    )
    right = _masked_copy(
        unmasked,
        observation_key=str(row["observation_key"]),
        label="right",
    )
    minimum_support = []
    for start, stop in context_crop_bounds(MINIMUM_CONTEXT_CADENCES):
        effective_flux_valid = (
            unmasked["flux_valid"][:, start:stop]
            & unmasked["time_valid"][None, start:stop]
        )
        minimum_support.append(
            {
                "start": start,
                "stop": stop,
                "n_time_valid": int(
                    np.count_nonzero(unmasked["time_valid"][start:stop])
                ),
                "n_flux_valid_by_view": [
                    int(value)
                    for value in np.count_nonzero(effective_flux_valid, axis=1)
                ],
            }
        )
    return {
        "unmasked": unmasked,
        "left": left,
        "right": right,
        "binding": {
            "observation_key": str(row["observation_key"]),
            "leakage_component_id": str(row["leakage_component_id"]),
            "sector": int(row["sector"]),
            "shard_sha256": observed_hash,
            "segment_id": int(spec.segment_id),
            "start_offset": int(spec.start_offset),
            "n_observed": int(spec.n_observed),
            "n_padded": int(spec.n_padded),
            "minimum_tile_support": minimum_support,
        },
    }


def partition_base_sample(
    sample: Mapping[str, np.ndarray],
    *,
    context_length: int,
) -> list[dict[str, np.ndarray]]:
    """Slice one prepared base into an exact tiling of model-ready crops."""

    bounds = context_crop_bounds(context_length)
    missing = [name for name in (*_CADENCE_KEYS, "view_present") if name not in sample]
    if missing:
        raise ValueError(f"base sample lacks required fields: {missing}")
    for name in _CADENCE_KEYS:
        if np.asarray(sample[name]).shape[-1] != FIXED_EXPOSURE_CADENCES:
            raise ValueError(f"base sample field {name} is not 2048 cadences")

    crops: list[dict[str, np.ndarray]] = []
    for start, stop in bounds:
        crop: dict[str, np.ndarray] = {}
        for name, value in sample.items():
            array = np.asarray(value)
            crop[name] = (
                array[..., start:stop].copy() if name in _CADENCE_KEYS else array.copy()
            )
        time_valid = np.asarray(crop["time_valid"], dtype=bool)
        local = np.asarray(crop["local_time_cadences"])
        rebased = np.zeros_like(local)
        valid_indices = np.flatnonzero(time_valid)
        if valid_indices.size:
            rebased[time_valid] = local[time_valid] - local[valid_indices[0]]
        crop["local_time_cadences"] = rebased
        crops.append(crop)

    if len(crops) * int(context_length) != FIXED_EXPOSURE_CADENCES:
        raise ValueError("context crops do not close to the fixed exposure")
    for name in _CADENCE_KEYS:
        if name == "local_time_cadences":
            continue
        reconstructed = np.concatenate([crop[name] for crop in crops], axis=-1)
        if not np.array_equal(reconstructed, np.asarray(sample[name])):
            raise ValueError(f"context tiling changed base field {name}")
    return crops


def average_crop_embeddings(
    embeddings: np.ndarray,
    *,
    n_visits: int,
    crops_per_visit: int,
) -> np.ndarray:
    """Average ordered crop embeddings to one equally weighted visit row."""

    array = np.asarray(embeddings)
    if (
        n_visits <= 0
        or crops_per_visit <= 0
        or array.ndim != 2
        or array.shape[0] != n_visits * crops_per_visit
        or not np.all(np.isfinite(array))
    ):
        raise ValueError("crop embeddings do not form a rectangular visit grid")
    return np.mean(array.reshape(n_visits, crops_per_visit, array.shape[1]), axis=1)


def _flatten_visit_crops(
    visits: Sequence[Mapping[str, Any]],
    *,
    flavor: str,
    context_length: int,
) -> tuple[list[dict[str, np.ndarray]], int]:
    if flavor not in {"left", "right", "unmasked"}:
        raise ValueError("unknown context diagnostic sample flavor")
    crops_per_visit = FIXED_EXPOSURE_CADENCES // int(context_length)
    flattened: list[dict[str, np.ndarray]] = []
    for visit in visits:
        crops = partition_base_sample(visit[flavor], context_length=context_length)
        if len(crops) != crops_per_visit:
            raise ValueError("context diagnostic produced a ragged crop grid")
        flattened.extend(crops)
    return flattened, crops_per_visit


def _crop_set_sha256(
    rows: Sequence[Mapping[str, str]],
    samples: Sequence[Mapping[str, np.ndarray]],
    *,
    crops_per_visit: int,
) -> str:
    expanded = [
        {"observation_key": (f"{row['observation_key']}:crop:{crop_index:02d}")}
        for row in rows
        for crop_index in range(crops_per_visit)
    ]
    return _sample_set_sha256(expanded, samples)


def _encode_visit_crops(
    *,
    model: Any,
    samples: Sequence[Mapping[str, np.ndarray]],
    n_visits: int,
    crops_per_visit: int,
    context_length: int,
    device: Any,
    batch_size: int,
) -> dict[str, np.ndarray]:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is optional locally
        raise RuntimeError("PyTorch is required for the context diagnostic") from exc
    if batch_size <= 0 or not samples:
        raise ValueError("context diagnostic batch or samples are empty")
    model.eval()
    collected: dict[str, list[np.ndarray]] = {name: [] for name in _REPRESENTATIONS}
    with torch.no_grad():
        for start in range(0, len(samples), batch_size):
            batch_samples = samples[start : start + batch_size]
            batch = move_batch_to_device(collate_fm0_samples(batch_samples), device)
            output = (
                model(batch)
                if int(context_length) == FIXED_EXPOSURE_CADENCES
                else model.forward_short_context(batch)
            )
            for name in _REPRESENTATIONS:
                values = output[name].detach().float().cpu().numpy()
                if (
                    values.ndim != 2
                    or values.shape[0] != len(batch_samples)
                    or not np.all(np.isfinite(values))
                ):
                    raise ValueError(f"context encoder output {name} is invalid")
                collected[name].append(values)
    return {
        name: average_crop_embeddings(
            np.concatenate(values, axis=0),
            n_visits=n_visits,
            crops_per_visit=crops_per_visit,
        )
        for name, values in collected.items()
    }


def _representation_summary(
    *,
    left: Mapping[str, np.ndarray],
    right: Mapping[str, np.ndarray],
    unmasked: Mapping[str, np.ndarray],
    component_ids: Sequence[str],
    sectors: Sequence[int],
) -> dict[str, Any]:
    if (
        set(left) != set(right)
        or set(left) != set(unmasked)
        or set(left) != set(_REPRESENTATIONS)
    ):
        raise ValueError("context representation stages differ")
    return {
        name: {
            "unmasked_visit_embedding_health": summarize_embedding_matrix(
                unmasked[name]
            ),
            "safe_mask_pair_separation": paired_similarity_summary(
                left[name], right[name], cluster_ids=component_ids
            ),
            "query_sector_excluded_cross_visit_retrieval": (
                query_sector_excluded_retrieval_summary(
                    unmasked[name], component_ids, sectors
                )
            ),
        }
        for name in _REPRESENTATIONS
    }


def _checkpoint_delta(
    *,
    step0: Mapping[str, Mapping[str, np.ndarray]],
    step2000: Mapping[str, Mapping[str, np.ndarray]],
    component_ids: Sequence[str],
    sectors: Sequence[int],
) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for name in _REPRESENTATIONS:
        initial_health = summarize_embedding_matrix(step0["unmasked"][name])
        trained_health = summarize_embedding_matrix(step2000["unmasked"][name])
        output[name] = {
            "unmasked_visit_effective_rank": {
                "step0": initial_health["effective_rank"],
                "step2000": trained_health["effective_rank"],
                "step2000_minus_step0": float(trained_health["effective_rank"])
                - float(initial_health["effective_rank"]),
            },
            "safe_mask_pair_separation": _paired_separation_delta(
                step0["left"][name],
                step0["right"][name],
                step2000["left"][name],
                step2000["right"][name],
                component_ids=component_ids,
            ),
            "query_sector_excluded_cross_visit_retrieval": (
                _paired_retrieval_delta(
                    step0["unmasked"][name],
                    step2000["unmasked"][name],
                    component_ids=component_ids,
                    sectors=sectors,
                )
            ),
        }
    return output


def evaluate_context_window_diagnostic(
    *,
    run_dir: str | Path,
    step0_checkpoint_path: str | Path,
    step2000_checkpoint_path: str | Path,
    temporal_panel_dir: str | Path,
    temporal_panel_receipt_sha256: str,
    max_repeated_components: int = 64,
    max_new_components: int = 64,
    batch_size: int = 32,
) -> dict[str, Any]:
    """Compare exact FM0 checkpoints across a fixed 2,048-cadence exposure."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is optional locally
        raise RuntimeError("PyTorch is required for the context diagnostic") from exc
    if batch_size <= 0:
        raise ValueError("context diagnostic batch size must be positive")
    root = Path(run_dir).resolve(strict=True)
    device = torch.device("cpu")
    panel_rows, panel_summary = load_temporal_panel(
        temporal_panel_dir, receipt_sha256=temporal_panel_receipt_sha256
    )
    model0, contract0, validation0 = _load_trusted_model(
        root,
        device=device,
        checkpoint_path=Path(step0_checkpoint_path),
    )
    model2000, contract2000, validation2000 = _load_trusted_model(
        root,
        device=device,
        checkpoint_path=Path(step2000_checkpoint_path),
    )
    if validation0.get("global_step") != EXPECTED_CHECKPOINT_STEPS[0]:
        raise ValueError("context diagnostic initial checkpoint is not exact step 0")
    if validation2000.get("global_step") != EXPECTED_CHECKPOINT_STEPS[1]:
        raise ValueError("context diagnostic trained checkpoint is not exact step 2000")
    if contract0 != contract2000:
        raise ValueError("context diagnostic checkpoints have different contracts")
    variant = str(contract0.get("variant", ""))
    if not variant:
        raise ValueError("context diagnostic checkpoint lacks its variant")
    selected = select_temporal_rows(
        panel_rows,
        max_repeated_components=max_repeated_components,
        max_new_components=max_new_components,
        required_view_indices=variant_view_indices(variant),
    )
    for model in (model0, model2000):
        if int(model.config.window_length) != FIXED_EXPOSURE_CADENCES or not callable(
            getattr(model, "forward_short_context", None)
        ):
            raise ValueError("checkpoint model lacks exact short-context support")

    cohort_data: dict[str, dict[str, Any]] = {}
    for cohort in TEMPORAL_COHORTS:
        rows = selected[cohort]
        visits = [_load_visit_base(row, variant=variant) for row in rows]
        base_samples = [visit["unmasked"] for visit in visits]
        cohort_data[cohort] = {
            "rows": rows,
            "visits": visits,
            "component_ids": [row["leakage_component_id"] for row in rows],
            "sectors": [int(row["sector"]) for row in rows],
            "base_interval_sha256": _sample_set_sha256(rows, base_samples),
        }

    checkpoint_contexts: dict[int, dict[str, Any]] = {
        step: {} for step in EXPECTED_CHECKPOINT_STEPS
    }
    deltas: dict[str, Any] = {}
    crop_bindings: dict[str, Any] = {}
    for context_length in CONTEXT_LENGTHS:
        context_key = str(context_length)
        deltas[context_key] = {}
        crop_bindings[context_key] = {
            "context_length": context_length,
            "crops_per_visit": FIXED_EXPOSURE_CADENCES // context_length,
            "bounds": [list(value) for value in context_crop_bounds(context_length)],
            "cohorts": {},
        }
        for cohort in TEMPORAL_COHORTS:
            data = cohort_data[cohort]
            samples_by_flavor: dict[str, list[dict[str, np.ndarray]]] = {}
            crops_per_visit = FIXED_EXPOSURE_CADENCES // context_length
            hashes: dict[str, str] = {}
            target_counts: dict[str, int] = {}
            for flavor in ("left", "right", "unmasked"):
                samples, observed_crops = _flatten_visit_crops(
                    data["visits"],
                    flavor=flavor,
                    context_length=context_length,
                )
                if observed_crops != crops_per_visit:
                    raise ValueError("context crop count differs from its grid")
                samples_by_flavor[flavor] = samples
                hashes[flavor] = _crop_set_sha256(
                    data["rows"],
                    samples,
                    crops_per_visit=crops_per_visit,
                )
                target_counts[flavor] = int(
                    sum(
                        np.count_nonzero(sample["reconstruction_mask"])
                        for sample in samples
                    )
                )
            crop_bindings[context_key]["cohorts"][cohort] = {
                "sample_set_sha256": hashes,
                "reconstruction_target_counts": target_counts,
            }

            encoded: dict[int, dict[str, dict[str, np.ndarray]]] = {}
            for step, model in (
                (EXPECTED_CHECKPOINT_STEPS[0], model0),
                (EXPECTED_CHECKPOINT_STEPS[1], model2000),
            ):
                encoded[step] = {
                    flavor: _encode_visit_crops(
                        model=model,
                        samples=samples_by_flavor[flavor],
                        n_visits=len(data["rows"]),
                        crops_per_visit=crops_per_visit,
                        context_length=context_length,
                        device=device,
                        batch_size=batch_size,
                    )
                    for flavor in ("left", "right", "unmasked")
                }
                checkpoint_contexts[step].setdefault(context_key, {})[cohort] = (
                    _representation_summary(
                        left=encoded[step]["left"],
                        right=encoded[step]["right"],
                        unmasked=encoded[step]["unmasked"],
                        component_ids=data["component_ids"],
                        sectors=data["sectors"],
                    )
                )
            deltas[context_key][cohort] = _checkpoint_delta(
                step0=encoded[EXPECTED_CHECKPOINT_STEPS[0]],
                step2000=encoded[EXPECTED_CHECKPOINT_STEPS[1]],
                component_ids=data["component_ids"],
                sectors=data["sectors"],
            )

    for cohort in TEMPORAL_COHORTS:
        for flavor in ("left", "right", "unmasked"):
            counts = {
                int(
                    crop_bindings[str(length)]["cohorts"][cohort][
                        "reconstruction_target_counts"
                    ][flavor]
                )
                for length in CONTEXT_LENGTHS
            }
            if len(counts) != 1:
                raise ValueError("context tilings changed reconstruction-mask exposure")

    population: dict[str, Any] = {}
    base_intervals: dict[str, Any] = {}
    for cohort in TEMPORAL_COHORTS:
        data = cohort_data[cohort]
        components = tuple(sorted({str(value) for value in data["component_ids"]}))
        population[cohort] = {
            "selected_leakage_components": len(components),
            "selected_observation_visits": len(data["rows"]),
            "selected_leakage_components_sha256": _identity_sha256(components),
            "selected_observation_keys_sha256": _identity_sha256(
                tuple(row["observation_key"] for row in data["rows"])
            ),
            "all_visits_retained_for_selected_components": True,
        }
        base_intervals[cohort] = {
            "sample_set_sha256": data["base_interval_sha256"],
            "visits": [visit["binding"] for visit in data["visits"]],
        }

    run_contract_path = root / "run_contract.json"
    return {
        "schema_version": CONTEXT_WINDOW_DIAGNOSTIC_SCHEMA_VERSION,
        "status": "descriptive_development_diagnostic_complete",
        "formal_gate_applied": False,
        "claim_limit": (
            "CPU-only, label-blind S66--S77 development context diagnostic; "
            "not training, event retention, sealed testing, or a model gate"
        ),
        "run": {
            "run_dir": str(root),
            "variant": variant,
            "architecture": contract0.get("architecture"),
            "run_git_sha": contract0.get("expected_git_sha"),
            "run_contract_path": str(run_contract_path),
            "run_contract_sha256": hashlib.sha256(
                run_contract_path.read_bytes()
            ).hexdigest(),
            "device": "cpu",
            "exact_checkpoint_steps": list(EXPECTED_CHECKPOINT_STEPS),
        },
        "temporal_panel": panel_summary,
        "evaluation_population": population,
        "fixed_exposure": {
            "cadences_per_visit": FIXED_EXPOSURE_CADENCES,
            "context_lengths": list(CONTEXT_LENGTHS),
            "crop_counts": {
                str(length): FIXED_EXPOSURE_CADENCES // length
                for length in CONTEXT_LENGTHS
            },
            "nonoverlapping_complete_tilings": True,
            "one_unpadded_single_segment_base_interval_per_visit": True,
            "base_intervals": base_intervals,
            "crop_bindings": crop_bindings,
            "crop_embeddings_averaged_to_visits_before_metrics": True,
            "tile_local_time_rebased": True,
            "reconstruction_mask_exposure_preserved_across_contexts": True,
        },
        "checkpoints": {
            "step0": {
                "global_step": EXPECTED_CHECKPOINT_STEPS[0],
                "checkpoint_path": validation0["selected_checkpoint_path"],
                "checkpoint_sha256": validation0["selected_checkpoint_sha256"],
                "projection_spectrum": _model_projection_spectrum(model0),
                "contexts": checkpoint_contexts[EXPECTED_CHECKPOINT_STEPS[0]],
            },
            "step2000": {
                "global_step": EXPECTED_CHECKPOINT_STEPS[1],
                "checkpoint_path": validation2000["selected_checkpoint_path"],
                "checkpoint_sha256": validation2000["selected_checkpoint_sha256"],
                "projection_spectrum": _model_projection_spectrum(model2000),
                "contexts": checkpoint_contexts[EXPECTED_CHECKPOINT_STEPS[1]],
            },
        },
        "paired_checkpoint_deltas": deltas,
        "pairing_contract": {
            "step0_and_step2000_use_identical_in_memory_crop_objects": True,
            "sources_base_intervals_masks_and_crop_bounds_are_identical": True,
            "all_metrics_have_one_row_per_visit": True,
            "cohorts_reported_separately": list(TEMPORAL_COHORTS),
            "query_sector_excluded_retrieval": True,
        },
        "data_access": {
            "development_shards_opened": sum(
                len(cohort_data[cohort]["rows"]) for cohort in TEMPORAL_COHORTS
            ),
            "each_development_shard_opened_once": True,
            "sealed_shards_opened": 0,
            "labels_candidates_or_events_opened": False,
            "artificial_events_injected": False,
            "optimizer_or_training_created": False,
        },
        "not_evaluated": {
            "event_retention": "not part of this context diagnostic",
            "labels_or_morphology": "not authorized or opened",
            "sealed_test": "not authorized or opened",
            "formal_model_gate": "not applied",
            "model_training": "not performed",
        },
    }


__all__ = [
    "CONTEXT_LENGTHS",
    "CONTEXT_WINDOW_DIAGNOSTIC_SCHEMA_VERSION",
    "FIXED_EXPOSURE_CADENCES",
    "average_crop_embeddings",
    "context_crop_bounds",
    "evaluate_context_window_diagnostic",
    "partition_base_sample",
]
