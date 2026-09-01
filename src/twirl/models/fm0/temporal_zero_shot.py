"""Paired, label-blind temporal zero-shot evaluation for TWIRL-FM0.

This evaluator opens only the frozen S66--S77 development panel and the
already-frozen S56--S64 FM input release.  It evaluates exact step-0 and
step-2000 checkpoints on the same in-memory rows, windows, and masks.  It
does not read labels, event products, candidate products, or sealed shards,
and it never creates an optimizer or performs training.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from .dataset import VIEW_NAMES, variant_view_indices
from .registry import sha256_file
from .representation_health import (
    BASELINE_TIME_BINS,
    DEFAULT_PCA_COMPONENTS,
    _encode,
    _feature_matrix,
    _identity_sha256,
    _load_trusted_model,
    _model_projection_spectrum,
    _model_window,
    _paired_similarity_values,
    _read_bound_manifest_rows,
    _read_pca_training_representatives,
    _stable_seed,
    binned_flux_baseline_features,
    fit_train_pca_baseline,
    fit_train_scalar_standardizer,
    paired_similarity_summary,
    query_sector_excluded_retrieval_summary,
    robust_scalar_baseline_features,
    source_clustered_mean_interval,
    summarize_embedding_matrix,
    transform_train_pca_baseline,
    transform_train_scalar_baseline,
)
from .temporal_panel import (
    BASELINE_SECTORS,
    DEVELOPMENT_PARTITION,
    TEMPORAL_PANEL_FIELDS,
    TEMPORAL_PANEL_READY_STATE,
    TEMPORAL_PANEL_RECEIPT_SCHEMA_VERSION,
    TEMPORAL_PANEL_SCHEMA_VERSION,
    TEMPORAL_PANEL_SECTOR_BINDING_FIELDS,
    TEMPORAL_PANEL_SECTOR_BINDING_SCHEMA_VERSION,
    TEMPORAL_PANEL_SECTORS,
)
from .validation import read_json

TEMPORAL_ZERO_SHOT_SCHEMA_VERSION = "twirl_fm0_temporal_zero_shot_v1"
TEMPORAL_ZERO_SHOT_SELECTION_SALT = "twirl_fm0_temporal_zero_shot_v1:bounded-components"
EXPECTED_CHECKPOINT_STEPS = (0, 2_000)
TEMPORAL_COHORTS = ("repeated", "new")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _exact_sha256(value: str, *, label: str) -> str:
    digest = str(value).strip().lower()
    if _SHA256.fullmatch(digest) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return digest


def _exact_int(value: Any, *, label: str, minimum: int = 1) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{label} must be an integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be an integer") from exc
    if str(value).strip() != str(parsed) or parsed < minimum:
        raise ValueError(f"{label} is outside its allowed range")
    return parsed


def _serialized_bool(value: Any, *, label: str) -> bool:
    text = str(value).strip().lower()
    if text == "true":
        return True
    if text == "false":
        return False
    raise ValueError(f"{label} must be serialized as true or false")


def _read_exact_csv(
    path: Path,
    fields: Sequence[str],
    *,
    label: str,
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ValueError(f"{label} schema differs from its frozen contract")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError(f"{label} is empty")
    return rows


def _require_read_only(path: Path, *, label: str) -> None:
    if path.stat().st_mode & 0o222:
        raise ValueError(f"{label} must be read-only")


def load_temporal_panel(
    panel_dir: str | Path,
    *,
    receipt_sha256: str,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Validate and load the immutable identity-only temporal panel."""

    raw_root = Path(panel_dir).expanduser()
    if raw_root.is_symlink():
        raise ValueError("temporal-panel root must be materialized")
    root = raw_root.resolve(strict=True)
    expected_receipt = _exact_sha256(
        receipt_sha256, label="temporal-panel receipt SHA-256"
    )
    if not root.is_dir():
        raise ValueError("temporal-panel root must be a materialized directory")
    _require_read_only(root, label="temporal-panel root")
    receipt_raw = root / "receipt.json"
    if receipt_raw.is_symlink():
        raise ValueError("temporal-panel receipt must be materialized")
    receipt_path = receipt_raw.resolve(strict=True)
    if receipt_path.parent != root:
        raise ValueError("temporal-panel receipt is not materialized in its root")
    _require_read_only(receipt_path, label="temporal-panel receipt")
    observed_receipt = sha256_file(receipt_path)
    if observed_receipt != expected_receipt:
        raise ValueError("temporal-panel receipt differs from its exact SHA-256")
    ready_raw = root / "READY"
    if ready_raw.is_symlink():
        raise ValueError("temporal-panel READY must be materialized")
    ready_path = ready_raw.resolve(strict=True)
    if ready_path.parent != root:
        raise ValueError("temporal-panel READY is not materialized in its root")
    _require_read_only(ready_path, label="temporal-panel READY")
    if ready_path.read_text(encoding="utf-8").strip() != observed_receipt:
        raise ValueError("temporal-panel READY does not bind its receipt")

    receipt = read_json(receipt_path)
    false_boundaries = (
        "scientific_training_eligible",
        "model_training_authorized",
        "sealed_test_access_authorized",
        "event_retention_authorized",
        "formal_model_gate_authorized",
        "production_model_claim",
        "prospective_test_claim",
    )
    if (
        receipt.get("schema_version") != TEMPORAL_PANEL_RECEIPT_SCHEMA_VERSION
        or receipt.get("ready_state") != TEMPORAL_PANEL_READY_STATE
        or receipt.get("passed") is not True
        or tuple(receipt.get("panel_sectors", ())) != TEMPORAL_PANEL_SECTORS
        or tuple(receipt.get("baseline_manifest_sectors", ())) != BASELINE_SECTORS
        or receipt.get("excluded_sectors") != [65]
        or receipt.get("label_blind") is not True
        or receipt.get("identity_only") is not True
        or receipt.get("shard_payloads_opened") is not False
        or receipt.get("temporal_panel_frozen") is not True
        or receipt.get("development_evaluation_eligible") is not True
        or any(receipt.get(field) is not False for field in false_boundaries)
    ):
        raise ValueError("temporal-panel receipt boundary differs")

    panel_name = str(receipt.get("panel_path", ""))
    bindings_name = str(receipt.get("sector_bindings_path", ""))
    if panel_name != "panel.csv" or bindings_name != "sector_bindings.csv":
        raise ValueError("temporal-panel artifact paths differ")
    panel_raw = root / panel_name
    bindings_raw = root / bindings_name
    if panel_raw.is_symlink() or bindings_raw.is_symlink():
        raise ValueError("temporal-panel artifacts must be materialized")
    panel_path = panel_raw.resolve(strict=True)
    bindings_path = bindings_raw.resolve(strict=True)
    if panel_path.parent != root or bindings_path.parent != root:
        raise ValueError("temporal-panel artifact escaped its root")
    _require_read_only(panel_path, label="temporal-panel CSV")
    _require_read_only(bindings_path, label="temporal-panel sector bindings")
    panel_hash = sha256_file(panel_path)
    bindings_hash = sha256_file(bindings_path)
    if panel_hash != _exact_sha256(
        receipt.get("panel_sha256", ""), label="temporal-panel CSV SHA-256"
    ) or bindings_hash != _exact_sha256(
        receipt.get("sector_bindings_sha256", ""),
        label="temporal-panel sector-binding SHA-256",
    ):
        raise ValueError("temporal-panel artifact hash differs")

    rows = _read_exact_csv(panel_path, TEMPORAL_PANEL_FIELDS, label="temporal panel")
    bindings = _read_exact_csv(
        bindings_path,
        TEMPORAL_PANEL_SECTOR_BINDING_FIELDS,
        label="temporal-panel sector bindings",
    )
    if tuple(_exact_int(row["sector"], label="sector binding") for row in bindings) != (
        TEMPORAL_PANEL_SECTORS
    ):
        raise ValueError("temporal-panel sector bindings are not exactly S66--S77")
    release_roots: dict[int, str] = {}
    receipt_hashes: dict[int, str] = {}
    manifest_hashes: dict[int, str] = {}
    expected_development_rows: dict[int, int] = {}
    for binding in bindings:
        sector = int(binding["sector"])
        if (
            binding["binding_schema_version"]
            != TEMPORAL_PANEL_SECTOR_BINDING_SCHEMA_VERSION
            or _exact_int(
                binding["n_sealed_rows_emitted"], label="sealed rows", minimum=0
            )
            != 0
            or _serialized_bool(
                binding["scientific_training_eligible"],
                label="sector scientific-training eligibility",
            )
            is not False
            or _serialized_bool(
                binding["development_evaluation_eligible"],
                label="sector development-evaluation eligibility",
            )
            is not True
        ):
            raise ValueError("temporal-panel sector binding boundary differs")
        release_roots[sector] = str(
            Path(binding["sector_release_root"]).expanduser().resolve(strict=True)
        )
        receipt_hashes[sector] = _exact_sha256(
            binding["sector_receipt_sha256"], label="sector receipt SHA-256"
        )
        manifest_hashes[sector] = _exact_sha256(
            binding["sector_manifest_sha256"], label="sector manifest SHA-256"
        )
        expected_development_rows[sector] = _exact_int(
            binding["n_development_rows"],
            label="sector development row count",
            minimum=0,
        )

    observations: set[str] = set()
    components_to_cohort: dict[str, str] = {}
    cohort_counts: dict[str, int] = defaultdict(int)
    sector_counts: dict[int, int] = defaultdict(int)
    for row in rows:
        observation = str(row["observation_key"]).strip()
        component = str(row["leakage_component_id"]).strip()
        cohort = str(row["temporal_cohort"]).strip()
        sector = _exact_int(row["sector"], label="panel sector")
        if (
            row["panel_schema_version"] != TEMPORAL_PANEL_SCHEMA_VERSION
            or sector not in TEMPORAL_PANEL_SECTORS
            or not observation
            or not component
            or cohort not in TEMPORAL_COHORTS
            or row["source_partition"] != DEVELOPMENT_PARTITION
            or _serialized_bool(row["full_visit_shard"], label="full_visit_shard")
            is not True
            or _serialized_bool(
                row["model_context_length_bound"],
                label="model_context_length_bound",
            )
            is not False
            or _serialized_bool(
                row["scientific_training_eligible"],
                label="scientific_training_eligible",
            )
            is not False
            or _serialized_bool(
                row["development_evaluation_eligible"],
                label="development_evaluation_eligible",
            )
            is not True
            or row["sector_release_root"] != release_roots[sector]
            or row["sector_receipt_sha256"] != receipt_hashes[sector]
            or row["sector_manifest_sha256"] != manifest_hashes[sector]
        ):
            raise ValueError("temporal-panel row boundary differs")
        if observation in observations:
            raise ValueError("temporal panel contains duplicate observations")
        observations.add(observation)
        prior = components_to_cohort.setdefault(component, cohort)
        if prior != cohort:
            raise ValueError("one component crosses temporal cohorts")
        _exact_sha256(row["sha256"], label="panel shard SHA-256")
        cohort_counts[cohort] += 1
        sector_counts[sector] += 1

    if len(rows) != _exact_int(receipt.get("n_panel_rows"), label="panel row count"):
        raise ValueError("temporal-panel row count differs from its receipt")
    if any(
        sector_counts[sector] != expected_development_rows[sector]
        for sector in TEMPORAL_PANEL_SECTORS
    ):
        raise ValueError("temporal-panel rows do not close to sector bindings")
    if set(cohort_counts) != set(TEMPORAL_COHORTS):
        raise ValueError("temporal panel must contain repeated and new cohorts")
    if cohort_counts["repeated"] != int(
        receipt.get("n_repeated_rows", -1)
    ) or cohort_counts["new"] != int(receipt.get("n_new_rows", -1)):
        raise ValueError("temporal-panel cohort counts differ from its receipt")
    component_counts = {
        cohort: sum(value == cohort for value in components_to_cohort.values())
        for cohort in TEMPORAL_COHORTS
    }
    if (
        len(components_to_cohort) != int(receipt.get("n_panel_components", -1))
        or component_counts["repeated"] != int(receipt.get("n_repeated_components", -1))
        or component_counts["new"] != int(receipt.get("n_new_components", -1))
        or int(receipt.get("n_sealed_rows_emitted", -1)) != 0
    ):
        raise ValueError("temporal-panel component or sealed count differs")
    return rows, {
        "root": str(root),
        "receipt_path": str(receipt_path),
        "receipt_sha256": observed_receipt,
        "ready_state": TEMPORAL_PANEL_READY_STATE,
        "panel_path": str(panel_path),
        "panel_sha256": panel_hash,
        "sector_bindings_path": str(bindings_path),
        "sector_bindings_sha256": bindings_hash,
        "config_path": receipt.get("config_path"),
        "config_sha256": receipt.get("config_sha256"),
        "admission_receipt_path": receipt.get("admission_receipt_path"),
        "admission_receipt_sha256": receipt.get("admission_receipt_sha256"),
        "baseline_manifest_path": receipt.get("baseline_manifest_path"),
        "baseline_manifest_sha256": receipt.get("baseline_manifest_sha256"),
        "n_panel_rows": len(rows),
        "n_panel_components": len(components_to_cohort),
        "cohort_row_counts": dict(sorted(cohort_counts.items())),
        "cohort_component_counts": component_counts,
        "label_blind": True,
        "identity_only": True,
        "sealed_rows": 0,
    }


def select_temporal_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    max_repeated_components: int,
    max_new_components: int,
    required_view_indices: Sequence[int] = (),
) -> dict[str, list[dict[str, str]]]:
    """Select bounded components whose every retained visit has required views."""

    limits = {
        "repeated": int(max_repeated_components),
        "new": int(max_new_components),
    }
    if any(value < 2 for value in limits.values()):
        raise ValueError("each temporal cohort requires at least two components")
    required_views = tuple(int(value) for value in required_view_indices)
    if any(value < 0 for value in required_views):
        raise ValueError("required temporal-panel view indices must be nonnegative")
    grouped: dict[str, dict[str, list[dict[str, str]]]] = {
        cohort: defaultdict(list) for cohort in TEMPORAL_COHORTS
    }
    seen_observations: set[str] = set()
    component_cohorts: dict[str, str] = {}
    for source_row in rows:
        row = dict(source_row)
        cohort = str(row.get("temporal_cohort", ""))
        component = str(row.get("leakage_component_id", ""))
        observation = str(row.get("observation_key", ""))
        if (
            cohort not in TEMPORAL_COHORTS
            or not component
            or not observation
            or row.get("source_partition") != DEVELOPMENT_PARTITION
        ):
            raise ValueError("temporal selection received an ineligible row")
        if observation in seen_observations:
            raise ValueError("temporal selection received duplicate observations")
        seen_observations.add(observation)
        prior = component_cohorts.setdefault(component, cohort)
        if prior != cohort:
            raise ValueError("temporal selection component crosses cohorts")
        grouped[cohort][component].append(row)

    selected: dict[str, list[dict[str, str]]] = {}
    for cohort in TEMPORAL_COHORTS:
        components: list[str] = []
        for component, component_rows in grouped[cohort].items():
            eligible = True
            if required_views:
                for row in component_rows:
                    try:
                        present = json.loads(row.get("view_present_json", "[]"))
                    except json.JSONDecodeError as exc:
                        raise ValueError(
                            "temporal selection received invalid view_present_json"
                        ) from exc
                    if (
                        not isinstance(present, list)
                        or len(present) != len(VIEW_NAMES)
                        or any(
                            type(value) is not int or value not in {0, 1}
                            for value in present
                        )
                    ):
                        raise ValueError(
                            "temporal selection received invalid view_present_json schema"
                        )
                    if any(not bool(present[index]) for index in required_views):
                        eligible = False
            if eligible:
                components.append(component)
        if len(components) < 2:
            raise ValueError(f"temporal cohort {cohort} has fewer than two components")
        ordered = sorted(
            components,
            key=lambda component: (
                hashlib.sha256(
                    f"{TEMPORAL_ZERO_SHOT_SELECTION_SALT}\x1f{cohort}\x1f{component}".encode()
                ).hexdigest(),
                component,
            ),
        )
        chosen = set(ordered[: min(limits[cohort], len(ordered))])
        selected[cohort] = sorted(
            (row for component in chosen for row in grouped[cohort][component]),
            key=lambda row: (int(row["sector"]), row["observation_key"]),
        )
    return selected


def _representatives(rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    for source_row in rows:
        row = dict(source_row)
        by_component[row["leakage_component_id"]].append(row)
    return [
        min(component_rows, key=lambda row: row["observation_key"])
        for _, component_rows in sorted(by_component.items())
    ]


def _sample_set_sha256(
    rows: Sequence[Mapping[str, str]],
    samples: Sequence[Mapping[str, np.ndarray]],
) -> str:
    if len(rows) != len(samples):
        raise ValueError("sample-set rows and tensors do not align")
    digest = hashlib.sha256()
    for row, sample in zip(rows, samples, strict=True):
        digest.update(str(row["observation_key"]).encode("utf-8"))
        digest.update(b"\0")
        for key in sorted(sample):
            array = np.ascontiguousarray(np.asarray(sample[key]))
            digest.update(key.encode("utf-8"))
            digest.update(b"\0")
            digest.update(array.dtype.str.encode("ascii"))
            digest.update(b"\0")
            digest.update(repr(array.shape).encode("ascii"))
            digest.update(b"\0")
            digest.update(array.tobytes())
    return digest.hexdigest()


def _build_model_samples(
    rows: Sequence[Mapping[str, str]],
    *,
    variant: str,
    mask_label: str,
    unmasked: bool,
) -> list[dict[str, np.ndarray]]:
    return [
        _model_window(
            release_root=Path(row["sector_release_root"]).resolve(strict=True),
            row=row,
            variant=variant,
            mask_label=mask_label,
            unmasked=unmasked,
        )
        for row in rows
    ]


def _retrieval_hits(
    embeddings: np.ndarray,
    component_ids: Sequence[str],
    sectors: Sequence[int],
) -> dict[int, float]:
    """Return per-query sector-excluded hits for paired checkpoint deltas."""

    array = np.asarray(embeddings, dtype=np.float64)
    groups = tuple(str(value) for value in component_ids)
    sector_values = np.asarray(sectors, dtype=np.int64)
    if (
        array.ndim != 2
        or array.shape[0] != len(groups)
        or sector_values.shape != (array.shape[0],)
        or not np.all(np.isfinite(array))
    ):
        raise ValueError("paired retrieval inputs do not align")
    component_sectors: dict[str, set[int]] = defaultdict(set)
    for component, sector in zip(groups, sector_values, strict=True):
        component_sectors[component].add(int(sector))
    norms = np.linalg.norm(array, axis=1)
    normalized = np.divide(
        array,
        norms[:, None],
        out=np.zeros_like(array),
        where=norms[:, None] > 0.0,
    )
    hits: dict[int, float] = {}
    for query, component in enumerate(groups):
        if len(component_sectors[component]) < 2 or norms[query] <= 0.0:
            continue
        candidates = (sector_values != sector_values[query]) & (norms > 0.0)
        if not np.any(candidates):
            continue
        indices = np.flatnonzero(candidates)
        nearest = int(indices[np.argmax(normalized[indices] @ normalized[query])])
        hits[query] = float(groups[nearest] == component)
    return hits


def _paired_separation_delta(
    step0_left: np.ndarray,
    step0_right: np.ndarray,
    step2000_left: np.ndarray,
    step2000_right: np.ndarray,
    *,
    component_ids: Sequence[str],
) -> dict[str, Any]:
    _, _, initial, initial_valid = _paired_similarity_values(step0_left, step0_right)
    _, _, trained, trained_valid = _paired_similarity_values(
        step2000_left, step2000_right
    )
    identifiers = tuple(str(value) for value in component_ids)
    if len(identifiers) != initial.shape[0]:
        raise ValueError("paired separation identifiers do not align")
    valid = initial_valid & trained_valid
    if not np.any(valid):
        raise ValueError("checkpoint separation delta has no finite pairs")
    values = trained[valid] - initial[valid]
    valid_ids = tuple(
        component
        for component, include in zip(identifiers, valid, strict=True)
        if include
    )
    return {
        "n_paired_source_rows": int(values.size),
        "step2000_minus_step0_mean": float(np.mean(values)),
        "step2000_minus_step0_median": float(np.median(values)),
        "step2000_minus_step0_source_clustered_95_interval": (
            source_clustered_mean_interval(
                values,
                valid_ids,
                seed=_stable_seed("fm0-temporal-checkpoint-separation-delta"),
            )
        ),
    }


def _paired_retrieval_delta(
    step0: np.ndarray,
    step2000: np.ndarray,
    *,
    component_ids: Sequence[str],
    sectors: Sequence[int],
) -> dict[str, Any]:
    initial = _retrieval_hits(step0, component_ids, sectors)
    trained = _retrieval_hits(step2000, component_ids, sectors)
    queries = sorted(set(initial) & set(trained))
    if not queries:
        return {
            "status": "not_available",
            "reason": "no_queries_evaluable_at_both_checkpoints",
            "query_sector_excluded": True,
        }
    delta = np.asarray([trained[index] - initial[index] for index in queries])
    identifiers = tuple(str(component_ids[index]) for index in queries)
    return {
        "status": "available",
        "query_sector_excluded": True,
        "n_paired_queries": len(queries),
        "step0_top1_same_component_retrieval": float(
            np.mean([initial[index] for index in queries])
        ),
        "step2000_top1_same_component_retrieval": float(
            np.mean([trained[index] for index in queries])
        ),
        "step2000_minus_step0_mean": float(np.mean(delta)),
        "step2000_minus_step0_source_clustered_95_interval": (
            source_clustered_mean_interval(
                delta,
                identifiers,
                seed=_stable_seed("fm0-temporal-checkpoint-retrieval-delta"),
            )
        ),
    }


def _representation_summary(
    left: Mapping[str, np.ndarray],
    right: Mapping[str, np.ndarray],
    visits: Mapping[str, np.ndarray],
    *,
    representative_components: Sequence[str],
    visit_components: Sequence[str],
    visit_sectors: Sequence[int],
) -> dict[str, Any]:
    if (
        set(left) != {"h_window", "z_window"}
        or set(right) != set(left)
        or set(visits) != set(left)
    ):
        raise ValueError("temporal checkpoint representation stages differ")
    return {
        name: {
            "embedding_health": summarize_embedding_matrix(left[name]),
            "safe_mask_pair_separation": paired_similarity_summary(
                left[name], right[name], cluster_ids=representative_components
            ),
            "query_sector_excluded_cross_visit_retrieval": (
                query_sector_excluded_retrieval_summary(
                    visits[name], visit_components, visit_sectors
                )
            ),
        }
        for name in ("h_window", "z_window")
    }


def _baseline_summary(
    left: np.ndarray,
    right: np.ndarray,
    visits: np.ndarray,
    *,
    representative_components: Sequence[str],
    visit_components: Sequence[str],
    visit_sectors: Sequence[int],
) -> dict[str, Any]:
    return {
        "embedding_health": summarize_embedding_matrix(left),
        "safe_mask_pair_separation": paired_similarity_summary(
            left, right, cluster_ids=representative_components
        ),
        "query_sector_excluded_cross_visit_retrieval": (
            query_sector_excluded_retrieval_summary(
                visits, visit_components, visit_sectors
            )
        ),
    }


def paired_checkpoint_delta(
    step0: Mapping[str, Any],
    step2000: Mapping[str, Any],
    *,
    representative_components: Sequence[str],
    visit_components: Sequence[str],
    visit_sectors: Sequence[int],
) -> dict[str, Any]:
    """Compute source/query-paired step-2000 minus step-0 diagnostics."""

    output: dict[str, Any] = {"delta_definition": "step2000_minus_step0"}
    representations: dict[str, Any] = {}
    for name in ("h_window", "z_window"):
        initial_health = summarize_embedding_matrix(step0["left"][name])
        trained_health = summarize_embedding_matrix(step2000["left"][name])
        representations[name] = {
            "effective_rank": {
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
                component_ids=representative_components,
            ),
            "query_sector_excluded_cross_visit_retrieval": (
                _paired_retrieval_delta(
                    step0["visits"][name],
                    step2000["visits"][name],
                    component_ids=visit_components,
                    sectors=visit_sectors,
                )
            ),
        }
    output["representations"] = representations

    initial_reconstruction = step0["reconstruction"]
    trained_reconstruction = step2000["reconstruction"]
    if initial_reconstruction["masked_valid_target_count"] != trained_reconstruction[
        "masked_valid_target_count"
    ] or not math.isclose(
        float(initial_reconstruction["zero_prediction_masked_huber_mean"]),
        float(trained_reconstruction["zero_prediction_masked_huber_mean"]),
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        raise ValueError("checkpoint reconstruction samples or masks differ")
    output["masked_reconstruction"] = {
        "masked_valid_target_count": initial_reconstruction[
            "masked_valid_target_count"
        ],
        "step2000_minus_step0_masked_huber_mean": float(
            trained_reconstruction["masked_huber_mean"]
        )
        - float(initial_reconstruction["masked_huber_mean"]),
        "step2000_minus_step0_model_to_zero_baseline_ratio": float(
            trained_reconstruction["model_to_zero_baseline_ratio"]
        )
        - float(initial_reconstruction["model_to_zero_baseline_ratio"]),
        "zero_prediction_masked_huber_mean": initial_reconstruction[
            "zero_prediction_masked_huber_mean"
        ],
    }
    return output


def _validate_baseline_binding(
    contract: Mapping[str, Any],
    *,
    manifest_path: Path,
    manifest_sha256: str,
    panel_summary: Mapping[str, Any],
) -> tuple[Path, str]:
    release = contract.get("input_release")
    if not isinstance(release, Mapping):
        raise TypeError("checkpoint run contract lacks its input release")
    root = manifest_path.parent.resolve(strict=True)
    observed = sha256_file(manifest_path)
    expected = _exact_sha256(manifest_sha256, label="baseline manifest SHA-256")
    contract_hash = _exact_sha256(
        release.get("manifest_sha256", ""),
        label="run-contract input manifest SHA-256",
    )
    panel_hash = _exact_sha256(
        panel_summary.get("baseline_manifest_sha256", ""),
        label="temporal-panel baseline manifest SHA-256",
    )
    contract_root = Path(str(release.get("release_root", ""))).resolve(strict=True)
    panel_path = Path(str(panel_summary.get("baseline_manifest_path", ""))).resolve(
        strict=True
    )
    if (
        manifest_path.name != "manifest.csv"
        or observed != expected
        or expected != contract_hash
        or expected != panel_hash
        or root != contract_root
        or manifest_path != panel_path
    ):
        raise ValueError(
            "baseline PCA/scalar manifest is not the exact frozen S56--S64 authority"
        )
    _read_bound_manifest_rows(root, manifest_sha256=expected)
    return root, observed


def evaluate_temporal_zero_shot(
    *,
    run_dir: str | Path,
    step0_checkpoint_path: str | Path,
    step2000_checkpoint_path: str | Path,
    temporal_panel_dir: str | Path,
    temporal_panel_receipt_sha256: str,
    baseline_manifest_path: str | Path,
    baseline_manifest_sha256: str,
    device: str = "cpu",
    max_repeated_components: int = 256,
    max_new_components: int = 256,
    max_baseline_train_components: int = 512,
    pca_components: int = DEFAULT_PCA_COMPONENTS,
    batch_size: int = 8,
) -> dict[str, Any]:
    """Evaluate exact FM0 step-0/step-2000 checkpoints on S66--S77."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is optional locally
        raise RuntimeError(
            "PyTorch is required for temporal zero-shot evaluation"
        ) from exc
    if max_baseline_train_components < 2 or pca_components <= 0 or batch_size <= 0:
        raise ValueError("temporal evaluator component and batch bounds are invalid")

    root = Path(run_dir).resolve(strict=True)
    target_device = torch.device(device)
    panel_rows, panel_summary = load_temporal_panel(
        temporal_panel_dir, receipt_sha256=temporal_panel_receipt_sha256
    )

    model0, contract0, validation0 = _load_trusted_model(
        root,
        device=target_device,
        checkpoint_path=Path(step0_checkpoint_path),
    )
    if validation0.get("global_step") != EXPECTED_CHECKPOINT_STEPS[0]:
        raise ValueError("initial checkpoint is not exact immutable step 0")
    model2000, contract2000, validation2000 = _load_trusted_model(
        root,
        device=target_device,
        checkpoint_path=Path(step2000_checkpoint_path),
    )
    if validation2000.get("global_step") != EXPECTED_CHECKPOINT_STEPS[1]:
        raise ValueError("trained checkpoint is not exact immutable step 2000")
    if contract0 != contract2000:
        raise ValueError("step-0 and step-2000 checkpoints have different contracts")
    variant = str(contract0.get("variant", ""))
    if not variant:
        raise ValueError("checkpoint run contract lacks its model variant")
    selected = select_temporal_rows(
        panel_rows,
        max_repeated_components=max_repeated_components,
        max_new_components=max_new_components,
        required_view_indices=variant_view_indices(variant),
    )
    huber_delta = float(contract0["optimization"]["huber_delta"])

    baseline_manifest = Path(baseline_manifest_path).resolve(strict=True)
    baseline_root, observed_baseline_hash = _validate_baseline_binding(
        contract0,
        manifest_path=baseline_manifest,
        manifest_sha256=baseline_manifest_sha256,
        panel_summary=panel_summary,
    )
    baseline_rows = _read_pca_training_representatives(
        baseline_root,
        manifest_sha256=observed_baseline_hash,
        max_components=max_baseline_train_components,
    )
    if len(baseline_rows) < 2 or any(
        row["source_partition"] != "poc_train" for row in baseline_rows
    ):
        raise ValueError("baseline fit is not restricted to frozen training rows")
    baseline_samples = [
        _model_window(
            release_root=baseline_root,
            row=row,
            variant=variant,
            mask_label="temporal-baseline-train-unmasked",
            unmasked=True,
        )
        for row in baseline_rows
    ]

    cohort_data: dict[str, dict[str, Any]] = {}
    for cohort in TEMPORAL_COHORTS:
        rows = selected[cohort]
        representatives = _representatives(rows)
        paired_left = _build_model_samples(
            representatives,
            variant=variant,
            mask_label=f"temporal-{cohort}-paired-left",
            unmasked=False,
        )
        paired_right = _build_model_samples(
            representatives,
            variant=variant,
            mask_label=f"temporal-{cohort}-paired-right",
            unmasked=False,
        )
        visits = _build_model_samples(
            rows,
            variant=variant,
            mask_label=f"temporal-{cohort}-unmasked",
            unmasked=True,
        )
        cohort_data[cohort] = {
            "rows": rows,
            "representatives": representatives,
            "left_samples": paired_left,
            "right_samples": paired_right,
            "visit_samples": visits,
            "representative_components": [
                row["leakage_component_id"] for row in representatives
            ],
            "visit_components": [row["leakage_component_id"] for row in rows],
            "visit_sectors": [int(row["sector"]) for row in rows],
            "sample_sha256": {
                "paired_left": _sample_set_sha256(representatives, paired_left),
                "paired_right": _sample_set_sha256(representatives, paired_right),
                "unmasked_visits": _sample_set_sha256(rows, visits),
            },
        }

    raw_checkpoints: dict[int, dict[str, Any]] = {}
    checkpoint_summaries: dict[str, Any] = {}
    for step, model, validation in (
        (0, model0, validation0),
        (2_000, model2000, validation2000),
    ):
        raw_checkpoints[step] = {}
        cohort_summaries: dict[str, Any] = {}
        for cohort in TEMPORAL_COHORTS:
            data = cohort_data[cohort]
            left, reconstruction = _encode(
                model=model,
                samples=data["left_samples"],
                device=target_device,
                batch_size=batch_size,
                with_reconstruction=True,
                huber_delta=huber_delta,
            )
            right, _ = _encode(
                model=model,
                samples=data["right_samples"],
                device=target_device,
                batch_size=batch_size,
                with_reconstruction=False,
                huber_delta=huber_delta,
            )
            visits, _ = _encode(
                model=model,
                samples=data["visit_samples"],
                device=target_device,
                batch_size=batch_size,
                with_reconstruction=False,
                huber_delta=huber_delta,
            )
            if reconstruction is None:
                raise ValueError("temporal checkpoint omitted masked reconstruction")
            raw = {
                "left": left,
                "right": right,
                "visits": visits,
                "reconstruction": reconstruction,
            }
            raw_checkpoints[step][cohort] = raw
            cohort_summaries[cohort] = {
                "sample_set_sha256": data["sample_sha256"],
                "masked_reconstruction": reconstruction,
                "representations": _representation_summary(
                    left,
                    right,
                    visits,
                    representative_components=data["representative_components"],
                    visit_components=data["visit_components"],
                    visit_sectors=data["visit_sectors"],
                ),
            }
        checkpoint_summaries[f"step{step}"] = {
            "global_step": step,
            "checkpoint_path": validation["selected_checkpoint_path"],
            "checkpoint_sha256": validation["selected_checkpoint_sha256"],
            "checkpoint_is_immutable_milestone": validation[
                "selected_checkpoint_is_immutable_milestone"
            ],
            "projection_spectrum": _model_projection_spectrum(model),
            "cohorts": cohort_summaries,
        }

    pca_train = _feature_matrix(baseline_samples, binned_flux_baseline_features)
    pca_fit, pca_metadata = fit_train_pca_baseline(
        pca_train, n_components=pca_components
    )
    scalar_train = _feature_matrix(baseline_samples, robust_scalar_baseline_features)
    scalar_fit, scalar_metadata = fit_train_scalar_standardizer(scalar_train)
    baseline_summaries: dict[str, Any] = {}
    for cohort in TEMPORAL_COHORTS:
        data = cohort_data[cohort]
        baseline_summaries[cohort] = {}
        for name, feature, fit, transform in (
            (
                "pca",
                binned_flux_baseline_features,
                pca_fit,
                transform_train_pca_baseline,
            ),
            (
                "robust_scalar",
                robust_scalar_baseline_features,
                scalar_fit,
                transform_train_scalar_baseline,
            ),
        ):
            left = transform(_feature_matrix(data["left_samples"], feature), fit)
            right = transform(_feature_matrix(data["right_samples"], feature), fit)
            visits = transform(_feature_matrix(data["visit_samples"], feature), fit)
            baseline_summaries[cohort][name] = _baseline_summary(
                left,
                right,
                visits,
                representative_components=data["representative_components"],
                visit_components=data["visit_components"],
                visit_sectors=data["visit_sectors"],
            )

    paired_deltas = {
        cohort: paired_checkpoint_delta(
            raw_checkpoints[0][cohort],
            raw_checkpoints[2_000][cohort],
            representative_components=cohort_data[cohort]["representative_components"],
            visit_components=cohort_data[cohort]["visit_components"],
            visit_sectors=cohort_data[cohort]["visit_sectors"],
        )
        for cohort in TEMPORAL_COHORTS
    }

    population: dict[str, Any] = {}
    for cohort in TEMPORAL_COHORTS:
        data = cohort_data[cohort]
        components = tuple({row["leakage_component_id"] for row in data["rows"]})
        population[cohort] = {
            "selected_leakage_components": len(components),
            "selected_observation_visits": len(data["rows"]),
            "selected_leakage_components_sha256": _identity_sha256(components),
            "selected_observation_keys_sha256": _identity_sha256(
                tuple(row["observation_key"] for row in data["rows"])
            ),
            "all_visits_retained_for_selected_components": True,
            "one_representative_visit_per_component": True,
            "sample_set_sha256": data["sample_sha256"],
        }

    train_components = tuple(row["leakage_component_id"] for row in baseline_rows)
    train_observations = tuple(row["observation_key"] for row in baseline_rows)
    pca_metadata.update(
        {
            "baseline_manifest_sha256": observed_baseline_hash,
            "baseline_manifest_sectors": list(BASELINE_SECTORS),
            "selected_training_components_sha256": _identity_sha256(train_components),
            "selected_training_observation_keys_sha256": _identity_sha256(
                train_observations
            ),
            "sample_set_sha256": _sample_set_sha256(baseline_rows, baseline_samples),
            "feature_definition": (
                f"per-view visible-only means/support in {BASELINE_TIME_BINS} bins"
            ),
        }
    )
    scalar_metadata.update(
        {
            "baseline_manifest_sha256": observed_baseline_hash,
            "baseline_manifest_sectors": list(BASELINE_SECTORS),
            "selected_training_components_sha256": _identity_sha256(train_components),
            "selected_training_observation_keys_sha256": _identity_sha256(
                train_observations
            ),
            "sample_set_sha256": _sample_set_sha256(baseline_rows, baseline_samples),
            "feature_definition": (
                "per-view visible-only robust scalars and support fraction"
            ),
        }
    )

    run_contract_path = root / "run_contract.json"
    return {
        "schema_version": TEMPORAL_ZERO_SHOT_SCHEMA_VERSION,
        "passed": True,
        "passed_definition": (
            "all checksum, pairing, sealed-access, numerical, and output checks "
            "completed; this is descriptive development evidence, not a model gate"
        ),
        "claim_limit": (
            "label-blind S66--S77 temporal development evaluation only; no event, "
            "sealed-test, prospective-test, production, or training claim"
        ),
        "run": {
            "run_dir": str(root),
            "variant": variant,
            "architecture": contract0.get("architecture"),
            "run_git_sha": contract0.get("expected_git_sha"),
            "run_contract_path": str(run_contract_path),
            "run_contract_sha256": sha256_file(run_contract_path),
            "exact_checkpoint_steps": list(EXPECTED_CHECKPOINT_STEPS),
        },
        "temporal_panel": panel_summary,
        "evaluation_population": population,
        "pairing_contract": {
            "step0_and_step2000_use_same_in_memory_sample_objects": True,
            "identical_rows_windows_and_masks_bound_by_sample_sha256": True,
            "selection": (
                "stable SHA-256 component ordering separately within repeated/new; "
                "all visits retained for selected components; every retained visit "
                "declares the checkpoint-required flux views"
            ),
            "query_sector_excluded_retrieval": True,
        },
        "baseline_fit_population": {
            "source_partition": "poc_train",
            "sectors": list(BASELINE_SECTORS),
            "manifest_path": str(baseline_manifest),
            "manifest_sha256": observed_baseline_hash,
            "selected_leakage_components": len(baseline_rows),
            "selected_observation_visits": len(baseline_rows),
            "development_fit_rows": 0,
            "later_sector_fit_rows": 0,
            "sealed_test_fit_rows": 0,
            "pca": pca_metadata,
            "robust_scalar": scalar_metadata,
        },
        "checkpoints": checkpoint_summaries,
        "paired_checkpoint_deltas": paired_deltas,
        "baselines": {
            "fit_only_on_frozen_s56_s64_train_representatives": True,
            "model_visible_identifiers_or_sectors": False,
            "cohorts": baseline_summaries,
        },
        "data_access": {
            "poc_train_window_reads": len(baseline_rows),
            "poc_development_window_reads": sum(
                2 * len(cohort_data[cohort]["representatives"])
                + len(cohort_data[cohort]["rows"])
                for cohort in TEMPORAL_COHORTS
            ),
            "poc_sealed_test_window_reads": 0,
            "sealed_test_access_count": 0,
            "labels_or_events_opened": False,
            "training_or_optimizer_created": False,
        },
        "not_evaluated": {
            "labels_or_morphology": "not authorized by the temporal panel",
            "event_retention": "not authorized by the temporal panel",
            "sealed_test": "not opened by this development-only evaluator",
            "formal_model_gate": "requires a separately frozen decision contract",
            "model_training": "not performed or authorized",
        },
    }


__all__ = [
    "EXPECTED_CHECKPOINT_STEPS",
    "TEMPORAL_COHORTS",
    "TEMPORAL_ZERO_SHOT_SCHEMA_VERSION",
    "evaluate_temporal_zero_shot",
    "load_temporal_panel",
    "paired_checkpoint_delta",
    "select_temporal_rows",
]
