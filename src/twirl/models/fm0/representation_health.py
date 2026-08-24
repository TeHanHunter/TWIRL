"""Development-only representation-health checks for TWIRL-FM0.

These checks intentionally use only the immutable FM input release.  Learned
encoders are evaluated on ``poc_development``; PCA and robust-scalar controls
are fit on ``poc_train`` and then evaluated on that same development panel.
They measure whether an encoder has a non-collapsed representation and whether
repeated observations of a held-out Gaia component are closer than unrelated
observations.  They do *not* use labels, BLS products, candidate tables,
injection truth, or the sealed source test.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
import csv
import hashlib
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from .dataset import (
    WINDOW_LENGTH,
    collate_fm0_samples,
    move_batch_to_device,
    prepare_model_window,
    variant_view_indices,
)
from .input_release import (
    MANIFEST_COLUMNS,
    ObservationRelease,
    WindowSpec,
    deterministic_training_window,
    extract_window,
    load_input_release,
)
from .registry import sha256_file
from .validation import (
    read_json,
    validate_real_run_release,
)


LEGACY_REPRESENTATION_HEALTH_SCHEMA_VERSION = (
    "twirl_fm0_1_representation_health_v1"
)
REPRESENTATION_HEALTH_SCHEMA_VERSION = "twirl_fm0_representation_health_v2"
OBSERVATION_SECTOR_AUTHORITY_COLUMNS = (
    "observation_key",
    "leakage_component_id",
    "source_partition",
    "sector",
)
BASELINE_TIME_BINS = 64
DEFAULT_PCA_COMPONENTS = 32
TRAINED_RANDOM_BOOTSTRAP_REPLICATES = 2_000
TRAINED_RANDOM_BOOTSTRAP_SEED = 56_006_701


def _stable_seed(*parts: object) -> int:
    digest = hashlib.sha256(
        "\x1f".join(str(part) for part in parts).encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def _finite_float(value: float | np.floating[Any]) -> float | None:
    candidate = float(value)
    return candidate if math.isfinite(candidate) else None


def _identity_sha256(values: Sequence[str]) -> str:
    return hashlib.sha256(
        "\n".join(sorted(str(value) for value in values)).encode("utf-8")
    ).hexdigest()


def load_observation_sector_authority(
    path: str | Path,
    *,
    expected_sha256: str,
) -> tuple[dict[str, dict[str, str]], dict[str, Any]]:
    """Load a checksum-bound evaluator-only observation-to-sector mapping.

    The authority deliberately excludes the sealed partition.  Its fields are
    evaluator metadata and are never added to a model batch.
    """

    source = Path(path).resolve(strict=True)
    expected = str(expected_sha256).strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("observation-sector authority SHA-256 is invalid")
    observed = sha256_file(source)
    if observed != expected:
        raise ValueError("observation-sector authority differs from its SHA-256")
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != OBSERVATION_SECTOR_AUTHORITY_COLUMNS:
            raise ValueError("observation-sector authority schema differs")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("observation-sector authority is empty")
    result: dict[str, dict[str, str]] = {}
    partitions: dict[str, int] = defaultdict(int)
    sectors: set[int] = set()
    for row in rows:
        key = str(row["observation_key"]).strip()
        component = str(row["leakage_component_id"]).strip()
        partition = str(row["source_partition"]).strip()
        if not key or key in result or not component:
            raise ValueError("observation-sector authority has an invalid identity")
        if partition not in {"poc_train", "poc_development"}:
            raise ValueError(
                "observation-sector authority must contain zero sealed-test rows"
            )
        try:
            sector = int(row["sector"])
        except (TypeError, ValueError) as exc:
            raise ValueError("observation-sector authority has an invalid sector") from exc
        if sector <= 0 or str(sector) != str(row["sector"]).strip():
            raise ValueError("observation-sector authority has an invalid sector")
        result[key] = {
            "observation_key": key,
            "leakage_component_id": component,
            "source_partition": partition,
            "sector": str(sector),
        }
        partitions[partition] += 1
        sectors.add(sector)
    return result, {
        "path": str(source),
        "sha256": observed,
        "schema_columns": list(OBSERVATION_SECTOR_AUTHORITY_COLUMNS),
        "n_rows": len(result),
        "partition_counts": dict(sorted(partitions.items())),
        "sectors": sorted(sectors),
        "sealed_test_rows": 0,
        "model_visible": False,
    }


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


def _paired_similarity_values(
    paired_left: np.ndarray,
    paired_right: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return paired, unrelated, difference, and finite-row masks."""

    left = np.asarray(paired_left, dtype=np.float64)
    right = np.asarray(paired_right, dtype=np.float64)
    if left.ndim != 2 or right.shape != left.shape or left.shape[0] < 2:
        raise ValueError("paired similarity requires at least two equal-shaped rows")
    paired = cosine_similarity_rows(left, right)
    unrelated = cosine_similarity_rows(left, np.roll(right, shift=1, axis=0))
    valid = np.isfinite(paired) & np.isfinite(unrelated)
    difference = np.full(paired.shape, np.nan, dtype=np.float64)
    difference[valid] = paired[valid] - unrelated[valid]
    return paired, unrelated, difference, valid


def paired_similarity_summary(
    paired_left: np.ndarray,
    paired_right: np.ndarray,
    *,
    cluster_ids: Sequence[str],
) -> dict[str, Any]:
    """Compare same-window masked pairs with deterministically unrelated rows."""

    left = np.asarray(paired_left, dtype=np.float64)
    paired, unrelated, difference, valid = _paired_similarity_values(
        left, paired_right
    )
    if not np.any(valid):
        raise ValueError("paired similarity has no finite comparisons")
    identifiers = tuple(str(value) for value in cluster_ids)
    if len(identifiers) != left.shape[0]:
        raise ValueError("paired similarity cluster IDs must align with rows")
    return {
        "n_pairs": int(np.count_nonzero(valid)),
        "paired_cosine_mean": _finite_float(np.mean(paired[valid])),
        "unrelated_cosine_mean": _finite_float(np.mean(unrelated[valid])),
        "paired_minus_unrelated_mean": _finite_float(np.mean(difference[valid])),
        "paired_minus_unrelated_median": _finite_float(
            np.median(difference[valid])
        ),
        "paired_minus_unrelated_source_clustered_95_interval": (
            source_clustered_mean_interval(
                difference[valid],
                tuple(
                    identifier
                    for identifier, include in zip(identifiers, valid)
                    if include
                ),
                seed=_stable_seed("fm0-paired-separation-bootstrap"),
            )
        ),
    }


def source_paired_trained_minus_random_summary(
    trained_left: np.ndarray,
    trained_right: np.ndarray,
    random_left: np.ndarray,
    random_right: np.ndarray,
    *,
    cluster_ids: Sequence[str],
    bootstrap_seed: int = TRAINED_RANDOM_BOOTSTRAP_SEED,
    n_bootstrap: int = TRAINED_RANDOM_BOOTSTRAP_REPLICATES,
) -> dict[str, Any]:
    """Compare trained and random separation on identical source rows."""

    trained_left_array = np.asarray(trained_left, dtype=np.float64)
    arrays = tuple(
        np.asarray(value, dtype=np.float64)
        for value in (trained_right, random_left, random_right)
    )
    if any(value.shape != trained_left_array.shape for value in arrays):
        raise ValueError("trained and random paired embeddings must align")
    _, _, trained_difference, trained_valid = _paired_similarity_values(
        trained_left_array, arrays[0]
    )
    _, _, random_difference, random_valid = _paired_similarity_values(
        arrays[1], arrays[2]
    )
    identifiers = tuple(str(value) for value in cluster_ids)
    if len(identifiers) != trained_left_array.shape[0]:
        raise ValueError("trained-minus-random cluster IDs must align with rows")
    valid = trained_valid & random_valid
    if not np.any(valid):
        raise ValueError("trained-minus-random comparison has no finite rows")
    delta = trained_difference[valid] - random_difference[valid]
    valid_identifiers = tuple(
        identifier for identifier, include in zip(identifiers, valid) if include
    )
    return {
        "n_source_pairs": int(np.count_nonzero(valid)),
        "bootstrap_seed": int(bootstrap_seed),
        "trained_paired_minus_unrelated_mean": _finite_float(
            np.mean(trained_difference[valid])
        ),
        "random_paired_minus_unrelated_mean": _finite_float(
            np.mean(random_difference[valid])
        ),
        "trained_minus_random_mean": _finite_float(np.mean(delta)),
        "trained_minus_random_median": _finite_float(np.median(delta)),
        "trained_minus_random_source_clustered_95_interval": (
            source_clustered_mean_interval(
                delta,
                valid_identifiers,
                seed=int(bootstrap_seed),
                n_bootstrap=int(n_bootstrap),
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


def query_sector_excluded_retrieval_summary(
    embeddings: np.ndarray,
    component_ids: Sequence[str],
    sectors: Sequence[int],
) -> dict[str, Any]:
    """Measure cross-visit retrieval after excluding the query's whole sector.

    A query is eligible only when its physical component has another visit in a
    different sector.  All visits from the query sector are removed from its
    candidate set, including visits of unrelated sources.  Sector and source
    identities remain evaluator-only metadata.
    """

    array = np.asarray(embeddings, dtype=np.float64)
    groups = tuple(str(value) for value in component_ids)
    sector_values = np.asarray(sectors)
    if (
        array.ndim != 2
        or array.shape[0] != len(groups)
        or sector_values.shape != (array.shape[0],)
    ):
        raise ValueError("sector-excluded retrieval inputs must align")
    if not np.all(np.isfinite(array)):
        raise ValueError("sector-excluded retrieval embeddings are non-finite")
    parsed: list[int] = []
    for value in sector_values:
        if isinstance(value, (bool, np.bool_)):
            raise ValueError("sector-excluded retrieval sectors must be integers")
        try:
            numeric = float(value)
            candidate = int(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(
                "sector-excluded retrieval sectors must be integers"
            ) from exc
        if not math.isfinite(numeric) or numeric != float(candidate):
            raise ValueError("sector-excluded retrieval sectors must be integers")
        parsed.append(candidate)
    sector_values = np.asarray(parsed, dtype=np.int64)
    if np.any(sector_values <= 0):
        raise ValueError("sector-excluded retrieval sectors must be positive")
    if array.shape[0] < 2:
        return {
            "status": "not_available",
            "reason": "fewer_than_two_visits",
            "query_sector_excluded": True,
        }

    component_sectors: dict[str, set[int]] = defaultdict(set)
    for component, sector in zip(groups, sector_values):
        component_sectors[component].add(int(sector))
    structurally_eligible = np.asarray(
        [len(component_sectors[component]) >= 2 for component in groups],
        dtype=bool,
    )
    if not np.any(structurally_eligible):
        return {
            "status": "not_available",
            "reason": "no_components_repeated_across_sectors",
            "query_sector_excluded": True,
            "n_visit_embeddings": int(array.shape[0]),
        }

    norms = np.linalg.norm(array, axis=1)
    normalized = np.divide(
        array,
        norms[:, None],
        out=np.zeros_like(array),
        where=norms[:, None] > 0.0,
    )
    recovered: list[float] = []
    evaluated_groups: list[str] = []
    zero_norm_queries = 0
    no_nonzero_candidate_queries = 0
    for query in np.flatnonzero(structurally_eligible):
        if norms[query] <= 0.0:
            zero_norm_queries += 1
            continue
        candidates = (sector_values != sector_values[query]) & (norms > 0.0)
        if not np.any(candidates):
            no_nonzero_candidate_queries += 1
            continue
        candidate_indices = np.flatnonzero(candidates)
        nearest = int(
            candidate_indices[
                np.argmax(normalized[candidate_indices] @ normalized[query])
            ]
        )
        recovered.append(float(groups[query] == groups[nearest]))
        evaluated_groups.append(groups[query])
    if not recovered:
        return {
            "status": "not_available",
            "reason": "no_nonzero_sector_excluded_queries",
            "query_sector_excluded": True,
            "n_visit_embeddings": int(array.shape[0]),
            "n_structurally_eligible_queries": int(
                np.count_nonzero(structurally_eligible)
            ),
            "zero_norm_queries": int(zero_norm_queries),
            "queries_without_nonzero_candidate": int(no_nonzero_candidate_queries),
        }
    recovered_array = np.asarray(recovered, dtype=np.float64)
    return {
        "status": "available",
        "query_sector_excluded": True,
        "candidate_definition": "all nonzero visits outside the query sector",
        "n_visit_embeddings": int(array.shape[0]),
        "n_sectors": int(np.unique(sector_values).size),
        "n_cross_sector_repeated_components": int(
            sum(len(values) >= 2 for values in component_sectors.values())
        ),
        "n_structurally_eligible_queries": int(
            np.count_nonzero(structurally_eligible)
        ),
        "n_evaluated_queries": int(recovered_array.size),
        "zero_norm_queries": int(zero_norm_queries),
        "queries_without_nonzero_candidate": int(no_nonzero_candidate_queries),
        "top1_same_component_retrieval": _finite_float(
            np.mean(recovered_array)
        ),
        "top1_source_clustered_95_interval": source_clustered_mean_interval(
            recovered_array,
            evaluated_groups,
            seed=_stable_seed("fm0-sector-excluded-retrieval-bootstrap"),
        ),
    }


def summarize_projection_spectrum(
    weight: np.ndarray,
    bias: np.ndarray | None = None,
) -> dict[str, Any]:
    """Summarize the full singular spectrum of a linear projection."""

    matrix = np.asarray(weight, dtype=np.float64)
    if matrix.ndim != 2 or min(matrix.shape) < 1:
        raise ValueError("projection weight must be a nonempty matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("projection weight contains non-finite values")
    bias_array = None if bias is None else np.asarray(bias, dtype=np.float64)
    if bias_array is not None:
        if bias_array.shape != (matrix.shape[0],) or not np.all(
            np.isfinite(bias_array)
        ):
            raise ValueError("projection bias is invalid")
    singular = np.linalg.svd(matrix, compute_uv=False)
    maximum = float(singular[0]) if singular.size else 0.0
    tolerance = float(max(matrix.shape) * np.finfo(np.float64).eps * maximum)
    numerical_rank = int(np.count_nonzero(singular > tolerance))
    energy = np.square(singular)
    total = float(np.sum(energy))
    if total > 0.0:
        weights = energy / total
        positive = weights[weights > 0.0]
        effective_rank = float(np.exp(-np.sum(positive * np.log(positive))))
        leading_share = float(np.max(weights))
    else:
        effective_rank = 0.0
        leading_share = 1.0
    retained = singular[singular > tolerance]
    condition = (
        float(retained[0] / retained[-1]) if retained.size else None
    )
    return {
        "input_dim": int(matrix.shape[1]),
        "output_dim": int(matrix.shape[0]),
        "numerical_rank": numerical_rank,
        "rank_tolerance": _finite_float(tolerance),
        "effective_rank": _finite_float(effective_rank),
        "leading_singular_energy_share": _finite_float(leading_share),
        "largest_singular_value": _finite_float(maximum),
        "smallest_singular_value": _finite_float(singular[-1]),
        "condition_number_on_numerical_support": _finite_float(condition)
        if condition is not None
        else None,
        "frobenius_norm": _finite_float(np.linalg.norm(matrix)),
        "bias_l2_norm": _finite_float(np.linalg.norm(bias_array))
        if bias_array is not None
        else None,
        "singular_values": [float(value) for value in singular],
    }


def _read_bound_manifest_rows(
    release_root: Path,
    *,
    manifest_sha256: str,
) -> list[dict[str, str]]:
    manifest = release_root / "manifest.csv"
    if sha256_file(manifest) != manifest_sha256:
        raise ValueError("FM0 release manifest differs from the run binding")
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise ValueError("FM0 release manifest schema differs from the contract")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("FM0 release manifest is empty")
    return rows


def _read_partition_rows(
    release_root: Path,
    *,
    manifest_sha256: str,
    partition: str,
    max_components: int,
    selection_salt: str,
) -> list[dict[str, str]]:
    rows = _read_bound_manifest_rows(
        release_root, manifest_sha256=manifest_sha256
    )
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["source_partition"] != partition:
            continue
        if row["scientific_training_eligible"] != "True":
            raise ValueError(
                f"{partition} selection contains a training-ineligible row"
            )
        grouped[row["leakage_component_id"]].append(row)
    if not grouped:
        raise ValueError(f"FM0 release has no {partition} rows")
    if max_components <= 1:
        raise ValueError("max_components must be at least two")
    ordered_components = sorted(
        grouped,
        key=lambda component: hashlib.sha256(
            f"{selection_salt}:{component}".encode("utf-8")
        ).digest(),
    )
    selected = ordered_components[:max_components]
    result: list[dict[str, str]] = []
    for component in selected:
        result.extend(sorted(grouped[component], key=lambda row: row["observation_key"]))
    return result


def _read_development_rows(
    release_root: Path,
    *,
    manifest_sha256: str,
    max_components: int,
) -> list[dict[str, str]]:
    # Preserve the exact FM0.1 selection salt so evaluator v2 is matched to the
    # already frozen 256-component development panel byte-for-byte.
    return _read_partition_rows(
        release_root,
        manifest_sha256=manifest_sha256,
        partition="poc_development",
        max_components=max_components,
        selection_salt=LEGACY_REPRESENTATION_HEALTH_SCHEMA_VERSION,
    )


def _read_pca_training_representatives(
    release_root: Path,
    *,
    manifest_sha256: str,
    max_components: int,
) -> list[dict[str, str]]:
    rows = _read_partition_rows(
        release_root,
        manifest_sha256=manifest_sha256,
        partition="poc_train",
        max_components=max_components,
        selection_salt="twirl_fm0_representation_health_v2:pca_train",
    )
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["leakage_component_id"]].append(row)
    return [
        sorted(component_rows, key=lambda row: row["observation_key"])[0]
        for _, component_rows in sorted(grouped.items())
    ]


def _validate_authority_against_manifest(
    manifest_rows: Sequence[Mapping[str, str]],
    authority: Mapping[str, Mapping[str, str]],
) -> dict[str, Any]:
    """Require complete authority coverage of nonsealed manifest rows only."""

    allowed = {"poc_train", "poc_development"}
    manifest_nonsealed = {
        str(row["observation_key"]): row
        for row in manifest_rows
        if row["source_partition"] in allowed
    }
    if set(manifest_nonsealed) != set(authority):
        raise ValueError(
            "observation-sector authority does not exactly cover train/development "
            "manifest rows"
        )
    for key, row in manifest_nonsealed.items():
        bound = authority[key]
        if (
            bound["leakage_component_id"] != row["leakage_component_id"]
            or bound["source_partition"] != row["source_partition"]
        ):
            raise ValueError(
                "observation-sector authority identity differs from the manifest"
            )
    return {
        "manifest_train_development_rows": int(len(manifest_nonsealed)),
        "authority_exactly_covers_manifest_train_development": True,
        "sealed_manifest_rows_not_in_authority": int(
            sum(
                row["source_partition"] == "poc_sealed_test"
                for row in manifest_rows
            )
        ),
    }


def _authority_sectors_for_rows(
    rows: Sequence[Mapping[str, str]],
    authority: Mapping[str, Mapping[str, str]],
    *,
    partition: str,
) -> list[int]:
    sectors: list[int] = []
    for row in rows:
        if row["source_partition"] != partition:
            raise ValueError("selected row is in an unexpected source partition")
        key = str(row["observation_key"])
        bound = authority.get(key)
        if bound is None:
            raise ValueError("selected observation is absent from sector authority")
        if (
            bound["leakage_component_id"] != row["leakage_component_id"]
            or bound["source_partition"] != partition
        ):
            raise ValueError("selected observation differs from sector authority")
        sectors.append(int(bound["sector"]))
    return sectors


def _eligible_model_window_spec(
    release: ObservationRelease,
    *,
    observation_key: str,
    variant: str,
) -> WindowSpec:
    """Return a deterministic window containing every required flux view.

    The ordinary data-independent window remains the first choice.  A valid
    visit can, however, contain a quality-masked interval longer than one
    2,048-cadence window.  If the ordinary choice lands entirely inside such
    an interval, search only within the same visit and select a deterministic
    eligible start.  This preserves the frozen component/visit selection and
    never substitutes a more convenient source.
    """

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
    for segment_id in ordered_segments:
        indices = np.flatnonzero(release.segment_id == segment_id)
        if indices.size == 0 or np.any(np.diff(indices) != 1):
            raise ValueError("FM0 release segment is empty or non-contiguous")
        valid = (
            release.time_valid[indices, None]
            & release.flux_valid[indices][:, required]
        )
        initial_start = initial.start_offset if segment_id == initial.segment_id else None
        if initial_start is not None:
            initial_stop = min(indices.size, initial_start + WINDOW_LENGTH)
            if np.all(np.any(valid[initial_start:initial_stop], axis=0)):
                return initial

        starts = np.arange(max(1, indices.size - WINDOW_LENGTH + 1))
        stops = np.minimum(starts + WINDOW_LENGTH, indices.size)
        cumulative = np.vstack(
            [np.zeros((1, required.size), dtype=np.int64), np.cumsum(valid, axis=0)]
        )
        counts = cumulative[stops] - cumulative[starts]
        eligible = starts[np.all(counts > 0, axis=1)]
        if eligible.size:
            offset = int(
                _stable_seed(
                    "fm0-representation-eligible-window",
                    observation_key,
                    tuple(int(value) for value in required),
                    segment_id,
                )
                % eligible.size
            )
            start = int(eligible[offset])
            observed = min(WINDOW_LENGTH, int(indices.size) - start)
            return WindowSpec(
                segment_id=segment_id,
                start_offset=start,
                n_observed=observed,
                n_padded=WINDOW_LENGTH - observed,
            )
    raise ValueError(
        "FM0 development visit has no window with valid cadences for every "
        f"required view: {observation_key}"
    )


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
    spec = _eligible_model_window_spec(
        release,
        observation_key=row["observation_key"],
        variant=variant,
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


def _baseline_visible_flux(
    sample: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    flux = np.asarray(sample["flux"], dtype=np.float64)
    valid = np.asarray(sample["flux_valid"], dtype=bool)
    reconstruction = np.asarray(sample["reconstruction_mask"], dtype=bool)
    present = np.asarray(sample["view_present"], dtype=bool)
    time_valid = np.asarray(sample["time_valid"], dtype=bool)
    if (
        flux.ndim != 2
        or valid.shape != flux.shape
        or reconstruction.shape != flux.shape
        or present.shape != (flux.shape[0],)
        or time_valid.shape != (flux.shape[1],)
    ):
        raise ValueError("FM0 baseline sample tensor shapes are invalid")
    visible = (
        valid
        & ~reconstruction
        & present[:, None]
        & time_valid[None, :]
        & np.isfinite(flux)
    )
    return flux, visible


def binned_flux_baseline_features(
    sample: Mapping[str, np.ndarray],
    *,
    n_time_bins: int = BASELINE_TIME_BINS,
) -> np.ndarray:
    """Return visible-only binned flux means and support fractions for PCA."""

    flux, visible = _baseline_visible_flux(sample)
    if n_time_bins <= 0 or n_time_bins > flux.shape[1]:
        raise ValueError("baseline time-bin count is invalid")
    features: list[float] = []
    for view in range(flux.shape[0]):
        for indices in np.array_split(np.arange(flux.shape[1]), n_time_bins):
            selected = visible[view, indices]
            features.append(
                float(np.mean(flux[view, indices][selected]))
                if np.any(selected)
                else 0.0
            )
        for indices in np.array_split(np.arange(flux.shape[1]), n_time_bins):
            features.append(float(np.mean(visible[view, indices])))
    result = np.asarray(features, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise ValueError("binned PCA baseline produced non-finite features")
    return result


def robust_scalar_baseline_features(
    sample: Mapping[str, np.ndarray],
) -> np.ndarray:
    """Return robust visible-only light-curve summaries for each flux view."""

    flux, visible = _baseline_visible_flux(sample)
    features: list[float] = []
    for view in range(flux.shape[0]):
        values = flux[view, visible[view]]
        if values.size:
            median = float(np.median(values))
            mad = float(np.median(np.abs(values - median)))
            quantiles = np.quantile(values, (0.05, 0.25, 0.75, 0.95))
        else:
            median = 0.0
            mad = 0.0
            quantiles = np.zeros(4, dtype=np.float64)
        adjacent = visible[view, :-1] & visible[view, 1:]
        if np.any(adjacent):
            differences = np.diff(flux[view])[adjacent]
            difference_median = float(np.median(differences))
            difference_mad = float(
                np.median(np.abs(differences - difference_median))
            )
        else:
            difference_mad = 0.0
        features.extend(
            (
                median,
                mad,
                *(float(value) for value in quantiles),
                difference_mad,
                float(np.mean(visible[view])),
            )
        )
    result = np.asarray(features, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise ValueError("robust scalar baseline produced non-finite features")
    return result


def _feature_matrix(
    samples: Sequence[Mapping[str, np.ndarray]],
    feature_function: Any,
) -> np.ndarray:
    if not samples:
        raise ValueError("baseline feature extraction received no samples")
    matrix = np.stack([feature_function(sample) for sample in samples], axis=0)
    if matrix.ndim != 2 or not np.all(np.isfinite(matrix)):
        raise ValueError("baseline feature matrix is invalid")
    return matrix


def _array_sha256(array: np.ndarray) -> str:
    values = np.ascontiguousarray(np.asarray(array, dtype="<f8"))
    header = f"shape={values.shape};dtype=<f8;".encode("ascii")
    return hashlib.sha256(header + values.tobytes()).hexdigest()


def fit_train_pca_baseline(
    train_features: np.ndarray,
    *,
    n_components: int = DEFAULT_PCA_COMPONENTS,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Fit a standardized PCA using training-partition features only."""

    matrix = np.asarray(train_features, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 1:
        raise ValueError("PCA fit requires at least two rows and one feature")
    if not np.all(np.isfinite(matrix)) or n_components <= 0:
        raise ValueError("PCA fit inputs are invalid")
    mean = np.mean(matrix, axis=0)
    raw_scale = np.std(matrix, axis=0, ddof=0)
    scale = np.where(raw_scale > 1.0e-12, raw_scale, 1.0)
    standardized = (matrix - mean) / scale
    _, singular, right = np.linalg.svd(standardized, full_matrices=False)
    retained = min(int(n_components), matrix.shape[0] - 1, matrix.shape[1])
    if retained <= 0:
        raise ValueError("PCA fit has no retainable components")
    components = right[:retained].copy()
    # SVD component signs are algebraically arbitrary.  Fix each sign by its
    # largest-magnitude loading so provenance hashes remain deterministic.
    pivots = np.argmax(np.abs(components), axis=1)
    signs = np.sign(components[np.arange(retained), pivots])
    signs[signs == 0.0] = 1.0
    components *= signs[:, None]
    energy = np.square(singular)
    total_energy = float(np.sum(energy))
    explained = (
        float(np.sum(energy[:retained]) / total_energy)
        if total_energy > 0.0
        else 0.0
    )
    fit = {"mean": mean, "scale": scale, "components": components}
    metadata = {
        "fit_partition": "poc_train",
        "n_fit_rows": int(matrix.shape[0]),
        "input_feature_dim": int(matrix.shape[1]),
        "requested_components": int(n_components),
        "retained_components": int(retained),
        "constant_input_features": int(np.count_nonzero(raw_scale <= 1.0e-12)),
        "explained_variance_fraction": _finite_float(explained),
        "fit_feature_matrix_sha256": _array_sha256(matrix),
        "fit_mean_sha256": _array_sha256(mean),
        "fit_scale_sha256": _array_sha256(scale),
        "fit_components_sha256": _array_sha256(components),
        "development_fit_rows": 0,
        "sealed_test_fit_rows": 0,
    }
    return fit, metadata


def transform_train_pca_baseline(
    features: np.ndarray,
    fit: Mapping[str, np.ndarray],
) -> np.ndarray:
    matrix = np.asarray(features, dtype=np.float64)
    mean = np.asarray(fit["mean"], dtype=np.float64)
    scale = np.asarray(fit["scale"], dtype=np.float64)
    components = np.asarray(fit["components"], dtype=np.float64)
    if (
        matrix.ndim != 2
        or mean.shape != (matrix.shape[1],)
        or scale.shape != mean.shape
        or components.ndim != 2
        or components.shape[1] != matrix.shape[1]
        or np.any(scale <= 0.0)
    ):
        raise ValueError("PCA transform inputs differ from the fit")
    transformed = ((matrix - mean) / scale) @ components.T
    if not np.all(np.isfinite(transformed)):
        raise ValueError("PCA transform produced non-finite features")
    return transformed


def fit_train_scalar_standardizer(
    train_features: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Fit training-only centering/scaling for the robust scalar baseline."""

    matrix = np.asarray(train_features, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 1:
        raise ValueError("scalar baseline fit requires at least two rows")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("scalar baseline fit contains non-finite features")
    median = np.median(matrix, axis=0)
    mad = np.median(np.abs(matrix - median), axis=0)
    scale = np.where(mad > 1.0e-12, 1.4826 * mad, 1.0)
    fit = {"center": median, "scale": scale}
    metadata = {
        "fit_partition": "poc_train",
        "n_fit_rows": int(matrix.shape[0]),
        "feature_dim": int(matrix.shape[1]),
        "constant_or_near_constant_features": int(
            np.count_nonzero(mad <= 1.0e-12)
        ),
        "fit_feature_matrix_sha256": _array_sha256(matrix),
        "fit_center_sha256": _array_sha256(median),
        "fit_scale_sha256": _array_sha256(scale),
        "development_fit_rows": 0,
        "sealed_test_fit_rows": 0,
    }
    return fit, metadata


def transform_train_scalar_baseline(
    features: np.ndarray,
    fit: Mapping[str, np.ndarray],
) -> np.ndarray:
    matrix = np.asarray(features, dtype=np.float64)
    center = np.asarray(fit["center"], dtype=np.float64)
    scale = np.asarray(fit["scale"], dtype=np.float64)
    if (
        matrix.ndim != 2
        or center.shape != (matrix.shape[1],)
        or scale.shape != center.shape
        or np.any(scale <= 0.0)
    ):
        raise ValueError("scalar baseline transform inputs differ from the fit")
    transformed = (matrix - center) / scale
    if not np.all(np.isfinite(transformed)):
        raise ValueError("scalar baseline transform produced non-finite features")
    return transformed


def _encode(
    *,
    model: Any,
    samples: Sequence[Mapping[str, np.ndarray]],
    device: Any,
    batch_size: int,
    with_reconstruction: bool,
    huber_delta: float,
) -> tuple[dict[str, np.ndarray], dict[str, Any] | None]:
    try:
        import torch
        import torch.nn.functional as functional
    except ImportError as exc:  # pragma: no cover - base package keeps Torch optional
        raise RuntimeError("PyTorch is required for FM0 representation evaluation") from exc
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    embeddings: dict[str, list[np.ndarray]] = {
        "h_window": [],
        "z_window": [],
    }
    losses: list[np.ndarray] = []
    baseline_losses: list[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(samples), batch_size):
            batch = collate_fm0_samples(samples[start : start + batch_size])
            batch = move_batch_to_device(batch, device)
            output = model(batch)
            for name in embeddings:
                if name not in output:
                    raise ValueError(f"FM0 encoder output lacks {name}")
                value = output[name].detach().float().cpu().numpy()
                if value.ndim != 2 or value.shape[0] != len(
                    samples[start : start + batch_size]
                ):
                    raise ValueError(f"FM0 encoder output {name} has invalid shape")
                if not np.all(np.isfinite(value)):
                    raise ValueError(f"FM0 encoder output {name} is non-finite")
                embeddings[name].append(value)
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
                baseline_loss = functional.huber_loss(
                    torch.zeros_like(output["reconstruction"]),
                    batch["flux"],
                    delta=float(huber_delta),
                    reduction="none",
                )
                for item, baseline_item, item_mask in zip(loss, baseline_loss, mask):
                    valid = item[item_mask]
                    if valid.numel():
                        losses.append(valid.detach().cpu().numpy())
                        baseline_losses.append(
                            baseline_item[item_mask].detach().cpu().numpy()
                        )
    merged = {
        name: np.concatenate(values, axis=0)
        for name, values in embeddings.items()
    }
    if not with_reconstruction:
        return merged, None
    all_losses = np.concatenate(losses) if losses else np.empty(0, dtype=np.float32)
    all_baseline_losses = (
        np.concatenate(baseline_losses)
        if baseline_losses
        else np.empty(0, dtype=np.float32)
    )
    model_mean = _finite_float(np.mean(all_losses)) if all_losses.size else None
    baseline_mean = (
        _finite_float(np.mean(all_baseline_losses))
        if all_baseline_losses.size
        else None
    )
    return merged, {
        "masked_valid_target_count": int(all_losses.size),
        "masked_huber_mean": model_mean,
        "zero_prediction_masked_huber_mean": baseline_mean,
        "model_to_zero_baseline_ratio": (
            _finite_float(float(model_mean) / float(baseline_mean))
            if model_mean is not None
            and baseline_mean is not None
            and baseline_mean > 0.0
            else None
        ),
    }


def _load_trusted_model(
    run_dir: Path,
    *,
    device: Any,
    checkpoint_path: Path | None = None,
) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - base package keeps Torch optional
        raise RuntimeError("PyTorch is required for FM0 representation evaluation") from exc
    from .model import FM0ModelConfig, TWIRLFM0

    validation = validate_real_run_release(run_dir, inspect_checkpoint=True)
    contract = read_json(run_dir / "run_contract.json")
    selected = (
        run_dir / "checkpoint.pt"
        if checkpoint_path is None
        else Path(checkpoint_path).resolve(strict=True)
    )
    if selected.parent != run_dir:
        raise ValueError("representation checkpoint must be directly inside run_dir")
    if selected.name != "checkpoint.pt":
        match = re.fullmatch(r"checkpoint_step_([0-9]{8})\.pt", selected.name)
        if match is None:
            raise ValueError("representation checkpoint name is not an immutable milestone")
        sidecar = selected.with_name(selected.name + ".sha256")
        if not sidecar.is_file():
            raise ValueError("representation milestone lacks its SHA-256 sidecar")
        fields = sidecar.read_text(encoding="utf-8").strip().split()
        digest = sha256_file(selected)
        if len(fields) != 2 or fields[0] != digest or fields[1] != selected.name:
            raise ValueError("representation milestone SHA-256 sidecar differs")
    else:
        digest = str(validation["artifact_sha256"]["checkpoint.pt"])
    checkpoint = torch.load(selected, map_location=device, weights_only=False)
    if not isinstance(checkpoint, Mapping) or checkpoint.get("run_contract") != contract:
        raise ValueError("representation checkpoint run contract differs")
    progress = checkpoint.get("progress")
    if not isinstance(progress, Mapping):
        raise ValueError("representation checkpoint progress is malformed")
    global_step = progress.get("global_step")
    if isinstance(global_step, bool) or not isinstance(global_step, int) or global_step < 0:
        raise ValueError("representation checkpoint step is invalid")
    if selected.name != "checkpoint.pt":
        filename_step = int(selected.name.removeprefix("checkpoint_step_").removesuffix(".pt"))
        if global_step != filename_step:
            raise ValueError("representation milestone filename and state differ")
        allowed = contract.get("immutable_milestone_steps")
        if not isinstance(allowed, list) or global_step not in allowed:
            raise ValueError("representation milestone was not predeclared")
    model = TWIRLFM0(FM0ModelConfig(**dict(checkpoint["model_config"]))).to(device)
    model.load_state_dict(checkpoint["model_state"], strict=True)
    if any(not torch.all(torch.isfinite(value)) for value in model.state_dict().values()):
        raise ValueError("representation checkpoint contains non-finite model state")
    model.eval()
    selected_validation = dict(validation)
    selected_validation["global_step"] = int(global_step)
    selected_validation["selected_checkpoint_path"] = str(selected)
    selected_validation["selected_checkpoint_sha256"] = digest
    selected_validation["selected_checkpoint_is_immutable_milestone"] = (
        selected.name != "checkpoint.pt"
    )
    return model, contract, selected_validation


def _model_projection_spectrum(model: Any) -> dict[str, Any]:
    layer = getattr(model, "embedding_projection", None)
    if layer is None or not hasattr(layer, "weight"):
        raise ValueError("FM0 model lacks its embedding projection")
    weight = layer.weight.detach().float().cpu().numpy()
    bias = (
        layer.bias.detach().float().cpu().numpy()
        if getattr(layer, "bias", None) is not None
        else None
    )
    return summarize_projection_spectrum(weight, bias)


def _representation_stage_summary(
    left: Mapping[str, np.ndarray],
    right: Mapping[str, np.ndarray],
    visits: Mapping[str, np.ndarray],
    *,
    representative_components: Sequence[str],
    visit_components: Sequence[str],
    visit_sectors: Sequence[int],
) -> dict[str, Any]:
    expected = {"h_window", "z_window"}
    if set(left) != expected or set(right) != expected or set(visits) != expected:
        raise ValueError("FM0 representation stages differ from evaluator v2")
    return {
        name: {
            "embedding_health": summarize_embedding_matrix(left[name]),
            "safe_mask_pair_separation": paired_similarity_summary(
                left[name],
                right[name],
                cluster_ids=representative_components,
            ),
            "query_sector_excluded_cross_visit_retrieval": (
                query_sector_excluded_retrieval_summary(
                    visits[name], visit_components, visit_sectors
                )
            ),
        }
        for name in ("h_window", "z_window")
    }


def _baseline_evaluation_summary(
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


def evaluate_representation_health(
    *,
    run_dir: str | Path,
    checkpoint_path: str | Path | None = None,
    observation_sector_authority_path: str | Path,
    observation_sector_authority_sha256: str,
    device: str = "cpu",
    max_components: int = 256,
    max_pca_train_components: int = 512,
    pca_components: int = DEFAULT_PCA_COMPONENTS,
    batch_size: int = 8,
    random_control_seed: int = 0,
) -> dict[str, Any]:
    """Evaluate one encoder with train-fit baselines and development metrics."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - base package keeps Torch optional
        raise RuntimeError("PyTorch is required for FM0 representation evaluation") from exc
    root = Path(run_dir).resolve(strict=True)
    target_device = torch.device(device)
    selected_checkpoint = (
        None if checkpoint_path is None else Path(checkpoint_path)
    )
    model, contract, run_validation = _load_trusted_model(
        root,
        device=target_device,
        checkpoint_path=selected_checkpoint,
    )
    release = contract.get("input_release")
    if not isinstance(release, Mapping):
        raise ValueError("real run contract lacks its input-release binding")
    release_root = Path(str(release.get("release_root", ""))).resolve(strict=True)
    manifest_sha256 = str(release.get("manifest_sha256", ""))
    variant = str(contract.get("variant", ""))
    if not variant:
        raise ValueError("real run contract lacks its variant")
    authority, authority_summary = load_observation_sector_authority(
        observation_sector_authority_path,
        expected_sha256=observation_sector_authority_sha256,
    )
    manifest_rows = _read_bound_manifest_rows(
        release_root, manifest_sha256=manifest_sha256
    )
    authority_summary.update(
        _validate_authority_against_manifest(manifest_rows, authority)
    )
    rows = _read_development_rows(
        release_root,
        manifest_sha256=manifest_sha256,
        max_components=max_components,
    )
    train_representatives = _read_pca_training_representatives(
        release_root,
        manifest_sha256=manifest_sha256,
        max_components=max_pca_train_components,
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
    component_ids = [row["leakage_component_id"] for row in rows]
    visit_sectors = _authority_sectors_for_rows(
        rows, authority, partition="poc_development"
    )
    _authority_sectors_for_rows(
        train_representatives, authority, partition="poc_train"
    )
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
    baseline_training = [
        _model_window(
            release_root=release_root,
            row=row,
            variant=variant,
            mask_label="baseline-train-unmasked",
            unmasked=True,
        )
        for row in train_representatives
    ]
    huber_delta = float(contract["optimization"]["huber_delta"])
    left_representations, reconstruction = _encode(
        model=model,
        samples=paired_left,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=True,
        huber_delta=huber_delta,
    )
    right_representations, _ = _encode(
        model=model,
        samples=paired_right,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    visit_representations, _ = _encode(
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
    random_left_representations, _ = _encode(
        model=random_model,
        samples=paired_left,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    random_right_representations, _ = _encode(
        model=random_model,
        samples=paired_right,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    random_visit_representations, _ = _encode(
        model=random_model,
        samples=unmasked_visits,
        device=target_device,
        batch_size=batch_size,
        with_reconstruction=False,
        huber_delta=huber_delta,
    )
    selected_component_sha256 = _identity_sha256(tuple(by_component))
    selected_observation_sha256 = _identity_sha256(
        tuple(row["observation_key"] for row in rows)
    )
    train_component_sha256 = _identity_sha256(
        tuple(row["leakage_component_id"] for row in train_representatives)
    )
    train_observation_sha256 = _identity_sha256(
        tuple(row["observation_key"] for row in train_representatives)
    )

    pca_train_features = _feature_matrix(
        baseline_training, binned_flux_baseline_features
    )
    pca_fit, pca_fit_metadata = fit_train_pca_baseline(
        pca_train_features, n_components=pca_components
    )
    pca_fit_metadata.update(
        {
            "selected_training_components_sha256": train_component_sha256,
            "selected_training_observation_keys_sha256": train_observation_sha256,
            "feature_definition": (
                f"per-view visible-only means and support fractions in "
                f"{BASELINE_TIME_BINS} nonoverlapping time bins"
            ),
        }
    )
    pca_left = transform_train_pca_baseline(
        _feature_matrix(paired_left, binned_flux_baseline_features), pca_fit
    )
    pca_right = transform_train_pca_baseline(
        _feature_matrix(paired_right, binned_flux_baseline_features), pca_fit
    )
    pca_visits = transform_train_pca_baseline(
        _feature_matrix(unmasked_visits, binned_flux_baseline_features), pca_fit
    )

    scalar_train_features = _feature_matrix(
        baseline_training, robust_scalar_baseline_features
    )
    scalar_fit, scalar_fit_metadata = fit_train_scalar_standardizer(
        scalar_train_features
    )
    scalar_fit_metadata.update(
        {
            "selected_training_components_sha256": train_component_sha256,
            "selected_training_observation_keys_sha256": train_observation_sha256,
            "feature_definition": (
                "per-view visible-only median, MAD, q05, q25, q75, q95, "
                "adjacent-difference MAD, and support fraction"
            ),
        }
    )
    scalar_left = transform_train_scalar_baseline(
        _feature_matrix(paired_left, robust_scalar_baseline_features), scalar_fit
    )
    scalar_right = transform_train_scalar_baseline(
        _feature_matrix(paired_right, robust_scalar_baseline_features), scalar_fit
    )
    scalar_visits = transform_train_scalar_baseline(
        _feature_matrix(unmasked_visits, robust_scalar_baseline_features), scalar_fit
    )

    trained_representations = _representation_stage_summary(
        left_representations,
        right_representations,
        visit_representations,
        representative_components=representative_components,
        visit_components=component_ids,
        visit_sectors=visit_sectors,
    )
    random_representations = _representation_stage_summary(
        random_left_representations,
        random_right_representations,
        random_visit_representations,
        representative_components=representative_components,
        visit_components=component_ids,
        visit_sectors=visit_sectors,
    )
    trained_minus_random = {
        name: source_paired_trained_minus_random_summary(
            left_representations[name],
            right_representations[name],
            random_left_representations[name],
            random_right_representations[name],
            cluster_ids=representative_components,
        )
        for name in ("h_window", "z_window")
    }
    return {
        "schema_version": REPRESENTATION_HEALTH_SCHEMA_VERSION,
        "passed": True,
        "passed_definition": (
            "all evaluator-v2 provenance, sealed-access, numerical, and output "
            "checks completed; scientific go/no-go thresholds are applied by a "
            "separate frozen gate"
        ),
        "claim_limit": (
            "development-only representation-health evidence; no morphology, "
            "injection-retention, sealed-test, or foundation-model claim"
        ),
        "run": {
            "run_dir": str(root),
            "variant": variant,
            "architecture": contract.get("architecture"),
            "run_git_sha": contract.get("expected_git_sha"),
            "global_step": int(run_validation["global_step"]),
            "checkpoint_path": run_validation["selected_checkpoint_path"],
            "checkpoint_sha256": run_validation["selected_checkpoint_sha256"],
            "checkpoint_is_immutable_milestone": run_validation[
                "selected_checkpoint_is_immutable_milestone"
            ],
        },
        "evaluation_population": {
            "source_partition": "poc_development",
            "selected_leakage_components": int(len(by_component)),
            "selected_observation_visits": int(len(rows)),
            "selected_leakage_components_sha256": selected_component_sha256,
            "selected_observation_keys_sha256": selected_observation_sha256,
            "selection": "stable SHA-256 ordering of leakage components",
            "one_deterministic_window_per_selected_visit": True,
            "quality_mask_fallback": (
                "retain the selected visit and deterministically choose an eligible "
                "window only when the ordinary window has no required-view cadences"
            ),
            "model_visible_identifiers": False,
            "observation_sector_authority": authority_summary,
            "query_sector_excluded_retrieval": True,
        },
        "baseline_fit_population": {
            "source_partition": "poc_train",
            "selected_leakage_components": int(len(train_representatives)),
            "selected_observation_visits": int(len(train_representatives)),
            "selected_leakage_components_sha256": train_component_sha256,
            "selected_observation_keys_sha256": train_observation_sha256,
            "one_representative_visit_per_component": True,
            "development_rows": 0,
            "sealed_test_rows": 0,
        },
        "data_access": {
            "poc_train_window_reads": int(len(train_representatives)),
            "poc_development_window_reads": int(
                2 * len(representatives) + len(rows)
            ),
            "poc_train_unique_shards": int(
                len({row["relative_path"] for row in train_representatives})
            ),
            "poc_development_unique_shards": int(
                len({row["relative_path"] for row in rows})
            ),
            "poc_sealed_test_window_reads": 0,
            "poc_sealed_test_unique_shards": 0,
            "sealed_test_access_count": 0,
        },
        "masked_reconstruction": reconstruction,
        "trained_encoder": {
            "projection_spectrum": _model_projection_spectrum(model),
            "representations": trained_representations,
        },
        "random_encoder_control": {
            "seed": int(random_control_seed),
            "projection_spectrum": _model_projection_spectrum(random_model),
            "representations": random_representations,
        },
        "source_paired_trained_minus_random": trained_minus_random,
        "baselines": {
            "model_visible_identifiers_or_sector": False,
            "pca": {
                "fit": pca_fit_metadata,
                "development_evaluation": _baseline_evaluation_summary(
                    pca_left,
                    pca_right,
                    pca_visits,
                    representative_components=representative_components,
                    visit_components=component_ids,
                    visit_sectors=visit_sectors,
                ),
            },
            "robust_scalar": {
                "fit": scalar_fit_metadata,
                "development_evaluation": _baseline_evaluation_summary(
                    scalar_left,
                    scalar_right,
                    scalar_visits,
                    representative_components=representative_components,
                    visit_components=component_ids,
                    visit_sectors=visit_sectors,
                ),
            },
        },
        "not_evaluated": {
            "frozen_linear_morphology_probe": "requires separately frozen label registry",
            "injection_event_retention": "requires separately frozen development injection canary",
            "shortcut_probes_beyond_sector_exclusion": (
                "requires separately frozen audit metadata registry"
            ),
            "sealed_test": "not opened by this development-only evaluator",
        },
    }
