"""Cadence-local collapse diagnostics for TWIRL-FM0.3.

The diagnostic consumes one token per input cadence.  It never bins, merges,
resamples, or averages token representations along the cadence axis.  Adjacent
differences are formed only across valid, contiguous within-segment edges, and
are centered across the batch separately at each fixed edge position.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

G_DELTA_RETENTION_MINIMUM = 0.5
CADENCE_ENERGY_RETENTION_MINIMUM = 0.5
DIFFERENCE_EFFECTIVE_RANK_MINIMUM = 26.0
MAXIMUM_CONTIGUOUS_DELTA_CADENCES = 1.5
_ENERGY_TOLERANCE = 1.0e-12


def _effective_rank(centered_differences: np.ndarray, *, degrees_freedom: int) -> float:
    if degrees_freedom <= 0 or centered_differences.size == 0:
        return 0.0
    covariance = centered_differences.T @ centered_differences / float(degrees_freedom)
    eigenvalues = np.clip(np.linalg.eigvalsh(covariance), 0.0, None)
    if eigenvalues.size == 0:
        return 0.0
    maximum = float(eigenvalues[-1])
    numerical_tolerance = max(covariance.shape) * np.finfo(np.float64).eps * maximum
    eigenvalues[eigenvalues <= numerical_tolerance] = 0.0
    total = float(np.sum(eigenvalues))
    if total <= 0.0:
        return 0.0
    weights = eigenvalues[eigenvalues > 0.0] / total
    return float(np.exp(-np.sum(weights * np.log(weights))))


def summarize_cadence_difference_geometry(
    cadence_tokens: np.ndarray,
    valid: np.ndarray,
    segment_boundary: np.ndarray,
    delta_time_cadences: np.ndarray,
    *,
    maximum_contiguous_delta_cadences: float = MAXIMUM_CONTIGUOUS_DELTA_CADENCES,
) -> dict[str, Any]:
    """Measure data-dependent adjacent-token geometry without temporal pooling.

    For every eligible edge ``t -> t+1``, let
    ``delta_z[b,t] = z[b,t+1] - z[b,t]``.  Differences and both endpoint-token
    arrays are centered across the batch independently for each fixed edge.
    Positions with fewer than two eligible batch rows do not contribute.

    ``G_delta`` is the centered difference energy divided by the mean centered
    energy of the two endpoint-token arrays.  Effective rank is the entropy
    rank of the pooled covariance of the already position-centered difference
    rows.  Pooling those scalar/covariance sufficient statistics does not
    create or replace any cadence representation.
    """

    tokens = np.asarray(cadence_tokens, dtype=np.float64)
    valid_array = np.asarray(valid)
    boundary_array = np.asarray(segment_boundary)
    delta_time = np.asarray(delta_time_cadences, dtype=np.float64)
    if tokens.ndim != 3 or min(tokens.shape) < 1:
        raise ValueError("cadence_tokens must have shape [batch, cadence, feature]")
    batch_size, cadence_length, embedding_dim = tokens.shape
    if batch_size < 2 or cadence_length < 2:
        raise ValueError(
            "cadence collapse diagnostics require at least two batch rows and cadences"
        )
    expected_shape = (batch_size, cadence_length)
    if valid_array.shape != expected_shape:
        raise ValueError("valid must have shape [batch, cadence]")
    if boundary_array.shape != expected_shape:
        raise ValueError("segment_boundary must have shape [batch, cadence]")
    if delta_time.shape != expected_shape:
        raise ValueError("delta_time_cadences must have shape [batch, cadence]")
    if valid_array.dtype != np.bool_ or boundary_array.dtype != np.bool_:
        raise ValueError("valid and segment_boundary must be boolean")
    maximum_delta = float(maximum_contiguous_delta_cadences)
    if not math.isfinite(maximum_delta) or maximum_delta <= 0.0:
        raise ValueError(
            "maximum_contiguous_delta_cadences must be finite and positive"
        )

    right_delta = delta_time[:, 1:]
    edge_eligible = (
        valid_array[:, :-1]
        & valid_array[:, 1:]
        & ~boundary_array[:, 1:]
        & np.isfinite(right_delta)
        & (right_delta > 0.0)
        & (right_delta <= maximum_delta)
    )
    finite_endpoints = np.all(np.isfinite(tokens[:, :-1, :]), axis=-1) & np.all(
        np.isfinite(tokens[:, 1:, :]), axis=-1
    )
    if np.any(edge_eligible & ~finite_endpoints):
        raise ValueError("eligible cadence edge contains a non-finite endpoint token")

    edge_counts = np.sum(edge_eligible, axis=0, dtype=np.int64)
    statistical_positions = edge_counts >= 2
    statistical_mask = edge_eligible & statistical_positions[None, :]
    degrees_freedom = int(np.sum(np.maximum(edge_counts - 1, 0)))

    difference_energy = 0.0
    endpoint_energy = 0.0
    centered_difference_rows: list[np.ndarray] = []
    for edge_index in np.flatnonzero(statistical_positions):
        rows = edge_eligible[:, edge_index]
        left = tokens[rows, edge_index, :]
        right = tokens[rows, edge_index + 1, :]
        differences = right - left
        centered_differences = differences - np.mean(differences, axis=0, keepdims=True)
        centered_left = left - np.mean(left, axis=0, keepdims=True)
        centered_right = right - np.mean(right, axis=0, keepdims=True)
        difference_energy += float(np.sum(np.square(centered_differences)))
        endpoint_energy += 0.5 * float(
            np.sum(np.square(centered_left)) + np.sum(np.square(centered_right))
        )
        centered_difference_rows.append(centered_differences)

    endpoint_tolerance = _ENERGY_TOLERANCE * max(1.0, endpoint_energy)
    difference_tolerance = _ENERGY_TOLERANCE * max(1.0, difference_energy)
    endpoint_energy_positive = endpoint_energy > endpoint_tolerance
    difference_energy_positive = difference_energy > difference_tolerance
    g_delta = (
        float(difference_energy / endpoint_energy)
        if endpoint_energy_positive and difference_energy_positive
        else 0.0
    )
    difference_matrix = (
        np.concatenate(centered_difference_rows, axis=0)
        if centered_difference_rows
        else np.empty((0, embedding_dim), dtype=np.float64)
    )
    effective_rank = (
        _effective_rank(difference_matrix, degrees_freedom=degrees_freedom)
        if difference_energy_positive
        else 0.0
    )
    diagnostic_valid = bool(
        np.any(statistical_positions)
        and degrees_freedom > 0
        and endpoint_energy_positive
        and difference_energy_positive
        and math.isfinite(g_delta)
        and math.isfinite(effective_rank)
    )
    return {
        "batch_size": int(batch_size),
        "cadence_length": int(cadence_length),
        "embedding_dim": int(embedding_dim),
        "maximum_contiguous_delta_cadences": maximum_delta,
        "candidate_adjacent_edges": int(batch_size * (cadence_length - 1)),
        "eligible_contiguous_edges": int(np.count_nonzero(edge_eligible)),
        "statistical_edges": int(np.count_nonzero(statistical_mask)),
        "statistical_edge_positions": int(np.count_nonzero(statistical_positions)),
        "statistical_degrees_of_freedom": degrees_freedom,
        "eligible_edges_by_position": [int(value) for value in edge_counts],
        "centered_difference_energy": difference_energy,
        "centered_endpoint_energy": endpoint_energy,
        "g_delta": g_delta,
        "difference_effective_rank": effective_rank,
        "endpoint_energy_positive": bool(endpoint_energy_positive),
        "difference_energy_positive": bool(difference_energy_positive),
        "diagnostic_valid": diagnostic_valid,
        "cadence_representations_averaged_or_merged": False,
    }


def _evaluate_cadence_collapse_math(
    step0_tokens: np.ndarray,
    step2000_tokens: np.ndarray,
    valid: np.ndarray,
    segment_boundary: np.ndarray,
    delta_time_cadences: np.ndarray,
) -> dict[str, Any]:
    """Compute the frozen collapse criteria without making an artifact claim.

    This private numerical kernel accepts arrays only and therefore cannot
    authenticate their scientific identity.  It deliberately emits no
    authoritative ``passed`` field.  The artifact runner is the only public
    path that may bind these criteria to real panel/checkpoint bytes and turn
    them into a gate decision.
    """

    step0_array = np.asarray(step0_tokens)
    step2000_array = np.asarray(step2000_tokens)
    if step0_array.ndim != 3 or step2000_array.ndim != 3:
        raise ValueError("step-0 and step-2000 tokens must have shape [B,L,E]")
    initial = summarize_cadence_difference_geometry(
        step0_array,
        valid,
        segment_boundary,
        delta_time_cadences,
        maximum_contiguous_delta_cadences=MAXIMUM_CONTIGUOUS_DELTA_CADENCES,
    )
    trained = summarize_cadence_difference_geometry(
        step2000_array,
        valid,
        segment_boundary,
        delta_time_cadences,
        maximum_contiguous_delta_cadences=MAXIMUM_CONTIGUOUS_DELTA_CADENCES,
    )
    if (
        initial["batch_size"],
        initial["cadence_length"],
        initial["embedding_dim"],
    ) != (
        trained["batch_size"],
        trained["cadence_length"],
        trained["embedding_dim"],
    ):
        raise ValueError("step-0 and step-2000 cadence-token shapes must match")

    g_delta_minimum = G_DELTA_RETENTION_MINIMUM * float(initial["g_delta"])
    g_delta_passed = bool(
        initial["diagnostic_valid"]
        and trained["diagnostic_valid"]
        and float(trained["g_delta"]) >= g_delta_minimum
    )
    cadence_energy_retention = (
        float(trained["centered_endpoint_energy"])
        / float(initial["centered_endpoint_energy"])
        if initial["diagnostic_valid"]
        and float(initial["centered_endpoint_energy"]) > 0.0
        else 0.0
    )
    cadence_energy_passed = bool(
        initial["diagnostic_valid"]
        and trained["diagnostic_valid"]
        and math.isfinite(cadence_energy_retention)
        and cadence_energy_retention >= CADENCE_ENERGY_RETENTION_MINIMUM
    )
    rank_passed = bool(
        trained["diagnostic_valid"]
        and float(trained["difference_effective_rank"])
        >= DIFFERENCE_EFFECTIVE_RANK_MINIMUM
    )
    criteria_satisfied = bool(g_delta_passed and cadence_energy_passed and rank_passed)
    return {
        "schema_version": "twirl_fm0_cadence_collapse_math_v1",
        "authoritative_artifact_gate": False,
        "step0": initial,
        "step2000": trained,
        "thresholds": {
            "g_delta_step2000_over_step0_minimum": G_DELTA_RETENTION_MINIMUM,
            "cadence_energy_step2000_over_step0_minimum": (
                CADENCE_ENERGY_RETENTION_MINIMUM
            ),
            "difference_effective_rank_minimum": (DIFFERENCE_EFFECTIVE_RANK_MINIMUM),
            "maximum_contiguous_delta_cadences": (MAXIMUM_CONTIGUOUS_DELTA_CADENCES),
        },
        "g_delta_step2000_minimum": g_delta_minimum,
        "g_delta_passed": g_delta_passed,
        "cadence_energy_step2000_over_step0": cadence_energy_retention,
        "cadence_energy_passed": cadence_energy_passed,
        "difference_effective_rank_passed": rank_passed,
        "criteria_satisfied": criteria_satisfied,
        "cadence_representations_averaged_or_merged": False,
    }


__all__ = [
    "CADENCE_ENERGY_RETENTION_MINIMUM",
    "DIFFERENCE_EFFECTIVE_RANK_MINIMUM",
    "G_DELTA_RETENTION_MINIMUM",
    "MAXIMUM_CONTIGUOUS_DELTA_CADENCES",
    "summarize_cadence_difference_geometry",
]
