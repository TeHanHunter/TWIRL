"""Locked BLS definition for the S56 A2v1 enrichment candidate table."""
from __future__ import annotations

import hashlib
import json
from dataclasses import fields
from typing import Any, Mapping

from twirl.search.bls import BLSConfig


A2V1_TEACHER_BLS_SEARCH_CONTRACT = "s56_a2v1_teacher_bls_search_v1"
A2V1_TEACHER_MIN_CADENCES = 200


def approved_a2v1_teacher_bls_config() -> dict[str, Any]:
    """Return a fresh JSON-native copy of the approved production search."""

    return {
        "apertures": ["DET_FLUX_ADP_SML", "DET_FLUX_ADP"],
        "n_periods": 50_000,
        "n_peaks": 10,
        "p_min_d": 0.12,
        "p_max_cap_d": 15.0,
        "max_period_fraction": 0.45,
        "durations_min": [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 30.0],
        "period_mask_frac": 0.005,
        "period_bin_edges": [],
        "max_peaks_per_period_bin": 0,
        "min_cadences": A2V1_TEACHER_MIN_CADENCES,
        "sigma_clip": 5.0,
        "orbit_edge_trim_d": 0.0,
        "source_product_tag": "A2v1",
    }


def approved_a2v1_teacher_bls_runtime_config() -> BLSConfig:
    """Return the runtime config represented by the locked JSON contract."""

    payload = approved_a2v1_teacher_bls_config()
    metadata_fields = {"source_product_tag"}
    runtime_fields = {field.name for field in fields(BLSConfig)}
    observed_runtime_fields = set(payload) - metadata_fields
    if observed_runtime_fields != runtime_fields:
        missing = sorted(runtime_fields - observed_runtime_fields)
        unexpected = sorted(observed_runtime_fields - runtime_fields)
        raise RuntimeError(
            "locked A2v1 BLS JSON/runtime fields disagree: "
            f"missing={missing}, unexpected={unexpected}"
        )
    runtime_payload = {
        key: value for key, value in payload.items() if key in runtime_fields
    }
    for key in ("apertures", "durations_min", "period_bin_edges"):
        runtime_payload[key] = tuple(runtime_payload[key])
    return BLSConfig(**runtime_payload)


def bls_config_sha256(config: Mapping[str, Any]) -> str:
    """Fingerprint a JSON-native BLS configuration without path dependence."""

    encoded = json.dumps(
        dict(config), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


__all__ = [
    "A2V1_TEACHER_BLS_SEARCH_CONTRACT",
    "A2V1_TEACHER_MIN_CADENCES",
    "approved_a2v1_teacher_bls_config",
    "approved_a2v1_teacher_bls_runtime_config",
    "bls_config_sha256",
]
