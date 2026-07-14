"""BLS-to-Teacher-v2 recovery accounting for locked and external injections."""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Mapping

import numpy as np
import pandas as pd


TMAG_SUPPORT_LABELS: tuple[str, ...] = (
    "Tmag < 17",
    "17 <= Tmag < 18",
    "18 <= Tmag < 19",
    "Tmag >= 19",
)


def _as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.fillna("").astype(str).str.lower().isin(
        {"1", "1.0", "true", "t", "yes", "y"}
    )


def aggregate_compact_recovery(
    manifest: pd.DataFrame,
    scored_candidates: pd.DataFrame,
    *,
    threshold: float,
    score_column: str = "p_compact_transit",
    max_peak_rank: int = 5,
    outcome_prefix: str = "teacher_v2",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return one BLS and compact-retention outcome per scheduled injection."""

    required = {
        "injection_id",
        "rep_peak_rank",
        "is_injected_signal_peak",
        score_column,
    }
    missing = sorted(required - set(scored_candidates.columns))
    if missing:
        raise KeyError(f"Teacher recovery scores are missing columns: {missing}")
    if manifest["injection_id"].fillna("").astype(str).duplicated().any():
        raise ValueError("recovery manifest injection IDs are not unique")
    scored = scored_candidates.copy()
    scored = scored.loc[
        pd.to_numeric(scored["rep_peak_rank"], errors="coerce").between(
            1, int(max_peak_rank)
        )
    ].copy()
    scored["_truth_match"] = _as_bool(scored["is_injected_signal_peak"])
    scored["_compact_pass"] = pd.to_numeric(
        scored[score_column], errors="coerce"
    ).ge(float(threshold))
    rows: list[dict[str, Any]] = []
    for injection_id, group in scored.groupby("injection_id", sort=False):
        matched = group.loc[group["_truth_match"]]
        rows.append(
            {
                "injection_id": str(injection_id),
                "n_scored_candidates": int(len(group)),
                "n_truth_matched_candidates": int(len(matched)),
                "bls_top5_recovered": bool(len(matched)),
                f"{outcome_prefix}_compact_recovered": bool(
                    len(matched) and matched["_compact_pass"].any()
                ),
                "matched_compact_score_max": (
                    float(pd.to_numeric(matched[score_column], errors="coerce").max())
                    if len(matched)
                    else np.nan
                ),
            }
        )
    outcome = pd.DataFrame(rows)
    result = manifest.merge(outcome, on="injection_id", how="left", validate="one_to_one")
    compact_column = f"{outcome_prefix}_compact_recovered"
    for column in ("bls_top5_recovered", compact_column):
        result[column] = result[column].fillna(False).astype(bool)
    result["n_scored_candidates"] = result["n_scored_candidates"].fillna(0).astype(int)
    result["n_truth_matched_candidates"] = result["n_truth_matched_candidates"].fillna(0).astype(int)
    bls = result["bls_top5_recovered"]
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "threshold": float(threshold),
        "score_column": score_column,
        "max_peak_rank": int(max_peak_rank),
        "outcome_prefix": outcome_prefix,
        "n_injections": int(len(result)),
        "n_bls_top5_recovered": int(bls.sum()),
        "bls_top5_fraction": float(bls.mean()) if len(result) else np.nan,
        "n_compact_recovered": int(result[compact_column].sum()),
        "compact_end_to_end_fraction": float(result[compact_column].mean()) if len(result) else np.nan,
        "compact_conditional_on_bls": (
            float(result.loc[bls, compact_column].mean()) if bls.any() else np.nan
        ),
    }
    return result, summary


def bls_topk_recovery_table(
    manifest: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    max_k: int = 5,
) -> pd.DataFrame:
    """Report exact top-1/3/5 recovery on one immutable candidate set."""

    truth = candidates.copy()
    truth["_truth_match"] = _as_bool(truth["is_injected_signal_peak"])
    rank = pd.to_numeric(truth["rep_peak_rank"], errors="coerce")
    rows: list[dict[str, Any]] = []
    for k in (1, 3, int(max_k)):
        recovered_ids = set(
            truth.loc[truth["_truth_match"] & rank.le(k), "injection_id"].astype(str)
        )
        recovered = manifest["injection_id"].astype(str).isin(recovered_ids)
        rows.append(
            {
                "top_k": int(k),
                "n_injections": int(len(manifest)),
                "n_recovered": int(recovered.sum()),
                "recovery_fraction": float(recovered.mean()) if len(manifest) else np.nan,
            }
        )
    return pd.DataFrame(rows).drop_duplicates("top_k").sort_values("top_k")


def period_radius_tmag_support(manifest: pd.DataFrame) -> pd.DataFrame:
    """Summarize host and grid support without smoothing or host reuse."""

    required = {"tic", "tmag", "grid_cell_id", "injection_id"}
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise KeyError(f"injection support manifest is missing columns: {missing}")
    work = manifest.copy()
    work["tic"] = pd.to_numeric(work["tic"], errors="raise").astype(np.int64)
    work["tmag"] = pd.to_numeric(work["tmag"], errors="coerce")
    if work["tmag"].isna().any():
        raise ValueError("injection support manifest contains invalid Tmag values")
    work["tmag_bin"] = pd.cut(
        work["tmag"],
        [-np.inf, 17.0, 18.0, 19.0, np.inf],
        right=False,
        labels=TMAG_SUPPORT_LABELS,
    )
    cells = pd.Index(sorted(work["grid_cell_id"].astype(str).unique()))
    rows: list[dict[str, Any]] = []
    for label in TMAG_SUPPORT_LABELS:
        subset = work.loc[work["tmag_bin"].astype(str).eq(label)].copy()
        support = (
            subset.assign(grid_cell_id=subset["grid_cell_id"].astype(str))
            .groupby("grid_cell_id")
            .size()
            .reindex(cells, fill_value=0)
            .to_numpy(dtype=int)
        )
        represented = int(np.count_nonzero(support))
        n_injections = int(len(subset))
        n_hosts = int(subset["tic"].nunique())
        coverage = float(represented / len(cells)) if len(cells) else 0.0
        rows.append(
            {
                "tmag_bin": label,
                "n_injections": n_injections,
                "n_unique_hosts": n_hosts,
                "host_reuse_factor": (
                    float(n_injections / n_hosts) if n_hosts else np.nan
                ),
                "n_expected_cells": int(len(cells)),
                "n_represented_cells": represented,
                "cell_coverage_fraction": coverage,
                "cell_support_min": int(support.min()) if len(support) else 0,
                "cell_support_p10": (
                    float(np.quantile(support, 0.10)) if len(support) else 0.0
                ),
                "cell_support_median": (
                    float(np.median(support)) if len(support) else 0.0
                ),
                "cell_support_max": int(support.max()) if len(support) else 0,
                "bright_enrichment_recommended": bool(
                    label != "Tmag >= 19" and coverage < 0.90
                ),
                "sampling_role": "primary_host_disjoint",
            }
        )
    return pd.DataFrame(rows)


def compare_recovery_models(
    outcomes: Mapping[str, pd.DataFrame],
    *,
    outcome_columns: Mapping[str, str],
) -> pd.DataFrame:
    """Compare model retention on the exact intersection of injection IDs."""

    names = list(outcomes)
    if set(names) != set(outcome_columns):
        raise ValueError("recovery outcomes and outcome columns must have identical keys")
    common: set[str] | None = None
    indexed: dict[str, pd.DataFrame] = {}
    for name in names:
        frame = outcomes[name].copy()
        ids = set(frame["injection_id"].astype(str))
        common = ids if common is None else common & ids
        indexed[name] = frame.set_index(frame["injection_id"].astype(str))
    common_ids = sorted(common or set())
    rows: list[dict[str, Any]] = []
    for name in names:
        frame = indexed[name].loc[common_ids]
        column = outcome_columns[name]
        bls = _as_bool(frame["bls_top5_recovered"])
        recovered = _as_bool(frame[column])
        rows.append(
            {
                "model": name,
                "n_common_injections": len(common_ids),
                "n_bls_recovered": int(bls.sum()),
                "n_model_recovered": int(recovered.sum()),
                "end_to_end_fraction": float(recovered.mean()) if len(frame) else np.nan,
                "conditional_on_bls": float(recovered[bls].mean()) if bls.any() else np.nan,
            }
        )
    return pd.DataFrame(rows)


__all__ = [
    "aggregate_compact_recovery",
    "bls_topk_recovery_table",
    "compare_recovery_models",
    "period_radius_tmag_support",
]
