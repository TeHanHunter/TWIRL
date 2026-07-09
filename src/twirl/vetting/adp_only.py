"""Hard ADP-only data contract for S56 vetting and machine learning."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


ADP_ONLY_CONTRACT_VERSION = "s56_adp_pair_v1"
ADP_ONLY_APERTURES: tuple[str, str] = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
ADP_ONLY_APERTURE_SIGNATURE = ",".join(ADP_ONLY_APERTURES)
ADP_ONLY_INJECTION_H5_BY_ORIGIN: dict[str, str] = {
    "s56_recovery50_teacher_queue": (
        "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
        "recovery50_adp_pair_subset/injected_lightcurves.h5"
    ),
    "s56_recovery50_teacher_queue_next4k": (
        "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
        "recovery50_adp_pair_subset_next4k/injected_lightcurves.h5"
    ),
}

# Only explicitly ADP-derived search/vetting quantities are allowed into the
# scalar branch. Generic queue BLS values are intentionally absent because
# historical real queues populated them from canonical DET_FLUX products.
ADP_ONLY_METADATA_COLUMNS: tuple[str, ...] = (
    "tmag",
    "anchor_period_d",
    "anchor_duration_min",
    "anchor_sde",
    "adp_sml_peak_rank",
    "adp_sml_period_d",
    "adp_sml_duration_min",
    "adp_sml_depth",
    "adp_sml_depth_snr",
    "adp_sml_sde",
    "adp_sml_log_power",
    "adp_sml_own_even_n_in",
    "adp_sml_own_even_depth",
    "adp_sml_own_odd_n_in",
    "adp_sml_own_odd_depth",
    "adp_sml_own_even_odd_depth_delta",
    "adp_sml_own_even_odd_sigma_delta",
    "adp_sml_anchor_depth",
    "adp_sml_anchor_snr",
    "adp_sml_anchor_n_in",
    "adp_sml_anchor_even_n_in",
    "adp_sml_anchor_even_depth",
    "adp_sml_anchor_odd_n_in",
    "adp_sml_anchor_odd_depth",
    "adp_sml_anchor_even_odd_depth_delta",
    "adp_sml_anchor_even_odd_sigma_delta",
    "adp_sml_trend_ptp",
    "adp_peak_rank",
    "adp_period_d",
    "adp_duration_min",
    "adp_depth",
    "adp_depth_snr",
    "adp_sde",
    "adp_log_power",
    "adp_own_even_n_in",
    "adp_own_even_depth",
    "adp_own_odd_n_in",
    "adp_own_odd_depth",
    "adp_own_even_odd_depth_delta",
    "adp_own_even_odd_sigma_delta",
    "adp_anchor_depth",
    "adp_anchor_snr",
    "adp_anchor_n_in",
    "adp_anchor_even_n_in",
    "adp_anchor_even_depth",
    "adp_anchor_odd_n_in",
    "adp_anchor_odd_depth",
    "adp_anchor_even_odd_depth_delta",
    "adp_anchor_even_odd_sigma_delta",
    "adp_trend_ptp",
    "aperture_period_rel_delta",
    "aperture_depth_ratio_primary_over_small",
    "aperture_disagreement_flag",
)


_LEGACY_GENERIC_RESULT_COLUMNS = frozenset(
    {
        "apertures_agree",
        "baseline_d",
        "blind_rank",
        "bls_run_id",
        "both_apertures_top1_recovered",
        "cadence_alias_delta",
        "cadence_alias_flag",
        "cadence_period_n",
        "cadence_period_nearest_n",
        "candidate_key",
        "centroid_delta_pix",
        "centroid_dx_pix",
        "centroid_dy_pix",
        "centroid_pass",
        "centroid_sigma_oot_pix",
        "centroid_status",
        "centroid_z",
        "class_rank",
        "display_ephemeris_period_rel_delta",
        "dur_envelope_min",
        "dur_envelope_pass",
        "dur_grid_ceiling_hit",
        "human_candidate_key",
        "min_twoap_sde",
        "n_apertures_agree",
        "n_cad_kept_max",
        "n_harmonic_apertures_agree",
        "n_in_transit",
        "n_oot_band",
        "n_topn_apertures_agree",
        "p_alias_pass",
        "p_cluster_pass",
        "period_cluster_count",
        "period_rel_err",
        "recovery_status",
        "rep_aperture",
        "rep_peak_rank",
        "roche_pass",
        "roche_period_d",
        "strict_top1_recovered",
        "t0_phase_err_min",
        "topn_apertures_agree",
        "topn_recovery_status",
        "truth_n_good_in_transit",
        "truth_n_in_transit",
        "vet_class",
        "vetted_rank",
        "any_exact_or_harmonic_recovered",
    }
)


def validate_adp_only_apertures(apertures: Sequence[str]) -> tuple[str, str]:
    """Require the ordered small-primary ADP pair used by active S56 models."""

    normalized = tuple(str(value) for value in apertures)
    if normalized != ADP_ONLY_APERTURES:
        raise ValueError(
            f"active S56 ML requires apertures={ADP_ONLY_APERTURES}; got {normalized}"
        )
    return ADP_ONLY_APERTURES


def canonical_det_flux_columns(columns: Iterable[str]) -> list[str]:
    """Return columns that explicitly encode a non-ADP DET_FLUX product."""

    bad: list[str] = []
    for column in columns:
        text = str(column)
        upper = text.upper()
        if "DET_FLUX" in upper and "DET_FLUX_ADP" not in upper:
            bad.append(text)
        elif text.startswith("sml_"):
            # Historical next4k injection metrics used DET_FLUX_SML and were
            # normalized to the ambiguous ``sml_*`` prefix.
            bad.append(text)
    return sorted(set(bad))


def classify_period_relation(
    reference: pd.Series,
    comparison: pd.Series,
    *,
    tolerance: float = 0.01,
    harmonic_factors: tuple[float, ...] = (0.25, 1.0 / 3.0, 0.5, 2.0, 3.0, 4.0),
) -> pd.Series:
    """Classify two period columns as exact, harmonic, unrelated, or missing."""

    left = pd.to_numeric(reference, errors="coerce").to_numpy(dtype=float)
    right = pd.to_numeric(comparison, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(left), "missing", dtype=object)
    finite = np.isfinite(left) & np.isfinite(right) & (left > 0) & (right > 0)
    ratio = np.full(len(left), np.nan, dtype=float)
    ratio[finite] = left[finite] / right[finite]
    exact = finite & (np.abs(ratio - 1.0) <= float(tolerance))
    result[exact] = "exact"
    remaining = finite & ~exact
    if np.any(remaining):
        factors = np.asarray(harmonic_factors, dtype=float)
        errors = np.min(np.abs(ratio[remaining, None] / factors[None, :] - 1.0), axis=1)
        harmonic_idx = np.flatnonzero(remaining)[errors <= float(tolerance)]
        result[harmonic_idx] = "harmonic"
    result[finite & (result == "missing")] = "unrelated"
    return pd.Series(result, index=reference.index, dtype=object)


def _as_bool(series: pd.Series, index: pd.Index) -> pd.Series:
    if series is None:
        return pd.Series(False, index=index, dtype=bool)
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def _combined_metrics(paths: Sequence[Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    wanted_prefixes = ("anchor_", "adp_sml_", "adp_", "aperture_")
    for path in paths:
        path = Path(path)
        if not path.exists():
            continue
        frame = pd.read_csv(path, low_memory=False)
        if "review_id" not in frame:
            continue
        keep = [
            column
            for column in frame.columns
            if column == "review_id" or column.startswith(wanted_prefixes)
        ]
        frames.append(frame.loc[:, keep].copy())
    if not frames:
        return pd.DataFrame(columns=["review_id"])
    metrics = pd.concat(frames, ignore_index=True, sort=False)
    metrics["review_id"] = metrics["review_id"].fillna("").astype(str)
    return metrics.drop_duplicates("review_id", keep="last")


def build_adp_only_training_frame(
    frame: pd.DataFrame,
    *,
    metrics_tables: Sequence[Path],
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Return a sanitized active training table with ADP-only observables."""

    out = frame.copy()
    out["review_id"] = out["review_id"].fillna("").astype(str)
    metrics = _combined_metrics(metrics_tables)
    metric_columns = [column for column in metrics.columns if column != "review_id"]
    overlap = [column for column in metric_columns if column in out.columns]
    if overlap:
        out = out.drop(columns=overlap)
    out = out.merge(metrics, on="review_id", how="left", validate="one_to_one")

    anchor_period = pd.to_numeric(out.get("anchor_period_d"), errors="coerce")
    anchor_t0 = pd.to_numeric(out.get("anchor_t0_bjd"), errors="coerce")
    anchor_duration = pd.to_numeric(out.get("anchor_duration_min"), errors="coerce")
    small_period = pd.to_numeric(out.get("adp_sml_period_d"), errors="coerce")
    primary_period = pd.to_numeric(out.get("adp_period_d"), errors="coerce")
    display_period = pd.to_numeric(out.get("display_period_d"), errors="coerce")
    anchor_aperture = out.get("anchor_aperture", pd.Series("", index=out.index)).fillna("").astype(str)

    small_ok = (
        anchor_aperture.eq(ADP_ONLY_APERTURES[0])
        & anchor_period.gt(0)
        & anchor_t0.notna()
        & anchor_duration.gt(0)
        & small_period.gt(0)
    )
    display_matches_small = classify_period_relation(display_period, small_period).eq("exact")
    primary_ok = primary_period.gt(0)
    review_ok = small_ok & display_matches_small & primary_ok
    reasons = pd.Series("", index=out.index, dtype=object)
    reasons.loc[~small_ok] = "missing_adp_small_anchor"
    reasons.loc[small_ok & ~display_matches_small] = "display_not_adp_small_anchor"
    reasons.loc[small_ok & display_matches_small & ~primary_ok] = "missing_adp_primary_supplement"

    out["adp_only_contract_version"] = ADP_ONLY_CONTRACT_VERSION
    out["adp_only_review_ok"] = review_ok.astype(bool)
    out["adp_only_exclusion_reason"] = reasons

    if {"source_kind", "origin_queue", "source_h5"}.issubset(out.columns):
        injected = out["source_kind"].fillna("").astype(str).eq("injection_recovery")
        rebuilt_h5 = out["origin_queue"].fillna("").astype(str).map(ADP_ONLY_INJECTION_H5_BY_ORIGIN)
        out.loc[injected & rebuilt_h5.notna(), "source_h5"] = rebuilt_h5.loc[
            injected & rebuilt_h5.notna()
        ]

    # Generic active ephemeris/BLS fields are rewritten from ADP-small metrics.
    replacements = {
        "period_d": "anchor_period_d",
        "t0_bjd": "anchor_t0_bjd",
        "duration_min": "anchor_duration_min",
        "depth": "adp_sml_depth",
        "depth_snr": "adp_sml_depth_snr",
        "sde_max": "adp_sml_sde",
        "rep_peak_rank": "adp_sml_peak_rank",
    }
    for destination, source in replacements.items():
        out[destination] = pd.to_numeric(out.get(source), errors="coerce")
    out["rep_aperture"] = ADP_ONLY_APERTURES[0]
    period_delta = pd.to_numeric(out.get("aperture_period_rel_delta"), errors="coerce")
    agrees = primary_ok & period_delta.le(0.02)
    out["n_apertures_agree"] = np.where(primary_ok, np.where(agrees, 2, 1), 0)
    out["apertures_agree"] = np.where(
        agrees,
        ",".join(ADP_ONLY_APERTURES),
        np.where(small_ok, ADP_ONLY_APERTURES[0], ""),
    )
    out["min_twoap_sde"] = np.fmin(
        pd.to_numeric(out.get("adp_sml_sde"), errors="coerce"),
        pd.to_numeric(out.get("adp_sde"), errors="coerce"),
    )

    if "main_teacher_include" in out:
        out["main_teacher_include"] = _as_bool(out["main_teacher_include"], out.index) & review_ok
    if "compact_morphology_include" in out:
        out["compact_morphology_include"] = _as_bool(out["compact_morphology_include"], out.index) & review_ok
    if "training_split" in out:
        out.loc[~review_ok, "training_split"] = "unlabeled_or_audit"

    invalid_small = ~small_ok
    for column in ("display_period_d", "display_t0_bjd", "display_duration_min"):
        if column in out:
            out.loc[invalid_small, column] = np.nan
    if "display_ephemeris_source" in out:
        out.loc[invalid_small, "display_ephemeris_source"] = "invalid_adp"

    explicit_bad = canonical_det_flux_columns(out.columns)
    drop_columns = set(explicit_bad)
    drop_columns.update(column for column in _LEGACY_GENERIC_RESULT_COLUMNS if column in out.columns)
    drop_columns.update(column for column in out.columns if column.startswith(("queue_", "selection_")))
    out = out.drop(columns=sorted(drop_columns), errors="ignore")

    remaining_bad = canonical_det_flux_columns(out.columns)
    if remaining_bad:
        raise ValueError(f"canonical DET_FLUX columns survived ADP sanitization: {remaining_bad}")
    summary = {
        "contract_version": ADP_ONLY_CONTRACT_VERSION,
        "n_rows": int(len(out)),
        "n_adp_only_review_ok": int(review_ok.sum()),
        "n_excluded": int((~review_ok).sum()),
        "exclusion_reason_counts": {
            str(key): int(value)
            for key, value in reasons[~review_ok].value_counts().sort_index().items()
        },
        "dropped_columns": sorted(drop_columns),
        "active_apertures": list(ADP_ONLY_APERTURES),
    }
    return out, summary


def assert_adp_only_training_frame(frame: pd.DataFrame) -> None:
    """Reject active training tables that do not satisfy the ADP contract."""

    if "adp_only_contract_version" not in frame:
        raise ValueError("training table is missing adp_only_contract_version")
    versions = set(frame["adp_only_contract_version"].dropna().astype(str))
    if versions != {ADP_ONLY_CONTRACT_VERSION}:
        raise ValueError(f"unexpected ADP-only contract versions: {sorted(versions)}")
    bad = canonical_det_flux_columns(frame.columns)
    if bad:
        raise ValueError(f"training table contains canonical DET_FLUX columns: {bad}")


def assert_adp_only_tensor_rows(frame: pd.DataFrame) -> None:
    """Reject tensor indices built before the ordered two-ADP contract."""

    required = {"adp_only_contract_version", "tensor_apertures"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"tensor rows are missing ADP contract columns: {missing}")
    versions = set(frame["adp_only_contract_version"].dropna().astype(str))
    if versions != {ADP_ONLY_CONTRACT_VERSION}:
        raise ValueError(f"unexpected tensor ADP contract versions: {sorted(versions)}")
    signatures = set(frame["tensor_apertures"].dropna().astype(str))
    if signatures != {ADP_ONLY_APERTURE_SIGNATURE}:
        raise ValueError(f"unexpected tensor aperture signatures: {sorted(signatures)}")


def assert_adp_only_search_frame(frame: pd.DataFrame) -> None:
    """Reject candidate/BLS tables that are not from the active ADP search."""

    required = {
        "adp_only_contract_version",
        "aperture",
        "bls_search_branch",
        "period_d",
        "t0_bjd",
        "duration_min",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"ADP search table is missing required columns: {missing}")
    versions = set(frame["adp_only_contract_version"].dropna().astype(str))
    if versions != {ADP_ONLY_CONTRACT_VERSION}:
        raise ValueError(f"unexpected ADP-only search contract versions: {sorted(versions)}")
    branches = set(frame["bls_search_branch"].dropna().astype(str))
    if branches != {"current_adp"}:
        raise ValueError(f"unexpected BLS search branches: {sorted(branches)}")
    apertures = set(frame["aperture"].dropna().astype(str))
    unexpected = sorted(apertures - set(ADP_ONLY_APERTURES))
    if unexpected:
        raise ValueError(f"candidate/BLS table contains non-ADP apertures: {unexpected}")
    if not apertures:
        raise ValueError("candidate/BLS table contains no aperture values")
    bad = canonical_det_flux_columns(frame.columns)
    if bad:
        raise ValueError(f"candidate/BLS table contains canonical DET_FLUX columns: {bad}")


__all__ = [
    "ADP_ONLY_APERTURES",
    "ADP_ONLY_APERTURE_SIGNATURE",
    "ADP_ONLY_CONTRACT_VERSION",
    "ADP_ONLY_INJECTION_H5_BY_ORIGIN",
    "ADP_ONLY_METADATA_COLUMNS",
    "assert_adp_only_search_frame",
    "assert_adp_only_tensor_rows",
    "assert_adp_only_training_frame",
    "build_adp_only_training_frame",
    "canonical_det_flux_columns",
    "classify_period_relation",
    "validate_adp_only_apertures",
]
