"""CNN tensor and teacher helpers for the S56 recovery50 vetting queue."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from twirl.io.hlsp import BJDREFI, HLSPLightCurve
from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURE_SIGNATURE,
    ADP_ONLY_CONTRACT_VERSION,
    assert_adp_only_tensor_rows,
    assert_adp_only_training_frame,
    canonical_det_flux_columns,
    validate_adp_only_apertures,
)
from twirl.vetting.recovery50_teacher import (
    DEFAULT_APERTURES,
    json_default,
    leakage_columns,
    load_feature_table,
    metadata_feature_columns,
    read_table,
    select_bls_ephemeris,
    select_teacher_classes,
    select_model_ephemeris,
    write_table,
    _confusion_metrics,
    _read_lc_for_row,
    _safe_float,
)


DEFAULT_FOLDED_POINTS = 512
DEFAULT_CONTEXT_POINTS = 512
DEFAULT_EVENT_POINTS = 128
DEFAULT_TIMELINE_POINTS = 16384
DEFAULT_FOLDED_RAW_POINTS = 8192
DEFAULT_CONTEXT_RAW_POINTS = 8192
DEFAULT_MAX_EVENTS = 16
DEFAULT_FOLDED_WINDOW_DURATIONS = 4.0
DEFAULT_CONTEXT_WINDOW_DURATIONS = 12.0
DEFAULT_EVENT_WINDOW_DURATIONS = 4.0
DEFAULT_PSEUDO_MIN_CONFIDENCE = 0.98
DEFAULT_PSEUDO_MIN_MARGIN = 0.50
DEFAULT_PSEUDO_WEIGHT = 0.25


@dataclass(frozen=True)
class TensorConfig:
    apertures: tuple[str, ...] = DEFAULT_APERTURES
    folded_points: int = DEFAULT_FOLDED_POINTS
    context_points: int = DEFAULT_CONTEXT_POINTS
    event_points: int = DEFAULT_EVENT_POINTS
    timeline_points: int = DEFAULT_TIMELINE_POINTS
    folded_raw_points: int = DEFAULT_FOLDED_RAW_POINTS
    context_raw_points: int = DEFAULT_CONTEXT_RAW_POINTS
    max_events: int = DEFAULT_MAX_EVENTS
    folded_window_durations: float = DEFAULT_FOLDED_WINDOW_DURATIONS
    context_window_durations: float = DEFAULT_CONTEXT_WINDOW_DURATIONS
    event_window_durations: float = DEFAULT_EVENT_WINDOW_DURATIONS
    min_event_points: int = 1


@dataclass(frozen=True)
class CnnTrainConfig:
    epochs: int = 80
    batch_size: int = 32
    learning_rate: float = 5.0e-4
    weight_decay: float = 1.0e-4
    dropout: float = 0.20
    early_stop_patience: int = 14
    min_class_count: int = 40
    validation_fraction: float = 0.20
    test_fraction: float = 0.20
    seed: int = 56
    require_cuda: bool = True
    pseudo_min_confidence: float = DEFAULT_PSEUDO_MIN_CONFIDENCE
    pseudo_min_margin: float = DEFAULT_PSEUDO_MIN_MARGIN
    pseudo_weight: float = DEFAULT_PSEUDO_WEIGHT


def _time_to_bjd(time: np.ndarray) -> np.ndarray:
    finite = time[np.isfinite(time)]
    if finite.size and np.nanmedian(finite) < 1.0e5:
        return time + BJDREFI
    return time


def _edges_from_centers(x: np.ndarray) -> np.ndarray:
    if len(x) < 2:
        step = 1.0
    else:
        step = float(np.nanmedian(np.diff(x)))
    return np.concatenate(([x[0] - 0.5 * step], 0.5 * (x[1:] + x[:-1]), [x[-1] + 0.5 * step]))


def _robust_baseline(x: np.ndarray, flux: np.ndarray, *, oot_min_duration: float = 1.5) -> float:
    oot = np.abs(x) > oot_min_duration
    values = flux[oot & np.isfinite(flux)]
    if values.size < 3:
        values = flux[np.isfinite(flux)]
    baseline = float(np.nanmedian(values)) if values.size else float("nan")
    if not np.isfinite(baseline) or abs(baseline) < 1.0e-8:
        return 1.0
    return baseline


def _bin_values(x: np.ndarray, y: np.ndarray, centers: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = _edges_from_centers(centers)
    bin_id = np.digitize(x, edges) - 1
    out = np.full(len(centers), np.nan, dtype=np.float32)
    counts = np.zeros(len(centers), dtype=np.int16)
    valid = (bin_id >= 0) & (bin_id < len(centers)) & np.isfinite(y)
    for idx in np.unique(bin_id[valid]):
        values = y[valid & (bin_id == idx)]
        values = values[np.isfinite(values)]
        if values.size:
            out[idx] = float(np.nanmedian(values))
            counts[idx] = min(int(values.size), np.iinfo(np.int16).max)
    mask = np.isfinite(out)
    return out, mask.astype(np.bool_), counts


def _phase_x_duration(time_bjd: np.ndarray, *, period_d: float, t0_bjd: float, duration_min: float) -> np.ndarray:
    duration_d = duration_min / 1440.0
    if duration_d <= 0 or period_d <= 0:
        return np.full_like(time_bjd, np.nan, dtype=np.float64)
    phase_d = ((time_bjd - t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d
    return phase_d / duration_d


def folded_view(
    *,
    time_bjd: np.ndarray,
    flux_by_aperture: Sequence[np.ndarray],
    quality: np.ndarray,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    n_points: int,
    window_durations: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return value, mask, count, and x-grid arrays for one folded view.

    Channels are small aperture, primary aperture, and primary-minus-small
    after per-aperture local normalization.
    """

    centers = np.linspace(-window_durations, window_durations, n_points, dtype=np.float32)
    x = _phase_x_duration(time_bjd, period_d=period_d, t0_bjd=t0_bjd, duration_min=duration_min)
    good_base = (quality == 0) & np.isfinite(time_bjd) & np.isfinite(x) & (np.abs(x) <= window_durations)

    values: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    counts: list[np.ndarray] = []
    normalized_samples: list[tuple[np.ndarray, np.ndarray]] = []
    for flux in flux_by_aperture:
        flux = np.asarray(flux, dtype=np.float64)
        good = good_base & np.isfinite(flux)
        if np.any(good):
            baseline = _robust_baseline(x[good], flux[good])
            y = flux[good] / baseline - 1.0
            val, mask, count = _bin_values(x[good], y, centers)
            normalized_samples.append((x[good], y))
        else:
            val = np.full(n_points, np.nan, dtype=np.float32)
            mask = np.zeros(n_points, dtype=np.bool_)
            count = np.zeros(n_points, dtype=np.int16)
            normalized_samples.append((np.empty(0), np.empty(0)))
        values.append(val)
        masks.append(mask)
        counts.append(count)

    if len(values) >= 2:
        diff = values[1] - values[0]
        diff_mask = masks[0] & masks[1] & np.isfinite(diff)
        diff_counts = np.minimum(counts[0], counts[1]).astype(np.int16)
        diff = np.where(diff_mask, diff, np.nan).astype(np.float32)
        values.append(diff)
        masks.append(diff_mask)
        counts.append(diff_counts)

    return (
        np.stack(values, axis=0).astype(np.float32),
        np.stack(masks, axis=0).astype(np.bool_),
        np.stack(counts, axis=0).astype(np.int16),
        centers,
    )


FOLD_VIEW_NAMES = ("model", "bls", "bls_half", "bls_double", "bls_secondary")


def fold_view_ephemerides(row: pd.Series) -> tuple[tuple[str, float, float, float], ...]:
    """Return fixed fold views for corrected-model and BLS harmonic context."""

    model_period, model_t0, model_duration, _ = select_model_ephemeris(row)
    bls_period, bls_t0, bls_duration, _ = select_bls_ephemeris(row)
    if np.isfinite(bls_period) and bls_period > 0:
        bls_half = bls_period / 2.0
        bls_double = bls_period * 2.0
        bls_secondary_t0 = bls_t0 + 0.5 * bls_period if np.isfinite(bls_t0) else float("nan")
    else:
        bls_half = float("nan")
        bls_double = float("nan")
        bls_secondary_t0 = float("nan")
    return (
        ("model", model_period, model_t0, model_duration),
        ("bls", bls_period, bls_t0, bls_duration),
        ("bls_half", bls_half, bls_t0, bls_duration),
        ("bls_double", bls_double, bls_t0, bls_duration),
        ("bls_secondary", bls_period, bls_secondary_t0, bls_duration),
    )


def multiview_folded_view(
    *,
    row: pd.Series,
    time_bjd: np.ndarray,
    flux_by_aperture: Sequence[np.ndarray],
    quality: np.ndarray,
    n_points: int,
    window_durations: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    counts: list[np.ndarray] = []
    x_grid = np.linspace(-window_durations, window_durations, n_points, dtype=np.float32)
    for _, period_d, t0_bjd, duration_min in fold_view_ephemerides(row):
        value, mask, count, x_grid = folded_view(
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            period_d=period_d,
            t0_bjd=t0_bjd,
            duration_min=duration_min,
            n_points=n_points,
            window_durations=window_durations,
        )
        values.append(value)
        masks.append(mask)
        counts.append(count)
    return (
        np.concatenate(values, axis=0).astype(np.float32),
        np.concatenate(masks, axis=0).astype(np.bool_),
        np.concatenate(counts, axis=0).astype(np.int16),
        x_grid,
    )


def raw_folded_view(
    *,
    row: pd.Series,
    time_bjd: np.ndarray,
    flux_by_aperture: Sequence[np.ndarray],
    quality: np.ndarray,
    n_points: int,
    window_durations: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return padded unbinned phase-folded point sequences for every fold view."""

    n_base_channels = len(flux_by_aperture) + (1 if len(flux_by_aperture) >= 2 else 0)
    n_channels_per_view = n_base_channels + 1
    values_by_view: list[np.ndarray] = []
    masks_by_view: list[np.ndarray] = []
    raw_counts = np.zeros(len(FOLD_VIEW_NAMES), dtype=np.int32)
    for view_idx, (_, period_d, t0_bjd, duration_min) in enumerate(fold_view_ephemerides(row)):
        values = np.full((n_channels_per_view, n_points), np.nan, dtype=np.float32)
        masks = np.zeros((n_channels_per_view, n_points), dtype=np.bool_)
        x = _phase_x_duration(time_bjd, period_d=period_d, t0_bjd=t0_bjd, duration_min=duration_min)
        good_time = (quality == 0) & np.isfinite(time_bjd) & np.isfinite(x) & (np.abs(x) <= window_durations)
        if not np.any(good_time):
            values_by_view.append(values)
            masks_by_view.append(masks)
            continue
        order_all = np.argsort(x[good_time], kind="mergesort")
        good_idx = np.where(good_time)[0][order_all]
        raw_counts[view_idx] = int(len(good_idx))
        if len(good_idx) > n_points:
            take = np.linspace(0, len(good_idx) - 1, n_points).round().astype(int)
            good_idx = good_idx[take]
        n = len(good_idx)
        x_sel = x[good_idx]

        aperture_values: list[np.ndarray] = []
        aperture_masks: list[np.ndarray] = []
        for flux in flux_by_aperture:
            flux = np.asarray(flux, dtype=np.float64)
            good_flux = good_time & np.isfinite(flux)
            channel = np.full(n, np.nan, dtype=np.float32)
            channel_mask = np.isfinite(flux[good_idx])
            if np.any(good_flux):
                baseline = _robust_baseline(x[good_flux], flux[good_flux])
                channel[channel_mask] = (flux[good_idx[channel_mask]] / baseline - 1.0).astype(np.float32)
            aperture_values.append(channel)
            aperture_masks.append(channel_mask)
        if len(aperture_values) >= 2:
            diff = aperture_values[1] - aperture_values[0]
            diff_mask = aperture_masks[0] & aperture_masks[1] & np.isfinite(diff)
            diff = np.where(diff_mask, diff, np.nan).astype(np.float32)
            aperture_values.append(diff)
            aperture_masks.append(diff_mask)

        for channel_idx, channel in enumerate(aperture_values):
            values[channel_idx, :n] = channel
            masks[channel_idx, :n] = aperture_masks[channel_idx]
        values[-1, :n] = (x_sel / max(window_durations, 1.0e-8)).astype(np.float32)
        masks[-1, :n] = True
        values_by_view.append(values)
        masks_by_view.append(masks)
    return (
        np.concatenate(values_by_view, axis=0).astype(np.float32),
        np.concatenate(masks_by_view, axis=0).astype(np.bool_),
        raw_counts,
    )


def timeline_view(
    *,
    time_bjd: np.ndarray,
    flux_by_aperture: Sequence[np.ndarray],
    quality: np.ndarray,
    n_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return fixed-bin chronological detrended-light-curve channels."""

    centers = np.linspace(0.0, 1.0, n_points, dtype=np.float32)
    n_channels = len(flux_by_aperture) + (1 if len(flux_by_aperture) >= 2 else 0)
    values = np.full((n_channels, n_points), np.nan, dtype=np.float32)
    masks = np.zeros((n_channels, n_points), dtype=np.bool_)
    counts = np.zeros((n_channels, n_points), dtype=np.int16)
    good_time = (quality == 0) & np.isfinite(time_bjd)
    if not np.any(good_time):
        return values, masks, counts, centers
    t_good = time_bjd[good_time]
    t_min = float(np.nanmin(t_good))
    t_span = float(np.nanmax(t_good) - t_min)
    if not np.isfinite(t_span) or t_span <= 0:
        return values, masks, counts, centers
    x_all = (time_bjd - t_min) / t_span

    aperture_values: list[np.ndarray] = []
    aperture_masks: list[np.ndarray] = []
    aperture_counts: list[np.ndarray] = []
    for flux in flux_by_aperture:
        flux = np.asarray(flux, dtype=np.float64)
        good = good_time & np.isfinite(flux)
        if np.any(good):
            baseline = float(np.nanmedian(flux[good]))
            if not np.isfinite(baseline) or abs(baseline) < 1.0e-8:
                baseline = 1.0
            y = flux[good] / baseline - 1.0
            val, mask, count = _bin_values(x_all[good], y, centers)
        else:
            val = np.full(n_points, np.nan, dtype=np.float32)
            mask = np.zeros(n_points, dtype=np.bool_)
            count = np.zeros(n_points, dtype=np.int16)
        aperture_values.append(val)
        aperture_masks.append(mask)
        aperture_counts.append(count)
    if len(aperture_values) >= 2:
        diff = aperture_values[1] - aperture_values[0]
        diff_mask = aperture_masks[0] & aperture_masks[1] & np.isfinite(diff)
        diff_count = np.minimum(aperture_counts[0], aperture_counts[1]).astype(np.int16)
        diff = np.where(diff_mask, diff, np.nan).astype(np.float32)
        aperture_values.append(diff)
        aperture_masks.append(diff_mask)
        aperture_counts.append(diff_count)
    for channel_idx, val in enumerate(aperture_values):
        values[channel_idx, :] = val
        masks[channel_idx, :] = aperture_masks[channel_idx]
        counts[channel_idx, :] = aperture_counts[channel_idx]
    return values, masks, counts, centers


def event_view(
    *,
    time_bjd: np.ndarray,
    flux_by_aperture: Sequence[np.ndarray],
    quality: np.ndarray,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    n_points: int,
    window_durations: float,
    max_events: int,
    min_event_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers = np.linspace(-window_durations, window_durations, n_points, dtype=np.float32)
    n_channels = len(flux_by_aperture) + (1 if len(flux_by_aperture) >= 2 else 0)
    values = np.full((n_channels, max_events, n_points), np.nan, dtype=np.float32)
    masks = np.zeros((n_channels, max_events, n_points), dtype=np.bool_)
    counts = np.zeros((n_channels, max_events, n_points), dtype=np.int16)
    event_centers = np.full(max_events, np.nan, dtype=np.float64)
    duration_d = duration_min / 1440.0
    if period_d <= 0 or duration_d <= 0 or len(time_bjd) == 0:
        return values, masks, counts, event_centers

    finite_time = time_bjd[np.isfinite(time_bjd)]
    if finite_time.size == 0:
        return values, masks, counts, event_centers
    k_min = int(np.floor((float(np.nanmin(finite_time)) - t0_bjd) / period_d)) - 1
    k_max = int(np.ceil((float(np.nanmax(finite_time)) - t0_bjd) / period_d)) + 1
    candidates: list[tuple[int, float, int]] = []
    for k in range(k_min, k_max + 1):
        center = t0_bjd + k * period_d
        x = (time_bjd - center) / duration_d
        n_good = int(np.sum((quality == 0) & np.isfinite(x) & (np.abs(x) <= window_durations)))
        if n_good >= min_event_points:
            candidates.append((k, center, n_good))
    if not candidates:
        return values, masks, counts, event_centers
    candidates = sorted(candidates, key=lambda item: item[1])[:max_events]

    for event_idx, (_, center, _) in enumerate(candidates):
        event_centers[event_idx] = center
        x = (time_bjd - center) / duration_d
        good_base = (quality == 0) & np.isfinite(x) & (np.abs(x) <= window_durations)
        aperture_values: list[np.ndarray] = []
        aperture_masks: list[np.ndarray] = []
        aperture_counts: list[np.ndarray] = []
        for flux in flux_by_aperture:
            flux = np.asarray(flux, dtype=np.float64)
            good = good_base & np.isfinite(flux)
            if np.any(good):
                baseline = _robust_baseline(x[good], flux[good])
                y = flux[good] / baseline - 1.0
                val, mask, count = _bin_values(x[good], y, centers)
            else:
                val = np.full(n_points, np.nan, dtype=np.float32)
                mask = np.zeros(n_points, dtype=np.bool_)
                count = np.zeros(n_points, dtype=np.int16)
            aperture_values.append(val)
            aperture_masks.append(mask)
            aperture_counts.append(count)
        if len(aperture_values) >= 2:
            diff = aperture_values[1] - aperture_values[0]
            diff_mask = aperture_masks[0] & aperture_masks[1] & np.isfinite(diff)
            diff_count = np.minimum(aperture_counts[0], aperture_counts[1]).astype(np.int16)
            diff = np.where(diff_mask, diff, np.nan).astype(np.float32)
            aperture_values.append(diff)
            aperture_masks.append(diff_mask)
            aperture_counts.append(diff_count)
        for channel_idx, val in enumerate(aperture_values):
            values[channel_idx, event_idx, :] = val
            masks[channel_idx, event_idx, :] = aperture_masks[channel_idx]
            counts[channel_idx, event_idx, :] = aperture_counts[channel_idx]
    return values, masks, counts, event_centers


def _tensor_status(mask: np.ndarray) -> str:
    if mask.size == 0 or not bool(mask.any()):
        return "no_observed_bins"
    per_channel = mask.reshape(mask.shape[0], -1).sum(axis=1)
    if per_channel[0] == 0:
        return "missing_primary_small_channel"
    return "ok"


def _read_tensor_lc(
    row: pd.Series,
    *,
    compact_lc_h5: Path | None,
    hlsp_root: Path | None,
    injection_h5_override: Path | None,
    apertures: Sequence[str],
) -> HLSPLightCurve | None:
    return _read_lc_for_row(
        row,
        apertures=apertures,
        compact_lc_h5=compact_lc_h5,
        hlsp_root=hlsp_root,
        injection_h5_override=injection_h5_override,
    )


def build_recovery50_cnn_tensors(
    *,
    training_table: Path,
    out_dir: Path,
    compact_lc_h5: Path | None,
    hlsp_root: Path | None,
    injection_h5_override: Path | None,
    config: TensorConfig = TensorConfig(),
    max_rows: int | None = None,
    progress_every: int = 50,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_table(training_table).copy()
    assert_adp_only_training_frame(rows)
    validate_adp_only_apertures(config.apertures)
    if "row_id" not in rows:
        rows.insert(0, "row_id", np.arange(len(rows), dtype=int))
    if max_rows is not None:
        rows = rows.head(max_rows).copy()

    folded_values: list[np.ndarray] = []
    folded_masks: list[np.ndarray] = []
    folded_counts: list[np.ndarray] = []
    context_values: list[np.ndarray] = []
    context_masks: list[np.ndarray] = []
    context_counts: list[np.ndarray] = []
    raw_folded_values: list[np.ndarray] = []
    raw_folded_masks: list[np.ndarray] = []
    raw_folded_counts: list[np.ndarray] = []
    raw_context_values: list[np.ndarray] = []
    raw_context_masks: list[np.ndarray] = []
    raw_context_counts: list[np.ndarray] = []
    event_values: list[np.ndarray] = []
    event_masks: list[np.ndarray] = []
    event_counts: list[np.ndarray] = []
    event_centers: list[np.ndarray] = []
    timeline_values: list[np.ndarray] = []
    timeline_masks: list[np.ndarray] = []
    timeline_counts: list[np.ndarray] = []
    records: list[dict[str, Any]] = []
    status_counts: dict[str, int] = {}

    for processed, (_, row) in enumerate(rows.iterrows(), start=1):
        period_d, t0_bjd, duration_min, ephemeris_source = select_model_ephemeris(row)
        rec = {
            "tensor_index": len(records),
            "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            "tensor_apertures": ADP_ONLY_APERTURE_SIGNATURE,
            "row_id": int(row.get("row_id", -1)),
            "review_id": str(row.get("review_id", "")),
            "tic": int(row.get("tic", -1)),
            "source_kind": str(row.get("source_kind", "")),
            "human_label": str(row.get("human_label", "")),
            "main_teacher_target": str(row.get("main_teacher_target", "")),
            "main_teacher_include": bool(row.get("main_teacher_include", False)),
            "training_split": str(row.get("training_split", "")),
            "period_d": period_d,
            "t0_bjd": t0_bjd,
            "duration_min": duration_min,
            "model_period_d": period_d,
            "model_t0_bjd": t0_bjd,
            "model_duration_min": duration_min,
            "model_ephemeris_source": ephemeris_source,
            "harmonic_suspect": bool(row.get("harmonic_suspect", False)),
            "refold_factor": _safe_float(row.get("refold_factor")),
            "refold_period_d": _safe_float(row.get("refold_period_d")),
            "ephemeris_status": str(row.get("ephemeris_status", "")),
            "queue_period_d": _safe_float(row.get("period_d")),
            "queue_t0_bjd": _safe_float(row.get("t0_bjd")),
            "queue_duration_min": _safe_float(row.get("duration_min")),
            "display_period_d": _safe_float(row.get("display_period_d")),
            "display_t0_bjd": _safe_float(row.get("display_t0_bjd")),
            "display_duration_min": _safe_float(row.get("display_duration_min")),
            "display_ephemeris_source": str(row.get("display_ephemeris_source", "")),
            "tmag": _safe_float(row.get("tmag")),
        }
        if not all(np.isfinite([period_d, t0_bjd, duration_min])) or period_d <= 0 or duration_min <= 0:
            rec["tensor_status"] = "invalid_ephemeris"
            status_counts[rec["tensor_status"]] = status_counts.get(rec["tensor_status"], 0) + 1
            continue
        lc = _read_tensor_lc(
            row,
            compact_lc_h5=compact_lc_h5,
            hlsp_root=hlsp_root,
            injection_h5_override=injection_h5_override,
            apertures=config.apertures,
        )
        if lc is None:
            rec["tensor_status"] = "missing_light_curve"
            status_counts[rec["tensor_status"]] = status_counts.get(rec["tensor_status"], 0) + 1
            continue
        primary_aperture = config.apertures[0]
        if primary_aperture not in lc.flux:
            rec["tensor_status"] = "missing_primary_aperture:" + primary_aperture
            status_counts[rec["tensor_status"]] = status_counts.get(rec["tensor_status"], 0) + 1
            continue

        time_bjd = _time_to_bjd(np.asarray(lc.time, dtype=np.float64))
        quality = np.asarray(lc.quality, dtype=np.int32)
        missing_supplemental = [ap for ap in config.apertures[1:] if ap not in lc.flux]
        rec["missing_supplemental_apertures"] = ",".join(missing_supplemental)
        flux_by_aperture = [
            np.asarray(lc.flux[ap], dtype=np.float64)
            if ap in lc.flux
            else np.full_like(time_bjd, np.nan, dtype=np.float64)
            for ap in config.apertures
        ]

        folded, folded_mask, folded_count, folded_x = multiview_folded_view(
            row=row,
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            n_points=config.folded_points,
            window_durations=config.folded_window_durations,
        )
        context, context_mask, context_count, context_x = multiview_folded_view(
            row=row,
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            n_points=config.context_points,
            window_durations=config.context_window_durations,
        )
        raw_folded, raw_folded_mask, raw_folded_count = raw_folded_view(
            row=row,
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            n_points=config.folded_raw_points,
            window_durations=config.folded_window_durations,
        )
        raw_context, raw_context_mask, raw_context_count = raw_folded_view(
            row=row,
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            n_points=config.context_raw_points,
            window_durations=config.context_window_durations,
        )
        events, event_mask, event_count, centers = event_view(
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            period_d=period_d,
            t0_bjd=t0_bjd,
            duration_min=duration_min,
            n_points=config.event_points,
            window_durations=config.event_window_durations,
            max_events=config.max_events,
            min_event_points=config.min_event_points,
        )
        timeline, timeline_mask, timeline_count, timeline_x = timeline_view(
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            n_points=config.timeline_points,
        )
        status = _tensor_status(folded_mask)
        rec["tensor_status"] = status
        rec["n_event_windows"] = int(np.isfinite(centers).sum())
        rec["folded_observed_fraction"] = float(folded_mask.sum() / folded_mask.size) if folded_mask.size else 0.0
        rec["context_observed_fraction"] = float(context_mask.sum() / context_mask.size) if context_mask.size else 0.0
        rec["raw_folded_observed_fraction"] = (
            float(raw_folded_mask.sum() / raw_folded_mask.size) if raw_folded_mask.size else 0.0
        )
        rec["raw_context_observed_fraction"] = (
            float(raw_context_mask.sum() / raw_context_mask.size) if raw_context_mask.size else 0.0
        )
        rec["raw_folded_max_samples"] = int(raw_folded_count.max()) if raw_folded_count.size else 0
        rec["raw_context_max_samples"] = int(raw_context_count.max()) if raw_context_count.size else 0
        rec["raw_folded_truncated"] = bool(raw_folded_count.max() > config.folded_raw_points) if raw_folded_count.size else False
        rec["raw_context_truncated"] = bool(raw_context_count.max() > config.context_raw_points) if raw_context_count.size else False
        rec["event_observed_fraction"] = float(event_mask.sum() / event_mask.size) if event_mask.size else 0.0
        rec["timeline_observed_fraction"] = float(timeline_mask.sum() / timeline_mask.size) if timeline_mask.size else 0.0
        if status != "ok":
            status_counts[status] = status_counts.get(status, 0) + 1
            continue

        folded_values.append(folded)
        folded_masks.append(folded_mask)
        folded_counts.append(folded_count)
        context_values.append(context)
        context_masks.append(context_mask)
        context_counts.append(context_count)
        raw_folded_values.append(raw_folded)
        raw_folded_masks.append(raw_folded_mask)
        raw_folded_counts.append(raw_folded_count)
        raw_context_values.append(raw_context)
        raw_context_masks.append(raw_context_mask)
        raw_context_counts.append(raw_context_count)
        event_values.append(events)
        event_masks.append(event_mask)
        event_counts.append(event_count)
        event_centers.append(centers)
        timeline_values.append(timeline)
        timeline_masks.append(timeline_mask)
        timeline_counts.append(timeline_count)
        records.append(rec)
        status_counts[status] = status_counts.get(status, 0) + 1
        if progress_every > 0 and (processed % progress_every == 0 or processed == len(rows)):
            print(f"[cnn-tensor] processed={processed:,}/{len(rows):,}; accepted={len(records):,}", flush=True)

    n_base_channels = len(config.apertures) + (1 if len(config.apertures) >= 2 else 0)
    n_fold_channels = n_base_channels * len(FOLD_VIEW_NAMES)
    n_raw_fold_channels = (n_base_channels + 1) * len(FOLD_VIEW_NAMES)
    folded_arr = (
        np.stack(folded_values, axis=0).astype(np.float32)
        if folded_values
        else np.empty((0, n_fold_channels, config.folded_points), dtype=np.float32)
    )
    folded_mask_arr = (
        np.stack(folded_masks, axis=0).astype(np.bool_)
        if folded_masks
        else np.empty((0, n_fold_channels, config.folded_points), dtype=np.bool_)
    )
    folded_count_arr = (
        np.stack(folded_counts, axis=0).astype(np.int16)
        if folded_counts
        else np.empty((0, n_fold_channels, config.folded_points), dtype=np.int16)
    )
    context_arr = (
        np.stack(context_values, axis=0).astype(np.float32)
        if context_values
        else np.empty((0, n_fold_channels, config.context_points), dtype=np.float32)
    )
    context_mask_arr = (
        np.stack(context_masks, axis=0).astype(np.bool_)
        if context_masks
        else np.empty((0, n_fold_channels, config.context_points), dtype=np.bool_)
    )
    context_count_arr = (
        np.stack(context_counts, axis=0).astype(np.int16)
        if context_counts
        else np.empty((0, n_fold_channels, config.context_points), dtype=np.int16)
    )
    raw_folded_arr = (
        np.stack(raw_folded_values, axis=0).astype(np.float16)
        if raw_folded_values
        else np.empty((0, n_raw_fold_channels, config.folded_raw_points), dtype=np.float16)
    )
    raw_folded_mask_arr = (
        np.stack(raw_folded_masks, axis=0).astype(np.bool_)
        if raw_folded_masks
        else np.empty((0, n_raw_fold_channels, config.folded_raw_points), dtype=np.bool_)
    )
    raw_folded_count_arr = (
        np.stack(raw_folded_counts, axis=0).astype(np.int32)
        if raw_folded_counts
        else np.empty((0, len(FOLD_VIEW_NAMES)), dtype=np.int32)
    )
    raw_context_arr = (
        np.stack(raw_context_values, axis=0).astype(np.float16)
        if raw_context_values
        else np.empty((0, n_raw_fold_channels, config.context_raw_points), dtype=np.float16)
    )
    raw_context_mask_arr = (
        np.stack(raw_context_masks, axis=0).astype(np.bool_)
        if raw_context_masks
        else np.empty((0, n_raw_fold_channels, config.context_raw_points), dtype=np.bool_)
    )
    raw_context_count_arr = (
        np.stack(raw_context_counts, axis=0).astype(np.int32)
        if raw_context_counts
        else np.empty((0, len(FOLD_VIEW_NAMES)), dtype=np.int32)
    )
    event_arr = (
        np.stack(event_values, axis=0).astype(np.float32)
        if event_values
        else np.empty((0, n_base_channels, config.max_events, config.event_points), dtype=np.float32)
    )
    event_mask_arr = (
        np.stack(event_masks, axis=0).astype(np.bool_)
        if event_masks
        else np.empty((0, n_base_channels, config.max_events, config.event_points), dtype=np.bool_)
    )
    event_count_arr = (
        np.stack(event_counts, axis=0).astype(np.int16)
        if event_counts
        else np.empty((0, n_base_channels, config.max_events, config.event_points), dtype=np.int16)
    )
    event_centers_arr = (
        np.stack(event_centers, axis=0).astype(np.float64)
        if event_centers
        else np.empty((0, config.max_events), dtype=np.float64)
    )
    timeline_arr = (
        np.stack(timeline_values, axis=0).astype(np.float32)
        if timeline_values
        else np.empty((0, n_base_channels, config.timeline_points), dtype=np.float32)
    )
    timeline_mask_arr = (
        np.stack(timeline_masks, axis=0).astype(np.bool_)
        if timeline_masks
        else np.empty((0, n_base_channels, config.timeline_points), dtype=np.bool_)
    )
    timeline_count_arr = (
        np.stack(timeline_counts, axis=0).astype(np.int16)
        if timeline_counts
        else np.empty((0, n_base_channels, config.timeline_points), dtype=np.int16)
    )

    out_npz = out_dir / "recovery50_cnn_tensors.npz"
    np.savez_compressed(
        out_npz,
        folded=folded_arr,
        folded_mask=folded_mask_arr,
        folded_counts=folded_count_arr,
        context=context_arr,
        context_mask=context_mask_arr,
        context_counts=context_count_arr,
        raw_folded=raw_folded_arr,
        raw_folded_mask=raw_folded_mask_arr,
        raw_folded_counts=raw_folded_count_arr,
        raw_context=raw_context_arr,
        raw_context_mask=raw_context_mask_arr,
        raw_context_counts=raw_context_count_arr,
        events=event_arr,
        events_mask=event_mask_arr,
        events_counts=event_count_arr,
        event_centers_bjd=event_centers_arr,
        timeline=timeline_arr,
        timeline_mask=timeline_mask_arr,
        timeline_counts=timeline_count_arr,
        folded_x_duration=np.linspace(
            -config.folded_window_durations,
            config.folded_window_durations,
            config.folded_points,
            dtype=np.float32,
        ),
        context_x_duration=np.linspace(
            -config.context_window_durations,
            config.context_window_durations,
            config.context_points,
            dtype=np.float32,
        ),
        event_x_duration=np.linspace(
            -config.event_window_durations,
            config.event_window_durations,
            config.event_points,
            dtype=np.float32,
        ),
        timeline_x_fraction=np.linspace(0.0, 1.0, config.timeline_points, dtype=np.float32),
        apertures=np.asarray(config.apertures, dtype=object),
        derived_channels=np.asarray(["primary_minus_small"] if len(config.apertures) >= 2 else [], dtype=object),
        fold_view_names=np.asarray(FOLD_VIEW_NAMES, dtype=object),
    )
    row_table = pd.DataFrame(records)
    row_csv = out_dir / "recovery50_cnn_tensor_rows.csv"
    row_table.to_csv(row_csv, index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "training_table": str(training_table),
        "compact_lc_h5": str(compact_lc_h5) if compact_lc_h5 is not None else "",
        "hlsp_root": str(hlsp_root) if hlsp_root is not None else "",
        "injection_h5_override": str(injection_h5_override) if injection_h5_override is not None else "",
        "out_dir": str(out_dir),
        "n_input_rows": int(len(rows)),
        "n_tensor_rows": int(len(row_table)),
        "folded_shape": list(folded_arr.shape),
        "context_shape": list(context_arr.shape),
        "raw_folded_shape": list(raw_folded_arr.shape),
        "raw_context_shape": list(raw_context_arr.shape),
        "events_shape": list(event_arr.shape),
        "timeline_shape": list(timeline_arr.shape),
        "status_counts": {str(k): int(v) for k, v in sorted(status_counts.items())},
        "raw_folded_truncated_rows": int(row_table.get("raw_folded_truncated", pd.Series(dtype=bool)).sum()) if len(row_table) else 0,
        "raw_context_truncated_rows": int(row_table.get("raw_context_truncated", pd.Series(dtype=bool)).sum()) if len(row_table) else 0,
        "apertures": list(config.apertures),
        "fold_view_names": list(FOLD_VIEW_NAMES),
        "tensor_config": asdict(config),
        "observable_channel_policy": "primary small aperture is required; supplemental aperture and difference channels are masked when absent; injected original/pre-injection arrays excluded",
        "ephemeris_policy": "folded/context tensors concatenate model/corrected, BLS, BLS/2, 2*BLS, and BLS secondary views; events use model/corrected ephemeris; timeline is chronological detrended flux before folding",
        "outputs": {"npz": str(out_npz), "rows": str(row_csv), "summary": str(out_dir / "summary.json")},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    (out_dir / "summary.md").write_text(render_tensor_markdown(summary))
    return summary


def render_tensor_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 Recovery50 CNN Tensor Build",
        "",
        f"- input rows: `{summary['n_input_rows']}`",
        f"- tensor rows: `{summary['n_tensor_rows']}`",
        f"- folded shape: `{tuple(summary['folded_shape'])}`",
        f"- context shape: `{tuple(summary['context_shape'])}`",
        f"- raw folded shape: `{tuple(summary['raw_folded_shape'])}`",
        f"- raw context shape: `{tuple(summary['raw_context_shape'])}`",
        f"- events shape: `{tuple(summary['events_shape'])}`",
        f"- timeline shape: `{tuple(summary['timeline_shape'])}`",
        f"- apertures: `{', '.join(summary['apertures'])}`",
        f"- fold views: `{', '.join(summary['fold_view_names'])}`",
        f"- status counts: `{summary['status_counts']}`",
        f"- raw folded truncated rows: `{summary['raw_folded_truncated_rows']}`",
        f"- raw context truncated rows: `{summary['raw_context_truncated_rows']}`",
        f"- channel policy: `{summary['observable_channel_policy']}`",
        f"- ephemeris policy: `{summary['ephemeris_policy']}`",
        "",
    ]
    return "\n".join(lines)


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=True) as npz:
        return {key: np.asarray(npz[key]) for key in npz.files}


def _prepare_value_mask(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    value = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    mask_value = mask.astype(np.float32)
    return np.concatenate([value, mask_value], axis=1).astype(np.float32)


def _balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray, n_classes: int) -> float:
    recalls = []
    for cls in range(n_classes):
        idx = y_true == cls
        if np.any(idx):
            recalls.append(float(np.mean(y_pred[idx] == cls)))
    return float(np.mean(recalls)) if recalls else float("nan")


def _grouped_splits(
    rows: pd.DataFrame,
    *,
    classes: Sequence[str],
    cfg: CnnTrainConfig,
) -> pd.Series:
    rng = np.random.default_rng(cfg.seed)
    split_by_group: dict[Any, str] = {}
    group_frame = (
        rows.groupby("tic", dropna=False)
        .agg(main_teacher_target=("main_teacher_target", lambda s: s.value_counts().index[0]))
        .reset_index()
    )
    for cls in classes:
        groups = group_frame.loc[group_frame["main_teacher_target"].eq(cls), "tic"].to_numpy().copy()
        rng.shuffle(groups)
        n = len(groups)
        if n < 3:
            for group in groups:
                split_by_group[group] = "train"
            continue
        n_test = max(1, int(round(cfg.test_fraction * n))) if cfg.test_fraction > 0 else 0
        n_val = max(1, int(round(cfg.validation_fraction * n))) if cfg.validation_fraction > 0 and n - n_test >= 3 else 0
        n_test = min(n_test, max(0, n - 1))
        n_val = min(n_val, max(0, n - n_test - 1))
        for group in groups[:n_test]:
            split_by_group[group] = "test"
        for group in groups[n_test:n_test + n_val]:
            split_by_group[group] = "validation"
        for group in groups[n_test + n_val:]:
            split_by_group[group] = "train"
    return rows["tic"].map(split_by_group).fillna("train").astype(str)


def _metadata_matrix(
    *,
    rows: pd.DataFrame,
    feature_table: pd.DataFrame,
    train_mask: np.ndarray,
) -> tuple[np.ndarray, list[str], dict[str, list[float]]]:
    join_cols = (
        ["review_id"]
        if "review_id" in rows.columns and "review_id" in feature_table.columns
        else ["row_id"] if "row_id" in rows.columns and "row_id" in feature_table.columns else []
    )
    if not join_cols:
        return np.empty((len(rows), 0), dtype=np.float32), [], {"center": [], "scale": []}
    merged = rows.loc[:, join_cols].merge(feature_table, on=join_cols, how="left")
    columns = metadata_feature_columns(merged)
    columns = [col for col in columns if col in merged.columns]
    leaks = leakage_columns(columns)
    if leaks:
        raise ValueError(f"metadata columns contain leakage: {leaks}")
    if not columns:
        return np.empty((len(rows), 0), dtype=np.float32), [], {"center": [], "scale": []}
    x = merged.loc[:, columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    center = np.nanmedian(x[train_mask], axis=0)
    center = np.where(np.isfinite(center), center, 0.0)
    scale = np.nanstd(x[train_mask], axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
    x = np.where(np.isfinite(x), x, center)
    x = ((x - center) / scale).astype(np.float32)
    return x, columns, {"center": center.astype(float).tolist(), "scale": scale.astype(float).tolist()}


def _metadata_matrix_from_saved_norm(
    *,
    rows: pd.DataFrame,
    feature_table: pd.DataFrame,
    columns: Sequence[str],
    norm: dict[str, Sequence[float]],
) -> np.ndarray:
    if not columns:
        return np.empty((len(rows), 0), dtype=np.float32)
    join_cols = (
        ["review_id"]
        if "review_id" in rows.columns and "review_id" in feature_table.columns
        else ["row_id"] if "row_id" in rows.columns and "row_id" in feature_table.columns else []
    )
    if join_cols:
        merged = rows.loc[:, join_cols].merge(feature_table, on=join_cols, how="left")
    else:
        merged = rows.copy()
    leaks = leakage_columns(columns)
    if leaks:
        raise ValueError(f"metadata columns contain leakage: {leaks}")
    for column in columns:
        if column not in merged:
            merged[column] = np.nan
    x = merged.loc[:, list(columns)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    center = np.asarray(norm.get("center", []), dtype=np.float32)
    scale = np.asarray(norm.get("scale", []), dtype=np.float32)
    if center.size != len(columns):
        center = np.zeros(len(columns), dtype=np.float32)
    if scale.size != len(columns):
        scale = np.ones(len(columns), dtype=np.float32)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0).astype(np.float32)
    center = np.where(np.isfinite(center), center, 0.0).astype(np.float32)
    x = np.where(np.isfinite(x), x, center)
    return ((x - center) / scale).astype(np.float32)


def _build_recovery50_cnn_model(
    *,
    folded_channels: int,
    raw_folded_channels: int,
    raw_context_channels: int,
    event_channels: int,
    timeline_channels: int,
    metadata_dim: int,
    n_classes: int,
    dropout: float,
) -> Any:
    import torch
    from torch import nn

    class Recovery50CNN(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.folded_stem = nn.Sequential(
                nn.Conv1d(folded_channels, 32, kernel_size=9, padding=4),
                nn.BatchNorm1d(32),
                nn.ReLU(),
                nn.Conv1d(32, 32, kernel_size=7, padding=3),
                nn.BatchNorm1d(32),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.context_stem = nn.Sequential(
                nn.Conv1d(folded_channels, 16, kernel_size=11, padding=5),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=7, padding=3),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.raw_folded_stem = nn.Sequential(
                nn.Conv1d(raw_folded_channels, 16, kernel_size=15, padding=7),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=9, padding=4),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.raw_context_stem = nn.Sequential(
                nn.Conv1d(raw_context_channels, 16, kernel_size=21, padding=10),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=11, padding=5),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.event_stem = nn.Sequential(
                nn.Conv2d(event_channels, 16, kernel_size=(3, 9), padding=(1, 4)),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.Conv2d(16, 16, kernel_size=(3, 7), padding=(1, 3)),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool2d((1, 1)),
            )
            self.timeline_stem = nn.Sequential(
                nn.Conv1d(timeline_channels, 16, kernel_size=17, padding=8),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=9, padding=4),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            hidden_in = 32 + 16 + 16 + 16 + 16 + 16 + metadata_dim
            self.head = nn.Sequential(
                nn.Dropout(dropout),
                nn.Linear(hidden_in, 64),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(64, n_classes),
            )

        def forward(
            self,
            folded_x: Any,
            context_x: Any,
            raw_folded_x: Any,
            raw_context_x: Any,
            event_x: Any,
            timeline_x: Any,
            metadata: Any,
        ) -> Any:
            f = self.folded_stem(folded_x).flatten(1)
            c = self.context_stem(context_x).flatten(1)
            rf = self.raw_folded_stem(raw_folded_x).flatten(1)
            rc = self.raw_context_stem(raw_context_x).flatten(1)
            e = self.event_stem(event_x).flatten(1)
            t = self.timeline_stem(timeline_x).flatten(1)
            parts = [f, c, rf, rc, e, t]
            if metadata.shape[1] > 0:
                parts.append(metadata)
            return self.head(torch.cat(parts, dim=1))

    return Recovery50CNN()


def _classification_report_arrays(
    truth: np.ndarray,
    pred: np.ndarray,
    *,
    labels: Sequence[str],
) -> dict[str, Any]:
    truth_s = pd.Series([labels[int(idx)] for idx in truth])
    pred_s = pd.Series([labels[int(idx)] for idx in pred])
    return _confusion_metrics(truth_s, pred_s, labels)


def _ece(prob: np.ndarray, y_true: np.ndarray, n_bins: int = 10) -> dict[str, Any]:
    if prob.size == 0:
        return {"ece": float("nan"), "bins": []}
    conf = prob.max(axis=1)
    pred = prob.argmax(axis=1)
    correct = pred == y_true
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    bins: list[dict[str, Any]] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = (conf >= lo) & (conf < hi if hi < 1.0 else conf <= hi)
        n = int(idx.sum())
        if n == 0:
            bins.append({"lo": float(lo), "hi": float(hi), "n": 0})
            continue
        acc = float(np.mean(correct[idx]))
        avg_conf = float(np.mean(conf[idx]))
        ece += (n / len(conf)) * abs(acc - avg_conf)
        bins.append({"lo": float(lo), "hi": float(hi), "n": n, "accuracy": acc, "confidence": avg_conf})
    return {"ece": float(ece), "bins": bins}


def train_recovery50_cnn_teacher(
    *,
    tensor_npz: Path,
    tensor_rows: Path,
    training_table: Path,
    out_dir: Path,
    metrics_tables: Sequence[Path] = (),
    cfg: CnnTrainConfig = CnnTrainConfig(),
) -> dict[str, Any]:
    import torch
    from torch import nn
    from torch.utils.data import DataLoader, TensorDataset

    out_dir.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    tensors = _load_npz(tensor_npz)
    rows = read_table(tensor_rows).copy()
    assert_adp_only_tensor_rows(rows)
    training = load_feature_table(training_table=training_table, metrics_tables=tuple(metrics_tables))
    assert_adp_only_training_frame(training)
    classes = select_teacher_classes(training, cfg.min_class_count)
    if len(classes) < 2:
        raise ValueError(f"need at least two teacher classes, got {classes}")

    if "review_id" in rows and "review_id" in training:
        join_columns = ["review_id"]
    elif "row_id" in rows and "row_id" in training:
        join_columns = ["row_id"]
    else:
        raise KeyError("tensor rows and training table have no shared row identifier")
    target_columns = [
        column
        for column in (
            *join_columns,
            "human_label",
            "main_teacher_target",
            "main_teacher_include",
            "training_split",
            "adp_only_review_ok",
        )
        if column in training
    ]
    active_targets = training.loc[:, target_columns].drop_duplicates(join_columns, keep="last")
    stale_target_columns = [
        column
        for column in target_columns
        if column not in join_columns and column in rows
    ]
    rows = rows.drop(columns=stale_target_columns).merge(
        active_targets,
        on=join_columns,
        how="left",
        validate="one_to_one",
    )
    active_include = (
        rows["main_teacher_include"]
        if rows["main_teacher_include"].dtype == bool
        else rows["main_teacher_include"].fillna(False).astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})
    )
    work = rows[active_include & rows["main_teacher_target"].isin(classes)].copy().reset_index(drop=True)
    if work.empty:
        raise ValueError("no tensor rows match selected teacher classes")
    split = _grouped_splits(work, classes=classes, cfg=cfg)
    work["cnn_training_split"] = split
    split_class_counts = pd.crosstab(work["main_teacher_target"], split)
    for cls in classes:
        if cfg.validation_fraction > 0 and int(split_class_counts.get("validation", pd.Series(dtype=int)).get(cls, 0)) == 0:
            raise ValueError(f"CNN validation split has no rows for class {cls}")
        if cfg.test_fraction > 0 and int(split_class_counts.get("test", pd.Series(dtype=int)).get(cls, 0)) == 0:
            raise ValueError(f"CNN test split has no rows for class {cls}")
    train_mask = split.to_numpy() == "train"
    if len(np.unique(work.loc[train_mask, "main_teacher_target"])) < 2:
        raise ValueError("CNN train split has fewer than two classes")

    tensor_idx = work["tensor_index"].to_numpy(dtype=int)
    folded = _prepare_value_mask(tensors["folded"][tensor_idx], tensors["folded_mask"][tensor_idx])
    context = _prepare_value_mask(tensors["context"][tensor_idx], tensors["context_mask"][tensor_idx])
    if "raw_folded" in tensors and "raw_folded_mask" in tensors:
        raw_folded = _prepare_value_mask(tensors["raw_folded"][tensor_idx], tensors["raw_folded_mask"][tensor_idx])
    else:
        raw_values = np.zeros((len(work), 1, DEFAULT_FOLDED_RAW_POINTS), dtype=np.float32)
        raw_masks = np.zeros_like(raw_values, dtype=np.bool_)
        raw_folded = _prepare_value_mask(raw_values, raw_masks)
    if "raw_context" in tensors and "raw_context_mask" in tensors:
        raw_context = _prepare_value_mask(tensors["raw_context"][tensor_idx], tensors["raw_context_mask"][tensor_idx])
    else:
        raw_values = np.zeros((len(work), 1, DEFAULT_CONTEXT_RAW_POINTS), dtype=np.float32)
        raw_masks = np.zeros_like(raw_values, dtype=np.bool_)
        raw_context = _prepare_value_mask(raw_values, raw_masks)
    events = _prepare_value_mask(tensors["events"][tensor_idx], tensors["events_mask"][tensor_idx])
    if "timeline" in tensors and "timeline_mask" in tensors:
        timeline = _prepare_value_mask(tensors["timeline"][tensor_idx], tensors["timeline_mask"][tensor_idx])
    else:
        timeline_values = np.zeros((len(work), 1, DEFAULT_TIMELINE_POINTS), dtype=np.float32)
        timeline_masks = np.zeros_like(timeline_values, dtype=np.bool_)
        timeline = _prepare_value_mask(timeline_values, timeline_masks)
    y_labels = work["main_teacher_target"].astype(str).to_numpy()
    class_to_id = {label: idx for idx, label in enumerate(classes)}
    y = np.asarray([class_to_id[label] for label in y_labels], dtype=np.int64)

    feature_table = load_feature_table(training_table=training_table, metrics_tables=tuple(metrics_tables))
    metadata_x, metadata_columns, metadata_norm = _metadata_matrix(
        rows=work,
        feature_table=feature_table,
        train_mask=train_mask,
    )

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if cfg.require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required but torch.cuda.is_available() is false")

    train_idx = np.where(split.to_numpy() == "train")[0]
    class_counts = np.bincount(y[train_idx], minlength=len(classes)).astype(np.float32)
    class_weights = class_counts.sum() / np.maximum(class_counts, 1.0)
    class_weights = class_weights / np.nanmean(class_weights)

    class Recovery50CNN(nn.Module):
        def __init__(
            self,
            folded_channels: int,
            raw_folded_channels: int,
            raw_context_channels: int,
            event_channels: int,
            timeline_channels: int,
            metadata_dim: int,
            n_classes: int,
        ):
            super().__init__()
            self.folded_stem = nn.Sequential(
                nn.Conv1d(folded_channels, 32, kernel_size=9, padding=4),
                nn.BatchNorm1d(32),
                nn.ReLU(),
                nn.Conv1d(32, 32, kernel_size=7, padding=3),
                nn.BatchNorm1d(32),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.context_stem = nn.Sequential(
                nn.Conv1d(folded_channels, 16, kernel_size=11, padding=5),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=7, padding=3),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.raw_folded_stem = nn.Sequential(
                nn.Conv1d(raw_folded_channels, 16, kernel_size=15, padding=7),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=9, padding=4),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.raw_context_stem = nn.Sequential(
                nn.Conv1d(raw_context_channels, 16, kernel_size=21, padding=10),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=11, padding=5),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            self.event_stem = nn.Sequential(
                nn.Conv2d(event_channels, 16, kernel_size=(3, 9), padding=(1, 4)),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.Conv2d(16, 16, kernel_size=(3, 7), padding=(1, 3)),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool2d((1, 1)),
            )
            self.timeline_stem = nn.Sequential(
                nn.Conv1d(timeline_channels, 16, kernel_size=17, padding=8),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.Conv1d(16, 16, kernel_size=9, padding=4),
                nn.BatchNorm1d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(1),
            )
            hidden_in = 32 + 16 + 16 + 16 + 16 + 16 + metadata_dim
            self.head = nn.Sequential(
                nn.Dropout(cfg.dropout),
                nn.Linear(hidden_in, 64),
                nn.ReLU(),
                nn.Dropout(cfg.dropout),
                nn.Linear(64, n_classes),
            )

        def forward(
            self,
            folded_x: Any,
            context_x: Any,
            raw_folded_x: Any,
            raw_context_x: Any,
            event_x: Any,
            timeline_x: Any,
            metadata: Any,
        ) -> Any:
            f = self.folded_stem(folded_x).flatten(1)
            c = self.context_stem(context_x).flatten(1)
            rf = self.raw_folded_stem(raw_folded_x).flatten(1)
            rc = self.raw_context_stem(raw_context_x).flatten(1)
            e = self.event_stem(event_x).flatten(1)
            t = self.timeline_stem(timeline_x).flatten(1)
            parts = [f, c, rf, rc, e, t]
            if metadata.shape[1] > 0:
                parts.append(metadata)
            return self.head(torch.cat(parts, dim=1))

    def run_profile(profile: str, use_metadata: bool) -> dict[str, Any]:
        profile_dir = out_dir / profile
        profile_dir.mkdir(parents=True, exist_ok=True)
        meta = metadata_x if use_metadata else np.empty((len(work), 0), dtype=np.float32)
        model = Recovery50CNN(
            folded_channels=folded.shape[1],
            raw_folded_channels=raw_folded.shape[1],
            raw_context_channels=raw_context.shape[1],
            event_channels=events.shape[1],
            timeline_channels=timeline.shape[1],
            metadata_dim=meta.shape[1],
            n_classes=len(classes),
        ).to(device)
        opt = torch.optim.AdamW(model.parameters(), lr=cfg.learning_rate, weight_decay=cfg.weight_decay)
        loss_fn = nn.CrossEntropyLoss(weight=torch.tensor(class_weights, dtype=torch.float32, device=device))
        ds = TensorDataset(
            torch.from_numpy(folded[train_idx]).float(),
            torch.from_numpy(context[train_idx]).float(),
            torch.from_numpy(raw_folded[train_idx]).float(),
            torch.from_numpy(raw_context[train_idx]).float(),
            torch.from_numpy(events[train_idx]).float(),
            torch.from_numpy(timeline[train_idx]).float(),
            torch.from_numpy(meta[train_idx]).float(),
            torch.from_numpy(y[train_idx]).long(),
        )
        loader = DataLoader(
            ds,
            batch_size=cfg.batch_size,
            shuffle=True,
            generator=torch.Generator().manual_seed(cfg.seed),
        )
        best_state = None
        best_val = -np.inf
        best_epoch = 0
        stale = 0
        history: list[dict[str, Any]] = []
        split_array = split.to_numpy()
        monitor_idx = np.flatnonzero(split_array != "test")
        for epoch in range(1, cfg.epochs + 1):
            model.train()
            total_loss = 0.0
            n_seen = 0
            for fb, cb, rfb, rcb, eb, tb, mb, yb in loader:
                fb = fb.to(device)
                cb = cb.to(device)
                rfb = rfb.to(device)
                rcb = rcb.to(device)
                eb = eb.to(device)
                tb = tb.to(device)
                mb = mb.to(device)
                yb = yb.to(device)
                opt.zero_grad(set_to_none=True)
                loss = loss_fn(model(fb, cb, rfb, rcb, eb, tb, mb), yb)
                loss.backward()
                opt.step()
                total_loss += float(loss.detach().cpu()) * len(yb)
                n_seen += len(yb)
            monitor_pred, _ = _predict_cnn(
                model,
                folded[monitor_idx],
                context[monitor_idx],
                raw_folded[monitor_idx],
                raw_context[monitor_idx],
                events[monitor_idx],
                timeline[monitor_idx],
                meta[monitor_idx],
                device=device,
                batch_size=cfg.batch_size,
            )
            row = {"epoch": epoch, "train_loss": total_loss / max(n_seen, 1)}
            monitor_split = split_array[monitor_idx]
            monitor_y = y[monitor_idx]
            for name in ("train", "validation"):
                idx = monitor_split == name
                row[f"{name}_n"] = int(idx.sum())
                row[f"{name}_accuracy"] = float(np.mean(monitor_pred[idx] == monitor_y[idx])) if idx.any() else float("nan")
                row[f"{name}_balanced_accuracy"] = _balanced_accuracy(monitor_y[idx], monitor_pred[idx], len(classes)) if idx.any() else float("nan")
            history.append(row)
            val = row.get("validation_balanced_accuracy", float("nan"))
            print(
                f"[cnn:{profile}] epoch={epoch:03d} loss={row['train_loss']:.5f} "
                f"val_bal={val:.3f}",
                flush=True,
            )
            if np.isfinite(val) and val > best_val:
                best_val = float(val)
                best_epoch = epoch
                best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
                stale = 0
            else:
                stale += 1
            if stale >= cfg.early_stop_patience:
                break
        if best_state is not None:
            model.load_state_dict(best_state)
        pred, prob = _predict_cnn(
            model,
            folded,
            context,
            raw_folded,
            raw_context,
            events,
            timeline,
            meta,
            device=device,
            batch_size=cfg.batch_size,
        )
        test_idx = split_array == "test"
        if history and test_idx.any():
            selected_history = next(
                (item for item in history if int(item["epoch"]) == int(best_epoch)),
                history[-1],
            )
            selected_history["test_n"] = int(test_idx.sum())
            selected_history["test_accuracy"] = float(np.mean(pred[test_idx] == y[test_idx]))
            selected_history["test_balanced_accuracy"] = _balanced_accuracy(
                y[test_idx], pred[test_idx], len(classes)
            )
        predictions = work.copy()
        predictions["cnn_profile"] = profile
        predictions["cnn_label"] = [classes[int(idx)] for idx in pred]
        predictions["cnn_confidence"] = prob.max(axis=1)
        sorted_prob = np.sort(prob, axis=1)
        predictions["cnn_margin"] = sorted_prob[:, -1] - sorted_prob[:, -2] if prob.shape[1] > 1 else 1.0
        predictions["cnn_correct"] = predictions["cnn_label"].eq(predictions["main_teacher_target"])
        for idx, label in enumerate(classes):
            predictions[f"cnn_p_{label}"] = prob[:, idx]
        pred_path = write_table(predictions, profile_dir / "cnn_predictions.parquet")
        pd.DataFrame(history).to_csv(profile_dir / "training_history.csv", index=False)
        torch.save(
            {
                "model_state_dict": model.state_dict(),
                "classes": list(classes),
                "config": asdict(cfg),
                "profile": profile,
                "metadata_columns": metadata_columns if use_metadata else [],
                "metadata_norm": metadata_norm if use_metadata else {"center": [], "scale": []},
                "tensor_npz": str(tensor_npz),
                "tensor_rows": str(tensor_rows),
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
                "test_evaluation_policy": "once_after_validation_epoch_selection",
                "tensor_shapes": {
                    "folded": list(folded.shape),
                    "context": list(context.shape),
                    "raw_folded": list(raw_folded.shape),
                    "raw_context": list(raw_context.shape),
                    "events": list(events.shape),
                    "timeline": list(timeline.shape),
                },
            },
            profile_dir / "cnn_teacher.pt",
        )
        metrics: dict[str, Any] = {}
        for name in ("train", "validation", "test"):
            idx = split_array == name
            metrics[name] = _classification_report_arrays(y[idx], pred[idx], labels=classes)
            metrics[name]["calibration"] = _ece(prob[idx], y[idx])
        for source in ("real_candidate", "injection_recovery"):
            idx = work["source_kind"].fillna("").astype(str).eq(source).to_numpy()
            metrics[f"source_{source}"] = _classification_report_arrays(y[idx], pred[idx], labels=classes)
            metrics[f"source_{source}"]["calibration"] = _ece(prob[idx], y[idx])
        summary = {
            "status": "ok",
            "profile": profile,
            "use_metadata": bool(use_metadata),
            "best_epoch": int(best_epoch),
            "best_validation_balanced_accuracy": best_val if np.isfinite(best_val) else None,
            "classes": list(classes),
            "n_rows": int(len(work)),
            "n_metadata_features": int(meta.shape[1]),
            "metadata_columns": metadata_columns if use_metadata else [],
            "test_evaluation_policy": "once_after_validation_epoch_selection",
            "tensor_shapes": {
                "folded": list(folded.shape),
                "context": list(context.shape),
                "raw_folded": list(raw_folded.shape),
                "raw_context": list(raw_context.shape),
                "events": list(events.shape),
                "timeline": list(timeline.shape),
            },
            "metrics": metrics,
            "outputs": {
                "model": str(profile_dir / "cnn_teacher.pt"),
                "predictions": str(pred_path),
                "history": str(profile_dir / "training_history.csv"),
                "summary": str(profile_dir / "summary.json"),
            },
        }
        (profile_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
        return summary

    profiles = {
        "cnn_shape_only": run_profile("cnn_shape_only", use_metadata=False),
        "cnn_shape_plus_bls": run_profile("cnn_shape_plus_bls", use_metadata=True),
    }
    best_profile = max(
        profiles,
        key=lambda name: profiles[name].get("best_validation_balanced_accuracy") or -np.inf,
    )
    top_summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tensor_npz": str(tensor_npz),
        "tensor_rows": str(tensor_rows),
        "training_table": str(training_table),
        "out_dir": str(out_dir),
        "classes": list(classes),
        "class_counts": {str(k): int(v) for k, v in work["main_teacher_target"].value_counts().sort_index().items()},
        "split_counts": {str(k): int(v) for k, v in work["cnn_training_split"].value_counts().sort_index().items()},
        "split_class_counts": {
            str(cls): {str(split_name): int(value) for split_name, value in row.items()}
            for cls, row in split_class_counts.to_dict(orient="index").items()
        },
        "group_split": "tic",
        "device": str(device),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "torch_device_count": int(torch.cuda.device_count()) if torch.cuda.is_available() else 0,
        "torch_version": str(torch.__version__),
        "torch_cuda_built": str(torch.version.cuda),
        "config": asdict(cfg),
        "profiles": profiles,
        "best_profile_by_validation_balanced_accuracy": best_profile,
        "observable_channel_policy": "only candidate-observable tensor channels; truth/recovery/source/selection columns excluded from inputs",
        "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
    }
    (out_dir / "summary.json").write_text(json.dumps(top_summary, indent=2, sort_keys=True, default=json_default) + "\n")
    (out_dir / "summary.md").write_text(render_cnn_training_markdown(top_summary))
    return top_summary


def _predict_cnn(
    model: Any,
    folded: np.ndarray,
    context: np.ndarray,
    raw_folded: np.ndarray,
    raw_context: np.ndarray,
    events: np.ndarray,
    timeline: np.ndarray,
    metadata: np.ndarray,
    *,
    device: Any,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    import torch

    model.eval()
    pred_parts: list[np.ndarray] = []
    prob_parts: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(folded), batch_size):
            end = min(start + batch_size, len(folded))
            logits = model(
                torch.from_numpy(folded[start:end]).float().to(device),
                torch.from_numpy(context[start:end]).float().to(device),
                torch.from_numpy(raw_folded[start:end]).float().to(device),
                torch.from_numpy(raw_context[start:end]).float().to(device),
                torch.from_numpy(events[start:end]).float().to(device),
                torch.from_numpy(timeline[start:end]).float().to(device),
                torch.from_numpy(metadata[start:end]).float().to(device),
            )
            prob = torch.softmax(logits, dim=1).detach().cpu().numpy()
            prob_parts.append(prob)
            pred_parts.append(np.argmax(prob, axis=1))
    return np.concatenate(pred_parts), np.concatenate(prob_parts)


def score_recovery50_cnn_teacher(
    *,
    tensor_npz: Path,
    tensor_rows: Path,
    model_path: Path,
    feature_table: Path | None,
    out_dir: Path,
    batch_size: int = 256,
    require_cuda: bool = False,
) -> dict[str, Any]:
    """Score tensor rows with a saved CNN teacher profile.

    This intentionally does not require labels in ``tensor_rows``. It is used
    by active-learning queues where the output probabilities are triage
    priorities, not final astrophysical labels.
    """

    import torch

    out_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = torch.load(model_path, map_location="cpu")
    if checkpoint.get("adp_only_contract_version") != ADP_ONLY_CONTRACT_VERSION:
        raise ValueError("refusing to score with a pre-ADP-only CNN checkpoint")
    classes = [str(label) for label in checkpoint["classes"]]
    cfg = dict(checkpoint.get("config", {}))
    tensors = _load_npz(tensor_npz)
    rows = read_table(tensor_rows).copy().reset_index(drop=True)
    assert_adp_only_tensor_rows(rows)
    if rows.empty:
        raise ValueError("no tensor rows to score")
    feature = read_table(feature_table) if feature_table is not None and Path(feature_table).exists() else rows.copy()

    tensor_idx = rows["tensor_index"].to_numpy(dtype=int)
    folded = _prepare_value_mask(tensors["folded"][tensor_idx], tensors["folded_mask"][tensor_idx])
    context = _prepare_value_mask(tensors["context"][tensor_idx], tensors["context_mask"][tensor_idx])
    raw_folded = _prepare_value_mask(tensors["raw_folded"][tensor_idx], tensors["raw_folded_mask"][tensor_idx])
    raw_context = _prepare_value_mask(tensors["raw_context"][tensor_idx], tensors["raw_context_mask"][tensor_idx])
    events = _prepare_value_mask(tensors["events"][tensor_idx], tensors["events_mask"][tensor_idx])
    timeline = _prepare_value_mask(tensors["timeline"][tensor_idx], tensors["timeline_mask"][tensor_idx])
    metadata_columns = [str(col) for col in checkpoint.get("metadata_columns", [])]
    bad_metadata = canonical_det_flux_columns(metadata_columns)
    if bad_metadata:
        raise ValueError(f"CNN checkpoint contains canonical DET_FLUX metadata: {bad_metadata}")
    metadata = _metadata_matrix_from_saved_norm(
        rows=rows,
        feature_table=feature,
        columns=metadata_columns,
        norm=checkpoint.get("metadata_norm", {"center": [], "scale": []}),
    )

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required but torch.cuda.is_available() is false")
    model = _build_recovery50_cnn_model(
        folded_channels=folded.shape[1],
        raw_folded_channels=raw_folded.shape[1],
        raw_context_channels=raw_context.shape[1],
        event_channels=events.shape[1],
        timeline_channels=timeline.shape[1],
        metadata_dim=metadata.shape[1],
        n_classes=len(classes),
        dropout=float(cfg.get("dropout", 0.0)),
    ).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    pred, prob = _predict_cnn(
        model,
        folded,
        context,
        raw_folded,
        raw_context,
        events,
        timeline,
        metadata,
        device=device,
        batch_size=batch_size,
    )

    predictions = rows.copy()
    predictions["cnn_profile"] = str(checkpoint.get("profile", model_path.parent.name))
    predictions["cnn_label"] = [classes[int(idx)] for idx in pred]
    predictions["cnn_confidence"] = prob.max(axis=1)
    sorted_prob = np.sort(prob, axis=1)
    predictions["cnn_margin"] = sorted_prob[:, -1] - sorted_prob[:, -2] if prob.shape[1] > 1 else 1.0
    for idx, label in enumerate(classes):
        predictions[f"cnn_p_{label}"] = prob[:, idx]
    pred_path = write_table(predictions, out_dir / "cnn_scored_candidates.parquet")
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tensor_npz": str(tensor_npz),
        "tensor_rows": str(tensor_rows),
        "model_path": str(model_path),
        "feature_table": str(feature_table) if feature_table else "",
        "out_dir": str(out_dir),
        "classes": classes,
        "profile": str(checkpoint.get("profile", model_path.parent.name)),
        "n_rows": int(len(predictions)),
        "n_metadata_features": int(metadata.shape[1]),
        "metadata_columns": metadata_columns,
        "device": str(device),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "outputs": {"predictions": str(pred_path), "summary": str(out_dir / "summary.json")},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    return summary


def render_cnn_training_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 Recovery50 CNN Teacher",
        "",
        f"- classes: `{', '.join(summary['classes'])}`",
        f"- class counts: `{summary['class_counts']}`",
        f"- split counts: `{summary['split_counts']}`",
        f"- group split: `{summary['group_split']}`",
        f"- device: `{summary['device']}`",
        f"- best profile: `{summary['best_profile_by_validation_balanced_accuracy']}`",
        "",
        "| Profile | Metadata | Best Val Balanced Accuracy | Test Balanced Accuracy | Real Balanced Accuracy | Injected Balanced Accuracy |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, payload in summary["profiles"].items():
        metrics = payload["metrics"]
        lines.append(
            f"| `{name}` | `{payload['use_metadata']}` | "
            f"{(payload.get('best_validation_balanced_accuracy') or float('nan')):.3f} | "
            f"{metrics['test']['balanced_accuracy']:.3f} | "
            f"{metrics['source_real_candidate']['balanced_accuracy']:.3f} | "
            f"{metrics['source_injection_recovery']['balanced_accuracy']:.3f} |"
        )
    lines.append("")
    return "\n".join(lines)
