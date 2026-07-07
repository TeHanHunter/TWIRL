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
from twirl.vetting.recovery50_teacher import (
    DEFAULT_APERTURES,
    json_default,
    leakage_columns,
    load_feature_table,
    metadata_feature_columns,
    read_table,
    select_teacher_classes,
    write_table,
    _confusion_metrics,
    _read_lc_for_row,
    _safe_float,
)


DEFAULT_FOLDED_POINTS = 512
DEFAULT_CONTEXT_POINTS = 512
DEFAULT_EVENT_POINTS = 128
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
    for idx in range(len(centers)):
        values = y[bin_id == idx]
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
    event_values: list[np.ndarray] = []
    event_masks: list[np.ndarray] = []
    event_counts: list[np.ndarray] = []
    event_centers: list[np.ndarray] = []
    records: list[dict[str, Any]] = []
    status_counts: dict[str, int] = {}

    for processed, (_, row) in enumerate(rows.iterrows(), start=1):
        period_d = _safe_float(row.get("period_d"))
        t0_bjd = _safe_float(row.get("t0_bjd"))
        duration_min = _safe_float(row.get("duration_min"))
        rec = {
            "tensor_index": len(records),
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
        missing = [ap for ap in config.apertures if ap not in lc.flux]
        if missing:
            rec["tensor_status"] = "missing_aperture:" + ",".join(missing)
            status_counts[rec["tensor_status"]] = status_counts.get(rec["tensor_status"], 0) + 1
            continue

        time_bjd = _time_to_bjd(np.asarray(lc.time, dtype=np.float64))
        quality = np.asarray(lc.quality, dtype=np.int32)
        flux_by_aperture = [np.asarray(lc.flux[ap], dtype=np.float64) for ap in config.apertures]

        folded, folded_mask, folded_count, folded_x = folded_view(
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            period_d=period_d,
            t0_bjd=t0_bjd,
            duration_min=duration_min,
            n_points=config.folded_points,
            window_durations=config.folded_window_durations,
        )
        context, context_mask, context_count, context_x = folded_view(
            time_bjd=time_bjd,
            flux_by_aperture=flux_by_aperture,
            quality=quality,
            period_d=period_d,
            t0_bjd=t0_bjd,
            duration_min=duration_min,
            n_points=config.context_points,
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
        status = _tensor_status(folded_mask)
        rec["tensor_status"] = status
        rec["n_event_windows"] = int(np.isfinite(centers).sum())
        rec["folded_observed_fraction"] = float(folded_mask.sum() / folded_mask.size) if folded_mask.size else 0.0
        rec["context_observed_fraction"] = float(context_mask.sum() / context_mask.size) if context_mask.size else 0.0
        rec["event_observed_fraction"] = float(event_mask.sum() / event_mask.size) if event_mask.size else 0.0
        if status != "ok":
            status_counts[status] = status_counts.get(status, 0) + 1
            continue

        folded_values.append(folded)
        folded_masks.append(folded_mask)
        folded_counts.append(folded_count)
        context_values.append(context)
        context_masks.append(context_mask)
        context_counts.append(context_count)
        event_values.append(events)
        event_masks.append(event_mask)
        event_counts.append(event_count)
        event_centers.append(centers)
        records.append(rec)
        status_counts[status] = status_counts.get(status, 0) + 1
        if progress_every > 0 and (processed % progress_every == 0 or processed == len(rows)):
            print(f"[cnn-tensor] processed={processed:,}/{len(rows):,}; accepted={len(records):,}", flush=True)

    n_channels = len(config.apertures) + (1 if len(config.apertures) >= 2 else 0)
    folded_arr = (
        np.stack(folded_values, axis=0).astype(np.float32)
        if folded_values
        else np.empty((0, n_channels, config.folded_points), dtype=np.float32)
    )
    folded_mask_arr = (
        np.stack(folded_masks, axis=0).astype(np.bool_)
        if folded_masks
        else np.empty((0, n_channels, config.folded_points), dtype=np.bool_)
    )
    folded_count_arr = (
        np.stack(folded_counts, axis=0).astype(np.int16)
        if folded_counts
        else np.empty((0, n_channels, config.folded_points), dtype=np.int16)
    )
    context_arr = (
        np.stack(context_values, axis=0).astype(np.float32)
        if context_values
        else np.empty((0, n_channels, config.context_points), dtype=np.float32)
    )
    context_mask_arr = (
        np.stack(context_masks, axis=0).astype(np.bool_)
        if context_masks
        else np.empty((0, n_channels, config.context_points), dtype=np.bool_)
    )
    context_count_arr = (
        np.stack(context_counts, axis=0).astype(np.int16)
        if context_counts
        else np.empty((0, n_channels, config.context_points), dtype=np.int16)
    )
    event_arr = (
        np.stack(event_values, axis=0).astype(np.float32)
        if event_values
        else np.empty((0, n_channels, config.max_events, config.event_points), dtype=np.float32)
    )
    event_mask_arr = (
        np.stack(event_masks, axis=0).astype(np.bool_)
        if event_masks
        else np.empty((0, n_channels, config.max_events, config.event_points), dtype=np.bool_)
    )
    event_count_arr = (
        np.stack(event_counts, axis=0).astype(np.int16)
        if event_counts
        else np.empty((0, n_channels, config.max_events, config.event_points), dtype=np.int16)
    )
    event_centers_arr = (
        np.stack(event_centers, axis=0).astype(np.float64)
        if event_centers
        else np.empty((0, config.max_events), dtype=np.float64)
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
        events=event_arr,
        events_mask=event_mask_arr,
        events_counts=event_count_arr,
        event_centers_bjd=event_centers_arr,
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
        apertures=np.asarray(config.apertures, dtype=object),
        derived_channels=np.asarray(["primary_minus_small"] if len(config.apertures) >= 2 else [], dtype=object),
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
        "events_shape": list(event_arr.shape),
        "status_counts": {str(k): int(v) for k, v in sorted(status_counts.items())},
        "apertures": list(config.apertures),
        "tensor_config": asdict(config),
        "observable_channel_policy": "only channels present for both real and injected rows; injected original/pre-injection arrays excluded",
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
        f"- events shape: `{tuple(summary['events_shape'])}`",
        f"- apertures: `{', '.join(summary['apertures'])}`",
        f"- status counts: `{summary['status_counts']}`",
        f"- channel policy: `{summary['observable_channel_policy']}`",
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
    join_cols = [col for col in ("row_id", "review_id") if col in rows.columns and col in feature_table.columns]
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
    training = load_feature_table(training_table=training_table, metrics_tables=tuple(metrics_tables))
    classes = select_teacher_classes(training, cfg.min_class_count)
    if len(classes) < 2:
        raise ValueError(f"need at least two teacher classes, got {classes}")

    work = rows[rows["main_teacher_target"].isin(classes)].copy().reset_index(drop=True)
    if work.empty:
        raise ValueError("no tensor rows match selected teacher classes")
    split = _grouped_splits(work, classes=classes, cfg=cfg)
    work["cnn_training_split"] = split
    train_mask = split.to_numpy() == "train"
    if len(np.unique(work.loc[train_mask, "main_teacher_target"])) < 2:
        raise ValueError("CNN train split has fewer than two classes")

    tensor_idx = work["tensor_index"].to_numpy(dtype=int)
    folded = _prepare_value_mask(tensors["folded"][tensor_idx], tensors["folded_mask"][tensor_idx])
    context = _prepare_value_mask(tensors["context"][tensor_idx], tensors["context_mask"][tensor_idx])
    events = _prepare_value_mask(tensors["events"][tensor_idx], tensors["events_mask"][tensor_idx])
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
        def __init__(self, folded_channels: int, event_channels: int, metadata_dim: int, n_classes: int):
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
            self.event_stem = nn.Sequential(
                nn.Conv2d(event_channels, 16, kernel_size=(3, 9), padding=(1, 4)),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.Conv2d(16, 16, kernel_size=(3, 7), padding=(1, 3)),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.AdaptiveAvgPool2d((1, 1)),
            )
            hidden_in = 32 + 16 + 16 + metadata_dim
            self.head = nn.Sequential(
                nn.Dropout(cfg.dropout),
                nn.Linear(hidden_in, 64),
                nn.ReLU(),
                nn.Dropout(cfg.dropout),
                nn.Linear(64, n_classes),
            )

        def forward(self, folded_x: Any, context_x: Any, event_x: Any, metadata: Any) -> Any:
            f = self.folded_stem(folded_x).flatten(1)
            c = self.context_stem(context_x).flatten(1)
            e = self.event_stem(event_x).flatten(1)
            parts = [f, c, e]
            if metadata.shape[1] > 0:
                parts.append(metadata)
            return self.head(torch.cat(parts, dim=1))

    def run_profile(profile: str, use_metadata: bool) -> dict[str, Any]:
        profile_dir = out_dir / profile
        profile_dir.mkdir(parents=True, exist_ok=True)
        meta = metadata_x if use_metadata else np.empty((len(work), 0), dtype=np.float32)
        model = Recovery50CNN(
            folded_channels=folded.shape[1],
            event_channels=events.shape[1],
            metadata_dim=meta.shape[1],
            n_classes=len(classes),
        ).to(device)
        opt = torch.optim.AdamW(model.parameters(), lr=cfg.learning_rate, weight_decay=cfg.weight_decay)
        loss_fn = nn.CrossEntropyLoss(weight=torch.tensor(class_weights, dtype=torch.float32, device=device))
        ds = TensorDataset(
            torch.from_numpy(folded[train_idx]).float(),
            torch.from_numpy(context[train_idx]).float(),
            torch.from_numpy(events[train_idx]).float(),
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
        for epoch in range(1, cfg.epochs + 1):
            model.train()
            total_loss = 0.0
            n_seen = 0
            for fb, cb, eb, mb, yb in loader:
                fb = fb.to(device)
                cb = cb.to(device)
                eb = eb.to(device)
                mb = mb.to(device)
                yb = yb.to(device)
                opt.zero_grad(set_to_none=True)
                loss = loss_fn(model(fb, cb, eb, mb), yb)
                loss.backward()
                opt.step()
                total_loss += float(loss.detach().cpu()) * len(yb)
                n_seen += len(yb)
            pred, prob = _predict_cnn(model, folded, context, events, meta, device=device, batch_size=cfg.batch_size)
            row = {"epoch": epoch, "train_loss": total_loss / max(n_seen, 1)}
            for name in ("train", "validation", "test"):
                idx = split_array == name
                row[f"{name}_n"] = int(idx.sum())
                row[f"{name}_accuracy"] = float(np.mean(pred[idx] == y[idx])) if idx.any() else float("nan")
                row[f"{name}_balanced_accuracy"] = _balanced_accuracy(y[idx], pred[idx], len(classes)) if idx.any() else float("nan")
            history.append(row)
            val = row.get("validation_balanced_accuracy", float("nan"))
            print(
                f"[cnn:{profile}] epoch={epoch:03d} loss={row['train_loss']:.5f} "
                f"val_bal={val:.3f} test_bal={row.get('test_balanced_accuracy', float('nan')):.3f}",
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
        pred, prob = _predict_cnn(model, folded, context, events, meta, device=device, batch_size=cfg.batch_size)
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
    }
    (out_dir / "summary.json").write_text(json.dumps(top_summary, indent=2, sort_keys=True, default=json_default) + "\n")
    (out_dir / "summary.md").write_text(render_cnn_training_markdown(top_summary))
    return top_summary


def _predict_cnn(
    model: Any,
    folded: np.ndarray,
    context: np.ndarray,
    events: np.ndarray,
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
                torch.from_numpy(events[start:end]).float().to(device),
                torch.from_numpy(metadata[start:end]).float().to(device),
            )
            prob = torch.softmax(logits, dim=1).detach().cpu().numpy()
            prob_parts.append(prob)
            pred_parts.append(np.argmax(prob, axis=1))
    return np.concatenate(pred_parts), np.concatenate(prob_parts)


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
