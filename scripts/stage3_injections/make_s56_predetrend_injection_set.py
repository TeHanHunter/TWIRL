#!/usr/bin/env python3
"""Build pre-detrend LC-level injections from raw TGLC HDF5 products.

This is the scientifically useful LC-level path for testing whether TWIRL-FS
detrending attenuates or reshapes short WD occultations:

1. read raw per-orbit TGLC HDF5 light curves,
2. inject a finite-exposure ``batman`` transit into each raw aperture flux,
3. rerun canonical and optional ADP TWIRL-FS detrending,
4. write the same injection HDF5 schema consumed by the Stage 5 review builder.

This is not a pixel-level injection. It starts after TGLC extraction, so it
does not measure extraction, deblending, aperture choice, or ePSF-fit failures.
It does measure the specific detrending failure mode that matters for the
current S56 vetting decision.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import importlib.util
import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import BJDREFI  # noqa: E402
from twirl.lightcurves.tglc_h5_reader import (  # noqa: E402
    APERTURE_KEYS,
    TGLCAperture,
    TGLCLightCurve,
    read_tglc_h5,
)
from twirl.search.injections import (  # noqa: E402
    batman_transit_model,
    box_transit_mask,
    choose_observed_epoch,
    inject_batman_transit,
)


DEFAULT_ORBIT_ROOTS = (
    Path("/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi"),
    Path("/pdo/users/tehan/tglc-gpu-production/orbit-120/ffi"),
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
    / "s56_predetrend_batman_smoke"
)

DEFAULT_OUTPUT_APERTURES = (
    "DET_FLUX_ADP",
    "DET_FLUX",
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP_LAG",
    "DET_FLUX_SML",
    "DET_FLUX_LAG",
)

DEFAULT_CADENCE_S = 200.0
TESS_FLUX_E_PER_S_AT_T10 = 1.5e4
BASELINE_SOURCES = ("raw_median", "inferred_aperture", "tessmag_total")

OUTPUT_TO_RAW_APERTURE = {
    "DET_FLUX": "Primary",
    "DET_FLUX_SML": "Small",
    "DET_FLUX_LAG": "Large",
    "DET_FLUX_ADP": "Primary",
    "DET_FLUX_ADP_SML": "Small",
    "DET_FLUX_ADP_LAG": "Large",
}

CANONICAL_OUTPUTS = {
    "DET_FLUX": "Primary",
    "DET_FLUX_SML": "Small",
    "DET_FLUX_LAG": "Large",
}
ADP_OUTPUTS = {
    "DET_FLUX_ADP": "Primary",
    "DET_FLUX_ADP_SML": "Small",
    "DET_FLUX_ADP_LAG": "Large",
}

PARAM_ATTR_KEYS = (
    "radius_rearth",
    "radius_rwd",
    "impact_b",
    "a_over_rwd",
    "inclination_deg",
    "target_depth",
    "geometric_depth",
    "model_depth",
    "sampled_model_depth",
    "duration_model",
    "injection_model",
    "sampling_mode",
    "grid_period_bin",
    "grid_radius_bin",
    "grid_depth_bin",
    "grid_cell_id",
    "period_bin_lo_d",
    "period_bin_hi_d",
    "radius_bin_lo_rearth",
    "radius_bin_hi_rearth",
    "depth_bin_lo",
    "depth_bin_hi",
)

MANIFEST_COLUMNS = (
    "injection_id",
    "tic",
    "sector",
    "camera",
    "ccd",
    "tessmag",
    "aperture",
    "apertures",
    "signal_family",
    "period_d",
    "t0_d",
    "t0_bjd",
    "duration_min",
    "depth",
    *PARAM_ATTR_KEYS,
    "n_good_in_transit",
    "baseline",
    "baseline_source",
    "expected_tess_flux_per_cadence",
    "split",
    "label",
    "label_source",
    "source_target_paths",
    "h5_group",
)


def _load_script_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_POSTDET = _load_script_module(
    "make_s56_lc_injection_training_set",
    REPO_ROOT / "scripts/stage3_injections/make_s56_lc_injection_training_set.py",
)
_TWIRLFS = _load_script_module(
    "build_twirl_hlsp",
    REPO_ROOT / "scripts/stage1_lightcurves/build_twirl_hlsp.py",
)


def _dataset(group, name: str, values, compression: str | None) -> None:
    kwargs = {"shuffle": True}
    if compression:
        kwargs["compression"] = compression
    group.create_dataset(name, data=values, **kwargs)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _scalar(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    return value


def _read_tic_filter(path: Path | None, tic_column: str) -> set[int] | None:
    if path is None:
        return None
    import pandas as pd

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        table = pd.read_parquet(path, columns=[tic_column])
    elif suffix == ".csv":
        table = pd.read_csv(path, usecols=[tic_column])
    elif suffix in {".json", ".jsonl"}:
        table = pd.read_json(path, lines=suffix == ".jsonl")
    else:
        raise ValueError(f"unsupported TIC table format: {path}")
    if tic_column not in table.columns:
        raise KeyError(f"TIC column {tic_column!r} not found in {path}")
    return {int(tic) for tic in table[tic_column].dropna().astype(int).tolist()}


def _read_target_tmag_table(path: Path | None, tic_column: str, tmag_column: str) -> dict[int, float] | None:
    if path is None:
        return None
    import pandas as pd

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        table = pd.read_parquet(path, columns=[tic_column, tmag_column])
    elif suffix == ".csv":
        table = pd.read_csv(path, usecols=[tic_column, tmag_column])
    elif suffix in {".json", ".jsonl"}:
        table = pd.read_json(path, lines=suffix == ".jsonl")
    else:
        raise ValueError(f"unsupported target Tmag table format: {path}")
    missing = [col for col in (tic_column, tmag_column) if col not in table.columns]
    if missing:
        raise KeyError(f"missing column(s) in {path}: {missing}")
    table = table[[tic_column, tmag_column]].dropna().drop_duplicates(tic_column)
    return {
        int(row[tic_column]): float(row[tmag_column])
        for _, row in table.iterrows()
    }


def _read_path_list(path: Path | None) -> tuple[Path, ...]:
    if path is None:
        return ()
    rows: list[Path] = []
    for raw in path.read_text().splitlines():
        value = raw.strip()
        if not value or value.startswith("#"):
            continue
        rows.append(Path(value))
    return tuple(rows)


def _parse_float_list(raw: str | None) -> tuple[float, ...] | None:
    if raw is None:
        return None
    values: list[float] = []
    for piece in str(raw).split(","):
        text = piece.strip().lower()
        if not text:
            continue
        if text in {"inf", "+inf", "infinity", "+infinity"}:
            values.append(float("inf"))
        elif text in {"-inf", "-infinity"}:
            values.append(float("-inf"))
        else:
            values.append(float(text))
    return tuple(values)


def _parse_tmag_bin_edges(raw: str | None) -> tuple[float, ...] | None:
    edges = _parse_float_list(raw)
    if edges is None:
        return None
    if len(edges) < 2:
        raise argparse.ArgumentTypeError("--target-tmag-bin-edges needs at least two comma-separated edges")
    for lo, hi in zip(edges[:-1], edges[1:]):
        if not hi > lo:
            raise argparse.ArgumentTypeError("--target-tmag-bin-edges must be strictly increasing")
    return edges


def _target_tmag_bin_label(lo: float, hi: float) -> str:
    if np.isneginf(lo):
        return f"Tmag < {hi:g}"
    if np.isposinf(hi):
        return f"Tmag > {lo:g}"
    return f"{lo:g} <= Tmag < {hi:g}"


def _read_target_tmag_map(target_paths: dict[int, list[Path]]) -> dict[int, float]:
    import h5py

    out: dict[int, float] = {}
    for tic, paths in target_paths.items():
        if not paths:
            continue
        try:
            with h5py.File(paths[0], "r") as h5:
                out[int(tic)] = float(h5.attrs.get("TessMag", np.nan))
        except Exception:
            out[int(tic)] = float("nan")
    return out


def _build_target_tmag_sampler(
    target_tmag_by_tic: dict[int, float],
    *,
    edges: tuple[float, ...] | None,
    weights: tuple[float, ...] | None,
) -> dict[str, Any] | None:
    if edges is None:
        return None
    n_bins = len(edges) - 1
    if weights is None:
        raw_weights = np.ones(n_bins, dtype=float)
    else:
        if len(weights) != n_bins:
            raise ValueError(
                f"target tmag weights length ({len(weights)}) must match "
                f"number of bins ({n_bins})"
            )
        raw_weights = np.asarray(weights, dtype=float)
        if np.any(~np.isfinite(raw_weights)) or np.any(raw_weights < 0):
            raise ValueError("target tmag weights must be finite non-negative values")
        if not np.any(raw_weights > 0):
            raise ValueError("at least one target tmag weight must be positive")

    bins: list[np.ndarray] = []
    labels: list[str] = []
    for idx, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        tics = sorted(
            tic
            for tic, tmag in target_tmag_by_tic.items()
            if np.isfinite(tmag) and tmag >= lo and tmag < hi
        )
        bins.append(np.asarray(tics, dtype=np.int64))
        labels.append(_target_tmag_bin_label(float(lo), float(hi)))
        if len(tics) == 0:
            raw_weights[idx] = 0.0
    if not np.any(raw_weights > 0):
        raise ValueError("no targets fall inside the requested target Tmag bins")
    probabilities = raw_weights / raw_weights.sum()
    return {
        "edges": tuple(float(edge) for edge in edges),
        "labels": tuple(labels),
        "bins": tuple(bins),
        "probabilities": probabilities,
        "n_targets_by_bin": tuple(int(len(tics)) for tics in bins),
    }


def _draw_target_tic(
    rng: np.random.Generator,
    target_tics: np.ndarray,
    target_tmag_sampler: dict[str, Any] | None,
) -> int:
    if target_tmag_sampler is None:
        return int(rng.choice(target_tics))
    probabilities = np.asarray(target_tmag_sampler["probabilities"], dtype=float)
    bin_index = int(rng.choice(np.arange(len(probabilities)), p=probabilities))
    tics = np.asarray(target_tmag_sampler["bins"][bin_index], dtype=np.int64)
    if tics.size == 0:
        return int(rng.choice(target_tics))
    return int(rng.choice(tics))


def _discover_target_paths(
    orbit_roots: tuple[Path, ...],
    *,
    target_h5_paths: tuple[Path, ...] | None,
    tic_filter: set[int] | None,
    limit_targets: int | None,
) -> dict[int, list[Path]]:
    if target_h5_paths:
        out: dict[int, list[Path]] = {}
        for path in target_h5_paths:
            try:
                tic = int(Path(path).stem)
            except ValueError:
                continue
            if tic_filter is not None and tic not in tic_filter:
                continue
            out.setdefault(tic, []).append(Path(path))
        if limit_targets is not None:
            keep = set(sorted(out)[: int(limit_targets)])
            out = {tic: paths for tic, paths in out.items() if tic in keep}
        return out

    per_orbit_index = [_TWIRLFS.discover_orbit_hdf5(root) for root in orbit_roots]
    all_tics = sorted(set().union(*(idx.keys() for idx in per_orbit_index)))
    if tic_filter is not None:
        all_tics = [tic for tic in all_tics if tic in tic_filter]
    if limit_targets is not None:
        all_tics = all_tics[: int(limit_targets)]
    out: dict[int, list[Path]] = {}
    for tic in all_tics:
        paths = [idx[tic] for idx in per_orbit_index if tic in idx]
        if paths:
            out[int(tic)] = paths
    return out


def _clone_with_raw_flux(
    lc: TGLCLightCurve,
    raw_flux_by_aperture: dict[str, np.ndarray],
) -> TGLCLightCurve:
    apertures: dict[str, TGLCAperture] = {}
    for ap_key in APERTURE_KEYS:
        old = lc.apertures[ap_key]
        raw_flux = np.asarray(raw_flux_by_aperture.get(ap_key, old.raw_flux), dtype=np.float64)
        apertures[ap_key] = TGLCAperture(
            name=old.name,
            raw_flux=raw_flux,
            raw_flux_err=np.asarray(old.raw_flux_err, dtype=np.float64),
            raw_magnitude=np.asarray(old.raw_magnitude, dtype=np.float64),
            raw_magnitude_err=np.asarray(old.raw_magnitude_err, dtype=np.float64),
            centroid_x=np.asarray(old.centroid_x, dtype=np.float64),
            centroid_y=np.asarray(old.centroid_y, dtype=np.float64),
            flux_was_synthesized=old.flux_was_synthesized,
        )
    cloned = TGLCLightCurve(
        tic=lc.tic,
        sector=lc.sector,
        orbit=lc.orbit,
        cam=lc.cam,
        ccd=lc.ccd,
        tmag=lc.tmag,
        ra=lc.ra,
        dec=lc.dec,
        bjd_offset=lc.bjd_offset,
        time=np.asarray(lc.time, dtype=np.float64),
        cadence=np.asarray(lc.cadence, dtype=np.int64),
        quality=np.asarray(lc.quality, dtype=np.int64),
        centroid_x=np.asarray(lc.centroid_x, dtype=np.float64),
        centroid_y=np.asarray(lc.centroid_y, dtype=np.float64),
        background=np.asarray(lc.background, dtype=np.float64),
        background_err=np.asarray(lc.background_err, dtype=np.float64),
        apertures=apertures,
        path=lc.path,
    )
    if hasattr(lc, "orbitid"):
        cloned.orbitid = np.asarray(getattr(lc, "orbitid"), dtype=np.int32)
    return cloned


def _positive_baseline(raw_flux: np.ndarray, quality: np.ndarray) -> float:
    raw_flux = np.asarray(raw_flux, dtype=np.float64)
    good = (np.asarray(quality) == 0) & np.isfinite(raw_flux)
    if not np.any(good):
        return float("nan")
    baseline = float(np.nanmedian(raw_flux[good]))
    if not np.isfinite(baseline) or baseline <= 0:
        return float("nan")
    return baseline


def tess_flux_per_cadence_from_tmag(tmag: float, cadence_s: float = DEFAULT_CADENCE_S) -> float:
    """Expected total TESS target flux per cadence from TESS magnitude."""
    tmag_arr = np.asarray(tmag, dtype=np.float64)
    out = TESS_FLUX_E_PER_S_AT_T10 * 10.0 ** ((10.0 - tmag_arr) / 2.5) * float(cadence_s)
    out = np.where(np.isfinite(tmag_arr), out, np.nan)
    if out.ndim == 0:
        return float(out)
    return out


def _infer_aperture_flux_baseline(
    *,
    raw_flux: np.ndarray,
    raw_magnitude: np.ndarray,
    quality: np.ndarray,
    tmag: float,
    cadence_s: float,
    min_points: int = 30,
) -> tuple[float, float]:
    """Infer target aperture flux and aperture fraction from TGLC raw columns.

    TGLC ``RawMagnitude`` is derived after correcting aperture flux back toward
    total target flux. Comparing positive ``RawFlux`` to the flux implied by
    ``RawMagnitude`` therefore estimates the aperture fraction. Multiplying
    that fraction by the target's catalog ``TessMag`` flux gives a physical
    injection baseline that is less fragile than a background-subtracted median.
    """
    raw_flux = np.asarray(raw_flux, dtype=np.float64)
    raw_magnitude = np.asarray(raw_magnitude, dtype=np.float64)
    quality = np.asarray(quality)
    mask = (
        (quality == 0)
        & np.isfinite(raw_flux)
        & (raw_flux > 0)
        & np.isfinite(raw_magnitude)
    )
    if np.count_nonzero(mask) < int(min_points):
        return float("nan"), float("nan")
    inferred_total = tess_flux_per_cadence_from_tmag(raw_magnitude[mask], cadence_s=cadence_s)
    with np.errstate(divide="ignore", invalid="ignore"):
        aperture_fraction = raw_flux[mask] / inferred_total
    aperture_fraction = aperture_fraction[
        np.isfinite(aperture_fraction) & (aperture_fraction > 0)
    ]
    if aperture_fraction.size < int(min_points):
        return float("nan"), float("nan")
    fraction = float(np.nanmedian(aperture_fraction))
    expected_total = tess_flux_per_cadence_from_tmag(tmag, cadence_s=cadence_s)
    baseline = float(expected_total * fraction) if np.isfinite(expected_total) else float("nan")
    return baseline, fraction


def _baseline_diagnostics(
    lc: TGLCLightCurve,
    *,
    quality: np.ndarray,
    baseline_source: str,
    cadence_s: float,
) -> dict[str, Any]:
    if baseline_source not in BASELINE_SOURCES:
        raise ValueError(f"unknown baseline source: {baseline_source!r}")
    expected_total = tess_flux_per_cadence_from_tmag(lc.tmag, cadence_s=cadence_s)
    raw_baseline_by_ap: dict[str, float] = {}
    inferred_baseline_by_ap: dict[str, float] = {}
    aperture_fraction_by_ap: dict[str, float] = {}
    injection_baseline_by_ap: dict[str, float] = {}

    for ap in APERTURE_KEYS:
        ap_data = lc.apertures[ap]
        raw_baseline = _positive_baseline(ap_data.raw_flux, quality)
        inferred_baseline, aperture_fraction = _infer_aperture_flux_baseline(
            raw_flux=ap_data.raw_flux,
            raw_magnitude=ap_data.raw_magnitude,
            quality=quality,
            tmag=lc.tmag,
            cadence_s=cadence_s,
        )
        raw_baseline_by_ap[ap] = raw_baseline
        inferred_baseline_by_ap[ap] = inferred_baseline
        aperture_fraction_by_ap[ap] = aperture_fraction

    primary_raw = raw_baseline_by_ap["Primary"]
    primary_inferred = inferred_baseline_by_ap["Primary"]
    for ap in APERTURE_KEYS:
        raw_baseline = raw_baseline_by_ap[ap]
        inferred_baseline = inferred_baseline_by_ap[ap]
        if baseline_source == "raw_median":
            baseline = raw_baseline
        elif baseline_source == "tessmag_total":
            baseline = expected_total
        else:
            baseline = inferred_baseline
        if not np.isfinite(baseline) or baseline <= 0:
            if np.isfinite(primary_inferred) and primary_inferred > 0:
                baseline = primary_inferred
            elif np.isfinite(primary_raw) and primary_raw > 0:
                baseline = primary_raw
        injection_baseline_by_ap[ap] = baseline

    return {
        "expected_tess_flux_per_cadence": expected_total,
        "raw_baseline_by_ap": raw_baseline_by_ap,
        "inferred_baseline_by_ap": inferred_baseline_by_ap,
        "aperture_fraction_by_ap": aperture_fraction_by_ap,
        "injection_baseline_by_ap": injection_baseline_by_ap,
    }


def _detrended_output_columns(
    lc: TGLCLightCurve,
    *,
    adaptive_bkspace: float,
    adaptive_gap_split: float,
) -> tuple[dict[str, np.ndarray], dict[str, dict[str, object]]]:
    columns: dict[str, np.ndarray] = {}
    diagnostics: dict[str, dict[str, object]] = {}
    canonical, canonical_diag = _TWIRLFS.detrend_apertures(
        lc,
        cfg=_TWIRLFS.CANONICAL_DETREND_CONFIG,
    )
    for out_name, ap_key in CANONICAL_OUTPUTS.items():
        columns[out_name] = np.asarray(canonical[ap_key][1], dtype=np.float32)
        diagnostics[out_name] = dict(canonical_diag.get(ap_key, {}))

    adp_cfg = _TWIRLFS.adaptive_detrend_config(adaptive_bkspace, adaptive_gap_split)
    adaptive, adaptive_diag = _TWIRLFS.detrend_apertures(lc, cfg=adp_cfg)
    for out_name, ap_key in ADP_OUTPUTS.items():
        columns[out_name] = np.asarray(adaptive[ap_key][1], dtype=np.float32)
        diagnostics[out_name] = dict(adaptive_diag.get(ap_key, {}))
    return columns, diagnostics


def _draw_params(
    injection_index: int,
    *,
    rng: np.random.Generator,
    sampling_mode: str,
    families: tuple[str, ...],
    grid_period_range_d: tuple[float, float],
    grid_radius_range_rearth: tuple[float, float],
    grid_depth_range: tuple[float, float],
    grid_period_bins: int,
    grid_radius_bins: int,
    grid_depth_bins: int,
    grid_depth_spacing: str,
) -> dict[str, Any]:
    if sampling_mode == "period_depth_grid":
        return _POSTDET._draw_depth_grid_params(
            injection_index,
            rng=rng,
            period_range=grid_period_range_d,
            depth_range=grid_depth_range,
            radius_range_rearth=grid_radius_range_rearth,
            period_bins=grid_period_bins,
            depth_bins=grid_depth_bins,
            depth_spacing=grid_depth_spacing,
        )
    if sampling_mode == "period_radius_grid":
        return _POSTDET._draw_grid_params(
            injection_index,
            rng=rng,
            period_range=grid_period_range_d,
            radius_range_rearth=grid_radius_range_rearth,
            period_bins=grid_period_bins,
            radius_bins=grid_radius_bins,
        )
    family = families[injection_index % len(families)]
    params = _POSTDET._draw_params(family, rng)
    params["sampling_mode"] = "family_cycle"
    if params.get("radius_rwd", "") == "":
        radius_rwd = float(np.sqrt(max(float(params["depth"]), 1.0e-6)))
        period_d = float(params["period_d"])
        a_over_rwd = _POSTDET._a_over_rwd_from_period(period_d)
        impact_b = float(rng.uniform(0.0, min(0.8 * (1.0 + radius_rwd), 0.8 * a_over_rwd)))
        params["radius_rwd"] = radius_rwd
        params["radius_rearth"] = radius_rwd * _POSTDET.WD_RADIUS_REARTH
        params["impact_b"] = impact_b
        params["a_over_rwd"] = a_over_rwd
        params["inclination_deg"] = float(
            np.degrees(np.arccos(np.clip(impact_b / a_over_rwd, 0.0, 1.0)))
        )
        params["geometric_depth"] = float(_POSTDET._circle_overlap_depth(radius_rwd, impact_b))
        params["target_depth"] = float(params["depth"])
        params["duration_model"] = "empirical_batman"
        params["injection_model"] = "batman_quadratic"
    return params


def _split_for_tic(tic: int) -> str:
    return _POSTDET._split_for_tic(tic)


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _scalar(row.get(key, "")) for key in fieldnames})


def make_predetrend_injections(
    *,
    orbit_roots: tuple[Path, ...],
    out_dir: Path,
    apertures: tuple[str, ...],
    n_injections: int,
    families: tuple[str, ...],
    sampling_mode: str,
    grid_period_range_d: tuple[float, float],
    grid_radius_range_rearth: tuple[float, float],
    grid_depth_range: tuple[float, float],
    grid_period_bins: int,
    grid_radius_bins: int,
    grid_depth_bins: int,
    grid_depth_spacing: str,
    random_state: int,
    min_in_transit: int,
    max_attempts_per_injection: int,
    adaptive_bkspace: float,
    adaptive_gap_split: float,
    baseline_source: str,
    cadence_s: float,
    target_tmag_bin_edges: tuple[float, ...] | None,
    target_tmag_bin_weights: tuple[float, ...] | None,
    target_tmag_by_tic: dict[int, float] | None,
    target_h5_paths: tuple[Path, ...] | None,
    tic_filter: set[int] | None,
    limit_targets: int | None,
    compression: str | None,
    store_original: bool,
    overwrite: bool,
    progress_every: int,
) -> dict[str, Any]:
    import h5py

    if not apertures:
        raise ValueError("at least one output aperture is required")
    unknown = [ap for ap in apertures if ap not in OUTPUT_TO_RAW_APERTURE]
    if unknown:
        raise ValueError(f"unknown output aperture(s): {unknown}")
    if not families:
        raise ValueError("at least one signal family is required")
    if baseline_source not in BASELINE_SOURCES:
        raise ValueError(f"unknown baseline source: {baseline_source!r}")

    target_paths = _discover_target_paths(
        tuple(Path(p) for p in orbit_roots),
        target_h5_paths=target_h5_paths,
        tic_filter=tic_filter,
        limit_targets=limit_targets,
    )
    if not target_paths:
        raise ValueError("no raw TGLC HDF5 targets found for requested inputs")
    print(f"[predetrend-injections] source targets: {len(target_paths):,}", flush=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_h5 = out_dir / "injected_lightcurves.h5"
    manifest_csv = out_dir / "injection_manifest.csv"
    labels_csv = out_dir / "injection_labels.csv"
    summary_json = out_dir / "summary.json"
    if out_h5.exists() and not overwrite:
        raise FileExistsError(f"output exists; pass --overwrite: {out_h5}")
    tmp_h5 = out_h5.with_suffix(out_h5.suffix + ".tmp")
    if tmp_h5.exists():
        tmp_h5.unlink()

    rng = np.random.default_rng(random_state)
    target_tics = np.asarray(sorted(target_paths.keys()), dtype=np.int64)
    target_tmag_sampler = None
    if target_tmag_bin_edges is not None:
        if target_tmag_by_tic is None:
            target_tmag_by_tic = _read_target_tmag_map(target_paths)
        target_tmag_sampler = _build_target_tmag_sampler(
            target_tmag_by_tic,
            edges=target_tmag_bin_edges,
            weights=target_tmag_bin_weights,
        )
        print(
            "[predetrend-injections] target Tmag sampler: "
            + ", ".join(
                f"{label}: {n_targets} targets, p={prob:.3f}"
                for label, n_targets, prob in zip(
                    target_tmag_sampler["labels"],
                    target_tmag_sampler["n_targets_by_bin"],
                    target_tmag_sampler["probabilities"],
                )
            ),
            flush=True,
        )
    rows: list[dict[str, Any]] = []
    skipped = {
        "read_failed": 0,
        "merge_failed": 0,
        "missing_raw_aperture": 0,
        "nonpositive_baseline": 0,
        "epoch_placement_failed": 0,
        "model_failed": 0,
        "detrend_failed": 0,
    }
    total_attempts = 0
    max_total_attempts = max(int(n_injections) * int(max_attempts_per_injection), int(n_injections))
    created_utc = datetime.now(timezone.utc).isoformat()

    with h5py.File(tmp_h5, "w") as out:
        out.attrs["created_utc"] = created_utc
        out.attrs["injection_level"] = "raw_flux_pre_detrend"
        out.attrs["injection_model"] = "batman_quadratic"
        out.attrs["detrend_methods"] = json.dumps(["twirl-fs-v2", "twirl-fs-v2-adp03q"])
        out.attrs["orbit_roots"] = json.dumps([str(p) for p in orbit_roots])
        out.attrs["aperture"] = apertures[0]
        out.attrs["apertures"] = json.dumps(list(apertures))
        out.attrs["random_state"] = int(random_state)
        out.attrs["time_unit"] = f"BJD - {BJDREFI}"
        out.attrs["families"] = json.dumps(list(families))
        out.attrs["sampling_mode"] = sampling_mode
        out.attrs["grid_period_range_d"] = json.dumps(list(grid_period_range_d))
        out.attrs["grid_radius_range_rearth"] = json.dumps(list(grid_radius_range_rearth))
        out.attrs["grid_depth_range"] = json.dumps(list(grid_depth_range))
        out.attrs["grid_period_bins"] = int(grid_period_bins)
        out.attrs["grid_radius_bins"] = int(grid_radius_bins)
        out.attrs["grid_depth_bins"] = int(grid_depth_bins)
        out.attrs["grid_depth_spacing"] = grid_depth_spacing
        out.attrs["adaptive_bkspace"] = float(adaptive_bkspace)
        out.attrs["adaptive_gap_split"] = float(adaptive_gap_split)
        out.attrs["baseline_source"] = str(baseline_source)
        out.attrs["cadence_s"] = float(cadence_s)
        out.attrs["target_tmag_bin_edges"] = json.dumps(list(target_tmag_bin_edges or ()))
        out.attrs["target_tmag_bin_weights"] = json.dumps(list(target_tmag_bin_weights or ()))
        if target_tmag_sampler is not None:
            out.attrs["target_tmag_bin_labels"] = json.dumps(list(target_tmag_sampler["labels"]))
            out.attrs["target_tmag_bin_target_counts"] = json.dumps(list(target_tmag_sampler["n_targets_by_bin"]))
            out.attrs["target_tmag_bin_probabilities"] = json.dumps(
                [float(value) for value in target_tmag_sampler["probabilities"]]
            )
        out.attrs["limb_dark_u1"] = float(_POSTDET.BATMAN_LIMB_DARKENING[0])
        out.attrs["limb_dark_u2"] = float(_POSTDET.BATMAN_LIMB_DARKENING[1])
        out.attrs["batman_supersample_factor"] = int(_POSTDET.BATMAN_SUPERSAMPLE_FACTOR)
        inj_group = out.create_group("injections")

        while len(rows) < n_injections and total_attempts < max_total_attempts:
            total_attempts += 1
            tic = _draw_target_tic(rng, target_tics, target_tmag_sampler)
            try:
                lcs = [read_tglc_h5(path) for path in target_paths[tic]]
                merged = _TWIRLFS.merge_orbits(lcs)
            except Exception:
                skipped["read_failed"] += 1
                continue
            if merged is None:
                skipped["merge_failed"] += 1
                continue
            if any(ap not in merged.apertures for ap in APERTURE_KEYS):
                skipped["missing_raw_aperture"] += 1
                continue

            params = _draw_params(
                len(rows),
                rng=rng,
                sampling_mode=sampling_mode,
                families=families,
                grid_period_range_d=grid_period_range_d,
                grid_radius_range_rearth=grid_radius_range_rearth,
                grid_depth_range=grid_depth_range,
                grid_period_bins=grid_period_bins,
                grid_radius_bins=grid_radius_bins,
                grid_depth_bins=grid_depth_bins,
                grid_depth_spacing=grid_depth_spacing,
            )
            time = np.asarray(merged.time, dtype=np.float64)
            quality = np.asarray(merged.quality, dtype=np.int32)
            orbitid = np.asarray(getattr(merged, "orbitid", np.zeros(len(time), dtype=np.int16)), dtype=np.int16)
            raw_by_ap = {
                ap: np.asarray(merged.apertures[ap].raw_flux, dtype=np.float64)
                for ap in APERTURE_KEYS
            }
            baseline_diag = _baseline_diagnostics(
                merged,
                quality=quality,
                baseline_source=baseline_source,
                cadence_s=cadence_s,
            )
            raw_baseline_by_ap = baseline_diag["raw_baseline_by_ap"]
            injection_baseline_by_ap = baseline_diag["injection_baseline_by_ap"]
            inferred_baseline_by_ap = baseline_diag["inferred_baseline_by_ap"]
            aperture_fraction_by_ap = baseline_diag["aperture_fraction_by_ap"]
            expected_tess_flux_per_cadence = baseline_diag["expected_tess_flux_per_cadence"]
            if not np.isfinite(injection_baseline_by_ap["Primary"]):
                skipped["nonpositive_baseline"] += 1
                continue
            for ap in APERTURE_KEYS:
                if not np.isfinite(injection_baseline_by_ap[ap]):
                    injection_baseline_by_ap[ap] = injection_baseline_by_ap["Primary"]
                if not np.isfinite(raw_baseline_by_ap[ap]):
                    raw_baseline_by_ap[ap] = raw_baseline_by_ap["Primary"]

            primary_raw = raw_by_ap["Primary"]
            primary_good = (quality == 0) & np.isfinite(time) & np.isfinite(primary_raw)
            try:
                t0_d, mask = choose_observed_epoch(
                    time,
                    period_d=float(params["period_d"]),
                    duration_min=float(params["duration_min"]),
                    rng=rng,
                    quality=quality,
                    finite_mask=np.isfinite(primary_raw),
                    min_in_transit=min_in_transit,
                    max_tries=500,
                )
            except ValueError:
                skipped["epoch_placement_failed"] += 1
                continue
            n_good_in_transit = int(np.count_nonzero(mask & primary_good))
            if n_good_in_transit < min_in_transit:
                skipped["epoch_placement_failed"] += 1
                continue

            injected_raw_by_ap: dict[str, np.ndarray] = {}
            transit_model = None
            try:
                for ap, raw in raw_by_ap.items():
                    baseline = float(injection_baseline_by_ap[ap])
                    if transit_model is None:
                        injected_raw, ap_mask, transit_model = inject_batman_transit(
                            time,
                            raw,
                            period_d=float(params["period_d"]),
                            t0_d=t0_d,
                            duration_min=float(params["duration_min"]),
                            radius_rstar=float(params["radius_rwd"]),
                            a_over_rstar=float(params["a_over_rwd"]),
                            impact_b=float(params["impact_b"]),
                            baseline=baseline,
                            limb_darkening=_POSTDET.BATMAN_LIMB_DARKENING,
                            supersample_factor=_POSTDET.BATMAN_SUPERSAMPLE_FACTOR,
                        )
                        mask = ap_mask
                    else:
                        injected_raw = np.asarray(raw, dtype=np.float64) + baseline * (transit_model - 1.0)
                    injected_raw_by_ap[ap] = injected_raw
            except Exception:
                skipped["model_failed"] += 1
                continue

            if transit_model is None:
                skipped["model_failed"] += 1
                continue
            sampled_model_depth = float(max(0.0, 1.0 - np.nanmin(transit_model)))
            duration_d = float(params["duration_min"]) / 1440.0
            dense_time = np.linspace(t0_d - 2.0 * duration_d, t0_d + 2.0 * duration_d, 2000)
            try:
                dense_model = batman_transit_model(
                    dense_time,
                    period_d=float(params["period_d"]),
                    t0_d=t0_d,
                    radius_rstar=float(params["radius_rwd"]),
                    a_over_rstar=float(params["a_over_rwd"]),
                    impact_b=float(params["impact_b"]),
                    limb_darkening=_POSTDET.BATMAN_LIMB_DARKENING,
                    supersample_factor=_POSTDET.BATMAN_SUPERSAMPLE_FACTOR,
                )
            except Exception:
                skipped["model_failed"] += 1
                continue
            params["model_depth"] = float(max(0.0, 1.0 - np.nanmin(dense_model)))
            params["sampled_model_depth"] = sampled_model_depth

            try:
                original_columns, original_diag = _detrended_output_columns(
                    merged,
                    adaptive_bkspace=adaptive_bkspace,
                    adaptive_gap_split=adaptive_gap_split,
                )
                injected_lc = _clone_with_raw_flux(merged, injected_raw_by_ap)
                injected_columns, injected_diag = _detrended_output_columns(
                    injected_lc,
                    adaptive_bkspace=adaptive_bkspace,
                    adaptive_gap_split=adaptive_gap_split,
                )
            except Exception:
                skipped["detrend_failed"] += 1
                continue

            injection_id = f"predet_{len(rows):06d}"
            group = inj_group.create_group(injection_id)
            group.attrs["injection_id"] = injection_id
            group.attrs["injection_level"] = "raw_flux_pre_detrend"
            group.attrs["tic"] = int(merged.tic)
            group.attrs["sector"] = int(merged.sector)
            group.attrs["camera"] = int(merged.cam)
            group.attrs["ccd"] = int(merged.ccd)
            group.attrs["tessmag"] = float(merged.tmag)
            group.attrs["aperture"] = apertures[0]
            group.attrs["apertures"] = json.dumps(list(apertures))
            group.attrs["raw_apertures"] = json.dumps(list(APERTURE_KEYS))
            group.attrs["signal_family"] = str(params["signal_family"])
            group.attrs["period_d"] = float(params["period_d"])
            group.attrs["t0_d"] = float(t0_d)
            group.attrs["t0_bjd"] = float(t0_d + BJDREFI)
            group.attrs["duration_min"] = float(params["duration_min"])
            group.attrs["depth"] = float(params["depth"])
            group.attrs["n_cadences"] = int(len(time))
            group.attrs["n_quality_zero"] = int(np.count_nonzero(quality == 0))
            group.attrs["n_in_transit"] = int(np.count_nonzero(mask))
            group.attrs["n_good_in_transit"] = n_good_in_transit
            for key in PARAM_ATTR_KEYS:
                if params.get(key, "") != "":
                    group.attrs[key] = _scalar(params[key])
            group.attrs["limb_dark_u1"] = float(_POSTDET.BATMAN_LIMB_DARKENING[0])
            group.attrs["limb_dark_u2"] = float(_POSTDET.BATMAN_LIMB_DARKENING[1])
            group.attrs["batman_supersample_factor"] = int(_POSTDET.BATMAN_SUPERSAMPLE_FACTOR)
            group.attrs["baseline"] = float(injection_baseline_by_ap["Primary"])
            group.attrs["baseline_source"] = str(baseline_source)
            group.attrs["cadence_s"] = float(cadence_s)
            group.attrs["expected_tess_flux_per_cadence"] = float(expected_tess_flux_per_cadence)
            for ap, baseline in raw_baseline_by_ap.items():
                group.attrs[f"raw_baseline_{ap}"] = float(baseline)
            for ap, baseline in inferred_baseline_by_ap.items():
                group.attrs[f"inferred_baseline_{ap}"] = float(baseline)
            for ap, baseline in injection_baseline_by_ap.items():
                group.attrs[f"injection_baseline_{ap}"] = float(baseline)
            for ap, fraction in aperture_fraction_by_ap.items():
                group.attrs[f"aperture_fraction_{ap}"] = float(fraction)
            for out_name, diag in injected_diag.items():
                group.attrs[f"{out_name}_cotrend_status"] = str(diag.get("cotrend_status", ""))
            group.attrs["label"] = "planet_like"
            group.attrs["label_source"] = "injection"
            group.attrs["split"] = _split_for_tic(int(merged.tic))
            group.attrs["source_target_paths"] = json.dumps([str(path) for path in target_paths[tic]])

            _dataset(group, "time", time.astype(np.float64), compression)
            _dataset(group, "quality", quality.astype(np.int32), compression)
            _dataset(group, "orbitid", orbitid.astype(np.int16), compression)
            _dataset(group, "in_transit", np.asarray(mask, dtype=np.bool_), compression)
            _dataset(group, "transit_model", transit_model.astype(np.float32), compression)
            for ap in APERTURE_KEYS:
                _dataset(group, f"RAW_FLUX_{ap}_injected", injected_raw_by_ap[ap].astype(np.float32), compression)
                if store_original:
                    _dataset(group, f"RAW_FLUX_{ap}_original", raw_by_ap[ap].astype(np.float32), compression)
            for aperture in apertures:
                _dataset(group, f"{aperture}_injected", injected_columns[aperture].astype(np.float32), compression)
                if store_original:
                    _dataset(group, f"{aperture}_original", original_columns[aperture].astype(np.float32), compression)
            _dataset(group, "flux_injected", injected_columns[apertures[0]].astype(np.float32), compression)
            if store_original:
                _dataset(group, "flux_original", original_columns[apertures[0]].astype(np.float32), compression)

            row = {
                "injection_id": injection_id,
                "tic": int(merged.tic),
                "sector": int(merged.sector),
                "camera": int(merged.cam),
                "ccd": int(merged.ccd),
                "tessmag": float(merged.tmag),
                "aperture": apertures[0],
                "apertures": "|".join(apertures),
                "signal_family": str(params["signal_family"]),
                "period_d": float(params["period_d"]),
                "t0_d": float(t0_d),
                "t0_bjd": float(t0_d + BJDREFI),
                "duration_min": float(params["duration_min"]),
                "depth": float(params["depth"]),
                "n_good_in_transit": n_good_in_transit,
                "baseline": float(injection_baseline_by_ap["Primary"]),
                "baseline_source": str(baseline_source),
                "expected_tess_flux_per_cadence": float(expected_tess_flux_per_cadence),
                "split": _split_for_tic(int(merged.tic)),
                "label": "planet_like",
                "label_source": "injection",
                "source_target_paths": "|".join(str(path) for path in target_paths[tic]),
                "h5_group": f"/injections/{injection_id}",
            }
            for key in PARAM_ATTR_KEYS:
                row[key] = _scalar(params.get(key, ""))
            rows.append(row)
            if progress_every > 0 and len(rows) % progress_every == 0:
                print(
                    f"  [predetrend-injections] accepted {len(rows):,}/{n_injections:,} "
                    f"after {total_attempts:,} attempts",
                    flush=True,
                )

        out.attrs["n_injections"] = int(len(rows))
        out.attrs["total_attempts"] = int(total_attempts)
        out.attrs["skipped"] = json.dumps(skipped, sort_keys=True)

    if len(rows) < n_injections:
        tmp_h5.unlink(missing_ok=True)
        raise RuntimeError(
            f"created {len(rows)} / {n_injections} injections after "
            f"{total_attempts} attempts; skipped={skipped}"
        )
    os.replace(tmp_h5, out_h5)
    _write_csv(manifest_csv, rows, MANIFEST_COLUMNS)
    _write_csv(labels_csv, rows, ("injection_id", "tic", "label", "label_source", "split", "signal_family"))
    summary = {
        "created_utc": created_utc,
        "injection_level": "raw_flux_pre_detrend",
        "out_h5": str(out_h5),
        "manifest_csv": str(manifest_csv),
        "labels_csv": str(labels_csv),
        "orbit_roots": [str(p) for p in orbit_roots],
        "n_source_targets": int(len(target_paths)),
        "n_injections": int(len(rows)),
        "apertures": list(apertures),
        "raw_apertures": list(APERTURE_KEYS),
        "sampling_mode": sampling_mode,
        "grid_period_range_d": list(grid_period_range_d),
        "grid_radius_range_rearth": list(grid_radius_range_rearth),
        "grid_depth_range": list(grid_depth_range),
        "grid_period_bins": int(grid_period_bins),
        "grid_radius_bins": int(grid_radius_bins),
        "grid_depth_bins": int(grid_depth_bins),
        "grid_depth_spacing": grid_depth_spacing,
        "adaptive_bkspace": float(adaptive_bkspace),
        "adaptive_gap_split": float(adaptive_gap_split),
        "baseline_source": str(baseline_source),
        "cadence_s": float(cadence_s),
        "target_tmag_bin_edges": list(target_tmag_bin_edges or ()),
        "target_tmag_bin_weights": list(target_tmag_bin_weights or ()),
        "target_tmag_sampler": (
            {
                "labels": list(target_tmag_sampler["labels"]),
                "n_targets_by_bin": list(target_tmag_sampler["n_targets_by_bin"]),
                "probabilities": [float(value) for value in target_tmag_sampler["probabilities"]],
            }
            if target_tmag_sampler is not None
            else None
        ),
        "random_state": int(random_state),
        "min_in_transit": int(min_in_transit),
        "total_attempts": int(total_attempts),
        "skipped": skipped,
    }
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def _parse_range(raw: str) -> tuple[float, float]:
    parts = [p.strip() for p in str(raw).split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("range must be 'lo,hi'")
    lo, hi = float(parts[0]), float(parts[1])
    if hi <= lo:
        raise argparse.ArgumentTypeError("range upper bound must exceed lower bound")
    return lo, hi


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--orbit-roots", type=Path, nargs="+", default=list(DEFAULT_ORBIT_ROOTS))
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--apertures", nargs="+", default=list(DEFAULT_OUTPUT_APERTURES))
    ap.add_argument("--n-injections", type=int, default=100)
    ap.add_argument("--families", nargs="+", default=list(_POSTDET.SIGNAL_FAMILIES))
    ap.add_argument(
        "--sampling-mode",
        choices=("family_cycle", "period_radius_grid", "period_depth_grid"),
        default="period_depth_grid",
    )
    ap.add_argument("--grid-period-range-d", type=_parse_range, default=_POSTDET.DEFAULT_GRID_PERIOD_RANGE_D)
    ap.add_argument("--grid-radius-range-rearth", type=_parse_range, default=_POSTDET.DEFAULT_GRID_RADIUS_RANGE_REARTH)
    ap.add_argument("--grid-depth-range", type=_parse_range, default=_POSTDET.DEFAULT_GRID_DEPTH_RANGE)
    ap.add_argument("--grid-period-bins", type=int, default=_POSTDET.DEFAULT_GRID_PERIOD_BINS)
    ap.add_argument("--grid-radius-bins", type=int, default=_POSTDET.DEFAULT_GRID_RADIUS_BINS)
    ap.add_argument("--grid-depth-bins", type=int, default=_POSTDET.DEFAULT_GRID_DEPTH_BINS)
    ap.add_argument("--grid-depth-spacing", choices=("linear", "log"), default="linear")
    ap.add_argument("--random-state", type=int, default=5605)
    ap.add_argument("--min-in-transit", type=int, default=2)
    ap.add_argument("--max-attempts-per-injection", type=int, default=50)
    ap.add_argument("--adaptive-bkspace", type=float, default=0.3)
    ap.add_argument("--adaptive-gap-split", type=float, default=0.2)
    ap.add_argument(
        "--baseline-source",
        choices=BASELINE_SOURCES,
        default="raw_median",
        help=(
            "Flux scale for BATMAN injection. raw_median preserves the first "
            "pre-detrend pilot behavior; inferred_aperture estimates the "
            "target aperture flux from RawMagnitude/RawFlux/TESSMAG; "
            "tessmag_total injects the full catalog target flux into every aperture."
        ),
    )
    ap.add_argument("--cadence-s", type=float, default=DEFAULT_CADENCE_S)
    ap.add_argument(
        "--target-tmag-bin-edges",
        type=_parse_tmag_bin_edges,
        default=None,
        help=(
            "Optional comma-separated Tmag bin edges for target sampling, e.g. "
            "'0,17,18,19,99'. When provided, targets are drawn by Tmag bin "
            "instead of uniformly over all available S56 targets."
        ),
    )
    ap.add_argument(
        "--target-tmag-bin-weights",
        type=_parse_float_list,
        default=None,
        help=(
            "Optional comma-separated non-negative target-bin weights. Length "
            "must be one fewer than --target-tmag-bin-edges. If omitted, "
            "nonempty Tmag bins are sampled equally."
        ),
    )
    ap.add_argument("--candidate-table", type=Path, default=None)
    ap.add_argument("--tic-column", default="tic")
    ap.add_argument(
        "--target-tmag-table",
        type=Path,
        default=None,
        help=(
            "Optional table with TIC and Tmag columns used by --target-tmag-bin-edges. "
            "This avoids opening every raw HDF5 file just to read TessMag."
        ),
    )
    ap.add_argument("--target-tmag-tic-column", default="tic")
    ap.add_argument("--target-tmag-column", default="tessmag")
    ap.add_argument(
        "--target-h5-paths",
        type=Path,
        nargs="+",
        default=None,
        help=(
            "Optional explicit raw TGLC HDF5 paths. Useful for smoke tests; "
            "paths are grouped by TIC filename stem instead of discovering the full tree."
        ),
    )
    ap.add_argument(
        "--target-h5-list",
        type=Path,
        default=None,
        help="Optional text file with one raw TGLC HDF5 path per line; combined with --target-h5-paths.",
    )
    ap.add_argument("--limit-targets", type=int, default=None)
    ap.add_argument("--compression", choices=("lzf", "gzip", "none"), default="lzf")
    ap.add_argument("--store-original", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--progress-every", type=int, default=50)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    compression = None if args.compression == "none" else args.compression
    tic_filter = _read_tic_filter(args.candidate_table, args.tic_column)
    target_tmag_by_tic = _read_target_tmag_table(
        args.target_tmag_table,
        args.target_tmag_tic_column,
        args.target_tmag_column,
    )
    target_h5_paths = tuple(args.target_h5_paths or ()) + _read_path_list(args.target_h5_list)
    summary = make_predetrend_injections(
        orbit_roots=tuple(args.orbit_roots),
        out_dir=args.out_dir,
        apertures=tuple(args.apertures),
        n_injections=args.n_injections,
        families=tuple(args.families),
        sampling_mode=args.sampling_mode,
        grid_period_range_d=tuple(args.grid_period_range_d),
        grid_radius_range_rearth=tuple(args.grid_radius_range_rearth),
        grid_depth_range=tuple(args.grid_depth_range),
        grid_period_bins=args.grid_period_bins,
        grid_radius_bins=args.grid_radius_bins,
        grid_depth_bins=args.grid_depth_bins,
        grid_depth_spacing=args.grid_depth_spacing,
        random_state=args.random_state,
        min_in_transit=args.min_in_transit,
        max_attempts_per_injection=args.max_attempts_per_injection,
        adaptive_bkspace=args.adaptive_bkspace,
        adaptive_gap_split=args.adaptive_gap_split,
        baseline_source=args.baseline_source,
        cadence_s=args.cadence_s,
        target_tmag_bin_edges=args.target_tmag_bin_edges,
        target_tmag_bin_weights=args.target_tmag_bin_weights,
        target_tmag_by_tic=target_tmag_by_tic,
        target_h5_paths=target_h5_paths if target_h5_paths else None,
        tic_filter=tic_filter,
        limit_targets=args.limit_targets,
        compression=compression,
        store_original=args.store_original,
        overwrite=args.overwrite,
        progress_every=args.progress_every,
    )
    print("[predetrend-injections] complete")
    print(f"  injections: {summary['n_injections']:,}")
    print(f"  source targets: {summary['n_source_targets']:,}")
    print(f"  apertures: {summary['apertures']}")
    print(f"  skipped: {summary['skipped']}")
    print(f"  out: {summary['out_h5']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
