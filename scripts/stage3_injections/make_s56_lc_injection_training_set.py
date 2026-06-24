#!/usr/bin/env python3
"""Build LC-level injected examples from a compact S56 TWIRL-FS export.

This creates positive examples for the S56 training program. These labels are
for injected light curves; they should be folded into candidate-level
self-training only after the injected products pass through search, vetting,
and candidate consolidation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import APERTURES, BJDREFI
from twirl.search.injections import batman_transit_model, choose_observed_epoch, inject_batman_transit

DEFAULT_EXPORT_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export.h5"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT / "data_local/stage3_injections/s56_twirlfs_v2_injection_training"
)

SIGNAL_FAMILIES = (
    "wd1856_like",
    "short_deep",
    "roche_boundary",
    "wd_small_body",
    "wd_earth_size",
    "wd_giant_or_bd",
    "wd_long_period",
    "wd_roche_edge",
)

WD_RADIUS_REARTH = 0.013 * 109.076
WD_DENSITY_G_CM3 = 3.85e5
G_SI = 6.67430e-11
DEFAULT_GRID_PERIOD_RANGE_D = (0.08, 13.0)
DEFAULT_GRID_RADIUS_RANGE_REARTH = (0.08, 16.8)
DEFAULT_GRID_DEPTH_RANGE = (0.003, 0.995)
DEFAULT_GRID_PERIOD_BINS = 10
DEFAULT_GRID_RADIUS_BINS = 10
DEFAULT_GRID_DEPTH_BINS = 10
MAX_TRANSIT_DUTY_CYCLE = 0.20
BATMAN_LIMB_DARKENING = (0.05, 0.05)
BATMAN_SUPERSAMPLE_FACTOR = 7


def _log_uniform(rng: np.random.Generator, lo: float, hi: float) -> float:
    return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))


def _a_over_rwd_from_period(period_d: float) -> float:
    rho_kg_m3 = WD_DENSITY_G_CM3 * 1000.0
    period_s = float(period_d) * 86400.0
    return float((G_SI * rho_kg_m3 * period_s**2 / (3.0 * np.pi)) ** (1.0 / 3.0))


def _duration_from_geometry_values(
    *,
    period_d: float,
    radius_rwd: float,
    impact_b: float,
    a_over_rwd: float,
) -> float:
    chord = float(np.sqrt(max((1.0 + radius_rwd) ** 2 - impact_b**2, 1.0e-6)))
    duration_d = float(period_d) / np.pi * chord / max(a_over_rwd, 1.0e-6)
    max_duration_min = min(180.0, MAX_TRANSIT_DUTY_CYCLE * float(period_d) * 1440.0)
    return float(np.clip(duration_d * 1440.0, 2.0, max(2.0, max_duration_min)))


def _duration_from_wd_geometry(
    *,
    period_d: float,
    radius_rearth: float,
    rng: np.random.Generator,
) -> tuple[float, float, float, float]:
    """Draw a physically motivated WD transit duration.

    This keeps the injected duration scale tied to WD density, period,
    companion size, and impact parameter instead of drawing arbitrary
    durations.
    """
    radius_rwd = float(radius_rearth) / WD_RADIUS_REARTH
    a_over_rwd = _a_over_rwd_from_period(period_d)
    max_b = max(0.05, 1.0 + radius_rwd)
    impact_b = float(rng.uniform(0.0, min(0.95 * max_b, 0.95 * a_over_rwd)))
    duration_min = _duration_from_geometry_values(
        period_d=period_d,
        radius_rwd=radius_rwd,
        impact_b=impact_b,
        a_over_rwd=a_over_rwd,
    )
    return duration_min, radius_rwd, impact_b, a_over_rwd


def _draw_physical_wd_params(
    *,
    family: str,
    rng: np.random.Generator,
    period_range: tuple[float, float],
    radius_range_rearth: tuple[float, float],
) -> dict[str, float | str]:
    radius_rearth = _log_uniform(rng, *radius_range_rearth)
    period_d = _log_uniform(rng, *period_range)
    duration_min, radius_rwd, impact_b, a_over_rwd = _duration_from_wd_geometry(
        period_d=period_d,
        radius_rearth=radius_rearth,
        rng=rng,
    )
    depth = float(min(radius_rwd**2, 0.995))
    return {
        "signal_family": family,
        "period_d": period_d,
        "duration_min": duration_min,
        "depth": depth,
        "radius_rearth": radius_rearth,
        "radius_rwd": radius_rwd,
        "impact_b": impact_b,
        "a_over_rwd": a_over_rwd,
        "inclination_deg": float(np.degrees(np.arccos(np.clip(impact_b / a_over_rwd, 0.0, 1.0)))),
        "duration_model": "wd_density_batman",
        "injection_model": "batman_quadratic",
    }


def _draw_params(family: str, rng: np.random.Generator) -> dict[str, float | str]:
    if family == "wd1856_like":
        return {
            "signal_family": family,
            "period_d": _log_uniform(rng, 0.8, 2.5),
            "duration_min": float(rng.uniform(4.0, 12.0)),
            "depth": float(rng.uniform(0.05, 0.55)),
            "radius_rearth": "",
            "radius_rwd": "",
            "impact_b": "",
            "duration_model": "empirical_uniform",
        }
    if family == "short_deep":
        return {
            "signal_family": family,
            "period_d": _log_uniform(rng, 0.12, 1.5),
            "duration_min": float(rng.uniform(3.0, 20.0)),
            "depth": float(rng.uniform(0.03, 0.50)),
            "radius_rearth": "",
            "radius_rwd": "",
            "impact_b": "",
            "duration_model": "empirical_uniform",
        }
    if family == "roche_boundary":
        return {
            "signal_family": family,
            "period_d": _log_uniform(rng, 0.18, 0.35),
            "duration_min": float(rng.uniform(3.0, 15.0)),
            "depth": float(rng.uniform(0.02, 0.30)),
            "radius_rearth": "",
            "radius_rwd": "",
            "impact_b": "",
            "duration_model": "empirical_uniform",
        }
    if family == "wd_small_body":
        return _draw_physical_wd_params(
            family=family,
            rng=rng,
            period_range=(0.08, 13.0),
            radius_range_rearth=(0.08, 0.35),
        )
    if family == "wd_earth_size":
        return _draw_physical_wd_params(
            family=family,
            rng=rng,
            period_range=(0.08, 13.0),
            radius_range_rearth=(0.35, 1.2),
        )
    if family == "wd_giant_or_bd":
        return _draw_physical_wd_params(
            family=family,
            rng=rng,
            period_range=(0.08, 13.0),
            radius_range_rearth=(1.2, DEFAULT_GRID_RADIUS_RANGE_REARTH[1]),
        )
    if family == "wd_long_period":
        return _draw_physical_wd_params(
            family=family,
            rng=rng,
            period_range=(2.5, 13.0),
            radius_range_rearth=(0.15, 6.0),
        )
    if family == "wd_roche_edge":
        return _draw_physical_wd_params(
            family=family,
            rng=rng,
            period_range=(0.08, 0.25),
            radius_range_rearth=(0.08, 3.0),
        )
    raise ValueError(f"unknown signal family: {family}")


def _log_edges(lo: float, hi: float, bins: int) -> np.ndarray:
    if lo <= 0 or hi <= lo:
        raise ValueError(f"invalid log range: {lo}, {hi}")
    if bins <= 0:
        raise ValueError("number of bins must be positive")
    return np.exp(np.linspace(np.log(lo), np.log(hi), bins + 1))


def _depth_edges(lo: float, hi: float, bins: int, spacing: str) -> np.ndarray:
    if spacing == "log":
        return _log_edges(lo, hi, bins)
    if lo < 0 or hi <= lo:
        raise ValueError(f"invalid depth range: {lo}, {hi}")
    if bins <= 0:
        raise ValueError("number of bins must be positive")
    return np.linspace(lo, hi, bins + 1)


def _circle_overlap_depth(radius_rwd: float, impact_b: float) -> float:
    """Uniform-brightness overlap depth for two circles, normalized by WD area."""
    r = float(radius_rwd)
    d = float(abs(impact_b))
    if r <= 0:
        return 0.0
    if d >= 1.0 + r:
        return 0.0
    if d <= abs(1.0 - r):
        return float(min(r * r, 1.0))
    cos1 = np.clip((d * d + 1.0 - r * r) / (2.0 * d), -1.0, 1.0)
    cos2 = np.clip((d * d + r * r - 1.0) / (2.0 * d * r), -1.0, 1.0)
    area = (
        np.arccos(cos1)
        + r * r * np.arccos(cos2)
        - 0.5 * np.sqrt(max((-d + 1.0 + r) * (d + 1.0 - r) * (d - 1.0 + r) * (d + 1.0 + r), 0.0))
    )
    return float(np.clip(area / np.pi, 0.0, 1.0))


def _impact_for_geometric_depth(radius_rwd: float, target_depth: float) -> float:
    central_depth = _circle_overlap_depth(radius_rwd, 0.0)
    target = float(np.clip(target_depth, 0.0, central_depth))
    if target >= central_depth:
        return 0.0
    if target <= 0:
        return float(1.0 + radius_rwd)
    lo = 0.0
    hi = float(1.0 + radius_rwd)
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        depth = _circle_overlap_depth(radius_rwd, mid)
        if depth > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _grid_family(period_d: float, radius_rearth: float) -> str:
    if period_d < 0.25:
        return "wd_roche_edge"
    if period_d > 2.5:
        return "wd_long_period"
    if radius_rearth < 0.35:
        return "wd_small_body"
    if radius_rearth < 1.2:
        return "wd_earth_size"
    return "wd_giant_or_bd"


def _draw_grid_params(
    injection_index: int,
    *,
    rng: np.random.Generator,
    period_range: tuple[float, float],
    radius_range_rearth: tuple[float, float],
    period_bins: int,
    radius_bins: int,
) -> dict[str, float | int | str]:
    period_edges = _log_edges(*period_range, period_bins)
    radius_edges = _log_edges(*radius_range_rearth, radius_bins)
    n_cells = period_bins * radius_bins
    cell = int(injection_index) % n_cells
    period_bin = cell // radius_bins
    radius_bin = cell % radius_bins
    period_d = _log_uniform(rng, float(period_edges[period_bin]), float(period_edges[period_bin + 1]))
    radius_rearth = _log_uniform(rng, float(radius_edges[radius_bin]), float(radius_edges[radius_bin + 1]))
    duration_min, radius_rwd, impact_b, a_over_rwd = _duration_from_wd_geometry(
        period_d=period_d,
        radius_rearth=radius_rearth,
        rng=rng,
    )
    depth = float(min(radius_rwd**2, 0.995))
    family = _grid_family(period_d, radius_rearth)
    return {
        "signal_family": family,
        "sampling_mode": "period_radius_grid",
        "grid_period_bin": period_bin,
        "grid_radius_bin": radius_bin,
        "grid_cell_id": f"p{period_bin:02d}_r{radius_bin:02d}",
        "period_bin_lo_d": float(period_edges[period_bin]),
        "period_bin_hi_d": float(period_edges[period_bin + 1]),
        "radius_bin_lo_rearth": float(radius_edges[radius_bin]),
        "radius_bin_hi_rearth": float(radius_edges[radius_bin + 1]),
        "period_d": period_d,
        "duration_min": duration_min,
        "depth": depth,
        "radius_rearth": radius_rearth,
        "radius_rwd": radius_rwd,
        "impact_b": impact_b,
        "a_over_rwd": a_over_rwd,
        "inclination_deg": float(np.degrees(np.arccos(np.clip(impact_b / a_over_rwd, 0.0, 1.0)))),
        "geometric_depth": depth,
        "target_depth": depth,
        "duration_model": "wd_density_batman",
        "injection_model": "batman_quadratic",
    }


def _draw_depth_grid_params(
    injection_index: int,
    *,
    rng: np.random.Generator,
    period_range: tuple[float, float],
    depth_range: tuple[float, float],
    radius_range_rearth: tuple[float, float],
    period_bins: int,
    depth_bins: int,
    depth_spacing: str,
) -> dict[str, float | int | str]:
    period_edges = _log_edges(*period_range, period_bins)
    depth_edges = _depth_edges(*depth_range, depth_bins, depth_spacing)
    n_cells = period_bins * depth_bins
    cell = int(injection_index) % n_cells
    period_bin = cell // depth_bins
    depth_bin = cell % depth_bins
    period_d = _log_uniform(rng, float(period_edges[period_bin]), float(period_edges[period_bin + 1]))
    if depth_spacing == "log":
        target_depth = _log_uniform(rng, float(depth_edges[depth_bin]), float(depth_edges[depth_bin + 1]))
    else:
        target_depth = float(rng.uniform(float(depth_edges[depth_bin]), float(depth_edges[depth_bin + 1])))

    min_radius_rwd = max(float(radius_range_rearth[0]) / WD_RADIUS_REARTH, np.sqrt(min(target_depth, 0.999)))
    max_radius_rwd = float(radius_range_rearth[1]) / WD_RADIUS_REARTH
    if min_radius_rwd > max_radius_rwd:
        min_radius_rwd = max_radius_rwd
        target_depth = min(target_depth, _circle_overlap_depth(min_radius_rwd, 0.0))
    radius_rwd = _log_uniform(rng, min_radius_rwd, max_radius_rwd)
    radius_rearth = float(radius_rwd * WD_RADIUS_REARTH)
    impact_b = _impact_for_geometric_depth(radius_rwd, target_depth)
    a_over_rwd = _a_over_rwd_from_period(period_d)
    impact_b = float(min(impact_b, 0.98 * a_over_rwd, 0.999 * (1.0 + radius_rwd)))
    geometric_depth = _circle_overlap_depth(radius_rwd, impact_b)
    duration_min = _duration_from_geometry_values(
        period_d=period_d,
        radius_rwd=radius_rwd,
        impact_b=impact_b,
        a_over_rwd=a_over_rwd,
    )
    family = _grid_family(period_d, radius_rearth)
    return {
        "signal_family": family,
        "sampling_mode": "period_depth_grid",
        "grid_period_bin": period_bin,
        "grid_depth_bin": depth_bin,
        "grid_cell_id": f"p{period_bin:02d}_d{depth_bin:02d}",
        "period_bin_lo_d": float(period_edges[period_bin]),
        "period_bin_hi_d": float(period_edges[period_bin + 1]),
        "depth_bin_lo": float(depth_edges[depth_bin]),
        "depth_bin_hi": float(depth_edges[depth_bin + 1]),
        "period_d": period_d,
        "duration_min": duration_min,
        "depth": geometric_depth,
        "target_depth": target_depth,
        "geometric_depth": geometric_depth,
        "radius_rearth": radius_rearth,
        "radius_rwd": radius_rwd,
        "impact_b": impact_b,
        "a_over_rwd": a_over_rwd,
        "inclination_deg": float(np.degrees(np.arccos(np.clip(impact_b / a_over_rwd, 0.0, 1.0)))),
        "duration_model": "wd_density_batman",
        "injection_model": "batman_quadratic",
    }


def _split_for_tic(tic: int) -> str:
    digest = hashlib.sha1(str(int(tic)).encode("ascii")).hexdigest()
    bucket = int(digest[:8], 16) % 100
    if bucket < 70:
        return "train"
    if bucket < 85:
        return "validation"
    return "test"


def _dataset(group, name: str, values, compression: str | None) -> None:
    kwargs = {"shuffle": True}
    if compression:
        kwargs["compression"] = compression
    group.create_dataset(name, data=values, **kwargs)


def _target_keys(h5) -> list[str]:
    if "targets" not in h5:
        raise KeyError("input export is missing /targets group")
    return sorted(h5["targets"].keys())


def _make_injections(
    *,
    export_h5: Path,
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
    compression: str | None,
    store_original: bool,
    overwrite: bool,
    progress_every: int,
) -> dict:
    import h5py
    import pandas as pd

    if not export_h5.exists():
        raise FileNotFoundError(f"missing compact LC export: {export_h5}")
    if not apertures:
        raise ValueError("at least one aperture is required")
    primary_aperture = apertures[0]
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
    created_utc = datetime.now(timezone.utc).isoformat()
    rows: list[dict] = []
    skipped = {
        "missing_aperture": 0,
        "no_good_flux": 0,
        "epoch_placement_failed": 0,
        "model_failed": 0,
    }
    total_attempts = 0
    max_total_attempts = max(n_injections * max_attempts_per_injection, n_injections)

    with h5py.File(export_h5, "r") as src, h5py.File(tmp_h5, "w") as out:
        keys = _target_keys(src)
        if not keys:
            raise ValueError("input export contains no targets")
        out.attrs["created_utc"] = created_utc
        out.attrs["source_export_h5"] = str(export_h5)
        out.attrs["aperture"] = primary_aperture
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
        out.attrs["injection_model"] = "batman_quadratic"
        out.attrs["limb_dark_u1"] = float(BATMAN_LIMB_DARKENING[0])
        out.attrs["limb_dark_u2"] = float(BATMAN_LIMB_DARKENING[1])
        out.attrs["batman_supersample_factor"] = int(BATMAN_SUPERSAMPLE_FACTOR)
        inj_group = out.create_group("injections")

        while len(rows) < n_injections and total_attempts < max_total_attempts:
            total_attempts += 1
            target_key = str(rng.choice(keys))
            target = src["targets"][target_key]
            missing = [ap for ap in apertures if ap not in target]
            if missing:
                skipped["missing_aperture"] += 1
                continue

            if sampling_mode == "period_depth_grid":
                params = _draw_depth_grid_params(
                    len(rows),
                    rng=rng,
                    period_range=grid_period_range_d,
                    depth_range=grid_depth_range,
                    radius_range_rearth=grid_radius_range_rearth,
                    period_bins=grid_period_bins,
                    depth_bins=grid_depth_bins,
                    depth_spacing=grid_depth_spacing,
                )
            elif sampling_mode == "period_radius_grid":
                params = _draw_grid_params(
                    len(rows),
                    rng=rng,
                    period_range=grid_period_range_d,
                    radius_range_rearth=grid_radius_range_rearth,
                    period_bins=grid_period_bins,
                    radius_bins=grid_radius_bins,
                )
            else:
                family = families[len(rows) % len(families)]
                params = _draw_params(family, rng)
                params["sampling_mode"] = "family_cycle"
                if params.get("radius_rwd", "") == "":
                    radius_rwd = float(np.sqrt(max(float(params["depth"]), 1.0e-6)))
                    period_d = float(params["period_d"])
                    a_over_rwd = _a_over_rwd_from_period(period_d)
                    impact_b = float(rng.uniform(0.0, min(0.8 * (1.0 + radius_rwd), 0.8 * a_over_rwd)))
                    params["radius_rwd"] = radius_rwd
                    params["radius_rearth"] = radius_rwd * WD_RADIUS_REARTH
                    params["impact_b"] = impact_b
                    params["a_over_rwd"] = a_over_rwd
                    params["inclination_deg"] = float(np.degrees(np.arccos(np.clip(impact_b / a_over_rwd, 0.0, 1.0))))
                    params["geometric_depth"] = float(_circle_overlap_depth(radius_rwd, impact_b))
                    params["target_depth"] = float(params["depth"])
                    params["duration_model"] = "empirical_batman"
                    params["injection_model"] = "batman_quadratic"
            time = np.asarray(target["time"], dtype=np.float64)
            quality = np.asarray(target["quality"], dtype=np.int32)
            orbitid = np.asarray(target["orbitid"], dtype=np.int16)
            flux_by_aperture = {
                ap: np.asarray(target[ap], dtype=np.float64)
                for ap in apertures
            }
            flux = flux_by_aperture[primary_aperture]
            good = (quality == 0) & np.isfinite(time) & np.isfinite(flux)
            if not np.any(good):
                skipped["no_good_flux"] += 1
                continue

            baseline = float(np.nanmedian(flux[good]))
            try:
                t0_d, mask = choose_observed_epoch(
                    time,
                    period_d=float(params["period_d"]),
                    duration_min=float(params["duration_min"]),
                    rng=rng,
                    quality=quality,
                    finite_mask=np.isfinite(flux),
                    min_in_transit=min_in_transit,
                    max_tries=500,
                )
            except ValueError:
                skipped["epoch_placement_failed"] += 1
                continue
            n_good_in_transit = int(np.count_nonzero(mask & good))
            if n_good_in_transit < min_in_transit:
                skipped["epoch_placement_failed"] += 1
                continue

            injected_by_aperture = {}
            baseline_by_aperture = {}
            mask = None
            transit_model = None
            try:
                for ap, ap_flux in flux_by_aperture.items():
                    ap_good = (quality == 0) & np.isfinite(time) & np.isfinite(ap_flux)
                    ap_baseline = float(np.nanmedian(ap_flux[ap_good])) if np.any(ap_good) else float("nan")
                    if not np.isfinite(ap_baseline):
                        ap_baseline = baseline
                    if transit_model is None:
                        injected_flux, ap_mask, transit_model = inject_batman_transit(
                            time,
                            ap_flux,
                            period_d=float(params["period_d"]),
                            t0_d=t0_d,
                            duration_min=float(params["duration_min"]),
                            radius_rstar=float(params["radius_rwd"]),
                            a_over_rstar=float(params["a_over_rwd"]),
                            impact_b=float(params["impact_b"]),
                            baseline=ap_baseline,
                            limb_darkening=BATMAN_LIMB_DARKENING,
                            supersample_factor=BATMAN_SUPERSAMPLE_FACTOR,
                        )
                    else:
                        injected_flux = np.asarray(ap_flux, dtype=np.float64) + float(ap_baseline) * (transit_model - 1.0)
                        ap_mask = mask
                    injected_by_aperture[ap] = injected_flux
                    baseline_by_aperture[ap] = ap_baseline
                    if mask is None:
                        mask = ap_mask
            except Exception:
                skipped["model_failed"] += 1
                continue
            if mask is None:
                skipped["no_good_flux"] += 1
                continue
            sampled_model_depth = float(max(0.0, 1.0 - np.nanmin(transit_model))) if transit_model is not None else float("nan")
            duration_d = float(params["duration_min"]) / 1440.0
            dense_time = np.linspace(t0_d - 2.0 * duration_d, t0_d + 2.0 * duration_d, 2000)
            dense_model = batman_transit_model(
                dense_time,
                period_d=float(params["period_d"]),
                t0_d=t0_d,
                radius_rstar=float(params["radius_rwd"]),
                a_over_rstar=float(params["a_over_rwd"]),
                impact_b=float(params["impact_b"]),
                limb_darkening=BATMAN_LIMB_DARKENING,
                supersample_factor=BATMAN_SUPERSAMPLE_FACTOR,
            )
            model_depth = float(max(0.0, 1.0 - np.nanmin(dense_model)))
            params["model_depth"] = model_depth
            params["sampled_model_depth"] = sampled_model_depth

            injection_id = f"inj_{len(rows):06d}"
            group = inj_group.create_group(injection_id)
            tic = int(target.attrs["tic"])
            group.attrs["injection_id"] = injection_id
            group.attrs["tic"] = tic
            group.attrs["sector"] = int(target.attrs["sector"])
            group.attrs["camera"] = int(target.attrs["camera"])
            group.attrs["ccd"] = int(target.attrs["ccd"])
            group.attrs["tessmag"] = float(target.attrs["tessmag"])
            group.attrs["aperture"] = primary_aperture
            group.attrs["apertures"] = json.dumps(list(apertures))
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
            for attr_key in (
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
            ):
                if params.get(attr_key, "") != "":
                    group.attrs[attr_key] = params[attr_key]
            group.attrs["limb_dark_u1"] = float(BATMAN_LIMB_DARKENING[0])
            group.attrs["limb_dark_u2"] = float(BATMAN_LIMB_DARKENING[1])
            group.attrs["batman_supersample_factor"] = int(BATMAN_SUPERSAMPLE_FACTOR)
            group.attrs["baseline"] = baseline
            for ap, ap_baseline in baseline_by_aperture.items():
                group.attrs[f"baseline_{ap}"] = float(ap_baseline)
            group.attrs["label"] = "planet_like"
            group.attrs["label_source"] = "injection"
            group.attrs["split"] = _split_for_tic(tic)
            group.attrs["source_target_group"] = f"/targets/{target_key}"

            _dataset(group, "time", time.astype(np.float64), compression)
            _dataset(group, "quality", quality.astype(np.int32), compression)
            _dataset(group, "orbitid", orbitid.astype(np.int16), compression)
            _dataset(group, "flux_injected", injected_by_aperture[primary_aperture].astype(np.float32), compression)
            _dataset(group, "in_transit", mask.astype(np.bool_), compression)
            if transit_model is not None:
                _dataset(group, "transit_model", transit_model.astype(np.float32), compression)
            for ap in apertures:
                _dataset(group, f"{ap}_injected", injected_by_aperture[ap].astype(np.float32), compression)
                if store_original:
                    _dataset(group, f"{ap}_original", flux_by_aperture[ap].astype(np.float32), compression)
            if store_original:
                _dataset(group, "flux_original", flux_by_aperture[primary_aperture].astype(np.float32), compression)

            row = {
                "injection_id": injection_id,
                "tic": tic,
                "sector": int(target.attrs["sector"]),
                "camera": int(target.attrs["camera"]),
                "ccd": int(target.attrs["ccd"]),
                "tessmag": float(target.attrs["tessmag"]),
                "aperture": primary_aperture,
                "apertures": "|".join(apertures),
                "signal_family": str(params["signal_family"]),
                "period_d": float(params["period_d"]),
                "t0_d": float(t0_d),
                "t0_bjd": float(t0_d + BJDREFI),
                "duration_min": float(params["duration_min"]),
                "depth": float(params["depth"]),
                "radius_rearth": params.get("radius_rearth", ""),
                "radius_rwd": params.get("radius_rwd", ""),
                "impact_b": params.get("impact_b", ""),
                "a_over_rwd": params.get("a_over_rwd", ""),
                "inclination_deg": params.get("inclination_deg", ""),
                "target_depth": params.get("target_depth", ""),
                "geometric_depth": params.get("geometric_depth", ""),
                "model_depth": params.get("model_depth", ""),
                "sampled_model_depth": params.get("sampled_model_depth", ""),
                "duration_model": params.get("duration_model", ""),
                "injection_model": params.get("injection_model", ""),
                "sampling_mode": params.get("sampling_mode", ""),
                "grid_period_bin": params.get("grid_period_bin", ""),
                "grid_radius_bin": params.get("grid_radius_bin", ""),
                "grid_depth_bin": params.get("grid_depth_bin", ""),
                "grid_cell_id": params.get("grid_cell_id", ""),
                "period_bin_lo_d": params.get("period_bin_lo_d", ""),
                "period_bin_hi_d": params.get("period_bin_hi_d", ""),
                "radius_bin_lo_rearth": params.get("radius_bin_lo_rearth", ""),
                "radius_bin_hi_rearth": params.get("radius_bin_hi_rearth", ""),
                "depth_bin_lo": params.get("depth_bin_lo", ""),
                "depth_bin_hi": params.get("depth_bin_hi", ""),
                "baseline": baseline,
                "n_cadences": int(len(time)),
                "n_quality_zero": int(np.count_nonzero(quality == 0)),
                "n_in_transit": int(np.count_nonzero(mask)),
                "n_good_in_transit": n_good_in_transit,
                "split": _split_for_tic(tic),
                "label": "planet_like",
                "label_source": "injection",
                "source_export_h5": str(export_h5),
                "source_target_group": f"/targets/{target_key}",
                "h5_group": f"/injections/{injection_id}",
            }
            rows.append(row)
            if progress_every > 0 and len(rows) % progress_every == 0:
                print(
                    f"  [make-injections] accepted {len(rows):,}/{n_injections:,} "
                    f"after {total_attempts:,} attempts",
                    flush=True,
                )

        out.attrs["n_injections"] = int(len(rows))
        out.attrs["total_attempts"] = int(total_attempts)

    if len(rows) < n_injections:
        tmp_h5.unlink(missing_ok=True)
        raise RuntimeError(
            f"created {len(rows)} / {n_injections} injections after "
            f"{total_attempts} attempts; skipped={skipped}"
        )

    os.replace(tmp_h5, out_h5)
    manifest = pd.DataFrame(rows)
    manifest.to_csv(manifest_csv, index=False)
    label_cols = [
        "injection_id",
        "tic",
        "split",
        "label",
        "label_source",
        "signal_family",
        "aperture",
        "apertures",
        "period_d",
        "t0_bjd",
        "duration_min",
        "depth",
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
        "n_good_in_transit",
        "h5_group",
    ]
    manifest.loc[:, label_cols].to_csv(labels_csv, index=False)
    summary = {
        "created_utc": created_utc,
        "source_export_h5": str(export_h5),
        "out_h5": str(out_h5),
        "manifest_csv": str(manifest_csv),
        "labels_csv": str(labels_csv),
        "aperture": primary_aperture,
        "apertures": list(apertures),
        "n_injections": int(len(rows)),
        "total_attempts": int(total_attempts),
        "skipped": skipped,
        "sampling_mode": sampling_mode,
        "grid_period_range_d": list(grid_period_range_d),
        "grid_radius_range_rearth": list(grid_radius_range_rearth),
        "grid_depth_range": list(grid_depth_range),
        "grid_period_bins": int(grid_period_bins),
        "grid_radius_bins": int(grid_radius_bins),
        "grid_depth_bins": int(grid_depth_bins),
        "grid_depth_spacing": grid_depth_spacing,
        "injection_model": "batman_quadratic",
        "limb_darkening": list(BATMAN_LIMB_DARKENING),
        "batman_supersample_factor": int(BATMAN_SUPERSAMPLE_FACTOR),
        "signal_family_counts": manifest["signal_family"].value_counts().to_dict(),
        "split_counts": manifest["split"].value_counts().to_dict(),
    }
    if "grid_cell_id" in manifest:
        summary["grid_cell_counts"] = manifest["grid_cell_id"].value_counts().sort_index().to_dict()
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--export-h5", type=Path, default=DEFAULT_EXPORT_H5)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument(
        "--apertures",
        nargs="+",
        default=list(APERTURES),
        help="Aperture columns to inject with the same signal. First aperture is the legacy primary.",
    )
    ap.add_argument(
        "--aperture",
        default=None,
        help="Deprecated single-aperture alias. Use --apertures for multi-channel injections.",
    )
    ap.add_argument("--n-injections", type=int, default=300)
    ap.add_argument("--families", nargs="+", choices=SIGNAL_FAMILIES,
                    default=list(SIGNAL_FAMILIES))
    ap.add_argument(
        "--sampling-mode",
        choices=("period_depth_grid", "period_radius_grid", "family_cycle"),
        default="period_depth_grid",
        help="Use an approximately even grid for injection-recovery maps.",
    )
    ap.add_argument("--grid-period-min-d", type=float, default=DEFAULT_GRID_PERIOD_RANGE_D[0])
    ap.add_argument("--grid-period-max-d", type=float, default=DEFAULT_GRID_PERIOD_RANGE_D[1])
    ap.add_argument("--grid-radius-min-rearth", type=float, default=DEFAULT_GRID_RADIUS_RANGE_REARTH[0])
    ap.add_argument("--grid-radius-max-rearth", type=float, default=DEFAULT_GRID_RADIUS_RANGE_REARTH[1])
    ap.add_argument("--grid-depth-min", type=float, default=DEFAULT_GRID_DEPTH_RANGE[0])
    ap.add_argument("--grid-depth-max", type=float, default=DEFAULT_GRID_DEPTH_RANGE[1])
    ap.add_argument("--grid-period-bins", type=int, default=DEFAULT_GRID_PERIOD_BINS)
    ap.add_argument("--grid-radius-bins", type=int, default=DEFAULT_GRID_RADIUS_BINS)
    ap.add_argument("--grid-depth-bins", type=int, default=DEFAULT_GRID_DEPTH_BINS)
    ap.add_argument("--grid-depth-spacing", choices=("linear", "log"), default="linear")
    ap.add_argument("--random-state", type=int, default=56)
    ap.add_argument("--min-in-transit", type=int, default=2)
    ap.add_argument("--max-attempts-per-injection", type=int, default=25)
    ap.add_argument("--compression", choices=("lzf", "gzip", "none"), default="lzf")
    ap.add_argument("--store-original", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--progress-every", type=int, default=100)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    compression = None if args.compression == "none" else args.compression
    apertures = tuple([args.aperture] if args.aperture else args.apertures)
    summary = _make_injections(
        export_h5=args.export_h5,
        out_dir=args.out_dir,
        apertures=apertures,
        n_injections=args.n_injections,
        families=tuple(args.families),
        sampling_mode=args.sampling_mode,
        grid_period_range_d=(args.grid_period_min_d, args.grid_period_max_d),
        grid_radius_range_rearth=(args.grid_radius_min_rearth, args.grid_radius_max_rearth),
        grid_depth_range=(args.grid_depth_min, args.grid_depth_max),
        grid_period_bins=args.grid_period_bins,
        grid_radius_bins=args.grid_radius_bins,
        grid_depth_bins=args.grid_depth_bins,
        grid_depth_spacing=args.grid_depth_spacing,
        random_state=args.random_state,
        min_in_transit=args.min_in_transit,
        max_attempts_per_injection=args.max_attempts_per_injection,
        compression=compression,
        store_original=args.store_original,
        overwrite=args.overwrite,
        progress_every=args.progress_every,
    )
    print("[make-injections] complete")
    print(f"  injections: {summary['n_injections']:,}")
    print(f"  apertures: {summary['apertures']}")
    print(f"  sampling: {summary['sampling_mode']}")
    print(f"  signal families: {summary['signal_family_counts']}")
    print(f"  splits: {summary['split_counts']}")
    print(f"  out: {summary['out_h5']}")
    print(f"  manifest: {summary['manifest_csv']}")
    print(f"  labels: {summary['labels_csv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
