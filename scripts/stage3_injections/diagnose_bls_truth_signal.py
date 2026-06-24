#!/usr/bin/env python3
"""Diagnose BLS recovery at the known truth ephemeris for one injected HLSP.

This is meant for injection-recovery debugging, especially cases where BLS
selects a broad wrong-period peak. It reports both the global BLS winner and
the nearest truth-period power/rank under the same periodogram.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
from astropy.timeseries import BoxLeastSquares

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import BJDREFI, quality_mask, read_hlsp  # noqa: E402
from twirl.search.bls import BLSConfig  # noqa: E402
from twirl.search.candidates import compute_sde  # noqa: E402
from twirl.search.grids import build_period_grid, duration_grid_days  # noqa: E402
from twirl.search.injections import box_transit_mask  # noqa: E402


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _resolve_repo_path(value: str | Path | None) -> Path | None:
    if value is None or str(value) == "":
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def _manifest_truth(manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(Path(manifest_path).read_text())
    return {
        "manifest": str(manifest_path),
        "injected_hlsp": _resolve_repo_path(manifest.get("twirlfs_injected")),
        "original_hlsp": _resolve_repo_path(manifest.get("twirlfs_original")),
        "truth_period_d": float(manifest["period_d"]),
        "truth_t0_bjd": float(manifest["t0_tjd"]) + float(BJDREFI),
        "truth_duration_min": float(manifest["duration_min"]),
        "n_good_in_transit_manifest": manifest.get("n_good_in_transit", ""),
        "model_depth": manifest.get("model_depth", ""),
        "target": manifest.get("target", {}),
    }


def _phase_error_min(t0_bjd: float, truth_t0_bjd: float, period_d: float) -> float:
    if not all(np.isfinite([t0_bjd, truth_t0_bjd, period_d])) or period_d <= 0:
        return float("nan")
    delta_d = ((t0_bjd - truth_t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d
    return float(delta_d * 1440.0)


def _truth_window_stats(
    time: np.ndarray,
    flux: np.ndarray,
    *,
    truth_period_d: float,
    truth_t0_bjd: float,
    truth_duration_min: float,
    quality_mask_values: np.ndarray,
) -> dict[str, Any]:
    truth_t0_tjd = float(truth_t0_bjd) - float(BJDREFI)
    in_tran = box_transit_mask(
        time,
        period_d=truth_period_d,
        t0_d=truth_t0_tjd,
        duration_min=truth_duration_min,
    )
    near = box_transit_mask(
        time,
        period_d=truth_period_d,
        t0_d=truth_t0_tjd,
        duration_min=max(5.0 * truth_duration_min, truth_duration_min + 20.0),
    )
    good = quality_mask_values & np.isfinite(flux)
    in_good = good & in_tran
    oot_good = good & ~near
    if np.any(in_good) and np.any(oot_good):
        in_med = float(np.nanmedian(flux[in_good]))
        oot_med = float(np.nanmedian(flux[oot_good]))
        depth = float((oot_med - in_med) / oot_med) if oot_med != 0 and np.isfinite(oot_med) else float("nan")
        scatter = float(1.4826 * np.nanmedian(np.abs(flux[oot_good] - oot_med)))
        snr = float((oot_med - in_med) / scatter * np.sqrt(np.count_nonzero(in_good))) if scatter > 0 else float("nan")
    else:
        in_med = oot_med = depth = scatter = snr = float("nan")
    return {
        "n_good": int(np.count_nonzero(good)),
        "n_truth_in_transit": int(np.count_nonzero(in_tran)),
        "n_good_truth_in_transit": int(np.count_nonzero(in_good)),
        "n_good_truth_oot": int(np.count_nonzero(oot_good)),
        "truth_in_median": in_med,
        "truth_oot_median": oot_med,
        "truth_window_depth": depth,
        "truth_oot_mad_sigma": scatter,
        "truth_window_snr_estimate": snr,
    }


def _diagnose_aperture(
    *,
    hlsp_path: Path,
    aperture: str,
    truth_period_d: float,
    truth_t0_bjd: float,
    truth_duration_min: float,
    durations_min: tuple[float, ...],
    n_periods: int,
    sigma_clip: float,
    p_min_d: float,
    max_period_fraction: float,
    p_max_cap_d: float,
    original_hlsp_path: Path | None,
) -> dict[str, Any]:
    columns = (aperture,)
    lc = read_hlsp(hlsp_path, columns=columns)
    if lc is None or aperture not in lc.flux:
        return {"status": "missing_or_unreadable", "aperture": aperture, "hlsp_path": str(hlsp_path)}
    cfg = BLSConfig(
        apertures=(aperture,),
        durations_min=durations_min,
        n_periods=int(n_periods),
        sigma_clip=float(sigma_clip),
        p_min_d=float(p_min_d),
        max_period_fraction=float(max_period_fraction),
        p_max_cap_d=float(p_max_cap_d),
    )
    mask = quality_mask(lc, aperture)
    flux = np.asarray(lc.flux[aperture], dtype=np.float64)
    truth_stats = _truth_window_stats(
        lc.time,
        flux,
        truth_period_d=truth_period_d,
        truth_t0_bjd=truth_t0_bjd,
        truth_duration_min=truth_duration_min,
        quality_mask_values=mask,
    )
    t = np.asarray(lc.time[mask], dtype=np.float64)
    f = np.asarray(flux[mask], dtype=np.float64)
    if len(t) < cfg.min_cadences:
        return {
            "status": "too_few_cadences",
            "aperture": aperture,
            "hlsp_path": str(hlsp_path),
            "truth_window": truth_stats,
        }
    med = float(np.nanmedian(f))
    y = f / med
    n_after_clip = len(y)
    if cfg.sigma_clip and cfg.sigma_clip > 0:
        mad = float(np.nanmedian(np.abs(y - 1.0)))
        if mad > 0 and np.isfinite(mad):
            sigma = 1.4826 * mad
            keep = (y - 1.0) <= cfg.sigma_clip * sigma
            if keep.sum() >= cfg.min_cadences:
                t = t[keep]
                y = y[keep]
                n_after_clip = int(keep.sum())

    durations_d = duration_grid_days(cfg.durations_min)
    periods = build_period_grid(
        baseline_d=float(np.nanmax(t) - np.nanmin(t)),
        p_min_d=cfg.p_min_d,
        max_period_fraction=cfg.max_period_fraction,
        p_max_cap_d=cfg.p_max_cap_d,
        n_periods=cfg.n_periods,
        durations_d=durations_d,
    )
    truth_period = np.array([truth_period_d], dtype=np.float64)
    if not np.any(np.isclose(periods, truth_period_d, rtol=0, atol=1.0e-10)):
        periods = np.sort(np.unique(np.concatenate([periods, truth_period])))
    bls = BoxLeastSquares(t, y)
    pg = bls.power(periods, durations_d, oversample=1)
    power = np.asarray(pg.power, dtype=np.float64)
    sde = compute_sde(power)
    order = np.argsort(sde)[::-1]
    truth_idx = int(np.argmin(np.abs(np.asarray(pg.period, dtype=np.float64) - truth_period_d)))
    truth_rank = int(np.flatnonzero(order == truth_idx)[0] + 1)
    best_idx = int(order[0])
    top = []
    for idx in order[:10]:
        top.append(
            {
                "period_d": float(pg.period[idx]),
                "duration_min": float(pg.duration[idx] * 1440.0),
                "t0_bjd": float(pg.transit_time[idx] + BJDREFI),
                "depth": float(pg.depth[idx]),
                "depth_snr": float(pg.depth_snr[idx]),
                "sde": float(sde[idx]),
                "t0_phase_err_min": _phase_error_min(
                    float(pg.transit_time[idx] + BJDREFI),
                    truth_t0_bjd,
                    truth_period_d,
                ),
            }
        )

    original_delta = {}
    if original_hlsp_path is not None and original_hlsp_path.exists():
        orig = read_hlsp(original_hlsp_path, columns=columns)
        if orig is not None and aperture in orig.flux and len(orig.time) == len(lc.time):
            delta = np.asarray(flux, dtype=np.float64) - np.asarray(orig.flux[aperture], dtype=np.float64)
            original_delta = _truth_window_stats(
                lc.time,
                delta,
                truth_period_d=truth_period_d,
                truth_t0_bjd=truth_t0_bjd,
                truth_duration_min=truth_duration_min,
                quality_mask_values=mask,
            )

    return {
        "status": "ok",
        "aperture": aperture,
        "hlsp_path": str(hlsp_path),
        "config": asdict(cfg),
        "n_periods_actual": int(len(periods)),
        "n_after_clip": int(n_after_clip),
        "truth_window": truth_stats,
        "injected_minus_original_window": original_delta,
        "best": {
            "period_d": float(pg.period[best_idx]),
            "duration_min": float(pg.duration[best_idx] * 1440.0),
            "t0_bjd": float(pg.transit_time[best_idx] + BJDREFI),
            "depth": float(pg.depth[best_idx]),
            "depth_snr": float(pg.depth_snr[best_idx]),
            "sde": float(sde[best_idx]),
            "t0_phase_err_min": _phase_error_min(
                float(pg.transit_time[best_idx] + BJDREFI),
                truth_t0_bjd,
                truth_period_d,
            ),
        },
        "truth_period": {
            "nearest_period_d": float(pg.period[truth_idx]),
            "duration_min": float(pg.duration[truth_idx] * 1440.0),
            "t0_bjd": float(pg.transit_time[truth_idx] + BJDREFI),
            "depth": float(pg.depth[truth_idx]),
            "depth_snr": float(pg.depth_snr[truth_idx]),
            "sde": float(sde[truth_idx]),
            "rank": truth_rank,
            "t0_phase_err_min": _phase_error_min(
                float(pg.transit_time[truth_idx] + BJDREFI),
                truth_t0_bjd,
                truth_period_d,
            ),
        },
        "top_peaks": top,
    }


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=None, help="Pixel-injection smoke manifest with truth and TWIRL-FS paths.")
    ap.add_argument("--hlsp", type=Path, default=None, help="Injected HLSP FITS path. Overrides manifest twirlfs_injected.")
    ap.add_argument("--original-hlsp", type=Path, default=None, help="Original HLSP FITS for injected-minus-original diagnostics.")
    ap.add_argument("--period-d", type=float, default=None)
    ap.add_argument("--t0-bjd", type=float, default=None)
    ap.add_argument("--duration-min", type=float, default=None)
    ap.add_argument("--apertures", nargs="+", default=["DET_FLUX_ADP", "DET_FLUX"])
    ap.add_argument("--durations-min", nargs="+", type=float, default=[1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 30.0])
    ap.add_argument("--n-periods", type=int, default=20000)
    ap.add_argument("--sigma-clip", type=float, default=5.0)
    ap.add_argument("--p-min-d", type=float, default=0.08)
    ap.add_argument("--max-period-fraction", type=float, default=0.45)
    ap.add_argument("--p-max-cap-d", type=float, default=15.0)
    ap.add_argument("--out-json", type=Path, default=None)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    manifest_info: dict[str, Any] = {}
    if args.manifest is not None:
        manifest_info = _manifest_truth(args.manifest)
    hlsp = args.hlsp or manifest_info.get("injected_hlsp")
    if hlsp is None:
        raise ValueError("provide --hlsp or --manifest with twirlfs_injected")
    original_hlsp = args.original_hlsp or manifest_info.get("original_hlsp")
    truth_period_d = float(args.period_d if args.period_d is not None else manifest_info["truth_period_d"])
    truth_t0_bjd = float(args.t0_bjd if args.t0_bjd is not None else manifest_info["truth_t0_bjd"])
    truth_duration_min = float(args.duration_min if args.duration_min is not None else manifest_info["truth_duration_min"])

    results = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(args.manifest) if args.manifest else "",
        "hlsp": str(hlsp),
        "original_hlsp": str(original_hlsp) if original_hlsp else "",
        "truth": {
            "period_d": truth_period_d,
            "t0_bjd": truth_t0_bjd,
            "duration_min": truth_duration_min,
            "n_good_in_transit_manifest": manifest_info.get("n_good_in_transit_manifest", ""),
            "model_depth": manifest_info.get("model_depth", ""),
            "target": manifest_info.get("target", {}),
        },
        "apertures": {},
    }
    for aperture in args.apertures:
        results["apertures"][aperture] = _diagnose_aperture(
            hlsp_path=Path(hlsp),
            aperture=str(aperture),
            truth_period_d=truth_period_d,
            truth_t0_bjd=truth_t0_bjd,
            truth_duration_min=truth_duration_min,
            durations_min=tuple(float(v) for v in args.durations_min),
            n_periods=int(args.n_periods),
            sigma_clip=float(args.sigma_clip),
            p_min_d=float(args.p_min_d),
            max_period_fraction=float(args.max_period_fraction),
            p_max_cap_d=float(args.p_max_cap_d),
            original_hlsp_path=Path(original_hlsp) if original_hlsp else None,
        )

    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(results, indent=2, sort_keys=True, default=_json_default) + "\n")
        print(f"[truth-bls] wrote {args.out_json}")
    else:
        print(json.dumps(results, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
