#!/usr/bin/env python3
"""Sweep flux-space detrending parameters against injected-event preservation.

This experiment asks a narrower question than the model-comparison script:
within the chosen robust-auto subtractive detrend family, how flexible can the
cotrend spline be before it starts absorbing transit/eclipse-like signals?

Outputs:

* config_case_metrics.csv      one row per light curve aperture and config
* injection_grid_summary.csv   aggregate retention by config x duration x depth
* config_summary.csv           aggregate trend/noise and worst retention by config
* summary.md                   compact readable report
* plots/*.png, *.pdf           retention and residual-trend diagnostics
"""
from __future__ import annotations

import argparse
import glob
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.flux_detrend import (  # noqa: E402
    FluxDetrendConfig,
    flux_space_detrend_result,
)
from twirl.lightcurves.tglc_h5_reader import APERTURE_KEYS, read_tglc_h5  # noqa: E402
from twirl.plotting.style import apply_twirl_style, get_ordered_palette  # noqa: E402


@dataclass(frozen=True)
class CurveCase:
    case_id: str
    tic: int
    sector: int
    orbit: int
    cam: int
    ccd: int
    tmag: float
    aperture: str
    time: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray | None
    quality: np.ndarray
    path: str


def robust_sigma(x: np.ndarray) -> float:
    finite = np.asarray(x, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    med = np.nanmedian(finite)
    mad = np.nanmedian(np.abs(finite - med))
    sigma = 1.4826 * float(mad)
    if np.isfinite(sigma) and sigma > 0:
        return sigma
    sigma = float(np.nanstd(finite))
    return sigma if np.isfinite(sigma) else np.nan


def robust_abs_scale(x: np.ndarray) -> float:
    finite = np.asarray(x, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    abs_flux = np.abs(finite)
    scale = float(np.nanpercentile(abs_flux, 68.0))
    if np.isfinite(scale) and scale > 0:
        return scale
    scale = float(np.nanmedian(abs_flux))
    return scale if np.isfinite(scale) else np.nan


def binned_trend_metrics(time: np.ndarray, flux: np.ndarray, quality: np.ndarray, bin_d: float) -> tuple[float, float, int]:
    good = np.isfinite(time) & np.isfinite(flux) & (quality == 0)
    if good.sum() < 100:
        return np.nan, np.nan, 0
    t = time[good]
    y = flux[good]
    start = float(np.nanmin(t))
    end = float(np.nanmax(t))
    edges = np.arange(start, end + bin_d, bin_d)
    if edges.size < 3:
        return np.nan, np.nan, 0
    medians: list[float] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = (t >= lo) & (t < hi)
        if in_bin.sum() >= 20:
            medians.append(float(np.nanmedian(y[in_bin])))
    if len(medians) < 3:
        return np.nan, np.nan, len(medians)
    arr = np.asarray(medians, dtype=np.float64)
    centered = arr - np.nanmedian(arr)
    return robust_sigma(centered), float(np.nanpercentile(arr, 95) - np.nanpercentile(arr, 5)), len(arr)


def injection_retention(
    case: CurveCase,
    cfg: FluxDetrendConfig,
    base_det: np.ndarray,
    *,
    depth: float,
    duration_min: float,
    period_d: float,
) -> dict:
    q0 = np.asarray(case.quality) == 0
    finite = np.isfinite(case.time) & np.isfinite(case.flux)
    good = finite & q0
    if good.sum() < 50:
        return {}
    flux_scale = robust_abs_scale(case.flux[good])
    if not np.isfinite(flux_scale) or flux_scale <= 0:
        return {}

    t0 = float(np.nanmin(case.time[finite]) + 0.37 * period_d)
    duration_d = duration_min / 60.0 / 24.0
    phase = ((case.time - t0 + 0.5 * period_d) % period_d) - 0.5 * period_d
    in_event = finite & q0 & (np.abs(phase) <= 0.5 * duration_d)
    out_event = finite & q0 & ~in_event
    if in_event.sum() < 1 or out_event.sum() < 50:
        return {}

    injected = case.flux.copy()
    injected[in_event] = injected[in_event] - depth * flux_scale
    inj = flux_space_detrend_result(
        case.time,
        injected,
        quality=case.quality,
        flux_err=case.flux_err,
        cfg=cfg,
    )
    delta = inj.det_flux - base_det
    if not np.isfinite(delta[in_event]).any():
        return {}
    local_offset = np.nanmedian(delta[out_event]) if np.isfinite(delta[out_event]).any() else 0.0
    recovered_depth = -(float(np.nanmedian(delta[in_event])) - float(local_offset))
    return {
        "depth": float(depth),
        "duration_min": float(duration_min),
        "period_d": float(period_d),
        "injection_flux_scale": float(flux_scale),
        "injection_cadences": int(in_event.sum()),
        "recovered_depth": recovered_depth,
        "depth_retention_frac": recovered_depth / depth if depth else np.nan,
    }


def collect_h5_paths(args: argparse.Namespace) -> list[Path]:
    paths: list[Path] = [Path(p) for p in args.h5]
    for pat in args.h5_glob:
        paths.extend(Path(p) for p in sorted(glob.glob(pat)))
    out: list[Path] = []
    seen: set[Path] = set()
    for p in paths:
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def parse_apertures(values: list[str]) -> tuple[str, ...]:
    apertures: list[str] = []
    for value in values:
        if value.lower() == "all":
            apertures.extend(APERTURE_KEYS)
        else:
            apertures.append(value)
    bad = [a for a in apertures if a not in APERTURE_KEYS]
    if bad:
        raise SystemExit(f"unknown aperture(s): {', '.join(bad)}; valid: {', '.join(APERTURE_KEYS)}")
    return tuple(dict.fromkeys(apertures))


def load_cases(paths: list[Path], apertures: tuple[str, ...], limit: int | None) -> list[CurveCase]:
    cases: list[CurveCase] = []
    n_files = 0
    for path in paths:
        try:
            lc = read_tglc_h5(path)
        except Exception as exc:
            print(f"[warn] failed to read {path}: {type(exc).__name__}: {exc}", file=sys.stderr)
            continue
        n_files += 1
        for ap in apertures:
            if ap not in lc.apertures:
                continue
            ap_data = lc.apertures[ap]
            cases.append(
                CurveCase(
                    case_id=f"tic{lc.tic}_s{lc.sector:04d}_o{lc.orbit}_c{lc.cam}-{lc.ccd}_{ap}",
                    tic=lc.tic,
                    sector=lc.sector,
                    orbit=lc.orbit,
                    cam=lc.cam,
                    ccd=lc.ccd,
                    tmag=lc.tmag,
                    aperture=ap,
                    time=lc.time,
                    flux=ap_data.raw_flux,
                    flux_err=ap_data.raw_flux_err,
                    quality=lc.quality,
                    path=str(path),
                )
            )
        if limit is not None and n_files >= limit:
            break
    return cases


def config_label(bkspace_d: float, sigma_clip: float) -> str:
    return f"bk{bkspace_d:g}_sig{sigma_clip:g}"


def run_sweep(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    apertures = parse_apertures(args.aperture)
    h5_paths = collect_h5_paths(args)
    cases = load_cases(h5_paths, apertures, args.limit)
    if not cases:
        raise SystemExit("[sweep-detrend] no cases selected")

    metric_rows: list[dict] = []
    injection_rows: list[dict] = []
    configs = [
        (bkspace, sigma_clip, FluxDetrendConfig(
            bkspace_d=bkspace,
            sigma_clip=sigma_clip,
            output_mode="subtractive",
            scale_strategy="auto",
        ))
        for bkspace in args.bkspace_d
        for sigma_clip in args.sigma_clip
    ]
    total = len(cases) * len(configs)
    i = 0
    for case in cases:
        for bkspace, sigma_clip, cfg in configs:
            i += 1
            label = config_label(bkspace, sigma_clip)
            print(f"[sweep-detrend] {i}/{total} {case.case_id} {label}")
            result = flux_space_detrend_result(
                case.time,
                case.flux,
                quality=case.quality,
                flux_err=case.flux_err,
                cfg=cfg,
            )
            q0 = case.quality == 0
            raw_negative_q0 = np.isfinite(case.flux) & (case.flux < 0) & q0
            good_det = np.isfinite(result.det_flux) & q0
            det_good = result.det_flux[good_det]
            det_mad = robust_sigma(det_good)
            trend_mad, trend_range, n_bins = binned_trend_metrics(
                case.time,
                result.det_flux,
                case.quality,
                args.trend_bin_d,
            )
            failed = (
                good_det.sum() < 50
                or not np.isfinite(det_mad)
                or det_mad > args.max_reasonable_mad
            )
            metric_rows.append({
                "config": label,
                "bkspace_d": bkspace,
                "sigma_clip": sigma_clip,
                "case_id": case.case_id,
                "tic": case.tic,
                "tmag": case.tmag,
                "aperture": case.aperture,
                "scale": result.scale,
                "scale_source": result.scale_source,
                "fit_count": result.fit_count,
                "n_quality0": int(q0.sum()),
                "n_det_finite_quality0": int(good_det.sum()),
                "raw_negative_quality0": int(raw_negative_q0.sum()),
                "raw_negative_quality0_retained": int((raw_negative_q0 & np.isfinite(result.det_flux)).sum()),
                "det_mad_quality0": det_mad,
                "trend_mad_binned": trend_mad,
                "trend_range_binned": trend_range,
                "trend_n_bins": n_bins,
                "failed": bool(failed),
                "path": case.path,
            })
            for duration_min in args.duration_min:
                for depth in args.depth:
                    inj = injection_retention(
                        case,
                        cfg,
                        result.det_flux,
                        depth=depth,
                        duration_min=duration_min,
                        period_d=args.period_d,
                    )
                    if inj:
                        inj.update({
                            "config": label,
                            "bkspace_d": bkspace,
                            "sigma_clip": sigma_clip,
                            "case_id": case.case_id,
                            "tic": case.tic,
                            "tmag": case.tmag,
                            "aperture": case.aperture,
                        })
                        injection_rows.append(inj)
    return pd.DataFrame(metric_rows), pd.DataFrame(injection_rows)


def summarize(metric_df: pd.DataFrame, inj_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    inj_summary = (
        inj_df.groupby(["config", "bkspace_d", "sigma_clip", "duration_min", "depth"], dropna=False)
        .agg(
            injection_cases=("case_id", "nunique"),
            median_retention=("depth_retention_frac", "median"),
            p16_retention=("depth_retention_frac", lambda x: np.nanpercentile(x, 16)),
            p84_retention=("depth_retention_frac", lambda x: np.nanpercentile(x, 84)),
            min_retention=("depth_retention_frac", "min"),
            max_retention=("depth_retention_frac", "max"),
            median_cadences=("injection_cadences", "median"),
        )
        .reset_index()
    )
    worst = (
        inj_summary.groupby(["config", "bkspace_d", "sigma_clip"], dropna=False)
        .agg(
            worst_p16_retention=("p16_retention", "min"),
            worst_median_retention=("median_retention", "min"),
            max_p84_retention=("p84_retention", "max"),
        )
        .reset_index()
    )
    base = (
        metric_df.groupby(["config", "bkspace_d", "sigma_clip"], dropna=False)
        .agg(
            cases=("case_id", "nunique"),
            failures=("failed", "sum"),
            median_det_mad=("det_mad_quality0", "median"),
            median_trend_mad=("trend_mad_binned", "median"),
            median_trend_range=("trend_range_binned", "median"),
            raw_neg_q0=("raw_negative_quality0", "sum"),
            retained_neg_q0=("raw_negative_quality0_retained", "sum"),
        )
        .reset_index()
    )
    config_summary = base.merge(worst, on=["config", "bkspace_d", "sigma_clip"], how="left")
    config_summary["retained_neg_q0_frac"] = (
        config_summary["retained_neg_q0"] / config_summary["raw_neg_q0"].replace(0, np.nan)
    )
    config_summary["pass_retention_90_110"] = (
        (config_summary["worst_p16_retention"] >= 0.90)
        & (config_summary["max_p84_retention"] <= 1.10)
    )
    return config_summary, inj_summary


def frame_to_markdown(df: pd.DataFrame, *, floatfmt: str = ".5g") -> str:
    if df.empty:
        return "_No rows._"
    headers = [str(c) for c in df.columns]
    body: list[list[str]] = []
    for _, row in df.iterrows():
        values: list[str] = []
        for value in row:
            if isinstance(value, float):
                values.append(format(value, floatfmt) if np.isfinite(value) else "")
            else:
                values.append(str(value))
        body.append(values)
    widths = [max(len(headers[i]), *(len(row[i]) for row in body)) for i in range(len(headers))]
    lines = [
        "| " + " | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    for row in body:
        lines.append("| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(headers))) + " |")
    return "\n".join(lines)


def write_summary(
    out_path: Path,
    config_summary: pd.DataFrame,
    inj_summary: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    cols = [
        "config",
        "bkspace_d",
        "sigma_clip",
        "failures",
        "median_det_mad",
        "median_trend_mad",
        "worst_p16_retention",
        "worst_median_retention",
        "max_p84_retention",
        "pass_retention_90_110",
    ]
    ranked = config_summary.sort_values(
        ["pass_retention_90_110", "median_trend_mad", "median_det_mad"],
        ascending=[False, True, True],
    )
    lines = [
        "# Flux Detrend Parameter Sweep",
        "",
        f"Durations tested: `{', '.join(str(x) for x in args.duration_min)} min`.",
        f"Depths tested: `{', '.join(str(x) for x in args.depth)}`.",
        f"Trend metric: MAD of `{args.trend_bin_d:g} d` binned detrended medians.",
        "",
        "## Config Summary",
        "",
        frame_to_markdown(ranked[cols], floatfmt=".5g"),
        "",
        "## Worst Retention Rows",
        "",
    ]
    worst_rows = inj_summary.sort_values(["p16_retention", "median_retention"]).head(20)
    inj_cols = [
        "config",
        "duration_min",
        "depth",
        "median_retention",
        "p16_retention",
        "p84_retention",
        "median_cadences",
    ]
    lines.append(frame_to_markdown(worst_rows[inj_cols], floatfmt=".5g"))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_outputs(config_summary: pd.DataFrame, inj_summary: pd.DataFrame, out_dir: Path) -> None:
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    apply_twirl_style("full_page")

    for sigma_clip in sorted(inj_summary["sigma_clip"].dropna().unique()):
        sigma_df = inj_summary[inj_summary["sigma_clip"] == sigma_clip]
        if sigma_df.empty:
            continue
        depth_for_plot = min(sigma_df["depth"].unique(), key=lambda x: abs(x - 0.10))
        ddf = sigma_df[sigma_df["depth"] == depth_for_plot]
        fig, ax = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)
        colors = get_ordered_palette(ddf["bkspace_d"].nunique(), "viridis")
        for color, (bkspace, g) in zip(colors, ddf.groupby("bkspace_d")):
            g = g.sort_values("duration_min")
            ax.plot(g["duration_min"], g["median_retention"], marker="o", lw=1.4, color=color, label=f"{bkspace:g} d")
            ax.fill_between(g["duration_min"], g["p16_retention"], g["p84_retention"], color=color, alpha=0.12, lw=0)
        ax.axhline(1.0, color="black", lw=0.8, ls=":")
        ax.axhline(0.9, color="0.4", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_xlabel("Injected duration (min)")
        ax.set_ylabel("Depth retention fraction")
        ax.legend(title=f"bkspace, sigma={sigma_clip:g}", loc="lower left", ncol=2)
        stem = f"retention_depth{depth_for_plot:g}_sigma{sigma_clip:g}".replace(".", "p")
        fig.savefig(plot_dir / f"{stem}.png", dpi=180)
        fig.savefig(plot_dir / f"{stem}.pdf")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.0), constrained_layout=True)
    for sigma_clip, g in config_summary.groupby("sigma_clip"):
        g = g.sort_values("bkspace_d")
        ax.plot(g["bkspace_d"], g["median_trend_mad"], marker="o", lw=1.4, label=f"sigma={sigma_clip:g}")
    ax.set_xlabel("Spline knot spacing (d)")
    ax.set_ylabel("Median binned-trend MAD")
    ax.legend(loc="best")
    fig.savefig(plot_dir / "trend_metric_vs_bkspace.png", dpi=180)
    fig.savefig(plot_dir / "trend_metric_vs_bkspace.pdf")
    plt.close(fig)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5", action="append", default=[], help="TGLC HDF5 file; repeatable.")
    ap.add_argument("--h5-glob", action="append", default=[], help="Quoted glob for TGLC HDF5 files.")
    ap.add_argument("--aperture", action="append", default=["Primary"], help="'Small', 'Primary', 'Large', or 'all'.")
    ap.add_argument("--limit", type=int, default=None, help="Maximum HDF5 files to read.")
    ap.add_argument("--bkspace-d", type=float, nargs="+", default=[0.3, 0.5, 0.8, 1.2, 2.0])
    ap.add_argument("--sigma-clip", type=float, nargs="+", default=[3.0, 5.0])
    ap.add_argument("--duration-min", type=float, nargs="+", default=[5.0, 10.0, 30.0, 60.0, 120.0, 180.0])
    ap.add_argument("--depth", type=float, nargs="+", default=[0.03, 0.10, 0.30, 0.50])
    ap.add_argument("--period-d", type=float, default=1.40793903)
    ap.add_argument("--trend-bin-d", type=float, default=0.5)
    ap.add_argument("--max-reasonable-mad", type=float, default=10.0)
    ap.add_argument("--out-dir", type=Path, default=Path("reports/stage1_lightcurves/detrend_experiments/param_sweep"))
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    metric_df, inj_df = run_sweep(args)
    config_summary, inj_summary = summarize(metric_df, inj_df)

    metric_path = args.out_dir / "config_case_metrics.csv"
    inj_path = args.out_dir / "injection_grid_summary.csv"
    config_path = args.out_dir / "config_summary.csv"
    metric_df.to_csv(metric_path, index=False)
    inj_summary.to_csv(inj_path, index=False)
    config_summary.to_csv(config_path, index=False)
    write_summary(args.out_dir / "summary.md", config_summary, inj_summary, args)
    plot_outputs(config_summary, inj_summary, args.out_dir)
    print(f"[sweep-detrend] wrote {metric_path}")
    print(f"[sweep-detrend] wrote {inj_path}")
    print(f"[sweep-detrend] wrote {config_path}")
    print(f"[sweep-detrend] wrote {args.out_dir / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
