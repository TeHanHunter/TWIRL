#!/usr/bin/env python3
"""Compare flux-space detrending variants on faint TGLC light curves.

This is the reproducible Stage 1.6 experiment for dim WD light curves whose
background-subtracted linear flux can be negative or near zero. It accepts real
TGLC HDF5 files from PDO, and also has a synthetic mode so the numerical
failure modes are testable without staging survey products locally.

Example real-data run on PDO:

    python scripts/stage1_lightcurves/compare_flux_detrend_models.py \\
        --h5-glob '/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi/cam*/ccd*/LC/*.h5' \\
        --h5-glob '/pdo/users/tehan/tglc-gpu-production/orbit-120/ffi/cam*/ccd*/LC/*.h5' \\
        --tmag-min 18.5 --limit 50 \\
        --out-dir reports/stage1_lightcurves/detrend_experiments/s56_faint

Local smoke test:

    python scripts/stage1_lightcurves/compare_flux_detrend_models.py \\
        --synthetic 8 \\
        --out-dir reports/stage1_lightcurves/detrend_experiments/synthetic_smoke
"""
from __future__ import annotations

import argparse
import glob
import json
import math
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
    FluxDetrendResult,
    flux_space_detrend_result,
)
from twirl.plotting.style import apply_twirl_style, get_ordered_palette  # noqa: E402


APERTURE_KEYS: tuple[str, ...] = ("Small", "Primary", "Large")

VARIANTS: dict[str, FluxDetrendConfig] = {
    "divisive": FluxDetrendConfig(output_mode="divisive"),
    "sub_median": FluxDetrendConfig(output_mode="subtractive", scale_strategy="median"),
    "sub_median_abs": FluxDetrendConfig(output_mode="subtractive", scale_strategy="median_abs"),
    "sub_auto": FluxDetrendConfig(output_mode="subtractive", scale_strategy="auto"),
}

MAX_REASONABLE_REL_MAD = 10.0


@dataclass
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
    raw_flux_source: str


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


def rms_about_median(x: np.ndarray) -> float:
    finite = np.asarray(x, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    med = np.nanmedian(finite)
    return float(np.sqrt(np.nanmean((finite - med) ** 2)))


def summarize_result(case: CurveCase, variant: str, result: FluxDetrendResult) -> dict:
    q0 = np.asarray(case.quality) == 0
    finite_raw = np.isfinite(case.flux)
    finite_det = np.isfinite(result.det_flux)
    raw_negative = finite_raw & (case.flux < 0)
    raw_negative_q0 = raw_negative & q0
    raw_negative_flagged = raw_negative & ~q0
    retained_negative_q0 = raw_negative_q0 & finite_det
    good_det = q0 & finite_det
    det_good = result.det_flux[good_det]

    if det_good.size:
        p05, p50, p95 = np.nanpercentile(det_good, [5, 50, 95])
    else:
        p05 = p50 = p95 = np.nan
    mad = robust_sigma(det_good)
    rms = rms_about_median(det_good)
    failure_reasons: list[str] = []
    if int(good_det.sum()) < 50:
        failure_reasons.append("too_few_quality0_cadences")
    if not np.isfinite(mad):
        failure_reasons.append("nonfinite_mad")
    if not np.isfinite(p50):
        failure_reasons.append("nonfinite_median")
    if np.isfinite(p50) and abs(float(p50)) > 1.0e6:
        failure_reasons.append("extreme_median")
    if np.isfinite(mad) and mad > MAX_REASONABLE_REL_MAD:
        failure_reasons.append("unstable_relative_mad")
    failed = bool(failure_reasons)

    return {
        "case_id": case.case_id,
        "tic": case.tic,
        "sector": case.sector,
        "orbit": case.orbit,
        "cam": case.cam,
        "ccd": case.ccd,
        "tmag": case.tmag,
        "aperture": case.aperture,
        "variant": variant,
        "output_mode": result.output_mode,
        "scale_strategy": result.scale_strategy,
        "scale": result.scale,
        "scale_source": result.scale_source,
        "fit_count": result.fit_count,
        "n_cadence": int(case.flux.size),
        "n_raw_finite": int(finite_raw.sum()),
        "n_quality0": int(q0.sum()),
        "n_det_finite": int(finite_det.sum()),
        "n_det_finite_quality0": int(good_det.sum()),
        "raw_negative": int(raw_negative.sum()),
        "raw_negative_quality0": int(raw_negative_q0.sum()),
        "raw_negative_flagged": int(raw_negative_flagged.sum()),
        "raw_negative_quality0_retained": int(retained_negative_q0.sum()),
        "raw_negative_quality0_retained_frac": (
            float(retained_negative_q0.sum() / raw_negative_q0.sum())
            if raw_negative_q0.sum()
            else np.nan
        ),
        "det_median_quality0": float(p50) if np.isfinite(p50) else np.nan,
        "det_mad_quality0": float(mad) if np.isfinite(mad) else np.nan,
        "det_rms_quality0": float(rms) if np.isfinite(rms) else np.nan,
        "det_p05_quality0": float(p05) if np.isfinite(p05) else np.nan,
        "det_p95_quality0": float(p95) if np.isfinite(p95) else np.nan,
        "failed": bool(failed),
        "failure_reason": ",".join(failure_reasons),
        "path": case.path,
        "raw_flux_source": case.raw_flux_source,
    }


def injection_preservation(
    case: CurveCase,
    variant: str,
    cfg: FluxDetrendConfig,
    *,
    depth: float,
    duration_min: float,
    period_d: float,
) -> dict:
    q0 = np.asarray(case.quality) == 0
    finite = np.isfinite(case.time) & np.isfinite(case.flux)
    fit = finite & q0
    if fit.sum() < 50 or depth <= 0:
        return {}

    fit_flux = case.flux[fit]
    flux_scale = float(np.nanpercentile(np.abs(fit_flux), 68.0))
    if not np.isfinite(flux_scale) or flux_scale <= 0:
        flux_scale = abs(float(np.nanmedian(fit_flux)))
    if not np.isfinite(flux_scale) or flux_scale <= 0:
        return {}

    t0 = float(np.nanmin(case.time[finite]) + 0.37 * period_d)
    duration_d = duration_min / 60.0 / 24.0
    phase = ((case.time - t0 + 0.5 * period_d) % period_d) - 0.5 * period_d
    in_event = finite & q0 & (np.abs(phase) <= 0.5 * duration_d)
    out_event = finite & q0 & ~in_event
    if in_event.sum() < 2 or out_event.sum() < 50:
        return {}

    injected = case.flux.copy()
    injected[in_event] = injected[in_event] - depth * flux_scale
    base = flux_space_detrend_result(
        case.time, case.flux, quality=case.quality, flux_err=case.flux_err, cfg=cfg
    )
    inj = flux_space_detrend_result(
        case.time, injected, quality=case.quality, flux_err=case.flux_err, cfg=cfg
    )
    delta = inj.det_flux - base.det_flux
    if not np.isfinite(delta[in_event]).any():
        return {}
    local_offset = np.nanmedian(delta[out_event]) if np.isfinite(delta[out_event]).any() else 0.0
    recovered_depth = -(float(np.nanmedian(delta[in_event])) - float(local_offset))
    return {
        "case_id": case.case_id,
        "variant": variant,
        "injection_depth": depth,
        "injection_duration_min": duration_min,
        "injection_period_d": period_d,
        "injection_flux_scale": flux_scale,
        "injection_cadences": int(in_event.sum()),
        "recovered_depth": recovered_depth,
        "depth_retention_frac": recovered_depth / depth if depth else np.nan,
    }


def cases_from_h5(
    paths: list[Path],
    *,
    apertures: tuple[str, ...],
    tmag_min: float | None,
    limit: int | None,
) -> list[CurveCase]:
    from twirl.lightcurves.tglc_h5_reader import read_tglc_h5

    cases: list[CurveCase] = []
    for path in paths:
        try:
            lc = read_tglc_h5(path)
        except Exception as exc:
            print(f"[warn] failed to read {path}: {type(exc).__name__}: {exc}", file=sys.stderr)
            continue
        if tmag_min is not None and np.isfinite(lc.tmag) and lc.tmag < tmag_min:
            continue
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
                    raw_flux_source="synthesized" if ap_data.flux_was_synthesized else "raw",
                )
            )
        if limit is not None and len(cases) >= limit * len(apertures):
            break
    return cases


def synthetic_cases(n: int, *, seed: int, apertures: tuple[str, ...]) -> list[CurveCase]:
    rng = np.random.default_rng(seed)
    cadence_d = 200.0 / 86400.0
    orbit_a = np.arange(0.0, 12.8, cadence_d)
    orbit_b = np.arange(14.0, 26.8, cadence_d)
    time = np.concatenate([orbit_a, orbit_b])
    cases: list[CurveCase] = []
    star_flux_grid = np.geomspace(0.5, 600.0, max(n, 2))
    for i in range(n):
        star_flux = float(star_flux_grid[i])
        noise = 70.0 + 12.0 * rng.random()
        tmag = 20.0 - 2.5 * math.log10(max(star_flux, 0.5) / 10.0)
        slow = 0.12 * star_flux * np.sin(2.0 * np.pi * time / 3.7 + rng.uniform(0, 2 * np.pi))
        ramp = 0.04 * star_flux * (time - np.nanmedian(time)) / np.ptp(time)
        base_flux = star_flux + slow + ramp + rng.normal(0.0, noise, size=time.size)
        quality = np.zeros(time.size, dtype=np.int64)
        quality[rng.random(time.size) < 0.015] = 2
        nan_mask = rng.random(time.size) < 0.005
        base_flux[nan_mask] = np.nan
        flux_err = np.full(time.size, noise, dtype=np.float64)
        flux_err[nan_mask] = np.nan
        for ap_i, ap in enumerate(apertures):
            ap_scale = 1.0 + 0.25 * (ap_i - 1)
            flux = ap_scale * base_flux + rng.normal(0.0, 0.1 * noise, size=time.size)
            cases.append(
                CurveCase(
                    case_id=f"synthetic{i:03d}_{ap}",
                    tic=900000000 + i,
                    sector=56,
                    orbit=-1,
                    cam=-1,
                    ccd=-1,
                    tmag=float(tmag),
                    aperture=ap,
                    time=time.copy(),
                    flux=flux,
                    flux_err=flux_err.copy(),
                    quality=quality.copy(),
                    path="synthetic",
                    raw_flux_source="synthetic",
                )
            )
    return cases


def collect_h5_paths(args: argparse.Namespace) -> list[Path]:
    paths: list[Path] = []
    for p in args.h5:
        paths.append(Path(p))
    for pat in args.h5_glob:
        paths.extend(Path(p) for p in sorted(glob.glob(pat)))
    seen: set[Path] = set()
    out: list[Path] = []
    for p in paths:
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def write_markdown_summary(df: pd.DataFrame, inj_df: pd.DataFrame, out_path: Path) -> None:
    lines = ["# Flux Detrend Model Comparison", ""]
    if df.empty:
        lines.append("No cases were processed.")
        out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    grouped = (
        df.groupby("variant", dropna=False)
        .agg(
            cases=("case_id", "nunique"),
            failures=("failed", "sum"),
            median_mad=("det_mad_quality0", "median"),
            median_rms=("det_rms_quality0", "median"),
            median_finite_q0=("n_det_finite_quality0", "median"),
            raw_neg_q0=("raw_negative_quality0", "sum"),
            retained_neg_q0=("raw_negative_quality0_retained", "sum"),
        )
        .reset_index()
    )
    grouped["retained_neg_q0_frac"] = grouped["retained_neg_q0"] / grouped["raw_neg_q0"].replace(0, np.nan)
    lines.append("## Aggregate Metrics")
    lines.append("")
    lines.append(frame_to_markdown(grouped, floatfmt=".5g"))
    lines.append("")

    if not inj_df.empty:
        inj_grouped = (
            inj_df.groupby("variant", dropna=False)
            .agg(
                injection_cases=("case_id", "nunique"),
                median_depth_retention=("depth_retention_frac", "median"),
                p16_depth_retention=("depth_retention_frac", lambda x: np.nanpercentile(x, 16)),
                p84_depth_retention=("depth_retention_frac", lambda x: np.nanpercentile(x, 84)),
            )
            .reset_index()
        )
        lines.append("## Injection Preservation")
        lines.append("")
        lines.append(frame_to_markdown(inj_grouped, floatfmt=".5g"))
        lines.append("")

    worst = df.sort_values(["failed", "det_mad_quality0"], ascending=[False, False]).head(12)
    lines.append("## Highest-Risk Case Rows")
    lines.append("")
    cols = [
        "case_id",
        "variant",
        "tmag",
        "scale_source",
        "scale",
        "raw_negative_quality0",
        "n_det_finite_quality0",
        "det_mad_quality0",
        "failed",
        "failure_reason",
    ]
    lines.append(frame_to_markdown(worst[cols], floatfmt=".5g"))
    while lines and lines[-1] == "":
        lines.pop()
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def frame_to_markdown(df: pd.DataFrame, *, floatfmt: str = ".5g") -> str:
    """Small dependency-free markdown table writer.

    Pandas' built-in ``to_markdown`` requires the optional ``tabulate``
    package, which is not part of TWIRL's runtime dependencies.
    """
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
    widths = [
        max(len(headers[i]), *(len(row[i]) for row in body))
        for i in range(len(headers))
    ]
    lines = [
        "| " + " | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    for row in body:
        lines.append("| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(headers))) + " |")
    return "\n".join(lines)


def choose_plot_cases(df: pd.DataFrame, n: int) -> list[str]:
    if n <= 0 or df.empty:
        return []
    pivot = df.pivot_table(
        index="case_id",
        values=["failed", "raw_negative_quality0", "det_mad_quality0", "tmag"],
        aggfunc={
            "failed": "sum",
            "raw_negative_quality0": "max",
            "det_mad_quality0": "max",
            "tmag": "first",
        },
    )
    pivot = pivot.sort_values(
        ["failed", "raw_negative_quality0", "det_mad_quality0", "tmag"],
        ascending=[False, False, False, False],
    )
    return list(pivot.head(n).index)


def plot_case(case: CurveCase, results: dict[str, FluxDetrendResult], out_dir: Path) -> None:
    template = apply_twirl_style("full_page")
    colors = dict(zip(results, get_ordered_palette(len(results), "viridis")))
    q0 = np.asarray(case.quality) == 0
    finite = np.isfinite(case.flux)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(template["figsize"][0], template["figsize"][1] * 1.25),
        sharex=True,
        constrained_layout=True,
    )
    ax = axes[0]
    ax.scatter(case.time[finite & ~q0], case.flux[finite & ~q0], s=3, c="0.78", lw=0, label="flagged")
    ax.scatter(case.time[finite & q0], case.flux[finite & q0], s=3, c="0.25", lw=0, label="quality=0")
    for name, result in results.items():
        ax.plot(case.time, result.cotrend, lw=0.9, color=colors[name], label=f"{name} cotrend")
    ax.set_ylabel("raw flux")
    ax.legend(loc="upper right", ncol=2)

    ax = axes[1]
    for name, result in results.items():
        y = result.det_flux
        keep = np.isfinite(case.time) & np.isfinite(y)
        label = f"{name} ({result.scale_source})"
        ax.plot(case.time[keep], y[keep], lw=0.65, color=colors[name], alpha=0.9, label=label)
    ax.axhline(1.0, color="black", lw=0.7, ls=":")
    ax.set_xlabel("BTJD")
    ax.set_ylabel("relative flux")
    ax.legend(loc="upper right", ncol=2)
    stem = safe_stem(case.case_id)
    fig.savefig(out_dir / f"{stem}.png", dpi=180)
    fig.savefig(out_dir / f"{stem}.pdf")
    plt.close(fig)


def safe_stem(text: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in text)


def parse_apertures(values: list[str]) -> tuple[str, ...]:
    apertures: list[str] = []
    for v in values:
        if v.lower() == "all":
            apertures.extend(APERTURE_KEYS)
        else:
            apertures.append(v)
    bad = [a for a in apertures if a not in APERTURE_KEYS]
    if bad:
        raise SystemExit(f"unknown aperture(s): {', '.join(bad)}; valid: {', '.join(APERTURE_KEYS)}")
    return tuple(dict.fromkeys(apertures))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5", action="append", default=[], help="TGLC HDF5 file; repeatable.")
    ap.add_argument("--h5-glob", action="append", default=[], help="Quoted glob for TGLC HDF5 files.")
    ap.add_argument("--synthetic", type=int, default=0, help="Generate N synthetic faint light curves.")
    ap.add_argument("--seed", type=int, default=20260517)
    ap.add_argument("--aperture", action="append", default=["Primary"], help="'Small', 'Primary', 'Large', or 'all'.")
    ap.add_argument("--tmag-min", type=float, default=None)
    ap.add_argument("--limit", type=int, default=None, help="Maximum HDF5 files to read after filtering.")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("reports/stage1_lightcurves/detrend_experiments/latest"),
    )
    ap.add_argument("--plot-cases", type=int, default=6)
    ap.add_argument("--inject-depth", type=float, default=0.10)
    ap.add_argument("--inject-duration-min", type=float, default=10.0)
    ap.add_argument("--inject-period-d", type=float, default=1.40793903)
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    apertures = parse_apertures(args.aperture)
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    cases: list[CurveCase] = []
    h5_paths = collect_h5_paths(args)
    if h5_paths:
        cases.extend(
            cases_from_h5(
                h5_paths,
                apertures=apertures,
                tmag_min=args.tmag_min,
                limit=args.limit,
            )
        )
    if args.synthetic:
        cases.extend(synthetic_cases(args.synthetic, seed=args.seed, apertures=apertures))

    if not cases:
        print("[compare-detrend] no cases selected", file=sys.stderr)
        return 2

    rows: list[dict] = []
    inj_rows: list[dict] = []
    result_cache: dict[str, dict[str, FluxDetrendResult]] = {}
    for i, case in enumerate(cases, 1):
        print(f"[compare-detrend] {i}/{len(cases)} {case.case_id}")
        per_case: dict[str, FluxDetrendResult] = {}
        for variant, cfg in VARIANTS.items():
            result = flux_space_detrend_result(
                case.time,
                case.flux,
                quality=case.quality,
                flux_err=case.flux_err,
                cfg=cfg,
            )
            per_case[variant] = result
            rows.append(summarize_result(case, variant, result))
            inj = injection_preservation(
                case,
                variant,
                cfg,
                depth=args.inject_depth,
                duration_min=args.inject_duration_min,
                period_d=args.inject_period_d,
            )
            if inj:
                inj_rows.append(inj)
        result_cache[case.case_id] = per_case

    df = pd.DataFrame(rows)
    inj_df = pd.DataFrame(inj_rows)
    summary_csv = out_dir / "detrend_model_summary.csv"
    injection_csv = out_dir / "detrend_injection_summary.csv"
    df.to_csv(summary_csv, index=False)
    inj_df.to_csv(injection_csv, index=False)
    write_markdown_summary(df, inj_df, out_dir / "summary.md")

    selected = set(choose_plot_cases(df, args.plot_cases))
    cases_by_id = {case.case_id: case for case in cases}
    for case_id in selected:
        plot_case(cases_by_id[case_id], result_cache[case_id], plot_dir)

    manifest = {
        "n_cases": len(cases),
        "n_rows": len(rows),
        "n_injection_rows": len(inj_rows),
        "variants": {
            name: {
                "output_mode": cfg.output_mode,
                "scale_strategy": cfg.scale_strategy,
                "bkspace_d": cfg.bkspace_d,
                "sigma_clip": cfg.sigma_clip,
            }
            for name, cfg in VARIANTS.items()
        },
        "inputs": {
            "h5": [str(p) for p in h5_paths],
            "synthetic": args.synthetic,
            "apertures": apertures,
            "tmag_min": args.tmag_min,
            "limit": args.limit,
        },
        "outputs": {
            "summary_csv": str(summary_csv),
            "injection_csv": str(injection_csv),
            "summary_md": str(out_dir / "summary.md"),
            "plot_dir": str(plot_dir),
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"[compare-detrend] wrote {summary_csv}")
    print(f"[compare-detrend] wrote {out_dir / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
