#!/usr/bin/env python3
"""Duration-aware S56 BLS recovery visualizations."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style  # noqa: E402


DEFAULT_BLS_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/"
    / "small_pair_200k/injection_bls_recoveries.csv"
)
DEFAULT_LEO_QUEUE_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/"
    / "review_queue.csv"
)
DEFAULT_LEO_METRICS_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/"
    / "leo_metrics.csv"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/"
    / "duration_aware"
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    text = series.astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y"})


def classify_recovery(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["strict_top1_recovered"] = out["recovery_status"].astype(str).eq("bls_recovered")
    out["topn_exact_recovered"] = out["topn_recovery_status"].astype(str).eq("bls_topn_recovered")
    out["topn_harmonic_recovered"] = out["topn_recovery_status"].astype(str).eq("bls_topn_harmonic_match")
    exact_cols = [c for c in out.columns if c.startswith("topn_exact_recovered_")]
    harmonic_cols = [c for c in out.columns if c.startswith("topn_harmonic_match_")]
    if exact_cols:
        out["topn_exact_recovered"] |= pd.concat([_bool_series(out[c]) for c in exact_cols], axis=1).any(axis=1)
    if harmonic_cols:
        out["topn_harmonic_recovered"] |= pd.concat([_bool_series(out[c]) for c in harmonic_cols], axis=1).any(axis=1)
    out["any_exact_or_harmonic_recovered"] = (
        out["strict_top1_recovered"] | out["topn_exact_recovered"] | out["topn_harmonic_recovered"]
    )
    return out


def prepare_bls_frame(path: Path) -> pd.DataFrame:
    df = classify_recovery(pd.read_csv(path))
    for col in (
        "truth_period_d",
        "truth_duration_min",
        "truth_sampled_model_depth",
        "truth_model_depth",
        "truth_depth",
        "truth_radius_rearth",
        "tmag",
    ):
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    depth = df.get("truth_sampled_model_depth", pd.Series(np.nan, index=df.index)).copy()
    for fallback in ("truth_model_depth", "truth_depth"):
        missing = ~np.isfinite(depth) | (depth <= 0)
        if missing.any() and fallback in df:
            depth.loc[missing] = pd.to_numeric(df.loc[missing, fallback], errors="coerce")
    df["plot_depth_frac"] = depth.clip(lower=1e-6, upper=1.0)
    df["plot_depth_pct"] = 100.0 * df["plot_depth_frac"]
    df["plot_radius_rearth"] = df["truth_radius_rearth"]
    good = (
        np.isfinite(df["truth_period_d"])
        & (df["truth_period_d"] > 0)
        & np.isfinite(df["truth_duration_min"])
        & (df["truth_duration_min"] > 0)
        & np.isfinite(df["plot_depth_pct"])
        & (df["plot_depth_pct"] > 0)
        & np.isfinite(df["plot_radius_rearth"])
        & (df["plot_radius_rearth"] > 0)
        & np.isfinite(df["tmag"])
    )
    return df[good].copy()


def _standardized_design(raw_features: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = np.nanmean(raw_features, axis=0)
    scale = np.nanstd(raw_features, axis=0)
    scale[scale == 0] = 1.0
    x = (raw_features - center) / scale
    x = np.column_stack([np.ones(len(x)), x])
    return x, center, scale


def fit_logistic_boundary(df: pd.DataFrame, *, ridge: float = 1.0) -> dict[str, Any]:
    log_period = np.log10(df["truth_period_d"].to_numpy(dtype=float))
    log_duration = np.log10(df["truth_duration_min"].to_numpy(dtype=float))
    log_radius = np.log10(df["plot_radius_rearth"].to_numpy(dtype=float))
    # BLS SNR scales approximately as depth * sqrt(total in-transit cadence
    # count). With constant WD radius in this injected grid, depth ~ R_p^2 and
    # total in-transit cadence count scales like duration / period.
    bls_radius_proxy = 2.0 * log_radius + 0.5 * log_duration - 0.5 * log_period
    raw = np.column_stack([bls_radius_proxy, df["tmag"].to_numpy(dtype=float)])
    y = df["any_exact_or_harmonic_recovered"].fillna(False).astype(bool).to_numpy(dtype=float)
    x, center, scale = _standardized_design(raw)
    beta = np.zeros(x.shape[1], dtype=float)
    penalty = np.eye(x.shape[1], dtype=float) * float(ridge)
    penalty[0, 0] = 0.0
    converged = False
    for iteration in range(80):
        eta = np.clip(x @ beta, -35.0, 35.0)
        p = 1.0 / (1.0 + np.exp(-eta))
        grad = x.T @ (p - y) + penalty @ beta
        weights = p * (1.0 - p)
        hessian = (x.T * weights) @ x + penalty
        try:
            step = np.linalg.solve(hessian, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(hessian, grad, rcond=None)[0]
        beta -= step
        if float(np.max(np.abs(step))) < 1e-6:
            converged = True
            break

    eta = np.clip(x @ beta, -35.0, 35.0)
    prob = 1.0 / (1.0 + np.exp(-eta))
    eps = 1e-12
    logloss = -float(np.mean(y * np.log(prob + eps) + (1.0 - y) * np.log(1.0 - prob + eps)))
    order = np.argsort(prob)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(len(prob))
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos and n_neg:
        auc = float((ranks[y.astype(bool)].sum() - n_pos * (n_pos - 1) / 2.0) / (n_pos * n_neg))
    else:
        auc = float("nan")

    # Convert from standardized coefficients to raw-feature coefficients:
    # logit = raw_intercept + sum(raw_coef_j * raw_feature_j).
    raw_coef = beta[1:] / scale
    raw_intercept = beta[0] - float(np.sum(beta[1:] * center / scale))
    return {
        "beta_standardized": beta,
        "center": center,
        "scale": scale,
        "raw_intercept": raw_intercept,
        "raw_coef": raw_coef,
        "feature_names": ["log10_radius2_sqrt_duration_over_period", "tmag"],
        "model_type": "physical_radius_proxy",
        "converged": converged,
        "iterations": iteration + 1,
        "ridge": float(ridge),
        "logloss": logloss,
        "auc": auc,
        "n_rows": int(len(df)),
        "n_recovered": n_pos,
    }


def boundary_radius_rearth(
    model: dict[str, Any],
    period_d: np.ndarray,
    duration_min: np.ndarray,
    tmag: float,
    *,
    probability: float = 0.5,
) -> np.ndarray:
    raw_coef = np.asarray(model["raw_coef"], dtype=float)
    intercept = float(model["raw_intercept"])
    logit = np.log(probability / (1.0 - probability))
    log_period = np.log10(period_d)
    log_duration = np.log10(duration_min)
    denom = raw_coef[0]
    if abs(float(denom)) < 1e-8:
        return np.full_like(log_period, np.nan, dtype=float)
    proxy_boundary = (logit - intercept - raw_coef[1] * tmag) / denom
    log_radius = 0.5 * (proxy_boundary - 0.5 * log_duration + 0.5 * log_period)
    return np.power(10.0, log_radius)


def _fraction_grid(
    x: pd.Series | np.ndarray,
    y: pd.Series | np.ndarray,
    recovered: pd.Series | np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    recovered_arr = np.asarray(recovered, dtype=float)
    good = np.isfinite(x_arr) & np.isfinite(y_arr) & np.isfinite(recovered_arr)
    count, _, _ = np.histogram2d(x_arr[good], y_arr[good], bins=(x_edges, y_edges))
    hit, _, _ = np.histogram2d(
        x_arr[good],
        y_arr[good],
        bins=(x_edges, y_edges),
        weights=recovered_arr[good],
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        fraction = hit / count
    return fraction, hit, count


def _gaussian_kernel(size: int = 5, sigma: float = 1.0) -> np.ndarray:
    axis = np.arange(size, dtype=float) - (size - 1.0) / 2.0
    xx, yy = np.meshgrid(axis, axis, indexing="xy")
    kernel = np.exp(-(xx**2 + yy**2) / (2.0 * sigma**2))
    return kernel / np.sum(kernel)


def _smooth_grid(grid: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    pad_y = kernel.shape[0] // 2
    pad_x = kernel.shape[1] // 2
    padded = np.pad(grid, ((pad_y, pad_y), (pad_x, pad_x)), mode="constant", constant_values=0.0)
    out = np.zeros_like(grid, dtype=float)
    for iy in range(out.shape[0]):
        for ix in range(out.shape[1]):
            out[iy, ix] = float(np.sum(padded[iy : iy + kernel.shape[0], ix : ix + kernel.shape[1]] * kernel))
    return out


def _smoothed_fraction(hit: np.ndarray, count: np.ndarray, *, min_effective_count: float = 2.0) -> np.ndarray:
    kernel = _gaussian_kernel(size=5, sigma=1.0)
    smooth_hit = _smooth_grid(hit, kernel)
    smooth_count = _smooth_grid(count, kernel)
    with np.errstate(invalid="ignore", divide="ignore"):
        fraction = smooth_hit / smooth_count
    fraction[smooth_count < min_effective_count] = np.nan
    return fraction


def plot_period_radius_boundary_lines(df: pd.DataFrame, model: dict[str, Any], out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    period_grid = np.geomspace(0.12, 13.0, 240)
    duration_values = [3.0, 5.0, 9.0, 16.0]
    duration_labels = ["3 min", "5 min", "9 min", "16 min"]
    tmag_values = [17.5, 18.5, 19.0, 19.5]
    colors = ["#2b6cb0", "#38a169", "#dd6b20", "#c53030"]
    max_radius = float(np.nanmax(df["plot_radius_rearth"]))
    min_radius = max(0.12, float(np.nanpercentile(df["plot_radius_rearth"], 0.5)))

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.8), sharex=True, sharey=True)
    recovered = df["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
    for ax, duration, label in zip(axes.ravel(), duration_values, duration_labels):
        width = 0.35 if duration < 4.0 else 0.45
        nearby = df["truth_duration_min"].between(duration / (1.0 + width), duration * (1.0 + width))
        sample = df[nearby].sample(n=min(850, int(nearby.sum())), random_state=int(duration * 1000)) if nearby.any() else df.iloc[0:0]
        if len(sample):
            sample_recovered = sample["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
            ax.scatter(
                sample.loc[~sample_recovered, "truth_period_d"],
                sample.loc[~sample_recovered, "plot_radius_rearth"],
                s=4,
                color="0.65",
                alpha=0.20,
                linewidths=0,
                rasterized=True,
            )
            ax.scatter(
                sample.loc[sample_recovered, "truth_period_d"],
                sample.loc[sample_recovered, "plot_radius_rearth"],
                s=6,
                color="#111827",
                alpha=0.38,
                linewidths=0,
                rasterized=True,
            )
        for tmag, color in zip(tmag_values, colors):
            radius = boundary_radius_rearth(model, period_grid, np.full_like(period_grid, duration), tmag)
            radius = np.where((radius >= min_radius) & (radius <= max_radius * 1.08), radius, np.nan)
            ax.plot(period_grid, radius, color=color, linewidth=2.0, label=f"Tmag {tmag:g}")
        ax.axvline(0.216, color="0.35", linestyle=":", linewidth=1.0)
        ax.set_title(f"Duration {label}")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.12, 13.0)
        ax.set_ylim(0.12, max_radius * 1.12)
        ax.set_xticks([0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
        ax.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10"])
        ax.grid(True, which="major", color="0.88", linewidth=0.7)
        n_nearby = int(nearby.sum())
        n_recovered = int((nearby & recovered).sum())
        ax.text(
            0.03,
            0.96,
            f"local injections: {n_nearby}\nBLS recovered: {n_recovered}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            color="0.25",
        )
    for ax in axes[-1, :]:
        ax.set_xlabel("Injected period [days]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Planet radius [R_earth]")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.01))
    fig.subplots_adjust(top=0.88, bottom=0.14, left=0.10, right=0.98, wspace=0.10, hspace=0.20)
    fig.text(
        0.5,
        0.02,
        "Lines show the fitted 50% BLS recovery radius for fixed transit durations; grey/black points are nearby injected misses/recoveries.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "period_radius_50pct_boundary_by_duration.png"
    pdf = out_dir / "period_radius_50pct_boundary_by_duration.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"period_radius_boundary_png": str(png), "period_radius_boundary_pdf": str(pdf)}


def plot_period_radius_recovery_maps(df: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    period_edges = np.geomspace(0.12, 13.0, 27)
    radius_edges = np.geomspace(0.12, 18.0, 25)
    period_centers = np.sqrt(period_edges[:-1] * period_edges[1:])
    radius_centers = np.sqrt(radius_edges[:-1] * radius_edges[1:])
    duration_bins = [
        ("2-4 min", 2.0, 4.0),
        ("4-8 min", 4.0, 8.0),
        ("8-25 min", 8.0, 25.0),
    ]
    tmag_bins = [
        ("Tmag < 18", -np.inf, 18.0),
        ("18 <= Tmag < 19", 18.0, 19.0),
        ("19 <= Tmag < 20", 19.0, 20.0),
    ]
    recovered = df["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#f2f2f2")
    fig, axes = plt.subplots(
        len(tmag_bins),
        len(duration_bins),
        figsize=(12.2, 9.2),
        sharex=True,
        sharey=True,
    )
    rows: list[dict[str, Any]] = []
    mesh = None
    for row_idx, (tmag_label, tmag_lo, tmag_hi) in enumerate(tmag_bins):
        for col_idx, (duration_label, duration_lo, duration_hi) in enumerate(duration_bins):
            ax = axes[row_idx, col_idx]
            mask = (
                df["truth_duration_min"].between(duration_lo, duration_hi, inclusive="left")
                & (df["tmag"] >= tmag_lo)
                & (df["tmag"] < tmag_hi)
            )
            subset = df[mask]
            fraction, hit, count = _fraction_grid(
                subset["truth_period_d"],
                subset["plot_radius_rearth"],
                subset["any_exact_or_harmonic_recovered"].astype(float),
                period_edges,
                radius_edges,
            )
            smooth_fraction = _smoothed_fraction(hit, count, min_effective_count=1.6)
            mesh = ax.pcolormesh(
                period_edges,
                radius_edges,
                smooth_fraction.T,
                cmap=cmap,
                vmin=0.0,
                vmax=1.0,
                shading="auto",
            )
            support = np.isfinite(smooth_fraction)
            if np.any(support) and np.nanmin(smooth_fraction) <= 0.5 <= np.nanmax(smooth_fraction):
                ax.contour(
                    period_centers,
                    radius_centers,
                    smooth_fraction.T,
                    levels=[0.5],
                    colors=["white"],
                    linewidths=1.3,
                )
            if len(subset):
                sample = subset.sample(n=min(500, len(subset)), random_state=1000 + row_idx * 10 + col_idx)
                sample_recovered = sample["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
                ax.scatter(
                    sample.loc[~sample_recovered, "truth_period_d"],
                    sample.loc[~sample_recovered, "plot_radius_rearth"],
                    s=2,
                    color="black",
                    alpha=0.10,
                    linewidths=0,
                    rasterized=True,
                )
                ax.scatter(
                    sample.loc[sample_recovered, "truth_period_d"],
                    sample.loc[sample_recovered, "plot_radius_rearth"],
                    s=3,
                    color="white",
                    alpha=0.22,
                    linewidths=0,
                    rasterized=True,
                )
            ax.axvline(0.216, color="white", linestyle=":", linewidth=0.9, alpha=0.9)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(0.12, 13.0)
            ax.set_ylim(0.12, 18.0)
            ax.grid(False)
            if row_idx == 0:
                ax.set_title(duration_label)
            if col_idx == 0:
                ax.set_ylabel(f"{tmag_label}\nPlanet radius [R_earth]")
            if row_idx == len(tmag_bins) - 1:
                ax.set_xlabel("Injected period [days]")
            ax.text(
                0.04,
                0.94,
                f"n={len(subset)}\nrec={int(recovered[mask].sum())}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7.5,
                color="white",
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "black", "alpha": 0.35, "edgecolor": "none"},
            )
            for ix in range(len(period_edges) - 1):
                for iy in range(len(radius_edges) - 1):
                    if count[ix, iy] > 0:
                        rows.append(
                            {
                                "tmag_bin": tmag_label,
                                "duration_bin": duration_label,
                                "period_lo_d": period_edges[ix],
                                "period_hi_d": period_edges[ix + 1],
                                "radius_lo_rearth": radius_edges[iy],
                                "radius_hi_rearth": radius_edges[iy + 1],
                                "count": int(count[ix, iy]),
                                "recovered": int(hit[ix, iy]),
                                "recovery_fraction": float(fraction[ix, iy]) if np.isfinite(fraction[ix, iy]) else np.nan,
                                "smoothed_recovery_fraction": (
                                    float(smooth_fraction[ix, iy]) if np.isfinite(smooth_fraction[ix, iy]) else np.nan
                                ),
                            }
                        )
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), fraction=0.030, pad=0.018)
        cbar.set_label("BLS recovery fraction")
    fig.subplots_adjust(left=0.10, right=0.89, top=0.93, bottom=0.11, wspace=0.06, hspace=0.08)
    fig.text(
        0.5,
        0.02,
        "White contours mark 50% empirical recovery after light count-weighted smoothing; blank cells have little local injection support.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "period_radius_recovery_fraction_by_duration_tmag.png"
    pdf = out_dir / "period_radius_recovery_fraction_by_duration_tmag.pdf"
    csv = out_dir / "period_radius_recovery_fraction_by_duration_tmag.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {
        "period_radius_recovery_map_png": str(png),
        "period_radius_recovery_map_pdf": str(pdf),
        "period_radius_recovery_map_csv": str(csv),
    }


def plot_snr_proxy_recovery_map(df: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    local = df.copy()
    if "truth_n_good_in_transit" in local:
        n_in = pd.to_numeric(local["truth_n_good_in_transit"], errors="coerce").clip(lower=1)
    else:
        n_in = 24.0 * 60.0 / 200.0 * local["truth_duration_min"].clip(lower=0) * (
            27.0 / local["truth_period_d"].clip(lower=0.01)
        )
    local["snr_radius_proxy"] = local["plot_radius_rearth"] ** 2 * np.sqrt(n_in)
    local = local[np.isfinite(local["snr_radius_proxy"]) & (local["snr_radius_proxy"] > 0)].copy()
    x_edges = np.geomspace(np.nanpercentile(local["snr_radius_proxy"], 0.5), np.nanpercentile(local["snr_radius_proxy"], 99.5), 32)
    y_edges = np.linspace(15.0, 20.2, 27)
    fraction, hit, count = _fraction_grid(
        local["snr_radius_proxy"],
        local["tmag"],
        local["any_exact_or_harmonic_recovered"].astype(float),
        x_edges,
        y_edges,
    )
    smooth_fraction = _smoothed_fraction(hit, count, min_effective_count=3.0)

    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#f2f2f2")
    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    mesh = ax.pcolormesh(x_edges, y_edges, smooth_fraction.T, cmap=cmap, vmin=0.0, vmax=1.0, shading="auto")
    ax.set_xscale("log")
    ax.set_xlabel(r"Radius-SNR proxy: $R_p^2 \sqrt{N_\mathrm{in}}$")
    ax.set_ylabel("Tmag")
    ax.invert_yaxis()
    x_centers = np.sqrt(x_edges[:-1] * x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    if np.nanmin(smooth_fraction) <= 0.5 <= np.nanmax(smooth_fraction):
        ax.contour(x_centers, y_centers, smooth_fraction.T, levels=[0.5], colors=["white"], linewidths=1.4)
    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("BLS recovery fraction")
    ax.grid(True, which="major", color="white", alpha=0.18)
    ax.set_title("BLS recovery collapses toward a radius-SNR boundary")
    fig.subplots_adjust(bottom=0.18, left=0.12, right=0.96, top=0.90)
    fig.text(
        0.5,
        0.03,
        "The proxy uses planet radius for depth and the in-transit cadence count for duration/period sampling; it is diagnostic, not a calibrated MES.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "radius_snr_proxy_recovery_by_tmag.png"
    pdf = out_dir / "radius_snr_proxy_recovery_by_tmag.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"radius_snr_proxy_png": str(png), "radius_snr_proxy_pdf": str(pdf)}


def plot_boundary_floor_maps(df: pd.DataFrame, model: dict[str, Any], out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    apply_twirl_style("full_page")
    period_grid = np.geomspace(0.12, 13.0, 120)
    duration_grid = np.geomspace(2.0, 25.0, 100)
    pp, dd = np.meshgrid(period_grid, duration_grid, indexing="xy")
    tmag_values = [17.5, 18.5, 19.0, 19.5]
    labels = ["Tmag=17.5", "Tmag=18.5", "Tmag=19.0", "Tmag=19.5"]
    max_injected_radius = float(np.nanmax(df["plot_radius_rearth"]))
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.1), sharex=True, sharey=True)
    mesh = None
    for ax, tmag, label in zip(axes.ravel(), tmag_values, labels):
        radius = boundary_radius_rearth(model, pp, dd, tmag)
        radius = np.where((radius >= 0.1) & (radius <= max_injected_radius), radius, np.nan)
        mesh = ax.pcolormesh(
            period_grid,
            duration_grid,
            radius,
            cmap="viridis",
            norm=LogNorm(vmin=0.3, vmax=max_injected_radius),
            shading="auto",
        )
        ax.scatter(
            df.loc[(df["tmag"] >= tmag - 0.5) & (df["tmag"] < tmag + 0.5), "truth_period_d"],
            df.loc[(df["tmag"] >= tmag - 0.5) & (df["tmag"] < tmag + 0.5), "truth_duration_min"],
            s=2,
            c="white",
            alpha=0.16,
            linewidths=0,
        )
        ax.axvline(0.216, color="white", linestyle=":", linewidth=1.0)
        ax.set_title(label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.12, 13.0)
        ax.set_ylim(2.0, 25.0)
        ax.set_xticks([0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
        ax.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10"])
        ax.set_yticks([2.0, 5.0, 10.0, 20.0])
        ax.set_yticklabels(["2", "5", "10", "20"])
        ax.grid(True, which="major", color="white", alpha=0.18)
        if np.count_nonzero(np.isfinite(radius)) < 0.03 * radius.size:
            ax.text(
                0.5,
                0.5,
                "No 50% boundary\nwithin injected radii",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9,
                color="0.25",
            )
    for ax in axes[-1, :]:
        ax.set_xlabel("Injected period [days]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Injected duration [min]")
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), fraction=0.032, pad=0.02)
        cbar.set_label("Companion radius at 50% BLS recovery [R_earth]")
    fig.suptitle("Duration-aware BLS radius recovery boundary", y=0.995)
    fig.text(
        0.5,
        0.005,
        "Color is the 50% BLS radius cutoff from a physically constrained R_p^2 sqrt(duration/period) model; white means the cutoff is beyond the injected radius range.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "duration_radius_50pct_boundary_floor_maps.png"
    pdf = out_dir / "duration_radius_50pct_boundary_floor_maps.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"duration_radius_floor_png": str(png), "duration_radius_floor_pdf": str(pdf)}


def plot_boundary_3d_planes(df: pd.DataFrame, model: dict[str, Any], out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    apply_twirl_style("full_page")
    period_grid = np.geomspace(0.12, 13.0, 34)
    duration_grid = np.geomspace(2.0, 25.0, 28)
    pp, dd = np.meshgrid(period_grid, duration_grid, indexing="xy")
    fig = plt.figure(figsize=(9.2, 6.8))
    ax = fig.add_subplot(111, projection="3d")
    colors = {17.5: "#1f77b4", 18.5: "#2ca02c", 19.0: "#ff7f0e", 19.5: "#d62728"}
    max_injected_radius = float(np.nanmax(df["plot_radius_rearth"]))
    for tmag, color in colors.items():
        radius = boundary_radius_rearth(model, pp, dd, tmag)
        radius = np.where((radius >= 0.1) & (radius <= max_injected_radius), radius, np.nan)
        z = np.log10(radius)
        ax.plot_surface(
            np.log10(pp),
            np.log10(dd),
            z,
            color=color,
            alpha=0.35,
            linewidth=0,
            antialiased=True,
        )
        # Matplotlib does not propagate surface labels reliably; add matching edge line.
        mid = len(duration_grid) // 2
        ax.plot(
            np.log10(period_grid),
            np.full_like(period_grid, np.log10(duration_grid[mid])),
            z[mid, :],
            color=color,
            linewidth=2.0,
            label=f"Tmag={tmag:g}",
        )
    recovered = df["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
    sample = df.sample(n=min(1400, len(df)), random_state=56)
    sample_recovered = sample["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
    ax.scatter(
        np.log10(sample.loc[~sample_recovered, "truth_period_d"]),
        np.log10(sample.loc[~sample_recovered, "truth_duration_min"]),
        np.log10(sample.loc[~sample_recovered, "plot_radius_rearth"]),
        c="0.6",
        s=4,
        alpha=0.10,
        marker="x",
        depthshade=False,
    )
    ax.scatter(
        np.log10(sample.loc[sample_recovered, "truth_period_d"]),
        np.log10(sample.loc[sample_recovered, "truth_duration_min"]),
        np.log10(sample.loc[sample_recovered, "plot_radius_rearth"]),
        c=sample.loc[sample_recovered, "tmag"],
        cmap="viridis_r",
        s=9,
        alpha=0.55,
        depthshade=False,
    )
    ax.set_xlabel("log10 period [days]")
    ax.set_ylabel("log10 duration [min]")
    ax.set_zlabel("log radius [R_earth]", labelpad=5)
    ax.text2D(0.015, 0.48, "log radius [R_earth]", transform=ax.transAxes, rotation=90, fontsize=9)
    ax.set_title("BLS 50% radius recovery boundary as magnitude-shifted planes")
    ax.view_init(elev=24, azim=-128)
    ax.legend(loc="upper left", fontsize=8)
    fig.text(
        0.5,
        0.02,
        "Surfaces are logistic 50% boundaries; points are a random audit subset of dense injections.",
        ha="center",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "duration_radius_3d_boundary_planes.png"
    pdf = out_dir / "duration_radius_3d_boundary_planes.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"duration_radius_3d_png": str(png), "duration_radius_3d_pdf": str(pdf)}


def leo_bls_contingency(queue_csv: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    q = classify_recovery(pd.read_csv(queue_csv))
    if "source_kind" in q:
        q = q[q["source_kind"].astype(str).eq("injection_recovery")].copy()
    q["leo_recovered_pc"] = q["leo_class"].astype(str).eq("PC")
    q["leo_recovered_pc_or_fp"] = q["leo_class"].astype(str).isin(["PC", "FP"])
    q["bls_recovered"] = q["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
    table = pd.crosstab(q["bls_recovered"], q["leo_class"], margins=True)
    rows: dict[str, Any] = {
        "n": int(len(q)),
        "leo_class_counts": q["leo_class"].value_counts(dropna=False).to_dict(),
        "bls_any_recovered_n": int(q["bls_recovered"].sum()),
        "bls_any_recovered_frac": float(q["bls_recovered"].mean()),
    }
    for label, pred_col in [("leo_pc", "leo_recovered_pc"), ("leo_pc_or_fp", "leo_recovered_pc_or_fp")]:
        pred = q[pred_col].fillna(False).astype(bool)
        bls = q["bls_recovered"]
        tp = int((pred & bls).sum())
        fp = int((pred & ~bls).sum())
        fn = int((~pred & bls).sum())
        tn = int((~pred & ~bls).sum())
        rows[label] = {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "tn": tn,
            "precision_vs_bls": tp / (tp + fp) if tp + fp else float("nan"),
            "recall_vs_bls": tp / (tp + fn) if tp + fn else float("nan"),
            "specificity_vs_bls": tn / (tn + fp) if tn + fp else float("nan"),
        }
    return table, rows


def plot_leo_bls_contingency(queue_csv: Path, out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    table, summary = leo_bls_contingency(queue_csv)
    core = table.drop(index="All", errors="ignore").drop(columns="All", errors="ignore")
    classes = [c for c in ["FA", "FP", "PC"] if c in core.columns]
    core = core[classes]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), gridspec_kw={"width_ratios": [1.15, 1.0]})
    ax = axes[0]
    image = ax.imshow(core.to_numpy(dtype=float), cmap="Blues")
    ax.set_xticks(np.arange(len(core.columns)), labels=core.columns)
    ax.set_yticks(np.arange(len(core.index)), labels=["BLS missed", "BLS recovered"])
    for i in range(core.shape[0]):
        for j in range(core.shape[1]):
            value = int(core.iloc[i, j])
            ax.text(j, i, str(value), ha="center", va="center", color="black", fontsize=11)
    ax.set_title("Injected rows: BLS vs LEO class")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="count")

    ax = axes[1]
    metrics = summary["leo_pc_or_fp"]
    bars = [
        metrics["precision_vs_bls"],
        metrics["recall_vs_bls"],
        metrics["specificity_vs_bls"],
    ]
    labels = ["precision", "recall", "specificity"]
    ax.bar(labels, bars, color=["#2ca02c", "#ff7f0e", "#1f77b4"])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("LEO PC/FP vs BLS match")
    ax.set_title("LEO PC/FP is high precision, low recall", pad=10)
    for idx, value in enumerate(bars):
        ax.text(idx, min(value + 0.035, 1.05), f"{value:.0%}", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, 1.12)
    fig.subplots_adjust(bottom=0.22, wspace=0.34)
    fig.text(
        0.5,
        0.005,
        "Here 'BLS recovered' means exact/top-N/harmonic period match; 'LEO recovered' means PC or FP rather than FA.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "bls_vs_leo_recovery_contingency.png"
    pdf = out_dir / "bls_vs_leo_recovery_contingency.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    table.to_csv(out_dir / "bls_vs_leo_recovery_contingency.csv")
    return {"bls_vs_leo_png": str(png), "bls_vs_leo_pdf": str(pdf)}


def merged_leo_metric_frame(queue_csv: Path, metrics_csv: Path) -> pd.DataFrame:
    queue = classify_recovery(pd.read_csv(queue_csv))
    if "source_kind" in queue:
        queue = queue[queue["source_kind"].astype(str).eq("injection_recovery")].copy()
    metrics = pd.read_csv(metrics_csv)
    if "source_kind" in metrics:
        metrics = metrics[metrics["source_kind"].astype(str).eq("injection_recovery")].copy()
    keep_cols = [
        "review_id",
        "tic",
        "source_kind",
        "per",
        "dur",
        "dep",
        "MES",
        "new_MES",
        "SHP",
        "N_transit",
        "new_N_transit",
        "n_in",
        "q",
        "FA1",
        "FA2",
        "leo_class",
        "leo_FA",
        "leo_FP",
        "leo_PC",
    ]
    keep_cols = [col for col in keep_cols if col in metrics.columns]
    metrics = metrics[keep_cols].drop_duplicates("review_id")
    merged = queue.merge(metrics, on="review_id", how="left", suffixes=("", "_metric"))
    leo_class = merged.get("leo_class_metric", merged.get("leo_class"))
    if leo_class is None:
        leo_class = pd.Series("", index=merged.index)
    merged["leo_class_effective"] = leo_class.fillna(merged.get("leo_class", "")).astype(str)
    merged["bls_recovered"] = merged["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
    merged["leo_recovered_pc_or_fp"] = merged["leo_class_effective"].isin(["PC", "FP"])
    merged["leo_diagnostic_group"] = "BLS missed, LEO FA"
    merged.loc[merged["bls_recovered"] & ~merged["leo_recovered_pc_or_fp"], "leo_diagnostic_group"] = (
        "BLS recovered, LEO FA"
    )
    merged.loc[merged["bls_recovered"] & merged["leo_recovered_pc_or_fp"], "leo_diagnostic_group"] = (
        "BLS recovered, LEO PC/FP"
    )
    return merged


def summarize_leo_metric_diagnostics(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["MES", "new_MES", "SHP", "N_transit", "new_N_transit", "n_in", "q", "dep", "FA1", "FA2"]
    rows: list[dict[str, Any]] = []
    for group, group_df in df.groupby("leo_diagnostic_group", dropna=False):
        for metric in metric_cols:
            if metric not in group_df:
                continue
            values = pd.to_numeric(group_df[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            if len(values) == 0:
                continue
            rows.append(
                {
                    "group": group,
                    "metric": metric,
                    "n": int(len(values)),
                    "median": float(values.median()),
                    "p16": float(values.quantile(0.16)),
                    "p84": float(values.quantile(0.84)),
                    "mean": float(values.mean()),
                }
            )
    return pd.DataFrame(rows)


def plot_leo_metric_diagnostics(queue_csv: Path, metrics_csv: Path, out_dir: Path) -> tuple[dict[str, str], dict[str, Any]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    df = merged_leo_metric_frame(queue_csv, metrics_csv)
    summary = summarize_leo_metric_diagnostics(df)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = out_dir / "leo_bls_metric_summary.csv"
    summary.to_csv(summary_csv, index=False)

    groups = ["BLS missed, LEO FA", "BLS recovered, LEO FA", "BLS recovered, LEO PC/FP"]
    labels = ["BLS missed\nLEO FA", "BLS recovered\nLEO FA", "BLS recovered\nLEO PC/FP"]
    metrics = [
        ("MES", "MES"),
        ("new_MES", "new MES"),
        ("SHP", "SHP"),
        ("N_transit", "N transit"),
        ("new_N_transit", "new N transit"),
        ("n_in", "in-transit points"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.2))
    for ax, (metric, label) in zip(axes.ravel(), metrics):
        box_values = []
        for group in groups:
            values = pd.to_numeric(df.loc[df["leo_diagnostic_group"].eq(group), metric], errors="coerce")
            values = values.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
            box_values.append(values)
        ax.boxplot(
            box_values,
            tick_labels=labels,
            showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": "#dbeafe", "edgecolor": "#1f2937", "linewidth": 0.9},
            medianprops={"color": "#b91c1c", "linewidth": 1.4},
            whiskerprops={"color": "#1f2937", "linewidth": 0.9},
            capprops={"color": "#1f2937", "linewidth": 0.9},
        )
        ax.set_title(label)
        ax.tick_params(axis="x", labelrotation=0, labelsize=7.5)
        ax.grid(True, axis="y", color="0.90", linewidth=0.7)
        if metric in {"MES", "new_MES", "N_transit", "new_N_transit", "n_in"}:
            ymax = max((np.nanpercentile(values, 95) if len(values) else 0.0) for values in box_values)
            if np.isfinite(ymax) and ymax > 0:
                ax.set_ylim(bottom=0, top=ymax * 1.20)
    counts = df["leo_diagnostic_group"].value_counts().reindex(groups, fill_value=0)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.20, wspace=0.30, hspace=0.40)
    fig.text(
        0.5,
        0.955,
        "Why BLS-recovered injections still fail LEO",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.04,
        (
            f"Group counts: missed/FA={int(counts.iloc[0])}, recovered/FA={int(counts.iloc[1])}, "
            f"recovered/PC-FP={int(counts.iloc[2])}. LEO misses are mostly lower-MES, fewer-point, higher-SHP cases."
        ),
        ha="center",
        va="bottom",
        fontsize=8,
    )
    png = out_dir / "leo_bls_metric_diagnostics.png"
    pdf = out_dir / "leo_bls_metric_diagnostics.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    key_summary: dict[str, Any] = {
        "group_counts": {group: int(counts.loc[group]) for group in groups},
        "summary_csv": str(summary_csv),
    }
    for metric in ["MES", "new_MES", "SHP", "N_transit", "n_in"]:
        rows = summary[summary["metric"].eq(metric)]
        key_summary[metric] = {
            str(row["group"]): {
                "median": float(row["median"]),
                "p16": float(row["p16"]),
                "p84": float(row["p84"]),
                "n": int(row["n"]),
            }
            for _, row in rows.iterrows()
        }
    return (
        {
            "leo_metric_diagnostics_png": str(png),
            "leo_metric_diagnostics_pdf": str(pdf),
            "leo_metric_summary_csv": str(summary_csv),
        },
        key_summary,
    )


def write_summary(
    out_dir: Path,
    *,
    bls_csv: Path,
    leo_queue_csv: Path,
    model: dict[str, Any],
    leo_summary: dict[str, Any],
    leo_metric_summary: dict[str, Any] | None,
    paths: dict[str, str],
) -> None:
    raw_coef = np.asarray(model["raw_coef"], dtype=float)
    feature_names = model["feature_names"]
    lines = [
        "# S56 Duration-Aware Recovery Visualizations",
        "",
        f"BLS input: `{bls_csv}`",
        f"LEO comparison input: `{leo_queue_csv}`",
        "",
        "## Logistic Boundary Model",
        "",
        f"- Rows: `{model['n_rows']}`",
        f"- Recovered exact/top-N/harmonic: `{model['n_recovered']}`",
        f"- AUC: `{model['auc']:.3f}`",
        f"- Log loss: `{model['logloss']:.3f}`",
        f"- Converged: `{model['converged']}` in `{model['iterations']}` iterations",
        "",
        "| feature | raw coefficient |",
        "| --- | ---: |",
    ]
    lines.extend(f"| `{name}` | `{coef:.4g}` |" for name, coef in zip(feature_names, raw_coef))
    lines.extend(
        [
            "",
            "The fitted 50% boundary uses a physically constrained BLS proxy, `R_p^2 * sqrt(duration / period)`, plus `Tmag`. This gives a monotonic radius cutoff in period-duration space and avoids over-interpreting the correlated period-duration-radius sampling as independent physics.",
            "",
            "## BLS vs LEO",
            "",
            f"- Injected rows in LEO queue: `{leo_summary['n']}`",
            f"- BLS exact/top-N/harmonic recovered: `{leo_summary['bls_any_recovered_n']}/{leo_summary['n']}` = `{leo_summary['bls_any_recovered_frac']:.1%}`",
            f"- LEO PC/FP precision relative to BLS recovery: `{leo_summary['leo_pc_or_fp']['precision_vs_bls']:.1%}`",
            f"- LEO PC/FP recall relative to BLS recovery: `{leo_summary['leo_pc_or_fp']['recall_vs_bls']:.1%}`",
            "",
            "LEO is therefore a high-confidence vetting subset of the BLS recoveries, not yet a near-complete proxy for recovery. BLS recovery is the better statement for whether the injected signal is tracked by the search; LEO recovery is a stricter downstream pass/fail.",
            "",
            "The metric diagnostic supports targeted WD retuning, not a blind loosening of every LEO cut. BLS-recovered LEO-FA rows sit between the BLS misses and the LEO PC/FP rows in MES/new-MES, in-transit cadence count, and shape diagnostics.",
            "",
            "## Artifacts",
            "",
        ]
    )
    lines.extend(f"- `{key}`: `{path}`" for key, path in paths.items())
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "bls_csv": str(bls_csv),
        "leo_queue_csv": str(leo_queue_csv),
        "model": {k: v for k, v in model.items() if k not in {"beta_standardized", "center", "scale", "raw_coef"}},
        "model_arrays": {
            "beta_standardized": np.asarray(model["beta_standardized"]).tolist(),
            "center": np.asarray(model["center"]).tolist(),
            "scale": np.asarray(model["scale"]).tolist(),
            "raw_coef": np.asarray(model["raw_coef"]).tolist(),
        },
        "leo_summary": leo_summary,
        "leo_metric_summary": leo_metric_summary,
        "paths": paths,
    }
    (out_dir / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bls-csv", type=Path, default=DEFAULT_BLS_CSV)
    parser.add_argument("--leo-queue-csv", type=Path, default=DEFAULT_LEO_QUEUE_CSV)
    parser.add_argument("--leo-metrics-csv", type=Path, default=DEFAULT_LEO_METRICS_CSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--ridge", type=float, default=1.0)
    parser.add_argument(
        "--write-3d",
        action="store_true",
        help="Also write the legacy perspective 3D boundary plot. Disabled by default because the faceted views are clearer.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = prepare_bls_frame(args.bls_csv)
    df.to_csv(args.out_dir / "duration_model_input.csv", index=False)
    model = fit_logistic_boundary(df, ridge=args.ridge)
    paths: dict[str, str] = {}
    paths.update(plot_period_radius_boundary_lines(df, model, args.out_dir))
    paths.update(plot_period_radius_recovery_maps(df, args.out_dir))
    paths.update(plot_snr_proxy_recovery_map(df, args.out_dir))
    paths.update(plot_boundary_floor_maps(df, model, args.out_dir))
    if args.write_3d:
        paths.update(plot_boundary_3d_planes(df, model, args.out_dir))
    table, leo_summary = leo_bls_contingency(args.leo_queue_csv)
    paths.update(plot_leo_bls_contingency(args.leo_queue_csv, args.out_dir))
    leo_metric_summary: dict[str, Any] | None = None
    if args.leo_metrics_csv.exists():
        metric_paths, leo_metric_summary = plot_leo_metric_diagnostics(
            args.leo_queue_csv,
            args.leo_metrics_csv,
            args.out_dir,
        )
        paths.update(metric_paths)
    write_summary(
        args.out_dir,
        bls_csv=args.bls_csv,
        leo_queue_csv=args.leo_queue_csv,
        model=model,
        leo_summary=leo_summary,
        leo_metric_summary=leo_metric_summary,
        paths=paths,
    )
    print(json.dumps({"out_dir": str(args.out_dir), "n": int(len(df)), "paths": paths}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
