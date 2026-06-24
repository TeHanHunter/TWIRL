#!/usr/bin/env python3
"""Plot S56 injected-signal recovery in period/depth/radius space."""
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


DEFAULT_INPUT = (
    REPO_ROOT
    / "reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/"
    / "small_pair_200k/injection_bls_recoveries.csv"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/"
    / "parameter_space"
)
DEFAULT_PERIOD_BINS = 50
DEFAULT_DEPTH_BINS = 50
DEFAULT_SMOOTH_GRID_SIZE = 90
DEFAULT_SMOOTH_SIGMA_LOGP = 0.13
DEFAULT_SMOOTH_SIGMA_DEPTH_PCT = 7.0


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


def prepare_frame(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = classify_recovery(df)
    for col in (
        "truth_period_d",
        "truth_sampled_model_depth",
        "truth_model_depth",
        "truth_depth",
        "truth_radius_rearth",
        "tmag",
        "sde_max",
        "truth_n_good_in_transit",
    ):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    depth = df.get("truth_sampled_model_depth", pd.Series(np.nan, index=df.index)).copy()
    for fallback in ("truth_model_depth", "truth_depth"):
        missing = ~np.isfinite(depth) | (depth <= 0)
        if missing.any() and fallback in df:
            depth.loc[missing] = pd.to_numeric(df.loc[missing, fallback], errors="coerce")
    df["plot_depth_frac"] = depth.clip(lower=0, upper=1.0)
    df["plot_depth_pct"] = 100.0 * df["plot_depth_frac"]
    df = df[np.isfinite(df["truth_period_d"]) & np.isfinite(df["plot_depth_pct"]) & np.isfinite(df["tmag"])]
    return df


def binned_summary(df: pd.DataFrame) -> pd.DataFrame:
    tmag_bins = [0, 16, 17, 18, 19, 20, 30]
    period_bins = [0.10, 0.216, 0.50, 1.00, 2.00, 5.00, 15.0]
    depth_bins = [0, 5, 10, 20, 40, 60, 80, 100]
    rows: list[dict[str, Any]] = []
    for tlo, thi in zip(tmag_bins[:-1], tmag_bins[1:]):
        tmask = (df["tmag"] >= tlo) & (df["tmag"] < thi)
        sub = df[tmask]
        rows.append(_summarize_subset(sub, f"Tmag [{tlo:g},{thi:g})", "tmag_bin", tlo, thi))
    for plo, phi in zip(period_bins[:-1], period_bins[1:]):
        pmask = (df["truth_period_d"] >= plo) & (df["truth_period_d"] < phi)
        rows.append(_summarize_subset(df[pmask], f"P [{plo:g},{phi:g}) d", "period_bin", plo, phi))
    for dlo, dhi in zip(depth_bins[:-1], depth_bins[1:]):
        dmask = (df["plot_depth_pct"] >= dlo) & (df["plot_depth_pct"] < dhi)
        rows.append(_summarize_subset(df[dmask], f"depth [{dlo:g},{dhi:g}) pct", "depth_bin", dlo, dhi))
    return pd.DataFrame(rows)


def _summarize_subset(sub: pd.DataFrame, label: str, axis: str, lo: float, hi: float) -> dict[str, Any]:
    n = int(len(sub))
    strict = int(sub["strict_top1_recovered"].sum()) if n else 0
    any_match = int(sub["any_exact_or_harmonic_recovered"].sum()) if n else 0
    return {
        "axis": axis,
        "label": label,
        "lo": lo,
        "hi": hi,
        "n": n,
        "strict_top1_n": strict,
        "strict_top1_frac": strict / n if n else np.nan,
        "any_exact_or_harmonic_n": any_match,
        "any_exact_or_harmonic_frac": any_match / n if n else np.nan,
        "median_tmag": float(np.nanmedian(sub["tmag"])) if n else np.nan,
        "median_period_d": float(np.nanmedian(sub["truth_period_d"])) if n else np.nan,
        "median_depth_pct": float(np.nanmedian(sub["plot_depth_pct"])) if n else np.nan,
        "median_radius_rearth": float(np.nanmedian(sub["truth_radius_rearth"])) if n and "truth_radius_rearth" in sub else np.nan,
    }


def plot_recovery_map(df: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.1), sharex=True, sharey=True)
    panels = [
        ("strict_top1_recovered", "Top-rank BLS recovery"),
        ("any_exact_or_harmonic_recovered", "Exact/top-N/harmonic period match"),
    ]
    norm = plt.Normalize(vmin=max(14.0, float(np.nanpercentile(df["tmag"], 1))), vmax=min(20.5, float(np.nanpercentile(df["tmag"], 99))))
    cmap = "viridis_r"
    dense = len(df) >= 3000
    missed_size = 7 if dense else 17
    recovered_size = 13 if dense else 33
    missed_alpha = 0.17 if dense else 0.24
    recovered_alpha = 0.72 if dense else 0.92
    scatter = None
    for ax, (col, title) in zip(axes, panels):
        recovered = df[col].fillna(False).astype(bool)
        ax.scatter(
            df.loc[~recovered, "truth_period_d"],
            df.loc[~recovered, "plot_depth_pct"],
            c=df.loc[~recovered, "tmag"],
            cmap=cmap,
            norm=norm,
            s=missed_size,
            alpha=missed_alpha,
            marker="x",
            linewidths=0.35,
            label="missed",
        )
        scatter = ax.scatter(
            df.loc[recovered, "truth_period_d"],
            df.loc[recovered, "plot_depth_pct"],
            c=df.loc[recovered, "tmag"],
            cmap=cmap,
            norm=norm,
            s=recovered_size,
            alpha=recovered_alpha,
            edgecolors="black",
            linewidths=0.15,
            label="recovered",
        )
        ax.axvline(0.216, color="black", linestyle=":", linewidth=1.4)
        ax.text(0.224, 8, "Roche limit", rotation=90, va="bottom", ha="left", fontsize=8)
        ax.axvspan(0.4, 1.3, facecolor="0.4", alpha=0.08, hatch="//", edgecolor="0.35", linewidth=0.0)
        ax.text(0.64, 93, "Typical HZ", ha="center", va="center", fontsize=8, color="0.2")
        ax.set_title(title)
        ax.set_xscale("log")
        ax.set_xlim(0.10, 15.5)
        ax.set_ylim(-1, 102)
        ax.grid(True, which="major", alpha=0.28)
        ax.grid(True, which="minor", axis="x", alpha=0.12)
        n = len(df)
        rec_n = int(recovered.sum())
        ax.text(
            0.97,
            0.04,
            f"{rec_n}/{n} = {rec_n / n:.1%}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.9, "boxstyle": "round,pad=0.25"},
        )
    axes[0].set_ylabel("Finite-exposure transit depth [%]")
    for ax in axes:
        ax.set_xlabel("Injected orbital period [days]")
    if scatter is not None:
        cbar = fig.colorbar(scatter, ax=axes, fraction=0.035, pad=0.02)
        cbar.set_label("TESS magnitude")
    fig.suptitle("S56 pre-detrend BATMAN injection recovery, small-aperture BLS", y=0.99)
    fig.text(
        0.5,
        0.005,
        "Transparent crosses are misses; filled circles are recovered under each panel's definition. "
        "BLS uses DET_FLUX_ADP_SML + DET_FLUX_SML, 200k period grid.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "period_depth_recovery_by_tmag.png"
    pdf = out_dir / "period_depth_recovery_by_tmag.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"period_depth_png": str(png), "period_depth_pdf": str(pdf)}


def _period_depth_edges(
    df: pd.DataFrame,
    *,
    period_bins: int,
    depth_bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    finite_period = df["truth_period_d"].to_numpy(dtype=float)
    finite_period = finite_period[np.isfinite(finite_period) & (finite_period > 0)]
    if finite_period.size == 0:
        raise ValueError("no finite periods available for binned plot")
    lo = max(0.10, float(np.nanpercentile(finite_period, 0.1)))
    hi = min(15.0, float(np.nanpercentile(finite_period, 99.9)))
    if hi <= lo:
        lo, hi = 0.12, 13.0
    period_edges = np.geomspace(lo, hi, int(period_bins) + 1)
    depth_edges = np.linspace(0.0, 100.0, int(depth_bins) + 1)
    return period_edges, depth_edges


def _recovery_fraction_grid(
    df: pd.DataFrame,
    recovery_col: str,
    *,
    period_edges: np.ndarray,
    depth_edges: np.ndarray,
    min_count: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    period = df["truth_period_d"].to_numpy(dtype=float)
    depth = df["plot_depth_pct"].to_numpy(dtype=float)
    recovered = df[recovery_col].fillna(False).astype(bool).to_numpy(dtype=float)
    count, _, _ = np.histogram2d(period, depth, bins=(period_edges, depth_edges))
    rec, _, _ = np.histogram2d(period, depth, bins=(period_edges, depth_edges), weights=recovered)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = rec / count
    frac[count < int(min_count)] = np.nan
    return frac, count


def _draw_fraction_panel(
    ax,
    df: pd.DataFrame,
    recovery_col: str,
    *,
    title: str,
    period_edges: np.ndarray,
    depth_edges: np.ndarray,
    min_count: int,
):
    import matplotlib.pyplot as plt

    frac, count = _recovery_fraction_grid(
        df,
        recovery_col,
        period_edges=period_edges,
        depth_edges=depth_edges,
        min_count=min_count,
    )
    mesh = ax.pcolormesh(period_edges, depth_edges, frac.T, vmin=0.0, vmax=1.0, cmap="magma", shading="auto")
    x_centers = np.sqrt(period_edges[:-1] * period_edges[1:])
    y_centers = 0.5 * (depth_edges[:-1] + depth_edges[1:])
    finite = np.isfinite(frac)
    if np.count_nonzero(finite) >= 9:
        levels = [level for level in (0.5, 0.8) if np.nanmin(frac) <= level <= np.nanmax(frac)]
        if levels:
            contours = ax.contour(
                x_centers,
                y_centers,
                frac.T,
                levels=levels,
                colors=["white", "cyan"][: len(levels)],
                linewidths=1.1,
            )
            ax.clabel(contours, fmt={0.5: "50%", 0.8: "80%"}, fontsize=7)
    ax.axvline(0.216, color="white", linestyle=":", linewidth=1.2)
    ax.axvspan(0.4, 1.3, facecolor="white", alpha=0.08, hatch="//", edgecolor="white", linewidth=0.0)
    ax.set_title(title)
    ax.set_xscale("log")
    ax.set_xlim(float(period_edges[0]), float(period_edges[-1]))
    ax.set_ylim(float(depth_edges[0]), float(depth_edges[-1]))
    ax.grid(False)
    ax.text(
        0.98,
        0.04,
        f"n={len(df):,}; bins >= {min_count}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color="white",
        bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.45, "boxstyle": "round,pad=0.2"},
    )
    return mesh, count


def plot_recovery_fraction_maps(
    df: pd.DataFrame,
    out_dir: Path,
    *,
    period_bins: int,
    depth_bins: int,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    period_edges, depth_edges = _period_depth_edges(df, period_bins=period_bins, depth_bins=depth_bins)
    min_count = 2 if len(df) >= 3000 else 1

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.1), sharex=True, sharey=True)
    panels = [
        ("strict_top1_recovered", "Top-rank BLS recovery fraction"),
        ("any_exact_or_harmonic_recovered", "Exact/top-N/harmonic recovery fraction"),
    ]
    mesh = None
    for ax, (col, title) in zip(axes, panels):
        mesh, _ = _draw_fraction_panel(
            ax,
            df,
            col,
            title=title,
            period_edges=period_edges,
            depth_edges=depth_edges,
            min_count=min_count,
        )
        ax.set_xlabel("Injected orbital period [days]")
    axes[0].set_ylabel("Finite-exposure transit depth [%]")
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes, fraction=0.035, pad=0.02)
        cbar.set_label("Recovered fraction")
    fig.suptitle("S56 pre-detrend BATMAN injection BLS sensitivity", y=0.99)
    fig.text(
        0.5,
        0.005,
        "Color shows recovered fraction in period-depth bins; contours mark 50% and 80% recovery where sampled.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "period_depth_recovery_fraction_grid.png"
    pdf = out_dir / "period_depth_recovery_fraction_grid.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    return {"period_depth_fraction_png": str(png), "period_depth_fraction_pdf": str(pdf)}


def plot_tmag_sliced_recovery_fraction(
    df: pd.DataFrame,
    out_dir: Path,
    *,
    period_bins: int,
    depth_bins: int,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    period_edges, depth_edges = _period_depth_edges(df, period_bins=period_bins, depth_bins=depth_bins)
    slices = [
        ("Tmag < 18", df["tmag"] < 18.0),
        ("18 <= Tmag < 19", (df["tmag"] >= 18.0) & (df["tmag"] < 19.0)),
        ("19 <= Tmag < 20", (df["tmag"] >= 19.0) & (df["tmag"] < 20.0)),
        ("Tmag >= 20", df["tmag"] >= 20.0),
    ]
    min_count = 1
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.0), sharex=True, sharey=True)
    mesh = None
    for ax, (title, mask) in zip(axes.ravel(), slices):
        sub = df[mask].copy()
        if sub.empty:
            ax.set_title(f"{title}; n=0")
            ax.set_xscale("log")
            ax.set_xlim(float(period_edges[0]), float(period_edges[-1]))
            ax.set_ylim(0, 100)
            continue
        mesh, _ = _draw_fraction_panel(
            ax,
            sub,
            "any_exact_or_harmonic_recovered",
            title=title,
            period_edges=period_edges,
            depth_edges=depth_edges,
            min_count=min_count,
        )
    for ax in axes[-1, :]:
        ax.set_xlabel("Injected orbital period [days]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Finite-exposure transit depth [%]")
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), fraction=0.032, pad=0.02)
        cbar.set_label("Recovered fraction")
    fig.suptitle("Exact/top-N/harmonic BLS recovery by TESS magnitude", y=0.995)
    png = out_dir / "period_depth_recovery_fraction_by_tmag.png"
    pdf = out_dir / "period_depth_recovery_fraction_by_tmag.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"period_depth_fraction_tmag_png": str(png), "period_depth_fraction_tmag_pdf": str(pdf)}


def _tmag_slices(df: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    return [
        ("Tmag < 18", df["tmag"] < 18.0),
        ("18 <= Tmag < 19", (df["tmag"] >= 18.0) & (df["tmag"] < 19.0)),
        ("19 <= Tmag < 20", (df["tmag"] >= 19.0) & (df["tmag"] < 20.0)),
        ("Tmag >= 20", df["tmag"] >= 20.0),
    ]


def _smoothed_recovery_surface(
    sub: pd.DataFrame,
    recovery_col: str,
    *,
    period_grid: np.ndarray,
    depth_grid: np.ndarray,
    sigma_logp: float,
    sigma_depth_pct: float,
    batch_size: int = 384,
) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian-kernel estimate of recovery fraction on a period-depth grid."""

    x = np.log10(sub["truth_period_d"].to_numpy(dtype=float))
    y = sub["plot_depth_pct"].to_numpy(dtype=float)
    recovered = sub[recovery_col].fillna(False).astype(bool).to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    recovered = recovered[finite]
    if x.size == 0:
        shape = (len(period_grid), len(depth_grid))
        return np.full(shape, np.nan), np.zeros(shape)

    gx, gy = np.meshgrid(np.log10(period_grid), depth_grid, indexing="ij")
    points = np.column_stack([gx.ravel(), gy.ravel()])
    frac = np.full(points.shape[0], np.nan)
    effective_n = np.zeros(points.shape[0])
    inv_logp = 1.0 / max(float(sigma_logp), 1e-6)
    inv_depth = 1.0 / max(float(sigma_depth_pct), 1e-6)

    for start in range(0, len(points), batch_size):
        stop = min(start + batch_size, len(points))
        dx = (points[start:stop, 0, None] - x[None, :]) * inv_logp
        dy = (points[start:stop, 1, None] - y[None, :]) * inv_depth
        weights = np.exp(-0.5 * (dx * dx + dy * dy))
        denom = weights.sum(axis=1)
        numerator = weights @ recovered
        with np.errstate(invalid="ignore", divide="ignore"):
            frac[start:stop] = numerator / denom
            effective_n[start:stop] = denom * denom / np.square(weights).sum(axis=1)

    shape = (len(period_grid), len(depth_grid))
    return frac.reshape(shape), effective_n.reshape(shape)


def _smoothed_boundary_rows(
    surface: np.ndarray,
    effective_n: np.ndarray,
    *,
    label: str,
    period_grid: np.ndarray,
    depth_grid: np.ndarray,
    threshold: float,
    min_effective_n: float,
    n_slice: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for p_idx, period_d in enumerate(period_grid):
        frac_col = surface[p_idx, :].copy()
        eff_col = effective_n[p_idx, :]
        frac_col[eff_col < min_effective_n] = np.nan
        valid = np.isfinite(frac_col)
        if valid.any():
            monotonic = np.maximum.accumulate(np.where(valid, frac_col, -np.inf))
            above = np.where(monotonic >= threshold)[0]
        else:
            above = np.array([], dtype=int)
        if above.size:
            idx = int(above[0])
            if idx == 0 or not np.isfinite(monotonic[idx - 1]):
                boundary_depth = float(depth_grid[idx])
            else:
                low_y = float(depth_grid[idx - 1])
                high_y = float(depth_grid[idx])
                low_f = float(monotonic[idx - 1])
                high_f = float(monotonic[idx])
                if high_f > low_f:
                    boundary_depth = low_y + (threshold - low_f) * (high_y - low_y) / (high_f - low_f)
                else:
                    boundary_depth = high_y
            recovery_frac = float(frac_col[idx])
            eff_n = float(eff_col[idx])
        else:
            boundary_depth = float("nan")
            recovery_frac = float("nan")
            eff_n = float("nan")
        rows.append(
            {
                "tmag_bin": label,
                "threshold": float(threshold),
                "period_d": float(period_d),
                "boundary_depth_pct": boundary_depth,
                "smoothed_recovery_frac_at_boundary": recovery_frac,
                "effective_n_at_boundary": eff_n,
                "min_effective_n": float(min_effective_n),
                "n_slice": int(n_slice),
            }
        )
    return rows


def plot_smoothed_tmag_recovery(
    df: pd.DataFrame,
    out_dir: Path,
    *,
    grid_size: int = DEFAULT_SMOOTH_GRID_SIZE,
    sigma_logp: float = DEFAULT_SMOOTH_SIGMA_LOGP,
    sigma_depth_pct: float = DEFAULT_SMOOTH_SIGMA_DEPTH_PCT,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    period_grid = np.geomspace(0.12, 13.0, int(grid_size))
    depth_grid = np.linspace(0.0, 100.0, int(grid_size))
    slices = _tmag_slices(df)
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.0), sharex=True, sharey=True)
    mesh = None
    boundary_rows: list[dict[str, Any]] = []
    colors = {0.5: "white", 0.8: "cyan"}

    for ax, (title, mask) in zip(axes.ravel(), slices):
        sub = df[mask].copy()
        ax.set_title(f"{title}; n={len(sub):,}")
        ax.set_xscale("log")
        ax.set_xlim(float(period_grid[0]), float(period_grid[-1]))
        ax.set_ylim(float(depth_grid[0]), float(depth_grid[-1]))
        ax.axvline(0.216, color="white", linestyle=":", linewidth=1.2)
        ax.axvspan(0.4, 1.3, facecolor="white", alpha=0.08, hatch="//", edgecolor="white", linewidth=0.0)
        if sub.empty:
            continue
        surface, effective_n = _smoothed_recovery_surface(
            sub,
            "any_exact_or_harmonic_recovered",
            period_grid=period_grid,
            depth_grid=depth_grid,
            sigma_logp=sigma_logp,
            sigma_depth_pct=sigma_depth_pct,
        )
        min_effective_n = 8.0 if len(sub) >= 800 else 5.0
        masked = np.where(effective_n >= min_effective_n, surface, np.nan)
        mesh = ax.pcolormesh(period_grid, depth_grid, masked.T, vmin=0.0, vmax=1.0, cmap="magma", shading="auto")
        for threshold in (0.5, 0.8):
            boundary_rows.extend(
                _smoothed_boundary_rows(
                    surface,
                    effective_n,
                    label=title,
                    period_grid=period_grid,
                    depth_grid=depth_grid,
                    threshold=threshold,
                    min_effective_n=min_effective_n,
                    n_slice=len(sub),
                )
            )
            if np.isfinite(masked).any() and np.nanmin(masked) <= threshold <= np.nanmax(masked):
                contour = ax.contour(
                    period_grid,
                    depth_grid,
                    masked.T,
                    levels=[threshold],
                    colors=[colors[threshold]],
                    linewidths=1.25,
                )
                ax.clabel(contour, fmt={threshold: f"{threshold:.0%}"}, fontsize=7)
        ax.text(
            0.98,
            0.04,
            f"kernel: {sigma_logp:.2f} dex, {sigma_depth_pct:.0f}%",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7,
            color="white",
            bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.45, "boxstyle": "round,pad=0.2"},
        )

    for ax in axes[-1, :]:
        ax.set_xlabel("Injected orbital period [days]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Finite-exposure transit depth [%]")
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), fraction=0.032, pad=0.02)
        cbar.set_label("Kernel-smoothed recovered fraction")
    fig.suptitle("Smoothed exact/top-N/harmonic BLS recovery by TESS magnitude", y=0.995)
    fig.text(
        0.5,
        0.005,
        "Gaussian-kernel smoothing is for visualizing the sensitivity boundary; raw binned fractions remain in the companion plots.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    png = out_dir / "period_depth_smoothed_recovery_by_tmag.png"
    pdf = out_dir / "period_depth_smoothed_recovery_by_tmag.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    boundary = pd.DataFrame(boundary_rows)
    boundary_csv = out_dir / "period_depth_smoothed_boundary_by_tmag.csv"
    boundary.to_csv(boundary_csv, index=False)
    boundary_png = out_dir / "period_depth_smoothed_50pct_boundary_by_tmag.png"
    boundary_pdf = out_dir / "period_depth_smoothed_50pct_boundary_by_tmag.pdf"
    plot_smoothed_boundary(boundary, boundary_png, boundary_pdf)

    return {
        "period_depth_smoothed_tmag_png": str(png),
        "period_depth_smoothed_tmag_pdf": str(pdf),
        "period_depth_smoothed_boundary_csv": str(boundary_csv),
        "period_depth_smoothed_50pct_boundary_png": str(boundary_png),
        "period_depth_smoothed_50pct_boundary_pdf": str(boundary_pdf),
    }


def plot_smoothed_boundary(boundary: pd.DataFrame, png: Path, pdf: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {
        "Tmag < 18": "#1f77b4",
        "18 <= Tmag < 19": "#2ca02c",
        "19 <= Tmag < 20": "#ff7f0e",
        "Tmag >= 20": "#d62728",
    }
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    subset = boundary[np.isclose(boundary["threshold"], 0.5)].copy()
    for label, group in subset.groupby("tmag_bin", sort=False):
        group = group.sort_values("period_d")
        valid = np.isfinite(group["boundary_depth_pct"])
        if not valid.any():
            continue
        ax.plot(
            group.loc[valid, "period_d"],
            group.loc[valid, "boundary_depth_pct"],
            linewidth=1.8,
            color=colors.get(label),
            label=f"{label} (n={int(group['n_slice'].iloc[0])})",
        )
    ax.axvline(0.216, color="black", linestyle=":", linewidth=1.3)
    ax.text(0.224, 8, "Roche limit", rotation=90, va="bottom", ha="left", fontsize=8)
    ax.axvspan(0.4, 1.3, facecolor="0.4", alpha=0.08, hatch="//", edgecolor="0.35", linewidth=0.0)
    ax.set_xscale("log")
    ax.set_xlim(0.10, 15.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Injected orbital period [days]")
    ax.set_ylabel("Smoothed depth at 50% recovery [%]")
    ax.set_title("Kernel-smoothed BLS recovery boundary")
    ax.grid(True, which="major", alpha=0.28)
    ax.grid(True, which="minor", axis="x", alpha=0.12)
    ax.legend(fontsize=8, loc="upper left", frameon=True)
    fig.text(
        0.5,
        0.005,
        "Boundary uses a monotonic cumulative recovery curve in depth after Gaussian-kernel smoothing.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def recovery_boundary_by_tmag(
    df: pd.DataFrame,
    *,
    recovery_col: str = "any_exact_or_harmonic_recovered",
    period_bins: int = 14,
    depth_bins: int = 24,
    threshold: float = 0.5,
    min_count: int = 3,
) -> pd.DataFrame:
    period_edges, depth_edges = _period_depth_edges(df, period_bins=period_bins, depth_bins=depth_bins)
    period_centers = np.sqrt(period_edges[:-1] * period_edges[1:])
    depth_centers = 0.5 * (depth_edges[:-1] + depth_edges[1:])
    slices = [
        ("Tmag < 18", df["tmag"] < 18.0),
        ("18 <= Tmag < 19", (df["tmag"] >= 18.0) & (df["tmag"] < 19.0)),
        ("19 <= Tmag < 20", (df["tmag"] >= 19.0) & (df["tmag"] < 20.0)),
        ("Tmag >= 20", df["tmag"] >= 20.0),
    ]
    rows: list[dict[str, Any]] = []
    for label, mask in slices:
        sub = df[mask].copy()
        if sub.empty:
            continue
        frac, count = _recovery_fraction_grid(
            sub,
            recovery_col,
            period_edges=period_edges,
            depth_edges=depth_edges,
            min_count=min_count,
        )
        for p_idx, period_d in enumerate(period_centers):
            column = frac[p_idx, :]
            count_col = count[p_idx, :]
            valid = np.isfinite(column) & (count_col >= min_count)
            above = np.where(valid & (column >= threshold))[0]
            if above.size:
                first = int(above[0])
                boundary_depth = float(depth_centers[first])
                recovery_frac = float(column[first])
                n_bin = int(count_col[first])
            else:
                boundary_depth = float("nan")
                recovery_frac = float("nan")
                n_bin = 0
            rows.append(
                {
                    "tmag_bin": label,
                    "threshold": float(threshold),
                    "period_d": float(period_d),
                    "period_bin_lo_d": float(period_edges[p_idx]),
                    "period_bin_hi_d": float(period_edges[p_idx + 1]),
                    "boundary_depth_pct": boundary_depth,
                    "recovery_frac_at_boundary": recovery_frac,
                    "n_at_boundary_bin": n_bin,
                    "n_slice": int(len(sub)),
                }
            )
    return pd.DataFrame(rows)


def plot_recovery_boundary_by_tmag(
    df: pd.DataFrame,
    out_dir: Path,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    boundary = recovery_boundary_by_tmag(df)
    boundary_csv = out_dir / "period_depth_50pct_boundary_by_tmag.csv"
    boundary.to_csv(boundary_csv, index=False)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    colors = {
        "Tmag < 18": "#1f77b4",
        "18 <= Tmag < 19": "#2ca02c",
        "19 <= Tmag < 20": "#ff7f0e",
        "Tmag >= 20": "#d62728",
    }
    for label, group in boundary.groupby("tmag_bin", sort=False):
        group = group.sort_values("period_d")
        valid = np.isfinite(group["boundary_depth_pct"])
        if not valid.any():
            continue
        ax.plot(
            group.loc[valid, "period_d"],
            group.loc[valid, "boundary_depth_pct"],
            marker="o",
            markersize=4.0,
            linewidth=1.6,
            color=colors.get(label),
            label=f"{label} (n={int(group['n_slice'].iloc[0])})",
        )
    ax.axvline(0.216, color="black", linestyle=":", linewidth=1.3)
    ax.text(0.224, 8, "Roche limit", rotation=90, va="bottom", ha="left", fontsize=8)
    ax.axvspan(0.4, 1.3, facecolor="0.4", alpha=0.08, hatch="//", edgecolor="0.35", linewidth=0.0)
    ax.set_xscale("log")
    ax.set_xlim(0.10, 15.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Injected orbital period [days]")
    ax.set_ylabel("Lowest sampled depth with >=50% recovery [%]")
    ax.set_title("BLS recovery boundary by TESS magnitude")
    ax.grid(True, which="major", alpha=0.28)
    ax.grid(True, which="minor", axis="x", alpha=0.12)
    ax.legend(fontsize=8, loc="upper left", frameon=True)
    fig.text(
        0.5,
        0.005,
        "Segments are omitted where no sampled depth in that period bin reaches 50% exact/top-N/harmonic recovery.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    png = out_dir / "period_depth_50pct_boundary_by_tmag.png"
    pdf = out_dir / "period_depth_50pct_boundary_by_tmag.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {
        "period_depth_50pct_boundary_csv": str(boundary_csv),
        "period_depth_50pct_boundary_png": str(png),
        "period_depth_50pct_boundary_pdf": str(pdf),
    }


def plot_radius_map(df: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    if "truth_radius_rearth" not in df or not np.isfinite(df["truth_radius_rearth"]).any():
        return {}
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    fig, ax = plt.subplots(figsize=(6.7, 4.9))
    recovered = df["any_exact_or_harmonic_recovered"].fillna(False).astype(bool)
    ax.scatter(
        df.loc[~recovered, "truth_period_d"],
        df.loc[~recovered, "truth_radius_rearth"],
        c=df.loc[~recovered, "tmag"],
        cmap="viridis_r",
        s=17,
        alpha=0.25,
        marker="x",
        linewidths=0.55,
    )
    scatter = ax.scatter(
        df.loc[recovered, "truth_period_d"],
        df.loc[recovered, "truth_radius_rearth"],
        c=df.loc[recovered, "tmag"],
        cmap="viridis_r",
        s=34,
        alpha=0.9,
        edgecolors="black",
        linewidths=0.25,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.10, 15.5)
    ax.axvline(0.216, color="black", linestyle=":", linewidth=1.4)
    ax.axvspan(0.4, 1.3, facecolor="0.4", alpha=0.08, hatch="//", edgecolor="0.35", linewidth=0.0)
    ax.set_xlabel("Injected orbital period [days]")
    ax.set_ylabel("Injected companion radius [R_earth]")
    ax.set_title("Recovered exact/top-N/harmonic matches")
    ax.grid(True, which="major", alpha=0.28)
    ax.grid(True, which="minor", alpha=0.12)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("TESS magnitude")
    fig.tight_layout()
    png = out_dir / "period_radius_recovery_by_tmag.png"
    pdf = out_dir / "period_radius_recovery_by_tmag.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)
    return {"period_radius_png": str(png), "period_radius_pdf": str(pdf)}


def write_summary(df: pd.DataFrame, binned: pd.DataFrame, out_dir: Path, paths: dict[str, str], input_csv: Path) -> None:
    n = len(df)
    strict = int(df["strict_top1_recovered"].sum())
    any_match = int(df["any_exact_or_harmonic_recovered"].sum())
    tmag_rows = binned[binned["axis"].eq("tmag_bin")].copy()
    lines = [
        "# S56 Injection Recovery Parameter Space",
        "",
        f"Input recovery table: `{input_csv}`",
        "",
        f"- Rows plotted: `{n}`",
        f"- Strict top-rank BLS recoveries: `{strict}/{n}` = `{strict / n:.1%}`",
        f"- Exact/top-N/harmonic matches: `{any_match}/{n}` = `{any_match / n:.1%}`",
        "",
        "## Tmag Bins",
        "",
        markdown_table(
            tmag_rows[
                [
                    "label",
                    "n",
                    "strict_top1_n",
                    "strict_top1_frac",
                    "any_exact_or_harmonic_n",
                    "any_exact_or_harmonic_frac",
                    "median_period_d",
                    "median_depth_pct",
                ]
            ]
        ),
        "",
        "## Interpretation",
        "",
        "The recovery boundary should be read conditionally on Tmag and cadence coverage. A single period-depth panel mixes bright and faint WDs, so the Tmag-sliced recovery-fraction plot is the better diagnostic for the gradual BLS boundary.",
        "",
        "## Artifacts",
        "",
    ]
    for key, path in paths.items():
        lines.append(f"- `{key}`: `{path}`")
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_csv": str(input_csv),
        "n": n,
        "strict_top1_n": strict,
        "strict_top1_frac": strict / n,
        "any_exact_or_harmonic_n": any_match,
        "any_exact_or_harmonic_frac": any_match / n,
        "paths": paths,
        "tmag_bins": tmag_rows.to_dict(orient="records"),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    headers = [str(c) for c in df.columns]
    rows: list[list[str]] = []
    for _, row in df.iterrows():
        values: list[str] = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                values.append(format(float(value), ".3g") if np.isfinite(value) else "")
            else:
                values.append(str(value))
        rows.append(values)
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) for i in range(len(headers))]
    lines = [
        "| " + " | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(row))) + " |")
    return "\n".join(lines)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--period-bins", type=int, default=DEFAULT_PERIOD_BINS)
    parser.add_argument("--depth-bins", type=int, default=DEFAULT_DEPTH_BINS)
    parser.add_argument("--smooth-grid-size", type=int, default=DEFAULT_SMOOTH_GRID_SIZE)
    parser.add_argument("--smooth-sigma-logp", type=float, default=DEFAULT_SMOOTH_SIGMA_LOGP)
    parser.add_argument("--smooth-sigma-depth-pct", type=float, default=DEFAULT_SMOOTH_SIGMA_DEPTH_PCT)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    df = prepare_frame(args.input_csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_dir / "plot_input.csv", index=False)
    binned = binned_summary(df)
    binned.to_csv(args.out_dir / "recovery_binned_summary.csv", index=False)
    paths: dict[str, str] = {}
    paths.update(plot_recovery_map(df, args.out_dir))
    paths.update(
        plot_recovery_fraction_maps(
            df,
            args.out_dir,
            period_bins=args.period_bins,
            depth_bins=args.depth_bins,
        )
    )
    paths.update(
        plot_tmag_sliced_recovery_fraction(
            df,
            args.out_dir,
            period_bins=args.period_bins,
            depth_bins=args.depth_bins,
        )
    )
    paths.update(plot_recovery_boundary_by_tmag(df, args.out_dir))
    paths.update(
        plot_smoothed_tmag_recovery(
            df,
            args.out_dir,
            grid_size=args.smooth_grid_size,
            sigma_logp=args.smooth_sigma_logp,
            sigma_depth_pct=args.smooth_sigma_depth_pct,
        )
    )
    paths.update(plot_radius_map(df, args.out_dir))
    write_summary(df, binned, args.out_dir, paths, args.input_csv)
    print(json.dumps({"out_dir": str(args.out_dir), "n": int(len(df)), "paths": paths}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
