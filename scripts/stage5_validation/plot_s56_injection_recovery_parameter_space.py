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
    paths.update(plot_radius_map(df, args.out_dir))
    write_summary(df, binned, args.out_dir, paths, args.input_csv)
    print(json.dumps({"out_dir": str(args.out_dir), "n": int(len(df)), "paths": paths}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
