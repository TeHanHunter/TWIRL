#!/usr/bin/env python3
"""Measure injected-signal survival in pre-detrend injection HDF5 products.

This diagnostic answers whether BLS failures are plausibly low-SNR misses or
whether the injected signal is being removed/distorted before search. It reads
the stored injected-minus-original arrays, measures the signal at the truth
ephemeris, and optionally joins the review-builder BLS recovery CSV.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


DEFAULT_APERTURES = (
    "RAW_FLUX_Primary",
    "RAW_FLUX_Small",
    "RAW_FLUX_Large",
    "DET_FLUX",
    "DET_FLUX_SML",
    "DET_FLUX_LAG",
    "DET_FLUX_ADP",
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP_LAG",
)


def _robust_sigma(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan")
    med = float(np.nanmedian(finite))
    mad = float(np.nanmedian(np.abs(finite - med)))
    sigma = 1.4826 * mad
    if np.isfinite(sigma) and sigma > 0:
        return sigma
    sigma = float(np.nanstd(finite))
    return sigma if np.isfinite(sigma) else float("nan")


def _read_recovery(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "injection_id" not in df.columns:
        raise KeyError(f"recovery table has no injection_id column: {path}")
    df["injection_id"] = df["injection_id"].astype(str)
    keep = [
        c
        for c in (
            "injection_id",
            "recovery_status",
            "period_d",
            "duration_min",
            "sde_max",
            "depth_snr",
            "period_rel_err",
            "t0_phase_err_min",
            "rep_aperture",
        )
        if c in df.columns
    ]
    return df.loc[:, keep].copy()


def _float_attr(attrs: Any, key: str) -> float:
    try:
        return float(attrs.get(key, np.nan))
    except (TypeError, ValueError):
        return float("nan")


def _int_attr(attrs: Any, key: str) -> int:
    try:
        value = int(attrs.get(key, 0))
    except (TypeError, ValueError):
        return 0
    return max(value, 0)


def _measure_aperture(group: h5py.Group, aperture: str, base: dict[str, Any]) -> dict[str, Any] | None:
    original_name = f"{aperture}_original"
    injected_name = f"{aperture}_injected"
    if original_name not in group or injected_name not in group:
        return None

    original = np.asarray(group[original_name], dtype=np.float64)
    injected = np.asarray(group[injected_name], dtype=np.float64)
    in_transit = np.asarray(group["in_transit"], dtype=bool)
    quality = np.asarray(group["quality"]) if "quality" in group else np.zeros(in_transit.size, dtype=int)
    finite = np.isfinite(original) & np.isfinite(injected)
    good = finite & (quality == 0)
    in_good = good & in_transit
    oot_good = good & ~in_transit
    if np.count_nonzero(in_good) < 1 or np.count_nonzero(oot_good) < 20:
        return None

    delta = injected - original
    original_oot_median = float(np.nanmedian(original[oot_good]))
    original_oot_sigma = _robust_sigma(original[oot_good])
    delta_in_median = float(np.nanmedian(delta[in_good]))
    delta_oot_median = float(np.nanmedian(delta[oot_good]))
    delta_depth_abs = -(delta_in_median - delta_oot_median)
    denom = abs(original_oot_median) if np.isfinite(original_oot_median) and original_oot_median != 0 else np.nan
    measured_delta_depth = delta_depth_abs / denom if np.isfinite(denom) else np.nan

    truth_sampled = float(base["truth_sampled_model_depth"])
    truth_model = float(base["truth_model_depth"])
    depth_retention_sampled = (
        measured_delta_depth / truth_sampled
        if np.isfinite(measured_delta_depth) and np.isfinite(truth_sampled) and truth_sampled > 0
        else np.nan
    )
    depth_retention_model = (
        measured_delta_depth / truth_model
        if np.isfinite(measured_delta_depth) and np.isfinite(truth_model) and truth_model > 0
        else np.nan
    )
    per_cad_snr = delta_depth_abs / original_oot_sigma if np.isfinite(original_oot_sigma) and original_oot_sigma > 0 else np.nan
    multi_snr = per_cad_snr * np.sqrt(np.count_nonzero(in_good)) if np.isfinite(per_cad_snr) else np.nan

    return {
        **base,
        "aperture": aperture,
        "n_good_ap_in_transit": int(np.count_nonzero(in_good)),
        "n_good_ap_oot": int(np.count_nonzero(oot_good)),
        "original_oot_median": original_oot_median,
        "original_oot_sigma_mad": original_oot_sigma,
        "delta_in_median": delta_in_median,
        "delta_oot_median": delta_oot_median,
        "delta_depth_abs": float(delta_depth_abs),
        "measured_delta_depth": float(measured_delta_depth),
        "depth_retention_frac": float(depth_retention_sampled),
        "depth_retention_vs_model_frac": float(depth_retention_model),
        "delta_signal_snr_mad": float(per_cad_snr),
        "delta_signal_multi_snr_mad": float(multi_snr),
    }


def measure_survival(
    *,
    injection_h5: Path,
    recovery_csv: Path | None,
    apertures: tuple[str, ...],
) -> pd.DataFrame:
    recovery = _read_recovery(recovery_csv)
    recovery_by_id = recovery.set_index("injection_id", drop=False) if not recovery.empty else pd.DataFrame()
    rows: list[dict[str, Any]] = []
    with h5py.File(injection_h5, "r") as h5:
        for key in sorted(h5["injections"].keys()):
            group = h5["injections"][key]
            attrs = group.attrs
            injection_id = str(attrs.get("injection_id", key))
            in_transit = np.asarray(group["in_transit"], dtype=bool) if "in_transit" in group else np.array([], dtype=bool)
            quality = np.asarray(group["quality"]) if "quality" in group else np.zeros(in_transit.size, dtype=int)
            n_good_truth = int(np.count_nonzero(in_transit & (quality == 0))) if in_transit.size else _int_attr(attrs, "n_good_in_transit")
            base = {
                "injection_id": injection_id,
                "tic": int(attrs.get("tic", -1)),
                "sector": int(attrs.get("sector", 56)),
                "cam": int(attrs.get("camera", -1)),
                "ccd": int(attrs.get("ccd", -1)),
                "tmag": _float_attr(attrs, "tessmag"),
                "truth_period_d": _float_attr(attrs, "period_d"),
                "truth_t0_bjd": _float_attr(attrs, "t0_bjd"),
                "truth_duration_min": _float_attr(attrs, "duration_min"),
                "truth_depth": _float_attr(attrs, "depth"),
                "truth_model_depth": _float_attr(attrs, "model_depth"),
                "truth_sampled_model_depth": _float_attr(attrs, "sampled_model_depth"),
                "truth_radius_rearth": _float_attr(attrs, "radius_rearth"),
                "truth_radius_rwd": _float_attr(attrs, "radius_rwd"),
                "truth_impact_b": _float_attr(attrs, "impact_b"),
                "baseline": _float_attr(attrs, "baseline"),
                "baseline_source": str(attrs.get("baseline_source", "")),
                "expected_tess_flux_per_cadence": _float_attr(attrs, "expected_tess_flux_per_cadence"),
                "n_truth_in_transit": int(np.count_nonzero(in_transit)) if in_transit.size else _int_attr(attrs, "n_in_transit"),
                "n_good_truth_in_transit": n_good_truth,
            }
            for suffix in ("Primary", "Small", "Large"):
                base[f"raw_baseline_{suffix}"] = _float_attr(attrs, f"raw_baseline_{suffix}")
                base[f"injection_baseline_{suffix}"] = _float_attr(attrs, f"injection_baseline_{suffix}")
                base[f"inferred_baseline_{suffix}"] = _float_attr(attrs, f"inferred_baseline_{suffix}")
                base[f"aperture_fraction_{suffix}"] = _float_attr(attrs, f"aperture_fraction_{suffix}")
            if not recovery_by_id.empty and injection_id in recovery_by_id.index:
                rec = recovery_by_id.loc[injection_id]
                if isinstance(rec, pd.DataFrame):
                    rec = rec.iloc[0]
                for col in recovery.columns:
                    if col != "injection_id":
                        base[col if col == "recovery_status" else f"bls_{col}"] = rec[col]
            for aperture in apertures:
                measured = _measure_aperture(group, aperture, base)
                if measured is not None:
                    rows.append(measured)
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> dict[str, Any]:
    out: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_rows": int(len(df)),
        "n_injections": int(df["injection_id"].nunique()) if "injection_id" in df else 0,
        "apertures": sorted(df["aperture"].dropna().unique().tolist()) if "aperture" in df else [],
        "by_aperture": {},
    }
    if df.empty:
        return out
    for aperture, group in df.groupby("aperture"):
        rec: dict[str, Any] = {"n": int(len(group))}
        if "recovery_status" in group:
            recovered = group["recovery_status"].astype(str).eq("bls_recovered")
            rec["bls_recovered"] = int(recovered.sum())
            rec["bls_recovered_frac"] = float(recovered.mean())
        for col in (
            "depth_retention_frac",
            "delta_signal_snr_mad",
            "delta_signal_multi_snr_mad",
            "tmag",
            "truth_period_d",
            "truth_sampled_model_depth",
            "n_good_truth_in_transit",
        ):
            if col in group:
                finite = pd.to_numeric(group[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
                if len(finite):
                    rec[f"{col}_p16"] = float(np.nanpercentile(finite, 16))
                    rec[f"{col}_median"] = float(np.nanmedian(finite))
                    rec[f"{col}_p84"] = float(np.nanpercentile(finite, 84))
        out["by_aperture"][str(aperture)] = rec
    return out


def plot_diagnostics(df: pd.DataFrame, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    adp = df[df["aperture"].eq("DET_FLUX_ADP")].copy()
    if adp.empty:
        return
    recovered = (
        adp["recovery_status"].astype(str).eq("bls_recovered")
        if "recovery_status" in adp
        else pd.Series(False, index=adp.index)
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)

    ax = axes[0, 0]
    fail = adp[~recovered]
    succ = adp[recovered]
    sc = ax.scatter(
        fail["truth_period_d"],
        fail["truth_sampled_model_depth"],
        c=fail["tmag"],
        s=18,
        alpha=0.35,
        cmap="viridis_r",
        vmin=15,
        vmax=20.5,
        linewidths=0,
    )
    if len(succ):
        ax.scatter(
            succ["truth_period_d"],
            succ["truth_sampled_model_depth"],
            c=succ["tmag"],
            s=42,
            cmap="viridis_r",
            vmin=15,
            vmax=20.5,
            edgecolor="black",
            linewidth=0.6,
        )
    ax.set_xscale("log")
    ax.set_xlabel("Injected period [d]")
    ax.set_ylabel("BATMAN sampled depth")
    ax.set_title("Strict BLS recovery in DET_FLUX_ADP")
    fig.colorbar(sc, ax=ax, label="Tmag")

    ax = axes[0, 1]
    ax.scatter(adp["truth_period_d"], adp["depth_retention_frac"], c=adp["tmag"], s=18, alpha=0.4, cmap="viridis_r", vmin=15, vmax=20.5, linewidths=0)
    ax.axhline(1.0, color="0.2", lw=1.0, ls="--")
    ax.axhline(float(np.nanmedian(adp["depth_retention_frac"])), color="#c4452f", lw=1.2)
    ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.2)
    ax.set_xlabel("Injected period [d]")
    ax.set_ylabel("Depth retained")
    ax.set_title("Injected-minus-original depth retention")

    ax = axes[1, 0]
    ax.scatter(adp["delta_signal_multi_snr_mad"], adp["depth_retention_frac"], c=adp["tmag"], s=18, alpha=0.4, cmap="viridis_r", vmin=15, vmax=20.5, linewidths=0)
    ax.axvline(7.0, color="0.25", lw=1.0, ls="--")
    ax.set_xlabel("Truth-window multi-cadence SNR estimate")
    ax.set_ylabel("Depth retained")
    ax.set_title("Retention vs empirical injected-signal SNR")

    ax = axes[1, 1]
    bins = [-np.inf, 1, 3, 5, 7, 10, 20, np.inf]
    labels = ["<1", "1-3", "3-5", "5-7", "7-10", "10-20", ">20"]
    tmp = adp.copy()
    tmp["snr_bin"] = pd.cut(tmp["delta_signal_multi_snr_mad"], bins=bins, labels=labels)
    tab = tmp.groupby("snr_bin", observed=False).apply(
        lambda g: pd.Series(
            {
                "count": len(g),
                "recovered": int(g["recovery_status"].astype(str).eq("bls_recovered").sum()) if "recovery_status" in g else 0,
                "frac": float(g["recovery_status"].astype(str).eq("bls_recovered").mean()) if "recovery_status" in g and len(g) else np.nan,
            }
        )
    )
    x = np.arange(len(tab))
    ax.bar(x, tab["frac"].fillna(0.0), color="#4f7db8")
    for idx, row in enumerate(tab.itertuples()):
        ax.text(idx, float(row.frac if np.isfinite(row.frac) else 0.0) + 0.015, f"{int(row.recovered)}/{int(row.count)}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.set_ylim(0, max(0.2, float(np.nanmax(tab["frac"])) * 1.25 if np.isfinite(tab["frac"]).any() else 0.2))
    ax.set_xlabel("Truth-window SNR bin")
    ax.set_ylabel("Recovered fraction")
    ax.set_title("Recovery vs empirical injected-signal SNR")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--injection-h5", type=Path, required=True)
    ap.add_argument("--recovery-csv", type=Path, default=None)
    ap.add_argument("--apertures", nargs="+", default=list(DEFAULT_APERTURES))
    ap.add_argument("--out-csv", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--plot", type=Path, default=None)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    df = measure_survival(
        injection_h5=args.injection_h5,
        recovery_csv=args.recovery_csv,
        apertures=tuple(args.apertures),
    )
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    summary = summarize(df)
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if args.plot is not None:
        plot_diagnostics(df, args.plot)
    print(f"[predetrend-survival] wrote {args.out_csv} rows={len(df):,}")
    if args.out_json is not None:
        print(f"[predetrend-survival] wrote {args.out_json}")
    if args.plot is not None:
        print(f"[predetrend-survival] wrote {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
