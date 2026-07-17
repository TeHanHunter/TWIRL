"""Two-aperture TWIRL vet sheets for S56 human triage."""
from __future__ import annotations

from dataclasses import asdict
import os
from pathlib import Path
import tempfile
from typing import Any

import numpy as np

_MPL_CACHE = Path(tempfile.gettempdir()) / "twirl_mplconfig"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

from twirl.io.hlsp import BJDREFI, HLSPLightCurve, quality_mask
from twirl.lightcurves.detrend_presets import (
    ADP015Q_COLUMN_TAG,
    TWIRL_FS_V2_ADP015Q_BRANCH,
    compare_column_names,
)
from twirl.plotting.style import apply_twirl_style
from twirl.search.bls import BLSConfig, run_bls_on_lc
from twirl.vetting.adpplus import (
    ADP_VET_APERTURES,
    ADPPlusBranch,
    binned_trend_ptp,
    branch_by_name,
    branched_light_curve,
)


_ADP015_COLUMNS = compare_column_names(ADP015Q_COLUMN_TAG)
DEFAULT_TWO_APERTURE_BRANCH = TWIRL_FS_V2_ADP015Q_BRANCH
DEFAULT_TWO_APERTURE_APERTURES: tuple[str, str] = (
    _ADP015_COLUMNS["small"],
    _ADP015_COLUMNS["primary"],
)


def _aperture_prefix(aperture: str) -> str:
    return aperture.replace("DET_FLUX_", "").lower()


def _phase_min(time_btjd: np.ndarray, *, period_d: float, t0_bjd: float) -> np.ndarray:
    t0_d = float(t0_bjd) - float(BJDREFI)
    phase_d = ((time_btjd - t0_d + 0.5 * period_d) % period_d) - 0.5 * period_d
    return phase_d * 1440.0


def _bin_xy(x: np.ndarray, y: np.ndarray, *, n_bins: int = 80) -> tuple[np.ndarray, np.ndarray]:
    finite = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(finite) < 3:
        return np.array([], dtype=float), np.array([], dtype=float)
    xg = x[finite]
    yg = y[finite]
    edges = np.linspace(float(np.nanmin(xg)), float(np.nanmax(xg)), n_bins + 1)
    if not np.all(np.isfinite(edges)) or np.unique(edges).size < 3:
        return np.array([], dtype=float), np.array([], dtype=float)
    idx = np.digitize(xg, edges) - 1
    centers = 0.5 * (edges[:-1] + edges[1:])
    meds = np.full(n_bins, np.nan, dtype=float)
    for i in range(n_bins):
        sel = idx == i
        if np.count_nonzero(sel) >= 2:
            meds[i] = float(np.nanmedian(yg[sel]))
    ok = np.isfinite(meds)
    return centers[ok], meds[ok]


def _norm_flux(lc: HLSPLightCurve, aperture: str) -> np.ndarray:
    flux = np.asarray(lc.flux[aperture], dtype=float)
    keep = quality_mask(lc, aperture)
    med = float(np.nanmedian(flux[keep])) if np.any(keep) else float(np.nanmedian(flux[np.isfinite(flux)]))
    if not np.isfinite(med) or med == 0:
        med = 1.0
    return flux / med


def _run_branch_bls(
    lc: HLSPLightCurve,
    *,
    aperture: str,
    branch: ADPPlusBranch,
    cfg: BLSConfig,
    period_d: float | None = None,
    t0_bjd: float | None = None,
    duration_min: float | None = None,
    return_periodogram: bool = True,
) -> tuple[Any, dict[str, np.ndarray] | None, HLSPLightCurve, dict[str, Any]]:
    branch_lc, meta = branched_light_curve(
        lc,
        aperture,
        branch,
        period_d=period_d,
        t0_bjd=t0_bjd,
        duration_min=duration_min,
    )
    result = run_bls_on_lc(branch_lc, cfg, aperture=aperture, return_periodogram=return_periodogram)
    if return_periodogram:
        res, spectrum = result
        return res, spectrum, branch_lc, dict(meta)
    return result, None, branch_lc, dict(meta)


def _best_peak_dict(result: Any) -> dict[str, float]:
    if not getattr(result, "peaks", None):
        return {
            "period_d": float("nan"),
            "t0_bjd": float("nan"),
            "duration_min": float("nan"),
            "depth": float("nan"),
            "depth_snr": float("nan"),
            "sde": float("nan"),
        }
    return asdict(result.peaks[0])


def _depth_at_ephemeris(
    lc: HLSPLightCurve,
    aperture: str,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
) -> tuple[float, float, int]:
    y = _norm_flux(lc, aperture)
    keep = quality_mask(lc, aperture)
    phase = _phase_min(lc.time, period_d=period_d, t0_bjd=t0_bjd)
    half = 0.5 * float(duration_min)
    in_tr = keep & np.isfinite(y) & (np.abs(phase) <= half)
    oot = keep & np.isfinite(y) & (np.abs(phase) > max(2.0 * half, 20.0)) & (np.abs(phase) < 90.0)
    if np.count_nonzero(in_tr) < 1 or np.count_nonzero(oot) < 10:
        return float("nan"), float("nan"), int(np.count_nonzero(in_tr))
    oot_med = float(np.nanmedian(y[oot]))
    in_med = float(np.nanmedian(y[in_tr]))
    sigma = 1.4826 * float(np.nanmedian(np.abs(y[oot] - oot_med)))
    depth = oot_med - in_med
    snr = depth / sigma * np.sqrt(np.count_nonzero(in_tr)) if sigma > 0 else float("nan")
    return float(depth), float(snr), int(np.count_nonzero(in_tr))


def build_two_aperture_metrics(
    lc: HLSPLightCurve,
    *,
    branch: ADPPlusBranch,
    cfg: BLSConfig,
    apertures: tuple[str, ...] = DEFAULT_TWO_APERTURE_APERTURES,
    anchor_aperture: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Run BLS for both apertures and return metrics plus plotting payload."""

    results: dict[str, Any] = {}
    spectra: dict[str, dict[str, np.ndarray]] = {}
    lcs: dict[str, HLSPLightCurve] = {}
    branch_meta: dict[str, Any] = {}
    for ap in apertures:
        if ap not in lc.flux:
            continue
        if branch.is_identity:
            res, spec, branch_lc, meta = _run_branch_bls(lc, aperture=ap, branch=branch, cfg=cfg)
        else:
            current = branch_by_name("current_adp")
            current_res, _, _, _ = _run_branch_bls(lc, aperture=ap, branch=current, cfg=cfg, return_periodogram=False)
            peak = _best_peak_dict(current_res)
            res, spec, branch_lc, meta = _run_branch_bls(
                lc,
                aperture=ap,
                branch=branch,
                cfg=cfg,
                period_d=peak["period_d"],
                t0_bjd=peak["t0_bjd"],
                duration_min=peak["duration_min"],
            )
        results[ap] = res
        if spec is not None:
            spectra[ap] = spec
        lcs[ap] = branch_lc
        branch_meta[ap] = meta

    if anchor_aperture and anchor_aperture in results:
        anchor_ap = anchor_aperture
    else:
        anchor_ap = apertures[0] if apertures and apertures[0] in results else next(iter(results), "")
    anchor_peak = _best_peak_dict(results[anchor_ap]) if anchor_ap else {}
    metrics: dict[str, Any] = {
        "tic": lc.tic,
        "sector": lc.sector,
        "cam": lc.cam,
        "ccd": lc.ccd,
        "tmag": lc.tmag,
        "vet_branch": branch.name,
        "anchor_aperture": anchor_ap,
        "anchor_period_d": anchor_peak.get("period_d", np.nan),
        "anchor_t0_bjd": anchor_peak.get("t0_bjd", np.nan),
        "anchor_duration_min": anchor_peak.get("duration_min", np.nan),
        "anchor_sde": anchor_peak.get("sde", np.nan),
    }
    for ap, res in results.items():
        peak = _best_peak_dict(res)
        prefix = _aperture_prefix(ap)
        metrics[f"{prefix}_status"] = res.status
        for key, value in peak.items():
            metrics[f"{prefix}_{key}"] = value
        if np.isfinite(metrics["anchor_period_d"]):
            depth, snr, n_in = _depth_at_ephemeris(
                lcs[ap],
                ap,
                period_d=float(metrics["anchor_period_d"]),
                t0_bjd=float(metrics["anchor_t0_bjd"]),
                duration_min=float(metrics["anchor_duration_min"]),
            )
            metrics[f"{prefix}_anchor_depth"] = depth
            metrics[f"{prefix}_anchor_snr"] = snr
            metrics[f"{prefix}_anchor_n_in"] = n_in
        metrics[f"{prefix}_trend_ptp"] = binned_trend_ptp(lcs[ap].time, lcs[ap].flux[ap], lcs[ap].quality)
        metrics[f"{prefix}_branch_status"] = branch_meta[ap].get("status", "")
        metrics[f"{prefix}_branch_window_d"] = branch_meta[ap].get("window_d", np.nan)

    small_prefix = _aperture_prefix(apertures[0]) if len(apertures) >= 1 else ""
    primary_prefix = _aperture_prefix(apertures[1]) if len(apertures) >= 2 else ""
    small_p = metrics.get(f"{small_prefix}_period_d", np.nan)
    primary_p = metrics.get(f"{primary_prefix}_period_d", np.nan)
    if np.isfinite(small_p) and np.isfinite(primary_p) and small_p > 0:
        metrics["aperture_period_rel_delta"] = abs(float(primary_p) - float(small_p)) / float(small_p)
    else:
        metrics["aperture_period_rel_delta"] = np.nan
    small_depth = metrics.get(f"{small_prefix}_anchor_depth", np.nan)
    primary_depth = metrics.get(f"{primary_prefix}_anchor_depth", np.nan)
    metrics["aperture_depth_ratio_primary_over_small"] = (
        float(primary_depth) / float(small_depth)
        if np.isfinite(small_depth) and abs(float(small_depth)) > 0
        else np.nan
    )
    metrics["aperture_disagreement_flag"] = bool(
        np.isfinite(metrics["aperture_period_rel_delta"]) and metrics["aperture_period_rel_delta"] > 0.02
    )
    payload = {"results": results, "spectra": spectra, "lcs": lcs, "metrics": metrics}
    return metrics, payload


def render_two_aperture_sheet(
    lc: HLSPLightCurve,
    out_path: Path,
    *,
    branch_name: str = DEFAULT_TWO_APERTURE_BRANCH,
    cfg: BLSConfig | None = None,
    apertures: tuple[str, ...] = DEFAULT_TWO_APERTURE_APERTURES,
    anchor_aperture: str | None = None,
    row_metadata: dict[str, Any] | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Render a PNG/PDF two-aperture vet sheet and return metrics."""

    import matplotlib

    mpl_cache = Path(tempfile.gettempdir()) / "twirl_mplconfig"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    branch = branch_by_name(branch_name)
    cfg = cfg or BLSConfig(apertures=apertures, n_periods=20_000, n_peaks=10)
    metrics, payload = build_two_aperture_metrics(
        lc,
        branch=branch,
        cfg=cfg,
        apertures=apertures,
        anchor_aperture=anchor_aperture,
    )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(12.0, 11.0))
    gs = fig.add_gridspec(5, 2, height_ratios=[1.05, 0.85, 1.0, 1.0, 0.62], hspace=0.36, wspace=0.24)
    apertures = [ap for ap in apertures if ap in payload["lcs"]]
    anchor_period = metrics.get("anchor_period_d", np.nan)
    anchor_t0 = metrics.get("anchor_t0_bjd", np.nan)
    anchor_dur = metrics.get("anchor_duration_min", np.nan)

    for col, ap in enumerate(apertures):
        blc: HLSPLightCurve = payload["lcs"][ap]
        y = _norm_flux(blc, ap)
        keep = quality_mask(blc, ap)
        peak = _best_peak_dict(payload["results"][ap])
        ax = fig.add_subplot(gs[0, col])
        ax.scatter(blc.time[keep], y[keep], s=2.2, c="0.15", alpha=0.60, rasterized=True)
        if np.isfinite(peak["period_d"]):
            n_lo = int(np.ceil((blc.time[keep].min() + BJDREFI - peak["t0_bjd"]) / peak["period_d"]))
            n_hi = int(np.floor((blc.time[keep].max() + BJDREFI - peak["t0_bjd"]) / peak["period_d"]))
            for n in range(n_lo, n_hi + 1):
                ax.axvline(peak["t0_bjd"] + n * peak["period_d"] - BJDREFI, color="#c4452f", lw=0.35, alpha=0.35)
        if keep.any():
            lo, hi = np.nanpercentile(y[keep], [1, 99.5])
            pad = 0.05 * max(hi - lo, 0.05)
            ax.set_ylim(lo - pad, hi + pad)
            ax.set_xlim(np.nanmin(blc.time[keep]), np.nanmax(blc.time[keep]))
        ax.set_title(f"{ap} full LC, {branch.name}", loc="left", fontsize=9)
        ax.set_xlabel(f"BJD - {BJDREFI}")
        ax.set_ylabel("relative flux")

        ax = fig.add_subplot(gs[1, col])
        spec = payload["spectra"].get(ap)
        if spec is not None:
            ax.semilogx(spec["period"], spec["sde"], color="0.2", lw=0.5, rasterized=True)
            ax.axvline(peak["period_d"], color="#c4452f", lw=0.9, ls="--")
            ax.set_xlim(float(np.nanmin(spec["period"])), float(np.nanmax(spec["period"])))
        else:
            ax.text(0.5, 0.5, "No periodogram", ha="center", va="center")
        ax.set_title(f"{ap} own peak P={peak['period_d']:.5g} d, SDE={peak['sde']:.1f}", loc="left", fontsize=9)
        ax.set_xlabel("period (d)")
        ax.set_ylabel("SDE")

        ax = fig.add_subplot(gs[2, col])
        if np.isfinite(peak["period_d"]):
            phase = _phase_min(blc.time, period_d=peak["period_d"], t0_bjd=peak["t0_bjd"])
            win = keep & (np.abs(phase) <= 90.0)
            ax.scatter(phase[win], y[win], s=4, c="0.35", alpha=0.42, rasterized=True)
            bx, by = _bin_xy(phase[win], y[win], n_bins=70)
            ax.plot(bx, by, color="#1f4e79", lw=1.4)
            half = 0.5 * peak["duration_min"]
            ax.plot([-90, -half, -half, half, half, 90], [1, 1, 1 - peak["depth"], 1 - peak["depth"], 1, 1], color="#c4452f", lw=1.2)
        ax.set_xlim(-90, 90)
        ax.set_title(f"{ap} folded at own peak", loc="left", fontsize=9)
        ax.set_xlabel("minutes from own BLS epoch")
        ax.set_ylabel("relative flux")

        ax = fig.add_subplot(gs[3, col])
        if np.isfinite(anchor_period):
            phase = _phase_min(blc.time, period_d=float(anchor_period), t0_bjd=float(anchor_t0))
            win = keep & (np.abs(phase) <= 90.0)
            ax.scatter(phase[win], y[win], s=4, c="0.35", alpha=0.42, rasterized=True)
            bx, by = _bin_xy(phase[win], y[win], n_bins=70)
            ax.plot(bx, by, color="#1f4e79", lw=1.4)
            half = 0.5 * float(anchor_dur)
            anchor_depth = metrics.get(f"{_aperture_prefix(ap)}_anchor_depth", np.nan)
            depth = float(anchor_depth) if np.isfinite(anchor_depth) and anchor_depth > 0 else peak["depth"]
            ax.plot([-90, -half, -half, half, half, 90], [1, 1, 1 - depth, 1 - depth, 1, 1], color="#c4452f", lw=1.2)
        ax.set_xlim(-90, 90)
        ax.set_title(f"{ap} folded at shared small-aperture anchor", loc="left", fontsize=9)
        ax.set_xlabel("minutes from anchor epoch")
        ax.set_ylabel("relative flux")

    ax_txt = fig.add_subplot(gs[4, :])
    ax_txt.axis("off")
    lines = [
        f"TIC {lc.tic}  S{lc.sector} cam{lc.cam}/ccd{lc.ccd}  T={lc.tmag:.2f}",
        (
            f"anchor={metrics.get('anchor_aperture', '')}  "
            f"P={metrics.get('anchor_period_d', np.nan):.6g} d  "
            f"T0={metrics.get('anchor_t0_bjd', np.nan):.6f}  "
            f"dur={metrics.get('anchor_duration_min', np.nan):.3g} min  "
            f"SDE={metrics.get('anchor_sde', np.nan):.3g}"
        ),
        (
            f"aperture period delta={metrics.get('aperture_period_rel_delta', np.nan):.3g}; "
            f"primary/small depth ratio={metrics.get('aperture_depth_ratio_primary_over_small', np.nan):.3g}; "
            f"aperture_disagreement={metrics.get('aperture_disagreement_flag')}"
        ),
    ]
    if row_metadata:
        lines.append(
            "input row: "
            + "  ".join(
                f"{k}={row_metadata.get(k)}"
                for k in ("review_id", "selection_bucket", "source_kind")
                if k in row_metadata
            )
        )
    ax_txt.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=9)
    fig.suptitle(f"TIC {lc.tic} two-aperture TWIRL vet sheet ({branch.name})", fontsize=13, y=0.99)
    fig.subplots_adjust(top=0.955, bottom=0.045, left=0.07, right=0.985)
    fig.savefig(out_path, dpi=160)
    try:
        fig.savefig(out_path.with_suffix(".pdf"))
    except Exception:
        pass
    plt.close(fig)
    metrics["twirl_vet_sheet_name"] = out_path.name
    metrics["twirl_vet_sheet_pdf_name"] = out_path.with_suffix(".pdf").name
    return out_path, metrics
