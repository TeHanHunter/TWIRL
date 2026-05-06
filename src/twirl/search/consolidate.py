"""Cross-aperture candidate consolidation.

Reads a per-sector candidates.parquet (one row per BLS peak, keyed by
TIC × aperture × peak_rank) and groups peaks across apertures into
clusters. A "cluster" is a set of peaks within `period_match_rel` in
period and within ±`t0_match_min` of mid-transit (modulo period) of
each other.

The output is one row per (TIC × cluster), with:
    - representative period, t0, duration, depth, sde from the highest-SDE
      member peak
    - n_apertures_agree: number of distinct apertures contributing peaks
      to this cluster (1, 2, or 3)
    - apertures_agree: comma-joined aperture names
    - per-aperture peak summaries (sde_<APERTURE>, period_<APERTURE>, ...)
    - best_peak_rank_in_aperture: lowest rank seen across contributors

A peak that lands in multiple apertures with consistent (P, T0) is much
more likely to be a real transit; a peak that only appears in one aperture
is almost always a per-aperture systematic (scattered light, detrending
residual, contamination quirk on that CCD).

Used downstream by candidate vetting / ranking.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


@dataclass
class ConsolidateConfig:
    period_match_rel: float = 0.005        # ±0.5% in period
    t0_match_min: float = 5.0              # ±5 minutes in mid-transit phase
    sde_min: float = 7.0                   # ignore peaks below this SDE
    apertures: tuple[str, ...] = ("DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG")


def _t0_phase_offset_min(t0a: float, t0b: float, period_d: float) -> float:
    """Minimum |t0a - (t0b + k * period)| for any integer k, in minutes."""
    dt = (t0a - t0b) % period_d
    if dt > period_d / 2:
        dt -= period_d
    return abs(dt) * 1440.0


def _cluster_peaks_for_target(
    target_rows: list[dict],
    cfg: ConsolidateConfig,
) -> list[dict]:
    """Cluster peaks for a single TIC.

    `target_rows` is the list of all BLS peak rows for this (TIC, sector).
    Returns one dict per cluster.
    """
    # Filter status==ok and SDE >= sde_min, and finite period.
    ok = []
    for r in target_rows:
        if r.get("status") != "ok":
            continue
        if not np.isfinite(r.get("period_d", np.nan)):
            continue
        if r.get("sde", -np.inf) < cfg.sde_min:
            continue
        ok.append(r)
    if not ok:
        return []

    # Sort by SDE descending so the strongest peak seeds each cluster.
    ok.sort(key=lambda r: r["sde"], reverse=True)

    clusters: list[list[dict]] = []
    used = [False] * len(ok)
    for i, seed in enumerate(ok):
        if used[i]:
            continue
        cluster = [seed]
        used[i] = True
        for j in range(i + 1, len(ok)):
            if used[j]:
                continue
            cand = ok[j]
            # Period match (allow harmonics 1/2 and 2 for now).
            ratio_match = False
            for mult in (1.0, 0.5, 2.0):
                target_p = seed["period_d"] * mult
                if abs(cand["period_d"] - target_p) / target_p < cfg.period_match_rel:
                    ratio_match = True
                    matched_period_d = target_p
                    break
            if not ratio_match:
                continue
            # T0 phase match against the seed's actual period (not harmonic).
            # Use seed period for phase folding because that's the cluster's
            # canonical period.
            dt = _t0_phase_offset_min(cand["t0_bjd"], seed["t0_bjd"],
                                       seed["period_d"])
            if dt > cfg.t0_match_min:
                continue
            cluster.append(cand)
            used[j] = True
        clusters.append(cluster)

    return [_summarize_cluster(c, cfg) for c in clusters]


def _summarize_cluster(cluster: list[dict], cfg: ConsolidateConfig) -> dict:
    """Reduce a list of clustered peaks to one summary row."""
    # Representative = highest-SDE member.
    cluster.sort(key=lambda r: r["sde"], reverse=True)
    rep = cluster[0]
    apertures_in = sorted({c["aperture"] for c in cluster})

    # Per-aperture best (highest SDE) within this cluster.
    by_ap: dict[str, dict] = {}
    for c in cluster:
        ap = c["aperture"]
        if ap not in by_ap or c["sde"] > by_ap[ap]["sde"]:
            by_ap[ap] = c

    out = {
        # Identity
        "tic": int(rep["tic"]),
        "sector": int(rep["sector"]),
        "cam": int(rep["cam"]),
        "ccd": int(rep["ccd"]),
        "tmag": float(rep["tmag"]),
        # Cluster representative
        "period_d": float(rep["period_d"]),
        "t0_bjd": float(rep["t0_bjd"]),
        "duration_min": float(rep["duration_min"]),
        "depth": float(rep["depth"]),
        "depth_snr": float(rep["depth_snr"]),
        "sde_max": float(rep["sde"]),
        "rep_aperture": str(rep["aperture"]),
        "rep_peak_rank": int(rep["peak_rank"]),
        # Cross-aperture coincidence
        "n_apertures_agree": len(apertures_in),
        "apertures_agree": ",".join(apertures_in),
        # Per-aperture best
        "sde_DET_FLUX_SML": float(by_ap["DET_FLUX_SML"]["sde"]) if "DET_FLUX_SML" in by_ap else float("nan"),
        "sde_DET_FLUX":     float(by_ap["DET_FLUX"]["sde"])     if "DET_FLUX"     in by_ap else float("nan"),
        "sde_DET_FLUX_LAG": float(by_ap["DET_FLUX_LAG"]["sde"]) if "DET_FLUX_LAG" in by_ap else float("nan"),
        "rank_DET_FLUX_SML": int(by_ap["DET_FLUX_SML"]["peak_rank"]) if "DET_FLUX_SML" in by_ap else 0,
        "rank_DET_FLUX":     int(by_ap["DET_FLUX"]["peak_rank"])     if "DET_FLUX"     in by_ap else 0,
        "rank_DET_FLUX_LAG": int(by_ap["DET_FLUX_LAG"]["peak_rank"]) if "DET_FLUX_LAG" in by_ap else 0,
        # Provenance
        "n_cad_kept_max": int(max(c["n_cad_kept"] for c in cluster)),
        "baseline_d": float(rep["baseline_d"]),
        "bls_run_id": str(rep["bls_run_id"]),
    }
    return out


def consolidate_candidates(
    candidates_parquet: Path,
    cfg: ConsolidateConfig | None = None,
) -> pa.Table:
    """Read a per-sector BLS candidates table and emit a consolidated table.

    The consolidated rows are sorted by `n_apertures_agree` desc, then
    `sde_max` desc. The most interesting candidates are at the top.
    """
    cfg = cfg or ConsolidateConfig()
    table = pq.read_table(candidates_parquet)
    df = table.to_pandas()

    rows: list[dict] = []
    # Group by (tic, sector) — typically a sector has one TIC per file.
    for (tic, sector), sub in df.groupby(["tic", "sector"], sort=False):
        if tic < 0:
            continue  # read_fail rows
        target_rows = sub.to_dict("records")
        rows.extend(_cluster_peaks_for_target(target_rows, cfg))

    # Sort by coincidence then SDE.
    rows.sort(key=lambda r: (-r["n_apertures_agree"], -r["sde_max"]))

    if not rows:
        return pa.table({})

    cols: dict[str, list] = {k: [] for k in rows[0].keys()}
    for r in rows:
        for k in cols:
            cols[k].append(r.get(k))
    return pa.table(cols)


__all__ = ["ConsolidateConfig", "consolidate_candidates"]
