"""Peak walk + dataclasses + Parquet schema for the per-sector BLS first pass.

A "candidate" is a row in the per-sector candidate table, keyed by
(tic, sector, aperture, peak_rank). Peaks are extracted by walking the BLS
SDE spectrum: argmax → mask the peak ± a fractional period window plus its
1/2× and 2× harmonics → repeat N times.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class BLSPeak:
    peak_rank: int
    period_d: float
    t0_bjd: float           # absolute BJD = bls.transit_time + 2457000
    duration_min: float
    depth: float            # fractional, positive = dip
    depth_snr: float
    sde: float
    log_power: float


@dataclass
class BLSResult:
    tic: int
    sector: int
    cam: int
    ccd: int
    tmag: float
    aperture: str
    n_cad_total: int
    n_cad_kept: int
    dropout_frac: float
    n_cad_quality: int | None = None
    n_cad_edge_trimmed: int = 0
    n_cad_sigma_clipped: int = 0
    quality_dropout_frac: float | None = None
    n_orbits: int = 0
    baseline_d: float = 0.0
    status: str = "ok"
    hlsp_path: str = ""
    peaks: list[BLSPeak] = field(default_factory=list)


def compute_sde(power: np.ndarray) -> np.ndarray:
    """SDE = (power - median(power)) / robust_std(power) using MAD scaling.

    Returns an array of the same shape as `power`. `power` should be the BLS
    power / log-likelihood spectrum returned by astropy. SDE >= ~7 is the
    classical detection threshold; for WD 1856 we expect SDE much higher.
    """
    med = np.nanmedian(power)
    mad = np.nanmedian(np.abs(power - med))
    if mad <= 0 or not np.isfinite(mad):
        std = np.nanstd(power)
        if std <= 0 or not np.isfinite(std):
            return np.zeros_like(power)
        return (power - med) / std
    return (power - med) / (1.4826 * mad)


def walk_peaks(
    period: np.ndarray,
    sde: np.ndarray,
    extra: dict[str, np.ndarray],
    n_peaks: int = 10,
    period_mask_frac: float = 0.005,
    harmonics: tuple[float, ...] = (0.5, 2.0, 3.0),
    period_bin_edges: tuple[float, ...] = (),
    max_peaks_per_period_bin: int = 0,
) -> list[dict[str, float]]:
    """Iteratively pick the top-N SDE peaks, masking each peak + harmonics.

    `extra` is a dict of period-indexed arrays (e.g. depth, duration, t0) that
    are sampled at the peak index for downstream record assembly.

    If `period_bin_edges` and `max_peaks_per_period_bin` are supplied, the walk
    still ranks by global SDE, but masks a period bin once its quota is used.
    This preserves strong peaks while preventing one crowded period range from
    consuming the entire top-N list.
    """
    if len(period) == 0:
        return []
    available = np.isfinite(sde)
    period_bins = np.full(len(period), -1, dtype=int)
    bin_counts: dict[int, int] = {}
    use_period_quota = bool(period_bin_edges) and max_peaks_per_period_bin > 0
    if use_period_quota:
        edges = np.asarray(period_bin_edges, dtype=float)
        edges = edges[np.isfinite(edges)]
        if edges.size >= 2:
            edges = np.unique(np.sort(edges))
            period_bins = np.searchsorted(edges, period, side="right") - 1
            inside = (period_bins >= 0) & (period < edges[-1])
            available &= inside
        else:
            use_period_quota = False
    peaks: list[dict[str, float]] = []
    sde_work = np.where(available, sde, -np.inf)
    for rank in range(1, n_peaks + 1):
        idx = int(np.argmax(sde_work))
        if not np.isfinite(sde_work[idx]):
            break
        p = float(period[idx])
        bin_idx = int(period_bins[idx]) if use_period_quota else -1
        rec: dict[str, float] = {
            "peak_rank": rank,
            "period_d": p,
            "sde": float(sde_work[idx]),
        }
        for k, arr in extra.items():
            rec[k] = float(arr[idx])
        peaks.append(rec)
        if use_period_quota:
            bin_counts[bin_idx] = bin_counts.get(bin_idx, 0) + 1
            if bin_counts[bin_idx] >= max_peaks_per_period_bin:
                sde_work[period_bins == bin_idx] = -np.inf
        # Mask the peak and harmonics in fractional period space.
        for mult in (1.0, *harmonics):
            target = p * mult
            if target <= 0:
                continue
            lo = target * (1.0 - period_mask_frac)
            hi = target * (1.0 + period_mask_frac)
            sde_work[(period >= lo) & (period <= hi)] = -np.inf
    return peaks


def _json_or_empty(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    try:
        return json.dumps(value)
    except TypeError:
        return str(value)


def _config_int(config: dict[str, Any], key: str, default: int = 0) -> int:
    try:
        return int(config.get(key, default))
    except (TypeError, ValueError):
        return int(default)


def _config_float(config: dict[str, Any], key: str, default: float = np.nan) -> float:
    try:
        value = float(config.get(key, default))
    except (TypeError, ValueError):
        return float(default)
    return value if np.isfinite(value) else float(default)


def result_to_rows(
    res: BLSResult,
    run_id: str,
    run_config: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Flatten a BLSResult into one or more candidate-table rows.

    Always emits at least one row. Status rows ("too_few_cadences", "read_fail",
    "all_nan") emit a single row with peak_rank=0 and NaN peak fields so the
    candidate table has full coverage of the input target list.
    """
    n_cad_quality = res.n_cad_kept if res.n_cad_quality is None else res.n_cad_quality
    quality_dropout_frac = (
        res.quality_dropout_frac
        if res.quality_dropout_frac is not None
        else 0.0 if res.n_cad_total == 0
        else (res.n_cad_total - n_cad_quality) / res.n_cad_total
    )
    cfg = dict(run_config or {})
    base = {
        "tic": res.tic,
        "sector": res.sector,
        "cam": res.cam,
        "ccd": res.ccd,
        "tmag": res.tmag,
        "aperture": res.aperture,
        "n_cad_total": res.n_cad_total,
        "n_cad_kept": res.n_cad_kept,
        "dropout_frac": res.dropout_frac,
        "n_cad_quality": n_cad_quality,
        "n_cad_edge_trimmed": res.n_cad_edge_trimmed,
        "n_cad_sigma_clipped": res.n_cad_sigma_clipped,
        "quality_dropout_frac": quality_dropout_frac,
        "n_orbits": res.n_orbits,
        "baseline_d": res.baseline_d,
        "status": res.status,
        "hlsp_path": res.hlsp_path,
        "bls_run_id": run_id,
        "bls_search_branch": str(cfg.get("search_branch", "")),
        "bls_n_periods": _config_int(cfg, "n_periods"),
        "bls_n_peaks": _config_int(cfg, "n_peaks"),
        "bls_p_min_d": _config_float(cfg, "p_min_d"),
        "bls_p_max_cap_d": _config_float(cfg, "p_max_cap_d"),
        "bls_max_period_fraction": _config_float(cfg, "max_period_fraction"),
        "bls_period_mask_frac": _config_float(cfg, "period_mask_frac"),
        "bls_period_bin_edges": _json_or_empty(cfg.get("period_bin_edges")),
        "bls_max_peaks_per_period_bin": _config_int(cfg, "max_peaks_per_period_bin"),
        "bls_sigma_clip": _config_float(cfg, "sigma_clip"),
        "bls_orbit_edge_trim_d": _config_float(cfg, "orbit_edge_trim_d"),
        "bls_durations_min": _json_or_empty(cfg.get("durations_min")),
    }
    if not res.peaks:
        empty = {
            "peak_rank": 0,
            "period_d": np.nan,
            "t0_bjd": np.nan,
            "duration_min": np.nan,
            "depth": np.nan,
            "depth_snr": np.nan,
            "sde": np.nan,
            "log_power": np.nan,
        }
        return [{**base, **empty}]
    rows = []
    for pk in res.peaks:
        rows.append({**base, **asdict(pk)})
    return rows


PARQUET_FIELDS: tuple[str, ...] = (
    "tic", "sector", "cam", "ccd", "tmag", "aperture",
    "peak_rank", "period_d", "t0_bjd", "duration_min", "depth",
    "depth_snr", "sde", "log_power",
    "n_cad_total", "n_cad_quality", "n_cad_kept",
    "n_cad_edge_trimmed", "n_cad_sigma_clipped",
    "dropout_frac", "quality_dropout_frac", "n_orbits", "baseline_d",
    "status", "hlsp_path", "bls_run_id", "bls_search_branch",
    "bls_n_periods", "bls_n_peaks", "bls_p_min_d", "bls_p_max_cap_d",
    "bls_max_period_fraction", "bls_period_mask_frac",
    "bls_period_bin_edges", "bls_max_peaks_per_period_bin",
    "bls_sigma_clip", "bls_orbit_edge_trim_d", "bls_durations_min",
)


def write_parquet(rows: list[dict[str, Any]], out_path: Path) -> None:
    """Write candidate rows to Parquet with a stable column order."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cols: dict[str, list[Any]] = {k: [] for k in PARQUET_FIELDS}
    for r in rows:
        for k in PARQUET_FIELDS:
            cols[k].append(r.get(k))
    table = pa.table(cols)
    pq.write_table(table, out_path, compression="zstd")
