"""Peak walk + dataclasses + Parquet schema for the per-sector BLS first pass.

A "candidate" is a row in the per-sector candidate table, keyed by
(tic, sector, aperture, peak_rank). Peaks are extracted by walking the BLS
SDE spectrum: argmax → mask the peak ± a fractional period window plus its
1/2× and 2× harmonics → repeat N times.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
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
) -> list[dict[str, float]]:
    """Iteratively pick the top-N SDE peaks, masking each peak + harmonics.

    `extra` is a dict of period-indexed arrays (e.g. depth, duration, t0) that
    are sampled at the peak index for downstream record assembly.
    """
    if len(period) == 0:
        return []
    available = np.isfinite(sde)
    peaks: list[dict[str, float]] = []
    sde_work = np.where(available, sde, -np.inf)
    for rank in range(1, n_peaks + 1):
        idx = int(np.argmax(sde_work))
        if not np.isfinite(sde_work[idx]):
            break
        p = float(period[idx])
        rec: dict[str, float] = {
            "peak_rank": rank,
            "period_d": p,
            "sde": float(sde_work[idx]),
        }
        for k, arr in extra.items():
            rec[k] = float(arr[idx])
        peaks.append(rec)
        # Mask the peak and harmonics in fractional period space.
        for mult in (1.0, *harmonics):
            target = p * mult
            if target <= 0:
                continue
            lo = target * (1.0 - period_mask_frac)
            hi = target * (1.0 + period_mask_frac)
            sde_work[(period >= lo) & (period <= hi)] = -np.inf
    return peaks


def result_to_rows(res: BLSResult, run_id: str) -> list[dict[str, Any]]:
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
    "status", "hlsp_path", "bls_run_id",
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
