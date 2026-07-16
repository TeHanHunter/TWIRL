#!/usr/bin/env python3
"""Build a full-S56 real-candidate BLS peak table from the ADP pair only."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
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

from twirl.io.compact_export import read_compact_lc_export  # noqa: E402
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from twirl.vetting.adp_only import (  # noqa: E402
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
    validate_adp_only_apertures,
)
from twirl.vetting.recovery50_teacher import json_default, write_table  # noqa: E402


DEFAULT_COMPACT_LC = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/"
    "s56_twirlfs_v2_adp_lc_export_pdo.h5"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_adp_real_bls_peaks"


def _sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _target_tics(path: Path) -> list[int]:
    with h5py.File(path, "r") as h5:
        if "targets" not in h5:
            raise KeyError(f"compact export has no /targets group: {path}")
        return sorted(int(key) for key in h5["targets"].keys())


def _result_rows(result: Any) -> list[dict[str, Any]]:
    base = {
        "tic": int(result.tic),
        "sector": int(result.sector),
        "cam": int(result.cam),
        "ccd": int(result.ccd),
        "tmag": float(result.tmag),
        "aperture": str(result.aperture),
        "n_cad_total": int(result.n_cad_total),
        "n_cad_quality": int(result.n_cad_quality or 0),
        "n_cad_kept": int(result.n_cad_kept),
        "n_cad_edge_trimmed": int(result.n_cad_edge_trimmed),
        "n_cad_sigma_clipped": int(result.n_cad_sigma_clipped),
        "dropout_frac": float(result.dropout_frac),
        "quality_dropout_frac": float(result.quality_dropout_frac or 0.0),
        "n_orbits": int(result.n_orbits),
        "baseline_d": float(result.baseline_d),
        "status": str(result.status),
        "bls_search_branch": "current_adp",
        "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
    }
    if not result.peaks:
        return [
            {
                **base,
                "peak_rank": 0,
                "period_d": np.nan,
                "t0_bjd": np.nan,
                "duration_min": np.nan,
                "depth": np.nan,
                "depth_snr": np.nan,
                "sde": np.nan,
                "log_power": np.nan,
            }
        ]
    return [{**base, **asdict(peak)} for peak in result.peaks]


def _process_target(payload: tuple[int, str, dict[str, Any]]) -> list[dict[str, Any]]:
    tic, compact_lc_s, cfg_payload = payload
    compact_lc = Path(compact_lc_s)
    cfg = BLSConfig(
        apertures=ADP_ONLY_APERTURES,
        n_periods=int(cfg_payload["n_periods"]),
        n_peaks=int(cfg_payload["n_peaks"]),
        p_min_d=float(cfg_payload["p_min_d"]),
        p_max_cap_d=float(cfg_payload["p_max_cap_d"]),
        max_period_fraction=float(cfg_payload["max_period_fraction"]),
        sigma_clip=float(cfg_payload["sigma_clip"]),
        orbit_edge_trim_d=float(cfg_payload["orbit_edge_trim_d"]),
    )
    config_provenance = {
        "bls_n_periods": int(cfg.n_periods),
        "bls_n_peaks": int(cfg.n_peaks),
        "bls_p_min_d": float(cfg.p_min_d),
        "bls_p_max_cap_d": float(cfg.p_max_cap_d),
        "bls_max_period_fraction": float(cfg.max_period_fraction),
        "bls_sigma_clip": float(cfg.sigma_clip),
        "bls_orbit_edge_trim_d": float(cfg.orbit_edge_trim_d),
    }
    lc = read_compact_lc_export(compact_lc, tic=tic, columns=ADP_ONLY_APERTURES)
    if lc is None:
        return [
            {
                "tic": tic,
                "status": "missing_adp_pair",
                "bls_search_branch": "current_adp",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
                "source_product_tag": str(cfg_payload.get("source_product_tag", "")),
                **config_provenance,
            }
        ]
    rows: list[dict[str, Any]] = []
    for aperture in ADP_ONLY_APERTURES:
        result = run_bls_on_lc(lc, cfg, aperture=aperture)
        current = _result_rows(result)
        for row in current:
            row["source_product_tag"] = str(cfg_payload.get("source_product_tag", ""))
            row.update(config_provenance)
        rows.extend(current)
    return rows


def build_peak_table(
    *,
    compact_lc: Path,
    out_dir: Path,
    workers: int,
    n_periods: int,
    n_peaks: int,
    max_targets: int | None,
    progress_every: int,
    shard_index: int = 0,
    n_shards: int = 1,
    resume: bool = False,
    source_product_tag: str = "",
) -> dict[str, Any]:
    validate_adp_only_apertures(ADP_ONLY_APERTURES)
    out_dir.mkdir(parents=True, exist_ok=True)
    tics = _target_tics(compact_lc)
    if max_targets is not None:
        tics = tics[: max(0, int(max_targets))]
    if n_shards < 1 or shard_index < 0 or shard_index >= n_shards:
        raise ValueError("shard_index must satisfy 0 <= shard_index < n_shards")
    n_targets_total = len(tics)
    tics = tics[int(shard_index) :: int(n_shards)]
    suffix = f"_{int(shard_index):03d}" if int(n_shards) > 1 else ""
    output_path = out_dir / f"real_adp_bls_peaks{suffix}.parquet"
    summary_path = out_dir / f"summary{suffix}.json"
    cfg_payload = {
        "n_periods": int(n_periods),
        "n_peaks": int(n_peaks),
        "p_min_d": 0.12,
        "p_max_cap_d": 15.0,
        "max_period_fraction": 0.45,
        "sigma_clip": 5.0,
        "orbit_edge_trim_d": 0.0,
        "source_product_tag": str(source_product_tag),
    }
    if resume and output_path.exists() and summary_path.exists():
        summary = json.loads(summary_path.read_text())
        if (
            summary.get("contract_version") == ADP_ONLY_CONTRACT_VERSION
            and summary.get("compact_lc") == str(compact_lc)
            and int(summary.get("shard_index", -1)) == int(shard_index)
            and int(summary.get("n_shards", -1)) == int(n_shards)
            and int(summary.get("n_targets", -1)) == len(tics)
            and int(summary.get("n_targets_total", -1)) == n_targets_total
            and summary.get("source_product_tag") == str(source_product_tag)
            and summary.get("config") == cfg_payload
            and summary.get("peak_table_sha256") == _sha256(output_path)
        ):
            print(json.dumps(summary, indent=2, sort_keys=True))
            return summary
    payloads = [(tic, str(compact_lc), cfg_payload) for tic in tics]
    rows: list[dict[str, Any]] = []
    workers = max(1, int(workers))
    if workers == 1:
        iterator = map(_process_target, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        iterator = executor.map(_process_target, payloads, chunksize=1)
    try:
        for index, batch in enumerate(iterator, start=1):
            rows.extend(batch)
            if progress_every > 0 and index % int(progress_every) == 0:
                print(f"[adp-real-bls] processed {index:,}/{len(payloads):,}", flush=True)
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    peaks = pd.DataFrame(rows)
    output_path = write_table(peaks, output_path)
    output_sha256 = _sha256(output_path)
    status = peaks.get("status", pd.Series(dtype=str)).fillna("").astype(str)
    valid = status.eq("ok") & pd.to_numeric(peaks.get("peak_rank"), errors="coerce").gt(0)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": ADP_ONLY_CONTRACT_VERSION,
        "compact_lc": str(compact_lc),
        "out_dir": str(out_dir),
        "apertures": list(ADP_ONLY_APERTURES),
        "n_targets": int(len(tics)),
        "n_targets_total": int(n_targets_total),
        "n_rows": int(len(peaks)),
        "n_valid_peak_rows": int(valid.sum()),
        "n_periods": int(n_periods),
        "n_peaks": int(n_peaks),
        "workers": int(workers),
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
        "source_product_tag": str(source_product_tag),
        "peak_table_sha256": output_sha256,
        "config": cfg_payload,
        "status_counts": {str(k): int(v) for k, v in status.value_counts().sort_index().items()},
        "aperture_counts": {
            str(k): int(v)
            for k, v in peaks.get("aperture", pd.Series(dtype=str)).fillna("").astype(str).value_counts().sort_index().items()
        },
        "outputs": {
            "peak_table": str(output_path),
            "summary": str(summary_path),
        },
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact-lc", type=Path, default=DEFAULT_COMPACT_LC)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--n-periods", type=int, default=50_000)
    parser.add_argument("--n-peaks", type=int, default=10)
    parser.add_argument("--max-targets", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--source-product-tag",
        default="",
        help="Set to A2v1 when the compact export came from the production A2v1 FITS tree.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    build_peak_table(
        compact_lc=args.compact_lc,
        out_dir=args.out_dir,
        workers=args.workers,
        n_periods=args.n_periods,
        n_peaks=args.n_peaks,
        max_targets=args.max_targets,
        progress_every=args.progress_every,
        shard_index=args.shard_index,
        n_shards=args.n_shards,
        resume=args.resume,
        source_product_tag=args.source_product_tag,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
