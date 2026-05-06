"""Per-sector BLS orchestrator for the TWIRL Stage 2 first pass.

Inputs:
    --hlsp-root  /pdo/users/tehan/tglc-gpu-production/hlsp_s00NN/   (or override)
    --sector     56

Output:
    {out_dir}/sector_{NN:04d}/candidates.parquet
    {out_dir}/sector_{NN:04d}/meta.json
    {out_dir}/sector_{NN:04d}/periodograms/tic_{TIC}.npz   (only for --save-periodograms-for)

Designed to be invoked as:
    OMP_NUM_THREADS=1 python -m twirl.search.sector_run --sector 56 --workers 32
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import platform
import socket
import subprocess
import time
import uuid
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from twirl.io.hlsp import discover_sector_targets, read_hlsp
from twirl.search.bls import BLSConfig, run_bls_on_lc, save_periodogram
from twirl.search.candidates import result_to_rows, write_parquet

DEFAULT_HLSP_ROOT_FMT = "/pdo/users/tehan/tglc-gpu-production/hlsp_s{sector:04d}"
WD_1856_TIC = 267574918
DEFAULT_OUT_DIR = "data_local/stage2/bls_first_pass"

# Module-level config + state for pool workers (set by _init_worker via initializer).
_WORKER_CFG: BLSConfig | None = None
_WORKER_PERIODOGRAM_TICS: set[int] = set()
_WORKER_PERIODOGRAM_DIR: Path | None = None
_WORKER_RUN_ID: str = ""


def _init_worker(cfg: BLSConfig, periodogram_tics: set[int],
                 periodogram_dir: Path | None, run_id: str) -> None:
    global _WORKER_CFG, _WORKER_PERIODOGRAM_TICS, _WORKER_PERIODOGRAM_DIR, _WORKER_RUN_ID
    _WORKER_CFG = cfg
    _WORKER_PERIODOGRAM_TICS = set(periodogram_tics)
    _WORKER_PERIODOGRAM_DIR = periodogram_dir
    _WORKER_RUN_ID = run_id


def _process_one(path: Path) -> list[dict[str, Any]]:
    """Read + run BLS on one HLSP path. Returns candidate-table rows."""
    cfg = _WORKER_CFG
    assert cfg is not None
    lc = read_hlsp(path)
    if lc is None:
        # Fabricate a minimal status row so the input file is accounted for.
        from twirl.search.candidates import BLSResult
        res = BLSResult(
            tic=-1, sector=-1, cam=-1, ccd=-1, tmag=float("nan"),
            aperture="DET_FLUX", n_cad_total=0, n_cad_kept=0,
            dropout_frac=0.0, n_orbits=0, baseline_d=0.0,
            status="read_fail", hlsp_path=str(path), peaks=[],
        )
        return result_to_rows(res, _WORKER_RUN_ID)

    rows: list[dict[str, Any]] = []
    save_pg = lc.tic in _WORKER_PERIODOGRAM_TICS and _WORKER_PERIODOGRAM_DIR is not None
    for ap in cfg.apertures:
        out = run_bls_on_lc(lc, cfg, aperture=ap, return_periodogram=save_pg)
        if save_pg:
            res, spec = out
            pg_path = _WORKER_PERIODOGRAM_DIR / f"tic_{lc.tic}_{ap}.npz"
            try:
                save_periodogram(spec, pg_path)
            except Exception:
                pass
        else:
            res = out
        rows.extend(result_to_rows(res, _WORKER_RUN_ID))
    return rows


def _git_short_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True
        ).strip()
    except Exception:
        return "unknown"


def resolve_hlsp_root(sector: int, override: Path | None) -> Path:
    if override is not None:
        return Path(override)
    return Path(DEFAULT_HLSP_ROOT_FMT.format(sector=sector))


def run_sector(
    sector: int,
    hlsp_root: Path,
    out_dir: Path,
    cfg: BLSConfig,
    workers: int,
    limit: int | None,
    periodogram_tics: list[int],
    out_suffix: str = "",
) -> Path:
    """Execute BLS over every HLSP target in `hlsp_root` for `sector`."""
    hlsp_root = Path(hlsp_root)
    if not hlsp_root.exists():
        raise FileNotFoundError(f"HLSP root not found: {hlsp_root}")

    paths = discover_sector_targets(hlsp_root, sector)
    if limit is not None:
        paths = paths[:limit]
    if not paths:
        raise RuntimeError(f"No HLSP files for sector {sector} under {hlsp_root}")

    sector_out = Path(out_dir) / f"sector_{sector:04d}{out_suffix}"
    sector_out.mkdir(parents=True, exist_ok=True)
    pg_dir = sector_out / "periodograms"
    pg_tics = set(int(t) for t in periodogram_tics)
    if pg_tics:
        pg_dir.mkdir(parents=True, exist_ok=True)

    run_id = f"{_git_short_sha()}-{uuid.uuid4().hex[:8]}"
    t_start = time.time()
    rows: list[dict[str, Any]] = []
    n_attempted = len(paths)
    n_ok = 0
    n_failed = 0

    print(f"[stage2-bls] sector={sector} hlsp_root={hlsp_root} n_targets={n_attempted}"
          f" workers={workers} run_id={run_id}", flush=True)

    if workers <= 1:
        _init_worker(cfg, pg_tics, pg_dir if pg_tics else None, run_id)
        for i, p in enumerate(paths):
            try:
                target_rows = _process_one(p)
                rows.extend(target_rows)
                if any(r.get("status") == "ok" for r in target_rows):
                    n_ok += 1
                else:
                    n_failed += 1
            except Exception:
                n_failed += 1
            if (i + 1) % 200 == 0:
                print(f"  [stage2-bls] {i+1}/{n_attempted} processed "
                      f"(ok={n_ok} fail={n_failed})", flush=True)
    else:
        with mp.Pool(
            processes=workers,
            initializer=_init_worker,
            initargs=(cfg, pg_tics, pg_dir if pg_tics else None, run_id),
        ) as pool:
            for i, target_rows in enumerate(
                pool.imap_unordered(_process_one, paths, chunksize=8)
            ):
                rows.extend(target_rows)
                if any(r.get("status") == "ok" for r in target_rows):
                    n_ok += 1
                else:
                    n_failed += 1
                if (i + 1) % 500 == 0:
                    print(f"  [stage2-bls] {i+1}/{n_attempted} processed "
                          f"(ok={n_ok} fail={n_failed})", flush=True)

    candidates_path = sector_out / "candidates.parquet"
    write_parquet(rows, candidates_path)

    wall_time = time.time() - t_start
    meta = {
        "sector": int(sector),
        "hlsp_root": str(hlsp_root),
        "out_dir": str(sector_out),
        "n_attempted": int(n_attempted),
        "n_ok": int(n_ok),
        "n_failed": int(n_failed),
        "n_rows": int(len(rows)),
        "workers": int(workers),
        "wall_time_s": float(wall_time),
        "run_id": run_id,
        "git_sha": _git_short_sha(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "config": asdict(cfg),
        "periodogram_tics": sorted(int(t) for t in pg_tics),
    }
    with open(sector_out / "meta.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=str)

    print(f"[stage2-bls] done sector={sector} ok={n_ok} fail={n_failed} "
          f"rows={len(rows)} wall={wall_time:.1f}s -> {candidates_path}",
          flush=True)
    return candidates_path


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--hlsp-root", type=Path, default=None,
                    help="Override default /pdo/users/tehan/tglc-gpu-production/hlsp_s00NN/.")
    ap.add_argument("--out-dir", type=Path, default=Path(DEFAULT_OUT_DIR))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) // 2))
    ap.add_argument("--limit", type=int, default=None,
                    help="Process only the first N HLSP files (smoke testing).")
    ap.add_argument("--apertures", type=str,
                    default="DET_FLUX_SML,DET_FLUX,DET_FLUX_LAG",
                    help="Comma-separated list of aperture flux columns to search. "
                         "All three by default — different CCDs/fields favor different "
                         "apertures (see WD 1856 S59 recovery: SML works, medium fails).")
    ap.add_argument("--n-periods", type=int, default=100_000)
    ap.add_argument("--n-peaks", type=int, default=10)
    ap.add_argument("--orbit-edge-trim-d", type=float, default=0.0,
                    help="Drop cadences within this many days of each orbit's "
                         "first/last cadence (default 0 = no trim).")
    ap.add_argument("--out-suffix", type=str, default="",
                    help="Append to sector_NNNN dir name (e.g. '_trim1d') so a "
                         "two-pass run does not clobber the first-pass output.")
    ap.add_argument("--save-periodograms-for", type=str, default="",
                    help="Comma-separated TIC IDs whose full periodograms to save. "
                         "WD 1856 (267574918) is added automatically.")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)

    apertures = tuple(s.strip() for s in args.apertures.split(",") if s.strip())
    pg_tics = [int(s) for s in args.save_periodograms_for.split(",") if s.strip()]
    if WD_1856_TIC not in pg_tics:
        pg_tics.append(WD_1856_TIC)  # always preserve WD 1856 spectrum if present

    cfg = BLSConfig(
        apertures=apertures,
        n_periods=int(args.n_periods),
        n_peaks=int(args.n_peaks),
        orbit_edge_trim_d=float(args.orbit_edge_trim_d),
    )
    hlsp_root = resolve_hlsp_root(args.sector, args.hlsp_root)

    run_sector(
        sector=int(args.sector),
        hlsp_root=hlsp_root,
        out_dir=Path(args.out_dir),
        cfg=cfg,
        workers=int(args.workers),
        limit=args.limit,
        periodogram_tics=pg_tics,
        out_suffix=args.out_suffix,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
