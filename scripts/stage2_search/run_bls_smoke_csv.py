#!/usr/bin/env python3
"""Run a small BLS smoke test and write CSV output without pyarrow."""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import read_hlsp  # noqa: E402
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402


def _hlsp_path(root: Path, sector: int, tic: int) -> Path:
    s = f"{tic:016d}"
    return (
        root / s[:4] / s[4:8] / s[8:12] / s[12:]
        / f"hlsp_twirlfs_tess_ffi_s{sector:04d}-{s}_tess_v01_llc.fits"
    )


def _parse_tics(raw: str) -> list[int]:
    return [int(tok.strip()) for tok in raw.split(",") if tok.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlsp-root", type=Path, required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--tics", type=str, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--apertures", default="DET_FLUX_SML,DET_FLUX,DET_FLUX_LAG")
    parser.add_argument("--n-periods", type=int, default=50_000)
    parser.add_argument("--n-peaks", type=int, default=3)
    args = parser.parse_args(argv)

    tics = _parse_tics(args.tics)
    apertures = tuple(s.strip() for s in args.apertures.split(",") if s.strip())
    cfg = BLSConfig(apertures=apertures, n_periods=args.n_periods, n_peaks=args.n_peaks)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    t0 = time.time()
    for tic in tics:
        path = _hlsp_path(args.hlsp_root, args.sector, tic)
        lc = read_hlsp(path)
        if lc is None:
            rows.append({"tic": tic, "status": "read_fail", "hlsp_path": str(path)})
            print(f"READ_FAIL tic={tic} path={path}", flush=True)
            continue
        print(f"RUN tic={tic} T={lc.tmag:.3f} cad={len(lc.time)}", flush=True)
        for aperture in apertures:
            result = run_bls_on_lc(lc, cfg, aperture=aperture)
            if not result.peaks:
                rows.append(
                    {
                        "tic": tic,
                        "tmag": lc.tmag,
                        "sector": lc.sector,
                        "cam": lc.cam,
                        "ccd": lc.ccd,
                        "aperture": aperture,
                        "status": result.status,
                        "peak_rank": 0,
                        "n_cad_kept": result.n_cad_kept,
                        "dropout_frac": result.dropout_frac,
                        "hlsp_path": str(path),
                    }
                )
                print(f"  {aperture}: status={result.status} no peaks", flush=True)
                continue
            best = result.peaks[0]
            print(
                f"  {aperture}: P={best.period_d:.6f} d "
                f"dur={best.duration_min:.2f} min depth={best.depth:.4g} "
                f"snr={best.depth_snr:.2f} sde={best.sde:.2f}",
                flush=True,
            )
            for peak in result.peaks:
                rows.append(
                    {
                        "tic": tic,
                        "tmag": lc.tmag,
                        "sector": lc.sector,
                        "cam": lc.cam,
                        "ccd": lc.ccd,
                        "aperture": aperture,
                        "status": result.status,
                        "peak_rank": peak.peak_rank,
                        "period_d": peak.period_d,
                        "t0_bjd": peak.t0_bjd,
                        "duration_min": peak.duration_min,
                        "depth": peak.depth,
                        "depth_snr": peak.depth_snr,
                        "sde": peak.sde,
                        "log_power": peak.log_power,
                        "n_cad_kept": result.n_cad_kept,
                        "dropout_frac": result.dropout_frac,
                        "hlsp_path": str(path),
                    }
                )

    csv_path = args.out_dir / "bls_smoke_top_peaks.csv"
    fields = [
        "tic",
        "tmag",
        "sector",
        "cam",
        "ccd",
        "aperture",
        "status",
        "peak_rank",
        "period_d",
        "t0_bjd",
        "duration_min",
        "depth",
        "depth_snr",
        "sde",
        "log_power",
        "n_cad_kept",
        "dropout_frac",
        "hlsp_path",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    meta = {
        "hlsp_root": str(args.hlsp_root),
        "sector": args.sector,
        "tics": tics,
        "apertures": apertures,
        "n_periods": cfg.n_periods,
        "n_peaks": cfg.n_peaks,
        "wall_s": time.time() - t0,
        "csv": str(csv_path),
    }
    (args.out_dir / "meta.json").write_text(json.dumps(meta, indent=2))

    print(f"DONE rows={len(rows)} wall={meta['wall_s']:.1f}s csv={csv_path}", flush=True)
    valid = [r for r in rows if r.get("peak_rank") == 1 and r.get("status") == "ok"]
    valid.sort(key=lambda row: float(row.get("sde", float("-inf"))), reverse=True)
    print("TOP_BY_SDE", flush=True)
    for row in valid[:10]:
        print(
            f"tic={row['tic']} ap={row['aperture']} P={float(row['period_d']):.6f} "
            f"dur={float(row['duration_min']):.2f} depth={float(row['depth']):.4g} "
            f"snr={float(row['depth_snr']):.2f} sde={float(row['sde']):.2f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
