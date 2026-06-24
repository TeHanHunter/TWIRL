#!/usr/bin/env python3
"""Run the TGLC centroid on-target test on every TIC in a vetted parquet.

Appends columns ``centroid_dx_pix``, ``centroid_dy_pix``,
``centroid_delta_pix``, ``centroid_sigma_oot_pix``, ``centroid_z``,
``centroid_status``, ``centroid_pass`` to the input table and writes
``vetted_per_tic_centroid.parquet`` next to it.

The test is the TWIRL-canonical pixel vetting (see
``src/twirl/vetting/centroid_offset.py``): an on-target signal leaves
the flux-weighted centroid unchanged through transit; a background-EB
contaminator pulls the centroid by an amount that scales with the
contaminator's distance and apparent depth. Runs in milliseconds per
TIC; multiprocessing-parallel.

Paths via env vars so the same script runs locally and on PDO:
    TWIRL_HLSP_ROOT     -- directory holding hlsp_*.fits files
    TWIRL_INPUT_PARQUET -- vetted_per_tic.parquet to enrich
    TWIRL_OUTPUT_PARQUET -- output path (default: <input>_centroid.parquet)
    TWIRL_WORKERS       -- multiprocessing pool size (default 1)
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import read_hlsp  # noqa: E402
from twirl.vetting.centroid_offset import centroid_in_out_shift  # noqa: E402
from twirl.vetting.lightcurve_label_app import find_hlsp_path  # noqa: E402

WD_1856_TIC = 267574918


def _hlsp_path(root: Path, tic: int, sector: int) -> Path | None:
    s = f"{tic:016d}"
    legacy = (
        root / s[0:4] / s[4:8] / s[8:12] / s[12:16]
        / f"hlsp_qlp_tess_ffi_s{sector:04d}-{s}_tess_v01_llc.fits"
    )
    if legacy.exists():
        return legacy
    return find_hlsp_path(root, tic, sector)


def _check_one(args: tuple) -> dict:
    tic, sector, P, T0, dur, hlsp_root_str = args
    hlsp_root = Path(hlsp_root_str)
    p = _hlsp_path(hlsp_root, tic, sector)
    if p is None or not p.exists():
        return {"tic": int(tic), "centroid_status": "no_hlsp", "centroid_pass": False}
    try:
        lc = read_hlsp(p, columns=("SAP_FLUX",))
        if lc is None:
            return {"tic": int(tic), "centroid_status": "read_fail", "centroid_pass": False}
        # SAP_X / SAP_Y aren't in our read_hlsp HLSPLightCurve dataclass; read them
        # directly from the FITS via astropy.
        from astropy.io import fits
        with fits.open(p) as h:
            d = h[1].data
            sap_x = np.asarray(d["SAP_X"], dtype=float) if "SAP_X" in d.columns.names else None
            sap_y = np.asarray(d["SAP_Y"], dtype=float) if "SAP_Y" in d.columns.names else None
            time_arr = np.asarray(d["TIME"], dtype=float) + 2457000.0
            quality = np.asarray(d["QUALITY"], dtype=int)
        if sap_x is None or sap_y is None:
            return {"tic": int(tic), "centroid_status": "no_sap_xy", "centroid_pass": False}
        r = centroid_in_out_shift(
            time_bjd=time_arr, sap_x=sap_x, sap_y=sap_y, quality=quality,
            period_d=float(P), t0_bjd=float(T0), duration_min=float(dur),
        )
        return {
            "tic": int(tic),
            "centroid_dx_pix": r.delta_x_pix,
            "centroid_dy_pix": r.delta_y_pix,
            "centroid_delta_pix": r.delta_pix,
            "centroid_sigma_oot_pix": r.sigma_oot_pix,
            "centroid_z": r.z_score,
            "centroid_status": r.status,
            "centroid_pass": r.centroid_pass,
            "n_in_transit": r.n_in_transit,
            "n_oot_band": r.n_oot,
        }
    except Exception as e:
        return {
            "tic": int(tic),
            "centroid_status": f"error: {type(e).__name__}: {e}"[:80],
            "centroid_pass": False,
        }


def main() -> int:
    hlsp_root = Path(os.environ.get(
        "TWIRL_HLSP_ROOT",
        "/Users/tehan/PycharmProjects/TWIRL/data_local/stage2/leo_inputs/sector_0056/pdo/users/tehan/tglc-gpu-production/hlsp_s0056",
    ))
    in_parquet = Path(os.environ.get(
        "TWIRL_INPUT_PARQUET",
        "/Users/tehan/PycharmProjects/TWIRL/data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic.parquet",
    ))
    out_parquet = Path(os.environ.get(
        "TWIRL_OUTPUT_PARQUET",
        str(in_parquet.parent / (in_parquet.stem + "_centroid.parquet")),
    ))
    workers = int(os.environ.get("TWIRL_WORKERS", "1"))

    print(f"[centroid] reading {in_parquet}")
    df = pd.read_parquet(in_parquet)
    print(f"[centroid] {len(df):,} TICs to test; hlsp_root={hlsp_root}; workers={workers}")

    work = [
        (int(r["tic"]), int(r["sector"]), float(r["period_d"]),
         float(r["t0_bjd"]), float(r["duration_min"]), str(hlsp_root))
        for _, r in df.iterrows()
    ]

    results: list[dict] = []
    t0 = time.time()
    if workers <= 1:
        for i, w in enumerate(work, 1):
            results.append(_check_one(w))
            if i % 500 == 0:
                rate = i / (time.time() - t0)
                print(f"  [centroid] {i:,}/{len(work):,}  rate={rate:.1f}/s "
                      f"eta={(len(work) - i) / rate / 60:.1f}min")
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_check_one, w): w[0] for w in work}
            for i, fut in enumerate(as_completed(futs), 1):
                results.append(fut.result())
                if i % 500 == 0:
                    rate = i / (time.time() - t0)
                    print(f"  [centroid] {i:,}/{len(work):,}  rate={rate:.1f}/s "
                          f"eta={(len(work) - i) / rate / 60:.1f}min")

    res_df = pd.DataFrame(results).set_index("tic")
    merged = df.set_index("tic").join(res_df, how="left").reset_index()

    # Update vet_class for off-target signals (don't touch on-target classes)
    if "vet_class" in merged.columns:
        off = (merged["centroid_status"] == "off_target")
        n_off = int(off.sum())
        if n_off > 0:
            merged.loc[off, "vet_class"] = "off_target_suspect"
        print(f"[centroid] reclassified {n_off} candidates as off_target_suspect")

    merged.to_parquet(out_parquet, compression="zstd")
    out_csv = out_parquet.with_suffix(".csv")
    merged.to_csv(out_csv, index=False)
    print(f"[centroid] wrote {out_parquet}  ({len(merged):,} rows)")
    print(f"[centroid] wrote {out_csv}  ({len(merged):,} rows)")

    status_counts = merged["centroid_status"].value_counts()
    print(f"\n[centroid] status distribution:")
    for k, v in status_counts.items():
        print(f"  {k:<25} {v:>6,}")

    # WD 1856 detail
    wd = merged[merged["tic"] == WD_1856_TIC]
    if not wd.empty:
        r = wd.iloc[0]
        print(f"\n[wd1856] centroid result:")
        print(f"  status        : {r.get('centroid_status')}")
        print(f"  delta_pix     : {r.get('centroid_delta_pix'):.4f}" if pd.notna(r.get("centroid_delta_pix")) else "  delta_pix     : nan")
        print(f"  sigma_oot_pix : {r.get('centroid_sigma_oot_pix'):.4f}" if pd.notna(r.get("centroid_sigma_oot_pix")) else "  sigma_oot_pix : nan")
        print(f"  z             : {r.get('centroid_z'):.2f}" if pd.notna(r.get("centroid_z")) else "  z             : nan")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
