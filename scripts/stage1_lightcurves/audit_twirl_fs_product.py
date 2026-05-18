#!/usr/bin/env python3
"""Audit selected TICs in a TWIRL-FS HLSP tree."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits


def _hlsp_path(root: Path, sector: int, tic: int) -> Path:
    s = f"{tic:016d}"
    return (
        root / s[0:4] / s[4:8] / s[8:12] / s[12:16]
        / f"hlsp_twirlfs_tess_ffi_s{sector:04d}-{s}_tess_v01_llc.fits"
    )


def _mad_sigma(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    med = np.nanmedian(x)
    return float(1.4826 * np.nanmedian(np.abs(x - med)))


def audit_one(path: Path, tic: int) -> dict[str, object]:
    row: dict[str, object] = {"tic": int(tic), "path": str(path), "exists": path.exists()}
    if not path.exists():
        return row
    with fits.open(path, memmap=False) as hdul:
        header = hdul[0].header
        data = hdul[1].data
        q0 = data["QUALITY"] == 0
        sap = np.asarray(data["SAP_FLUX"], dtype=float)
        det = np.asarray(data["DET_FLUX"], dtype=float)
        neg_q0 = q0 & np.isfinite(sap) & (sap < 0)
        det_q0 = det[q0]
        row.update(
            pipeline=str(header.get("PIPELINE", "")),
            method=str(header.get("METHOD", "")),
            bkspace=float(header.get("BKSPACE", np.nan)),
            sigclip=float(header.get("SIGCLIP", np.nan)),
            outmode=str(header.get("OUTMODE", "")),
            scale=str(header.get("SCALE", "")),
            n_cadences=int(len(data)),
            n_q0=int(q0.sum()),
            finite_sap_q0=int(np.isfinite(sap[q0]).sum()),
            finite_det_q0=int(np.isfinite(det_q0).sum()),
            neg_sap_q0=int(neg_q0.sum()),
            finite_det_for_neg_sap_q0=int(np.isfinite(det[neg_q0]).sum()),
            det_median_q0=float(np.nanmedian(det_q0)),
            det_mad_q0=_mad_sigma(det_q0),
        )
    return row


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--tic", type=int, action="append", required=True)
    args = parser.parse_args(argv)

    rows = [audit_one(_hlsp_path(args.root, args.sector, tic), tic) for tic in args.tic]
    summary = {
        "root": str(args.root),
        "sector": int(args.sector),
        "n_targets": len(rows),
        "n_existing": sum(1 for row in rows if row["exists"]),
        "pipelines": sorted({row.get("pipeline", "") for row in rows if row["exists"]}),
        "methods": sorted({row.get("method", "") for row in rows if row["exists"]}),
        "bkspace_values": sorted({row.get("bkspace", None) for row in rows if row["exists"]}),
        "sigclip_values": sorted({row.get("sigclip", None) for row in rows if row["exists"]}),
        "total_q0": int(sum(row.get("n_q0", 0) for row in rows)),
        "total_finite_det_q0": int(sum(row.get("finite_det_q0", 0) for row in rows)),
        "total_neg_sap_q0": int(sum(row.get("neg_sap_q0", 0) for row in rows)),
        "total_finite_det_for_neg_sap_q0": int(
            sum(row.get("finite_det_for_neg_sap_q0", 0) for row in rows)
        ),
        "targets": rows,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
