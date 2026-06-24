#!/usr/bin/env python3
"""Export compact S56 TWIRL-FS light curves for ORCD/H200 training runs.

This is the transfer boundary for the first S56 injection and ML pilot: keep
raw TGLC/TICA products on PDO, but move compact downstream HLSP-derived arrays
plus a manifest to ORCD.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import APERTURES, BJDREFI, read_hlsp

DEFAULT_HLSP_ROOT = (
    REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export.h5"
)


def _read_tic_filter(path: Path | None, tic_column: str) -> set[int] | None:
    if path is None:
        return None
    import pandas as pd

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        table = pd.read_parquet(path, columns=[tic_column])
    elif suffix == ".csv":
        table = pd.read_csv(path, usecols=[tic_column])
    elif suffix in {".json", ".jsonl"}:
        table = pd.read_json(path, lines=suffix == ".jsonl")
    else:
        raise ValueError(f"unsupported TIC table format: {path}")
    if tic_column not in table.columns:
        raise KeyError(f"TIC column {tic_column!r} not found in {path}")
    return {int(tic) for tic in table[tic_column].dropna().astype(int).tolist()}


def _dataset(group, name: str, values, compression: str | None) -> None:
    kwargs = {"shuffle": True}
    if compression:
        kwargs["compression"] = compression
    group.create_dataset(name, data=values, **kwargs)


def _iter_sector_targets(hlsp_root: Path, sector: int, limit: int | None):
    pat = f"hlsp_*_tess_ffi_s{sector:04d}-*.fits"
    count = 0
    for path in Path(hlsp_root).rglob(pat):
        yield path
        count += 1
        if limit is not None and count >= limit:
            break


def _write_export(
    *,
    hlsp_root: Path,
    out_h5: Path,
    sector: int,
    columns: tuple[str, ...],
    tic_filter: set[int] | None,
    limit: int | None,
    compression: str | None,
    progress_every: int,
    overwrite: bool,
) -> dict:
    import h5py

    if tic_filter is not None:
        # Keep path discovery cheap by filtering after FITS header read; file
        # names encode TIC, but the header is authoritative for the export.
        target_count = len(tic_filter)
    else:
        target_count = None
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    if out_h5.exists() and not overwrite:
        raise FileExistsError(f"output exists; pass --overwrite: {out_h5}")

    tmp_path = out_h5.with_suffix(out_h5.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    manifest_records = []
    skipped = {"read_failed": 0, "tic_filter": 0, "duplicate_tic": 0, "no_flux_columns": 0}
    discovered = 0
    exported = 0
    created_utc = datetime.now(timezone.utc).isoformat()

    with h5py.File(tmp_path, "w") as h5:
        h5.attrs["created_utc"] = created_utc
        h5.attrs["sector"] = int(sector)
        h5.attrs["hlsp_root"] = str(hlsp_root)
        h5.attrs["time_column"] = "TIME"
        h5.attrs["time_unit"] = f"BJD - {BJDREFI}"
        h5.attrs["flux_columns"] = json.dumps(list(columns))
        h5.attrs["target_filter_size"] = -1 if target_count is None else int(target_count)
        targets = h5.create_group("targets")

        for idx, path in enumerate(_iter_sector_targets(hlsp_root, sector, limit), start=1):
            discovered += 1
            lc = read_hlsp(path, columns=columns)
            if lc is None:
                skipped["read_failed"] += 1
                continue
            if tic_filter is not None and lc.tic not in tic_filter:
                skipped["tic_filter"] += 1
                continue
            key = f"{lc.tic:016d}"
            if key in targets:
                skipped["duplicate_tic"] += 1
                continue
            flux_cols = [col for col in columns if col in lc.flux]
            if not flux_cols:
                skipped["no_flux_columns"] += 1
                continue

            group = targets.create_group(key)
            group.attrs["tic"] = int(lc.tic)
            group.attrs["sector"] = int(lc.sector)
            group.attrs["camera"] = int(lc.cam)
            group.attrs["ccd"] = int(lc.ccd)
            group.attrs["tessmag"] = float(lc.tmag)
            group.attrs["ra_obj"] = float(getattr(lc, "ra", np.nan))
            group.attrs["dec_obj"] = float(getattr(lc, "dec", np.nan))
            group.attrs["source_fits"] = str(path)

            _dataset(group, "time", np.asarray(lc.time, dtype=np.float64), compression)
            _dataset(group, "cadenceno", np.asarray(lc.cadenceno, dtype=np.int32), compression)
            _dataset(group, "quality", np.asarray(lc.quality, dtype=np.int32), compression)
            _dataset(group, "orbitid", np.asarray(lc.orbitid, dtype=np.int16), compression)
            for col in flux_cols:
                _dataset(group, col, np.asarray(lc.flux[col], dtype=np.float32), compression)

            exported += 1
            manifest_records.append(
                {
                    "tic": int(lc.tic),
                    "sector": int(lc.sector),
                    "camera": int(lc.cam),
                    "ccd": int(lc.ccd),
                    "tessmag": float(lc.tmag),
                    "n_cadences": int(len(lc.time)),
                    "flux_columns": flux_cols,
                    "source_fits": str(path),
                }
            )
            if progress_every > 0 and exported % progress_every == 0:
                print(f"[export-lc] exported {exported:,} targets", flush=True)

        h5.attrs["n_targets"] = int(exported)

    if discovered == 0:
        tmp_path.unlink(missing_ok=True)
        raise FileNotFoundError(f"no S{sector:04d} HLSP FITS found under {hlsp_root}")

    os.replace(tmp_path, out_h5)
    manifest = {
        "created_utc": created_utc,
        "sector": int(sector),
        "hlsp_root": str(hlsp_root),
        "out_h5": str(out_h5),
        "time_unit": f"BJD - {BJDREFI}",
        "requested_columns": list(columns),
        "n_discovered_files": int(discovered),
        "n_exported_targets": int(exported),
        "skipped": skipped,
        "records": manifest_records,
    }
    manifest_path = out_h5.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    ap.add_argument("--out-h5", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--sector", type=int, default=56)
    ap.add_argument("--candidate-table", type=Path, default=None,
                    help="Optional CSV/Parquet/JSON table used to restrict exported TICs.")
    ap.add_argument("--tic-column", default="tic")
    ap.add_argument("--columns", nargs="+", default=list(APERTURES),
                    help="Flux columns to export from the LIGHTCURVE extension.")
    ap.add_argument("--limit", type=int, default=None,
                    help="Optional file-count limit for smoke tests.")
    ap.add_argument("--compression", choices=("lzf", "gzip", "none"), default="lzf")
    ap.add_argument("--progress-every", type=int, default=1000)
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    compression = None if args.compression == "none" else args.compression
    tic_filter = _read_tic_filter(args.candidate_table, args.tic_column)
    manifest = _write_export(
        hlsp_root=args.hlsp_root,
        out_h5=args.out_h5,
        sector=args.sector,
        columns=tuple(args.columns),
        tic_filter=tic_filter,
        limit=args.limit,
        compression=compression,
        progress_every=args.progress_every,
        overwrite=args.overwrite,
    )
    print("[export-lc] complete")
    print(f"  exported targets: {manifest['n_exported_targets']:,}")
    print(f"  skipped: {manifest['skipped']}")
    print(f"  out: {manifest['out_h5']}")
    print(f"  manifest: {Path(manifest['out_h5']).with_suffix('.manifest.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
