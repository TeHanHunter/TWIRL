#!/usr/bin/env python3
"""Generate faint TIC-backed TGLC light curves from existing source/ePSF files.

This is a TWIRL helper for production tests like S94 where the shared PDO tree
already has source cutouts and ePSFs, but the current light curves were emitted
with a bright target catalog.  The script writes only to a user-owned output
tree and treats the shared TGLC products as read-only inputs.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.tglc_faint_allstar import (  # noqa: E402
    PDO_USER_ROOT,
    enforce_pdo_user_output,
    normalize_tic_catalog,
    source_target_table,
)


_WORKER_TIC_CATALOG: Table | None = None
_WORKER_TMAG_MAX: float | None = None
_WORKER_MAX_TARGETS_PER_SOURCE: int | None = None


def parse_ccd(value: str | None) -> list[tuple[int, int]]:
    if value is None or value.lower() == "all":
        return [(camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)]
    if "," in value and "-" not in value and ":" not in value:
        camera_s, ccd_s = value.split(",", 1)
        camera = int(camera_s)
        ccd = int(ccd_s)
        if not (1 <= camera <= 4 and 1 <= ccd <= 4):
            raise ValueError(f"Invalid camera/CCD: {camera}-{ccd}")
        return [(camera, ccd)]
    out = []
    for item in value.split(","):
        if "-" in item:
            camera_s, ccd_s = item.split("-", 1)
        elif ":" in item:
            camera_s, ccd_s = item.split(":", 1)
        else:
            raise ValueError("CCD entries must look like '1-1' or use --ccd all")
        camera = int(camera_s)
        ccd = int(ccd_s)
        if not (1 <= camera <= 4 and 1 <= ccd <= 4):
            raise ValueError(f"Invalid camera/CCD: {camera}-{ccd}")
        out.append((camera, ccd))
    return out


def parse_stages(value: str) -> set[str]:
    stages = {stage.strip().lower() for stage in value.split(",") if stage.strip()}
    allowed = {"catalogs", "lightcurves"}
    unknown = stages - allowed
    if unknown:
        raise ValueError(f"Unknown stages {sorted(unknown)}; allowed stages are {sorted(allowed)}")
    return stages


def catalog_path(out_root: Path, orbit: int, camera: int, ccd: int) -> Path:
    return out_root / f"orbit-{orbit}" / "ffi" / "catalogs" / f"TIC_cam{camera}_ccd{ccd}.ecsv"


def lc_dir(out_root: Path, orbit: int, camera: int, ccd: int) -> Path:
    return out_root / f"orbit-{orbit}" / "ffi" / f"cam{camera}" / f"ccd{ccd}" / "LC"


def source_dir(input_root: Path, orbit: int, camera: int, ccd: int) -> Path:
    return input_root / f"orbit-{orbit}" / "ffi" / f"cam{camera}" / f"ccd{ccd}" / "source"


def epsf_dir(input_root: Path, orbit: int, camera: int, ccd: int) -> Path:
    return input_root / f"orbit-{orbit}" / "ffi" / f"cam{camera}" / f"ccd{ccd}" / "epsf"


def write_large_ecsv(table: Table, path: Path, *, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table[:0].write(path, format="ascii.ecsv", overwrite=overwrite)
    with path.open("a") as handle:
        table.write(handle, format="ascii.fast_no_header", delimiter=" ", strip_whitespace=False)


def ensure_tic_catalog(
    *,
    path: Path,
    orbit: int,
    camera: int,
    ccd: int,
    tmag_max: float,
    nprocs: int,
    replace: bool,
) -> Table:
    if path.is_file() and not replace:
        table = Table.read(path)
        print(f"[catalog] reuse {path} ({len(table):,} rows)", flush=True)
        return table

    from tglc.scripts.catalogs import get_tic_catalog_data

    print(
        f"[catalog] query orbit={orbit} cam={camera} ccd={ccd} "
        f"Tmag<={tmag_max} -> {path}",
        flush=True,
    )
    table = get_tic_catalog_data(
        orbit,
        camera,
        ccd,
        magnitude_cutoff=float(tmag_max) + 1e-6,
        mdwarf_magnitude_cutoff=float(tmag_max) + 1e-6,
        nprocs=nprocs,
    )
    table = normalize_tic_catalog(table, tmag_max=tmag_max)
    write_large_ecsv(table, path, overwrite=True)
    print(f"[catalog] wrote {path} ({len(table):,} rows)", flush=True)
    return table


def source_epsf_pairs(
    *,
    input_root: Path,
    orbit: int,
    camera: int,
    ccd: int,
    limit_sources: int | None,
) -> list[tuple[Path, Path]]:
    sdir = source_dir(input_root, orbit, camera, ccd)
    edir = epsf_dir(input_root, orbit, camera, ccd)
    pairs: list[tuple[Path, Path]] = []
    for source_file in sorted(sdir.glob("source_*.pkl")):
        epsf_file = edir / f"epsf{source_file.stem.removeprefix('source')}.npy"
        if epsf_file.is_file():
            pairs.append((source_file, epsf_file))
    if limit_sources is not None:
        pairs = pairs[: int(limit_sources)]
    return pairs


def init_worker(tic_catalog_path: str, tmag_max: float, max_targets_per_source: int | None) -> None:
    global _WORKER_TIC_CATALOG, _WORKER_TMAG_MAX, _WORKER_MAX_TARGETS_PER_SOURCE
    _WORKER_TIC_CATALOG = Table.read(tic_catalog_path)
    _WORKER_TMAG_MAX = tmag_max
    _WORKER_MAX_TARGETS_PER_SOURCE = max_targets_per_source


def install_hdf5(tmp_path: Path, final_path: Path, *, replace: bool) -> bool:
    if replace:
        tmp_path.replace(final_path)
        return True
    try:
        os.link(tmp_path, final_path)
        tmp_path.unlink()
        return True
    except FileExistsError:
        tmp_path.unlink(missing_ok=True)
        return False


def process_source_job(job: tuple[str, str, str, bool, int, int]) -> dict[str, Any]:
    source_path_s, epsf_path_s, lc_dir_s, replace, psf_size, oversample = job
    if _WORKER_TIC_CATALOG is None:
        raise RuntimeError("worker TIC catalog was not initialized")

    source_path = Path(source_path_s)
    epsf_path = Path(epsf_path_s)
    out_dir = Path(lc_dir_s)
    out_dir.mkdir(parents=True, exist_ok=True)

    from tglc.light_curve import generate_light_curves

    with source_path.open("rb") as handle:
        source = pickle.load(handle)

    matched = source_target_table(
        source.gaia,
        _WORKER_TIC_CATALOG,
        tmag_max=_WORKER_TMAG_MAX,
        max_targets=_WORKER_MAX_TARGETS_PER_SOURCE,
    )
    n_matched = len(matched)
    if n_matched == 0:
        return {
            "source": source_path.name,
            "matched": 0,
            "requested": 0,
            "written": 0,
            "skipped_existing": 0,
            "errors": [],
        }

    keep = []
    skipped_existing = 0
    for tic in np.asarray(matched["TIC"], dtype=np.int64):
        if replace or not (out_dir / f"{int(tic)}.h5").is_file():
            keep.append(True)
        else:
            keep.append(False)
            skipped_existing += 1
    matched = matched[np.asarray(keep, dtype=bool)]
    if len(matched) == 0:
        return {
            "source": source_path.name,
            "matched": n_matched,
            "requested": 0,
            "written": 0,
            "skipped_existing": skipped_existing,
            "errors": [],
        }

    source.tic = matched
    epsf = np.load(epsf_path)
    written = 0
    skipped_race = 0
    errors: list[str] = []
    for light_curve in generate_light_curves(source, epsf, psf_size, oversample, tic_ids=None):
        tic = int(light_curve.meta["tic_id"])
        final_path = out_dir / f"{tic}.h5"
        tmp_path = out_dir / f".{tic}.h5.tmp.{os.getpid()}"
        try:
            light_curve.write_hdf5(tmp_path)
            if install_hdf5(tmp_path, final_path, replace=replace):
                written += 1
            else:
                skipped_race += 1
        except Exception as exc:  # noqa: BLE001 - capture per-target failures for production logs.
            tmp_path.unlink(missing_ok=True)
            errors.append(f"TIC {tic}: {type(exc).__name__}: {exc}")

    return {
        "source": source_path.name,
        "matched": n_matched,
        "requested": len(matched),
        "written": written,
        "skipped_existing": skipped_existing + skipped_race,
        "errors": errors,
    }


def write_manifest(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run_lightcurves_for_ccd(
    *,
    input_root: Path,
    out_root: Path,
    orbit: int,
    camera: int,
    ccd: int,
    tic_catalog_path: Path,
    workers: int,
    replace: bool,
    psf_size: int,
    oversample: int,
    tmag_max: float,
    limit_sources: int | None,
    max_targets_per_source: int | None,
) -> dict[str, Any]:
    pairs = source_epsf_pairs(
        input_root=input_root,
        orbit=orbit,
        camera=camera,
        ccd=ccd,
        limit_sources=limit_sources,
    )
    out_lc_dir = lc_dir(out_root, orbit, camera, ccd)
    out_lc_dir.mkdir(parents=True, exist_ok=True)
    print(
        f"[lightcurves] orbit={orbit} cam={camera} ccd={ccd} "
        f"sources={len(pairs):,} workers={workers} out={out_lc_dir}",
        flush=True,
    )
    if not pairs:
        return {"orbit": orbit, "camera": camera, "ccd": ccd, "sources": 0, "written": 0}

    jobs = [
        (str(source_file), str(epsf_file), str(out_lc_dir), replace, psf_size, oversample)
        for source_file, epsf_file in pairs
    ]
    totals = {
        "orbit": orbit,
        "camera": camera,
        "ccd": ccd,
        "sources": len(jobs),
        "matched": 0,
        "requested": 0,
        "written": 0,
        "skipped_existing": 0,
        "error_count": 0,
        "errors": [],
    }
    t0 = time.time()
    if workers <= 1:
        init_worker(str(tic_catalog_path), tmag_max, max_targets_per_source)
        iterator = map(process_source_job, jobs)
        for i, result in enumerate(iterator, start=1):
            accumulate_result(totals, result)
            report_progress(i, len(jobs), totals, t0)
    else:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=init_worker,
            initargs=(str(tic_catalog_path), tmag_max, max_targets_per_source),
        ) as executor:
            futures = [executor.submit(process_source_job, job) for job in jobs]
            for i, future in enumerate(as_completed(futures), start=1):
                accumulate_result(totals, future.result())
                report_progress(i, len(jobs), totals, t0)
    return totals


def accumulate_result(totals: dict[str, Any], result: dict[str, Any]) -> None:
    totals["matched"] += int(result.get("matched", 0))
    totals["requested"] += int(result.get("requested", 0))
    totals["written"] += int(result.get("written", 0))
    totals["skipped_existing"] += int(result.get("skipped_existing", 0))
    errors = result.get("errors", [])
    totals["error_count"] += len(errors)
    if errors and len(totals["errors"]) < 20:
        totals["errors"].extend(errors[: 20 - len(totals["errors"])])


def report_progress(i: int, total: int, totals: dict[str, Any], t0: float) -> None:
    if i == total or i % 10 == 0:
        elapsed = max(time.time() - t0, 1e-6)
        rate = i / elapsed
        eta_min = (total - i) / rate / 60.0 if rate > 0 else float("nan")
        print(
            f"  [lightcurves] {i:,}/{total:,} sources "
            f"written={totals['written']:,} skipped={totals['skipped_existing']:,} "
            f"errors={totals['error_count']:,} eta={eta_min:.1f} min",
            flush=True,
        )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--orbits", type=int, nargs="+", required=True)
    ap.add_argument("--input-root", type=Path, default=Path("/pdo/qlp-data"))
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--ccd", default="all", help="'all' or comma-separated camera-ccd entries like 1-1,1-2")
    ap.add_argument("--stages", default="catalogs,lightcurves")
    ap.add_argument("--tmag-max", type=float, default=18.5)
    ap.add_argument("--catalog-nprocs", type=int, default=16)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--psf-size", type=int, default=11)
    ap.add_argument("--oversample", type=int, default=2)
    ap.add_argument("--replace-catalogs", action="store_true")
    ap.add_argument("--replace-lightcurves", action="store_true")
    ap.add_argument("--limit-sources", type=int, default=None)
    ap.add_argument("--max-targets-per-source", type=int, default=None)
    args = ap.parse_args(argv)

    out_root = enforce_pdo_user_output(args.out_root)
    input_root = Path(args.input_root).expanduser().resolve(strict=False)
    ccds = parse_ccd(args.ccd)
    stages = parse_stages(args.stages)

    run_manifest: dict[str, Any] = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "sector": args.sector,
        "orbits": args.orbits,
        "input_root": str(input_root),
        "out_root": str(out_root),
        "pdo_user_root": str(PDO_USER_ROOT),
        "tmag_max": args.tmag_max,
        "ccds": [f"{camera}-{ccd}" for camera, ccd in ccds],
        "stages": sorted(stages),
        "catalog_nprocs": args.catalog_nprocs,
        "workers": args.workers,
        "psf_size": args.psf_size,
        "oversample": args.oversample,
        "results": [],
    }

    for orbit in args.orbits:
        for camera, ccd in ccds:
            tic_path = catalog_path(out_root, orbit, camera, ccd)
            if "catalogs" in stages:
                table = ensure_tic_catalog(
                    path=tic_path,
                    orbit=orbit,
                    camera=camera,
                    ccd=ccd,
                    tmag_max=args.tmag_max,
                    nprocs=args.catalog_nprocs,
                    replace=args.replace_catalogs,
                )
                finite_tmag = np.asarray(table["tmag"], dtype=float)
                n_faint = int(np.count_nonzero((finite_tmag > 18.0) & (finite_tmag <= args.tmag_max)))
                print(
                    f"[catalog] orbit={orbit} cam={camera} ccd={ccd} "
                    f"rows={len(table):,} 18<Tmag<={args.tmag_max}: {n_faint:,}",
                    flush=True,
                )
            if "lightcurves" in stages:
                if not tic_path.is_file():
                    raise FileNotFoundError(f"Missing TIC catalog for lightcurves: {tic_path}")
                result = run_lightcurves_for_ccd(
                    input_root=input_root,
                    out_root=out_root,
                    orbit=orbit,
                    camera=camera,
                    ccd=ccd,
                    tic_catalog_path=tic_path,
                    workers=args.workers,
                    replace=args.replace_lightcurves,
                    psf_size=args.psf_size,
                    oversample=args.oversample,
                    tmag_max=args.tmag_max,
                    limit_sources=args.limit_sources,
                    max_targets_per_source=args.max_targets_per_source,
                )
                run_manifest["results"].append(result)
                manifest_path = (
                    out_root
                    / "manifests"
                    / f"s{args.sector:04d}_orbit{orbit}_cam{camera}_ccd{ccd}_faint_allstar.json"
                )
                write_manifest(manifest_path, {**run_manifest, "last_result": result})
                print(f"[manifest] wrote {manifest_path}", flush=True)

    run_manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    summary_path = out_root / "manifests" / f"s{args.sector:04d}_faint_allstar_summary.json"
    write_manifest(summary_path, run_manifest)
    print(f"[manifest] wrote {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
