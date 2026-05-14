#!/usr/bin/env python3
"""Run one-orbit TGLC production for a selected CCD set.

This driver is designed for the Sector 56 benchmark on PDO:

- stage read-only TICA FFIs from `/pdo/qlp-data/...` into a user-owned
  orbit tree via symlinks,
- run `tglc catalogs -> cutouts -> epsfs -> lightcurves` sequentially per CCD,
- parallelize over multiple independent CCD jobs at the outer level,
- keep `epsfs --nprocs` fixed per CCD so outer/inner parallelism tests remain
  interpretable,
- extract only TWIRL WD TIC IDs at the `lightcurves` stage.

All writes are constrained to `/pdo/users/tehan/...`.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits


DEFAULT_TGLC_DATA_DIR = Path("/pdo/users/tehan/tglc-deep-catalogs")
DEFAULT_TICA_ROOT = Path("/pdo/qlp-data/tica-delivery")
DEFAULT_OBSERVATIONS_FITS = Path(
    "/pdo/users/tehan/TWIRL/data_local/catalogs/twirl_master_catalog/"
    "twirl_wd_tess_observations_v0.fits"
)
DEFAULT_PYTHON_BIN = Path("/sw/qlp-environment/.venv/bin/python")
DEFAULT_FORK_PATH = Path("/pdo/users/tehan/tess-gaia-light-curve-twirl")
DEFAULT_LD_PREFIX = (
    "/pdo/app/anaconda/anaconda2-4.4.0/lib:"
    "/pdo/app/python-versions/python-3.11.9/lib"
)
DEFAULT_LOG_DIR = Path("/pdo/users/tehan/tglc-deep-catalogs/twirl_logs")
USER_WRITE_ROOT = Path("/pdo/users/tehan").resolve()


@dataclass(frozen=True)
class CcdJob:
    camera: int
    ccd: int

    @property
    def label(self) -> str:
        return f"cam{self.camera}_ccd{self.ccd}"


def _parse_ccd(value: str) -> CcdJob:
    parts = value.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"CCD must have format 'camera,ccd'; got {value!r}"
        )
    try:
        camera = int(parts[0])
        ccd = int(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"CCD must contain integer camera and CCD values; got {value!r}"
        ) from exc
    if camera not in {1, 2, 3, 4} or ccd not in {1, 2, 3, 4}:
        raise argparse.ArgumentTypeError(
            f"Camera/CCD must be in 1..4; got {value!r}"
        )
    return CcdJob(camera=camera, ccd=ccd)


def _default_ccd_jobs() -> list[CcdJob]:
    return [CcdJob(camera=camera, ccd=ccd) for camera in range(1, 5) for ccd in range(1, 5)]


def _require_user_owned_write_path(path: Path) -> None:
    if path.exists() and not path.is_symlink():
        resolved = path.resolve()
    else:
        resolved = path.parent.resolve() / path.name
    if resolved != USER_WRITE_ROOT and USER_WRITE_ROOT not in resolved.parents:
        raise ValueError(
            f"Refusing to write outside {USER_WRITE_ROOT}: {path}"
        )


def _build_runtime_env(
    *,
    fork_path: Path,
    ld_library_prefix: str,
) -> dict[str, str]:
    env = dict(os.environ)
    existing_ld = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = (
        f"{ld_library_prefix}:{existing_ld}" if existing_ld else ld_library_prefix
    )
    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        f"{fork_path}:{existing_pythonpath}" if existing_pythonpath else str(fork_path)
    )
    return env


def _load_tic_lists_by_ccd(
    *,
    observations_fits: Path,
    orbit: int,
    ccd_jobs: list[CcdJob],
) -> dict[CcdJob, list[int]]:
    with fits.open(observations_fits, memmap=True) as hdul:
        rows = hdul[1].data
        orbit_col = np.asarray(rows["orbit"])
        camera_col = np.asarray(rows["camera"])
        ccd_col = np.asarray(rows["ccd"])
        tic_col = np.asarray(rows["tic_id"])

        tic_lists: dict[CcdJob, list[int]] = {}
        orbit_mask = orbit_col == orbit
        for job in ccd_jobs:
            mask = (
                orbit_mask
                & (camera_col == job.camera)
                & (ccd_col == job.ccd)
                & (tic_col > 0)
            )
            tic_lists[job] = np.unique(tic_col[mask]).astype(np.int64).tolist()
    return tic_lists


def _stage_ffi_symlinks(
    *,
    sector: int,
    orbit_tag: str,
    camera: int,
    ccd: int,
    tica_root: Path,
    tglc_data_dir: Path,
    orbit: int,
) -> tuple[Path, int]:
    source_dir = tica_root / f"s{sector:04d}" / f"cam{camera}-ccd{ccd}"
    if not source_dir.exists():
        raise FileNotFoundError(f"Missing TICA source directory: {source_dir}")

    ffi_dir = tglc_data_dir / f"orbit-{orbit}" / "ffi" / f"cam{camera}" / f"ccd{ccd}" / "ffi"
    _require_user_owned_write_path(ffi_dir)
    ffi_dir.mkdir(parents=True, exist_ok=True)

    pattern = f"*-{orbit_tag}-*-cam{camera}-ccd{ccd}_tess_v01_img.fits"
    source_files = sorted(source_dir.glob(pattern))
    if not source_files:
        raise FileNotFoundError(
            f"No FFI files found under {source_dir} with pattern {pattern!r}"
        )

    for source_file in source_files:
        link_path = ffi_dir / source_file.name
        _require_user_owned_write_path(link_path)
        if link_path.is_symlink() or link_path.exists():
            if link_path.resolve() == source_file.resolve():
                continue
            if not link_path.is_symlink():
                raise FileExistsError(
                    f"Refusing to replace non-symlink file {link_path}"
                )
            link_path.unlink()
        link_path.symlink_to(source_file)

    return ffi_dir, len(source_files)


def _count_files(directory: Path, pattern: str) -> int:
    if not directory.exists():
        return 0
    return sum(1 for _ in directory.glob(pattern))


def _run_stage(
    *,
    stage_name: str,
    command: list[str],
    env: dict[str, str],
    log_path: Path,
) -> dict[str, object]:
    _require_user_owned_write_path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"START stage={stage_name}\n")
        log_handle.write(f"COMMAND {shlex.join(command)}\n")
        log_handle.flush()
        proc = subprocess.run(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            env=env,
        )
        log_handle.write(f"RETURN_CODE {proc.returncode}\n")
    wall_seconds = time.perf_counter() - start

    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    return {
        "stage": stage_name,
        "return_code": proc.returncode,
        "wall_seconds": wall_seconds,
        "log_path": str(log_path),
        "invalid_ffi_count": log_text.count("Invalid FFI file"),
        "runtimewarning_count": log_text.count("RuntimeWarning"),
        "cadence_gap_warning_count": log_text.count("cadence gaps != 1 detected"),
    }


def _reuse_sibling_catalogs(
    *,
    job: CcdJob,
    orbit: int,
    reuse_from_orbit: int,
    tglc_data_dir: Path,
    catalogs_dir: Path,
) -> list[str]:
    """Symlink Gaia/TIC catalogs from a sibling orbit so the catalogs stage can skip the query.

    Within a TESS sector, orbits point at the same sky, and TGLC's catalogs stage keys its
    Gaia and TIC queries on (sector, camera, ccd). Sibling orbits in the same sector produce
    identical catalog files, so reusing them avoids rerunning very slow dense-field queries.

    The caller is responsible for passing a `reuse_from_orbit` that is in the same sector.
    """
    source_catalogs_dir = tglc_data_dir / f"orbit-{reuse_from_orbit}" / "ffi" / "catalogs"
    symlinked: list[str] = []
    for basename in (f"Gaia_cam{job.camera}_ccd{job.ccd}.ecsv", f"TIC_cam{job.camera}_ccd{job.ccd}.ecsv"):
        target = catalogs_dir / basename
        if target.exists() or target.is_symlink():
            continue
        source = source_catalogs_dir / basename
        if not source.exists():
            continue
        _require_user_owned_write_path(target)
        target.symlink_to(source.resolve())
        symlinked.append(basename)
    return symlinked


def _process_one_ccd(
    *,
    job: CcdJob,
    orbit: int,
    sector: int,
    orbit_tag: str,
    tglc_data_dir: Path,
    tica_root: Path,
    python_bin: Path,
    runtime_env: dict[str, str],
    log_root: Path,
    tic_list: list[int],
    max_magnitude: float,
    catalogs_nprocs: int,
    cutouts_nprocs: int,
    epsfs_nprocs: int,
    lightcurves_nprocs: int,
    stop_on_warning: bool,
    reuse_catalogs_from_orbit: int | None,
    no_gpu: bool = False,
    gpu_id: str | None = None,
) -> dict[str, object]:
    job_start = time.perf_counter()
    if gpu_id is not None:
        runtime_env = {**runtime_env, "CUDA_VISIBLE_DEVICES": gpu_id}
    job_root = tglc_data_dir / f"orbit-{orbit}" / "ffi" / f"cam{job.camera}" / f"ccd{job.ccd}"
    catalogs_dir = tglc_data_dir / f"orbit-{orbit}" / "ffi" / "catalogs"
    source_dir = job_root / "source"
    epsf_dir = job_root / "epsf"
    lc_dir = job_root / "LC"
    tic_list_path = log_root / "tic_lists" / f"orbit-{orbit}_{job.label}_wd_tics.txt"

    for path in [job_root, catalogs_dir, source_dir, epsf_dir, lc_dir, tic_list_path.parent, log_root]:
        _require_user_owned_write_path(path)

    if reuse_catalogs_from_orbit is not None and reuse_catalogs_from_orbit != orbit:
        reused = _reuse_sibling_catalogs(
            job=job,
            orbit=orbit,
            reuse_from_orbit=reuse_catalogs_from_orbit,
            tglc_data_dir=tglc_data_dir,
            catalogs_dir=catalogs_dir,
        )
        if reused:
            print(
                f"[orbit-pipeline] {job.label}: reusing catalogs from orbit "
                f"{reuse_catalogs_from_orbit}: {', '.join(reused)}",
                flush=True,
            )

    print(f"[orbit-pipeline] {job.label}: staging FFIs", flush=True)
    ffi_dir, n_symlinked = _stage_ffi_symlinks(
        sector=sector,
        orbit_tag=orbit_tag,
        camera=job.camera,
        ccd=job.ccd,
        tica_root=tica_root,
        tglc_data_dir=tglc_data_dir,
        orbit=orbit,
    )

    tic_list_path.parent.mkdir(parents=True, exist_ok=True)
    tic_list_path.write_text(
        "".join(f"{tic_id}\n" for tic_id in tic_list),
        encoding="utf-8",
    )

    stage_records: list[dict[str, object]] = []
    command_prefix = [str(python_bin), "-m", "tglc"]

    catalogs_record = _run_stage(
        stage_name="catalogs",
        command=[
            *command_prefix,
            "catalogs",
            "--orbit",
            str(orbit),
            "--ccd",
            f"{job.camera},{job.ccd}",
            "--max-magnitude",
            str(max_magnitude),
            "--nprocs",
            str(catalogs_nprocs),
            "--tglc-data-dir",
            str(tglc_data_dir),
        ],
        env=runtime_env,
        log_path=log_root / f"orbit-{orbit}_{job.label}_catalogs.log",
    )
    stage_records.append(catalogs_record)
    if catalogs_record["return_code"] != 0:
        raise RuntimeError(
            f"catalogs failed for orbit {orbit} {job.label}; see {catalogs_record['log_path']}"
        )

    cutouts_record = _run_stage(
        stage_name="cutouts",
        command=[
            *command_prefix,
            "cutouts",
            "--orbit",
            str(orbit),
            "--ccd",
            f"{job.camera},{job.ccd}",
            "--nprocs",
            str(cutouts_nprocs),
            "--replace",
            "--tglc-data-dir",
            str(tglc_data_dir),
        ],
        env=runtime_env,
        log_path=log_root / f"orbit-{orbit}_{job.label}_cutouts.log",
    )
    stage_records.append(cutouts_record)
    if cutouts_record["return_code"] != 0:
        raise RuntimeError(
            f"cutouts failed for orbit {orbit} {job.label}; see {cutouts_record['log_path']}"
        )
    if stop_on_warning and (
        cutouts_record["invalid_ffi_count"]
        or cutouts_record["runtimewarning_count"]
        or cutouts_record["cadence_gap_warning_count"]
    ):
        raise RuntimeError(
            f"cutouts warning pattern for orbit {orbit} {job.label}; see {cutouts_record['log_path']}"
        )

    source_count = _count_files(source_dir, "source_*.pkl")
    if source_count != 196:
        raise RuntimeError(
            f"Expected 196 source pickles for orbit {orbit} {job.label}, got {source_count}"
        )

    epsfs_command = [
        *command_prefix,
        "epsfs",
        "--orbit",
        str(orbit),
        "--ccd",
        f"{job.camera},{job.ccd}",
        "--nprocs",
        str(epsfs_nprocs),
        "--replace",
        "--tglc-data-dir",
        str(tglc_data_dir),
    ]
    if no_gpu:
        epsfs_command.append("--no-gpu")
    elif gpu_id is not None:
        print(
            f"[orbit-pipeline] {job.label}: epsfs pinned to CUDA_VISIBLE_DEVICES={gpu_id}",
            flush=True,
        )
    epsfs_record = _run_stage(
        stage_name="epsfs",
        command=epsfs_command,
        env=runtime_env,
        log_path=log_root / f"orbit-{orbit}_{job.label}_epsfs.log",
    )
    stage_records.append(epsfs_record)
    if epsfs_record["return_code"] != 0:
        raise RuntimeError(
            f"epsfs failed for orbit {orbit} {job.label}; see {epsfs_record['log_path']}"
        )
    if stop_on_warning and epsfs_record["runtimewarning_count"]:
        raise RuntimeError(
            f"epsfs RuntimeWarning for orbit {orbit} {job.label}; see {epsfs_record['log_path']}"
        )

    epsf_count = _count_files(epsf_dir, "epsf_*.npy")
    if epsf_count != 196:
        raise RuntimeError(
            f"Expected 196 ePSFs for orbit {orbit} {job.label}, got {epsf_count}"
        )

    if tic_list:
        lightcurve_command = [
            *command_prefix,
            "lightcurves",
            "--orbit",
            str(orbit),
            "--ccd",
            f"{job.camera},{job.ccd}",
            "--nprocs",
            str(lightcurves_nprocs),
            "--replace",
            "--tglc-data-dir",
            str(tglc_data_dir),
            "--tic",
            *[str(tic_id) for tic_id in tic_list],
        ]
        lightcurves_record = _run_stage(
            stage_name="lightcurves",
            command=lightcurve_command,
            env=runtime_env,
            log_path=log_root / f"orbit-{orbit}_{job.label}_lightcurves.log",
        )
        stage_records.append(lightcurves_record)
        if lightcurves_record["return_code"] != 0:
            raise RuntimeError(
                f"lightcurves failed for orbit {orbit} {job.label}; see {lightcurves_record['log_path']}"
            )

    lc_count = _count_files(lc_dir, "*.h5")
    job_summary = {
        "orbit": orbit,
        "sector": sector,
        "orbit_tag": orbit_tag,
        "camera": job.camera,
        "ccd": job.ccd,
        "job_label": job.label,
        "ffi_dir": str(ffi_dir),
        "n_symlinked_ffi": n_symlinked,
        "n_requested_wd_tics": len(tic_list),
        "tic_list_path": str(tic_list_path),
        "n_source_pickles": source_count,
        "n_epsf_files": epsf_count,
        "n_lightcurve_h5": lc_count,
        "stage_records": stage_records,
        "job_wall_seconds": time.perf_counter() - job_start,
    }
    summary_path = log_root / f"orbit-{orbit}_{job.label}_summary.json"
    _require_user_owned_write_path(summary_path)
    summary_path.write_text(
        json.dumps(job_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(
        "[orbit-pipeline] "
        f"{job.label}: DONE "
        f"ffi={n_symlinked} tic={len(tic_list)} "
        f"source={source_count} epsf={epsf_count} lc={lc_count} "
        f"wall={job_summary['job_wall_seconds'] / 3600:.2f}h",
        flush=True,
    )
    return job_summary


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one-orbit TGLC catalogs->cutouts->epsfs->lightcurves jobs."
    )
    parser.add_argument("--orbit", type=int, required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument(
        "--orbit-tag",
        required=True,
        help="Orbit-within-sector tag in TICA filenames, e.g. o1 or o2.",
    )
    parser.add_argument(
        "--ccd",
        type=_parse_ccd,
        action="append",
        default=[],
        help="CCD to process, format 'camera,ccd'. Repeat to restrict the set.",
    )
    parser.add_argument(
        "--skip-ccd",
        type=_parse_ccd,
        action="append",
        default=[],
        help="CCD to skip, format 'camera,ccd'. Repeat as needed.",
    )
    parser.add_argument(
        "--max-parallel-ccd-jobs",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--catalogs-nprocs",
        type=int,
        default=16,
    )
    parser.add_argument(
        "--cutouts-nprocs",
        type=int,
        default=16,
    )
    parser.add_argument(
        "--epsfs-nprocs",
        type=int,
        default=16,
    )
    parser.add_argument(
        "--lightcurves-nprocs",
        type=int,
        default=16,
    )
    parser.add_argument(
        "--max-magnitude",
        type=float,
        default=20.0,
    )
    parser.add_argument(
        "--tglc-data-dir",
        type=Path,
        default=DEFAULT_TGLC_DATA_DIR,
    )
    parser.add_argument(
        "--tica-root",
        type=Path,
        default=DEFAULT_TICA_ROOT,
    )
    parser.add_argument(
        "--observations-fits",
        type=Path,
        default=DEFAULT_OBSERVATIONS_FITS,
    )
    parser.add_argument(
        "--python-bin",
        type=Path,
        default=DEFAULT_PYTHON_BIN,
    )
    parser.add_argument(
        "--fork-path",
        type=Path,
        default=DEFAULT_FORK_PATH,
    )
    parser.add_argument(
        "--ld-library-prefix",
        default=DEFAULT_LD_PREFIX,
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
    )
    parser.add_argument(
        "--run-label",
        default="",
    )
    parser.add_argument(
        "--stop-on-warning",
        action="store_true",
        help="Treat warning markers in stage logs as fatal.",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help=(
            "Pass --no-gpu to `tglc epsfs`. Use when CuPy is unavailable or when "
            "benchmarking CPU baseline."
        ),
    )
    parser.add_argument(
        "--gpu-list",
        default="",
        help=(
            "Comma-separated CUDA device IDs (e.g. '0,1,2,3,4,6') to cycle across "
            "concurrent CCD jobs. Each CCD subprocess gets CUDA_VISIBLE_DEVICES set "
            "to one entry. Empty (default) means no pinning (CuPy picks device 0 "
            "for every concurrent job — only safe with --max-parallel-ccd-jobs 1 "
            "on GPU, or with --no-gpu). Ignored when --no-gpu is set."
        ),
    )
    parser.add_argument(
        "--reuse-catalogs-from-orbit",
        type=int,
        default=None,
        help=(
            "Before the catalogs stage, symlink Gaia/TIC ecsv files from this sibling "
            "orbit's catalogs directory into the current orbit's catalogs directory. "
            "The catalogs stage will then skip the query via its native skip-if-exists "
            "check. ONLY pass an orbit in the same TESS sector as --orbit; cross-sector "
            "reuse would pull the wrong sky patch. Useful for dense fields where the "
            "Gaia query is flaky (e.g. sector 56 cam3/ccd2)."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()

    ccd_jobs = args.ccd if args.ccd else _default_ccd_jobs()
    skip_ccd_jobs = set(args.skip_ccd)
    ccd_jobs = [job for job in ccd_jobs if job not in skip_ccd_jobs]
    if not ccd_jobs:
        raise ValueError("No CCD jobs selected after applying --skip-ccd.")

    for path in [args.tglc_data_dir, args.log_dir]:
        _require_user_owned_write_path(path)
    if not args.python_bin.exists():
        raise FileNotFoundError(f"Missing Python executable: {args.python_bin}")
    if not args.fork_path.exists():
        raise FileNotFoundError(f"Missing TGLC fork checkout: {args.fork_path}")
    if not args.observations_fits.exists():
        raise FileNotFoundError(
            f"Missing TWIRL observation FITS: {args.observations_fits}"
        )

    run_label = args.run_label or f"s{args.sector:04d}_orbit{args.orbit}_p{args.max_parallel_ccd_jobs}"
    run_log_root = args.log_dir / run_label
    _require_user_owned_write_path(run_log_root)
    run_log_root.mkdir(parents=True, exist_ok=True)

    print(
        "[orbit-pipeline] "
        f"orbit={args.orbit} sector={args.sector} orbit_tag={args.orbit_tag} "
        f"ccd_jobs={len(ccd_jobs)} max_parallel={args.max_parallel_ccd_jobs} "
        f"catalogs_nprocs={args.catalogs_nprocs} cutouts_nprocs={args.cutouts_nprocs} "
        f"epsfs_nprocs={args.epsfs_nprocs} lightcurves_nprocs={args.lightcurves_nprocs}",
        flush=True,
    )
    print(
        f"[orbit-pipeline] loading WD TIC lists from {args.observations_fits}",
        flush=True,
    )
    tic_lists_by_ccd = _load_tic_lists_by_ccd(
        observations_fits=args.observations_fits,
        orbit=args.orbit,
        ccd_jobs=ccd_jobs,
    )

    manifest_path = run_log_root / "run_manifest.json"
    _require_user_owned_write_path(manifest_path)
    manifest_path.write_text(
        json.dumps(
            {
                "orbit": args.orbit,
                "sector": args.sector,
                "orbit_tag": args.orbit_tag,
                "max_parallel_ccd_jobs": args.max_parallel_ccd_jobs,
                "catalogs_nprocs": args.catalogs_nprocs,
                "cutouts_nprocs": args.cutouts_nprocs,
                "epsfs_nprocs": args.epsfs_nprocs,
                "lightcurves_nprocs": args.lightcurves_nprocs,
                "max_magnitude": args.max_magnitude,
                "tglc_data_dir": str(args.tglc_data_dir),
                "tica_root": str(args.tica_root),
                "observations_fits": str(args.observations_fits),
                "python_bin": str(args.python_bin),
                "fork_path": str(args.fork_path),
                "ld_library_prefix": args.ld_library_prefix,
                "ccd_jobs": [
                    {
                        "camera": job.camera,
                        "ccd": job.ccd,
                        "n_requested_wd_tics": len(tic_lists_by_ccd[job]),
                    }
                    for job in ccd_jobs
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    runtime_env = _build_runtime_env(
        fork_path=args.fork_path,
        ld_library_prefix=args.ld_library_prefix,
    )

    gpu_ids: list[str] = [g.strip() for g in args.gpu_list.split(",") if g.strip()]
    if gpu_ids and args.no_gpu:
        print(
            "[orbit-pipeline] --gpu-list ignored because --no-gpu is set",
            flush=True,
        )
        gpu_ids = []

    run_start = time.perf_counter()
    summaries: list[dict[str, object]] = []
    with cf.ThreadPoolExecutor(max_workers=args.max_parallel_ccd_jobs) as executor:
        futures = {
            executor.submit(
                _process_one_ccd,
                job=job,
                gpu_id=(gpu_ids[idx % len(gpu_ids)] if gpu_ids else None),
                no_gpu=args.no_gpu,
                orbit=args.orbit,
                sector=args.sector,
                orbit_tag=args.orbit_tag,
                tglc_data_dir=args.tglc_data_dir,
                tica_root=args.tica_root,
                python_bin=args.python_bin,
                runtime_env=runtime_env,
                log_root=run_log_root,
                tic_list=tic_lists_by_ccd[job],
                max_magnitude=args.max_magnitude,
                catalogs_nprocs=args.catalogs_nprocs,
                cutouts_nprocs=args.cutouts_nprocs,
                epsfs_nprocs=args.epsfs_nprocs,
                lightcurves_nprocs=args.lightcurves_nprocs,
                stop_on_warning=args.stop_on_warning,
                reuse_catalogs_from_orbit=args.reuse_catalogs_from_orbit,
            ): job
            for idx, job in enumerate(ccd_jobs)
        }

        for future in cf.as_completed(futures):
            job = futures[future]
            try:
                summaries.append(future.result())
            except Exception as exc:
                print(
                    f"[orbit-pipeline] FAILED {job.label}: {exc}",
                    flush=True,
                )
                raise

    run_wall_seconds = time.perf_counter() - run_start
    final_summary_path = run_log_root / "run_summary.json"
    _require_user_owned_write_path(final_summary_path)
    final_summary_path.write_text(
        json.dumps(
            {
                "orbit": args.orbit,
                "sector": args.sector,
                "orbit_tag": args.orbit_tag,
                "n_ccd_jobs": len(ccd_jobs),
                "max_parallel_ccd_jobs": args.max_parallel_ccd_jobs,
                "run_wall_seconds": run_wall_seconds,
                "summaries": sorted(
                    summaries,
                    key=lambda row: (int(row["camera"]), int(row["ccd"])),
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        "[orbit-pipeline] "
        f"ALL DONE n_ccd_jobs={len(ccd_jobs)} wall={run_wall_seconds / 3600:.2f}h "
        f"summary={final_summary_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
