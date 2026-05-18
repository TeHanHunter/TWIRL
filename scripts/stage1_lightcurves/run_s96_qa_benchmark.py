#!/usr/bin/env python3
"""Run the Sector 96 TGLC-TICA vs QLP QA benchmark on PDO.

For a short list of known transit hosts (reports/stage1_lightcurves/s96_qa/targets.json), this
driver:

  1. Resolves each target's S96 (camera, ccd) via tess-point.
  2. Scans /pdo/qlp-data/tica-delivery/s0096/camC-ccdD/ to determine the
     S96 orbit numbers and their cadence ranges, and writes an orbit-map
     JSON consumed by detrend_wrapper.py and hlsp_wrapper.py.
  3. For each unique (orbit, cam, ccd) needed, runs tglc catalogs ->
     cutouts -> epsfs -> lightcurves with lightcurves restricted to the
     target TIC list. Reuses _process_one_ccd from run_tglc_orbit_pipeline.
  4. Runs qlp lctools detrend per orbit (via detrend_wrapper.py).
  5. Runs qlp lctools hlsp at sector level (via hlsp_wrapper.py) with
     --flag-type tica (S >= 67).
  6. Copies the resulting HLSP FITS into the user-provided output dir
     so the local benchmark plotter can pick them up.

All writes are constrained to /pdo/users/tehan/...

Typical invocation (inside tmux on pdogpu1):

    cd /pdo/users/tehan/TWIRL
    source scripts/activate_qlp_env.sh    # only needed if detrend/hlsp share env
    python scripts/stage1_lightcurves/run_s96_qa_benchmark.py \
        --targets-json reports/stage1_lightcurves/s96_qa/targets.json \
        --tglc-data-dir /pdo/users/tehan/tglc-deep-catalogs-s96qa \
        --hlsp-out-dir /pdo/users/tehan/tglc-deep-catalogs-s96qa/hlsp_s0096 \
        --max-parallel-ccd-jobs 2 --epsfs-nprocs 16

The script does NOT download QLP HLSPs from MAST or generate plots; that
is done locally by reports/exploratory/_scripts/plot_tglc_tica_vs_qlp_s96.py after scping the
produced HLSP FITS back.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_tglc_orbit_pipeline import (  # noqa: E402
    CcdJob,
    USER_WRITE_ROOT,
    _build_runtime_env,
    _process_one_ccd,
    _require_user_owned_write_path,
    DEFAULT_FORK_PATH,
    DEFAULT_LD_PREFIX,
    DEFAULT_PYTHON_BIN,
    DEFAULT_TICA_ROOT,
)

DEFAULT_QLP_PYTHON = Path("/pdo/app/qlp-environment/.venv/bin/python")
DEFAULT_TWIRL_ROOT = Path("/pdo/users/tehan/TWIRL")


@dataclass(frozen=True)
class Target:
    slug: str
    display: str
    tic: int
    ra_deg: float
    dec_deg: float
    tmag: float
    camera: int
    ccd: int


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets-json", type=Path, required=True)
    parser.add_argument("--sector", type=int, default=96)
    parser.add_argument(
        "--tglc-data-dir",
        type=Path,
        default=Path("/pdo/users/tehan/tglc-deep-catalogs-s96qa"),
    )
    parser.add_argument("--tica-root", type=Path, default=DEFAULT_TICA_ROOT)
    parser.add_argument("--python-bin", type=Path, default=DEFAULT_PYTHON_BIN)
    parser.add_argument("--qlp-python-bin", type=Path, default=DEFAULT_QLP_PYTHON)
    parser.add_argument("--fork-path", type=Path, default=DEFAULT_FORK_PATH)
    parser.add_argument("--ld-library-prefix", default=DEFAULT_LD_PREFIX)
    parser.add_argument("--twirl-root", type=Path, default=DEFAULT_TWIRL_ROOT)
    parser.add_argument(
        "--hlsp-out-dir",
        type=Path,
        default=Path("/pdo/users/tehan/tglc-deep-catalogs-s96qa/hlsp_s0096"),
    )
    parser.add_argument("--max-parallel-ccd-jobs", type=int, default=2)
    parser.add_argument("--catalogs-nprocs", type=int, default=16)
    parser.add_argument("--cutouts-nprocs", type=int, default=16)
    parser.add_argument("--epsfs-nprocs", type=int, default=16)
    parser.add_argument("--lightcurves-nprocs", type=int, default=16)
    parser.add_argument("--detrend-nprocs", type=int, default=32)
    parser.add_argument("--hlsp-nprocs", type=int, default=32)
    parser.add_argument("--max-magnitude", type=float, default=16.5)
    parser.add_argument(
        "--flag-type",
        default="tica",
        help="qlp lctools hlsp --flag-type. S >= 67 uses 'tica'; S < 67 'spoc'.",
    )
    parser.add_argument(
        "--flag-source",
        default="fits",
        help="qlp lctools hlsp --flag-source.",
    )
    parser.add_argument("--no-gpu", action="store_true")
    parser.add_argument(
        "--skip-tglc", action="store_true", help="Skip tglc stages (debug)."
    )
    parser.add_argument(
        "--skip-detrend", action="store_true", help="Skip qlp detrend stage."
    )
    parser.add_argument(
        "--skip-hlsp", action="store_true", help="Skip qlp hlsp stage."
    )
    return parser.parse_args()


def _resolve_ccds_with_tesspoint(
    targets_raw: list[dict], sector: int
) -> list[Target]:
    """Resolve each target's (cam, ccd) for the requested sector via tess-point."""
    from tess_stars2px import tess_stars2px_function_entry

    tics = [int(t["tic"]) for t in targets_raw]
    ras = [float(t["ra_deg"]) for t in targets_raw]
    decs = [float(t["dec_deg"]) for t in targets_raw]
    (
        out_id,
        out_ra,
        out_dec,
        out_sector,
        out_cam,
        out_ccd,
        *_,
    ) = tess_stars2px_function_entry(tics, ras, decs)

    resolved: list[Target] = []
    for raw in targets_raw:
        tic = int(raw["tic"])
        match = None
        for oid, osec, ocam, occd in zip(out_id, out_sector, out_cam, out_ccd):
            if int(oid) == tic and int(osec) == sector:
                match = (int(ocam), int(occd))
                break
        if match is None:
            raise RuntimeError(
                f"tess-point: TIC {tic} not observed in S{sector}"
            )
        resolved.append(
            Target(
                slug=raw["slug"],
                display=raw["display"],
                tic=tic,
                ra_deg=float(raw["ra_deg"]),
                dec_deg=float(raw["dec_deg"]),
                tmag=float(raw["tmag"]),
                camera=match[0],
                ccd=match[1],
            )
        )
    return resolved


_FFI_RE = re.compile(
    r"hlsp_tica_tess_ffi_s(?P<sector>\d{4})-o(?P<orbit_tag>\d+)-"
    r"(?P<cadence>\d+)-cam\d-ccd\d_tess_v\d+_img\.fits$"
)


def _scan_s96_orbits(
    tica_root: Path, sector: int
) -> dict[int, dict]:
    """Return {orbit_tag_int: {'cadence_range': (min, max), 'orbit_tag': 'o<N>'}}.

    Scans one CCD dir (cam1-ccd1) — all CCDs in a sector share orbit_tag and
    cadence range."""
    probe_dir = tica_root / f"s{sector:04d}" / "cam1-ccd1"
    if not probe_dir.exists():
        raise FileNotFoundError(f"Missing probe dir: {probe_dir}")
    by_tag: dict[int, list[int]] = {}
    for f in probe_dir.glob("hlsp_tica_tess_ffi_*_img.fits"):
        m = _FFI_RE.search(f.name)
        if not m:
            continue
        if int(m.group("sector")) != sector:
            continue
        tag = int(m.group("orbit_tag"))
        by_tag.setdefault(tag, []).append(int(m.group("cadence")))
    if not by_tag:
        raise RuntimeError(
            f"No TICA FFI files matched under {probe_dir} for S{sector}"
        )
    return {
        tag: {
            "orbit_tag": f"o{tag}",
            "cadence_range": (min(cads), max(cads)),
        }
        for tag, cads in by_tag.items()
    }


def _compute_orbit_numbers(sector: int, orbit_tags: list[int]) -> dict[int, int]:
    """Map o1/o2/... to global orbit numbers.

    Based on S56 -> (119, 120) (i.e. orbit = 2*sector + 7 for o1). S97/S98 are
    known 4-orbit exceptions (not relevant for S96). For S96 this gives
    (199, 200)."""
    base = 2 * sector + 7
    return {tag: base + (tag - 1) for tag in sorted(orbit_tags)}


def _write_orbit_map_json(
    path: Path, sector: int, orbits: dict[int, dict], orbit_numbers: dict[int, int]
) -> None:
    entries = [
        {
            "orbit_number": orbit_numbers[tag],
            "cadence_range": list(orbits[tag]["cadence_range"]),
            "mjd_range": [0.0, 0.0],
        }
        for tag in sorted(orbits.keys())
    ]
    path.write_text(json.dumps({str(sector): entries}, indent=2) + "\n")


def _activate_qlp_env(twirl_root: Path) -> dict[str, str]:
    """Replicate scripts/activate_qlp_env.sh in env form."""
    env = dict(os.environ)
    # PYTHONPATH must point to TGLC fork only — FFITools' qlp 0.1 on the system
    # PYTHONPATH otherwise shadows the pip-installed qlp 0.13.2.
    env["PYTHONPATH"] = str(DEFAULT_FORK_PATH)
    existing_ld = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = (
        f"{DEFAULT_LD_PREFIX}:{existing_ld}" if existing_ld else DEFAULT_LD_PREFIX
    )
    return env


def _run_detrend(
    *,
    orbit: int,
    twirl_root: Path,
    tglc_data_dir: Path,
    qlp_python: Path,
    env: dict[str, str],
    ccd_jobs: list[CcdJob],
    nprocs: int,
    log_path: Path,
) -> None:
    _require_user_owned_write_path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path = tglc_data_dir / f"orbit-{orbit}" / "ffi" / "run" / "qlp.cfg"
    autolist = _autolist_arg(ccd_jobs)
    cmd = [
        str(qlp_python),
        str(twirl_root / "scripts" / "detrend_wrapper.py"),
        "lctools", "detrend",
        "-c", str(cfg_path),
        "--autolist", *autolist,
        "--best",
        "-n", str(nprocs),
        "--logfile", str(log_path),
    ]
    print(f"[s96-qa] detrend orbit={orbit}: {shlex.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, env=env, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"detrend failed for orbit {orbit}; see {log_path}"
        )


def _autolist_arg(ccd_jobs: list[CcdJob]) -> list[str]:
    by_cam: dict[int, list[int]] = {}
    for job in ccd_jobs:
        by_cam.setdefault(job.camera, []).append(job.ccd)
    return [f"{cam}:{','.join(str(c) for c in sorted(ccds))}"
            for cam, ccds in sorted(by_cam.items())]


def _run_hlsp(
    *,
    sector: int,
    orbit: int,
    twirl_root: Path,
    tglc_data_dir: Path,
    qlp_python: Path,
    env: dict[str, str],
    nprocs: int,
    hlsp_out_dir: Path,
    flag_type: str,
    flag_source: str,
    log_path: Path,
) -> None:
    _require_user_owned_write_path(log_path)
    _require_user_owned_write_path(hlsp_out_dir)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    hlsp_out_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = tglc_data_dir / f"orbit-{orbit}" / "ffi" / "run" / "qlp.cfg"
    cmd = [
        str(qlp_python),
        str(twirl_root / "scripts" / "hlsp_wrapper.py"),
        "lctools", "hlsp",
        "-c", str(cfg_path),
        "-s", str(sector),
        "--autolist", "all",
        "--flag-type", flag_type,
        "--flag-source", flag_source,
        "--basedir", str(tglc_data_dir) + "/",
        "-o", str(hlsp_out_dir),
        "-n", str(nprocs),
    ]
    print(f"[s96-qa] hlsp sector={sector}: {shlex.join(cmd)}", flush=True)
    with log_path.open("w") as f:
        proc = subprocess.run(
            cmd, env=env, stdout=f, stderr=subprocess.STDOUT, check=False
        )
    if proc.returncode != 0:
        raise RuntimeError(f"hlsp failed; see {log_path}")


def _ensure_qlp_cfg(
    tglc_data_dir: Path, orbit: int, sector: int, orbit_tag: str
) -> Path:
    """Minimal qlp.cfg writer so detrend/hlsp can find the data tree.

    We don't have a DB; the wrappers patch orbit/sector lookup. The cfg only
    needs to point IOSettings at the correct datadir."""
    run_dir = tglc_data_dir / f"orbit-{orbit}" / "ffi" / "run"
    _require_user_owned_write_path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = run_dir / "qlp.cfg"
    if cfg_path.exists():
        return cfg_path
    cfg_path.write_text(
        "[IOSettings]\n"
        f"datadir = {tglc_data_dir}/\n"
        f"sector = {sector}\n"
        f"orbit_id = {orbit}\n"
        f"orbit_tag = {orbit_tag}\n"
    )
    return cfg_path


def main() -> int:
    args = _parse_args()
    for path in [args.tglc_data_dir, args.hlsp_out_dir]:
        _require_user_owned_write_path(path)
    args.tglc_data_dir.mkdir(parents=True, exist_ok=True)

    targets_raw = json.loads(args.targets_json.read_text())["targets"]
    print(f"[s96-qa] resolving {len(targets_raw)} targets via tess-point", flush=True)
    targets = _resolve_ccds_with_tesspoint(targets_raw, args.sector)
    for t in targets:
        print(f"  {t.slug:10s} TIC {t.tic:11d} T={t.tmag:5.2f} "
              f"-> cam{t.camera}/ccd{t.ccd}", flush=True)

    orbits = _scan_s96_orbits(args.tica_root, args.sector)
    orbit_numbers = _compute_orbit_numbers(args.sector, list(orbits.keys()))
    print(f"[s96-qa] S{args.sector} orbits: "
          + ", ".join(f"{orbit_numbers[t]} (o{t}, cadence {orbits[t]['cadence_range']})"
                      for t in sorted(orbits.keys())),
          flush=True)

    orbit_map_path = args.tglc_data_dir / "s96_orbit_map.json"
    _require_user_owned_write_path(orbit_map_path)
    _write_orbit_map_json(orbit_map_path, args.sector, orbits, orbit_numbers)

    ccd_groups: dict[tuple[int, int], list[Target]] = {}
    for t in targets:
        ccd_groups.setdefault((t.camera, t.ccd), []).append(t)
    ccd_jobs = [CcdJob(camera=c, ccd=d) for (c, d) in sorted(ccd_groups.keys())]

    runtime_env = _build_runtime_env(
        fork_path=args.fork_path,
        ld_library_prefix=args.ld_library_prefix,
    )
    runtime_env["TWIRL_SECTOR_ORBIT_JSON"] = str(orbit_map_path)

    log_root = args.tglc_data_dir / "twirl_logs" / f"s{args.sector:04d}_qa"
    log_root.mkdir(parents=True, exist_ok=True)

    if not args.skip_tglc:
        for tag_int in sorted(orbits.keys()):
            orbit = orbit_numbers[tag_int]
            orbit_tag = orbits[tag_int]["orbit_tag"]
            print(f"\n[s96-qa] === TGLC orbit {orbit} ({orbit_tag}) ===",
                  flush=True)
            for job in ccd_jobs:
                tic_list = [t.tic for t in ccd_groups[(job.camera, job.ccd)]]
                _process_one_ccd(
                    job=job,
                    orbit=orbit,
                    sector=args.sector,
                    orbit_tag=orbit_tag,
                    tglc_data_dir=args.tglc_data_dir,
                    tica_root=args.tica_root,
                    python_bin=args.python_bin,
                    runtime_env=runtime_env,
                    log_root=log_root,
                    tic_list=tic_list,
                    max_magnitude=args.max_magnitude,
                    catalogs_nprocs=args.catalogs_nprocs,
                    cutouts_nprocs=args.cutouts_nprocs,
                    epsfs_nprocs=args.epsfs_nprocs,
                    lightcurves_nprocs=args.lightcurves_nprocs,
                    stop_on_warning=False,
                    reuse_catalogs_from_orbit=None,
                    no_gpu=args.no_gpu,
                )

    # qlp.cfg per orbit and quat symlinks are prerequisites for detrend/hlsp.
    for tag_int in sorted(orbits.keys()):
        orbit = orbit_numbers[tag_int]
        orbit_tag = orbits[tag_int]["orbit_tag"]
        _ensure_qlp_cfg(args.tglc_data_dir, orbit, args.sector, orbit_tag)
        # Link cam<C>_quat.txt for each cam used (QLP expects them in run/).
        cams_needed = sorted({job.camera for job in ccd_jobs})
        run_dir = args.tglc_data_dir / f"orbit-{orbit}" / "ffi" / "run"
        for cam in cams_needed:
            src = Path(f"/pdo/qlp-data/orbit-{orbit}/ffi/run/cam{cam}_quat.txt")
            dst = run_dir / f"cam{cam}_quat.txt"
            if not dst.exists() and src.exists():
                dst.symlink_to(src)
            elif not src.exists():
                print(f"[s96-qa] WARNING: missing {src}", flush=True)

    qlp_env = _activate_qlp_env(args.twirl_root)
    qlp_env["TWIRL_SECTOR_ORBIT_JSON"] = str(orbit_map_path)

    if not args.skip_detrend:
        for tag_int in sorted(orbits.keys()):
            orbit = orbit_numbers[tag_int]
            _run_detrend(
                orbit=orbit,
                twirl_root=args.twirl_root,
                tglc_data_dir=args.tglc_data_dir,
                qlp_python=args.qlp_python_bin,
                env=qlp_env,
                ccd_jobs=ccd_jobs,
                nprocs=args.detrend_nprocs,
                log_path=log_root / f"detrend_orbit-{orbit}.log",
            )

    if not args.skip_hlsp:
        # HLSP runs at sector level; pick the first orbit's qlp.cfg (they are
        # equivalent for our layout).
        first_orbit = orbit_numbers[min(orbits.keys())]
        _run_hlsp(
            sector=args.sector,
            orbit=first_orbit,
            twirl_root=args.twirl_root,
            tglc_data_dir=args.tglc_data_dir,
            qlp_python=args.qlp_python_bin,
            env=qlp_env,
            nprocs=args.hlsp_nprocs,
            hlsp_out_dir=args.hlsp_out_dir,
            flag_type=args.flag_type,
            flag_source=args.flag_source,
            log_path=log_root / f"hlsp_s{args.sector:04d}.log",
        )

    summary = {
        "sector": args.sector,
        "targets": [
            {
                "slug": t.slug, "tic": t.tic, "tmag": t.tmag,
                "camera": t.camera, "ccd": t.ccd,
            }
            for t in targets
        ],
        "orbit_numbers": {str(tag): orbit_numbers[tag] for tag in orbits},
        "hlsp_out_dir": str(args.hlsp_out_dir),
    }
    summary_path = log_root / "s96_qa_summary.json"
    _require_user_owned_write_path(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"\n[s96-qa] DONE. Summary: {summary_path}", flush=True)
    print(f"[s96-qa] HLSP FITS: {args.hlsp_out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
