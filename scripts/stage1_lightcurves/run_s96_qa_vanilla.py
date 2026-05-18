#!/usr/bin/env python3
"""Run vanilla TGLC (TeHanHunter/TESS_Gaia_Light_Curve, mitpdo branch) on S96
TICA FFIs for a small list of transit hosts, to compare against QLP HLSPs.

Entirely independent of the MIT QLP-integrated TGLC fork. Uses:

    tglc.ffi.ffi_tica(ccd, camera, sector, size=150, local_directory)
        -> stacks all TICA FFIs under /pdo/qlp-data/tica-delivery/s{NNNN}/camC-ccdD/
           into 196 source_XX_YY.pkl patches under local_directory/source/C-D/

    tglc.run.lc_per_ccd(local_directory, cam, ccd, cores)
        -> iterates lc_per_cut over those 196 patches, writing one FITS per
           Gaia DR3 source into local_directory/lc/C-D/hlsp_tglc_*.fits

The caller must import the mitpdo TGLC checkout via PYTHONPATH; this script
does NOT add it to sys.path so the operator stays in control of which tglc
package is used.

Runs from /pdo/users/tehan/TWIRL on pdogpu1 after:

    export PYTHONPATH=/pdo/users/tehan/TESS_Gaia_Light_Curve
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

No QLP lctools. No detrend. No HLSP stitching. Just vanilla TGLC photometry.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time as _time
from pathlib import Path

USER_WRITE_ROOT = Path("/pdo/users/tehan")


def _install_gaia_retry_shim(max_attempts: int = 6, base_delay: float = 5.0) -> None:
    """Wrap astroquery Gaia.launch_job and launch_job_async to retry on
    transient HTTP errors. Vanilla tglc.ffi.Source / convert_gaia_id make
    raw Gaia TAP calls with no retry logic, so a single HTTP 500 from the
    Gaia archive crashes ffi_tica mid-way through the 196-patch loop.
    """
    from astroquery.gaia import Gaia
    from requests.exceptions import HTTPError, ConnectionError, Timeout

    retryable = (HTTPError, ConnectionError, Timeout)

    for name in ("launch_job_async", "launch_job"):
        original = getattr(Gaia, name)

        def make_wrapper(orig, fname):
            def wrapper(*args, **kwargs):
                last_exc = None
                for attempt in range(1, max_attempts + 1):
                    try:
                        return orig(*args, **kwargs)
                    except retryable as exc:
                        last_exc = exc
                        delay = base_delay * (2 ** (attempt - 1))
                        print(
                            f"[gaia-retry] {fname} attempt {attempt}/{max_attempts} "
                            f"failed ({type(exc).__name__}: {exc}); sleeping {delay:.0f}s",
                            flush=True,
                        )
                        _time.sleep(delay)
                raise last_exc

            return wrapper

        setattr(Gaia, name, make_wrapper(original, name))
    print(
        f"[gaia-retry] installed retry shim on Gaia.launch_job[_async] "
        f"(max_attempts={max_attempts}, base_delay={base_delay}s)",
        flush=True,
    )


def _require_user_owned_write_path(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if USER_WRITE_ROOT not in resolved.parents and resolved != USER_WRITE_ROOT:
        raise SystemExit(
            f"refusing to write outside {USER_WRITE_ROOT}: {resolved}"
        )
    return resolved


def _resolve_ccds(targets: list[dict]) -> dict[tuple[int, int], list[dict]]:
    """Group targets by (camera, ccd) using tess-point for S96."""
    from tess_stars2px import tess_stars2px_function_entry as s2p

    tics = [int(t["tic"]) for t in targets]
    ras = [float(t["ra_deg"]) for t in targets]
    decs = [float(t["dec_deg"]) for t in targets]
    (out_id, _, _, sectors, cams, ccds, *_rest) = s2p(tics, ras, decs)

    by_ccd: dict[tuple[int, int], list[dict]] = {}
    for tic, sec, cam, ccd in zip(out_id, sectors, cams, ccds):
        if int(sec) != 96:
            continue
        tgt = next(t for t in targets if int(t["tic"]) == int(tic))
        tgt = {**tgt, "camera": int(cam), "ccd": int(ccd)}
        by_ccd.setdefault((int(cam), int(ccd)), []).append(tgt)
    return by_ccd


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--targets-json", type=Path, required=True)
    ap.add_argument(
        "--local-directory",
        type=Path,
        default=Path("/pdo/users/tehan/tglc-vanilla-s96qa"),
        help="Vanilla TGLC working directory. Must be under /pdo/users/tehan.",
    )
    ap.add_argument("--sector", type=int, default=96)
    ap.add_argument("--size", type=int, default=150)
    ap.add_argument(
        "--cores",
        type=int,
        default=16,
        help="Processes for lc_per_ccd. Cap BLAS threads to 1 via env.",
    )
    ap.add_argument(
        "--skip-ffi-tica",
        action="store_true",
        help="Skip the ffi_tica stacking (useful when source/ pkls exist).",
    )
    ap.add_argument(
        "--skip-lc",
        action="store_true",
        help="Skip lc_per_ccd photometry (useful to regenerate source/ only).",
    )
    ap.add_argument(
        "--only-cam-ccd",
        action="append",
        default=None,
        help='Restrict to a single "C,D" pair; may be repeated.',
    )
    args = ap.parse_args()

    local_directory = _require_user_owned_write_path(args.local_directory)
    local_directory.mkdir(parents=True, exist_ok=True)
    # TGLC's code paths append f'...{local_directory}lc/...' with no separator,
    # so local_directory must end with a trailing slash for its f-strings.
    local_dir_str = str(local_directory).rstrip("/") + "/"

    targets = json.loads(args.targets_json.read_text())["targets"]
    by_ccd = _resolve_ccds(targets)
    if not by_ccd:
        raise SystemExit("no targets fall in sector 96 per tess-point")

    only = None
    if args.only_cam_ccd:
        only = {tuple(int(x) for x in p.split(",")) for p in args.only_cam_ccd}

    print(
        f"[s96-qa-vanilla] sector={args.sector} size={args.size} "
        f"local_directory={local_directory}",
        flush=True,
    )
    for (cam, ccd), group in sorted(by_ccd.items()):
        if only is not None and (cam, ccd) not in only:
            continue
        names = ", ".join(
            f"{t['slug']} (TIC {t['tic']}, T={t['tmag']})" for t in group
        )
        print(f"[s96-qa-vanilla] cam{cam}/ccd{ccd}: {names}", flush=True)

    # Install retry shim BEFORE importing tglc.ffi — Source() / convert_gaia_id
    # fire Gaia.launch_job_async with no retry, so a transient 500 would crash
    # ffi_tica mid-loop and waste the full ~2h of source-pkl work.
    _install_gaia_retry_shim()

    # Late imports so the vanilla TGLC PYTHONPATH is in effect.
    from tglc.ffi import ffi_tica
    from tglc.run import lc_per_ccd

    for (cam, ccd), group in sorted(by_ccd.items()):
        if only is not None and (cam, ccd) not in only:
            continue
        print(
            f"[s96-qa-vanilla] === cam{cam}/ccd{ccd} (S{args.sector}) ===",
            flush=True,
        )
        if not args.skip_ffi_tica:
            print(
                f"[s96-qa-vanilla] ffi_tica cam={cam} ccd={ccd} -> "
                f"{local_dir_str}source/{cam}-{ccd}/",
                flush=True,
            )
            # ffi_tica needs source/cam-ccd/ to exist.
            (local_directory / "source" / f"{cam}-{ccd}").mkdir(
                parents=True, exist_ok=True
            )
            (local_directory / "ffi").mkdir(parents=True, exist_ok=True)
            (local_directory / "lc" / f"{cam}-{ccd}").mkdir(
                parents=True, exist_ok=True
            )
            (local_directory / "epsf" / f"{cam}-{ccd}").mkdir(
                parents=True, exist_ok=True
            )
            (local_directory / "log").mkdir(parents=True, exist_ok=True)
            ffi_tica(
                ccd=ccd,
                camera=cam,
                sector=args.sector,
                size=args.size,
                local_directory=local_dir_str,
            )
        if not args.skip_lc:
            print(
                f"[s96-qa-vanilla] lc_per_ccd cam={cam} ccd={ccd} cores={args.cores}",
                flush=True,
            )
            lc_per_ccd(
                local_directory=local_dir_str,
                cam=cam,
                ccd=ccd,
                cores=args.cores,
            )
        print(f"[s96-qa-vanilla] cam{cam}/ccd{ccd} done.", flush=True)

    print("[s96-qa-vanilla] all requested CCDs complete.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
