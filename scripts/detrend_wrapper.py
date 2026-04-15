#!/usr/bin/env python3
"""Wrapper for qlp lctools detrend that patches orbit info to bypass DB lookup.

lightcurvedb 3.0.0 requires tess_orbit table not populated in stelvar.
This wrapper patches both call sites so detrend runs without a live database.

It also:
  * caps BLAS/OpenMP threads to 1 per worker (the upstream code silently
    spawns threads equal to core count, so `-n 8` becomes ~1000 threads on
    a 128-core box and everything thrashes)
  * wraps the per-task entry in try/except so one bad LC doesn't abort the
    pool
  * short-circuits LCs that already have `bestdmagkey` written, so reruns
    are O(1) instead of re-fitting every target

Usage on PDO (from /pdo/users/tehan/TWIRL/):
    LD_LIBRARY_PATH=/pdo/app/anaconda/anaconda2-4.4.0/lib:/pdo/app/python-versions/python-3.11.9/lib \\
    PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl \\
    /pdo/app/qlp-environment/.venv/bin/python detrend_wrapper.py \\
        lctools detrend \\
        -c /pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/run/qlp.cfg \\
        --autolist 1:1,2,3,4 2:1,2,3,4 3:1,2,3,4 4:1,2,3,4 \\
        --best -n 32 \\
        --logfile /pdo/users/tehan/TWIRL/logs/detrend_orbit119.log
"""

# Thread caps MUST be set before numpy/scipy are imported, so do this first.
import os as _os

for _v in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    _os.environ.setdefault(_v, "1")

import sys
from qlp.io.datastructures import OrbitalInfo

ORBIT_SECTOR_MAP = {119: 56, 120: 56}


def _mock_read_orbit_info(orbit_number: int) -> OrbitalInfo:
    sector = ORBIT_SECTOR_MAP.get(int(orbit_number), 56)
    return OrbitalInfo(
        orbit_number=int(orbit_number),
        sector=sector,
        cadence_range=(0, 0),
        mjd_range=(0.0, 0.0),
        basename="",
    )


# Patch 1: qlp.util.util._get_orbit_info (used by get_default_keys)
import qlp.util.util
qlp.util.util._get_orbit_info = _mock_read_orbit_info

# Patch 2: default_io_backend.read_orbit_info (used directly in detrend_qsp_h5)
import qlp.io.backends
qlp.io.backends.default_io_backend.read_orbit_info = _mock_read_orbit_info

# Patch 3: fault-tolerant + fast-skip per-task wrapper.
# - try/except: one bad LC (e.g. an off-by-one cadence mismatch between BJD
#   and quaternion features) would otherwise abort the whole pool.
# - fast-skip: detrend_qsp_h5 doesn't short-circuit when bestdmagkey is
#   already written; it still runs the expensive spline lstsq. Check here
#   and return immediately if the LC is already detrended.
import h5py
import qlp.lctools.bin.detrend as _qlp_detrend
_orig_detrend_wrapper = _qlp_detrend.detrend_wrapper


def _already_detrended(lc_file) -> bool:
    # A partial/crashed detrend can leave bestdmagkey set with no corresponding
    # magnitude dataset on disk; treat that as not-done so we redo the fit.
    try:
        with h5py.File(lc_file, "r") as f:
            ap = f["LightCurve"]["AperturePhotometry"]
            key = ap.attrs.get("bestdmagkey", None)
            if key is None:
                return False
            if isinstance(key, bytes):
                key = key.decode()
            primary = ap.get("PrimaryAperture") or ap.get("Aperture_000")
            if primary is None:
                return False
            return key in primary
    except Exception:
        return False


def _safe_detrend_wrapper(task):
    lc_file = task[0] if task else None
    if lc_file is not None and _already_detrended(lc_file):
        return None
    try:
        return _orig_detrend_wrapper(task)
    except Exception as e:
        print(f"[detrend skip] {lc_file}: {type(e).__name__}: {e}", flush=True)
        return None


_qlp_detrend.detrend_wrapper = _safe_detrend_wrapper

from qlp.__main__ import qlp_main
sys.exit(qlp_main())
