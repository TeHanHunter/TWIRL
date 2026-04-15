#!/usr/bin/env python3
"""Wrapper for qlp lctools hlsp that patches sector/orbit lookup to bypass DB.

lightcurvedb 3.0.0's hlsp stage calls `read_sector_info(sector)`, which queries
the `tess_orbit` table (missing on PDO stelvar). This wrapper replaces the
postgres/default backend's `read_sector_info` with a hard-coded lookup for the
TWIRL benchmark sectors.

Companion to `detrend_wrapper.py`, which patches the orbit-level lookup used by
`qlp lctools detrend`.

Usage (via activate_qlp_env.sh env vars):
    source /pdo/users/tehan/TWIRL/scripts/activate_qlp_env.sh
    "$TWIRL_QLP_PYTHON" /pdo/users/tehan/TWIRL/scripts/hlsp_wrapper.py \
        lctools hlsp \
        -c /pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/run/qlp.cfg \
        -s 56 --autolist all \
        --flag-type spoc --flag-source fits \
        -o /pdo/users/tehan/tglc-deep-catalogs/hlsp_s0056 \
        -n 8 -r
"""

# Thread caps MUST be set before numpy/scipy are imported. Without them, each
# worker spawns BLAS threads equal to core count, so `-n 8` becomes ~1000
# threads on a 128-core box and everything thrashes (see detrend_wrapper.py).
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

from qlp.io.datastructures import OrbitalInfo, SectorInfo


# Hard-coded TWIRL benchmark sector-orbit metadata. Cadence ranges come from
# the TICA FFI filenames under /pdo/qlp-data/tica-delivery/s00XX/camC-ccdD/.
# mjd_range is a placeholder (hlsp only iterates orbit_number at the call sites
# we hit). Extend this table when adding new benchmark sectors.
SECTOR_ORBIT_TABLE: dict[int, list[OrbitalInfo]] = {
    56: [
        OrbitalInfo(
            orbit_number=119,
            sector=56,
            cadence_range=(690247, 695975),
            mjd_range=(0.0, 0.0),
            basename="",
        ),
        OrbitalInfo(
            orbit_number=120,
            sector=56,
            cadence_range=(696067, 702293),
            mjd_range=(0.0, 0.0),
            basename="",
        ),
    ],
}

ORBIT_SECTOR_MAP = {
    orbit.orbit_number: sector
    for sector, orbits in SECTOR_ORBIT_TABLE.items()
    for orbit in orbits
}


def _mock_read_sector_info(sector: int) -> SectorInfo:
    orbits = SECTOR_ORBIT_TABLE.get(int(sector))
    if not orbits:
        raise ValueError(
            f"hlsp_wrapper: sector {sector} not in SECTOR_ORBIT_TABLE; "
            f"add its orbit list before running."
        )
    min_cadence = orbits[0].cadence_range[0]
    max_cadence = orbits[-1].cadence_range[1]
    min_mid_tjd = orbits[0].mjd_range[0]
    max_mid_tjd = orbits[-1].mjd_range[1]
    return SectorInfo(
        sector=int(sector),
        orbit_info=orbits,
        cadence_range=(min_cadence, max_cadence),
        mjd_range=(min_mid_tjd, max_mid_tjd),
    )


def _mock_read_orbit_info(orbit_number: int) -> OrbitalInfo:
    sector = ORBIT_SECTOR_MAP.get(int(orbit_number))
    if sector is None:
        return OrbitalInfo(
            orbit_number=int(orbit_number),
            sector=56,
            cadence_range=(0, 0),
            mjd_range=(0.0, 0.0),
            basename="",
        )
    for orbit in SECTOR_ORBIT_TABLE[sector]:
        if orbit.orbit_number == int(orbit_number):
            return orbit
    raise RuntimeError("unreachable")


# Patch the default IO backend used by hlsp and friends.
import qlp.io.backends

qlp.io.backends.default_io_backend.read_sector_info = _mock_read_sector_info
qlp.io.backends.default_io_backend.read_orbit_info = _mock_read_orbit_info

# Also patch util helper in case any other code path calls it.
import qlp.util.util

qlp.util.util._get_orbit_info = _mock_read_orbit_info

# Make the per-target HLSP worker fault-tolerant. Some h5s lack `bestdmagkey`
# (detrend failed on them due to an upstream cadence mismatch) or have other
# per-target defects. Without this, one bad target aborts the entire pool.
import qlp.lctools.bin.hlsp as _qlp_hlsp
_orig_hlsp_worker = _qlp_hlsp.generate_qlp_hlsp_fits_file_wrapper


def _safe_hlsp_worker(tic_id_camera_ccd_shared_flags, *args, **kwargs):
    try:
        return _orig_hlsp_worker(tic_id_camera_ccd_shared_flags, *args, **kwargs)
    except Exception as e:
        tic = cam = ccd = "?"
        try:
            tic, cam, ccd, *_ = tic_id_camera_ccd_shared_flags
        except Exception:
            pass
        print(
            f"[hlsp skip] tic={tic} cam={cam} ccd={ccd}: "
            f"{type(e).__name__}: {e}",
            flush=True,
        )
        return False


_qlp_hlsp.generate_qlp_hlsp_fits_file_wrapper = _safe_hlsp_worker

from qlp.__main__ import qlp_main

sys.exit(qlp_main())
