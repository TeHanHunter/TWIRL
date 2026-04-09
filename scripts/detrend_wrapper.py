#!/usr/bin/env python3
"""Wrapper for qlp lctools detrend that patches orbit info to bypass DB lookup.

lightcurvedb 3.0.0 requires tess_orbit table not populated in stelvar.
This wrapper patches both call sites so detrend runs without a live database.

Usage on PDO (from /pdo/users/tehan/TWIRL/):
    LD_LIBRARY_PATH=/pdo/app/anaconda/anaconda2-4.4.0/lib:/pdo/app/python-versions/python-3.11.9/lib \\
    PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl \\
    /pdo/app/qlp-environment/.venv/bin/python detrend_wrapper.py \\
        lctools detrend \\
        -c /pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/run/qlp.cfg \\
        --autolist 1:1,2,3,4 2:1,2,3,4 3:1,2,3,4 4:1,2,3,4 \\
        --best -n 8 \\
        --logfile /pdo/users/tehan/TWIRL/logs/detrend_orbit119.log
"""

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

from qlp.__main__ import qlp_main
sys.exit(qlp_main())
