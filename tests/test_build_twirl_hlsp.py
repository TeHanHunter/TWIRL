from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "scripts" / "stage1_lightcurves"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import build_twirl_hlsp as module  # noqa: E402
from twirl.lightcurves.tglc_h5_reader import (  # noqa: E402
    APERTURE_KEYS,
    TGLCAperture,
    TGLCLightCurve,
)


def _light_curve(*, orbit: int, tic: int = 123, cam: int = 1, ccd: int = 1):
    values = np.array([1.0, 2.0]) + orbit
    apertures = {
        name: TGLCAperture(
            name=name,
            raw_flux=values.copy(),
            raw_flux_err=np.ones(2),
            raw_magnitude=np.ones(2),
            raw_magnitude_err=np.ones(2),
            centroid_x=np.ones(2),
            centroid_y=np.ones(2),
            flux_was_synthesized=False,
        )
        for name in APERTURE_KEYS
    }
    return TGLCLightCurve(
        tic=tic,
        sector=56,
        orbit=orbit,
        cam=cam,
        ccd=ccd,
        tmag=15.0,
        ra=1.0,
        dec=2.0,
        bjd_offset=2_457_000.0,
        time=np.array([orbit + 0.1, orbit + 0.2]),
        cadence=np.array([orbit * 10 + 1, orbit * 10 + 2]),
        quality=np.zeros(2, dtype=np.int32),
        centroid_x=np.ones(2),
        centroid_y=np.ones(2),
        background=np.ones(2),
        background_err=np.ones(2),
        apertures=apertures,
        path=None,
    )


def test_discover_orbit_hdf5_rejects_cross_detector_duplicate(tmp_path: Path):
    first = tmp_path / "ffi" / "cam1" / "ccd1" / "LC" / "123.h5"
    second = tmp_path / "ffi" / "cam2" / "ccd2" / "LC" / "123.h5"
    first.parent.mkdir(parents=True)
    second.parent.mkdir(parents=True)
    first.touch()
    second.touch()

    with pytest.raises(ValueError, match="detector selection must be explicit"):
        module.discover_orbit_hdf5(tmp_path)


@pytest.mark.parametrize(
    ("changed", "value"),
    [("tic", 456), ("sector", 57), ("cam", 2), ("ccd", 2)],
)
def test_merge_orbits_rejects_identity_or_detector_mismatch(changed, value):
    first = _light_curve(orbit=119)
    second = _light_curve(orbit=120)
    setattr(second, changed, value)

    with pytest.raises(ValueError, match=f"different {changed}"):
        module.merge_orbits([first, second])


def test_merge_orbits_preserves_detector_identity():
    merged = module.merge_orbits(
        [_light_curve(orbit=120), _light_curve(orbit=119)]
    )

    assert merged is not None
    assert (merged.tic, merged.sector, merged.cam, merged.ccd) == (123, 56, 1, 1)
    assert merged.orbitid.tolist() == [119, 119, 120, 120]
