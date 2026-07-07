from __future__ import annotations

import numpy as np
from astropy.io import fits
from astropy.table import Table

from scripts.stage1_lightcurves.build_source_tic_overlays import (
    load_observation_targets,
    normalize_tic_catalog_with_positions,
    overlay_for_source_bounds,
    parse_source_grid,
    source_bounds,
)


def test_parse_source_grid_and_bounds_match_tglc_stride() -> None:
    assert parse_source_grid(type("P", (), {"name": "source_7_4.pkl"})()) == (7, 4)
    assert source_bounds(7, 4) == (1066.0, 1216.0, 584.0, 734.0)


def test_normalize_tic_catalog_with_positions_filters_requested_tics() -> None:
    table = Table()
    table["id"] = [10, 11, 12, 13]
    table["gaia3"] = [100, 101, 102, 103]
    table["tmag"] = [19.0, 20.1, 18.0, 17.0]
    table["ra"] = [1.0, 2.0, np.nan, 4.0]
    table["dec"] = [5.0, 6.0, 7.0, 8.0]

    got = normalize_tic_catalog_with_positions(
        table,
        tmag_max=99,
        allowed_tics={10, 11, 12},
    )

    assert list(got["TIC"]) == [10, 11]
    assert list(got["gaia3"]) == [100, 101]
    assert list(got["ra"]) == [1.0, 2.0]


def test_overlay_for_source_bounds_writes_tglc_source_tic_shape() -> None:
    projected = Table()
    projected["TIC"] = [201, 202, 203, 204]
    projected["gaia3"] = [1001, 1002, 1003, 1004]
    projected["ccd_x"] = [10.0, 20.0, 30.0, 40.0]
    projected["ccd_y"] = [10.0, 25.0, 35.0, 50.0]

    got = overlay_for_source_bounds(projected, (15.0, 35.0, 20.0, 40.0))

    assert got.colnames == ["TIC", "gaia3"]
    assert list(got["TIC"]) == [202, 203]
    assert list(got["gaia3"]) == [1002, 1003]


def test_load_observation_targets_uses_detector_coordinates(tmp_path) -> None:
    path = tmp_path / "observations.fits"
    rows = Table()
    rows["source_id"] = [1001, 1002, 1003, 1002]
    rows["tic_id"] = [201, 202, 203, 202]
    rows["orbit"] = [119, 119, 120, 119]
    rows["camera"] = [4, 4, 4, 4]
    rows["ccd"] = [4, 4, 4, 4]
    rows["colpix"] = [10.0, 20.0, 30.0, 21.0]
    rows["rowpix"] = [11.0, 22.0, 33.0, 23.0]
    rows["tmag"] = [19.0, 20.5, 18.0, 20.5]
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(rows)]).writeto(path)

    got = load_observation_targets(path, orbit=119, camera=4, ccd=4, tmag_max=99)

    assert got.colnames == ["TIC", "gaia3", "ccd_x", "ccd_y"]
    assert list(got["TIC"]) == [201, 202]
    assert list(got["gaia3"]) == [1001, 1002]
    assert list(got["ccd_x"]) == [10.0, 20.0]
