from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "stage1_lightcurves" / "plot_a2v1_pixel_mask_maps.py"
SPEC = importlib.util.spec_from_file_location("a2v1_pixel_mask_maps", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_central_epsf_index_uses_the_center_of_the_psf_block() -> None:
    assert MODULE.central_epsf_index(535) == 264
    assert MODULE.central_epsf_index(447) == 220


def test_central_epsf_index_rejects_non_square_or_even_psf_blocks() -> None:
    with pytest.raises(ValueError, match="odd square block"):
        MODULE.central_epsf_index(536)
    with pytest.raises(ValueError, match="odd square block"):
        MODULE.central_epsf_index(10)


def test_parse_sector_orbit() -> None:
    assert MODULE.parse_sector_orbit("64:135") == (64, 135)
    with pytest.raises(Exception, match="sector orbit"):
        MODULE.parse_sector_orbit("64")


def test_science_pixel_slices_match_the_tglc_source_coordinate_frame() -> None:
    y_slice, x_slice, extent = MODULE.science_pixel_slices("[45:2092,1:2048]")

    assert (y_slice.start, y_slice.stop) == (0, 2048)
    assert (x_slice.start, x_slice.stop) == (44, 2092)
    assert extent == (44.0, 2092.0, 0.0, 2048.0)


def test_bad_pixel_proxy_expands_a_saturated_pixel_to_cross_neighbors() -> None:
    image = np.ones((5, 5), dtype=float)
    image[2, 2] = 100.0

    mask = MODULE.bad_pixel_proxy(image)

    assert int(mask.sum()) == 5
    assert mask[2, 2]
    assert mask[1, 2]
    assert mask[3, 2]
    assert mask[2, 1]
    assert mask[2, 3]
