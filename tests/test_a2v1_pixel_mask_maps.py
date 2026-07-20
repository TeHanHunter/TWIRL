from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

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
