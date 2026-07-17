from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from twirl.lightcurves.tglc_faint_allstar import (
    enforce_pdo_user_output,
    normalize_tic_catalog,
    source_target_table,
)


def test_enforce_pdo_user_output_rejects_shared_tree() -> None:
    with pytest.raises(ValueError):
        enforce_pdo_user_output(Path("/pdo/qlp-data/orbit-195/ffi"))

    got = enforce_pdo_user_output(Path("/pdo/users/tehan/tglc-s94-faint-allstar"))
    assert str(got).startswith("/pdo/users/tehan/")


def test_normalize_tic_catalog_filters_tmag_and_bad_gaia() -> None:
    table = Table()
    table["id"] = [101, 102, 103, 104]
    table["gaia3"] = [1001, 1002, None, 1004]
    table["tmag"] = [17.9, 18.5, 18.2, 18.6]

    got = normalize_tic_catalog(table, tmag_max=18.5)

    assert list(got["TIC"]) == [101, 102]
    assert list(got["gaia3"]) == [1001, 1002]
    assert np.all(np.asarray(got["tmag"]) <= 18.5)


def test_source_target_table_matches_only_gaia_ids_inside_source() -> None:
    source_gaia = Table()
    source_gaia["designation"] = [
        "Gaia DR3 1001",
        "Gaia DR3 1003",
        "Gaia DR3 9999",
    ]

    tic_catalog = Table()
    tic_catalog["id"] = [201, 202, 203, 204]
    tic_catalog["gaia3"] = [1001, 1002, 1003, 1004]
    tic_catalog["tmag"] = [17.5, 17.0, 18.4, 18.7]

    got = source_target_table(source_gaia, tic_catalog, tmag_max=18.5)

    assert got.colnames == ["TIC", "gaia3"]
    assert list(got["TIC"]) == [201, 203]
    assert list(got["gaia3"]) == [1001, 1003]
