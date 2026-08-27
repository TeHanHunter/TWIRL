from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from twirl.models.fm0.hdf5_quality_admission import (
    Hdf5QualityAuthority,
    check_hdf5_product,
)
from twirl.models.fm0.registry import FM0ContractError


def _authority(*, exclusions: frozenset[int] = frozenset()) -> Hdf5QualityAuthority:
    cadences = frozenset((100, 101))
    return Hdf5QualityAuthority(
        sector=67,
        quaternion={(141, 1): cadences},
        qflag={(141, 1, 1): {100: 0, 101: 1}},
        mission={(1, 1): {100: 4, 101: 0, 102: 0}},
        exclusions={(1, 1): exclusions},
    )


def _write_hdf5(path: Path, *, cadences: tuple[int, ...] = (100, 101)) -> None:
    with h5py.File(path, "w") as handle:
        handle.attrs["Sector"] = 67
        handle.attrs["Orbit"] = 141
        handle.attrs["Camera"] = 1
        handle.attrs["CCD"] = 1
        handle.attrs["TIC ID"] = 42
        group = handle.create_group("LightCurve")
        group.create_dataset("Cadence", data=np.asarray(cadences, dtype=np.int64))
        group.create_dataset(
            "QualityFlag", data=np.zeros(len(cadences), dtype=np.int64)
        )


def test_check_hdf5_product_joins_internal_and_external_quality(
    tmp_path: Path,
) -> None:
    path = tmp_path / "42.h5"
    _write_hdf5(path)
    result = check_hdf5_product(
        path=path,
        sector=67,
        orbit=141,
        camera=1,
        ccd=1,
        authority=_authority(),
    )
    assert result == {
        "n_cadences": 2,
        "n_internal_bad": 0,
        "n_external_bad": 2,
        "n_authority_excluded": 0,
    }


def test_check_hdf5_product_accepts_only_declared_authority_exclusion(
    tmp_path: Path,
) -> None:
    path = tmp_path / "42.h5"
    _write_hdf5(path, cadences=(100, 102))
    with pytest.raises(FM0ContractError, match="unexplained cadence 102"):
        check_hdf5_product(
            path=path,
            sector=67,
            orbit=141,
            camera=1,
            ccd=1,
            authority=_authority(),
        )
    result = check_hdf5_product(
        path=path,
        sector=67,
        orbit=141,
        camera=1,
        ccd=1,
        authority=_authority(exclusions=frozenset((102,))),
    )
    assert result["n_authority_excluded"] == 1
    assert result["n_external_bad"] == 2


def test_check_hdf5_product_rejects_identity_and_array_drift(tmp_path: Path) -> None:
    path = tmp_path / "42.h5"
    _write_hdf5(path)
    with h5py.File(path, "r+") as handle:
        handle.attrs["CCD"] = 2
    with pytest.raises(FM0ContractError, match="wrong CCD"):
        check_hdf5_product(
            path=path,
            sector=67,
            orbit=141,
            camera=1,
            ccd=1,
            authority=_authority(),
        )
