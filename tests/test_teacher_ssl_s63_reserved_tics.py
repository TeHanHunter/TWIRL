from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
from astropy.table import Table
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts"
    / "stage5_validation"
    / "build_teacher_ssl_s63_reserved_tics.py"
)
SPEC = importlib.util.spec_from_file_location("reserved_s63_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_reserved_s63_inventory_is_identity_only_and_deterministic(
    tmp_path: Path,
) -> None:
    observations = tmp_path / "observations.fits"
    Table(
        {
            "sector": np.asarray([63, 63, 63, 62, 63]),
            "orbit": np.asarray([134, 133, 133, 132, 135]),
            "tic_id": np.asarray([9, 7, 9, 5, 11]),
        }
    ).write(observations)
    output = tmp_path / "s63_reserved_tics.txt"

    summary = MODULE.build_reserved_tics(
        observations=observations,
        out_tics=output,
    )

    assert output.read_text() == "7\n9\n"
    assert summary["n_selected_observation_rows"] == 3
    assert summary["n_reserved_tics"] == 2
    assert summary["identity_only"] is True
    assert summary["light_curves_opened"] is False
    on_disk = json.loads(output.with_suffix(".summary.json").read_text())
    assert on_disk["reserved_tics_sha256"] == summary["reserved_tics_sha256"]


def test_reserved_s63_inventory_rejects_scope_drift(tmp_path: Path) -> None:
    observations = tmp_path / "observations.fits"
    Table(
        {
            "sector": [63],
            "orbit": [133],
            "tic_id": [7],
        }
    ).write(observations)

    with pytest.raises(ValueError, match="bounded to S63"):
        MODULE.build_reserved_tics(
            observations=observations,
            out_tics=tmp_path / "bad.txt",
            sector=62,
        )
    with pytest.raises(ValueError, match="exact orbits"):
        MODULE.build_reserved_tics(
            observations=observations,
            out_tics=tmp_path / "bad.txt",
            orbits=(133,),
        )
