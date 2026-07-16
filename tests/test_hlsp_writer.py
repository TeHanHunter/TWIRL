from pathlib import Path
import importlib.util
import json
import sys

import h5py
import numpy as np
import pytest
from astropy.table import Table
from astropy.io import fits

from twirl.lightcurves.flux_detrend import FluxDetrendConfig
from twirl.lightcurves.hlsp_writer import HLSPTarget, hlsp_path, write_twirl_hlsp


def test_write_twirl_hlsp_can_emit_a2v1_only_flux_columns(tmp_path: Path) -> None:
    n = 5
    target = HLSPTarget(
        tic=123456789,
        sector=56,
        cam=4,
        ccd=1,
        tmag=18.5,
        ra=1.0,
        dec=2.0,
    )
    out_path = tmp_path / "a2v1.fits"
    base = np.ones(n, dtype=np.float32)
    branch_columns = {
        "DET_FLUX_ADP": base,
        "DET_FLUX_ADP_ERR": base * 0.01,
        "DET_FLUX_ADP_SML": base + 0.1,
        "DET_FLUX_ADP_LAG": base + 0.2,
        "DET_FLUX_ADP015": base + 0.3,
        "DET_FLUX_ADP015_ERR": base * 0.02,
        "DET_FLUX_ADP015_SML": base + 0.4,
        "DET_FLUX_ADP015_LAG": base + 0.5,
    }

    write_twirl_hlsp(
        out_path,
        target,
        time_btjd=np.arange(n, dtype=float),
        cadenceno=np.arange(n, dtype=np.int32),
        sap_flux=base,
        det_flux=base,
        det_flux_err=base,
        quality=np.zeros(n, dtype=np.int32),
        orbitid=np.full(n, 119, dtype=np.int32),
        sap_x=base,
        sap_y=base,
        sap_bkg=base,
        sap_bkg_err=base,
        det_flux_sml=base,
        det_flux_lag=base,
        detrend_config=FluxDetrendConfig(),
        method_version="A2v1",
        extra_flux_columns=branch_columns,
        extra_header={"PRODTAG": ("A2v1", "TWIRL production product tag")},
        include_canonical_det_flux=False,
        include_sys_rm_flux=False,
    )

    with fits.open(out_path, memmap=False) as hdul:
        columns = set(hdul["LIGHTCURVE"].columns.names)
        assert hdul[0].header["METHOD"] == "A2v1"
        assert hdul[0].header["PRODTAG"] == "A2v1"

    assert "DET_FLUX" not in columns
    assert "DET_FLUX_ERR" not in columns
    assert "DET_FLUX_SML" not in columns
    assert "DET_FLUX_LAG" not in columns
    assert "SYS_RM_FLUX" not in columns
    assert set(branch_columns) <= columns


def _load_a2v1_validator():
    script = Path(__file__).parents[1] / "scripts" / "stage1_lightcurves" / "validate_a2v1_product.py"
    spec = importlib.util.spec_from_file_location("validate_a2v1_product", script)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_observation_table(
    path: Path,
    *,
    edge_warn: bool = False,
    orbits: tuple[int, ...] = (119, 120),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Table(
        {
            "sector": [56] * len(orbits),
            "orbit": list(orbits),
            "camera": [4] * len(orbits),
            "ccd": [1] * len(orbits),
            "tic_id": [123456789] * len(orbits),
            "edge_warn": [edge_warn] * len(orbits),
        }
    ).write(path)


def _write_placeholder_h5(root: Path, orbit: int) -> None:
    h5_path = root / f"orbit-{orbit}" / "ffi" / "cam4" / "ccd1" / "LC" / "123456789.h5"
    h5_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(h5_path, "w") as handle:
        handle.create_group("LightCurve")


def _write_a2v1_fits(root: Path, *, include_canonical: bool) -> None:
    n = 5
    target = HLSPTarget(
        tic=123456789,
        sector=56,
        cam=4,
        ccd=1,
        tmag=18.5,
        ra=1.0,
        dec=2.0,
    )
    base = np.ones(n, dtype=np.float32)
    branch_columns = {
        "DET_FLUX_ADP": base,
        "DET_FLUX_ADP_ERR": base * 0.01,
        "DET_FLUX_ADP_SML": base + 0.1,
        "DET_FLUX_ADP_LAG": base + 0.2,
        "DET_FLUX_ADP015": base + 0.3,
        "DET_FLUX_ADP015_ERR": base * 0.02,
        "DET_FLUX_ADP015_SML": base + 0.4,
        "DET_FLUX_ADP015_LAG": base + 0.5,
    }
    write_twirl_hlsp(
        hlsp_path(root, target),
        target,
        time_btjd=np.arange(n, dtype=float),
        cadenceno=np.arange(n, dtype=np.int32),
        sap_flux=base,
        det_flux=base,
        det_flux_err=base,
        quality=np.zeros(n, dtype=np.int32),
        orbitid=np.array([119, 119, 119, 120, 120], dtype=np.int32),
        sap_x=base,
        sap_y=base,
        sap_bkg=base,
        sap_bkg_err=base,
        det_flux_sml=base,
        det_flux_lag=base,
        detrend_config=FluxDetrendConfig(),
        method_version="A2v1",
        extra_flux_columns=branch_columns,
        extra_header={
            "PRODTAG": ("A2v1", "TWIRL production product tag"),
            "A2V1": (True, "A2v1 ADP/ADP015-only branch product"),
            "BRANCHES": ("ADP,ADP015", "A2v1 branch families"),
        },
        include_canonical_det_flux=include_canonical,
        include_sys_rm_flux=include_canonical,
    )


def test_validate_a2v1_product_accepts_adp_adp015_only_schema(tmp_path: Path) -> None:
    module = _load_a2v1_validator()
    root = tmp_path / "a2v1"
    hlsp_root = root / "hlsp_s0056_A2v1"
    observations = tmp_path / "observations.fits"
    summary_json = tmp_path / "summary.json"
    _write_observation_table(observations)
    _write_placeholder_h5(root, 119)
    _write_placeholder_h5(root, 120)
    _write_a2v1_fits(hlsp_root, include_canonical=False)

    rc = module.main(
        [
            "--a2v1-root", str(root),
            "--hlsp-root", str(hlsp_root),
            "--observations", str(observations),
            "--sector", "56",
            "--orbits", "119", "120",
            "--check-h5-open",
            "--summary-json", str(summary_json),
        ]
    )

    summary = json.loads(summary_json.read_text())
    assert rc == 0
    assert summary["ok"] is True
    assert summary["h5"]["n_missing_h5"] == 0
    assert summary["h5"]["n_unreadable_h5"] == 0
    assert summary["fits"]["n_bad_checked_fits"] == 0


def test_validate_a2v1_product_rejects_unreadable_h5(tmp_path: Path) -> None:
    module = _load_a2v1_validator()
    root = tmp_path / "a2v1"
    observations = tmp_path / "observations.fits"
    summary_json = tmp_path / "summary.json"
    _write_observation_table(observations)
    _write_placeholder_h5(root, 119)
    _write_placeholder_h5(root, 120)
    bad_path = root / "orbit-120" / "ffi" / "cam4" / "ccd1" / "LC" / "123456789.h5"
    bad_path.write_bytes(b"not an HDF5 file")

    rc = module.main(
        [
            "--a2v1-root", str(root),
            "--observations", str(observations),
            "--sector", "56",
            "--orbits", "119", "120",
            "--skip-fits",
            "--check-h5-open",
            "--summary-json", str(summary_json),
        ]
    )

    summary = json.loads(summary_json.read_text())
    assert rc == 1
    assert summary["ok"] is False
    assert summary["h5"]["n_unreadable_h5"] == 1
    assert summary["h5"]["unreadable_h5_examples"][0]["path"] == str(bad_path)
    assert summary["fits"]["skipped"] is True


def test_validate_a2v1_product_rejects_canonical_det_flux_columns(tmp_path: Path) -> None:
    module = _load_a2v1_validator()
    root = tmp_path / "a2v1"
    hlsp_root = root / "hlsp_s0056_A2v1"
    observations = tmp_path / "observations.fits"
    summary_json = tmp_path / "summary.json"
    _write_observation_table(observations)
    _write_placeholder_h5(root, 119)
    _write_placeholder_h5(root, 120)
    _write_a2v1_fits(hlsp_root, include_canonical=True)

    rc = module.main(
        [
            "--a2v1-root", str(root),
            "--hlsp-root", str(hlsp_root),
            "--observations", str(observations),
            "--sector", "56",
            "--orbits", "119", "120",
            "--allow-missing-fits",
            "--summary-json", str(summary_json),
        ]
    )

    summary = json.loads(summary_json.read_text())
    assert rc == 1
    assert summary["ok"] is False
    assert summary["fits"]["n_bad_checked_fits"] == 1
    bad = summary["fits"]["bad_fits_examples"][0]
    assert "DET_FLUX" in bad["forbidden_columns_present"]


def test_validate_a2v1_product_rejects_empty_expected_contract(tmp_path: Path) -> None:
    module = _load_a2v1_validator()
    root = tmp_path / "missing-a2v1"
    observations = tmp_path / "observations.fits"
    summary_json = tmp_path / "summary.json"
    _write_observation_table(observations, orbits=())

    rc = module.main(
        [
            "--a2v1-root", str(root),
            "--observations", str(observations),
            "--sector", "56",
            "--orbits", "119", "120",
            "--skip-fits",
            "--allow-missing-h5",
            "--summary-json", str(summary_json),
        ]
    )

    summary = json.loads(summary_json.read_text())
    assert rc == 1
    assert summary["ok"] is False
    assert summary["ok_h5"] is False
    assert summary["n_expected_h5"] == 0
    assert summary["expected_contract"]["ok"] is False
    assert summary["expected_contract"]["missing_requested_orbits"] == [119, 120]


def test_validate_a2v1_product_rejects_missing_requested_orbit_rows(tmp_path: Path) -> None:
    module = _load_a2v1_validator()
    root = tmp_path / "a2v1"
    observations = tmp_path / "observations.fits"
    summary_json = tmp_path / "summary.json"
    _write_observation_table(observations, orbits=(119,))
    _write_placeholder_h5(root, 119)

    rc = module.main(
        [
            "--a2v1-root", str(root),
            "--observations", str(observations),
            "--sector", "56",
            "--orbits", "119", "120",
            "--skip-fits",
            "--summary-json", str(summary_json),
        ]
    )

    summary = json.loads(summary_json.read_text())
    assert rc == 1
    assert summary["expected_contract"]["observed_orbits"] == [119]
    assert summary["expected_contract"]["missing_requested_orbits"] == [120]


def test_validate_a2v1_product_allows_missing_fits_only(tmp_path: Path) -> None:
    module = _load_a2v1_validator()
    root = tmp_path / "a2v1"
    observations = tmp_path / "observations.fits"
    summary_json = tmp_path / "summary.json"
    _write_observation_table(observations)
    _write_placeholder_h5(root, 119)
    _write_placeholder_h5(root, 120)

    rc = module.main(
        [
            "--a2v1-root", str(root),
            "--observations", str(observations),
            "--sector", "56",
            "--orbits", "119", "120",
            "--allow-missing-fits",
            "--summary-json", str(summary_json),
        ]
    )

    summary = json.loads(summary_json.read_text())
    assert rc == 0
    assert summary["ok"] is True
    assert summary["fits"]["n_missing_fits_tics"] == 1
    assert summary["fits"]["n_bad_checked_fits"] == 0


def test_validate_a2v1_product_can_allow_documented_edge_exclusions(tmp_path: Path) -> None:
    module = _load_a2v1_validator()
    root = tmp_path / "a2v1"
    hlsp_root = root / "hlsp_s0056_A2v1"
    observations = tmp_path / "observations.fits"
    summary_json = tmp_path / "summary.json"
    _write_observation_table(observations, edge_warn=True)

    rc = module.main(
        [
            "--a2v1-root", str(root),
            "--hlsp-root", str(hlsp_root),
            "--observations", str(observations),
            "--sector", "56",
            "--orbits", "119", "120",
            "--allow-edge-warn-missing",
            "--summary-json", str(summary_json),
        ]
    )

    summary = json.loads(summary_json.read_text())
    assert rc == 0
    assert summary["ok"] is True
    assert summary["h5"]["n_missing_h5_edge_warn"] == 2
    assert summary["h5"]["n_missing_h5_non_edge"] == 0
    assert summary["fits"]["n_missing_fits_edge_warn_tics"] == 1
    assert summary["fits"]["n_missing_fits_non_edge_tics"] == 0


def test_write_twirl_hlsp_rejects_misaligned_core_array_before_writing(
    tmp_path: Path,
) -> None:
    n = 5
    target = HLSPTarget(
        tic=123456789,
        sector=56,
        cam=4,
        ccd=1,
        tmag=18.5,
        ra=1.0,
        dec=2.0,
    )
    out_path = tmp_path / "nested" / "misaligned.fits"
    base = np.ones(n, dtype=np.float32)

    with pytest.raises(ValueError, match="'quality' has length 4; expected 5"):
        write_twirl_hlsp(
            out_path,
            target,
            time_btjd=np.arange(n, dtype=float),
            cadenceno=np.arange(n, dtype=np.int32),
            sap_flux=base,
            det_flux=base,
            det_flux_err=base,
            quality=np.zeros(n - 1, dtype=np.int32),
            orbitid=np.full(n, 119, dtype=np.int32),
            sap_x=base,
            sap_y=base,
            sap_bkg=base,
            sap_bkg_err=base,
            det_flux_sml=base,
            det_flux_lag=base,
        )

    assert not out_path.exists()
    assert not out_path.parent.exists()
