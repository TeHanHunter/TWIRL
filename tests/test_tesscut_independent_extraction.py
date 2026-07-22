from __future__ import annotations

import json
from pathlib import Path
import stat
import subprocess
import sys

from astropy.io import fits
import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.lightcurves.tesscut_independent_extraction as tesscut_module
from twirl.lightcurves.a2v1_cadence_reference import CADENCE_REFERENCE_COLUMNS
from twirl.lightcurves.a2v1_qa import WD1856_TIC, file_sha256
from twirl.lightcurves.tesscut_independent_extraction import (
    BUILDER_VERSION,
    WD1856_TESSCUT_DEC_DEG,
    WD1856_TESSCUT_RA_DEG,
    build_wd1856_tesscut_independent_extraction,
    centered_rolling_median_detrend,
)


SOURCE_URL = (
    "https://mast.stsci.edu/api/v0.1/Download/file/?uri="
    "mast:TESS/product/wd1856-s56-tesscut.fits"
)


def _write_tesscut(
    path: Path,
    *,
    n_rows: int = 9,
    size: int = 5,
    sector: int = 56,
    camera: int = 4,
    ccd: int = 1,
    ra_deg: float = WD1856_TESSCUT_RA_DEG,
    dec_deg: float = WD1856_TESSCUT_DEC_DEG,
    include_wcs: bool = True,
    nonfinite_flux: bool = False,
    duplicate_time: bool = False,
) -> dict[str, np.ndarray]:
    time = 2_500.0 + np.arange(n_rows, dtype=float) * 200.0 / 86_400.0
    if duplicate_time:
        time[1] = time[0]
    quality = np.zeros(n_rows, dtype=np.int32)
    quality[n_rows // 2] = 8
    cube = np.empty((n_rows, size, size), dtype=np.float32)
    for index in range(n_rows):
        cube[index] = 2.0 + 0.1 * index
        cube[index, size // 2, size // 2] = 100.0 + 2.5 * index
    if nonfinite_flux:
        cube[0, 0, 0] = np.nan

    columns = [
        fits.Column(name="TIME", format="D", array=time),
        fits.Column(name="QUALITY", format="J", array=quality),
        fits.Column(
            name="FLUX",
            format=f"{size * size}E",
            dim=f"({size},{size})",
            array=cube,
        ),
    ]
    pixels = fits.BinTableHDU.from_columns(columns, name="PIXELS")
    pixels.header["BJDREFI"] = 2457000
    pixels.header["BJDREFF"] = 0.0
    pixels.header["TIMEUNIT"] = "d"
    pixels.header["TIMESYS"] = "TDB"
    if include_wcs:
        flux_column = 3
        pixels.header[f"1CTYP{flux_column}"] = "RA---TAN"
        pixels.header[f"2CTYP{flux_column}"] = "DEC--TAN"
        pixels.header[f"1CRPX{flux_column}"] = float(size // 2 + 1)
        pixels.header[f"2CRPX{flux_column}"] = float(size // 2 + 1)
        pixels.header[f"1CRVL{flux_column}"] = WD1856_TESSCUT_RA_DEG
        pixels.header[f"2CRVL{flux_column}"] = WD1856_TESSCUT_DEC_DEG
        pixels.header[f"1CUNI{flux_column}"] = "deg"
        pixels.header[f"2CUNI{flux_column}"] = "deg"
        pixels.header[f"1CDLT{flux_column}"] = -0.006
        pixels.header[f"2CDLT{flux_column}"] = 0.006

    primary = fits.PrimaryHDU()
    primary.header["ORIGIN"] = "STScI/MAST"
    primary.header["CREATOR"] = "astrocut"
    primary.header["PROCVER"] = "test-astrocut-1.0"
    primary.header["SECTOR"] = sector
    primary.header["CAMERA"] = camera
    primary.header["CCD"] = ccd
    primary.header["RA_OBJ"] = ra_deg
    primary.header["DEC_OBJ"] = dec_deg
    fits.HDUList([primary, pixels]).writeto(path, checksum=True)
    small_raw = cube[:, size // 2, size // 2].astype(np.float64)
    primary_raw = np.sum(
        cube[
            :,
            size // 2 - 1 : size // 2 + 2,
            size // 2 - 1 : size // 2 + 2,
        ],
        axis=(1, 2),
        dtype=np.float64,
    )
    return {
        "time": time,
        "quality": quality,
        "small_raw": small_raw,
        "primary_raw": primary_raw,
    }


def _write_compact(
    path: Path,
    time: np.ndarray,
    *,
    offset_seconds: float = 17.0,
    duplicate_time: bool = False,
) -> np.ndarray:
    compact_time = np.asarray(time, dtype=float) + offset_seconds / 86_400.0
    if duplicate_time:
        compact_time[1] = compact_time[0]
    cadence = np.arange(1_250_000, 1_250_000 + len(time), dtype=np.int64)
    quality = np.zeros(len(time), dtype=np.int32)
    quality[0] = 4
    orbitid = np.where(np.arange(len(time)) < len(time) // 2, 119, 120)
    with h5py.File(path, "w") as h5:
        h5.attrs["sector"] = 56
        h5.attrs["time_unit"] = "BJD - 2457000"
        targets = h5.create_group("targets")
        group = targets.create_group(f"{WD1856_TIC:016d}")
        group.attrs["tic"] = WD1856_TIC
        group.attrs["sector"] = 56
        group.attrs["camera"] = 4
        group.attrs["ccd"] = 1
        group.create_dataset("time", data=compact_time)
        group.create_dataset("cadenceno", data=cadence)
        group.create_dataset("quality", data=quality)
        group.create_dataset("orbitid", data=orbitid)
    return cadence


def _write_quality_reference(root: Path, cadence: np.ndarray) -> tuple[Path, Path]:
    table_path = root / "cadence_reference.csv"
    manifest_path = root / "cadence_reference.json"
    orbitid = np.where(np.arange(len(cadence)) < len(cadence) // 2, 119, 120)
    spoc_quality = np.zeros(len(cadence), dtype=np.int64)
    qlp_quality = np.zeros(len(cadence), dtype=np.int64)
    spoc_quality[1] = 16
    qlp_quality[-2] = 1
    frame = pd.DataFrame(
        {
            "sector": 56,
            "orbitid": orbitid,
            "camera": 4,
            "ccd": 1,
            "cadenceno": cadence,
            "spoc_quality": spoc_quality,
            "qlp_quality": qlp_quality,
            "external_quality": spoc_quality | (qlp_quality << 30),
        },
        columns=CADENCE_REFERENCE_COLUMNS,
    )
    frame.to_csv(table_path, index=False)
    sources: list[dict[str, object]] = []

    def add_source(role: str, index: int, **metadata: int) -> None:
        sources.append(
            {
                "role": role,
                "path": str((root / f"quality_source_{index:02d}.txt").resolve()),
                "sha256": f"{index:064x}",
                **metadata,
            }
        )

    add_source("qlp_cam_quat", 1, orbitid=119, camera=4)
    add_source("qlp_cam_quat", 2, orbitid=120, camera=4)
    add_source("qlp_detector_qflag", 3, orbitid=119, camera=4, ccd=1)
    add_source("qlp_detector_qflag", 4, orbitid=120, camera=4, ccd=1)
    add_source("spoc_flag_file", 5, camera=4, ccd=1)
    add_source("spoc_quality_table", 6)
    add_source("spoc_quality_provenance", 7)
    manifest = {
        "contract_version": "s56_a2v1_cadence_reference_v1",
        "sector": 56,
        "cadence_authority": "qlp_cam_quat",
        "quality_authority": "spoc_and_qlp_quality_flags",
        "quality_composition": {
            "external_quality": "spoc_quality | (qlp_quality << 30)",
            "qlp_quality_raw_values": [0, 1],
            "qlp_quality_external_bit": 30,
        },
        "table_sha256": file_sha256(table_path),
        "table_columns": list(CADENCE_REFERENCE_COLUMNS),
        "n_rows": len(frame),
        "detectors": ["cam4_ccd1"],
        "orbits": [119, 120],
        "n_rows_by_detector": {"cam4_ccd1": len(frame)},
        "n_nonzero_spoc_quality": 1,
        "n_nonzero_qlp_quality": 1,
        "n_nonzero_external_quality": 2,
        "n_spoc_authority_files_verified": 1,
        "n_qlp_qflag_files_verified": 2,
        "source_file_sha256": {
            str(source["path"]): str(source["sha256"]) for source in sources
        },
        "sources": sources,
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return table_path, manifest_path


def _inputs(tmp_path: Path, **tesscut_options: object) -> tuple[Path, Path, dict, np.ndarray]:
    tesscut = tmp_path / "tesscut.fits"
    raw = _write_tesscut(tesscut, **tesscut_options)
    compact = tmp_path / "compact.h5"
    cadence = _write_compact(compact, raw["time"])
    _write_quality_reference(tmp_path, cadence)
    return tesscut, compact, raw, cadence


def _build(tmp_path: Path, tesscut: Path, compact: Path, **options: object) -> dict:
    return build_wd1856_tesscut_independent_extraction(
        tesscut_fits=tesscut,
        compact_lc=compact,
        cadence_reference_table=tmp_path / "cadence_reference.csv",
        cadence_reference_manifest=tmp_path / "cadence_reference.json",
        output_fits=tmp_path / "reference.fits",
        manifest_json=tmp_path / "reference.json",
        source_url=SOURCE_URL,
        rolling_window_cadences=5,
        **options,
    )


def test_build_extracts_centered_apertures_and_detrends_independently(
    tmp_path: Path,
) -> None:
    tesscut, compact, raw, cadence = _inputs(tmp_path)
    effective_quality = (raw["quality"] != 0).astype(np.int32)
    effective_quality[[1, len(effective_quality) - 2]] = 1
    expected_small, expected_small_trend, minimum = centered_rolling_median_detrend(
        raw["small_raw"], effective_quality, window_cadences=5
    )
    expected_primary, expected_primary_trend, _ = centered_rolling_median_detrend(
        raw["primary_raw"], effective_quality, window_cadences=5
    )

    manifest = _build(tmp_path, tesscut, compact)

    with fits.open(tmp_path / "reference.fits", memmap=False, checksum=True) as hdul:
        assert hdul[0].header["TICID"] == WD1856_TIC
        assert hdul[0].header["SECTOR"] == 56
        assert hdul[0].header["CAMERA"] == 4
        assert hdul[0].header["CCD"] == 1
        assert hdul[0].header["RAWSHA"] == file_sha256(tesscut)
        assert hdul[0].header["MAPSHA"] == manifest[
            "compact_cadence_authority"
        ]["cadence_map_sha256"]
        assert hdul[0].header["QREFSHA"] == file_sha256(
            tmp_path / "cadence_reference.csv"
        )
        assert hdul[0].header["QMANISHA"] == file_sha256(
            tmp_path / "cadence_reference.json"
        )
        assert hdul[0].header["QPOLICY"] == manifest[
            "external_quality_overlay"
        ]["policy_contract"]
        assert hdul[0].header["SRCURL"] == SOURCE_URL
        assert hdul[0].header["ROLLWIN"] == 5
        assert hdul[0].header["APERSML"] == "1x1 central pixel"
        assert hdul[0].header["APERBIG"] == "3x3 centered sum"
        table = hdul["LIGHTCURVE"].data
        assert list(table.names) == manifest["output"]["columns"]
        np.testing.assert_array_equal(table["TICID"], WD1856_TIC)
        np.testing.assert_array_equal(table["SECTOR"], 56)
        np.testing.assert_array_equal(table["CADENCENO"], cadence)
        np.testing.assert_allclose(table["TIME"], raw["time"], atol=0, rtol=0)
        np.testing.assert_array_equal(table["QUALITY"], effective_quality)
        np.testing.assert_allclose(
            table["TESSCUT_FLUX_SML"], expected_small, atol=0, rtol=1e-14
        )
        np.testing.assert_allclose(
            table["TESSCUT_FLUX"], expected_primary, atol=0, rtol=1e-14
        )

    assert manifest["builder_version"] == BUILDER_VERSION
    assert manifest["tesscut_source"]["procver"] == ["test-astrocut-1.0"]
    assert manifest["compact_cadence_authority"]["datasets_read"] == [
        "time",
        "cadenceno",
        "quality",
        "orbitid",
    ]
    assert manifest["compact_cadence_authority"]["flux_read_or_used"] is False
    overlay = manifest["external_quality_overlay"]
    assert overlay["cadence_reference_table_sha256"] == file_sha256(
        tmp_path / "cadence_reference.csv"
    )
    assert overlay["cadence_reference_manifest_sha256"] == file_sha256(
        tmp_path / "cadence_reference.json"
    )
    assert overlay["compact_full_audit_counts"] == {
        "n_cad_total": 9,
        "n_cad_internal_bad": 1,
        "n_cad_external_bad": 2,
        "n_cad_external_only_bad": 2,
        "n_cad_effective_bad": 3,
    }
    assert overlay["tesscut_mapped_audit_counts"]["n_cad_effective_bad"] == 3
    assert manifest["extraction"]["detrending"]["minimum_baseline_cadences"] == minimum
    assert manifest["diagnostics"]["trend_1x1"] == {
        "min": float(np.min(expected_small_trend)),
        "median": float(np.median(expected_small_trend)),
        "max": float(np.max(expected_small_trend)),
    }
    assert manifest["diagnostics"]["trend_3x3"] == {
        "min": float(np.min(expected_primary_trend)),
        "median": float(np.median(expected_primary_trend)),
        "max": float(np.max(expected_primary_trend)),
    }


def test_output_and_sidecar_are_hash_bound_and_overwrite_is_explicit(
    tmp_path: Path,
) -> None:
    tesscut, compact, _, _ = _inputs(tmp_path)
    manifest = _build(tmp_path, tesscut, compact)
    published = json.loads((tmp_path / "reference.json").read_text())
    assert published == manifest
    assert published["output"]["sha256"] == file_sha256(tmp_path / "reference.fits")
    assert published["tesscut_source"]["sha256"] == file_sha256(tesscut)
    assert published["compact_cadence_authority"]["sha256"] == file_sha256(compact)
    assert len(published["compact_cadence_authority"]["cadence_map_sha256"]) == 64
    assert stat.S_IMODE((tmp_path / "reference.fits").stat().st_mode) == 0o640
    assert stat.S_IMODE((tmp_path / "reference.json").stat().st_mode) == 0o640

    with pytest.raises(FileExistsError, match="overwrite=True"):
        _build(tmp_path, tesscut, compact)
    replaced = _build(tmp_path, tesscut, compact, overwrite=True)
    assert replaced["output"]["sha256"] == file_sha256(tmp_path / "reference.fits")


@pytest.mark.parametrize(
    ("tesscut_options", "message"),
    [
        ({"sector": 57}, "SECTOR must be 56"),
        ({"camera": 3}, "CAMERA must be 4"),
        ({"ccd": 2}, "CCD must be 1"),
        ({"ra_deg": WD1856_TESSCUT_RA_DEG + 0.01}, "not WD 1856"),
        ({"include_wcs": False}, "celestial WCS"),
        ({"nonfinite_flux": True}, "entirely finite"),
        ({"duplicate_time": True}, "duplicate timestamps"),
        ({"size": 4}, "odd square"),
    ],
)
def test_invalid_tesscut_identity_schema_or_arrays_fail_closed(
    tmp_path: Path, tesscut_options: dict[str, object], message: str
) -> None:
    tesscut, compact, _, _ = _inputs(tmp_path, **tesscut_options)
    with pytest.raises(ValueError, match=message):
        _build(tmp_path, tesscut, compact)
    assert not (tmp_path / "reference.fits").exists()
    assert not (tmp_path / "reference.json").exists()


def test_cadence_mapping_requires_unique_subminute_matches(tmp_path: Path) -> None:
    tesscut = tmp_path / "tesscut.fits"
    raw = _write_tesscut(tesscut)
    compact = tmp_path / "compact.h5"
    cadence = _write_compact(compact, raw["time"], duplicate_time=True)
    _write_quality_reference(tmp_path, cadence)
    with pytest.raises(ValueError, match="duplicate timestamps"):
        _build(tmp_path, tesscut, compact)

    compact.unlink()
    cadence = _write_compact(compact, raw["time"], offset_seconds=61.0)
    _write_quality_reference(tmp_path, cadence)
    with pytest.raises(ValueError, match="above 60"):
        _build(tmp_path, tesscut, compact)


def test_source_mutation_rolls_back_without_publishing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tesscut, compact, _, _ = _inputs(tmp_path)
    original_hash = tesscut_module.file_sha256
    raw_path = tesscut.resolve()
    raw_calls = 0

    def changing_hash(path: Path) -> str:
        nonlocal raw_calls
        if Path(path).resolve() == raw_path:
            raw_calls += 1
            if raw_calls >= 2:
                return "f" * 64
        return original_hash(path)

    monkeypatch.setattr(tesscut_module, "file_sha256", changing_hash)
    with pytest.raises(RuntimeError, match="source mutated"):
        _build(tmp_path, tesscut, compact)
    assert not (tmp_path / "reference.fits").exists()
    assert not (tmp_path / "reference.json").exists()


def test_cli_builds_the_same_reference_contract(tmp_path: Path) -> None:
    tesscut, compact, _, _ = _inputs(tmp_path)
    output = tmp_path / "cli_reference.fits"
    manifest = tmp_path / "cli_reference.json"
    script = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "stage1_lightcurves"
        / "build_wd1856_tesscut_independent_extraction.py"
    )
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--tesscut-fits",
            str(tesscut),
            "--compact-lc",
            str(compact),
            "--cadence-reference-table",
            str(tmp_path / "cadence_reference.csv"),
            "--cadence-reference-manifest",
            str(tmp_path / "cadence_reference.json"),
            "--output-fits",
            str(output),
            "--manifest-json",
            str(manifest),
            "--source-url",
            SOURCE_URL,
            "--rolling-window-cadences",
            "5",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)
    assert payload == json.loads(manifest.read_text())
    assert payload["output"]["sha256"] == file_sha256(output)
