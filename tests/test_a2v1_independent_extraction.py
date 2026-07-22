from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.lightcurves.a2v1_independent_extraction as independent_module
from twirl.lightcurves.a2v1_cadence_reference import CADENCE_REFERENCE_COLUMNS
from twirl.lightcurves.a2v1_independent_extraction import (
    IndependentExtractionProvenance,
    build_wd1856_independent_metrics,
    parse_reference_flux_column_mappings,
)
from twirl.lightcurves.a2v1_qa import (
    WD1856_PERIOD_D,
    WD1856_T0_BJD,
    WD1856_TIC,
    file_sha256,
)
from twirl.lightcurves.a2v1_tier1_qa import (
    Tier1QAConfig,
    _dataframe_content_sha256,
    evaluate_independent_extraction,
)
from twirl.vetting.adp_only import ADP_ONLY_APERTURES


def _phase_minutes(time_bjd: np.ndarray) -> np.ndarray:
    return (
        ((time_bjd - WD1856_T0_BJD + 0.5 * WD1856_PERIOD_D) % WD1856_PERIOD_D)
        - 0.5 * WD1856_PERIOD_D
    ) * 1440.0


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    rng = np.random.default_rng(1856)
    n = 8_000
    cadence = np.arange(300_000, 300_000 + n, dtype=np.int64)
    start = WD1856_T0_BJD + 740 * WD1856_PERIOD_D - 0.25
    time_bjd = start + np.arange(n) * 200.0 / 86_400.0
    in_event = np.abs(_phase_minutes(time_bjd)) <= 4.0
    quality = np.zeros(n, dtype=np.int32)
    quality[:3] = 1

    common_noise = rng.normal(0.0, 0.0017, n)
    small = 1.0 + common_noise + rng.normal(0.0, 0.0005, n) - 0.54 * in_event
    primary = 1.0 + common_noise + rng.normal(0.0, 0.0007, n) - 0.50 * in_event
    compact = tmp_path / "s56_a2v1_adp_pair.h5"
    with h5py.File(compact, "w") as h5:
        h5.attrs["sector"] = 56
        h5.attrs["time_unit"] = "BJD - 2457000"
        targets = h5.create_group("targets")
        group = targets.create_group(f"{WD1856_TIC:016d}")
        group.attrs["tic"] = WD1856_TIC
        group.attrs["sector"] = 56
        group.attrs["camera"] = 4
        group.attrs["ccd"] = 1
        group.attrs["tessmag"] = 16.0
        group.create_dataset("time", data=time_bjd - 2_457_000.0)
        group.create_dataset("cadenceno", data=cadence)
        group.create_dataset("quality", data=quality)
        group.create_dataset(
            "orbitid", data=np.where(np.arange(n) < n // 2, 119, 120)
        )
        group.create_dataset(ADP_ONLY_APERTURES[0], data=small)
        group.create_dataset(ADP_ONLY_APERTURES[1], data=primary)

    keep = np.ones(n, dtype=bool)
    keep[::997] = False
    reference_quality = quality.copy()
    reference_quality[10:13] = 2
    reference_flux_small = (
        1.0
        + common_noise
        + rng.normal(0.0, 0.0006, n)
        - 0.53 * in_event
    )
    reference_flux_primary = (
        1.0
        + common_noise
        + rng.normal(0.0, 0.0008, n)
        - 0.51 * in_event
    )
    reference = tmp_path / "external_wd1856.csv"
    pd.DataFrame(
        {
            "TICID": np.full(np.count_nonzero(keep), WD1856_TIC, dtype=np.int64),
            "SECTOR": np.full(np.count_nonzero(keep), 56, dtype=np.int16),
            "CADENCENO": cadence[keep],
            "TIME": (time_bjd - 2_457_000.0)[keep],
            "QUALITY": reference_quality[keep],
            "REFERENCE_FLUX_SML": reference_flux_small[keep],
            "REFERENCE_FLUX": reference_flux_primary[keep],
        }
    ).to_csv(reference, index=False)
    _write_quality_reference(tmp_path, cadence)
    return compact, reference


def _write_quality_reference(root: Path, cadence: np.ndarray) -> tuple[Path, Path]:
    table_path = root / "cadence_reference.csv"
    manifest_path = root / "cadence_reference.json"
    orbitid = np.where(np.arange(len(cadence)) < len(cadence) // 2, 119, 120)
    spoc_quality = np.zeros(len(cadence), dtype=np.int64)
    qlp_quality = np.zeros(len(cadence), dtype=np.int64)
    spoc_quality[20] = 16
    qlp_quality[-20] = 1
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


def _quality_kwargs(root: Path) -> dict[str, Path]:
    return {
        "cadence_reference_table": root / "cadence_reference.csv",
        "cadence_reference_manifest": root / "cadence_reference.json",
    }


def _reference_flux_columns() -> dict[str, str]:
    return {
        ADP_ONLY_APERTURES[0]: "REFERENCE_FLUX_SML",
        ADP_ONLY_APERTURES[1]: "REFERENCE_FLUX",
    }


def test_reference_flux_mappings_require_exact_active_aperture_coverage(
    tmp_path: Path,
) -> None:
    specifications = [
        f"{ADP_ONLY_APERTURES[0]}=REFERENCE_FLUX_SML",
        f"{ADP_ONLY_APERTURES[1]}=REFERENCE_FLUX",
    ]
    assert parse_reference_flux_column_mappings(specifications) == (
        _reference_flux_columns()
    )
    with pytest.raises(ValueError, match="CURRENT=REFERENCE"):
        parse_reference_flux_column_mappings(["REFERENCE_FLUX"])
    with pytest.raises(ValueError, match="duplicate"):
        parse_reference_flux_column_mappings(
            [specifications[0], specifications[0], specifications[1]]
        )
    with pytest.raises(ValueError, match="exactly cover"):
        parse_reference_flux_column_mappings(specifications[:1])
    with pytest.raises(ValueError, match="distinct external column"):
        parse_reference_flux_column_mappings(
            [
                f"{ADP_ONLY_APERTURES[0]}=REFERENCE_FLUX",
                f"{ADP_ONLY_APERTURES[1]}=REFERENCE_FLUX",
            ]
        )

    compact, reference = _write_inputs(tmp_path)
    with pytest.raises(ValueError, match="exactly cover"):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=reference,
            reference_product=reference,
            **_quality_kwargs(tmp_path),
            metrics_csv=tmp_path / "incomplete_metrics.csv",
            manifest_json=tmp_path / "incomplete_manifest.json",
            provenance=_provenance(),
            reference_time_system="BTJD",
            reference_flux_columns={
                ADP_ONLY_APERTURES[0]: "REFERENCE_FLUX_SML"
            },
        )


def _provenance(**changes: str) -> IndependentExtractionProvenance:
    values = {
        "current_repository": "mit-kavli-institute/tess-gaia-light-curve+TWIRL",
        "current_revision": "tglc-a2v1-revision",
        "reference_extractor_family": "NASA SPOC aperture photometry",
        "reference_repository": "NASA-SPOC/pipeline-release",
        "reference_revision": "sector-56-release",
        "reference_pixel_source": "official independently calibrated SPOC target pixels",
        "reference_product_aperture": "SAP_FLUX pipeline aperture",
        "reference_target_id": f"TIC {WD1856_TIC}",
        "independence_basis": (
            "independent calibration, background model, aperture, and extraction code"
        ),
    }
    values.update(changes)
    return IndependentExtractionProvenance(**values)


def test_builder_writes_exact_tier1_metrics_and_manifest(tmp_path: Path) -> None:
    compact, reference = _write_inputs(tmp_path)
    metrics_path = tmp_path / "independent_metrics.csv"
    manifest_path = tmp_path / "independent_manifest.json"
    metrics, manifest = build_wd1856_independent_metrics(
        compact_lc=compact,
        reference_table=reference,
        reference_product=reference,
        **_quality_kwargs(tmp_path),
        metrics_csv=metrics_path,
        manifest_json=manifest_path,
        provenance=_provenance(),
        reference_time_system="BTJD",
        reference_flux_columns=_reference_flux_columns(),
    )

    assert list(metrics["aperture"]) == list(ADP_ONLY_APERTURES)
    assert dict(
        zip(metrics["aperture"], metrics["reference_flux_column"], strict=True)
    ) == _reference_flux_columns()
    assert metrics["n_common_cadences"].eq(7_991).all()
    assert metrics["n_in_event_cadences"].min() >= 4
    assert metrics["n_out_of_event_cadences"].min() >= 100
    assert manifest["metrics_file_sha256"] == file_sha256(metrics_path)
    assert manifest["metrics_content_sha256"] == _dataframe_content_sha256(
        pd.read_csv(metrics_path, float_precision="round_trip")
    )
    assert manifest["reference_table_identity"]["tic"] == WD1856_TIC
    assert manifest["reference_table_identity"]["sector"] == 56
    assert manifest["reference_table_identity"]["format"] == "csv"
    assert manifest["reference_flux_columns"] == _reference_flux_columns()
    assert manifest["comparison_mode"] == "signal_timing_only"
    assert manifest["contract_version"] == "s56_a2v1_independent_extraction_v2"
    assert manifest["reference_columns"]["flux_by_current_aperture"] == (
        _reference_flux_columns()
    )
    overlay = manifest["external_quality_overlay"]
    assert overlay["applied_before_common_cadence_metrics"] is True
    assert overlay["cadence_reference_table_sha256"] == file_sha256(
        tmp_path / "cadence_reference.csv"
    )
    assert overlay["cadence_reference_manifest_sha256"] == file_sha256(
        tmp_path / "cadence_reference.json"
    )
    assert overlay["current_compact_full_audit_counts"] == {
        "n_cad_total": 8_000,
        "n_cad_internal_bad": 3,
        "n_cad_external_bad": 2,
        "n_cad_external_only_bad": 2,
        "n_cad_effective_bad": 5,
    }
    assert overlay["current_common_audit_counts"]["n_cad_total"] == 7_991
    assert overlay["reference_common_audit_counts"]["n_cad_external_bad"] == 2
    assert json.loads(manifest_path.read_text()) == manifest

    compact_sha256 = file_sha256(compact)
    config = replace(
        Tier1QAConfig(),
        independent_contract_template="s{sector}_a2v1_independent_extraction_v2",
        expected_compact_sha256=compact_sha256,
        expected_cadence_reference_sha256=file_sha256(
            tmp_path / "cadence_reference.csv"
        ),
        expected_cadence_reference_manifest_sha256=file_sha256(
            tmp_path / "cadence_reference.json"
        ),
        expected_independent_reference_product_sha256=manifest[
            "reference_product_sha256"
        ],
    )
    gate = evaluate_independent_extraction(
        pd.read_csv(metrics_path, float_precision="round_trip"),
        manifest,
        config,
        sector=56,
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        current_compact_sha256=compact_sha256,
    )
    assert gate["status"] == "pass"
    assert gate["wd1856"]["status"] == "pass"

    legacy_manifest = dict(manifest)
    legacy_manifest.pop("reference_flux_columns")
    legacy_metrics = metrics.drop(columns="reference_flux_column")
    legacy_gate = evaluate_independent_extraction(
        legacy_metrics,
        legacy_manifest,
        config,
        sector=56,
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        current_compact_sha256=compact_sha256,
    )
    assert legacy_gate["status"] == "fail"
    assert any("reference_flux_columns" in reason for reason in legacy_gate["reasons"])

    reused_manifest = json.loads(json.dumps(manifest))
    reused_manifest["reference_flux_columns"] = {
        aperture: "REFERENCE_FLUX" for aperture in ADP_ONLY_APERTURES
    }
    reused_manifest["reference_columns"]["flux_by_current_aperture"] = dict(
        reused_manifest["reference_flux_columns"]
    )
    reused_gate = evaluate_independent_extraction(
        metrics,
        reused_manifest,
        config,
        sector=56,
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        current_compact_sha256=compact_sha256,
    )
    assert reused_gate["status"] == "fail"
    assert any("reuses one external column" in reason for reason in reused_gate["reasons"])


def test_builder_rejects_nonindependent_or_misaligned_reference(tmp_path: Path) -> None:
    compact, reference = _write_inputs(tmp_path)
    with pytest.raises(ValueError, match="outside the TGLC/TWIRL family"):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=reference,
            reference_product=reference,
            **_quality_kwargs(tmp_path),
            metrics_csv=tmp_path / "same_family.csv",
            manifest_json=tmp_path / "same_family.json",
            provenance=_provenance(reference_extractor_family="another TGLC run"),
            reference_time_system="BTJD",
            reference_flux_columns=_reference_flux_columns(),
        )

    shifted = pd.read_csv(reference)
    shifted["TIME"] += 30.0 / 86_400.0
    shifted_path = tmp_path / "shifted_external.csv"
    shifted.to_csv(shifted_path, index=False)
    with pytest.raises(ValueError, match="timestamps differ"):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=shifted_path,
            reference_product=shifted_path,
            **_quality_kwargs(tmp_path),
            metrics_csv=tmp_path / "shifted.csv",
            manifest_json=tmp_path / "shifted.json",
            provenance=_provenance(),
            reference_time_system="BTJD",
            reference_flux_columns=_reference_flux_columns(),
        )
    assert not (tmp_path / "shifted.csv").exists()
    assert not (tmp_path / "shifted.json").exists()


def test_builder_rejects_reference_changed_during_measurement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    compact, reference = _write_inputs(tmp_path)
    original_file_sha256 = independent_module.file_sha256
    reference_hash_calls = 0

    def changing_file_sha256(path: Path) -> str:
        nonlocal reference_hash_calls
        if Path(path).resolve() == reference.resolve():
            reference_hash_calls += 1
            if reference_hash_calls == 2:
                reference.write_text(reference.read_text() + "\n")
        return original_file_sha256(path)

    monkeypatch.setattr(
        independent_module, "file_sha256", changing_file_sha256
    )
    metrics_path = tmp_path / "changed_metrics.csv"
    manifest_path = tmp_path / "changed_manifest.json"
    with pytest.raises(ValueError, match="changed during the build"):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=reference,
            reference_product=reference,
            **_quality_kwargs(tmp_path),
            metrics_csv=metrics_path,
            manifest_json=manifest_path,
            provenance=_provenance(),
            reference_time_system="BTJD",
            reference_flux_columns=_reference_flux_columns(),
        )
    assert not metrics_path.exists()
    assert not manifest_path.exists()


def test_builder_rejects_unbound_external_quality_evidence(tmp_path: Path) -> None:
    compact, reference = _write_inputs(tmp_path)
    cadence_reference = tmp_path / "cadence_reference.csv"
    cadence_reference.write_text(
        cadence_reference.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="table hash mismatch"):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=reference,
            reference_product=reference,
            **_quality_kwargs(tmp_path),
            metrics_csv=tmp_path / "unbound_metrics.csv",
            manifest_json=tmp_path / "unbound_manifest.json",
            provenance=_provenance(),
            reference_time_system="BTJD",
            reference_flux_columns=_reference_flux_columns(),
        )


def test_data_derived_wrong_reference_ephemeris_fails_live_gate(
    tmp_path: Path,
) -> None:
    compact, reference = _write_inputs(tmp_path)
    wrong = pd.read_csv(reference)
    time_bjd = wrong["TIME"].to_numpy(dtype=float) + 2_457_000.0
    wrong_period = WD1856_PERIOD_D * 1.003
    wrong_phase_min = (
        ((time_bjd - WD1856_T0_BJD + 0.5 * wrong_period) % wrong_period)
        - 0.5 * wrong_period
    ) * 1440.0
    rng = np.random.default_rng(56_1856)
    for flux_column, depth in (
        ("REFERENCE_FLUX_SML", 0.53),
        ("REFERENCE_FLUX", 0.51),
    ):
        wrong[flux_column] = (
            1.0
            + rng.normal(0.0, 0.0015, len(wrong))
            - depth * (np.abs(wrong_phase_min) <= 4.0)
        )
    wrong_reference = tmp_path / "wrong_ephemeris_external.csv"
    wrong.to_csv(wrong_reference, index=False)

    metrics_path = tmp_path / "wrong_ephemeris_metrics.csv"
    manifest_path = tmp_path / "wrong_ephemeris_manifest.json"
    metrics, manifest = build_wd1856_independent_metrics(
        compact_lc=compact,
        reference_table=wrong_reference,
        reference_product=wrong_reference,
        **_quality_kwargs(tmp_path),
        metrics_csv=metrics_path,
        manifest_json=manifest_path,
        provenance=_provenance(),
        reference_time_system="BTJD",
        reference_flux_columns=_reference_flux_columns(),
    )
    recovered_error = np.abs(metrics["reference_period_d"] / WD1856_PERIOD_D - 1.0)
    assert recovered_error.min() > 5.0e-4
    compact_sha256 = file_sha256(compact)
    config = replace(
        Tier1QAConfig(),
        independent_contract_template="s{sector}_a2v1_independent_extraction_v2",
        expected_compact_sha256=compact_sha256,
        expected_cadence_reference_sha256=file_sha256(
            tmp_path / "cadence_reference.csv"
        ),
        expected_cadence_reference_manifest_sha256=file_sha256(
            tmp_path / "cadence_reference.json"
        ),
        expected_independent_reference_product_sha256=manifest[
            "reference_product_sha256"
        ],
    )
    gate = evaluate_independent_extraction(
        metrics,
        manifest,
        config,
        sector=56,
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        current_compact_sha256=compact_sha256,
    )
    assert gate["status"] == "fail"
    assert gate["wd1856"]["status"] == "fail"


@pytest.mark.parametrize(
    ("column", "wrong_value", "message"),
    (
        ("TICID", 123456789, "WD 1856"),
        ("SECTOR", 57, "Sector 56"),
    ),
)
def test_nonfits_reference_identity_fails_closed(
    tmp_path: Path,
    column: str,
    wrong_value: int,
    message: str,
) -> None:
    compact, reference = _write_inputs(tmp_path)
    table = pd.read_csv(reference)
    table[column] = wrong_value
    wrong_identity = tmp_path / f"wrong_{column.lower()}.csv"
    table.to_csv(wrong_identity, index=False)
    with pytest.raises(ValueError, match=message):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=wrong_identity,
            reference_product=wrong_identity,
            **_quality_kwargs(tmp_path),
            metrics_csv=tmp_path / f"wrong_{column.lower()}_metrics.csv",
            manifest_json=tmp_path / f"wrong_{column.lower()}_manifest.json",
            provenance=_provenance(),
            reference_time_system="BTJD",
            reference_flux_columns=_reference_flux_columns(),
        )


def test_nonfits_reference_requires_identity_columns(tmp_path: Path) -> None:
    compact, reference = _write_inputs(tmp_path)
    table = pd.read_csv(reference).drop(columns=["TICID"])
    missing_identity = tmp_path / "missing_identity.csv"
    table.to_csv(missing_identity, index=False)
    with pytest.raises(ValueError, match="lacks explicit target/sector metadata"):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=missing_identity,
            reference_product=missing_identity,
            **_quality_kwargs(tmp_path),
            metrics_csv=tmp_path / "missing_identity_metrics.csv",
            manifest_json=tmp_path / "missing_identity_manifest.json",
            provenance=_provenance(),
            reference_time_system="BTJD",
            reference_flux_columns=_reference_flux_columns(),
        )


def test_separate_source_product_identity_is_inspected(tmp_path: Path) -> None:
    compact, reference = _write_inputs(tmp_path)
    wrong_product_table = pd.read_csv(reference)
    wrong_product_table["SECTOR"] = 57
    wrong_product = tmp_path / "wrong_source_product.csv"
    wrong_product_table.to_csv(wrong_product, index=False)
    with pytest.raises(ValueError, match="external source product.*Sector 56"):
        build_wd1856_independent_metrics(
            compact_lc=compact,
            reference_table=reference,
            reference_product=wrong_product,
            **_quality_kwargs(tmp_path),
            metrics_csv=tmp_path / "wrong_product_metrics.csv",
            manifest_json=tmp_path / "wrong_product_manifest.json",
            provenance=_provenance(),
            reference_time_system="BTJD",
            reference_flux_columns=_reference_flux_columns(),
        )


def test_builder_refuses_overwrite_and_cli_exposes_required_provenance(
    tmp_path: Path,
) -> None:
    compact, reference = _write_inputs(tmp_path)
    metrics_path = tmp_path / "metrics.csv"
    manifest_path = tmp_path / "manifest.json"
    kwargs = {
        "compact_lc": compact,
        "reference_table": reference,
        "reference_product": reference,
        **_quality_kwargs(tmp_path),
        "metrics_csv": metrics_path,
        "manifest_json": manifest_path,
        "provenance": _provenance(),
        "reference_time_system": "BTJD",
        "reference_flux_columns": _reference_flux_columns(),
    }
    build_wd1856_independent_metrics(**kwargs)
    with pytest.raises(FileExistsError, match="overwrite"):
        build_wd1856_independent_metrics(**kwargs)

    script = (
        Path(__file__).resolve().parents[1]
        / "scripts/stage1_lightcurves/build_a2v1_independent_extraction.py"
    )
    cli_metrics = tmp_path / "cli_metrics.csv"
    cli_manifest = tmp_path / "cli_manifest.json"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--compact-lc",
            str(compact),
            "--reference-table",
            str(reference),
            "--reference-product",
            str(reference),
            "--cadence-reference-table",
            str(tmp_path / "cadence_reference.csv"),
            "--cadence-reference-manifest",
            str(tmp_path / "cadence_reference.json"),
            "--reference-time-system",
            "BTJD",
            "--reference-flux-column",
            f"{ADP_ONLY_APERTURES[0]}=REFERENCE_FLUX_SML",
            "--reference-flux-column",
            f"{ADP_ONLY_APERTURES[1]}=REFERENCE_FLUX",
            "--current-repository",
            "mit-kavli-institute/tess-gaia-light-curve+TWIRL",
            "--current-revision",
            "tglc-a2v1-revision",
            "--reference-extractor-family",
            "NASA SPOC aperture photometry",
            "--reference-repository",
            "NASA-SPOC/pipeline-release",
            "--reference-revision",
            "sector-56-release",
            "--reference-pixel-source",
            "official independently calibrated SPOC target pixels",
            "--reference-product-aperture",
            "SAP_FLUX pipeline aperture",
            "--reference-target-id",
            f"TIC {WD1856_TIC}",
            "--independence-basis",
            "independent calibration background aperture and extraction code",
            "--metrics-csv",
            str(cli_metrics),
            "--manifest-json",
            str(cli_manifest),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert cli_metrics.is_file()
    assert cli_manifest.is_file()
    cli_payload = json.loads(cli_manifest.read_text())
    assert cli_payload["metrics_file_sha256"] == file_sha256(cli_metrics)
    assert cli_payload["independent"] is True
    assert cli_payload["reference_flux_columns"] == _reference_flux_columns()
