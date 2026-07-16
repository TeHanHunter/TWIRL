from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json

import h5py
import numpy as np
import pandas as pd
import pytest

from twirl.lightcurves.a2v1_qa import (
    A2V1_PHOTOMETRIC_QA_VERSION,
    WD1856_PERIOD_D,
    WD1856_T0_BJD,
    WD1856_TIC,
    file_sha256,
)
from twirl.lightcurves.a2v1_tier1_qa import (
    TIER1_QA_CONTRACT_VERSION,
    Tier1QAConfig,
    _dataframe_content_sha256,
    audit_compact_population,
    evaluate_aperture_outliers,
    evaluate_fixed_injections,
    evaluate_independent_extraction,
    evaluate_tier1_gates,
    injection_metadata_sha256,
    load_tier1_config,
    run_a2v1_tier1_qa,
    summarize_fixed_injection_shards,
    write_strict_json,
)
from twirl.vetting.adp_only import ADP_ONLY_APERTURES


TEST_COMPACT_SHA256 = "a" * 64
TEST_INJECTION_IDS = tuple(f"inj_00_{index:03d}" for index in range(8))
TEST_INJECTION_SELECTION_SHA256 = hashlib.sha256(
    ("\n".join(TEST_INJECTION_IDS) + "\n").encode()
).hexdigest()
TEST_INJECTION_SOURCE_SHA256 = "9" * 64
TEST_INJECTION_METADATA_SHA256 = injection_metadata_sha256(
    pd.DataFrame(
        {
            "injection_id": TEST_INJECTION_IDS,
            "tic": [5000 + index for index in range(8)],
            "tmag": [19.4 if index < 4 else 18.0 for index in range(8)],
            "camera": [1 + index % 4 for index in range(8)],
            "ccd": [1 + (index // 4) % 4 for index in range(8)],
            "period_d": np.geomspace(0.2, 2.0, 8),
            "duration_min": np.geomspace(2.0, 10.0, 8),
            "model_depth": np.geomspace(0.01, 0.1, 8),
            "shard_index": np.zeros(8, dtype=int),
        }
    )
)


def _write_compact(path, *, n_targets: int = 40, n_cadences: int = 200) -> None:
    rng = np.random.default_rng(56)
    with h5py.File(path, "w") as h5:
        h5.attrs["sector"] = 56
        targets = h5.create_group("targets")
        for index in range(n_targets):
            tic = WD1856_TIC if index == 0 else 1000 + index
            group = targets.create_group(f"{tic:016d}")
            group.attrs["tic"] = tic
            group.attrs["sector"] = 56
            group.attrs["camera"] = 4 if index == 0 else 1 + index % 2
            group.attrs["ccd"] = 1 if index == 0 else 1 + index % 4
            group.attrs["tessmag"] = 17.0 + 4.0 * index / max(n_targets - 1, 1)
            group.create_dataset("time", data=np.arange(n_cadences) / 432.0)
            group.create_dataset("cadenceno", data=np.arange(n_cadences))
            quality = np.zeros(n_cadences, dtype=np.int32)
            quality[:2] = 1
            group.create_dataset("quality", data=quality)
            group.create_dataset(
                "orbitid",
                data=np.where(np.arange(n_cadences) < n_cadences // 2, 119, 120),
            )
            sigma = 0.002 * 10 ** (0.12 * (float(group.attrs["tessmag"]) - 17.0))
            base = 1.0 + rng.normal(0.0, sigma, n_cadences)
            group.create_dataset(ADP_ONLY_APERTURES[0], data=base)
            group.create_dataset(
                ADP_ONLY_APERTURES[1], data=1.0 + 1.05 * (base - 1.0)
            )


def _cadence_reference(*, n_cadences: int = 200) -> pd.DataFrame:
    rows: list[dict] = []
    for camera, ccd in ((1, 1), (1, 3), (2, 2), (2, 4), (4, 1)):
        for cadence in range(n_cadences):
            rows.append(
                {
                    "sector": 56,
                    "orbitid": 119 if cadence < n_cadences // 2 else 120,
                    "camera": camera,
                    "ccd": ccd,
                    "cadenceno": cadence,
                    "quality": 1 if cadence < 2 else 0,
                }
            )
    return pd.DataFrame(rows)


def _write_cadence_evidence(path) -> tuple[pd.DataFrame, dict]:
    table = _cadence_reference()
    table.to_csv(path, index=False)
    manifest = {
        "contract_version": "s56_a2v1_cadence_reference_v1",
        "sector": 56,
        "cadence_authority": "qlp_cam_quat",
        "quality_authority": "spoc_quality_flags",
        "table_sha256": file_sha256(path),
        "n_rows": len(table),
        "detectors": [
            "cam1_ccd1",
            "cam1_ccd3",
            "cam2_ccd2",
            "cam2_ccd4",
            "cam4_ccd1",
        ],
        "orbits": [119, 120],
        "source_file_sha256": {"authoritative-input": "c" * 64},
    }
    return table, manifest


def _write_injection_shard(
    path,
    *,
    n_ids: int = 8,
    retention: float = 0.96,
    shard_index: int = 0,
) -> None:
    rng = np.random.default_rng(123 + shard_index)
    with h5py.File(path, "w") as h5:
        h5.attrs["contract_version"] = "s56_a2v1_fresh_injection_pair_v1"
        h5.attrs["shard_index"] = shard_index
        h5.attrs["n_shards"] = 40
        h5.attrs["n_injections"] = n_ids
        h5.attrs["selection_mode"] = "shard"
        h5.attrs["source_adp_h5"] = "/test/s56_A2v1_adp_pair_rebuilt.h5"
        groups = h5.create_group("injections")
        periods = np.geomspace(0.2, 2.0, n_ids)
        durations = np.geomspace(2.0, 10.0, n_ids)
        depths = np.geomspace(0.01, 0.1, n_ids)
        for index in range(n_ids):
            group = groups.create_group(f"inj_{shard_index:02d}_{index:03d}")
            group.attrs["tic"] = 5000 + shard_index * 1000 + index
            group.attrs["tmag"] = 19.4 if index < n_ids // 2 else 18.0
            group.attrs["camera"] = 1 + index % 4
            group.attrs["ccd"] = 1 + (index // 4) % 4
            group.attrs["period_d"] = periods[index]
            group.attrs["duration_min"] = durations[index]
            group.attrs["model_depth"] = depths[index]
            group.attrs["radius_rearth"] = 0.2 + index
            n = 200
            model = np.ones(n)
            model[40:45] = np.array([0.97, 0.93, 0.90, 0.93, 0.97])
            group.create_dataset("quality", data=np.zeros(n, dtype=np.int32))
            group.create_dataset("transit_model", data=model)
            for aperture in ADP_ONLY_APERTURES:
                original = 1.0 + rng.normal(0.0, 1.0e-4, n)
                injected = original - retention * (1.0 - model) + 2.0e-4
                group.create_dataset(f"{aperture}_original", data=original)
                group.create_dataset(f"{aperture}_injected", data=injected)


def _small_injection_config(**changes) -> Tier1QAConfig:
    base = {
        "expected_compact_sha256": TEST_COMPACT_SHA256,
        "injection_required_shard_indices": (0,),
        "injection_expected_ids": 8,
        "expected_injection_selection_sha256": TEST_INJECTION_SELECTION_SHA256,
        "expected_injection_metadata_sha256": TEST_INJECTION_METADATA_SHA256,
        "expected_injection_shard_sha256": ("0" * 64,),
        "expected_injection_source_adp_sha256": TEST_INJECTION_SOURCE_SHA256,
        "expected_injection_parity_report_sha256": "8" * 64,
        "expected_independent_reference_product_sha256": "e" * 64,
        "min_faint_injection_ids": 4,
        "min_injection_detectors": 4,
        "min_injection_period_ratio": 5.0,
        "min_injection_duration_ratio": 2.0,
        "min_injection_depth_ratio": 5.0,
    }
    base.update(changes)
    return replace(Tier1QAConfig(), **base)


def _passing_independent(
    *, current_compact_sha256: str = TEST_COMPACT_SHA256
) -> tuple[pd.DataFrame, dict]:
    metrics = pd.DataFrame(
        {
            "tic": [WD1856_TIC, WD1856_TIC],
            "sector": [56, 56],
            "aperture": list(ADP_ONLY_APERTURES),
            "detector": ["cam4_ccd1", "cam4_ccd1"],
            "n_current_cadences": [11_775, 11_775],
            "n_reference_cadences": [11_770, 11_770],
            "n_common_cadences": [11_760, 11_760],
            "current_scatter_ppm": [100_000.0, 105_000.0],
            "reference_scatter_ppm": [110_000.0, 112_000.0],
            "current_period_d": [WD1856_PERIOD_D * 1.00005] * 2,
            "reference_period_d": [WD1856_PERIOD_D * 0.99995] * 2,
            "current_t0_bjd": [WD1856_T0_BJD + 10 * WD1856_PERIOD_D] * 2,
            "reference_t0_bjd": [WD1856_T0_BJD + 10 * WD1856_PERIOD_D] * 2,
            "current_depth": [0.57, 0.55],
            "reference_depth": [0.52, 0.51],
            "current_bls_duration_min": [8.0, 8.0],
            "reference_bls_duration_min": [8.0, 8.0],
            "current_bls_depth": [0.57, 0.55],
            "reference_bls_depth": [0.52, 0.51],
            "current_bls_depth_snr": [20.0, 20.0],
            "reference_bls_depth_snr": [20.0, 20.0],
            "current_bls_power": [100.0, 100.0],
            "reference_bls_power": [95.0, 95.0],
            "n_common_quality0_finite": [11_700, 11_700],
            "n_in_event_cadences": [20, 20],
            "n_out_of_event_cadences": [11_000, 11_000],
            "max_abs_time_delta_s": [1.0, 1.0],
        }
    )
    manifest = {
        "contract_version": "s56_a2v1_independent_extraction_v1",
        "sector": 56,
        "tic": WD1856_TIC,
        "independent": True,
        "independence_basis": "independent pixels and extraction code",
        "current_extractor_family": "MIT TGLC A2v1",
        "reference_extractor_family": "MAST QLP",
        "current_repository": "TeHanHunter/TESS_Gaia_Light_Curve",
        "reference_repository": "mit-qlp",
        "current_revision": "abc123",
        "reference_revision": "release-1",
        "pixel_source": "independently calibrated TICA cutout",
        "cadence_match_policy": "exact common CADENCENO intersection",
        "scatter_definition": "out-of-transit robust MAD in ppm on common cadences",
        "depth_definition": "median in-event deficit on a fixed ephemeris",
        "ephemeris_recovery_definition": (
            "independent bounded BoxLeastSquares recovery on matched cadences"
        ),
        "current_apertures": list(ADP_ONLY_APERTURES),
        "reference_product_aperture": "official pipeline aperture",
        "reference_target_id": f"TIC {WD1856_TIC}",
        "current_compact_sha256": current_compact_sha256,
        "current_product_sha256": current_compact_sha256,
        "reference_product_sha256": "e" * 64,
        "reference_table_sha256": "1" * 64,
        "reference_table_identity": {
            "format": "csv",
            "tic": WD1856_TIC,
            "sector": 56,
            "identity_source": {"tic": "TICID", "sector": "SECTOR"},
        },
        "reference_product_identity": {
            "format": "fits",
            "tic": WD1856_TIC,
            "sector": 56,
            "identity_source": {"tic": "FITS:TICID", "sector": "FITS:SECTOR"},
        },
        "fixed_ephemeris": {
            "period_d": WD1856_PERIOD_D,
            "t0_bjd": WD1856_T0_BJD,
        },
        "bounded_bls": {
            "period_relative_half_width": 0.01,
            "n_periods": 4001,
            "durations_min": [6.0, 8.0, 10.0],
            "oversample": 20,
            "min_depth_snr": 5.0,
        },
        "max_time_delta_seconds_allowed": 5.0,
        "metrics_file_sha256": "f" * 64,
        "metrics_content_sha256": _dataframe_content_sha256(metrics),
    }
    return metrics, manifest


def _passing_tier0(compact_sha256: str) -> dict:
    return {
        "sector": 56,
        "passed": True,
        "contract_version": A2V1_PHOTOMETRIC_QA_VERSION,
        "qa_tier": "tier0_integrity_and_benchmark",
        "provenance": {
            "compact_lc_sha256": compact_sha256,
            "schema_summary_sha256": "b" * 64,
            "detrended_apertures": list(ADP_ONLY_APERTURES),
        },
    }


def _write_injection_parity(path, compact_sha256: str) -> None:
    payload = {
        "passed": True,
        "failures": [],
        "hashes": {
            "s56_A2v1_adp_pair.h5": compact_sha256,
            "s56_A2v1_adp_pair_rebuilt.h5": TEST_INJECTION_SOURCE_SHA256,
        },
        "compact_rebuild_parity": {
            "passed": True,
            "reference_h5": "/test/s56_A2v1_adp_pair.h5",
            "active_h5": "/test/s56_A2v1_adp_pair_rebuilt.h5",
            "n_dataset_comparisons": 188_700,
            "n_mismatched_targets": 0,
            "n_missing_active_targets": 0,
            "n_extra_active_targets": 0,
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def test_tier1_config_is_strict_and_scoped(tmp_path) -> None:
    payload = asdict(Tier1QAConfig())
    path = tmp_path / "config.json"
    path.write_text(json.dumps(payload))
    config = load_tier1_config(path)
    assert config.scope == "active_search_pair"
    assert config.contract_version == TIER1_QA_CONTRACT_VERSION
    assert config.injection_required_shard_indices == (0, 13, 26, 39)

    payload["unknown_threshold"] = 1
    path.write_text(json.dumps(payload))
    with pytest.raises(KeyError, match="unknown"):
        load_tier1_config(path)

    with pytest.raises(ValueError, match="promotion is disabled"):
        replace(Tier1QAConfig(), promotion_enabled=True).validate()


def test_compact_population_uses_authoritative_cadence_reference(tmp_path) -> None:
    compact = tmp_path / "compact.h5"
    _write_compact(compact)
    targets, apertures, pairs = audit_compact_population(
        compact,
        apertures=ADP_ONLY_APERTURES,
        cadence_reference=_cadence_reference(),
    )
    assert len(targets) == 40
    assert len(apertures) == 80
    assert len(pairs) == 40
    assert targets["missing_cadence_fraction"].max() == 0.0
    assert targets["n_quality_mismatches"].sum() == 0
    assert apertures["finite_quality0_fraction"].min() == 1.0
    assert pairs["mad_ratio"].between(0.9, 1.2).all()


def test_model_weighted_injection_retention_and_full_shard_gate(tmp_path) -> None:
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    metrics, manifest = summarize_fixed_injection_shards([shard])
    assert len(metrics) == 16
    assert metrics["status"].eq("ok").all()
    assert np.allclose(metrics["depth_retention_fraction"], 0.96, atol=1.0e-8)

    config = _small_injection_config(
        expected_injection_shard_sha256=(file_sha256(shard),)
    )
    gate = evaluate_fixed_injections(metrics, manifest, sector=56, config=config)
    assert gate["status"] == "pass"
    assert gate["n_faint_injection_ids"] == 4
    assert gate["period_support_ratio"] >= 5.0

    wrong_shard = dict(manifest)
    wrong_shard["shard_indices"] = [1]
    assert (
        evaluate_fixed_injections(metrics, wrong_shard, sector=56, config=config)[
            "status"
        ]
        == "fail"
    )

    substituted_shard = dict(manifest)
    substituted_shard["shard_index_sha256"] = {"0": "b" * 64}
    assert (
        evaluate_fixed_injections(
            metrics, substituted_shard, sector=56, config=config
        )["status"]
        == "fail"
    )

    missing_pair = metrics.loc[
        ~(
            metrics["injection_id"].eq("inj_00_000")
            & metrics["aperture"].eq(ADP_ONLY_APERTURES[1])
        )
    ]
    assert (
        evaluate_fixed_injections(
            missing_pair, manifest, sector=56, config=config
        )["status"]
        == "fail"
    )

    substituted = metrics.copy()
    substituted["injection_id"] = "other_" + substituted["injection_id"].astype(str)
    substituted_manifest = dict(manifest)
    substituted_manifest["metrics_content_sha256"] = _dataframe_content_sha256(
        substituted
    )
    assert (
        evaluate_fixed_injections(
            substituted, substituted_manifest, sector=56, config=config
        )["status"]
        == "fail"
    )

    remapped = metrics.copy()
    remapped["tic"] = pd.to_numeric(remapped["tic"]) + 1_000_000
    remapped["period_d"] = pd.to_numeric(remapped["period_d"]) * 1.5
    remapped_manifest = dict(manifest)
    remapped_manifest["metrics_content_sha256"] = _dataframe_content_sha256(remapped)
    remapped_manifest["metadata_sha256"] = injection_metadata_sha256(remapped)
    assert (
        evaluate_fixed_injections(
            remapped, remapped_manifest, sector=56, config=config
        )["status"]
        == "fail"
    )


def test_independent_gate_rejects_same_family_and_unphysical_depth() -> None:
    metrics, manifest = _passing_independent()
    passed = evaluate_independent_extraction(
        metrics,
        manifest,
        _small_injection_config(),
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
    )
    assert passed["status"] == "pass"
    assert (
        passed["wd1856"]["apertures"][0]["reference_epoch_residual_min"]
        < 1.0e-6
    )

    same_family = dict(manifest)
    same_family["reference_extractor_family"] = "MIT TGLC A2v1"
    failed = evaluate_independent_extraction(
        metrics,
        same_family,
        _small_injection_config(),
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
    )
    assert failed["status"] == "fail"
    assert any("not independent" in reason for reason in failed["reasons"])

    bad_depth = metrics.copy()
    bad_depth.loc[0, "reference_depth"] = 1.2
    bad_manifest = dict(manifest)
    bad_manifest["metrics_content_sha256"] = _dataframe_content_sha256(bad_depth)
    assert (
        evaluate_independent_extraction(
            bad_depth,
            bad_manifest,
            _small_injection_config(),
            catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        )["status"]
        == "fail"
    )

    bad_detector = metrics.copy()
    bad_detector["detector"] = "totally_not_a_detector"
    bad_detector_manifest = dict(manifest)
    bad_detector_manifest["metrics_content_sha256"] = _dataframe_content_sha256(
        bad_detector
    )
    assert (
        evaluate_independent_extraction(
            bad_detector,
            bad_detector_manifest,
            _small_injection_config(),
            catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        )["status"]
        == "fail"
    )


def test_aperture_gate_rejects_anticorrelated_channels() -> None:
    pairs = pd.DataFrame(
        {
            "tic": np.arange(100),
            "mad_ratio": np.ones(100),
            "correlation": np.full(100, -0.8),
        }
    )
    assert evaluate_aperture_outliers(pairs, Tier1QAConfig())["status"] == "fail"


def test_active_scope_can_be_enrichment_ready_but_not_science_ready(tmp_path) -> None:
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    config = _small_injection_config(
        min_population_targets=40,
        min_targets_per_magnitude_bin=3,
        expected_injection_shard_sha256=(file_sha256(shard),),
    )
    compact = tmp_path / "compact.h5"
    _write_compact(compact)
    targets, apertures, pairs = audit_compact_population(
        compact,
        apertures=config.apertures,
        cadence_reference=_cadence_reference(),
    )
    injections, injection_manifest = summarize_fixed_injection_shards([shard])
    independent, independent_manifest = _passing_independent()
    evaluated, bins, target_flags = evaluate_tier1_gates(
        sector=56,
        config=config,
        tier0_summary=_passing_tier0(TEST_COMPACT_SHA256),
        target_metrics=targets,
        aperture_metrics=apertures,
        pair_metrics=pairs,
        injection_metrics=injections,
        injection_manifest=injection_manifest,
        independent_metrics=independent,
        independent_manifest=independent_manifest,
        current_compact_sha256=TEST_COMPACT_SHA256,
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        cadence_reference_gate={"status": "pass", "reasons": []},
        injection_source_parity_gate={"status": "pass", "reasons": []},
    )
    assert evaluated["status"] == "pass"
    assert evaluated["enrichment_ready"] is True
    assert evaluated["science_ready"] is False
    assert len(bins) >= 6
    assert target_flags["relative_cadence_loss"].max() == 0.0
    assert target_flags["tier1_target_qa_pass"].all()

    old_tier0 = _passing_tier0(TEST_COMPACT_SHA256)
    old_tier0["contract_version"] = "a2v1_photometric_qa_v1"
    failed, _, _ = evaluate_tier1_gates(
        sector=56,
        config=config,
        tier0_summary=old_tier0,
        target_metrics=targets,
        aperture_metrics=apertures,
        pair_metrics=pairs,
        injection_metrics=injections,
        injection_manifest=injection_manifest,
        independent_metrics=independent,
        independent_manifest=independent_manifest,
        current_compact_sha256=TEST_COMPACT_SHA256,
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
        cadence_reference_gate={"status": "pass", "reasons": []},
        injection_source_parity_gate={"status": "pass", "reasons": []},
    )
    assert failed["status"] == "fail"
    assert failed["enrichment_ready"] is False


def test_strict_json_replaces_nonfinite_values_atomically(tmp_path) -> None:
    path = tmp_path / "summary.json"
    write_strict_json(path, {"finite": 1.0, "nan": np.nan, "inf": np.inf})
    payload = json.loads(path.read_text())
    assert payload == {"finite": 1.0, "inf": None, "nan": None}
    assert not path.with_suffix(".json.tmp").exists()


def test_end_to_end_runner_writes_a_scoped_pass(tmp_path) -> None:
    compact = tmp_path / "compact.h5"
    _write_compact(compact)
    compact_sha256 = file_sha256(compact)
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    injection_parity_path = tmp_path / "injection_parity.json"
    _write_injection_parity(injection_parity_path, compact_sha256)
    config = _small_injection_config(
        expected_compact_sha256=compact_sha256,
        expected_injection_shard_sha256=(file_sha256(shard),),
        expected_injection_parity_report_sha256=file_sha256(injection_parity_path),
        min_population_targets=40,
        min_targets_per_magnitude_bin=3,
    )
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(asdict(config)))
    tier0_path = tmp_path / "tier0.json"
    tier0_path.write_text(json.dumps(_passing_tier0(compact_sha256)))

    cadence_path = tmp_path / "cadence.csv"
    _, cadence_manifest = _write_cadence_evidence(cadence_path)
    cadence_manifest_path = tmp_path / "cadence.json"
    cadence_manifest_path.write_text(json.dumps(cadence_manifest))
    config = replace(
        config,
        expected_cadence_reference_sha256=file_sha256(cadence_path),
        expected_cadence_reference_manifest_sha256=file_sha256(
            cadence_manifest_path
        ),
    )
    config_path.write_text(json.dumps(asdict(config)))

    independent, independent_manifest = _passing_independent(
        current_compact_sha256=compact_sha256
    )
    independent_path = tmp_path / "independent.csv"
    independent.to_csv(independent_path, index=False)
    reloaded_independent = pd.read_csv(
        independent_path, low_memory=False, float_precision="round_trip"
    )
    independent_manifest["metrics_file_sha256"] = file_sha256(independent_path)
    independent_manifest["metrics_content_sha256"] = _dataframe_content_sha256(
        reloaded_independent
    )
    independent_manifest_path = tmp_path / "independent.json"
    independent_manifest_path.write_text(json.dumps(independent_manifest))
    config = replace(
        config,
        expected_independent_metrics_sha256=file_sha256(independent_path),
        expected_independent_manifest_sha256=file_sha256(
            independent_manifest_path
        ),
        expected_independent_reference_product_sha256=independent_manifest[
            "reference_product_sha256"
        ],
    )
    config_path.write_text(json.dumps(asdict(config)))

    out_dir = tmp_path / "out"
    gate_json = tmp_path / "gate.json"
    summary = run_a2v1_tier1_qa(
        sector=56,
        config_path=config_path,
        tier0_summary_path=tier0_path,
        compact_lc=compact,
        cadence_reference_path=cadence_path,
        cadence_reference_manifest_path=cadence_manifest_path,
        injection_source_parity_path=injection_parity_path,
        injection_shards=[shard],
        independent_metrics_path=independent_path,
        independent_manifest_path=independent_manifest_path,
        out_dir=out_dir,
        gate_json=gate_json,
    )
    assert summary["status"] == "pass"
    assert summary["enrichment_ready"] is True
    assert summary["science_ready"] is False
    assert summary["target_qa"]["n_pass"] == 40
    assert gate_json.exists()
    assert (out_dir / "target_metrics.parquet").exists()
    published_manifest_path = out_dir / "fixed_injection_manifest.json"
    assert published_manifest_path.exists()
    published_manifest = json.loads(published_manifest_path.read_text())
    assert published_manifest["shard_sha256"] == {
        str(shard): file_sha256(shard)
    }
    assert published_manifest["selection_sha256"] == config.expected_injection_selection_sha256
    assert published_manifest["metadata_sha256"] == config.expected_injection_metadata_sha256
    assert published_manifest["metrics_file_sha256"] == file_sha256(
        out_dir / "fixed_injection_metrics.csv"
    )
    assert summary["provenance"]["published_injection_manifest_sha256"] == file_sha256(
        published_manifest_path
    )
    assert json.loads(gate_json.read_text())["passed"] is True

    # The published pair must be reusable without losing the hashes of the
    # precomputed inputs that supplied the second run.
    precomputed_out = tmp_path / "precomputed_out"
    precomputed_gate = tmp_path / "precomputed_gate.json"
    precomputed = run_a2v1_tier1_qa(
        sector=56,
        config_path=config_path,
        tier0_summary_path=tier0_path,
        compact_lc=compact,
        cadence_reference_path=cadence_path,
        cadence_reference_manifest_path=cadence_manifest_path,
        injection_source_parity_path=injection_parity_path,
        injection_metrics_path=out_dir / "fixed_injection_metrics.csv",
        injection_manifest_path=published_manifest_path,
        independent_metrics_path=independent_path,
        independent_manifest_path=independent_manifest_path,
        out_dir=precomputed_out,
        gate_json=precomputed_gate,
    )
    assert precomputed["passed"] is True
    assert precomputed["provenance"]["injection_metrics_input_sha256"] == file_sha256(
        out_dir / "fixed_injection_metrics.csv"
    )
    assert precomputed["provenance"]["injection_manifest_input_sha256"] == file_sha256(
        published_manifest_path
    )


def test_runner_rejects_input_output_collision(tmp_path) -> None:
    compact = tmp_path / "compact.h5"
    _write_compact(compact)
    config = _small_injection_config(expected_compact_sha256=file_sha256(compact))
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(asdict(config)))
    with pytest.raises(ValueError, match="collide"):
        run_a2v1_tier1_qa(
            sector=56,
            config_path=config_path,
            tier0_summary_path=config_path,
            compact_lc=compact,
            cadence_reference_path=config_path,
            cadence_reference_manifest_path=config_path,
            injection_source_parity_path=config_path,
            injection_shards=[config_path],
            independent_metrics_path=config_path,
            independent_manifest_path=config_path,
            out_dir=tmp_path / "out",
            gate_json=config_path,
        )

    with pytest.raises(ValueError, match="output paths collide"):
        run_a2v1_tier1_qa(
            sector=56,
            config_path=config_path,
            tier0_summary_path=tmp_path / "tier0.json",
            compact_lc=compact,
            cadence_reference_path=tmp_path / "cadence.csv",
            cadence_reference_manifest_path=tmp_path / "cadence.json",
            injection_source_parity_path=tmp_path / "injection_parity.json",
            injection_shards=[tmp_path / "injections.h5"],
            independent_metrics_path=tmp_path / "independent.csv",
            independent_manifest_path=tmp_path / "independent.json",
            out_dir=tmp_path / "out",
            gate_json=tmp_path / "out" / "target_metrics.parquet",
        )
