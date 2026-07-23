from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json

import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.lightcurves.a2v1_tier1_qa as tier1_qa
from twirl.injections.a2v1_recovery import epoch_quality_provenance
from twirl.lightcurves.a2v1_cadence_reference import (
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    AUTHORITY_EXCLUSION_POLICY,
    AUTHORITY_EXCLUSION_POLICY_CONTRACT,
    CADENCE_REFERENCE_BUILDER_VERSION,
    authority_exclusions_sha256,
)
from twirl.lightcurves.a2v1_qa import (
    A2V1_PHOTOMETRIC_QA_VERSION,
    WD1856_PERIOD_D,
    WD1856_T0_BJD,
    WD1856_TIC,
    file_sha256,
)
from twirl.lightcurves.a2v1_tier1_qa import (
    TIER1_QA_CONTRACT_VERSION,
    TIER1_TARGET_REASON_CODES,
    Tier1QAConfig,
    _dataframe_content_sha256,
    attach_target_qa_flags,
    audit_compact_population,
    build_target_eligibility,
    evaluate_aperture_outliers,
    evaluate_cadence_quality,
    evaluate_fixed_injections,
    evaluate_independent_extraction,
    evaluate_tier0_prerequisite,
    evaluate_tier1_gates,
    injection_metadata_sha256,
    load_tier1_config,
    run_a2v1_tier1_qa,
    summarize_fixed_injection_shards,
    write_strict_json,
)
from twirl.lightcurves.external_quality import (
    EFFECTIVE_QUALITY_POLICY,
    EXTERNAL_QUALITY_POLICY_CONTRACT,
    load_external_quality_reference,
)
from twirl.vetting.adp_only import ADP_ONLY_APERTURES


TEST_COMPACT_SHA256 = "a" * 64
TEST_TIER0_SUMMARY_SHA256 = "6" * 64
TEST_TIER0_BLS_SHA256 = "7" * 64
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


def _write_compact(path, *, n_targets: int = 40, n_cadences: int = 202) -> None:
    rng = np.random.default_rng(56)
    orbit_boundary = min(100, n_cadences // 2)
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
            if index == 0:
                group.attrs["gaia_dr3_source_id"] = "2098419251571450880"
            group.create_dataset("time", data=np.arange(n_cadences) / 432.0)
            group.create_dataset("cadenceno", data=np.arange(n_cadences))
            quality = np.zeros(n_cadences, dtype=np.int32)
            quality[:2] = 1
            group.create_dataset("quality", data=quality)
            group.create_dataset(
                "orbitid",
                data=np.where(np.arange(n_cadences) < orbit_boundary, 119, 120),
            )
            sigma = 0.002 * 10 ** (0.12 * (float(group.attrs["tessmag"]) - 17.0))
            base = 1.0 + rng.normal(0.0, sigma, n_cadences)
            group.create_dataset(ADP_ONLY_APERTURES[0], data=base)
            group.create_dataset(
                ADP_ONLY_APERTURES[1], data=1.0 + 1.05 * (base - 1.0)
            )


def _cadence_reference(*, n_cadences: int = 202) -> pd.DataFrame:
    rows: list[dict] = []
    orbit_boundary = min(100, n_cadences // 2)
    for camera in range(1, 5):
        for ccd in range(1, 5):
            for cadence in range(n_cadences):
                rows.append(
                    {
                        "sector": 56,
                        "orbitid": 119 if cadence < orbit_boundary else 120,
                        "camera": camera,
                        "ccd": ccd,
                        "cadenceno": cadence,
                        "spoc_quality": 1 if cadence < 2 else 0,
                        "qlp_quality": 0,
                        "external_quality": 1 if cadence < 2 else 0,
                    }
                )
    return pd.DataFrame(rows)


def _write_cadence_evidence(path) -> tuple[pd.DataFrame, dict]:
    table = _cadence_reference()
    table.to_csv(path, index=False)
    sources: list[dict] = []
    source_hashes: dict[str, str] = {}

    def source(role: str, **metadata: int) -> None:
        index = len(sources) + 1
        source_path = str((path.parent / f"authority_{index:03d}.dat").resolve())
        digest = f"{index:064x}"
        sources.append(
            {"role": role, "path": source_path, "sha256": digest, **metadata}
        )
        source_hashes[source_path] = digest

    source("spoc_quality_table")
    source("spoc_quality_provenance")
    for camera in range(1, 5):
        for ccd in range(1, 5):
            source("spoc_flag_file", camera=camera, ccd=ccd)
    for orbitid in (119, 120):
        for camera in range(1, 5):
            source("qlp_cam_quat", orbitid=orbitid, camera=camera)
            for ccd in range(1, 5):
                source(
                    "qlp_detector_qflag",
                    orbitid=orbitid,
                    camera=camera,
                    ccd=ccd,
                )
    detectors = [
        f"cam{camera}_ccd{ccd}"
        for camera in range(1, 5)
        for ccd in range(1, 5)
    ]
    authority_exclusions = {
        "contract_version": AUTHORITY_EXCLUSION_POLICY_CONTRACT,
        "policy": AUTHORITY_EXCLUSION_POLICY,
        "external_bit": AUTHORITY_EXCLUSION_EXTERNAL_BIT,
        "n_rows": 0,
        "by_detector": {
            detector: {"n_rows": 0, "rows": []} for detector in detectors
        },
    }
    manifest = {
        "contract_version": "s56_a2v1_cadence_reference_v1",
        "builder_version": CADENCE_REFERENCE_BUILDER_VERSION,
        "sector": 56,
        "cadence_authority": "qlp_cam_quat",
        "quality_authority": "spoc_and_qlp_quality_flags",
        "quality_composition": {
            "external_quality": "spoc_quality | (qlp_quality << 30)",
            "qlp_quality_raw_values": [0, 1],
            "qlp_quality_external_bit": 30,
        },
        "table_sha256": file_sha256(path),
        "table_columns": list(table.columns),
        "n_rows": len(table),
        "detectors": detectors,
        "orbits": [119, 120],
        "n_rows_by_detector": {
            detector: int(len(table) / len(detectors)) for detector in detectors
        },
        "source_file_sha256": source_hashes,
        "sources": sources,
        "n_spoc_authority_files_verified": 16,
        "n_qlp_qflag_files_verified": 32,
        "n_nonzero_spoc_quality": int((table["spoc_quality"] != 0).sum()),
        "n_nonzero_qlp_quality": int((table["qlp_quality"] != 0).sum()),
        "n_nonzero_external_quality": int((table["external_quality"] != 0).sum()),
        "n_spoc_rows_excluded_by_quat": 0,
        "authority_exclusions": authority_exclusions,
        "authority_exclusions_sha256": authority_exclusions_sha256(
            authority_exclusions
        ),
    }
    return table, manifest


def _write_injection_shard(
    path,
    *,
    n_ids: int = 8,
    retention: float = 1.0,
    shard_index: int = 0,
) -> None:
    rng = np.random.default_rng(123 + shard_index)
    with h5py.File(path, "w") as h5:
        h5.attrs["contract_version"] = "s56_a2v1_fresh_injection_pair_v2"
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
            group.attrs["sector"] = 56
            group.attrs["period_d"] = periods[index]
            group.attrs["duration_min"] = durations[index]
            group.attrs["model_depth"] = depths[index]
            group.attrs["radius_rearth"] = 0.2 + index
            group.attrs["injection_baseline_Small"] = 120.0
            group.attrs["DET_FLUX_ADP_SML_scale"] = 240.0
            group.attrs["injection_baseline_Primary"] = 150.0
            group.attrs["DET_FLUX_ADP_scale"] = 500.0
            n = 200
            model = np.ones(n)
            model[40:45] = np.array([0.97, 0.93, 0.90, 0.93, 0.97])
            group.create_dataset("quality", data=np.zeros(n, dtype=np.int32))
            group.create_dataset("cadenceno", data=np.arange(n, dtype=np.int64))
            group.create_dataset(
                "orbitid",
                data=np.where(np.arange(n) < n // 2, 119, 120),
            )
            group.create_dataset("transit_model", data=model)
            for aperture in ADP_ONLY_APERTURES:
                baseline_over_scale = (
                    float(group.attrs["injection_baseline_Small"])
                    / float(group.attrs["DET_FLUX_ADP_SML_scale"])
                    if aperture == ADP_ONLY_APERTURES[0]
                    else float(group.attrs["injection_baseline_Primary"])
                    / float(group.attrs["DET_FLUX_ADP_scale"])
                )
                original = 1.0 + rng.normal(0.0, 1.0e-4, n)
                injected = (
                    original
                    - retention * baseline_over_scale * (1.0 - model)
                    + 2.0e-4
                )
                group.create_dataset(f"{aperture}_original", data=original)
                group.create_dataset(f"{aperture}_injected", data=injected)


def _bind_injection_shard_to_cadence_reference(
    shard_path,
    cadence_path,
    cadence_manifest_path,
) -> None:
    reference = load_external_quality_reference(
        table_path=cadence_path,
        manifest_path=cadence_manifest_path,
        sector=56,
    )
    provenance = epoch_quality_provenance(reference)
    with h5py.File(shard_path, "r+") as h5:
        for name, value in provenance.items():
            h5.attrs[name] = value
        for injection_id, group in h5["injections"].items():
            overlay = reference.apply(
                sector=int(group.attrs["sector"]),
                camera=int(group.attrs["camera"]),
                ccd=int(group.attrs["ccd"]),
                cadenceno=np.asarray(group["cadenceno"]),
                orbitid=np.asarray(group["orbitid"]),
                internal_quality=np.asarray(group["quality"]),
                context=f"test injection {injection_id}",
            )
            for name in ("external_quality", "effective_quality"):
                if name in group:
                    del group[name]
            group.create_dataset(
                "external_quality",
                data=np.asarray(overlay.external_quality, dtype=np.int64),
            )
            group.create_dataset(
                "effective_quality",
                data=np.asarray(overlay.quality, dtype=np.int32),
            )
            for name, value in provenance.items():
                group.attrs[name] = value
            for name, value in overlay.counts.items():
                group.attrs[f"epoch_quality_{name}"] = int(value)
    reference.assert_unchanged()


def _small_injection_config(**changes) -> Tier1QAConfig:
    base = {
        "expected_compact_sha256": TEST_COMPACT_SHA256,
        "expected_tier0_summary_sha256": TEST_TIER0_SUMMARY_SHA256,
        "expected_tier0_bls_peaks_sha256": TEST_TIER0_BLS_SHA256,
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
    *,
    current_compact_sha256: str = TEST_COMPACT_SHA256,
    cadence_reference_sha256: str = "0" * 64,
    cadence_reference_manifest_sha256: str = "0" * 64,
) -> tuple[pd.DataFrame, dict]:
    metrics = pd.DataFrame(
        {
            "tic": [WD1856_TIC, WD1856_TIC],
            "sector": [56, 56],
            "aperture": list(ADP_ONLY_APERTURES),
            "reference_flux_column": ["REFERENCE_FLUX_SML", "REFERENCE_FLUX"],
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
        "contract_version": "s56_a2v1_independent_extraction_v2",
        "sector": 56,
        "tic": WD1856_TIC,
        "independent": True,
        "comparison_mode": "signal_timing_only",
        "independence_basis": "independent pixels and extraction code",
        "current_extractor_family": "MIT TGLC A2v1",
        "reference_extractor_family": "MAST QLP",
        "current_repository": "TeHanHunter/TESS_Gaia_Light_Curve",
        "reference_repository": "mit-qlp",
        "current_revision": "abc123",
        "reference_revision": "release-1",
        "pixel_source": "independently calibrated TICA cutout",
        "cadence_match_policy": (
            "exact common CADENCENO intersection with effective quality from "
            "native internal OR authoritative external flags"
        ),
        "scatter_definition": "out-of-transit robust MAD in ppm on common cadences",
        "depth_definition": "median in-event deficit on a fixed ephemeris",
        "ephemeris_recovery_definition": (
            "independent bounded BoxLeastSquares recovery on matched cadences"
        ),
        "current_apertures": list(ADP_ONLY_APERTURES),
        "reference_flux_columns": {
            ADP_ONLY_APERTURES[0]: "REFERENCE_FLUX_SML",
            ADP_ONLY_APERTURES[1]: "REFERENCE_FLUX",
        },
        "reference_columns": {
            "flux_by_current_aperture": {
                ADP_ONLY_APERTURES[0]: "REFERENCE_FLUX_SML",
                ADP_ONLY_APERTURES[1]: "REFERENCE_FLUX",
            }
        },
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
        "n_current_cadences": 11_775,
        "n_reference_cadences": 11_770,
        "n_common_cadences": 11_760,
        "external_quality_overlay": {
            "policy_contract": EXTERNAL_QUALITY_POLICY_CONTRACT,
            "effective_quality_policy": EFFECTIVE_QUALITY_POLICY,
            "sector": 56,
            "cadence_reference_contract_version": (
                "s56_a2v1_cadence_reference_v1"
            ),
            "cadence_reference_cadence_authority": "qlp_cam_quat",
            "cadence_reference_quality_authority": (
                "spoc_and_qlp_quality_flags"
            ),
            "cadence_reference_table_sha256": cadence_reference_sha256,
            "cadence_reference_manifest_sha256": (
                cadence_reference_manifest_sha256
            ),
            "cadence_reference_source_declaration_sha256": "1" * 64,
            "applied_before_common_cadence_metrics": True,
            "current_compact_full_audit_counts": {
                "n_cad_total": 11_775,
                "n_cad_internal_bad": 100,
                "n_cad_external_bad": 80,
                "n_cad_authority_excluded": 0,
                "n_cad_external_only_bad": 30,
                "n_cad_effective_bad": 130,
            },
            "current_common_audit_counts": {
                "n_cad_total": 11_760,
                "n_cad_internal_bad": 99,
                "n_cad_external_bad": 79,
                "n_cad_authority_excluded": 0,
                "n_cad_external_only_bad": 29,
                "n_cad_effective_bad": 128,
            },
            "reference_common_audit_counts": {
                "n_cad_total": 11_760,
                "n_cad_internal_bad": 50,
                "n_cad_external_bad": 79,
                "n_cad_authority_excluded": 0,
                "n_cad_external_only_bad": 40,
                "n_cad_effective_bad": 90,
            },
        },
        "metrics_file_sha256": "f" * 64,
        "metrics_content_sha256": _dataframe_content_sha256(metrics),
    }
    return metrics, manifest


def _passing_tier0(compact_sha256: str) -> dict:
    gates = {
        name: {"passed": True}
        for name in (
            "schema_and_completeness",
            "bls_target_aperture_coverage",
            "sampled_adp_photometry",
            "aperture_consistency",
            "prior_extraction_tree_comparison",
            "benchmark",
        )
    }
    return {
        "sector": 56,
        "passed": True,
        "contract_version": A2V1_PHOTOMETRIC_QA_VERSION,
        "qa_tier": "tier0_integrity_and_benchmark",
        "science_ready": False,
        "gates": gates,
        "benchmarks": {"wd1856": {"passed": True}},
        "provenance": {
            "compact_lc_sha256": compact_sha256,
            "schema_summary_sha256": "b" * 64,
            "bls_peaks_sha256": TEST_TIER0_BLS_SHA256,
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
    with pytest.raises(ValueError, match="unsupported Tier-1 contract_version"):
        replace(
            Tier1QAConfig(),
            contract_version="a2v1_tier1_science_qa_v1",
        ).validate()
    with pytest.raises(ValueError, match="locked A2v1 BLS minimum"):
        replace(Tier1QAConfig(), min_searchable_cadences=199).validate()


def test_tier0_prerequisite_is_hash_bound_and_checks_every_nested_gate() -> None:
    config = _small_injection_config()
    summary = _passing_tier0(TEST_COMPACT_SHA256)
    passed = evaluate_tier0_prerequisite(
        summary,
        sector=56,
        config=config,
        tier0_summary_sha256=TEST_TIER0_SUMMARY_SHA256,
        current_compact_sha256=TEST_COMPACT_SHA256,
    )
    assert passed["status"] == "pass"

    stale = _passing_tier0(TEST_COMPACT_SHA256)
    stale["gates"]["benchmark"]["passed"] = False
    stale["provenance"]["schema_summary_sha256"] = "not-even-a-digest"
    failed = evaluate_tier0_prerequisite(
        stale,
        sector=56,
        config=config,
        tier0_summary_sha256="5" * 64,
        current_compact_sha256=TEST_COMPACT_SHA256,
    )
    assert failed["status"] == "fail"
    assert "Tier-0 summary is not pinned by the locked configuration" in failed[
        "reasons"
    ]
    assert "Tier-0 gate benchmark did not pass" in failed["reasons"]
    assert "Tier-0 schema evidence hash is invalid" in failed["reasons"]


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
    assert targets["max_orbit_missing_cadence_fraction"].max() == 0.0
    assert targets["missing_orbit_ids"].eq("").all()
    assert targets.loc[targets["tic"].eq(WD1856_TIC), "gaia_dr3_source_id"].iloc[0] == (
        "2098419251571450880"
    )
    assert targets["n_quality_mask_disagreements"].sum() == 0
    assert apertures["finite_quality0_fraction"].min() == 1.0
    assert pairs["mad_ratio"].between(0.9, 1.2).all()

    external_mask = _cadence_reference()
    externally_bad = external_mask["cadenceno"].eq(3)
    external_mask.loc[externally_bad, "spoc_quality"] = 32
    external_mask.loc[externally_bad, "external_quality"] = 32
    masked_targets, masked_apertures, _ = audit_compact_population(
        compact,
        apertures=ADP_ONLY_APERTURES,
        cadence_reference=external_mask,
    )
    assert masked_targets["n_quality_mask_disagreements"].eq(1).all()
    assert masked_targets["internal_quality0_fraction"].eq(200 / 202).all()
    assert masked_targets["external_quality0_fraction"].eq(199 / 202).all()
    assert masked_targets["quality0_fraction"].eq(199 / 202).all()
    cadence_gate, _ = evaluate_cadence_quality(
        masked_targets, masked_apertures, Tier1QAConfig()
    )
    assert cadence_gate["status"] == "pass"


def test_compact_population_masks_only_declared_authority_exclusions(
    tmp_path,
) -> None:
    compact = tmp_path / "compact.h5"
    _write_compact(compact)
    cadence_reference = _cadence_reference()
    cadence_reference = cadence_reference.loc[
        ~cadence_reference["cadenceno"].eq(3)
    ].reset_index(drop=True)
    exclusions = {
        (56, camera, ccd): frozenset({3})
        for camera in range(1, 5)
        for ccd in range(1, 5)
    }

    targets, apertures, _ = audit_compact_population(
        compact,
        apertures=ADP_ONLY_APERTURES,
        cadence_reference=cadence_reference,
        authority_exclusions=exclusions,
    )

    assert targets["n_authority_excluded_cadences"].eq(1).all()
    assert targets["n_unexpected_cadences"].eq(0).all()
    assert targets["quality0_fraction"].eq(199 / 202).all()
    gate, _ = evaluate_cadence_quality(targets, apertures, Tier1QAConfig())
    assert gate["status"] == "pass"
    assert gate["authority_excluded_cadences_total"] == len(targets)

    unknown_reference = cadence_reference.loc[
        ~cadence_reference["cadenceno"].eq(4)
    ].reset_index(drop=True)
    unknown_targets, unknown_apertures, _ = audit_compact_population(
        compact,
        apertures=ADP_ONLY_APERTURES,
        cadence_reference=unknown_reference,
        authority_exclusions=exclusions,
    )
    assert unknown_targets["n_authority_excluded_cadences"].eq(1).all()
    assert unknown_targets["n_unexpected_cadences"].eq(1).all()
    unknown_gate, _ = evaluate_cadence_quality(
        unknown_targets, unknown_apertures, Tier1QAConfig()
    )
    assert unknown_gate["status"] == "fail"


def test_flagged_fraction_is_diagnostic_until_bls_cadence_floor(tmp_path) -> None:
    searchable_compact = tmp_path / "searchable.h5"
    _write_compact(searchable_compact, n_targets=2, n_cadences=500)
    with h5py.File(searchable_compact, "r+") as handle:
        for group in handle["targets"].values():
            quality = np.zeros(500, dtype=np.int32)
            quality[:250] = 1
            group["quality"][:] = quality
    targets, apertures, _ = audit_compact_population(
        searchable_compact,
        apertures=ADP_ONLY_APERTURES,
        cadence_reference=_cadence_reference(n_cadences=500),
    )
    gate, flagged = evaluate_cadence_quality(
        targets, apertures, Tier1QAConfig()
    )
    assert gate["status"] == "review"
    assert flagged["n_finite_quality0_min"].eq(250).all()
    assert flagged["usable_cadence_fraction_min"].eq(0.5).all()
    assert flagged["tier1_target_searchable"].all()
    assert flagged["tier1_target_searchability_reasons"].eq("").all()

    excluded_compact = tmp_path / "excluded.h5"
    _write_compact(excluded_compact, n_targets=2, n_cadences=500)
    with h5py.File(excluded_compact, "r+") as handle:
        for group in handle["targets"].values():
            quality = np.zeros(500, dtype=np.int32)
            quality[:301] = 1
            group["quality"][:] = quality
    targets, apertures, _ = audit_compact_population(
        excluded_compact,
        apertures=ADP_ONLY_APERTURES,
        cadence_reference=_cadence_reference(n_cadences=500),
    )
    _, flagged = evaluate_cadence_quality(targets, apertures, Tier1QAConfig())
    assert flagged["n_finite_quality0_min"].eq(199).all()
    assert not flagged["tier1_target_searchable"].any()
    assert flagged["tier1_target_searchability_reasons"].eq(
        "too_few_effective_cadences"
    ).all()


def test_model_weighted_injection_retention_and_full_shard_gate(tmp_path) -> None:
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    cadence_path = tmp_path / "cadence.csv"
    _, cadence_manifest = _write_cadence_evidence(cadence_path)
    cadence_manifest_path = tmp_path / "cadence.json"
    cadence_manifest_path.write_text(json.dumps(cadence_manifest) + "\n")
    _bind_injection_shard_to_cadence_reference(
        shard,
        cadence_path,
        cadence_manifest_path,
    )
    metrics, manifest = summarize_fixed_injection_shards(
        [shard],
        sector=56,
        cadence_reference_path=cadence_path,
        cadence_reference_manifest_path=cadence_manifest_path,
    )
    assert len(metrics) == 16
    assert metrics["status"].eq("ok").all()
    assert np.allclose(metrics["depth_retention_fraction"], 1.0, atol=1.0e-8)
    assert set(metrics["baseline_over_scale"].round(8)) == {0.3, 0.5}
    assert np.allclose(
        metrics["raw_depth_response_slope"],
        metrics["baseline_over_scale"],
        atol=1.0e-8,
    )
    assert (
        manifest["retention_metric"]["version"]
        == tier1_qa.INJECTION_RETENTION_METRIC_VERSION
    )

    config = _small_injection_config(
        expected_injection_shard_sha256=(file_sha256(shard),),
        expected_cadence_reference_sha256=file_sha256(cadence_path),
        expected_cadence_reference_manifest_sha256=file_sha256(
            cadence_manifest_path
        ),
    )
    gate = evaluate_fixed_injections(metrics, manifest, sector=56, config=config)
    assert gate["status"] == "pass"
    assert gate["n_faint_injection_ids"] == 4
    assert gate["period_support_ratio"] >= 5.0
    assert all(
        np.isclose(aperture["median_depth_retention"], 1.0)
        for aperture in gate["apertures"]
    )

    internal_only = dict(manifest)
    internal_only.pop("external_quality_overlay")
    rejected = evaluate_fixed_injections(
        metrics, internal_only, sector=56, config=config
    )
    assert rejected["status"] == "fail"
    assert any("overlay is missing" in value for value in rejected["reasons"])

    bad_counts = json.loads(json.dumps(manifest))
    bad_counts["external_quality_overlay"]["aggregate_counts"][
        "n_cad_effective_bad"
    ] += 1
    rejected = evaluate_fixed_injections(
        metrics, bad_counts, sector=56, config=config
    )
    assert rejected["status"] == "fail"
    assert any("arithmetic" in value for value in rejected["reasons"])

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


@pytest.mark.parametrize(
    ("attribute", "replacement"),
    (
        ("DET_FLUX_ADP_SML_scale", None),
        ("DET_FLUX_ADP_SML_scale", 0.0),
        ("injection_baseline_Small", np.nan),
    ),
)
def test_injection_retention_rejects_missing_or_invalid_normalization(
    tmp_path,
    attribute: str,
    replacement: float | None,
) -> None:
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    with h5py.File(shard, "r+") as h5:
        group = h5["injections/inj_00_000"]
        if replacement is None:
            del group.attrs[attribute]
        else:
            group.attrs[attribute] = replacement

    cadence_path = tmp_path / "cadence.csv"
    _, cadence_manifest = _write_cadence_evidence(cadence_path)
    cadence_manifest_path = tmp_path / "cadence.json"
    cadence_manifest_path.write_text(json.dumps(cadence_manifest) + "\n")
    _bind_injection_shard_to_cadence_reference(
        shard,
        cadence_path,
        cadence_manifest_path,
    )
    metrics, manifest = summarize_fixed_injection_shards(
        [shard],
        sector=56,
        cadence_reference_path=cadence_path,
        cadence_reference_manifest_path=cadence_manifest_path,
    )

    affected = metrics.loc[
        metrics["injection_id"].eq("inj_00_000")
        & metrics["aperture"].eq(ADP_ONLY_APERTURES[0])
    ].iloc[0]
    assert affected["status"] == "invalid_normalization_metadata"
    assert not np.isfinite(affected["baseline_over_scale"])
    assert manifest["malformed_shards"]

    config = _small_injection_config(
        expected_injection_shard_sha256=(file_sha256(shard),),
        expected_cadence_reference_sha256=file_sha256(cadence_path),
        expected_cadence_reference_manifest_sha256=file_sha256(
            cadence_manifest_path
        ),
    )
    gate = evaluate_fixed_injections(
        metrics,
        manifest,
        sector=56,
        config=config,
    )
    assert gate["status"] == "fail"
    assert any("malformed" in reason for reason in gate["reasons"])


def test_injection_retention_rejects_stale_epoch_quality_provenance_and_masks(
    tmp_path,
) -> None:
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    cadence_path = tmp_path / "cadence.csv"
    _, cadence_manifest = _write_cadence_evidence(cadence_path)
    cadence_manifest_path = tmp_path / "cadence.json"
    cadence_manifest_path.write_text(json.dumps(cadence_manifest) + "\n")
    _bind_injection_shard_to_cadence_reference(
        shard,
        cadence_path,
        cadence_manifest_path,
    )

    with h5py.File(shard, "r+") as h5:
        h5.attrs["cadence_reference_manifest_sha256"] = "f" * 64
        group = h5["injections/inj_00_000"]
        group["external_quality"][10] = 2048
        group["effective_quality"][10] = 1

    metrics, manifest = summarize_fixed_injection_shards(
        [shard],
        sector=56,
        cadence_reference_path=cadence_path,
        cadence_reference_manifest_path=cadence_manifest_path,
    )
    assert manifest["stored_epoch_quality_validation"]["status"] == "fail"
    assert (
        manifest["stored_epoch_quality_validation"]["n_provenance_failures"]
        >= 1
    )
    assert manifest["stored_epoch_quality_validation"]["n_mask_failures"] >= 1
    affected = metrics.loc[
        metrics["injection_id"].eq("inj_00_000")
    ]
    assert affected["status"].eq("quality_overlay_error").all()
    assert manifest["malformed_shards"]

    config = _small_injection_config(
        expected_injection_shard_sha256=(file_sha256(shard),),
        expected_cadence_reference_sha256=file_sha256(cadence_path),
        expected_cadence_reference_manifest_sha256=file_sha256(
            cadence_manifest_path
        ),
    )
    gate = evaluate_fixed_injections(
        metrics,
        manifest,
        sector=56,
        config=config,
    )
    assert gate["status"] == "fail"
    assert any(
        "epoch-quality" in reason or "malformed" in reason
        for reason in gate["reasons"]
    )


def test_injection_gate_propagates_nonzero_authority_exclusion(tmp_path) -> None:
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    cadence_path = tmp_path / "cadence.csv"
    table, cadence_manifest = _write_cadence_evidence(cadence_path)
    excluded_cadence = 10
    table = table.loc[
        ~(
            table["camera"].eq(1)
            & table["ccd"].eq(1)
            & table["cadenceno"].eq(excluded_cadence)
        )
    ].reset_index(drop=True)
    table.to_csv(cadence_path, index=False)
    exclusion = cadence_manifest["authority_exclusions"]
    exclusion["n_rows"] = 1
    exclusion["by_detector"]["cam1_ccd1"] = {
        "n_rows": 1,
        "rows": [
            {
                "cadenceno": excluded_cadence,
                "spoc_quality": 0,
            }
        ],
    }
    cadence_manifest["n_spoc_rows_excluded_by_quat"] = 1
    cadence_manifest["authority_exclusions_sha256"] = (
        authority_exclusions_sha256(exclusion)
    )
    cadence_manifest["table_sha256"] = file_sha256(cadence_path)
    cadence_manifest["n_rows"] = len(table)
    cadence_manifest["n_rows_by_detector"]["cam1_ccd1"] -= 1
    cadence_manifest_path = tmp_path / "cadence.json"
    cadence_manifest_path.write_text(
        json.dumps(cadence_manifest) + "\n",
        encoding="utf-8",
    )
    _bind_injection_shard_to_cadence_reference(
        shard,
        cadence_path,
        cadence_manifest_path,
    )

    metrics, manifest = summarize_fixed_injection_shards(
        [shard],
        sector=56,
        cadence_reference_path=cadence_path,
        cadence_reference_manifest_path=cadence_manifest_path,
    )

    overlay = manifest["external_quality_overlay"]
    assert overlay["aggregate_counts"]["n_cad_authority_excluded"] == 1
    assert overlay["aggregate_counts"]["n_cad_external_bad"] == 17
    assert overlay["aggregate_counts"]["n_cad_effective_bad"] == 17
    shard_counts = next(iter(overlay["counts_by_shard"].values()))
    assert shard_counts["n_cad_authority_excluded"] == 1
    config = _small_injection_config(
        expected_injection_shard_sha256=(file_sha256(shard),),
        expected_cadence_reference_sha256=file_sha256(cadence_path),
        expected_cadence_reference_manifest_sha256=file_sha256(
            cadence_manifest_path
        ),
    )
    gate = evaluate_fixed_injections(
        metrics,
        manifest,
        sector=56,
        config=config,
    )
    assert gate["status"] == "pass"


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

    # Raw-pixel and decontaminated extractions can have very different
    # dilution and noise.  Those ratios are retained as diagnostics, while
    # signal presence and ephemeris timing define this bounded gate.
    diluted = metrics.copy()
    diluted["reference_depth"] = [0.015, 0.0048]
    diluted["reference_scatter_ppm"] = [8_759.0, 4_090.0]
    diluted_manifest = dict(manifest)
    diluted_manifest["metrics_content_sha256"] = _dataframe_content_sha256(
        diluted
    )
    diluted_gate = evaluate_independent_extraction(
        diluted,
        diluted_manifest,
        _small_injection_config(),
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
    )
    assert diluted_gate["status"] == "pass"
    assert diluted_gate["scatter_ratio_is_diagnostic_only"] is True
    assert (
        diluted_gate["wd1856"]["apertures"][0][
            "depth_ratio_is_diagnostic_only"
        ]
        is True
    )

    truncated_reference = metrics.copy()
    truncated_reference["n_reference_cadences"] = 100
    truncated_reference["n_common_cadences"] = 100
    truncated_reference["n_common_quality0_finite"] = 90
    truncated_reference["n_in_event_cadences"] = 10
    truncated_reference["n_out_of_event_cadences"] = 80
    truncated_manifest = dict(manifest)
    truncated_manifest["metrics_content_sha256"] = _dataframe_content_sha256(
        truncated_reference
    )
    truncated_gate = evaluate_independent_extraction(
        truncated_reference,
        truncated_manifest,
        _small_injection_config(),
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
    )
    assert truncated_gate["status"] == "fail"
    assert truncated_gate["current_cadence_coverage_fraction_min"] < 0.01

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

    bad_quality = json.loads(json.dumps(manifest))
    bad_quality["external_quality_overlay"][
        "current_common_audit_counts"
    ]["n_cad_effective_bad"] += 1
    bad_quality_gate = evaluate_independent_extraction(
        metrics,
        bad_quality,
        _small_injection_config(),
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
    )
    assert bad_quality_gate["status"] == "fail"
    assert any(
        "quality audit current_common_audit_counts is inconsistent" in reason
        for reason in bad_quality_gate["reasons"]
    )

    loose_timing = dict(manifest)
    loose_timing["max_time_delta_seconds_allowed"] = 60.1
    loose_gate = evaluate_independent_extraction(
        metrics,
        loose_timing,
        _small_injection_config(),
        catalog_detectors={WD1856_TIC: "cam4_ccd1"},
    )
    assert loose_gate["status"] == "fail"
    assert any("tolerance" in reason for reason in loose_gate["reasons"])

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


def test_target_aperture_correlation_has_review_band() -> None:
    targets = pd.DataFrame(
        {
            "sector": [56],
            "tic": [123],
            "camera": [1],
            "ccd": [1],
            "detector": ["cam1_ccd1"],
            "tmag": [18.0],
            "tier1_cadence_qa_status": ["pass"],
            "tier1_cadence_qa_reasons": [""],
            "tier1_target_searchable": [True],
            "tier1_target_searchability_reasons": [""],
        }
    )
    apertures = pd.DataFrame(
        {
            "tic": [123, 123],
            "tmag": [18.0, 18.0],
            "aperture": list(ADP_ONLY_APERTURES),
            "mad_ppm": [10_000.0, 11_000.0],
            "rms5_ppm": [15_000.0, 16_000.0],
        }
    )
    pairs = pd.DataFrame({"tic": [123], "mad_ratio": [1.1], "correlation": [-0.1]})
    reviewed = attach_target_qa_flags(targets, apertures, pairs, Tier1QAConfig())
    assert reviewed.loc[0, "tier1_target_qa_status"] == "review"
    assert reviewed.loc[0, "tier1_target_qa_reasons"] == "aperture_correlation_low"

    pairs.loc[0, "correlation"] = -0.3
    failed = attach_target_qa_flags(targets, apertures, pairs, Tier1QAConfig())
    assert failed.loc[0, "tier1_target_qa_status"] == "fail"
    assert failed.loc[0, "tier1_target_qa_reasons"] == "aperture_anticorrelation"


def test_target_eligibility_rejects_coerced_or_inconsistent_pass_flags() -> None:
    target = pd.DataFrame(
        {
            "sector": [56],
            "tic": [123],
            "camera": [1],
            "ccd": [1],
            "detector": ["cam1_ccd1"],
            "tmag": [18.0],
            "tier1_cadence_qa_status": ["pass"],
            "tier1_scatter_qa_status": ["pass"],
            "tier1_aperture_pair_qa_status": ["pass"],
            "tier1_target_qa_status": ["fail"],
            "tier1_target_qa_reasons": ["finite_flux"],
            "tier1_target_qa_pass": ["False"],
            "tier1_target_searchable": [False],
            "tier1_target_searchability_reasons": [
                "too_few_effective_cadences"
            ],
            "n_finite_quality0_min": [0],
            "usable_cadence_fraction_min": [0.0],
        }
    )
    with pytest.raises(ValueError, match="non-null booleans"):
        build_target_eligibility(target, Tier1QAConfig())

    target["tier1_target_qa_pass"] = True
    with pytest.raises(ValueError, match="status and pass flag are inconsistent"):
        build_target_eligibility(target, Tier1QAConfig())


def test_missing_orbit_gets_explicit_target_reason(tmp_path) -> None:
    compact = tmp_path / "compact.h5"
    _write_compact(compact)
    with h5py.File(compact, "r+") as handle:
        group = handle[f"targets/{1001:016d}"]
        group["orbitid"][:] = 119
    targets, apertures, _ = audit_compact_population(
        compact,
        apertures=ADP_ONLY_APERTURES,
        cadence_reference=_cadence_reference(),
    )
    gate, flagged = evaluate_cadence_quality(targets, apertures, Tier1QAConfig())
    row = flagged.loc[flagged["tic"].eq(1001)].iloc[0]
    assert row["missing_orbit_ids"] == "120"
    assert row["max_orbit_missing_cadence_fraction"] == 1.0
    reasons = set(str(row["tier1_cadence_qa_reasons"]).split(";"))
    assert {"orbit_missing", "orbit_cadence_loss", "orbit_mismatch"} <= reasons
    assert reasons <= set(TIER1_TARGET_REASON_CODES)
    assert gate["status"] == "fail"


def test_active_scope_can_be_enrichment_ready_but_not_science_ready(tmp_path) -> None:
    shard = tmp_path / "injections.h5"
    _write_injection_shard(shard)
    cadence_path = tmp_path / "cadence.csv"
    _, cadence_manifest = _write_cadence_evidence(cadence_path)
    cadence_manifest_path = tmp_path / "cadence.json"
    cadence_manifest_path.write_text(json.dumps(cadence_manifest) + "\n")
    _bind_injection_shard_to_cadence_reference(
        shard,
        cadence_path,
        cadence_manifest_path,
    )
    config = _small_injection_config(
        min_population_targets=40,
        min_targets_per_magnitude_bin=3,
        expected_injection_shard_sha256=(file_sha256(shard),),
        expected_cadence_reference_sha256=file_sha256(cadence_path),
        expected_cadence_reference_manifest_sha256=file_sha256(
            cadence_manifest_path
        ),
    )
    compact = tmp_path / "compact.h5"
    _write_compact(compact)
    targets, apertures, pairs = audit_compact_population(
        compact,
        apertures=config.apertures,
        cadence_reference=_cadence_reference(),
    )
    injections, injection_manifest = summarize_fixed_injection_shards(
        [shard],
        sector=56,
        cadence_reference_path=cadence_path,
        cadence_reference_manifest_path=cadence_manifest_path,
    )
    independent, independent_manifest = _passing_independent(
        cadence_reference_sha256=file_sha256(cadence_path),
        cadence_reference_manifest_sha256=file_sha256(cadence_manifest_path),
    )
    evaluated, bins, target_flags = evaluate_tier1_gates(
        sector=56,
        config=config,
        tier0_summary=_passing_tier0(TEST_COMPACT_SHA256),
        tier0_summary_sha256=TEST_TIER0_SUMMARY_SHA256,
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
    assert target_flags["tier1_target_searchable"].all()

    review_targets = targets.copy()
    review_targets["quality0_fraction"] = 0.75
    reviewed, _, reviewed_target_flags = evaluate_tier1_gates(
        sector=56,
        config=config,
        tier0_summary=_passing_tier0(TEST_COMPACT_SHA256),
        tier0_summary_sha256=TEST_TIER0_SUMMARY_SHA256,
        target_metrics=review_targets,
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
    assert reviewed["status"] == "review"
    assert reviewed["passed"] is True
    assert reviewed["enrichment_ready"] is True
    assert reviewed["blocking_gates"] == []
    assert reviewed["warning_gates"] == ["cadence_and_finite_data"]
    assert reviewed_target_flags["tier1_target_searchable"].all()

    old_tier0 = _passing_tier0(TEST_COMPACT_SHA256)
    old_tier0["contract_version"] = "a2v1_photometric_qa_v1"
    failed, _, _ = evaluate_tier1_gates(
        sector=56,
        config=config,
        tier0_summary=old_tier0,
        tier0_summary_sha256=TEST_TIER0_SUMMARY_SHA256,
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


def test_end_to_end_runner_writes_a_scoped_pass(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
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
    config = replace(
        config,
        expected_tier0_summary_sha256=file_sha256(tier0_path),
    )

    cadence_path = tmp_path / "cadence.csv"
    _, cadence_manifest = _write_cadence_evidence(cadence_path)
    cadence_manifest_path = tmp_path / "cadence.json"
    cadence_manifest_path.write_text(json.dumps(cadence_manifest))
    _bind_injection_shard_to_cadence_reference(
        shard,
        cadence_path,
        cadence_manifest_path,
    )
    config = replace(
        config,
        expected_injection_shard_sha256=(file_sha256(shard),),
        expected_cadence_reference_sha256=file_sha256(cadence_path),
        expected_cadence_reference_manifest_sha256=file_sha256(
            cadence_manifest_path
        ),
    )
    config_path.write_text(json.dumps(asdict(config)))

    independent, independent_manifest = _passing_independent(
        current_compact_sha256=compact_sha256,
        cadence_reference_sha256=file_sha256(cadence_path),
        cadence_reference_manifest_sha256=file_sha256(cadence_manifest_path),
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
    assert summary["target_qa"]["n_searchable"] == 40
    assert summary["target_qa"]["n_excluded"] == 0
    assert summary["target_qa"]["observation_key"] == ["sector", "tic"]
    assert gate_json.exists()
    assert (out_dir / "target_metrics.parquet").exists()
    eligibility_path = out_dir / "target_eligibility.csv"
    eligibility = pd.read_csv(eligibility_path, dtype={"gaia_dr3_source_id": "string"})
    assert len(eligibility) == 40
    assert not eligibility.duplicated(["sector", "tic"]).any()
    assert eligibility.loc[
        eligibility["tic"].eq(WD1856_TIC), "gaia_dr3_source_id"
    ].iloc[0] == "2098419251571450880"
    assert set(
        reason
        for value in eligibility["tier1_target_qa_reasons"].fillna("").astype(str)
        for reason in value.split(";")
        if reason
    ) <= set(TIER1_TARGET_REASON_CODES)
    assert (out_dir / "detector_summary.csv").exists()
    for name in (
        "tier1_qa_diagnostics.png",
        "tier1_qa_diagnostics.pdf",
        "tier1_detector_eligibility.png",
        "tier1_detector_eligibility.pdf",
    ):
        assert (out_dir / name).stat().st_size > 0
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

    failing_shard = tmp_path / "failing_injections.h5"
    _write_injection_shard(failing_shard)
    _bind_injection_shard_to_cadence_reference(
        failing_shard,
        cadence_path,
        cadence_manifest_path,
    )
    with h5py.File(failing_shard, "r+") as h5:
        del h5["injections/inj_00_000"].attrs["DET_FLUX_ADP_SML_scale"]
    failure_config = replace(
        config,
        expected_injection_shard_sha256=(file_sha256(failing_shard),),
    )
    failure_config_path = tmp_path / "failure_config.json"
    failure_config_path.write_text(json.dumps(asdict(failure_config)))
    failure_out = tmp_path / "failure_out"
    failure_gate_path = tmp_path / "failure_gate.json"
    with pytest.raises(ValueError, match="fixed-injection preflight failed"):
        run_a2v1_tier1_qa(
            sector=56,
            config_path=failure_config_path,
            tier0_summary_path=tier0_path,
            compact_lc=compact,
            cadence_reference_path=cadence_path,
            cadence_reference_manifest_path=cadence_manifest_path,
            injection_source_parity_path=injection_parity_path,
            injection_shards=[failing_shard],
            independent_metrics_path=independent_path,
            independent_manifest_path=independent_manifest_path,
            out_dir=failure_out,
            gate_json=failure_gate_path,
        )
    assert (failure_out / "fixed_injection_metrics.csv").exists()
    failure_manifest_path = failure_out / "fixed_injection_manifest.json"
    assert failure_manifest_path.exists()
    failure_gate = json.loads(failure_gate_path.read_text())
    assert failure_gate["status"] == "fail"
    assert failure_gate["passed"] is False
    assert failure_gate["enrichment_ready"] is False
    assert failure_gate["science_ready"] is False
    assert failure_gate["population_scan_started"] is False
    assert failure_gate["failure_stage"] == "fixed_injection_preservation"
    failure_manifest = json.loads(failure_manifest_path.read_text())
    assert failure_manifest["preflight_failure"]["population_scan_started"] is False
    assert failure_manifest["metrics_file_sha256"] == file_sha256(
        failure_out / "fixed_injection_metrics.csv"
    )
    assert not (failure_out / "target_metrics.parquet").exists()

    original_plot = tier1_qa.plot_detector_eligibility

    def mutate_evidence_after_plot(*args, **kwargs) -> None:
        original_plot(*args, **kwargs)
        tier0_path.write_text(tier0_path.read_text() + "\n")

    monkeypatch.setattr(
        tier1_qa,
        "plot_detector_eligibility",
        mutate_evidence_after_plot,
    )
    race_out = tmp_path / "race_out"
    race_gate = tmp_path / "race_gate.json"
    with pytest.raises(ValueError, match="evidence inputs changed"):
        run_a2v1_tier1_qa(
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
            out_dir=race_out,
            gate_json=race_gate,
        )
    assert not race_gate.exists()
    assert not (race_out / "summary.json").exists()


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
