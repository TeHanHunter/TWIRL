from __future__ import annotations

import pandas as pd
import pytest

from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURES,
    ADP_ONLY_APERTURE_SIGNATURE,
    ADP_ONLY_CONTRACT_VERSION,
    ADP_ONLY_INJECTION_H5_BY_ORIGIN,
    assert_adp_only_training_frame,
    assert_adp_only_tensor_rows,
    build_adp_only_training_frame,
    canonical_det_flux_columns,
    classify_period_relation,
    validate_adp_only_apertures,
)


def test_period_relation_classifies_exact_harmonic_and_unrelated() -> None:
    left = pd.Series([1.0, 2.0, 0.73, float("nan")])
    right = pd.Series([1.001, 1.0, 1.0, 1.0])
    assert classify_period_relation(left, right).tolist() == [
        "exact",
        "harmonic",
        "unrelated",
        "missing",
    ]


def test_build_adp_only_training_frame_drops_canonical_results(tmp_path) -> None:
    base = pd.DataFrame(
        [
            {
                "review_id": "real:1",
                "source_kind": "real_candidate",
                "period_d": 7.0,
                "t0_bjd": 2459000.0,
                "duration_min": 30.0,
                "display_period_d": 1.5,
                "display_t0_bjd": 2459001.0,
                "display_duration_min": 12.0,
                "display_ephemeris_source": "twirl_vet_anchor",
                "main_teacher_include": True,
                "compact_morphology_include": True,
                "training_split": "train",
                "sde_DET_FLUX_SML": 99.0,
                "rank_DET_FLUX_SML": 1,
                "centroid_z": 8.0,
                "queue_period_d": 7.0,
            },
            {
                "review_id": "inj:2",
                "source_kind": "injection_recovery",
                "origin_queue": "s56_recovery50_teacher_queue_next4k",
                "source_h5": "old_canonical_pair.h5",
                "period_d": 4.0,
                "display_period_d": 2.0,
                "display_t0_bjd": 2459002.0,
                "display_duration_min": 6.0,
                "display_ephemeris_source": "twirl_vet_anchor",
                "main_teacher_include": True,
                "compact_morphology_include": True,
                "training_split": "train",
            },
        ]
    )
    metrics = tmp_path / "metrics.csv"
    pd.DataFrame(
        [
            {
                "review_id": "real:1",
                "anchor_aperture": ADP_ONLY_APERTURES[0],
                "anchor_period_d": 1.5,
                "anchor_t0_bjd": 2459001.0,
                "anchor_duration_min": 12.0,
                "anchor_sde": 10.0,
                "adp_sml_peak_rank": 1,
                "adp_sml_period_d": 1.5,
                "adp_sml_depth": 0.2,
                "adp_sml_depth_snr": 5.0,
                "adp_sml_sde": 10.0,
                "adp_period_d": 1.51,
                "adp_sde": 9.0,
                "aperture_period_rel_delta": 0.006,
            },
            {
                "review_id": "inj:2",
                "anchor_aperture": ADP_ONLY_APERTURES[0],
                "anchor_period_d": 2.0,
                "anchor_t0_bjd": 2459002.0,
                "anchor_duration_min": 6.0,
                "anchor_sde": 8.0,
                "adp_sml_peak_rank": 1,
                "adp_sml_period_d": 2.0,
                "adp_sml_sde": 8.0,
                "sml_period_d": 2.0,
            },
        ]
    ).to_csv(metrics, index=False)

    out, summary = build_adp_only_training_frame(base, metrics_tables=(metrics,))

    assert out.loc[0, "period_d"] == 1.5
    assert out.loc[0, "sde_max"] == 10.0
    assert bool(out.loc[0, "adp_only_review_ok"])
    assert not bool(out.loc[1, "adp_only_review_ok"])
    assert out.loc[1, "adp_only_exclusion_reason"] == "missing_adp_primary_supplement"
    assert not bool(out.loc[1, "main_teacher_include"])
    assert out.loc[1, "source_h5"] == ADP_ONLY_INJECTION_H5_BY_ORIGIN[
        "s56_recovery50_teacher_queue_next4k"
    ]
    assert summary["n_adp_only_review_ok"] == 1
    assert not canonical_det_flux_columns(out.columns)
    assert "centroid_z" not in out
    assert "queue_period_d" not in out
    assert_adp_only_training_frame(out)


def test_adp_only_aperture_validation_is_strict() -> None:
    assert validate_adp_only_apertures(ADP_ONLY_APERTURES) == ADP_ONLY_APERTURES
    with pytest.raises(ValueError, match="active S56 ML requires"):
        validate_adp_only_apertures(("DET_FLUX_ADP_SML", "DET_FLUX_SML"))


def test_adp_contract_rejects_historical_training_table() -> None:
    frame = pd.DataFrame(
        {
            "adp_only_contract_version": [ADP_ONLY_CONTRACT_VERSION],
            "sde_DET_FLUX_SML": [8.0],
        }
    )
    with pytest.raises(ValueError, match="canonical DET_FLUX"):
        assert_adp_only_training_frame(frame)


def test_tensor_contract_requires_ordered_adp_pair() -> None:
    valid = pd.DataFrame(
        {
            "adp_only_contract_version": [ADP_ONLY_CONTRACT_VERSION],
            "tensor_apertures": [ADP_ONLY_APERTURE_SIGNATURE],
        }
    )
    assert_adp_only_tensor_rows(valid)
    with pytest.raises(ValueError, match="tensor aperture signatures"):
        assert_adp_only_tensor_rows(valid.assign(tensor_apertures="DET_FLUX_ADP"))
