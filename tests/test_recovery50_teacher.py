from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.stage5_validation.build_s56_recovery50_teacher_queue import _parse_real_quotas
from twirl.vetting.recovery50_teacher import (
    add_deterministic_splits,
    add_label_roles,
    bin_folded_channel,
    leakage_columns,
    metadata_feature_columns,
    shape_feature_columns,
    train_teacher_student_smoke,
)


def test_parse_recovery50_real_quota_overrides() -> None:
    quotas = _parse_real_quotas(["real_eb_pceb_like=900", "real_high_sde_planet_like=700"])

    assert quotas == [("real_eb_pceb_like", 900), ("real_high_sde_planet_like", 700)]


def test_human_label_remains_main_target_for_injected_rows() -> None:
    df = pd.DataFrame(
        [
            {
                "source_kind": "injection_recovery",
                "human_label": "instrumental_or_systematic",
                "topn_recovery_status": "bls_top1_recovered",
            },
            {
                "source_kind": "injection_recovery",
                "human_label": "planet_like",
                "topn_recovery_status": "bls_top1_recovered",
            },
        ]
    )

    out = add_label_roles(df)

    assert out.loc[0, "truth_signal_present"]
    assert out.loc[0, "bls_truth_match"]
    assert out.loc[0, "main_teacher_target"] == "instrumental_or_systematic"
    assert out.loc[1, "main_teacher_target"] == "planet_like"


def test_uncertain_maps_to_negative_teacher_target_without_overwriting_raw_label() -> None:
    df = pd.DataFrame(
        [
            {
                "source_kind": "real_candidate",
                "human_label": "uncertain",
                "topn_recovery_status": "",
            },
        ]
    )

    out = add_label_roles(df)

    assert out.loc[0, "human_label"] == "uncertain"
    assert out.loc[0, "main_teacher_target"] == "instrumental_or_systematic"
    assert bool(out.loc[0, "main_teacher_include"])


def test_leakage_columns_are_detected_and_allowlist_is_clean() -> None:
    leaks = leakage_columns(
        [
            "truth_period_d",
            "recovery_status",
            "source_kind",
            "pseudo_label",
            "sde_max",
            "shape_flux_adp_sml_bin_000",
        ]
    )

    assert leaks == ["pseudo_label", "recovery_status", "source_kind", "truth_period_d"]
    df = pd.DataFrame({"sde_max": [1.0], "truth_period_d": [2.0], "source_kind": ["real_candidate"]})
    assert metadata_feature_columns(df) == ["sde_max"]
    assert shape_feature_columns(pd.DataFrame({"shape_flux_adp_sml_bin_000": [0.0]})) == [
        "shape_flux_adp_sml_bin_000"
    ]


def test_folded_channel_recovers_centered_dip() -> None:
    period = 1.0
    t0 = 2457000.5
    duration_min = 60.0
    time = np.linspace(t0 - 2.0, t0 + 2.0, 500)
    phase_hr = (((time - t0 + 0.5 * period) % period) - 0.5 * period) * 24.0
    flux = np.ones_like(time)
    flux[np.abs(phase_hr) < 0.5] -= 0.2
    quality = np.zeros_like(time, dtype=int)

    values, mask, counts, stats = bin_folded_channel(
        time=time,
        flux=flux,
        quality=quality,
        period_d=period,
        t0_bjd=t0,
        duration_min=duration_min,
        n_bins=21,
        window_durations=3.0,
    )

    assert stats["status"] == "ok"
    assert mask.sum() > 10
    assert counts.max() > 0
    assert values[len(values) // 2] < -0.1
    assert stats["folded_depth"] > 0.1


def test_teacher_smoke_trains_without_truth_features(tmp_path) -> None:
    rows = []
    for i in range(80):
        positive = i < 40
        rows.append(
            {
                "row_id": i,
                "review_id": f"row-{i}",
                "tic": 1000 + i,
                "sector": 56,
                "source_kind": "injection_recovery" if i % 3 == 0 else "real_candidate",
                "human_label": "planet_like" if positive else "instrumental_or_systematic",
                "main_teacher_target": "planet_like" if positive else "instrumental_or_systematic",
                "main_teacher_include": True,
                "truth_period_d": 1.0 if positive else 5.0,
                "sde_max": 30.0 if positive else 5.0,
                "depth": 0.5 if positive else 0.05,
                "shape_flux_adp_sml_bin_000": -0.2 if positive else 0.02,
                "shape_mask_adp_sml_bin_000": 1.0,
            }
        )
    table = add_deterministic_splits(pd.DataFrame(rows), validation_fraction=0.2, test_fraction=0.2)

    summary = train_teacher_student_smoke(feature_table=table, out_dir=tmp_path, min_class_count=40)

    assert summary["profiles"]["metadata_only"]["status"] == "ok"
    assert summary["profiles"]["shape_only"]["status"] == "ok"
    for profile in ("metadata_only", "shape_only", "combined"):
        features = summary["profiles"][profile]["feature_columns"]
        assert "truth_period_d" not in features
        assert not leakage_columns(features)
