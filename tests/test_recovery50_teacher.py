from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.stage5_validation.build_s56_recovery50_teacher_queue import _parse_real_quotas
from twirl.vetting.recovery50_teacher import (
    add_display_ephemeris,
    add_deterministic_splits,
    add_harmonic_ephemeris_annotations,
    add_label_roles,
    bin_folded_channel,
    infer_harmonic_refold,
    leakage_columns,
    metadata_feature_columns,
    select_model_ephemeris,
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


def test_bls_truth_match_accepts_aperture_specific_recovery_columns() -> None:
    out = add_label_roles(
        pd.DataFrame(
            {
                "source_kind": ["injection_recovery", "injection_recovery"],
                "human_label": ["planet_like", "uncertain"],
                "topn_exact_recovered_DET_FLUX_ADP_SML": [True, False],
                "topn_harmonic_match_DET_FLUX_ADP_SML": [False, True],
            }
        )
    )

    assert out["bls_truth_match"].tolist() == [True, True]


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


def test_display_ephemeris_prefers_vet_sheet_anchor(tmp_path) -> None:
    rows = pd.DataFrame(
        [
            {
                "review_id": "real:1",
                "period_d": 1.0,
                "t0_bjd": 2459000.0,
                "duration_min": 30.0,
            },
            {
                "review_id": "real:2",
                "period_d": 2.0,
                "t0_bjd": 2459001.0,
                "duration_min": 12.0,
            },
        ]
    )
    metrics_path = tmp_path / "metrics.csv"
    pd.DataFrame(
        [
            {
                "review_id": "real:1",
                "anchor_period_d": 1.5,
                "anchor_t0_bjd": 2459000.25,
                "anchor_duration_min": 21.0,
                "anchor_aperture": "DET_FLUX_ADP_SML",
                "anchor_sde": 9.0,
                "twirl_vet_sheet_name": "real_1.png",
                "twirl_vet_status": "ok",
            }
        ]
    ).to_csv(metrics_path, index=False)

    out = add_display_ephemeris(rows, metrics_tables=(metrics_path,))

    assert out.loc[0, "display_period_d"] == 1.5
    assert out.loc[0, "display_t0_bjd"] == 2459000.25
    assert out.loc[0, "display_duration_min"] == 21.0
    assert out.loc[0, "display_ephemeris_source"] == "twirl_vet_anchor"
    assert out.loc[1, "display_period_d"] == 2.0
    assert out.loc[1, "display_ephemeris_source"] == "queue"


def test_model_ephemeris_selection_order() -> None:
    row = pd.Series(
        {
            "period_d": 1.0,
            "t0_bjd": 2459000.0,
            "duration_min": 30.0,
            "display_period_d": 2.0,
            "display_t0_bjd": 2459001.0,
            "display_duration_min": 12.0,
            "display_ephemeris_source": "twirl_vet_anchor",
        }
    )

    period_d, t0_bjd, duration_min, source = select_model_ephemeris(row)

    assert period_d == 2.0
    assert t0_bjd == 2459001.0
    assert duration_min == 12.0
    assert source == "twirl_vet_anchor"


def test_harmonic_notes_set_refold_model_ephemeris_without_training_leakage() -> None:
    df = pd.DataFrame(
        [
            {
                "human_label": "skip",
                "human_notes": "this is transit like, but at half period",
                "display_period_d": 2.0,
                "display_t0_bjd": 2459001.0,
                "display_duration_min": 18.0,
                "display_ephemeris_source": "twirl_vet_anchor",
            }
        ]
    )

    out = add_harmonic_ephemeris_annotations(df)
    period_d, t0_bjd, duration_min, source = select_model_ephemeris(out.loc[0])

    assert infer_harmonic_refold("possible harmonic; refold at P/2") == (True, 0.5, "half_period_note")
    assert bool(out.loc[0, "harmonic_suspect"])
    assert out.loc[0, "refold_factor"] == 0.5
    assert out.loc[0, "refold_period_d"] == 1.0
    assert period_d == 1.0
    assert t0_bjd == 2459001.0
    assert duration_min == 18.0
    assert source.startswith("human_half_period_note")
    assert set(
        leakage_columns(
            [
                "harmonic_suspect",
                "refold_period_d",
                "ephemeris_status",
                "compact_morphology_target",
                "preserve_signal_target",
                "adp_sml_sde",
            ]
        )
    ) == {
        "compact_morphology_target",
        "ephemeris_status",
        "harmonic_suspect",
        "preserve_signal_target",
        "refold_period_d",
    }


def test_stellar_variability_is_preserved_as_compact_morphology() -> None:
    out = add_harmonic_ephemeris_annotations(
        pd.DataFrame(
            [
                {
                    "human_label": "stellar_variability",
                    "display_period_d": 1.0,
                    "display_t0_bjd": 2459000.0,
                    "display_duration_min": 20.0,
                }
            ]
        )
    )

    assert out.loc[0, "preserve_signal_target"] == "preserve_signal"
    assert out.loc[0, "compact_morphology_target"] == "stellar_variability"
    assert bool(out.loc[0, "compact_morphology_include"])


def test_leakage_columns_are_detected_and_allowlist_is_clean() -> None:
    leaks = leakage_columns(
        [
            "truth_period_d",
            "recovery_status",
            "source_kind",
            "pseudo_label",
            "adp_sml_sde",
            "shape_flux_adp_sml_bin_000",
        ]
    )

    assert leaks == ["pseudo_label", "recovery_status", "source_kind", "truth_period_d"]
    df = pd.DataFrame(
        {"adp_sml_sde": [1.0], "truth_period_d": [2.0], "source_kind": ["real_candidate"]}
    )
    assert metadata_feature_columns(df) == ["adp_sml_sde"]
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
                "adp_sml_sde": 30.0 if positive else 5.0,
                "adp_sml_depth": 0.5 if positive else 0.05,
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
