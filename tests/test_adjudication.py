from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from twirl.io.hlsp import HLSPLightCurve
from twirl.vetting.adjudication import (
    ADJUDICATION_VET_SHEET_VERSION,
    REPEAT_QUOTAS,
    _public_and_manifest,
    _select_repeat_sources,
    join_browser_labels,
    place_blinded_repeats,
    reanchor_deprecated_sources,
    verify_queue_contract,
)
from twirl.vetting.adjudication_audit import (
    add_harmonic_cnn_targets,
    harmonic_target_for_factor,
    load_adjudication_reviews,
    normalize_period_factor,
    parse_note_period_factor,
    repeat_agreement,
)
from twirl.vetting.label_io import canonicalize_candidate_key, normalize_review_queue
from twirl.vetting.recovery50_teacher import leakage_columns
from twirl.vetting.two_aperture import _plot_harmonic_fold


def _unique_source_rows() -> pd.DataFrame:
    signal_labels = (
        ["planet_like"] * 40
        + ["eclipsing_binary_or_pceb"] * 30
        + ["stellar_variability"] * 130
        + ["wide_transit_like"] * 23
    )
    control_labels = ["instrumental_or_systematic"] * 50 + ["uncertain"] * 50
    labels = signal_labels + control_labels
    rows = []
    for index, label in enumerate(labels):
        signal = index < len(signal_labels)
        rows.append(
            {
                "source_uid": f"source:{index}",
                "source_cohort": "signal" if signal else "teacher2k_control",
                "source_review_id": f"real:{index}",
                "origin_queue": "toy",
                "labeling_era": "recent",
                "first_human_label": label,
                "raw_human_label": label,
                "raw_human_label_source": "human",
                "raw_human_labeler": "tester",
                "raw_human_notes": "",
                "raw_human_updated_utc": "2026-07-10T00:00:00Z",
                "historical_signal_label": label if signal else "",
                "ephemeris_source": "current_adp_sml",
                "tic": 1_000_000 + index,
                "sector": 56,
                "cam": 1,
                "ccd": 1,
                "tmag": 17.0,
                "period_d": 1.0 + index / 1000.0,
                "t0_bjd": 2459825.0,
                "duration_min": 10.0,
                "depth": 0.1,
                "depth_snr": 5.0,
                "sde_max": 8.0,
                "aperture_period_rel_delta": 0.0,
                "aperture_disagreement_flag": False,
            }
        )
    return pd.DataFrame(rows)


def test_adjudication_sheet_version_is_stable() -> None:
    assert ADJUDICATION_VET_SHEET_VERSION == "S56-ADP-HV2"


def test_fixed_queue_layout_has_343_rows_and_hidden_repeats() -> None:
    sources = _unique_source_rows()
    repeats = _select_repeat_sources(sources, repeat_quotas=REPEAT_QUOTAS, seed=56)
    placed = place_blinded_repeats(sources, repeats, seed=56)
    public, manifest = _public_and_manifest(placed)
    result = verify_queue_contract(public, manifest)

    assert result["passed"], result["failures"]
    assert len(public) == 343
    assert manifest["source_uid"].nunique() == 323
    assert int(manifest["is_repeat"].sum()) == 20
    assert "raw_human_label" not in public
    assert "repeat_group" not in public
    assert "source_cohort" not in public
    assert public["twirl_vet_sheet_pdf_name"].fillna("").eq("").all()


def test_normalized_join_recovers_labels_when_upstream_ids_repeat(tmp_path: Path) -> None:
    queue = pd.DataFrame(
        [
            {"row_id": 8, "review_id": "a", "tic": 1, "sector": 56, "period_d": 1.0, "t0_bjd": 2.0, "source_bucket": "real"},
            {"row_id": 8, "review_id": "b", "tic": 2, "sector": 56, "period_d": 2.0, "t0_bjd": 3.0, "source_bucket": "real"},
            {"row_id": 9, "review_id": "c", "tic": 3, "sector": 56, "period_d": 3.0, "t0_bjd": 4.0, "source_bucket": "real"},
        ]
    )
    labels = pd.DataFrame(
        [
            {"row_id": 0, "candidate_key": "1|56|1.0|2.0|real", "label": "planet_like"},
            {"row_id": 1, "candidate_key": "2|56|2.0|3.0|real", "label": "uncertain"},
            {"row_id": 2, "candidate_key": "3|56|3.0|4.0|real", "label": "stellar_variability"},
        ]
    )
    queue_path = tmp_path / "queue.csv"
    labels_path = tmp_path / "labels.csv"
    queue.to_csv(queue_path, index=False)
    labels.to_csv(labels_path, index=False)

    joined = join_browser_labels(queue_path, labels_path)

    assert joined["row_id"].tolist() == [0, 1, 2]
    assert joined["browser_label"].tolist() == ["planet_like", "uncertain", "stellar_variability"]


def test_candidate_key_canonicalization_absorbs_parser_roundoff_only() -> None:
    legacy = "162462912||0.3326643093383712|2459825.535826442|real_eb_miner_priority"
    reparsed = "162462912||0.3326643093383712|2459825.5358264414|real_eb_miner_priority"

    assert canonicalize_candidate_key(legacy) == canonicalize_candidate_key(reparsed)
    assert canonicalize_candidate_key(legacy) != canonicalize_candidate_key(
        reparsed.replace("0.3326643093383712", "0.4")
    )


def test_deprecated_ephemeris_is_replaced_by_adp_small_top1() -> None:
    deprecated = pd.DataFrame(
        [{"tic": 10, "period_d": 9.0, "t0_bjd": 1.0, "duration_min": 30.0, "ephemeris_source": "stale"}]
    )
    peaks = pd.DataFrame(
        [
            {"tic": 10, "aperture": "DET_FLUX_ADP", "peak_rank": 1, "status": "ok", "period_d": 4.0, "t0_bjd": 5.0, "duration_min": 12.0, "sde": 8.0},
            {"tic": 10, "aperture": "DET_FLUX_ADP_SML", "peak_rank": 2, "status": "ok", "period_d": 3.0, "t0_bjd": 6.0, "duration_min": 13.0, "sde": 9.0},
            {"tic": 10, "aperture": "DET_FLUX_ADP_SML", "peak_rank": 1, "status": "ok", "period_d": 2.0, "t0_bjd": 7.0, "duration_min": 14.0, "sde": 10.0},
        ]
    )

    rebuilt = reanchor_deprecated_sources(deprecated, peaks)

    assert rebuilt.loc[0, "period_d"] == 2.0
    assert rebuilt.loc[0, "t0_bjd"] == 7.0
    assert rebuilt.loc[0, "duration_min"] == 14.0
    assert rebuilt.loc[0, "sde_max"] == 10.0
    assert rebuilt.loc[0, "ephemeris_source"] == "current_adp_sml_bls_top1"


def test_period_factor_and_repeat_agreement_are_separate_from_label() -> None:
    assert normalize_period_factor("0.5") == ("0.5", "refolded", 0.5)
    unresolved = normalize_period_factor("unresolved")
    assert unresolved[:2] == ("unresolved", "unresolved")
    rows = pd.DataFrame(
        [
            {
                "repeat_group": "repeat:00",
                "repeat_occurrence": "original",
                "source_uid": "source:1",
                "row_id": 1,
                "repeat_reference_label": "planet_like",
                "adjudicated_row_label": "planet_like",
                "adjudicated_row_period_factor": "1",
            },
            {
                "repeat_group": "repeat:00",
                "repeat_occurrence": "blind_repeat",
                "source_uid": "source:1",
                "row_id": 201,
                "repeat_reference_label": "planet_like",
                "adjudicated_row_label": "planet_like",
                "adjudicated_row_period_factor": "0.5",
            },
        ]
    )

    pairs, summary = repeat_agreement(rows)

    assert bool(pairs.loc[0, "label_agreement"])
    assert not bool(pairs.loc[0, "exact_agreement"])
    assert bool(pairs.loc[0, "discordant"])
    assert summary["n_discordant_pairs"] == 1


def test_harmonic_cnn_targets_use_note_override_and_broad_preserve_only() -> None:
    frame = pd.DataFrame(
        {
            "human_label": [
                "eclipsing_binary_or_pceb",
                "stellar_variability",
                "wide_transit_like",
                "planet_like",
                "instrumental_or_systematic",
                "uncertain",
                "skip",
            ],
            "human_notes": ["this should be 3P", "maybe P/14", "", "", "", "", ""],
            "period_d": [2.0] * 7,
            "adjudication_final": [True] * 7,
            "adjudicated_period_factor": ["1", "1", "0.25", "0.5", "1", "1", "1"],
            "adjudicated_period_status": [
                "review_period",
                "unresolved",
                "refolded",
                "refolded",
                "review_period",
                "review_period",
                "review_period",
            ],
            "refold_factor": [np.nan] * 7,
        }
    )

    out = add_harmonic_cnn_targets(frame)

    assert out.loc[0, "note_period_factor"] == 3.0
    assert out.loc[0, "effective_period_factor"] == 3.0
    assert out.loc[0, "effective_period_d"] == 6.0
    assert out.loc[0, "period_factor_source"] == "note"
    assert out.loc[0, "harmonic_target_v1"] == "3p"
    assert out.loc[0, "morphology_target_v1"] == "eclipse_contact"
    assert np.isclose(out.loc[1, "note_period_factor"], 1.0 / 14.0)
    assert out.loc[1, "morphology_target_v1"] == "smooth_variable"
    assert bool(out.loc[1, "morphology_include_v1"])
    assert out.loc[1, "period_task"] == "variable_period_refinement"
    assert not bool(out.loc[1, "harmonic_include_v1"])
    assert out.loc[2, "morphology_target_v1"] == ""
    assert out.loc[2, "preserve_target_v1"] == "preserve"
    assert bool(out.loc[2, "broad_preserve_only"])
    assert out.loc[2, "harmonic_target_v1"] == "p_over_4"
    assert out.loc[3, "morphology_target_v1"] == "planet_like"
    assert out.loc[4, "morphology_target_v1"] == "other"
    assert out.loc[5, "preserve_target_v1"] == "reject"
    assert not bool(out.loc[6, "morphology_include_v1"])
    assert not bool(out.loc[6, "preserve_include_v1"])


def test_note_period_factor_parser_is_strict_and_supports_thirds() -> None:
    assert parse_note_period_factor("period of 3P") == 3.0
    assert parse_note_period_factor("refold at P/3") == 1.0 / 3.0
    assert parse_note_period_factor("maybe p/17??") == 1.0 / 17.0
    assert parse_note_period_factor("this is at half period") == 0.5
    assert parse_note_period_factor("try double the period") == 2.0
    assert np.isnan(parse_note_period_factor("possible harmonic"))
    assert harmonic_target_for_factor(1.0 / 3.0) == "p_over_3"
    assert harmonic_target_for_factor(3.0) == "3p"
    assert harmonic_target_for_factor(1.0 / 14.0) == ""


def test_injected_harmonic_target_uses_truth_ratio_without_changing_morphology() -> None:
    frame = pd.DataFrame(
        {
            "human_label": ["planet_like", "planet_like"],
            "human_notes": ["", ""],
            "source_kind": ["injection_recovery", "injection_recovery"],
            "is_injected_row": [True, True],
            "period_d": [1.0, 1.0],
            "truth_period_d": [2.001, 1.37],
            "adjudication_final": [False, False],
            "adjudicated_period_factor": [np.nan, np.nan],
            "adjudicated_period_status": ["", ""],
            "refold_factor": [np.nan, np.nan],
        }
    )

    out = add_harmonic_cnn_targets(frame)

    assert out.loc[0, "morphology_target_v1"] == "planet_like"
    assert out.loc[0, "effective_period_factor"] == 2.0
    assert out.loc[0, "harmonic_target_v1"] == "2p"
    assert bool(out.loc[0, "harmonic_include_v1"])
    assert out.loc[0, "period_factor_source"] == "injection_truth_harmonic"
    assert out.loc[1, "morphology_target_v1"] == "planet_like"
    assert not bool(out.loc[1, "harmonic_include_v1"])
    assert out.loc[1, "period_factor_source"] == "injection_truth_unresolved"


def test_review_join_computes_corrected_period_and_unresolved(tmp_path: Path) -> None:
    queue = normalize_review_queue(
        pd.DataFrame(
            [
                {"review_id": "a", "tic": 1, "sector": 56, "period_d": 4.0, "t0_bjd": 5.0, "source_bucket": "blind"},
                {"review_id": "b", "tic": 2, "sector": 56, "period_d": 8.0, "t0_bjd": 6.0, "source_bucket": "blind"},
            ]
        )
    )
    manifest = queue.loc[:, ["row_id", "candidate_key"]].copy()
    manifest["source_uid"] = ["source:1", "source:2"]
    manifest["review_period_d"] = queue["period_d"]
    manifest["repeat_group"] = ""
    labels = pd.DataFrame(
        [
            {"row_id": 0, "candidate_key": queue.loc[0, "candidate_key"], "label": "planet_like", "period_factor": "0.5"},
            {"row_id": 1, "candidate_key": queue.loc[1, "candidate_key"], "label": "planet_like", "period_factor": "unresolved"},
        ]
    )
    queue_path = tmp_path / "queue.csv"
    manifest_path = tmp_path / "manifest.csv"
    labels_path = tmp_path / "labels.csv"
    queue.to_csv(queue_path, index=False)
    manifest.to_csv(manifest_path, index=False)
    labels.to_csv(labels_path, index=False)

    rows = load_adjudication_reviews(
        queue_csv=queue_path,
        labels_csv=labels_path,
        manifest_csv=manifest_path,
    )

    assert rows.loc[0, "adjudicated_row_period_d"] == 2.0
    assert rows.loc[0, "adjudicated_row_period_status"] == "refolded"
    assert np.isnan(rows.loc[1, "adjudicated_row_period_d"])
    assert rows.loc[1, "adjudicated_row_period_status"] == "unresolved"


def test_adjudication_columns_are_leakage() -> None:
    columns = [
        "human_label_adjudicated",
        "adjudicated_period_factor",
        "repeat_group",
        "raw_human_label",
        "period_factor",
        "morphology_target_v1",
        "harmonic_target_v1",
        "note_period_factor",
        "effective_period_factor",
        "period_task",
        "anchor_sde",
    ]
    leaks = leakage_columns(columns)

    assert set(columns[:-1]).issubset(set(leaks))
    assert "anchor_sde" not in leaks


def test_harmonic_fold_shows_raw_bins_and_both_phase_markers() -> None:
    time = np.linspace(2825.0, 2852.0, 500)
    phase = ((time - 2825.2 + 0.5) % 1.0) - 0.5
    flux = 1.0 - 0.08 * (np.abs(phase) < 0.03)
    aperture = "DET_FLUX_ADP_SML"
    lc = HLSPLightCurve(
        tic=1,
        tmag=17.0,
        sector=56,
        cam=1,
        ccd=1,
        ra=0.0,
        dec=0.0,
        time=time,
        cadenceno=np.arange(len(time)),
        orbitid=np.ones(len(time)),
        quality=np.zeros(len(time), dtype=int),
        flux={aperture: flux},
        path=Path("toy.fits"),
    )
    fig, ax = plt.subplots()

    _plot_harmonic_fold(
        ax,
        lc,
        aperture,
        period_d=1.0,
        t0_bjd=2459825.2,
        duration_min=60.0,
        depth=0.08,
        title="P",
        show_ylabel=True,
    )

    vertical = [float(line.get_xdata()[0]) for line in ax.lines if len(np.unique(line.get_xdata())) == 1]
    assert 0.0 in vertical
    assert 0.5 in vertical
    assert len(ax.collections) >= 2
    assert ax.get_xlim() == (-0.25, 0.75)
    plt.close(fig)
