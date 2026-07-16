from __future__ import annotations

import numpy as np
import pandas as pd

from twirl.vetting.adp_only import ADP_ONLY_CONTRACT_VERSION
from twirl.vetting.harmonic_export import _native_input_mask
from twirl.vetting.harmonic_inference import prepare_inference_rows, rank_planet_enrichment
from twirl.vetting.teacher_active_learning import (
    A2V1_TEACHER_INPUT_CONTRACT,
    audit_legacy_franklin_queue,
    build_a2v1_label_transfer_table,
    build_active_learning_batch,
    effective_human_label,
    evaluate_a2v1_transfer_gate,
    evaluate_a2v1_scored_transfer,
    sector_rollout_readiness,
    teacher_v2_readiness,
)
from twirl.vetting.teacher_candidates import normalize_a2v1_peak_candidates


def test_legacy_franklin_rows_are_quarantined_from_adp_morphology() -> None:
    queue = pd.DataFrame(
        {
            "review_id": ["legacy-a", "legacy-b"],
            "tic": [10, 20],
            "sector": [56, 56],
            "period_d": [1.0, 2.0],
            "t0_bjd": [2459000.0, 2459001.0],
            "duration_min": [10.0, 20.0],
            "rep_aperture": ["DET_FLUX", "DET_FLUX_LAG"],
        }
    )
    peaks = pd.DataFrame(
        {
            "tic": [10, 10, 20, 20],
            "aperture": ["DET_FLUX_ADP_SML"] * 4,
            "peak_rank": [1, 2, 1, 2],
            "period_d": [1.0, 3.0, 0.7, 1.0],
        }
    )
    rows, summary = audit_legacy_franklin_queue(queue, peaks)
    assert not rows["active_adp_morphology_eligible"].any()
    assert rows["adp_re_review_required"].all()
    assert summary["adp_rank1_direct_match_fraction"] == 0.5
    assert summary["adp_topn_harmonic_match_fraction"] == 1.0


def test_effective_label_only_overrides_rows_with_adjudication() -> None:
    frame = pd.DataFrame(
        {
            "human_label": ["planet_like", "stellar_variability"],
            "human_label_adjudicated": ["eclipsing_binary_or_pceb", np.nan],
        }
    )
    assert effective_human_label(frame).tolist() == [
        "eclipsing_binary_or_pceb",
        "stellar_variability",
    ]


def _a2v1_peaks() -> pd.DataFrame:
    rows: list[dict] = []
    for tic in (10, 20):
        for rank, period in ((1, 1.0 + tic / 1000), (2, 2.0 + tic / 1000)):
            rows.append(
                {
                    "tic": tic,
                    "sector": 56,
                    "cam": 1,
                    "ccd": 2,
                    "tmag": 17.0,
                    "aperture": "DET_FLUX_ADP_SML",
                    "peak_rank": rank,
                    "period_d": period,
                    "t0_bjd": 2459000.0,
                    "duration_min": 10.0,
                    "depth": 0.1,
                    "depth_snr": 5.0,
                    "sde": 12.0,
                    "log_power": 1.0,
                    "status": "ok",
                    "bls_search_branch": "current_adp",
                    "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
                    "source_product_tag": "A2v1",
                }
            )
        rows.append(
            {
                "tic": tic,
                "sector": 56,
                "cam": 1,
                "ccd": 2,
                "tmag": 17.0,
                "aperture": "DET_FLUX_ADP",
                "peak_rank": 1,
                "period_d": 1.0 + tic / 1000,
                "t0_bjd": 2459000.0,
                "duration_min": 10.0,
                "depth": 0.11,
                "depth_snr": 5.5,
                "sde": 13.0,
                "log_power": 1.1,
                "status": "ok",
                "bls_search_branch": "current_adp",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
                "source_product_tag": "A2v1",
            }
        )
    return pd.DataFrame(rows)


def test_a2v1_candidate_table_uses_small_peaks_and_primary_context() -> None:
    candidates = normalize_a2v1_peak_candidates(_a2v1_peaks(), small_peaks_per_tic=2)
    assert len(candidates) == 4
    assert candidates["rep_aperture"].eq("DET_FLUX_ADP_SML").all()
    assert candidates["adp_period_d"].notna().all()
    assert candidates["native_input_include"].all()
    assert candidates["input_contract_version"].eq(A2V1_TEACHER_INPUT_CONTRACT).all()
    assert candidates["review_id"].is_unique


def test_a2v1_candidate_table_preserves_s57_sector_provenance() -> None:
    peaks = _a2v1_peaks().assign(sector=57)
    candidates = normalize_a2v1_peak_candidates(peaks, small_peaks_per_tic=1)
    assert candidates["sector"].eq(57).all()
    assert candidates["review_id"].str.startswith("s0057-A2v1-").all()


def test_planet_enrichment_ranking_keeps_one_best_candidate_per_tic() -> None:
    scored = pd.DataFrame(
        {
            "review_id": ["a", "b", "c"],
            "tic": [10, 10, 20],
            "sector": [57, 57, 57],
            "p_planet_like": [0.80, 0.90, 0.85],
            "p_preserve": [0.95, 0.70, 0.90],
            "std_p_planet_like": [0.02, 0.05, 0.03],
            "sde_max": [12.0, 20.0, 15.0],
        }
    )
    ranked = rank_planet_enrichment(scored)
    assert ranked["review_id"].tolist() == ["b", "c"]
    assert ranked["planet_rank"].tolist() == [1, 2]
    assert ranked["tic"].is_unique


def test_a2v1_transfer_uses_human_labels_and_matching_new_peak() -> None:
    labels = pd.DataFrame(
        {
            "tic": [10, 20],
            "period_d": [1.01, 1.02],
            "human_label": ["planet_like", "stellar_variability"],
            "source_kind": ["real_candidate", "real_candidate"],
        }
    )
    transfer, summary = build_a2v1_label_transfer_table(labels, _a2v1_peaks())
    assert summary["top10_harmonic_match_fraction"] == 1.0
    assert transfer["a2v1_scoring_eligible"].all()
    scores = pd.DataFrame(
        {
            "tic": [10, 20],
            "rep_peak_rank": [1, 1],
            "p_planet_like": [0.8, 0.1],
            "p_eclipse_contact": [0.05, 0.1],
            "p_smooth_variable": [0.05, 0.7],
            "p_other": [0.1, 0.1],
        }
    )
    evaluated, scored = evaluate_a2v1_scored_transfer(transfer, scores)
    assert len(evaluated) == 2
    assert scored["truth_policy"].startswith("human morphology")
    assert scored["human_morphology_metrics"]["accuracy"] == 1.0


def _score_rows(n: int = 1300) -> pd.DataFrame:
    rng = np.random.default_rng(56)
    probability = rng.dirichlet(np.ones(4), size=n)
    return pd.DataFrame(
        {
            "review_id": [f"source-{index}" for index in range(n)],
            "tic": np.arange(100_000, 100_000 + n),
            "sector": 56,
            "period_d": rng.uniform(0.12, 12.0, n),
            "t0_bjd": rng.uniform(2459000.0, 2459020.0, n),
            "duration_min": rng.uniform(3.0, 90.0, n),
            "sde_max": rng.uniform(5.0, 100.0, n),
            "p_planet_like": probability[:, 0],
            "p_eclipse_contact": probability[:, 1],
            "p_smooth_variable": probability[:, 2],
            "p_other": probability[:, 3],
            "p_preserve": rng.uniform(0.0, 1.0, n),
            "p_harmonic_P_over_2": rng.uniform(0.0, 1.0, n),
            "morphology_entropy": rng.uniform(0.0, np.log(4), n),
            "morphology_margin": rng.uniform(0.0, 1.0, n),
            "ensemble_disagreement": rng.uniform(0.0, 0.2, n),
            "source_kind": "real_candidate",
        }
    )


def test_active_learning_batch_has_locked_quotas_and_hidden_scores() -> None:
    scores = _score_rows()
    excluded = scores.loc[:49, ["tic"]].copy()
    queue, overlap, hidden, summary = build_active_learning_batch(
        scores, batch_index=1, excluded_tables=(excluded,)
    )
    assert len(queue) == 1000
    assert queue["tic"].nunique() == 1000
    assert len(overlap) == 100
    assert set(overlap["candidate_key"]).issubset(set(queue["candidate_key"]))
    assert "p_planet_like" not in queue
    assert "selection_bucket" not in queue
    assert "source_kind" not in queue
    assert "input_contract_version" not in queue
    assert "selection_bucket" in hidden
    assert hidden["selection_weight"].gt(0).all()
    assert set(queue["tic"]).isdisjoint(set(excluded["tic"]))
    assert summary["excluded_tics"] == 50
    assert summary["selection_bucket_counts"] == {
        "broad_dip": 100,
        "disagreement_harmonic": 150,
        "eclipse_contact": 200,
        "planet_preserve": 300,
        "smooth_variable": 150,
        "stratified_control": 100,
    }


def test_transfer_and_student_gates_require_real_support() -> None:
    compatibility = pd.DataFrame({"a2v1_topn_harmonic_match": [True] * 91 + [False] * 9})
    transfer = evaluate_a2v1_transfer_gate(
        compatibility,
        real_macro_f1=0.72,
        predicted_class_counts={
            "planet_like": 2,
            "eclipse_contact": 3,
            "smooth_variable": 4,
            "other": 91,
        },
    )
    assert transfer["passed"]

    labels = pd.DataFrame(
        {
            "tic": np.arange(200),
            "human_label": (
                ["planet_like"] * 50
                + ["wide_transit_like"] * 50
                + ["eclipsing_binary_or_pceb"] * 50
                + ["stellar_variability"] * 50
            ),
            "source_kind": "real_candidate",
        }
    )
    metrics = {
        "balanced_accuracy": 0.80,
        "ece": 0.05,
        "per_class": {
            "planet_like": {"n": 10, "recall": 0.80},
            "eclipse_contact": {"n": 10, "recall": 0.70},
            "smooth_variable": {"n": 10, "recall": 0.70},
            "other": {"n": 10, "recall": 0.90},
        },
    }
    readiness = teacher_v2_readiness(labels, test_metrics=metrics)
    assert readiness["teacher_output_count"] == 5
    assert readiness["teacher_classes"][3] == "broad_dip"
    assert readiness["student_ready"]


def test_inference_rows_have_no_targets_and_native_export_accepts_inference_flag() -> None:
    candidates = pd.DataFrame(
        {
            "review_id": ["a"],
            "tic": [123],
            "period_d": [1.0],
            "t0_bjd": [2459000.0],
            "duration_min": [10.0],
        }
    )
    rows = prepare_inference_rows(candidates)
    assert rows.loc[0, "morphology_target_index"] == -1
    assert rows.loc[0, "native_group_path"] == "targets/0000000000000123"
    mask = _native_input_mask(pd.DataFrame({"native_input_include": [True, False]}))
    assert mask.tolist() == [True, False]


def test_s57_rollout_requires_completed_enriched_s56_batch() -> None:
    labels = pd.DataFrame(
        {
            "selection_bucket": (
                ["planet_preserve"] * 300
                + ["stratified_control"] * 100
                + ["other_bucket"] * 600
            ),
            "human_label": (
                ["planet_like"] * 12
                + ["uncertain"] * 288
                + ["planet_like"]
                + ["uncertain"] * 99
                + ["uncertain"] * 600
            ),
        }
    )
    gate = sector_rollout_readiness(
        labels, transfer_gate_passed=True, product_qa_passed=True
    )
    assert gate["s57_plus_ready"]
    assert gate["enrichment_ratio"] == 4.0
