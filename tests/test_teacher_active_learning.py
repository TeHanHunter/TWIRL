from __future__ import annotations

import numpy as np
import pandas as pd

from twirl.search.a2v1_bls_contract import (
    A2V1_TEACHER_BLS_SEARCH_CONTRACT,
    approved_a2v1_teacher_bls_config,
    bls_config_sha256,
)
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
from twirl.vetting.teacher_candidates import (
    filter_tier1_eligible_candidates,
    normalize_a2v1_peak_candidates,
    validate_quality_bound_bls_evidence,
    validate_tier1_enrichment_gate,
)


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


def _tier1_target_eligibility() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "tier1_contract_version": ["tier1-v1"] * 3,
            "tier1_config_name": ["active-pair-v1"] * 3,
            "tier1_scope": ["active_search_pair"] * 3,
            "sector_tic_key": [
                "s0056-tic0000000000000010",
                "s0056-tic0000000000000020",
                "s0057-tic0000000000000010",
            ],
            "sector": [56, 56, 57],
            "tic": [10, 20, 10],
            "tier1_target_qa_status": ["pass", "review", "fail"],
            "tier1_target_qa_reasons": [
                np.nan,
                "scatter_absolute_high",
                "finite_flux",
            ],
            "tier1_target_qa_pass": ["True", "False", "False"],
        }
    )


def _tier1_gate_summary(eligibility_sha256: str = "a" * 64) -> dict:
    gates = {
        name: {"status": "pass"}
        for name in (
            "cadence_reference_prerequisite",
            "injection_source_parity_prerequisite",
            "tier0_prerequisite",
            "population_scatter",
            "cadence_and_finite_data",
            "aperture_outliers",
            "fixed_injection_preservation",
            "independent_extraction",
        )
    }
    gates["cadence_reference_prerequisite"].update(
        {
            "table_sha256": "d" * 64,
            "manifest_sha256": "e" * 64,
            "cadence_authority": "qlp_cam_quat",
            "quality_authority": "spoc_and_qlp_quality_flags",
        }
    )
    return {
        "contract_version": "tier1-v1",
        "config_name": "active-pair-v1",
        "qa_tier": "tier1_bounded_enrichment_qa",
        "scope": "active_search_pair",
        "sector": 56,
        "status": "pass",
        "passed": True,
        "enrichment_ready": True,
        "science_ready": False,
        "promotion_enabled": False,
        "apertures": ["DET_FLUX_ADP_SML", "DET_FLUX_ADP"],
        "target_qa": {
            "observation_key": ["sector", "tic"],
            "candidate_teacher_filter": "tier1_target_qa_pass == True",
            "n_pass": 1,
            "n_review": 1,
            "n_fail": 0,
        },
        "provenance": {
            "compact_lc_sha256": "c" * 64,
            "cadence_reference_sha256": "d" * 64,
            "cadence_reference_manifest_sha256": "e" * 64,
            "output_sha256": {"target_eligibility": eligibility_sha256},
        },
        "gates": gates,
    }


def _tier1_gate_eligibility() -> pd.DataFrame:
    return _tier1_target_eligibility().loc[lambda frame: frame["sector"].eq(56)].copy()


def _quality_bound_bls_evidence() -> tuple[pd.DataFrame, dict]:
    config = approved_a2v1_teacher_bls_config()
    config_sha256 = bls_config_sha256(config)
    peaks = _a2v1_peaks().assign(
        external_quality_policy_contract=(
            "a2v1_bls_internal_or_authoritative_external_quality_v1"
        ),
        cadence_reference_sha256="d" * 64,
        cadence_reference_manifest_sha256="e" * 64,
        bls_search_contract_version=A2V1_TEACHER_BLS_SEARCH_CONTRACT,
        bls_config_sha256=config_sha256,
    )
    summary = {
        "sector": 56,
        "contract_version": ADP_ONLY_CONTRACT_VERSION,
        "bls_search_contract_version": A2V1_TEACHER_BLS_SEARCH_CONTRACT,
        "bls_config_sha256": config_sha256,
        "config": config,
        "external_quality_policy_contract": (
            "a2v1_bls_internal_or_authoritative_external_quality_v1"
        ),
        "compact_lc_sha256": "c" * 64,
        "cadence_reference_sha256": "d" * 64,
        "cadence_reference_manifest_sha256": "e" * 64,
        "cadence_reference_contract_version": "s56_a2v1_cadence_reference_v1",
        "cadence_reference_cadence_authority": "qlp_cam_quat",
        "cadence_reference_quality_authority": "spoc_and_qlp_quality_flags",
        "cadence_reference_source_hashes_sha256": "f" * 64,
        "apertures": ["DET_FLUX_ADP_SML", "DET_FLUX_ADP"],
        "source_product_tag": "A2v1",
        "n_shards": 1,
        "shard_index": 0,
        "n_targets": 2,
        "n_targets_total": 2,
        "n_unique_tics": 2,
        "n_rows": len(peaks),
        "peak_table_sha256": "b" * 64,
    }
    return peaks, summary


def test_quality_bound_bls_evidence_binds_tier1_and_target_coverage() -> None:
    peaks, summary = _quality_bound_bls_evidence()
    audit = validate_quality_bound_bls_evidence(
        peaks,
        summary,
        _tier1_gate_summary(),
        _tier1_gate_eligibility(),
        peak_table_sha256="b" * 64,
        compact_lc_sha256="c" * 64,
    )
    assert audit["status"] == "pass"
    assert audit["n_targets"] == 2
    assert audit["n_tier1_passing_targets"] == 1

    stale = dict(summary)
    stale["cadence_reference_sha256"] = "9" * 64
    with np.testing.assert_raises_regex(ValueError, "does not match the Tier-1"):
        validate_quality_bound_bls_evidence(
            peaks,
            stale,
            _tier1_gate_summary(),
            _tier1_gate_eligibility(),
            peak_table_sha256="b" * 64,
            compact_lc_sha256="c" * 64,
        )

    missing_aperture = peaks.loc[
        ~(
            peaks["tic"].eq(10)
            & peaks["aperture"].eq("DET_FLUX_ADP")
        )
    ].copy()
    incomplete_summary = dict(summary)
    incomplete_summary["n_rows"] = len(missing_aperture)
    with np.testing.assert_raises_regex(ValueError, "lack rank-1 BLS"):
        validate_quality_bound_bls_evidence(
            missing_aperture,
            incomplete_summary,
            _tier1_gate_summary(),
            _tier1_gate_eligibility(),
            peak_table_sha256="b" * 64,
            compact_lc_sha256="c" * 64,
        )

    custom_search = dict(summary)
    custom_search["config"] = dict(summary["config"])
    custom_search["config"]["n_periods"] = 1_000
    custom_search["bls_config_sha256"] = bls_config_sha256(
        custom_search["config"]
    )
    with np.testing.assert_raises_regex(ValueError, "approved search config"):
        validate_quality_bound_bls_evidence(
            peaks,
            custom_search,
            _tier1_gate_summary(),
            _tier1_gate_eligibility(),
            peak_table_sha256="b" * 64,
            compact_lc_sha256="c" * 64,
        )


def test_tier1_enrichment_gate_binds_passed_summary_to_eligibility() -> None:
    audit = validate_tier1_enrichment_gate(
        _tier1_gate_summary(),
        _tier1_gate_eligibility(),
        target_eligibility_sha256="a" * 64,
        compact_lc_sha256="c" * 64,
    )

    assert audit["passed"]
    assert audit["enrichment_ready"]
    assert audit["sector"] == 56
    assert audit["contract_version"] == "tier1-v1"
    assert audit["target_eligibility_sha256"] == "a" * 64
    assert audit["compact_lc_sha256"] == "c" * 64


def test_tier1_enrichment_gate_rejects_global_failure_or_nonready_run() -> None:
    failed = _tier1_gate_summary()
    failed["passed"] = False
    with np.testing.assert_raises_regex(ValueError, "passed is not true"):
        validate_tier1_enrichment_gate(
            failed,
            _tier1_gate_eligibility(),
            target_eligibility_sha256="a" * 64,
            compact_lc_sha256="c" * 64,
        )

    nonready = _tier1_gate_summary()
    nonready["enrichment_ready"] = False
    with np.testing.assert_raises_regex(ValueError, "enrichment_ready is not true"):
        validate_tier1_enrichment_gate(
            nonready,
            _tier1_gate_eligibility(),
            target_eligibility_sha256="a" * 64,
            compact_lc_sha256="c" * 64,
        )

    nested_failure = _tier1_gate_summary()
    nested_failure["gates"]["independent_extraction"]["status"] = "fail"
    with np.testing.assert_raises_regex(
        ValueError, "nested gate independent_extraction"
    ):
        validate_tier1_enrichment_gate(
            nested_failure,
            _tier1_gate_eligibility(),
            target_eligibility_sha256="a" * 64,
            compact_lc_sha256="c" * 64,
        )


def test_tier1_enrichment_gate_rejects_unbound_or_incompatible_eligibility() -> None:
    with np.testing.assert_raises_regex(ValueError, "SHA-256 does not match"):
        validate_tier1_enrichment_gate(
            _tier1_gate_summary(),
            _tier1_gate_eligibility(),
            target_eligibility_sha256="b" * 64,
            compact_lc_sha256="c" * 64,
        )

    with np.testing.assert_raises_regex(ValueError, "compact LC SHA-256"):
        validate_tier1_enrichment_gate(
            _tier1_gate_summary(),
            _tier1_gate_eligibility(),
            target_eligibility_sha256="a" * 64,
            compact_lc_sha256="d" * 64,
        )

    incompatible = _tier1_gate_eligibility()
    incompatible["tier1_config_name"] = "different-config"
    with np.testing.assert_raises_regex(ValueError, "config_name does not match"):
        validate_tier1_enrichment_gate(
            _tier1_gate_summary(),
            incompatible,
            target_eligibility_sha256="a" * 64,
            compact_lc_sha256="c" * 64,
        )


def test_tier1_filter_uses_exact_sector_tic_key_and_keeps_provenance() -> None:
    candidates = pd.DataFrame(
        {
            "review_id": ["a", "b", "c"],
            "sector": [56, 56, 56],
            "tic": [10, 10, 20],
        }
    )
    filtered, summary = filter_tier1_eligible_candidates(
        candidates, _tier1_target_eligibility()
    )

    assert filtered["review_id"].tolist() == ["a", "b"]
    assert filtered["tier1_target_qa_pass"].all()
    assert filtered["tier1_target_qa_status"].eq("pass").all()
    assert filtered["tier1_target_qa_reasons"].eq("").all()
    assert filtered["tier1_contract_version"].eq("tier1-v1").all()
    assert summary["join_key"] == ["sector", "tic"]
    assert summary["n_candidates_before"] == 3
    assert summary["n_candidates_after"] == 2
    assert summary["n_candidate_tics_before"] == 2
    assert summary["n_candidate_tics_after"] == 1
    assert summary["n_candidate_observations_before"] == 2
    assert summary["n_candidate_observations_after"] == 1
    assert summary["candidate_status_counts"] == {"pass": 2, "review": 1}


def test_tier1_candidate_filter_rejects_duplicate_eligibility_keys() -> None:
    candidates = pd.DataFrame({"sector": [56], "tic": [10]})
    eligibility = _tier1_target_eligibility()
    eligibility = pd.concat([eligibility, eligibility.iloc[[0]]], ignore_index=True)

    with np.testing.assert_raises_regex(ValueError, "unique \\(sector, TIC\\)"):
        filter_tier1_eligible_candidates(candidates, eligibility)


def test_tier1_candidate_filter_rejects_incomplete_exact_key_coverage() -> None:
    candidates = pd.DataFrame({"sector": [58], "tic": [10]})

    with np.testing.assert_raises_regex(ValueError, "does not completely cover"):
        filter_tier1_eligible_candidates(candidates, _tier1_target_eligibility())


def test_tier1_candidate_filter_rejects_inconsistent_status_and_pass_flag() -> None:
    candidates = pd.DataFrame({"sector": [56], "tic": [10]})
    eligibility = _tier1_target_eligibility()
    eligibility.loc[eligibility["tic"].eq(10), "tier1_target_qa_pass"] = "False"

    with np.testing.assert_raises_regex(ValueError, "status and pass flag"):
        filter_tier1_eligible_candidates(candidates, eligibility)


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
