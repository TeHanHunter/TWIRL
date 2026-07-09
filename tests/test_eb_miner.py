from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scripts.stage5_validation.train_score_s56_eb_miner import _resolve_score_profile
from twirl.vetting.adp_only import ADP_ONLY_CONTRACT_VERSION
from twirl.vetting.eb_miner import (
    EBTrainingConfig,
    EB_TARGET_NEGATIVE,
    EB_TARGET_POSITIVE,
    aggregate_ensemble_scores,
    build_candidate_scoring_pool,
    build_eb_miner_training_table,
    evaluate_eb_miner_release_gate,
    select_eb_priority_queue,
)
from twirl.vetting.recovery50_teacher import metadata_feature_columns, read_table


def _write_csv(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_eb_training_table_uses_latest_real_labels_and_excludes_injections(tmp_path: Path) -> None:
    queue = _write_csv(
        tmp_path / "queue.csv",
        [
            {
                "row_id": 0,
                "tic": 1,
                "sector": 56,
                "source_kind": "real_candidate",
                "review_id": "real:1",
                "period_d": 1.0,
                "t0_bjd": 2459000.0,
                "duration_min": 10,
                "source_bucket": "real_test",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            },
            {
                "row_id": 1,
                "tic": 2,
                "sector": 56,
                "source_kind": "real_candidate",
                "review_id": "real:2",
                "period_d": 2.0,
                "t0_bjd": 2459001.0,
                "duration_min": 12,
                "source_bucket": "real_test",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            },
            {
                "row_id": 2,
                "tic": 3,
                "sector": 56,
                "source_kind": "real_candidate",
                "review_id": "real:3",
                "period_d": 3.0,
                "t0_bjd": 2459002.0,
                "duration_min": 14,
                "source_bucket": "real_test",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            },
            {
                "row_id": 3,
                "tic": 4,
                "sector": 56,
                "source_kind": "injection_recovery",
                "review_id": "inj:4",
                "period_d": 4.0,
                "t0_bjd": 2459003.0,
                "duration_min": 16,
                "source_bucket": "inj_test",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            },
        ],
    )
    labels = _write_csv(
        tmp_path / "labels.csv",
        [
            {"row_id": 0, "label": "uncertain", "updated_utc": "2026-07-01T00:00:00Z"},
            {"row_id": 0, "label": "eclipsing_binary_or_pceb", "updated_utc": "2026-07-02T00:00:00Z"},
            {"row_id": 1, "label": "uncertain", "updated_utc": "2026-07-01T00:00:00Z"},
            {"row_id": 2, "label": "skip", "updated_utc": "2026-07-01T00:00:00Z"},
            {"row_id": 3, "label": "eclipsing_binary_or_pceb", "updated_utc": "2026-07-01T00:00:00Z"},
        ],
    )

    summary = build_eb_miner_training_table(
        joined_tables=(),
        queue_label_pairs=((queue, labels, "toy"),),
        out_dir=tmp_path / "out",
        cfg=EBTrainingConfig(min_positive_rows=1, max_uncertain_negatives=10, max_negatives_per_label=10),
    )
    audit = read_table(Path(summary["outputs"]["audit"]))
    training = read_table(Path(summary["outputs"]["training_table"]))

    assert audit["source_kind"].eq("real_candidate").all()
    assert "inj:4" not in set(audit["review_id"])
    latest = audit.set_index("review_id")
    assert latest.loc["real:1", "human_label"] == "eclipsing_binary_or_pceb"
    assert latest.loc["real:1", "eb_miner_target"] == EB_TARGET_POSITIVE
    assert latest.loc["real:2", "eb_miner_target"] == EB_TARGET_NEGATIVE
    assert latest.loc["real:3", "eb_miner_include"] in {False, 0}
    assert set(training["eb_miner_target"]) == {EB_TARGET_POSITIVE, EB_TARGET_NEGATIVE}


def test_candidate_scoring_pool_excludes_prior_tics_and_review_ids(tmp_path: Path) -> None:
    primary = _write_csv(
        tmp_path / "ranker.csv",
        [
            {
                "tic": 10,
                "peak_rank": 1,
                "sector": 56,
                "period_d": 1.0,
                "t0_bjd": 2459000.0,
                "duration_min": 10,
                "sde": 20,
                "aperture": "DET_FLUX_ADP_SML",
                "bls_search_branch": "current_adp",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            },
            {
                "tic": 11,
                "peak_rank": 1,
                "sector": 56,
                "period_d": 2.0,
                "t0_bjd": 2459001.0,
                "duration_min": 11,
                "sde": 19,
                "depth": 0.10,
                "depth_snr": 8.0,
                "log_power": 4.0,
                "aperture": "DET_FLUX_ADP_SML",
                "bls_search_branch": "current_adp",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            },
            {
                "tic": 11,
                "peak_rank": 1,
                "sector": 56,
                "period_d": 2.01,
                "t0_bjd": 2459001.0,
                "duration_min": 12,
                "sde": 17,
                "depth": 0.05,
                "depth_snr": 6.0,
                "log_power": 3.0,
                "aperture": "DET_FLUX_ADP",
                "bls_search_branch": "current_adp",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            },
        ],
    )
    fallback = _write_csv(
        tmp_path / "fallback.csv",
        [
            {
                "tic": 12,
                "sector": 56,
                "period_d": 3.0,
                "t0_bjd": 2459002.0,
                "duration_min": 12,
                "sde_max": 18,
                "sde": 18,
                "peak_rank": 1,
                "aperture": "DET_FLUX_ADP",
                "bls_search_branch": "current_adp",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            }
        ],
    )
    exclude = _write_csv(
        tmp_path / "exclude.csv",
        [
            {"tic": 10, "source_kind": "real_candidate", "review_id": "real:10"},
            {"tic": 12, "source_kind": "real_candidate", "review_id": "real:12"},
        ],
    )

    summary = build_candidate_scoring_pool(
        primary_candidates=primary,
        fallback_candidates=fallback,
        exclude_tables=(exclude,),
        out_dir=tmp_path / "pool",
    )
    pool = read_table(Path(summary["outputs"]["candidate_pool"]))

    assert set(pool["tic"]) == {11}
    small_row = pool.loc[pool["aperture"].eq("DET_FLUX_ADP_SML")].iloc[0]
    assert small_row["review_id"] == "real:11:adp_bls:adp_sml:peak:1"
    assert small_row["source_kind"] == "real_candidate"
    assert small_row["adp_sml_sde"] == 19
    assert small_row["adp_sde"] == 17
    assert small_row["aperture_period_relation"] == "exact"


def test_candidate_scoring_pool_rejects_canonical_bls_table(tmp_path: Path) -> None:
    bad = _write_csv(
        tmp_path / "canonical.csv",
        [
            {
                "tic": 10,
                "peak_rank": 1,
                "period_d": 1.0,
                "t0_bjd": 2459000.0,
                "duration_min": 10,
                "aperture": "DET_FLUX_SML",
                "bls_search_branch": "current_adp",
                "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
            }
        ],
    )
    with pytest.raises(ValueError, match="non-ADP apertures"):
        build_candidate_scoring_pool(
            primary_candidates=bad,
            exclude_tables=(),
            out_dir=tmp_path / "pool",
        )


def test_ensemble_aggregation_and_priority_queue_png_names(tmp_path: Path) -> None:
    base = pd.DataFrame(
        {
            "review_id": [f"real:{i}:ranker:1" for i in range(5)],
            "tic": [100, 101, 101, 102, 103],
            "period_d": [1.0, 2.0, 4.0, 3.0, 5.0],
            "t0_bjd": [2459000.0] * 5,
            "duration_min": [10.0] * 5,
            "sector": [56] * 5,
            "source_kind": ["real_candidate"] * 5,
            "source_pool": ["ranker_selected"] * 5,
            "sde_max": [9, 8, 7, 6, 5],
            "cnn_p_eb_pceb": [0.95, 0.90, 0.85, 0.80, 0.70],
        }
    )
    path1 = _write_csv(tmp_path / "score1.csv", base.to_dict(orient="records"))
    path2 = _write_csv(
        tmp_path / "score2.csv",
        base.assign(cnn_p_eb_pceb=base["cnn_p_eb_pceb"] - 0.05).to_dict(orient="records"),
    )
    scored = aggregate_ensemble_scores([path1, path2])
    queue, summary = select_eb_priority_queue(scored, n_review=3, harmonic_probability_floor=0.8)

    assert len(queue) == 3
    assert summary["n_review"] == 3
    assert queue["twirl_vet_sheet_name"].str.endswith(".png").all()
    assert queue["twirl_vet_sheet_pdf_name"].str.endswith(".pdf").all()
    assert queue["source_kind"].eq("real_candidate").all()


def test_eb_release_gate_requires_grouped_heldout_enrichment() -> None:
    good_member = {
        "eb_binary_eval": {
            "validation": {
                "n_positive": 2,
                "tp": 1,
                "fn": 1,
                "average_precision": 0.50,
            },
            "test": {
                "n_positive": 2,
                "tp": 1,
                "fn": 1,
                "average_precision": 0.40,
            },
        }
    }
    passed = evaluate_eb_miner_release_gate(
        {"member_summaries": [good_member, good_member, good_member]}
    )
    failed = evaluate_eb_miner_release_gate(
        {"member_summaries": [{"eb_binary_eval": {}}]}
    )

    assert passed["passed"]
    assert passed["heldout_recall"] == 0.5
    assert not failed["passed"]


def test_auto_score_profile_uses_validation_winner() -> None:
    summary = {"best_profile_by_validation_balanced_accuracy": "cnn_shape_only"}
    assert _resolve_score_profile(summary, "auto") == "cnn_shape_only"
    assert _resolve_score_profile(summary, "cnn_shape_plus_bls") == "cnn_shape_plus_bls"

def test_metadata_allowlist_rejects_truth_and_keeps_search_features() -> None:
    frame = pd.DataFrame(
        {
            "adp_sml_sde": [10.0],
            "adp_sml_depth": [0.1],
            "truth_period_d": [1.0],
            "human_label": ["eclipsing_binary_or_pceb"],
        }
    )

    assert metadata_feature_columns(frame) == ["adp_sml_depth", "adp_sml_sde"]
