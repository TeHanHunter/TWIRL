from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from twirl.vetting.teacher_v2_active_learning import (
    DEFAULT_ENRICHMENT_QUOTAS,
    build_existing_teacher_enrichment_batch,
    combine_existing_ranker_scores,
    verify_enrichment_batch,
)


def _score_tables(n_tics: int = 1200) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(56)
    n_rows = 2 * n_tics
    tic = np.repeat(np.arange(10_000, 10_000 + n_tics), 2)
    peak = np.tile([1, 2], n_tics)
    base = pd.DataFrame(
        {
            "review_id": [f"candidate-{value}-{rank}" for value, rank in zip(tic, peak)],
            "tic": tic,
            "sector": 56,
            "period_d": 0.2 + rng.random(n_rows) * 10.0,
            "t0_bjd": 2_459_000.0 + rng.random(n_rows),
            "duration_min": 3.0 + rng.random(n_rows) * 40.0,
            "tmag": 16.0 + rng.random(n_rows) * 5.0,
            "sde_max": 3.0 + rng.random(n_rows) * 50.0,
            "rep_peak_rank": peak,
            "cam": 1,
            "ccd": 1,
            "source_kind": "real_candidate",
            "is_injected_row": False,
        }
    )
    compact = base.copy()
    compact["p_compact_transit"] = rng.random(n_rows)
    compact["std_p_compact_transit"] = rng.random(n_rows) * 0.1
    compact["model_profile"] = "shape_plus_raw_chronology"
    compact["model_version"] = "teacher-v2"
    for label in ("planet_like", "eclipse_contact", "smooth_variable", "other"):
        compact[f"p_{label}"] = 0.25
        compact[f"std_p_{label}"] = 0.01
    compact["p_preserve"] = 0.5
    compact["std_p_preserve"] = 0.01

    morphology = base.copy()
    probabilities = rng.dirichlet(np.ones(4), size=n_rows)
    for index, label in enumerate(("planet_like", "eclipse_contact", "smooth_variable", "other")):
        morphology[f"p_{label}"] = probabilities[:, index]
        morphology[f"std_p_{label}"] = rng.random(n_rows) * 0.1
    morphology["p_preserve"] = rng.random(n_rows)
    morphology["std_p_preserve"] = rng.random(n_rows) * 0.1
    morphology["p_compact_transit"] = rng.random(n_rows)
    morphology["std_p_compact_transit"] = rng.random(n_rows) * 0.1
    morphology["model_profile"] = "shape_plus_periodogram_bls"
    morphology["model_version"] = "teacher-v2"
    return compact, morphology


def test_existing_rankers_join_exactly_and_keep_intended_outputs() -> None:
    compact, morphology = _score_tables(20)
    combined = combine_existing_ranker_scores(compact, morphology)
    assert len(combined) == len(compact)
    assert set(combined["compact_model_profile"]) == {"shape_plus_raw_chronology"}
    assert set(combined["morphology_model_profile"]) == {"shape_plus_periodogram_bls"}
    assert np.allclose(combined["p_eclipse_contact"], morphology["p_eclipse_contact"])
    assert np.allclose(combined["p_compact_transit"], compact["p_compact_transit"])
    assert combined["morphology_entropy"].between(0.0, 1.0).all()


def test_existing_ranker_join_rejects_ephemeris_mismatch() -> None:
    compact, morphology = _score_tables(20)
    morphology.loc[0, "period_d"] += 0.01
    with pytest.raises(ValueError, match="period_d"):
        combine_existing_ranker_scores(compact, morphology)


def test_enrichment_batch_is_deterministic_blinded_and_unique() -> None:
    compact, morphology = _score_tables()
    scores = combine_existing_ranker_scores(compact, morphology)
    first = build_existing_teacher_enrichment_batch(scores, sector=56, batch_index=0)
    second = build_existing_teacher_enrichment_batch(scores, sector=56, batch_index=0)
    queue, overlap, hidden, summary = first
    assert queue["candidate_key"].tolist() == second[0]["candidate_key"].tolist()
    assert len(queue) == DEFAULT_ENRICHMENT_QUOTAS.total == 1000
    assert queue["tic"].nunique() == 1000
    assert len(overlap) == 100
    assert summary["training_performed"] is False
    assert hidden["selection_bucket"].value_counts().to_dict() == {
        "compact_transit": 400,
        "eclipse_contact": 300,
        "smooth_variable": 100,
        "model_disagreement": 100,
        "stratified_control": 100,
    }
    assert not any(
        column.startswith(("p_", "std_p_", "member_", "selection_", "model_"))
        for column in queue
    )
    assert verify_enrichment_batch(queue, overlap, hidden)["passed"]


def test_enrichment_batch_excludes_prior_tics() -> None:
    compact, morphology = _score_tables(1300)
    scores = combine_existing_ranker_scores(compact, morphology)
    excluded = pd.DataFrame({"tic": np.arange(10_000, 10_050)})
    queue, _, _, summary = build_existing_teacher_enrichment_batch(
        scores,
        sector=56,
        batch_index=0,
        excluded_tables=[excluded],
    )
    assert set(queue["tic"]).isdisjoint(set(excluded["tic"]))
    assert summary["n_excluded_tics"] == 50
