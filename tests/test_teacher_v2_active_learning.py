from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from twirl.vetting.teacher_v2_active_learning import (
    DEFAULT_ENRICHMENT_QUOTAS,
    EnrichmentQuotas,
    build_existing_teacher_enrichment_batch,
    build_mixed_existing_teacher_enrichment_batch,
    combine_existing_ranker_scores,
    verify_enrichment_batch,
    write_existing_teacher_enrichment_batch,
    write_mixed_existing_teacher_enrichment_batch,
)


def _score_tables(
    n_tics: int = 1200,
    *,
    sector: int = 56,
    tic_start: int = 10_000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(56 + sector)
    n_rows = 2 * n_tics
    tic = np.repeat(np.arange(tic_start, tic_start + n_tics), 2)
    peak = np.tile([1, 2], n_tics)
    base = pd.DataFrame(
        {
            "review_id": [f"candidate-{value}-{rank}" for value, rank in zip(tic, peak)],
            "tic": tic,
            "sector": sector,
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


def test_written_queue_candidate_keys_survive_csv_roundtrip(tmp_path) -> None:
    compact, morphology = _score_tables()
    compact_path = tmp_path / "compact.parquet"
    morphology_path = tmp_path / "morphology.parquet"
    compact.to_parquet(compact_path, index=False)
    morphology.to_parquet(morphology_path, index=False)
    summary = write_existing_teacher_enrichment_batch(
        compact_scores_path=compact_path,
        morphology_scores_path=morphology_path,
        out_dir=tmp_path / "batch",
        sector=56,
        batch_index=0,
    )
    assert summary["n_rows"] == 1000


def test_mixed_batch_is_sector_balanced_and_globally_deduplicated() -> None:
    compact56, morphology56 = _score_tables(180, sector=56, tic_start=10_000)
    compact57, morphology57 = _score_tables(180, sector=57, tic_start=10_100)
    scores = {
        56: combine_existing_ranker_scores(compact56, morphology56),
        57: combine_existing_ranker_scores(compact57, morphology57),
    }
    quotas = EnrichmentQuotas(
        compact_transit=20,
        eclipse_contact=15,
        smooth_variable=10,
        model_disagreement=10,
        stratified_control=5,
    )
    queue, overlap, hidden, summary = build_mixed_existing_teacher_enrichment_batch(
        scores,
        batch_index=0,
        per_sector_quotas=quotas,
    )
    assert len(queue) == 120
    assert queue["tic"].nunique() == 120
    assert queue["sector"].value_counts().to_dict() == {56: 60, 57: 60}
    assert len(overlap) == 100
    assert summary["cross_sector_tic_deduplication"] is True
    assert hidden["selection_bucket"].value_counts().to_dict() == {
        "compact_transit": 40,
        "eclipse_contact": 30,
        "smooth_variable": 20,
        "model_disagreement": 20,
        "stratified_control": 10,
    }
    assert not any(
        column.startswith(("p_", "std_p_", "member_", "selection_", "model_"))
        for column in queue
    )


def test_written_mixed_queue_survives_csv_roundtrip(tmp_path) -> None:
    paths: dict[int, tuple[Path, Path]] = {}
    for sector, tic_start in ((56, 10_000), (57, 20_000)):
        compact, morphology = _score_tables(180, sector=sector, tic_start=tic_start)
        compact_path = tmp_path / f"compact_{sector}.parquet"
        morphology_path = tmp_path / f"morphology_{sector}.parquet"
        compact.to_parquet(compact_path, index=False)
        morphology.to_parquet(morphology_path, index=False)
        paths[sector] = (compact_path, morphology_path)
    summary = write_mixed_existing_teacher_enrichment_batch(
        sector_score_paths=paths,
        out_dir=tmp_path / "mixed",
        batch_index=0,
        per_sector_quotas=EnrichmentQuotas(
            compact_transit=20,
            eclipse_contact=15,
            smooth_variable=10,
            model_disagreement=10,
            stratified_control=5,
        ),
    )
    queue = pd.read_csv(tmp_path / "mixed/review_queue_1k.csv")
    assert summary["n_rows"] == 120
    assert queue["sector"].value_counts().to_dict() == {56: 60, 57: 60}
