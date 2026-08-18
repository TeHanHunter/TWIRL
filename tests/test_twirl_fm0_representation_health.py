from __future__ import annotations

import numpy as np

from twirl.models.fm0.representation_health import (
    paired_similarity_summary,
    same_component_retrieval_summary,
    source_clustered_mean_interval,
    summarize_embedding_matrix,
)


def test_embedding_summary_detects_rank_and_constant_dimension() -> None:
    embeddings = np.asarray(
        [[0.0, 0.0, 2.0], [1.0, 0.0, 2.0], [2.0, 0.0, 2.0]],
        dtype=np.float64,
    )
    summary = summarize_embedding_matrix(embeddings)
    assert summary["embedding_dim"] == 3
    assert summary["zero_or_constant_dimensions"] == 2
    assert summary["effective_rank"] == 1.0
    assert summary["leading_principal_component_share"] == 1.0


def test_paired_similarity_exceeds_deterministic_unrelated_control() -> None:
    left = np.asarray([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]])
    right = np.asarray([[0.99, 0.01], [0.01, 0.99], [-0.99, 0.01]])
    summary = paired_similarity_summary(left, right, cluster_ids=["a", "b", "c"])
    assert summary["n_pairs"] == 3
    assert summary["paired_cosine_mean"] > 0.99
    assert summary["paired_minus_unrelated_mean"] > 0.5
    assert summary["paired_minus_unrelated_source_clustered_95_interval"]["lower"] > 0.5


def test_same_component_retrieval_excludes_single_visit_components() -> None:
    embeddings = np.asarray(
        [[1.0, 0.0], [0.9, 0.1], [0.0, 1.0], [0.0, 0.9]], dtype=np.float64
    )
    summary = same_component_retrieval_summary(
        embeddings, ["a", "a", "b", "single"]
    )
    assert summary["status"] == "available"
    assert summary["n_repeated_component_queries"] == 2
    assert summary["top1_same_component_retrieval"] == 1.0
    assert summary["top1_source_clustered_95_interval"]["lower"] == 1.0


def test_source_clustered_interval_uses_components_not_raw_row_count() -> None:
    interval = source_clustered_mean_interval(
        np.asarray([1.0, 1.0, 0.0]), ["a", "a", "b"], seed=7, n_bootstrap=100
    )
    assert interval["n_source_clusters"] == 2
    assert interval["mean"] == 0.5
