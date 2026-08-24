from __future__ import annotations

import hashlib

import numpy as np
import pytest

from twirl.models.fm0.representation_health import (
    _encode,
    _validate_authority_against_manifest,
    binned_flux_baseline_features,
    fit_train_pca_baseline,
    fit_train_scalar_standardizer,
    load_observation_sector_authority,
    query_sector_excluded_retrieval_summary,
    robust_scalar_baseline_features,
    source_paired_trained_minus_random_summary,
    summarize_projection_spectrum,
    transform_train_pca_baseline,
    transform_train_scalar_baseline,
)


def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sample(*, length: int = 128) -> dict[str, np.ndarray]:
    flux = np.stack(
        (
            np.linspace(-1.0, 1.0, length),
            np.linspace(1.0, -1.0, length),
        )
    ).astype(np.float32)
    reconstruction = np.zeros_like(flux, dtype=bool)
    reconstruction[:, 12:20] = True
    return {
        "flux": flux,
        "flux_valid": np.ones_like(flux, dtype=bool),
        "flux_error": np.ones((2, length), dtype=np.float32),
        "error_valid": np.ones((2, length), dtype=bool),
        "local_time_cadences": np.arange(length, dtype=np.float32),
        "delta_time_cadences": np.r_[
            np.float32(0.0), np.ones(length - 1, dtype=np.float32)
        ],
        "time_valid": np.ones(length, dtype=bool),
        "segment_boundary": np.r_[True, np.zeros(length - 1, dtype=bool)],
        "view_present": np.ones(2, dtype=bool),
        "temporal_mask": np.any(reconstruction, axis=0),
        "reconstruction_mask": reconstruction,
    }


def test_sector_authority_is_checksum_bound_and_rejects_sealed_rows(tmp_path) -> None:
    authority = tmp_path / "observation_sector.csv"
    authority.write_text(
        "observation_key,leakage_component_id,source_partition,sector\n"
        "train-a,component-a,poc_train,56\n"
        "dev-a,component-a,poc_development,57\n",
        encoding="utf-8",
    )
    rows, summary = load_observation_sector_authority(
        authority, expected_sha256=_sha256(authority)
    )
    assert rows["dev-a"]["sector"] == "57"
    assert summary["sealed_test_rows"] == 0
    assert summary["model_visible"] is False
    with pytest.raises(ValueError, match="differs from its SHA-256"):
        load_observation_sector_authority(authority, expected_sha256="0" * 64)

    sealed = tmp_path / "sealed.csv"
    sealed.write_text(
        "observation_key,leakage_component_id,source_partition,sector\n"
        "sealed-a,component-z,poc_sealed_test,58\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="zero sealed-test rows"):
        load_observation_sector_authority(sealed, expected_sha256=_sha256(sealed))

    manifest = [
        {
            "observation_key": "train-a",
            "leakage_component_id": "component-a",
            "source_partition": "poc_train",
        },
        {
            "observation_key": "dev-a",
            "leakage_component_id": "component-a",
            "source_partition": "poc_development",
        },
        {
            "observation_key": "sealed-a",
            "leakage_component_id": "component-z",
            "source_partition": "poc_sealed_test",
        },
    ]
    binding = _validate_authority_against_manifest(manifest, rows)
    assert binding["authority_exactly_covers_manifest_train_development"] is True
    assert binding["sealed_manifest_rows_not_in_authority"] == 1
    with pytest.raises(ValueError, match="exactly cover"):
        _validate_authority_against_manifest(manifest, {"train-a": rows["train-a"]})


def test_query_sector_exclusion_removes_all_same_sector_candidates() -> None:
    embeddings = np.asarray(
        [
            [1.0, 0.0],  # source a, sector 56
            [0.8, 0.6],  # source a, sector 57
            [1.0, 0.0],  # unrelated but identical, sector 56
            [-1.0, 0.0],  # unrelated, sector 57
        ]
    )
    summary = query_sector_excluded_retrieval_summary(
        embeddings,
        ["a", "a", "b", "c"],
        [56, 57, 56, 57],
    )
    assert summary["status"] == "available"
    assert summary["query_sector_excluded"] is True
    assert summary["n_structurally_eligible_queries"] == 2
    assert summary["n_evaluated_queries"] == 2
    assert summary["top1_same_component_retrieval"] == 1.0
    assert summary["top1_source_clustered_95_interval"]["n_source_clusters"] == 1

    unavailable = query_sector_excluded_retrieval_summary(
        embeddings[:2], ["a", "b"], [56, 57]
    )
    assert unavailable["status"] == "not_available"
    assert unavailable["reason"] == "no_components_repeated_across_sectors"
    with pytest.raises(ValueError, match="must be integers"):
        query_sector_excluded_retrieval_summary(
            embeddings[:2], ["a", "a"], [56.5, 57]
        )


def test_source_paired_trained_minus_random_uses_identical_source_rows() -> None:
    trained = np.asarray([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]])
    random = np.ones((3, 2), dtype=np.float64)
    summary = source_paired_trained_minus_random_summary(
        trained,
        trained,
        random,
        random,
        cluster_ids=["a", "b", "c"],
    )
    assert summary["n_source_pairs"] == 3
    assert summary["trained_paired_minus_unrelated_mean"] > 1.0
    assert summary["random_paired_minus_unrelated_mean"] == 0.0
    interval = summary["trained_minus_random_source_clustered_95_interval"]
    assert interval["n_source_clusters"] == 3
    assert interval["lower"] > 0.0


def test_projection_spectrum_reports_rank_and_full_singular_values() -> None:
    summary = summarize_projection_spectrum(
        np.diag([3.0, 1.0, 0.0]), np.asarray([1.0, 2.0, 2.0])
    )
    assert summary["numerical_rank"] == 2
    assert summary["singular_values"] == [3.0, 1.0, 0.0]
    assert 1.0 < summary["effective_rank"] < 2.0
    assert summary["condition_number_on_numerical_support"] == 3.0
    assert summary["bias_l2_norm"] == 3.0


def test_baseline_features_never_read_artificially_masked_flux() -> None:
    first = _sample()
    second = {name: value.copy() for name, value in first.items()}
    first["flux"][:, 12:20] = 1.0e6
    second["flux"][:, 12:20] = -1.0e6
    np.testing.assert_array_equal(
        binned_flux_baseline_features(first, n_time_bins=8),
        binned_flux_baseline_features(second, n_time_bins=8),
    )
    np.testing.assert_array_equal(
        robust_scalar_baseline_features(first),
        robust_scalar_baseline_features(second),
    )


def test_train_fit_pca_and_scalar_baselines_are_deterministic_and_finite() -> None:
    rng = np.random.default_rng(17)
    training = rng.normal(size=(12, 8))
    evaluation = rng.normal(size=(5, 8))
    pca_fit_a, pca_metadata_a = fit_train_pca_baseline(
        training, n_components=4
    )
    pca_fit_b, pca_metadata_b = fit_train_pca_baseline(
        training, n_components=4
    )
    assert pca_metadata_a == pca_metadata_b
    transformed_a = transform_train_pca_baseline(evaluation, pca_fit_a)
    transformed_b = transform_train_pca_baseline(evaluation, pca_fit_b)
    np.testing.assert_array_equal(transformed_a, transformed_b)
    assert transformed_a.shape == (5, 4)
    assert np.all(np.isfinite(transformed_a))
    assert pca_metadata_a["fit_partition"] == "poc_train"
    assert pca_metadata_a["development_fit_rows"] == 0
    assert pca_metadata_a["sealed_test_fit_rows"] == 0

    scalar_fit, scalar_metadata = fit_train_scalar_standardizer(training)
    scalar = transform_train_scalar_baseline(evaluation, scalar_fit)
    assert scalar.shape == evaluation.shape
    assert np.all(np.isfinite(scalar))
    assert scalar_metadata["fit_partition"] == "poc_train"
    assert scalar_metadata["development_fit_rows"] == 0
    assert scalar_metadata["sealed_test_fit_rows"] == 0


def test_encode_collects_preprojection_and_postprojection_representations() -> None:
    torch = pytest.importorskip("torch")

    class TinyEncoder(torch.nn.Module):
        def forward(self, batch):
            hidden = batch["flux"].mean(dim=-1)
            return {
                "h_window": hidden,
                "z_window": 2.0 * hidden,
                "reconstruction": torch.zeros_like(batch["flux"]),
            }

    representations, reconstruction = _encode(
        model=TinyEncoder(),
        samples=[_sample(), _sample()],
        device=torch.device("cpu"),
        batch_size=2,
        with_reconstruction=True,
        huber_delta=0.01,
    )
    assert set(representations) == {"h_window", "z_window"}
    assert representations["h_window"].shape == (2, 2)
    np.testing.assert_allclose(
        representations["z_window"], 2.0 * representations["h_window"]
    )
    assert reconstruction is not None
    assert reconstruction["masked_valid_target_count"] > 0
    assert reconstruction["zero_prediction_masked_huber_mean"] is not None
