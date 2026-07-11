from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from twirl.vetting.harmonic_cnn import (
    HARMONIC_CLASSES,
    HarmonicModelConfig,
    build_grouped_test_and_cv_folds,
    build_harmonic_cnn,
    build_sample_loss_weights,
)


def test_loss_weights_enforce_source_and_other_caps() -> None:
    rows = pd.DataFrame(
        {
            "morphology_target_v1": (
                ["planet_like"] * 12
                + ["planet_like"] * 100
                + ["other"] * 30
                + ["other"] * 90
                + ["other"] * 80
            ),
            "morphology_include_v1": [True] * 312,
            "is_injected_row": [False] * 12 + [True] * 100 + [False] * 120 + [True] * 80,
            "human_label": (
                ["planet_like"] * 112
                + ["instrumental_or_systematic"] * 30
                + ["uncertain"] * 90
                + ["uncertain"] * 80
            ),
        }
    )

    weight = build_sample_loss_weights(rows)
    injected = rows["is_injected_row"].to_numpy()
    planet = rows["morphology_target_v1"].eq("planet_like").to_numpy()
    other = rows["morphology_target_v1"].eq("other").to_numpy()
    systematic = rows["human_label"].eq("instrumental_or_systematic").to_numpy()
    flat = rows["human_label"].eq("uncertain").to_numpy()

    assert np.isclose(weight[planet & injected].sum(), weight[planet & ~injected].sum())
    assert np.isclose(weight[other & ~injected & systematic].sum(), weight[other & ~injected & flat].sum())
    assert weight[other & injected].sum() <= 0.25 * weight[other & ~injected].sum() + 1.0e-5


def test_grouped_test_and_cv_splits_keep_host_tics_together() -> None:
    rows = pd.DataFrame(
        {
            "tic": np.repeat(np.arange(100, 160), 2),
            "split_stratum": np.tile(
                ["planet_like", "eclipse_contact", "smooth_variable", "other", "broad_preserve"],
                24,
            ),
        }
    )
    split = build_grouped_test_and_cv_folds(rows)

    joined = pd.concat([rows, split], axis=1)
    for _, group in joined.groupby("tic"):
        assert group["fixed_split"].nunique() == 1
        assert group["cv_fold"].nunique() == 1
    assert set(joined["fixed_split"]) == {"development", "test"}
    assert set(joined.loc[joined["fixed_split"].eq("development"), "cv_fold"]) == set(range(5))
    assert joined.loc[joined["fixed_split"].eq("test"), "cv_fold"].eq(-1).all()


def test_model_has_three_heads_and_native_ragged_masks() -> None:
    torch = pytest.importorskip("torch")
    model = build_harmonic_cnn(HarmonicModelConfig(metadata_dim=9), profile="full_combined")
    batch_size = 2
    batch = {
        "chronology_small": torch.randn(batch_size, 10, 37),
        "chronology_small_mask": torch.ones(batch_size, 10, 37, dtype=torch.bool),
        "chronology_supplemental": torch.randn(batch_size, 5, 37),
        "chronology_supplemental_mask": torch.ones(batch_size, 5, 37, dtype=torch.bool),
        "harmonic_values": torch.randn(batch_size, 7, 7, 37),
        "harmonic_mask": torch.ones(batch_size, 7, 7, 37, dtype=torch.bool),
        "local_values": torch.randn(batch_size, 14, 7, 19),
        "local_mask": torch.ones(batch_size, 14, 7, 19, dtype=torch.bool),
        "periodogram_values": torch.randn(batch_size, 4, 31),
        "periodogram_mask": torch.ones(batch_size, 4, 31, dtype=torch.bool),
        "metadata": torch.randn(batch_size, 9),
    }

    outputs = model(batch)

    assert outputs["morphology_logits"].shape == (batch_size, 4)
    assert outputs["preserve_logits"].shape == (batch_size, 2)
    assert outputs["harmonic_logits"].shape == (batch_size, len(HARMONIC_CLASSES))


def test_fully_masked_sequence_and_harmonic_branches_are_exactly_zero() -> None:
    torch = pytest.importorskip("torch")
    model = build_harmonic_cnn(HarmonicModelConfig(metadata_dim=0), profile="full_combined")
    model.eval()

    sequence_values = torch.randn(2, 5, 31)
    sequence_mask = torch.zeros_like(sequence_values, dtype=torch.bool)
    sequence_embedding = model.supp_encoder(sequence_values, sequence_mask)
    assert torch.count_nonzero(sequence_embedding) == 0

    harmonic_values = torch.randn(2, 7, 7, 31)
    harmonic_mask = torch.zeros_like(harmonic_values, dtype=torch.bool)
    harmonic_embedding = model.harmonic_encoder(harmonic_values, harmonic_mask)
    assert torch.count_nonzero(harmonic_embedding) == 0


def test_unknown_profile_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown harmonic CNN profile"):
        build_harmonic_cnn(HarmonicModelConfig(), profile="not_a_profile")
