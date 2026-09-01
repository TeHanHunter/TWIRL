from __future__ import annotations

from pathlib import Path

import numpy as np

from twirl.models.fm0.cadence_preserving_mechanics import (
    CANDIDATES,
    CONTEXT_LENGTH,
    deterministic_base_samples,
    load_mechanics_config,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/models/twirl_fm0_3_stride1_mechanics_v1.yaml"


def test_stride1_mechanics_config_is_exact_and_training_free() -> None:
    config, source, digest = load_mechanics_config(CONFIG)
    assert source == CONFIG.resolve()
    assert len(digest) == 64
    assert config["cadence_contract"]["patch_stride"] == 1
    assert config["cadence_contract"]["cadence_averaging"] is False
    assert config["authorization"]["optimizer_or_model_training"] is False
    assert set(config["candidates"]) == set(CANDIDATES)


def test_stride1_mechanics_inputs_are_deterministic_and_cadence_preserving() -> None:
    first = deterministic_base_samples(n_sources=3)
    second = deterministic_base_samples(n_sources=3)
    assert len(first) == len(second) == 3
    for left, right in zip(first, second, strict=True):
        assert left["flux"].shape == (2, CONTEXT_LENGTH)
        assert left["time_valid"].shape == (CONTEXT_LENGTH,)
        assert not np.any(left["reconstruction_mask"])
        for key in left:
            assert np.array_equal(left[key], right[key])
