from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np
import pytest

from twirl.models.fm0 import context_window_diagnostic as diagnostic


def _base_sample() -> dict[str, np.ndarray]:
    length = diagnostic.FIXED_EXPOSURE_CADENCES
    flux = np.vstack(
        (
            np.arange(length, dtype=np.float32),
            np.arange(length, dtype=np.float32) + 10_000.0,
        )
    )
    reconstruction = np.zeros_like(flux, dtype=bool)
    reconstruction[:, 11::17] = True
    temporal = np.any(reconstruction, axis=0)
    boundary = np.zeros(length, dtype=bool)
    boundary[0] = True
    return {
        "flux": flux,
        "flux_valid": np.ones_like(flux, dtype=bool),
        "flux_error": np.ones((2, length), dtype=np.float32),
        "error_valid": np.ones((2, length), dtype=bool),
        "local_time_cadences": np.arange(length, dtype=np.float32),
        "delta_time_cadences": np.ones(length, dtype=np.float32),
        "time_valid": np.ones(length, dtype=bool),
        "segment_boundary": boundary,
        "view_present": np.ones(2, dtype=bool),
        "temporal_mask": temporal,
        "reconstruction_mask": reconstruction,
    }


@pytest.mark.parametrize(
    ("context_length", "expected_crops"),
    ((256, 8), (512, 4), (1_024, 2), (2_048, 1)),
)
def test_partition_is_an_exact_fixed_exposure_tiling(
    context_length: int, expected_crops: int
) -> None:
    base = _base_sample()
    before = {name: value.copy() for name, value in base.items()}
    crops = diagnostic.partition_base_sample(base, context_length=context_length)

    assert len(crops) == expected_crops
    assert diagnostic.context_crop_bounds(context_length) == tuple(
        (index * context_length, (index + 1) * context_length)
        for index in range(expected_crops)
    )
    for name in diagnostic._CADENCE_KEYS:
        if name != "local_time_cadences":
            np.testing.assert_array_equal(
                np.concatenate([crop[name] for crop in crops], axis=-1),
                base[name],
            )
    assert all(crop["local_time_cadences"][0] == 0.0 for crop in crops)
    assert sum(
        np.count_nonzero(crop["reconstruction_mask"]) for crop in crops
    ) == np.count_nonzero(base["reconstruction_mask"])
    for name in base:
        np.testing.assert_array_equal(base[name], before[name])


def test_partition_rebases_to_the_first_valid_cadence() -> None:
    base = _base_sample()
    base["time_valid"][256:259] = False
    crop = diagnostic.partition_base_sample(base, context_length=256)[1]
    assert np.all(crop["local_time_cadences"][:3] == 0.0)
    assert crop["local_time_cadences"][3] == 0.0
    assert crop["local_time_cadences"][4] == 1.0


def test_full_coverage_window_is_deterministic_and_rejects_bad_minimum_tile() -> None:
    length = diagnostic.FIXED_EXPOSURE_CADENCES
    release = SimpleNamespace(
        segment_id=np.zeros(length, dtype=np.int32),
        time_valid=np.ones(length, dtype=bool),
        flux_valid=np.ones((length, 6), dtype=bool),
    )
    first = diagnostic._full_coverage_window_spec(
        release, observation_key="observation-a", variant="TWIRL-FM0.2.1"
    )
    second = diagnostic._full_coverage_window_spec(
        release, observation_key="observation-a", variant="TWIRL-FM0.2.1"
    )
    assert first == second
    assert first.n_observed == length
    assert first.n_padded == 0

    release.flux_valid[256:512, 2] = False
    with pytest.raises(ValueError, match="support in every 256-cadence tile"):
        diagnostic._full_coverage_window_spec(
            release, observation_key="observation-a", variant="TWIRL-FM0.2.1"
        )


def test_average_crop_embeddings_rejects_ragged_input() -> None:
    crops = np.arange(24, dtype=np.float64).reshape(6, 4)
    averaged = diagnostic.average_crop_embeddings(crops, n_visits=2, crops_per_visit=3)
    np.testing.assert_allclose(averaged, crops.reshape(2, 3, 4).mean(axis=1))
    with pytest.raises(ValueError, match="rectangular visit grid"):
        diagnostic.average_crop_embeddings(crops[:-1], n_visits=2, crops_per_visit=3)


class _FakeTensor:
    def __init__(self, values: np.ndarray) -> None:
        self.values = values

    def detach(self) -> _FakeTensor:
        return self

    def float(self) -> _FakeTensor:
        return self

    def cpu(self) -> _FakeTensor:
        return self

    def numpy(self) -> np.ndarray:
        return self.values


class _NoGrad:
    def __enter__(self) -> None:
        return None

    def __exit__(self, *_args: object) -> None:
        return None


class _FakeModel:
    def __init__(self) -> None:
        self.full_calls = 0
        self.short_calls = 0

    def eval(self) -> _FakeModel:
        return self

    @staticmethod
    def _output(samples: list[dict[str, np.ndarray]]) -> dict[str, _FakeTensor]:
        values = np.asarray(
            [[sample["marker"], 2.0 * sample["marker"]] for sample in samples],
            dtype=np.float32,
        )
        return {"h_window": _FakeTensor(values), "z_window": _FakeTensor(values + 1)}

    def __call__(self, samples: list[dict[str, np.ndarray]]) -> dict[str, _FakeTensor]:
        self.full_calls += 1
        return self._output(samples)

    def forward_short_context(
        self, samples: list[dict[str, np.ndarray]]
    ) -> dict[str, _FakeTensor]:
        self.short_calls += 1
        return self._output(samples)


def test_encoder_dispatches_short_path_and_averages_to_visits(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_torch = SimpleNamespace(no_grad=lambda: _NoGrad())
    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    monkeypatch.setattr(
        diagnostic, "collate_fm0_samples", lambda samples: list(samples)
    )
    monkeypatch.setattr(
        diagnostic, "move_batch_to_device", lambda samples, _device: samples
    )
    samples = [
        {"marker": np.asarray(value, dtype=np.float32)}
        for value in (1.0, 3.0, 5.0, 7.0)
    ]
    model = _FakeModel()
    short = diagnostic._encode_visit_crops(
        model=model,
        samples=samples,
        n_visits=2,
        crops_per_visit=2,
        context_length=1_024,
        device="cpu",
        batch_size=4,
    )
    np.testing.assert_allclose(short["h_window"], [[2.0, 4.0], [6.0, 12.0]])
    assert model.short_calls == 1
    assert model.full_calls == 0

    diagnostic._encode_visit_crops(
        model=model,
        samples=samples[:2],
        n_visits=2,
        crops_per_visit=1,
        context_length=2_048,
        device="cpu",
        batch_size=2,
    )
    assert model.full_calls == 1


def test_metrics_receive_one_row_per_visit() -> None:
    unmasked = np.asarray(
        [[1.0, 0.0], [1.0, 0.1], [0.0, 1.0], [0.1, 1.0]],
        dtype=np.float64,
    )
    left = unmasked + 0.01
    right = unmasked - 0.01
    result = diagnostic._representation_summary(
        left={"h_window": left, "z_window": left},
        right={"h_window": right, "z_window": right},
        unmasked={"h_window": unmasked, "z_window": unmasked},
        component_ids=("a", "a", "b", "b"),
        sectors=(66, 67, 66, 67),
    )
    health = result["z_window"]["unmasked_visit_embedding_health"]
    assert health["n_embeddings"] == 4
    retrieval = result["z_window"]["query_sector_excluded_cross_visit_retrieval"]
    assert retrieval["status"] == "available"
    assert retrieval["n_visit_embeddings"] == 4


def test_sealed_row_is_rejected_before_any_path_access() -> None:
    with pytest.raises(ValueError, match="only development rows"):
        diagnostic._load_visit_base(
            {"source_partition": "poc_sealed_test"},
            variant="TWIRL-FM0.2.1",
        )
