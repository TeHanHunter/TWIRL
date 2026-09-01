from __future__ import annotations

import hashlib
import inspect
import json
import sys
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from twirl.models.fm0 import centered_event_context_diagnostic as diagnostic

CONFIG = (
    Path(__file__).parents[1]
    / "configs/models/twirl_fm0_2_s66_s77_centered_event_context_v1.yaml"
)


def _base_sample() -> dict[str, np.ndarray]:
    length = diagnostic.BASE_INTERVAL_CADENCES
    boundary = np.zeros(length, dtype=bool)
    boundary[0] = True
    flux = np.vstack(
        (
            np.linspace(-0.01, 0.01, length, dtype=np.float32),
            np.linspace(0.02, -0.02, length, dtype=np.float32),
        )
    )
    return {
        "flux": flux,
        "flux_valid": np.ones_like(flux, dtype=bool),
        "flux_error": np.full((2, length), 0.01, dtype=np.float32),
        "error_valid": np.ones((2, length), dtype=bool),
        "local_time_cadences": np.arange(length, dtype=np.float32),
        "delta_time_cadences": np.ones(length, dtype=np.float32),
        "time_valid": np.ones(length, dtype=bool),
        "segment_boundary": boundary,
        "view_present": np.ones(2, dtype=bool),
        "temporal_mask": np.zeros(length, dtype=bool),
        "reconstruction_mask": np.zeros_like(flux, dtype=bool),
    }


def test_frozen_config_loads_and_binds_exact_bytes(tmp_path: Path) -> None:
    config, source, digest = diagnostic._load_frozen_config(CONFIG)
    assert source == CONFIG.resolve()
    assert digest == hashlib.sha256(CONFIG.read_bytes()).hexdigest()
    assert config["selection"]["max_components_per_cohort"] == 48
    assert config["contexts"]["effective_time_valid_counts_reported"] is True

    changed = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    changed["authorization"]["bls_or_period_features"] = True
    changed_path = tmp_path / "changed.yaml"
    changed_path.write_text(yaml.safe_dump(changed), encoding="utf-8")
    with pytest.raises(ValueError, match="bls_or_period_features"):
        diagnostic._load_frozen_config(changed_path)


def test_inference_loader_never_constructs_optimizer_or_scheduler(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    torch = pytest.importorskip("torch")
    from twirl.models.fm0.model import (
        TWIRLFM0,
        FM0ModelConfig,
        count_trainable_parameters,
    )

    config = FM0ModelConfig(
        architecture="tcn",
        n_flux_views=2,
        window_length=128,
        d_model=8,
        embedding_dim=8,
        stem_kernel=3,
        dropout=0.0,
        tcn_blocks=1,
        tcn_kernel=3,
        tcn_dilation_cycle=(1,),
        conformer_heads=2,
        conformer_conv_kernel=3,
        minimum_parameters=1,
        maximum_parameters=1_000_000,
    )
    model = TWIRLFM0(config)
    contract = {
        "variant": "TWIRL-FM0.2.1",
        "architecture": "tcn",
        "immutable_milestone_steps": [0, 2_000],
    }
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run_contract.json").write_text(json.dumps(contract), encoding="utf-8")
    (run_dir / "summary.json").write_text(
        json.dumps({"parameter_count": count_trainable_parameters(model)}),
        encoding="utf-8",
    )
    checkpoint_path = run_dir / "checkpoint_step_00000000.pt"
    torch.save(
        {
            "schema_version": diagnostic._TRUSTED_CHECKPOINT_SCHEMA_VERSION,
            "model_state": model.state_dict(),
            "model_config": asdict(config),
            "progress": {"global_step": 0},
            "run_contract": contract,
        },
        checkpoint_path,
    )
    checkpoint_hash = diagnostic.sha256_file(checkpoint_path)
    checkpoint_path.with_name(checkpoint_path.name + ".sha256").write_text(
        f"{checkpoint_hash}  {checkpoint_path.name}\n", encoding="utf-8"
    )
    monkeypatch.setattr(
        diagnostic,
        "validate_real_run_release",
        lambda *_args, **kwargs: {
            "artifact_sha256": {"checkpoint.pt": "0" * 64},
            "architecture": "tcn",
            "variant": "TWIRL-FM0.2.1",
            "global_step": 2_000,
            "checkpoint_inspected": kwargs["inspect_checkpoint"],
        },
    )

    def forbidden(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("optimizer or scheduler construction was touched")

    monkeypatch.setattr(torch.optim, "AdamW", forbidden)
    loaded, loaded_contract, validation = diagnostic._load_inference_only_trusted_model(
        run_dir,
        device=torch.device("cpu"),
        checkpoint_path=checkpoint_path,
    )
    assert loaded.training is False
    assert loaded_contract == contract
    assert validation["global_step"] == 0
    assert validation["checkpoint_model_state_inspected"] is True
    assert validation["checkpoint_inspection_mode"] == "model_only_inference"
    assert validation["optimizer_or_scheduler_constructed"] is False


def test_diagnostic_source_has_no_training_module_dependency() -> None:
    source = inspect.getsource(diagnostic)
    assert "from .training" not in source
    assert "make_optimizer_and_scheduler" not in source
    assert "validate_real_run_release(root, inspect_checkpoint=False)" in source


@pytest.mark.parametrize(
    ("length", "bounds"),
    (
        (128, (960, 1_088)),
        (256, (896, 1_152)),
        (512, (768, 1_280)),
        (2_048, (0, 2_048)),
    ),
)
def test_contexts_are_direct_nested_crops_with_one_designated_center(
    length: int, bounds: tuple[int, int]
) -> None:
    assert diagnostic.centered_context_bounds(length) == bounds
    start, stop = diagnostic.centered_event_bounds(length, 9)
    assert (start, stop) == (length // 2 - 4, length // 2 + 5)
    base_start, _ = bounds
    assert base_start + length // 2 == diagnostic.BASE_INTERVAL_CADENCES // 2


def test_cadence_center_trapezoids_are_symmetric_and_frozen() -> None:
    np.testing.assert_allclose(diagnostic.centered_trapezoid(1, 0.1), [0.1])
    np.testing.assert_allclose(diagnostic.centered_trapezoid(3, 0.3), [0.15, 0.3, 0.15])
    np.testing.assert_allclose(
        diagnostic.centered_trapezoid(9, 0.3) / 0.3,
        [1 / 6, 1 / 2, 5 / 6, 1, 1, 1, 5 / 6, 1 / 2, 1 / 6],
    )
    with pytest.raises(ValueError, match="duration"):
        diagnostic.centered_trapezoid(27, 0.1)


def test_injection_changes_only_flux_and_contains_one_centered_event() -> None:
    original = diagnostic.slice_centered_context(_base_sample(), context_length=256)
    before = {key: value.copy() for key, value in original.items()}
    injected, support = diagnostic.inject_centered_event(
        original, duration_cadences=3, fractional_depth=0.1
    )
    start, stop = diagnostic.centered_event_bounds(256, 3)
    assert np.flatnonzero(support).tolist() == list(range(start, stop))
    expected_profile = diagnostic.centered_trapezoid(3, 0.1)
    expected = (1.0 + original["flux"][:, start:stop]) * (
        1.0 - expected_profile[None, :]
    ) - 1.0
    np.testing.assert_allclose(
        injected["flux"][:, start:stop], expected, rtol=1.0e-6, atol=1.0e-8
    )
    np.testing.assert_array_equal(
        injected["flux"][:, :start], original["flux"][:, :start]
    )
    np.testing.assert_array_equal(
        injected["flux"][:, stop:], original["flux"][:, stop:]
    )
    for key in original:
        np.testing.assert_array_equal(original[key], before[key])
        if key != "flux":
            np.testing.assert_array_equal(injected[key], original[key])
    assert not np.any(injected["reconstruction_mask"])


def test_injection_rejects_hidden_contract_fields_and_invalid_support() -> None:
    sample = diagnostic.slice_centered_context(_base_sample(), context_length=128)
    with pytest.raises(ValueError, match="tensor allowlist"):
        diagnostic.inject_centered_event(
            {**sample, "period_days": np.asarray(1.0)},
            duration_cadences=1,
            fractional_depth=0.1,
        )
    sample["flux_valid"][0, 64] = False
    with pytest.raises(diagnostic.CenteredEventIneligibleError):
        diagnostic.inject_centered_event(
            sample, duration_cadences=1, fractional_depth=0.1
        )


def test_centered_slice_rebases_local_time_and_preserves_physical_fields() -> None:
    base = _base_sample()
    start, stop = diagnostic.centered_context_bounds(128)
    base["time_valid"][start : start + 3] = False
    crop = diagnostic.slice_centered_context(base, context_length=128)
    assert crop["local_time_cadences"][:3].tolist() == [0.0, 0.0, 0.0]
    assert crop["local_time_cadences"][3] == 0.0
    assert crop["local_time_cadences"][4] == 1.0
    np.testing.assert_array_equal(
        crop["delta_time_cadences"], base["delta_time_cadences"][start:stop]
    )
    np.testing.assert_array_equal(
        crop["segment_boundary"], base["segment_boundary"][start:stop]
    )


def test_centered_window_requires_fraction_and_full_event_support() -> None:
    length = diagnostic.BASE_INTERVAL_CADENCES
    release = SimpleNamespace(
        segment_id=np.zeros(length, dtype=np.int32),
        time_valid=np.ones(length, dtype=bool),
        flux_valid=np.ones((length, 6), dtype=bool),
    )
    spec, counts = diagnostic._centered_window_spec(
        release,
        observation_key="observation-a",
        variant="TWIRL-FM0.2.1",
    )
    assert spec == diagnostic.WindowSpec(0, 0, length, 0)
    assert counts[128] == (128, 128)

    release.flux_valid[length // 2, 2] = False
    with pytest.raises(diagnostic.CenteredEventIneligibleError, match="no unpadded"):
        diagnostic._centered_window_spec(
            release,
            observation_key="observation-a",
            variant="TWIRL-FM0.2.1",
        )


def _panel_row(
    component: str,
    observation: str,
    cohort: str,
    sector: int,
    *,
    present: str = "[1,1,1,1,1,1]",
) -> dict[str, str]:
    return {
        "leakage_component_id": component,
        "observation_key": observation,
        "temporal_cohort": cohort,
        "source_partition": "poc_development",
        "sector": str(sector),
        "view_present_json": present,
    }


def test_selection_keeps_one_structurally_eligible_visit_per_component(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = [
        _panel_row("repeat-a", "repeat-a-1", "repeated", 66),
        _panel_row("repeat-a", "repeat-a-2", "repeated", 67),
        _panel_row("repeat-b", "repeat-b", "repeated", 68),
        _panel_row("new-a", "new-a", "new", 69),
        _panel_row("new-b", "new-b", "new", 70),
    ]
    repeat_a = [row for row in rows if row["leakage_component_id"] == "repeat-a"]
    first = min(repeat_a, key=diagnostic._visit_order_key)["observation_key"]

    def load(row: dict[str, str], *, variant: str) -> dict[str, object]:
        assert variant == "TWIRL-FM0.2.1"
        if row["observation_key"] == first:
            raise diagnostic.CenteredEventIneligibleError("bad centered support")
        return {
            "sample": _base_sample(),
            "binding": {"observation_key": row["observation_key"]},
        }

    monkeypatch.setattr(diagnostic, "_load_centered_visit", load)
    selected, audit = diagnostic.select_centered_event_visits(
        rows,
        variant="TWIRL-FM0.2.1",
        max_components_per_cohort=2,
    )
    assert len(selected["repeated"]["rows"]) == 2
    assert len(selected["new"]["rows"]) == 2
    assert (
        len({row["leakage_component_id"] for row in selected["repeated"]["rows"]}) == 2
    )
    assert first not in {row["observation_key"] for row in selected["repeated"]["rows"]}
    assert audit["repeated"]["structurally_ineligible_visits_skipped"] == 1
    assert audit["repeated"]["screening_shards_opened"] == 3
    assert audit["new"]["screening_shards_opened"] == 2
    access = diagnostic._development_shard_access_summary(audit)
    assert access == {
        "screening_development_shards_opened": 5,
        "selected_development_shards": 4,
        "screening_shards_opened_but_not_selected": 1,
        "each_screened_development_shard_opened_at_most_once": True,
    }


def test_selection_skips_a_missing_view_visit_without_opening_it(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = [
        _panel_row(
            "repeat-a",
            "repeat-a-missing",
            "repeated",
            66,
            present="[1,1,0,1,1,1]",
        ),
        _panel_row("repeat-a", "repeat-a-good", "repeated", 67),
        _panel_row("repeat-b", "repeat-b", "repeated", 68),
        _panel_row("new-a", "new-a", "new", 69),
        _panel_row("new-b", "new-b", "new", 70),
    ]
    opened: list[str] = []

    def load(row: dict[str, str], *, variant: str) -> dict[str, object]:
        opened.append(row["observation_key"])
        return {
            "sample": _base_sample(),
            "binding": {"observation_key": row["observation_key"]},
        }

    monkeypatch.setattr(diagnostic, "_load_centered_visit", load)
    selected, audit = diagnostic.select_centered_event_visits(
        rows,
        variant="TWIRL-FM0.2.1",
        max_components_per_cohort=2,
    )
    assert "repeat-a-missing" not in opened
    assert "repeat-a-good" in {
        row["observation_key"] for row in selected["repeated"]["rows"]
    }
    assert audit["repeated"]["missing_required_view_visits_skipped"] == 1


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
        self.short_calls = 0
        self.full_calls = 0

    def eval(self) -> _FakeModel:
        return self

    @staticmethod
    def _output(samples: list[dict[str, np.ndarray]]) -> dict[str, _FakeTensor]:
        markers = np.asarray([sample["marker"] for sample in samples], dtype=np.float32)
        length = samples[0]["flux"].shape[-1]
        reconstruction = np.stack([sample["flux"] for sample in samples]) * 0.5
        return {
            "h_window": _FakeTensor(np.column_stack((markers, 2 * markers))),
            "z_window": _FakeTensor(np.column_stack((markers + 1, 2 * markers + 1))),
            "reconstruction": _FakeTensor(
                reconstruction.reshape(len(samples), 2, length)
            ),
        }

    def __call__(self, samples: list[dict[str, np.ndarray]]) -> dict[str, _FakeTensor]:
        self.full_calls += 1
        return self._output(samples)

    def forward_short_context(
        self, samples: list[dict[str, np.ndarray]]
    ) -> dict[str, _FakeTensor]:
        self.short_calls += 1
        return self._output(samples)


def test_encoder_uses_direct_short_path_and_keeps_reconstruction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(
        sys.modules, "torch", SimpleNamespace(no_grad=lambda: _NoGrad())
    )
    monkeypatch.setattr(
        diagnostic, "collate_fm0_samples", lambda samples: list(samples)
    )
    monkeypatch.setattr(
        diagnostic, "move_batch_to_device", lambda samples, _device: samples
    )
    samples = [
        {
            "flux": np.full((2, 128), value, dtype=np.float32),
            "marker": np.asarray(value, dtype=np.float32),
        }
        for value in (1.0, 2.0)
    ]
    model = _FakeModel()
    encoded = diagnostic._encode_samples(
        model=model,
        samples=samples,
        context_length=128,
        device="cpu",
        batch_size=2,
    )
    assert encoded["h_window"].shape == (2, 2)
    assert encoded["reconstruction"].shape == (2, 2, 128)
    assert model.short_calls == 1
    assert model.full_calls == 0

    full = [{"flux": np.ones((2, 2_048), dtype=np.float32), "marker": np.asarray(1.0)}]
    diagnostic._encode_samples(
        model=model,
        samples=full,
        context_length=2_048,
        device="cpu",
        batch_size=1,
    )
    assert model.full_calls == 1


def test_embedding_response_reports_scaled_displacement_and_coherence() -> None:
    original = np.asarray(
        [[1.0, 1.0], [2.0, 3.0], [3.0, 2.0], [4.0, 4.0]], dtype=np.float64
    )
    delta = np.asarray([0.5, -0.25], dtype=np.float64)
    injected = original + delta
    scale, audit = diagnostic._fit_original_robust_scale(original)
    result = diagnostic.summarize_embedding_response(
        original,
        injected,
        robust_scale=scale,
        component_ids=("a", "b", "c", "d"),
        condition_ids=("event",) * 4,
    )
    assert audit["positive_scale_dimensions"] == 2
    assert result["raw_l2_displacement"]["overall"]["mean"] == pytest.approx(
        np.linalg.norm(delta)
    )
    coherence = result["coherent_shift_and_direction"]["overall"]
    assert coherence["unit_delta_direction_coherence"] == pytest.approx(1.0)
    assert coherence["zero_displacement_pairs"] == 0


def test_visible_input_decoder_response_has_known_gain_and_residual() -> None:
    original = np.ones((2, 2, 9), dtype=np.float64)
    injected = original.copy()
    injected[:, :, 4] -= 0.2
    support = np.zeros((2, 9), dtype=bool)
    support[:, 4] = True
    valid = np.ones_like(original, dtype=bool)
    truth = injected - original

    perfect = diagnostic.visible_event_decoder_response(
        original,
        injected,
        np.zeros_like(original),
        truth,
        flux_valid=valid,
        event_support=support,
    )
    np.testing.assert_allclose(perfect["visible_event_projection_gain"], 1.0)
    np.testing.assert_allclose(perfect["event_support_normalized_residual"], 0.0)

    zero = diagnostic.visible_event_decoder_response(
        original,
        injected,
        np.zeros_like(original),
        np.zeros_like(original),
        flux_valid=valid,
        event_support=support,
    )
    np.testing.assert_allclose(zero["visible_event_projection_gain"], 0.0)
    np.testing.assert_allclose(zero["event_support_normalized_residual"], 1.0)


def test_condition_grid_applies_all_twelve_single_event_conditions() -> None:
    conditions = diagnostic._condition_grid()
    assert len(conditions) == 12
    assert {
        (item["duration_cadences"], item["fractional_depth"]) for item in conditions
    } == {
        (duration, depth)
        for duration in diagnostic.EVENT_DURATIONS_CADENCES
        for depth in diagnostic.EVENT_FRACTIONAL_DEPTHS
    }


def test_evaluator_exposes_the_thin_cli_signature() -> None:
    signature = inspect.signature(diagnostic.evaluate_centered_event_context)
    assert tuple(signature.parameters) == (
        "config_path",
        "run_dir",
        "step0_checkpoint_path",
        "step2000_checkpoint_path",
        "temporal_panel_dir",
        "temporal_panel_receipt_sha256",
        "batch_size",
    )
    assert all(
        parameter.kind is inspect.Parameter.KEYWORD_ONLY
        for parameter in signature.parameters.values()
    )


def test_cpu_wrapper_pins_config_and_forbids_gpu_training() -> None:
    root = Path(__file__).parents[1]
    wrapper = (
        root / "scripts/orcd/slurm_twirl_fm0_s66_s77_centered_event_context_cpu.sbatch"
    ).read_text(encoding="utf-8")
    config_hash = hashlib.sha256(CONFIG.read_bytes()).hexdigest()

    assert f'CONFIG_SHA256="{config_hash}"' in wrapper
    assert "#SBATCH --gres" not in wrapper
    assert "evaluate_twirl_fm0_centered_event_context.py" in wrapper
    assert "TWIRL_FM0_CENTERED_EVENT_OUTPUT" in wrapper
    assert "no BLS, period, labels, probe, training, or sealed access" in wrapper
