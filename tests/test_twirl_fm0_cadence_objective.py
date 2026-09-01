from __future__ import annotations

import copy
import math
from dataclasses import asdict
from pathlib import Path

import pytest

torch = pytest.importorskip("torch")

import twirl.models.fm0.training as fm0_training
import twirl.models.fm0.validation as fm0_validation
from twirl.models.fm0.cadence_objective import (
    position_centered_cadence_vicreg_loss,
)
from twirl.models.fm0.dataset import (
    SyntheticFM0Config,
    SyntheticFM0Dataset,
    collate_fm0_samples,
)
from twirl.models.fm0.model import (
    FM0ModelConfig,
    build_fm0_model,
)
from twirl.models.fm0.training import (
    CADENCE_VICREG_OBJECTIVE_IDENTITY,
    OBJECTIVE_STATE_SCHEMA_V2,
    FM0OptimizationConfig,
    run_synthetic_training,
    seed_everything,
)


def _tiny_fm03_model(*, dropout: float = 0.1):
    return build_fm0_model(
        "TWIRL-FM0.3.1",
        config_override=FM0ModelConfig(
            architecture="tcn",
            n_flux_views=2,
            window_length=128,
            d_model=8,
            embedding_dim=8,
            stem_kernel=5,
            dropout=dropout,
            tcn_blocks=1,
            tcn_dilation_cycle=(1,),
            conformer_heads=2,
            conformer_conv_kernel=7,
            patch_stride=1,
            minimum_parameters=1,
            maximum_parameters=1_000_000,
        ),
        enforce_parameter_budget=False,
    )


def _fm03_dataset_config(*, windows_per_epoch: int = 8) -> SyntheticFM0Config:
    return SyntheticFM0Config(
        variant="TWIRL-FM0.3.1",
        seed=560067,
        windows_per_epoch=windows_per_epoch,
        window_length=128,
        mask_target_fraction=0.15,
        mask_span_range=(1, 4),
    )


def _fm03_contract(
    *,
    optimizer_horizon: int = 1,
    milestones: list[int] | None = None,
) -> dict[str, object]:
    contract: dict[str, object] = {
        "schema_version": "test",
        "variant": "TWIRL-FM0.3.1",
        "objective": CADENCE_VICREG_OBJECTIVE_IDENTITY,
        "use_vicreg": True,
        "reconstruct_second_view": False,
        "mask_target_fraction": 0.15,
        "mask_span_range": [1, 4],
        "training_horizon_step": optimizer_horizon,
        "stop_after_step_is_execution_state_not_scientific_contract": True,
    }
    if milestones is not None:
        contract["immutable_milestone_steps"] = milestones
    return contract


def test_position_centering_does_not_treat_time_offsets_as_batch_variance() -> None:
    # Each fixed cadence position is constant across the batch, even though
    # values differ strongly between positions.  Pooling batch and time would
    # incorrectly treat this as a high-variance representation.
    cadence_values = torch.tensor([-10.0, 10.0]).view(1, 2, 1)
    first = cadence_values.expand(4, -1, -1).clone().requires_grad_()
    second = first.detach().clone().requires_grad_()
    visible = torch.ones((4, 2), dtype=torch.bool)

    loss, diagnostics = position_centered_cadence_vicreg_loss(
        first,
        second,
        visible,
        visible,
        invariance_weight=0.0,
        covariance_weight=0.0,
    )

    torch.testing.assert_close(diagnostics["variance"], torch.tensor(0.99))
    torch.testing.assert_close(loss, torch.tensor(24.75))
    assert diagnostics["invariance"] == 0
    assert diagnostics["statistical_cadence_positions"] == 2
    assert diagnostics["statistical_degrees_of_freedom"] == 6


def test_common_visible_mask_and_fp32_gradient_are_stable() -> None:
    first = torch.tensor(
        [
            [[0.0, 1.0], [float("nan"), float("nan")]],
            [[2.0, 3.0], [4.0, 5.0]],
            [[4.0, 5.0], [6.0, 7.0]],
        ],
        dtype=torch.bfloat16,
        requires_grad=True,
    )
    second = torch.tensor(
        [
            [[0.5, 1.5], [2.0, 3.0]],
            [[2.5, 3.5], [4.5, 5.5]],
            [[4.5, 5.5], [6.5, 7.5]],
        ],
        dtype=torch.bfloat16,
        requires_grad=True,
    )
    first_visible = torch.tensor(
        [[True, False], [True, True], [True, True]], dtype=torch.bool
    )
    second_visible = torch.tensor(
        [[True, True], [True, True], [True, True]], dtype=torch.bool
    )

    with torch.autocast(device_type="cpu", dtype=torch.bfloat16):
        loss, diagnostics = position_centered_cadence_vicreg_loss(
            first, second, first_visible, second_visible
        )
    loss.backward()

    assert loss.dtype == torch.float32
    assert torch.isfinite(loss)
    assert torch.isfinite(first.grad[first_visible]).all()
    assert torch.isfinite(second.grad[second_visible]).all()
    assert diagnostics["position_batch_counts"].tolist() == [3, 2]
    assert diagnostics["common_visible_tokens"] == 5
    assert diagnostics["statistical_degrees_of_freedom"] == 3


def test_pooled_statistics_use_global_position_centered_df() -> None:
    first = torch.tensor(
        [
            [[0.00, 0.00], [0.00, 0.50]],
            [[0.25, 0.50], [0.50, 0.00]],
            [[0.50, 1.00], [9.00, 9.00]],
        ]
    )
    second = first + 0.1
    first_visible = torch.tensor(
        [[True, True], [True, True], [True, False]], dtype=torch.bool
    )
    second_visible = torch.ones((3, 2), dtype=torch.bool)

    loss, diagnostics = position_centered_cadence_vicreg_loss(
        first, second, first_visible, second_visible
    )

    expected_invariance = 0.01
    expected_variance = 0.5 * (
        1.0 - math.sqrt(0.25 / 3.0 + 1.0e-4) + 1.0 - math.sqrt(0.625 / 3.0 + 1.0e-4)
    )
    expected_covariance = 2.0 * (0.125 / 3.0) ** 2
    expected_loss = (
        25.0 * expected_invariance + 25.0 * expected_variance + expected_covariance
    )
    torch.testing.assert_close(
        diagnostics["invariance"], torch.tensor(expected_invariance)
    )
    torch.testing.assert_close(diagnostics["variance"], torch.tensor(expected_variance))
    torch.testing.assert_close(
        diagnostics["covariance"], torch.tensor(expected_covariance)
    )
    torch.testing.assert_close(loss, torch.tensor(expected_loss))
    assert diagnostics["statistical_degrees_of_freedom"] == 3


@pytest.mark.parametrize(
    ("first_visible", "second_visible", "message"),
    [
        (
            [[True, False], [False, False]],
            [[True, False], [False, False]],
            "insufficient degrees of freedom",
        ),
        (
            [[False, False], [False, False]],
            [[False, False], [False, False]],
            "no common-visible tokens",
        ),
    ],
)
def test_cadence_vicreg_fails_closed_without_position_df(
    first_visible: list[list[bool]],
    second_visible: list[list[bool]],
    message: str,
) -> None:
    first = torch.zeros((2, 2, 3))
    second = torch.zeros_like(first)

    with pytest.raises(ValueError, match=message):
        position_centered_cadence_vicreg_loss(
            first,
            second,
            torch.tensor(first_visible, dtype=torch.bool),
            torch.tensor(second_visible, dtype=torch.bool),
        )


def test_singleton_position_is_ignored_by_pooled_statistics() -> None:
    first = torch.tensor(
        [
            [[100.0, -100.0], [0.0, 1.0]],
            [[0.0, 0.0], [2.0, 3.0]],
            [[0.0, 0.0], [4.0, 5.0]],
            [[0.0, 0.0], [6.0, 7.0]],
        ]
    )
    second = first.clone()
    visible = torch.tensor(
        [[True, True], [False, True], [False, True], [False, True]],
        dtype=torch.bool,
    )

    loss, diagnostics = position_centered_cadence_vicreg_loss(
        first, second, visible, visible
    )

    assert torch.isfinite(loss)
    assert diagnostics["position_batch_counts"].tolist() == [1, 4]
    assert diagnostics["statistical_cadence_positions"] == 1
    assert diagnostics["statistical_degrees_of_freedom"] == 3


def test_cadence_vicreg_rejects_nonfinite_common_visible_tokens() -> None:
    first = torch.zeros((3, 2, 3))
    first[0, 0, 0] = float("nan")
    second = torch.zeros_like(first)
    visible = torch.ones((3, 2), dtype=torch.bool)

    with pytest.raises(ValueError, match="nonfinite common-visible token"):
        position_centered_cadence_vicreg_loss(first, second, visible, visible)


def test_paired_visibility_excludes_either_mask_and_rejects_source_drift() -> None:
    flux = torch.zeros((1, 2, 5))
    flux_valid = torch.tensor(
        [[[True, True, True, True, False], [True, True, True, True, True]]]
    )
    common = {
        "flux": flux,
        "flux_valid": flux_valid,
        "flux_error": torch.ones((1, 2, 5)),
        "error_valid": torch.ones((1, 2, 5), dtype=torch.bool),
        "local_time_cadences": torch.arange(5).view(1, 5).float(),
        "delta_time_cadences": torch.ones((1, 5)),
        "time_valid": torch.tensor([[True, True, True, False, True]]),
        "segment_boundary": torch.tensor([[True, False, False, False, False]]),
        "view_present": torch.tensor([[True, False]]),
    }
    first = {
        **common,
        "reconstruction_mask": torch.tensor(
            [[[False] * 5, [True, False, False, False, False]]]
        ),
    }
    second = {
        **common,
        "reconstruction_mask": torch.tensor(
            [[[False, True, False, False, False], [False] * 5]]
        ),
    }

    first_visible, second_visible = fm0_training._paired_cadence_visibility_masks(
        first, second
    )
    assert first_visible.tolist() == [[False, True, True, False, False]]
    assert second_visible.tolist() == [[True, False, True, False, False]]
    assert (first_visible & second_visible).tolist() == [
        [False, False, True, False, False]
    ]

    for name in (
        "flux",
        "flux_valid",
        "flux_error",
        "error_valid",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
        "view_present",
    ):
        drifted_second = dict(second)
        drifted_second[name] = second[name].clone()
        flattened = drifted_second[name].reshape(-1)
        if flattened.dtype == torch.bool:
            flattened[0] = ~flattened[0]
        else:
            flattened[0] += 1
        with pytest.raises(ValueError, match=f"paired source tensor differs: {name}"):
            fm0_training._paired_cadence_visibility_masks(first, drifted_second)


def test_legacy_boolean_objective_state_never_infers_fm03_identity() -> None:
    legacy = fm0_training._validated_objective_state(
        {"use_vicreg": True, "reconstruct_second_view": False},
        label="checkpoint objective state",
    )
    cadence = fm0_training._objective_state(
        use_vicreg=True,
        reconstruct_second_view=False,
        objective_identity=CADENCE_VICREG_OBJECTIVE_IDENTITY,
    )
    assert legacy == {"use_vicreg": True, "reconstruct_second_view": False}
    assert cadence["identity"] == CADENCE_VICREG_OBJECTIVE_IDENTITY
    assert legacy != cadence


def _cadence_validation_skeleton(
    model_config: dict[str, object],
) -> tuple[dict, dict, dict]:
    contract = {
        "variant": "TWIRL-FM0.3.1",
        "architecture": "tcn",
        "objective": CADENCE_VICREG_OBJECTIVE_IDENTITY,
        "use_vicreg": True,
        "reconstruct_second_view": False,
        "mask_target_fraction": 0.15,
        "mask_span_range": [1, 4],
        "training_horizon_step": 1,
        "stop_after_step_is_execution_state_not_scientific_contract": True,
        "optimization": {"max_optimizer_steps": 1},
    }
    summary = {
        "global_step": 1,
        "requested_stop_after_step": 1,
        "architecture": "tcn",
        "variant": "TWIRL-FM0.3.1",
        "parameter_count": 0,
    }
    checkpoint = {
        "schema_version": fm0_training.CHECKPOINT_SCHEMA_VERSION,
        "model_state": {},
        "optimizer_state": {},
        "scheduler_state": {},
        "scaler_state": None,
        "rng_state": {},
        "progress": {"global_step": 1},
        "sampler_state": {},
        "model_config": model_config,
        "optimization_config": {},
        "synthetic_dataset_config": {},
        "dataset_contract": {},
        "run_contract": contract,
        "loss_history": [],
        "objective_state": {
            "schema_version": OBJECTIVE_STATE_SCHEMA_V2,
            "identity": CADENCE_VICREG_OBJECTIVE_IDENTITY,
            "use_vicreg": True,
            "reconstruct_second_view": False,
        },
    }
    return checkpoint, contract, summary


def test_release_validation_rejects_missing_fm03_objective_state(
    tmp_path: Path,
) -> None:
    checkpoint, contract, summary = _cadence_validation_skeleton(
        asdict(_tiny_fm03_model().config)
    )
    del checkpoint["objective_state"]
    path = tmp_path / "checkpoint.pt"
    torch.save(checkpoint, path)
    with pytest.raises(ValueError, match="cadence objective checkpoint contract"):
        fm0_validation._inspect_trusted_checkpoint(
            path,
            contract=contract,
            summary=summary,
        )


def test_release_validation_rejects_noncanonical_fm03_model_config(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from twirl.models.fm0 import model as fm0_model

    canonical = _tiny_fm03_model(dropout=0.1)
    drifted = _tiny_fm03_model(dropout=0.2)
    monkeypatch.setattr(
        fm0_model, "build_fm0_model", lambda *_args, **_kwargs: canonical
    )
    checkpoint, contract, summary = _cadence_validation_skeleton(asdict(drifted.config))
    path = tmp_path / "checkpoint.pt"
    torch.save(checkpoint, path)
    with pytest.raises(ValueError, match="canonical native-cadence model"):
        fm0_validation._inspect_trusted_checkpoint(
            path,
            contract=contract,
            summary=summary,
        )


def test_cadence_replay_matches_retained_graph_and_rng() -> None:
    dataset = SyntheticFM0Dataset(_fm03_dataset_config(windows_per_epoch=4))
    paired_batches = []
    for start in (0, 2):
        paired_batches.append(
            (
                collate_fm0_samples(
                    [
                        dataset.sample(index, mask_view=0)
                        for index in range(start, start + 2)
                    ]
                ),
                collate_fm0_samples(
                    [
                        dataset.sample(index, mask_view=1)
                        for index in range(start, start + 2)
                    ]
                ),
            )
        )
    optimization = FM0OptimizationConfig(
        effective_batch_windows=4,
        max_optimizer_steps=1,
        warmup_steps=0,
        vicreg_total_weight=2.0e-4,
    )
    seed_everything(12345)
    template = _tiny_fm03_model(dropout=0.1)
    retained_model = copy.deepcopy(template).train()
    replay_model = copy.deepcopy(template).train()

    torch.manual_seed(991)
    retained_first = []
    retained_second = []
    retained_first_visible = []
    retained_second_visible = []
    for first_batch, second_batch in paired_batches:
        first_output = retained_model(first_batch)
        retained_first.append(
            fm0_training._project_cadence_tokens(retained_model, first_output)
        )
        second_output = retained_model(second_batch)
        retained_second.append(
            fm0_training._project_cadence_tokens(retained_model, second_output)
        )
        first_visible, second_visible = fm0_training._paired_cadence_visibility_masks(
            first_batch, second_batch
        )
        retained_first_visible.append(first_visible)
        retained_second_visible.append(second_visible)
    retained_loss, retained_diagnostics = position_centered_cadence_vicreg_loss(
        torch.cat(retained_first),
        torch.cat(retained_second),
        torch.cat(retained_first_visible),
        torch.cat(retained_second_visible),
        invariance_weight=optimization.vicreg_invariance_weight,
        variance_weight=optimization.vicreg_variance_weight,
        covariance_weight=optimization.vicreg_covariance_weight,
    )
    (optimization.vicreg_total_weight * retained_loss).backward()

    torch.manual_seed(991)
    replay_records = []
    replay_first = []
    replay_second = []
    replay_first_visible = []
    replay_second_visible = []
    for first_batch, second_batch in paired_batches:
        first_rng = fm0_training._capture_forward_rng_state(torch.device("cpu"))
        with torch.no_grad():
            first_output = replay_model(first_batch)
            first_z = fm0_training._project_cadence_tokens(replay_model, first_output)
        second_rng = fm0_training._capture_forward_rng_state(torch.device("cpu"))
        with torch.no_grad():
            second_output = replay_model(second_batch)
            second_z = fm0_training._project_cadence_tokens(replay_model, second_output)
        first_visible, second_visible = fm0_training._paired_cadence_visibility_masks(
            first_batch, second_batch
        )
        replay_first.append(first_z)
        replay_second.append(second_z)
        replay_first_visible.append(first_visible)
        replay_second_visible.append(second_visible)
        replay_records.append(
            {
                "first_batch": first_batch,
                "second_batch": second_batch,
                "first_rng": first_rng,
                "second_rng": second_rng,
                "batch_windows": 2,
            }
        )
    post_forward_rng = fm0_training._capture_forward_rng_state(torch.device("cpu"))
    replay_loss, replay_diagnostics = (
        fm0_training._backprop_effective_batch_cadence_vicreg(
            model=replay_model,
            replay_records=replay_records,
            first_embeddings=replay_first,
            second_embeddings=replay_second,
            first_visible_masks=replay_first_visible,
            second_visible_masks=replay_second_visible,
            config=optimization,
            device=torch.device("cpu"),
            precision="fp32",
        )
    )

    assert torch.equal(torch.get_rng_state(), post_forward_rng["cpu"])
    torch.testing.assert_close(replay_loss, retained_loss.detach())
    for name in ("invariance", "variance", "covariance"):
        torch.testing.assert_close(
            replay_diagnostics[name], retained_diagnostics[name].detach()
        )
    for (retained_name, retained_parameter), (
        replay_name,
        replay_parameter,
    ) in zip(retained_model.named_parameters(), replay_model.named_parameters()):
        assert retained_name == replay_name
        if retained_parameter.grad is None or replay_parameter.grad is None:
            assert retained_parameter.grad is None
            assert replay_parameter.grad is None
            continue
        torch.testing.assert_close(
            replay_parameter.grad,
            retained_parameter.grad,
            rtol=1.0e-5,
            atol=1.0e-7,
            msg=retained_name,
        )


def test_fm03_training_uses_full_cadence_objective_and_first_reconstruction(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    original_cadence_loss = fm0_training.position_centered_cadence_vicreg_loss

    def recording_cadence_loss(first, second, first_visible, second_visible, **kwargs):
        observed.append((tuple(first.shape), tuple(first_visible.shape)))
        return original_cadence_loss(
            first, second, first_visible, second_visible, **kwargs
        )

    def reject_legacy_vicreg(*_args, **_kwargs):
        raise AssertionError("FM0.3 must not optimize pooled z_window VICReg")

    monkeypatch.setattr(
        fm0_training,
        "position_centered_cadence_vicreg_loss",
        recording_cadence_loss,
    )
    monkeypatch.setattr(fm0_training, "vicreg_loss", reject_legacy_vicreg)
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=0,
        max_optimizer_steps=1,
        effective_batch_windows=4,
        vicreg_total_weight=2.0e-4,
    )
    seed_everything(560067)
    model = _tiny_fm03_model()
    ordinary_forward = model.forward

    def poisoned_window_diagnostics(batch):
        output = ordinary_forward(batch)
        output["h_window"] = torch.full_like(output["h_window"], float("nan"))
        output["z_window"] = torch.full_like(output["z_window"], float("nan"))
        return output

    monkeypatch.setattr(model, "forward", poisoned_window_diagnostics)
    result = run_synthetic_training(
        model=model,
        dataset=SyntheticFM0Dataset(_fm03_dataset_config(windows_per_epoch=4)),
        output_dir=tmp_path,
        run_contract=_fm03_contract(),
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=True,
        objective_identity=CADENCE_VICREG_OBJECTIVE_IDENTITY,
    )

    assert observed == [((4, 128, 8), (4, 128))]
    metrics = result["final_metrics"]
    assert metrics["reconstruction"] == metrics["reconstruction_first"]
    assert metrics["embedding_projection_gradient_norm"] > 0.0
    for name in (
        "vicreg_invariance",
        "vicreg_variance",
        "vicreg_covariance",
        "vicreg_common_visible_tokens",
        "vicreg_statistical_cadence_positions",
        "vicreg_statistical_degrees_of_freedom",
    ):
        assert math.isfinite(metrics[name])
    checkpoint = torch.load(tmp_path / "checkpoint.pt", weights_only=False)
    assert checkpoint["objective_state"] == {
        "schema_version": OBJECTIVE_STATE_SCHEMA_V2,
        "identity": CADENCE_VICREG_OBJECTIVE_IDENTITY,
        "use_vicreg": True,
        "reconstruct_second_view": False,
    }


@pytest.mark.parametrize(
    ("contract_update", "dataset_update", "message"),
    [
        ({"mask_target_fraction": 0.2}, {}, "mask target fraction"),
        ({"mask_span_range": [1, 5]}, {}, "mask span range"),
        ({}, {"mask_span_range": (1, 5)}, "dataset mask geometry"),
    ],
)
def test_fm03_training_rejects_mask_contract_drift(
    tmp_path: Path,
    contract_update: dict[str, object],
    dataset_update: dict[str, object],
    message: str,
) -> None:
    contract = _fm03_contract()
    contract.update(contract_update)
    config_values = {
        "variant": "TWIRL-FM0.3.1",
        "seed": 560067,
        "windows_per_epoch": 4,
        "window_length": 128,
        "mask_target_fraction": 0.15,
        "mask_span_range": (1, 4),
        **dataset_update,
    }
    with pytest.raises(ValueError, match=message):
        run_synthetic_training(
            model=_tiny_fm03_model(),
            dataset=SyntheticFM0Dataset(SyntheticFM0Config(**config_values)),
            output_dir=tmp_path,
            run_contract=contract,
            optimization=FM0OptimizationConfig(
                max_optimizer_steps=1, effective_batch_windows=4
            ),
            target_step=1,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=True,
            objective_identity=CADENCE_VICREG_OBJECTIVE_IDENTITY,
        )


@pytest.mark.parametrize(
    ("contract_update", "message"),
    [
        ({"target_step": 1}, "invocation stop must remain runtime-only"),
        ({"training_horizon_step": 2}, "optimizer horizon differs"),
        (
            {"stop_after_step_is_execution_state_not_scientific_contract": False},
            "runtime-only stop declaration is missing",
        ),
    ],
)
def test_fm03_training_rejects_invocation_stop_in_invariant_contract(
    tmp_path: Path,
    contract_update: dict[str, object],
    message: str,
) -> None:
    contract = _fm03_contract()
    contract.update(contract_update)
    with pytest.raises(ValueError, match=message):
        run_synthetic_training(
            model=_tiny_fm03_model(),
            dataset=SyntheticFM0Dataset(_fm03_dataset_config(windows_per_epoch=4)),
            output_dir=tmp_path,
            run_contract=contract,
            optimization=FM0OptimizationConfig(
                max_optimizer_steps=1, effective_batch_windows=4
            ),
            target_step=1,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=True,
            objective_identity=CADENCE_VICREG_OBJECTIVE_IDENTITY,
        )


def test_fm03_checkpoint_resume_matches_uninterrupted_training(
    tmp_path: Path,
) -> None:
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=1,
        max_optimizer_steps=2,
        effective_batch_windows=2,
        vicreg_total_weight=2.0e-4,
    )
    contract = _fm03_contract(optimizer_horizon=2)
    dataset_config = _fm03_dataset_config(windows_per_epoch=8)

    full_dir = tmp_path / "full"
    full_dir.mkdir()
    seed_everything(560067)
    run_synthetic_training(
        model=_tiny_fm03_model(),
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=full_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=2,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=True,
        objective_identity=CADENCE_VICREG_OBJECTIVE_IDENTITY,
    )

    resumed_dir = tmp_path / "resumed"
    resumed_dir.mkdir()
    seed_everything(560067)
    run_synthetic_training(
        model=_tiny_fm03_model(),
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=resumed_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=True,
        objective_identity=CADENCE_VICREG_OBJECTIVE_IDENTITY,
    )
    seed_everything(999)
    run_synthetic_training(
        model=_tiny_fm03_model(),
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=resumed_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=2,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=True,
        objective_identity=CADENCE_VICREG_OBJECTIVE_IDENTITY,
        resume_checkpoint=resumed_dir / "checkpoint.pt",
    )

    full = torch.load(full_dir / "checkpoint.pt", weights_only=False)
    resumed = torch.load(resumed_dir / "checkpoint.pt", weights_only=False)
    assert full["loss_history"] == resumed["loss_history"]
    assert full["progress"] == resumed["progress"]
    assert full["sampler_state"] == resumed["sampler_state"]
    assert full["objective_state"] == resumed["objective_state"]
    assert torch.equal(
        full["rng_state"]["torch_cpu"], resumed["rng_state"]["torch_cpu"]
    )
    assert set(full["model_state"]) == set(resumed["model_state"])
    for name in full["model_state"]:
        torch.testing.assert_close(
            full["model_state"][name], resumed["model_state"][name]
        )
