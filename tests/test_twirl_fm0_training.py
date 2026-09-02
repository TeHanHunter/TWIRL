from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

torch = pytest.importorskip("torch")

import twirl.models.fm0.training as fm0_training  # noqa: E402
from twirl.models.fm0.dataset import (  # noqa: E402
    SyntheticFM0Config,
    SyntheticFM0Dataset,
    collate_fm0_samples,
)
from twirl.models.fm0.model import FM0ModelConfig, build_fm0_model  # noqa: E402
from twirl.models.fm0.training import (  # noqa: E402
    FM0OptimizationConfig,
    fm0_objective,
    masked_huber_reconstruction_loss,
    run_synthetic_training,
    seed_everything,
    vicreg_loss,
)
from twirl.models.fm0.validation import (  # noqa: E402
    RUN_CONTRACT_SCHEMA_VERSION,
    RUN_SUMMARY_SCHEMA_VERSION,
    read_json,
    validate_run_release,
    write_json_with_sha256,
    write_sha256_sidecar,
)


def _tiny_model():
    config = FM0ModelConfig(
        architecture="tcn",
        n_flux_views=2,
        window_length=32,
        d_model=8,
        embedding_dim=8,
        stem_kernel=5,
        dropout=0.1,
        tcn_blocks=1,
        tcn_dilation_cycle=(1,),
        conformer_heads=2,
        conformer_conv_kernel=7,
        minimum_parameters=1,
        maximum_parameters=1_000_000,
    )
    return build_fm0_model(
        "TWIRL-FM0.1.1", config_override=config, enforce_parameter_budget=False
    )


def _tiny_vicreg_model():
    config = FM0ModelConfig(
        architecture="tcn",
        n_flux_views=6,
        window_length=32,
        d_model=8,
        embedding_dim=8,
        stem_kernel=5,
        dropout=0.1,
        tcn_blocks=1,
        tcn_dilation_cycle=(1,),
        conformer_heads=2,
        conformer_conv_kernel=7,
        minimum_parameters=1,
        maximum_parameters=1_000_000,
    )
    return build_fm0_model(
        "TWIRL-FM0.1.5",
        development_winner="tcn",
        config_override=config,
        enforce_parameter_budget=False,
    )


def _write_structural_release(run_dir: Path) -> None:
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=1,
        max_optimizer_steps=1,
        effective_batch_windows=2,
    )
    contract = {
        "schema_version": RUN_CONTRACT_SCHEMA_VERSION,
        "variant": "TWIRL-FM0.1.1",
        "architecture": "tcn",
        "seed": 560067,
        "synthetic_smoke": True,
        "target_step": 1,
        "optimization": asdict(optimization),
    }
    run_dir.mkdir()
    contract_hash = write_json_with_sha256(
        run_dir / "run_contract.json", contract
    )
    seed_everything(560067)
    model = _tiny_model()
    result = run_synthetic_training(
        model=model,
        dataset=SyntheticFM0Dataset(
            SyntheticFM0Config(
                variant="TWIRL-FM0.1.1",
                seed=560067,
                windows_per_epoch=8,
                window_length=32,
            )
        ),
        output_dir=run_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=False,
    )
    checkpoint_hash = write_sha256_sidecar(run_dir / "checkpoint.pt")
    write_json_with_sha256(
        run_dir / "summary.json",
        {
            "schema_version": RUN_SUMMARY_SCHEMA_VERSION,
            "passed": True,
            "synthetic_only": True,
            "global_step": result["global_step"],
            "variant": "TWIRL-FM0.1.1",
            "architecture": "tcn",
            "parameter_count": sum(
                parameter.numel() for parameter in model.parameters()
            ),
            "final_metrics": result["final_metrics"],
            "run_contract_sha256": contract_hash,
            "checkpoint_sha256": checkpoint_hash,
        },
    )


def _rewrite_checkpoint_binding(run_dir: Path, checkpoint: dict) -> None:
    torch.save(checkpoint, run_dir / "checkpoint.pt")
    checkpoint_hash = write_sha256_sidecar(run_dir / "checkpoint.pt")
    summary = read_json(run_dir / "summary.json")
    summary["checkpoint_sha256"] = checkpoint_hash
    write_json_with_sha256(run_dir / "summary.json", summary)


def _assert_nested_equal(first, second, *, context: str = "root") -> None:
    if isinstance(first, torch.Tensor):
        assert isinstance(second, torch.Tensor), context
        assert torch.equal(first, second), context
        return
    if isinstance(first, np.ndarray):
        assert isinstance(second, np.ndarray), context
        assert np.array_equal(first, second), context
        return
    if isinstance(first, Mapping):
        assert isinstance(second, Mapping), context
        assert set(first) == set(second), context
        for key in first:
            _assert_nested_equal(
                first[key], second[key], context=f"{context}.{key}"
            )
        return
    if isinstance(first, (list, tuple)):
        assert isinstance(second, type(first)), context
        assert len(first) == len(second), context
        for index, (left, right) in enumerate(zip(first, second)):
            _assert_nested_equal(
                left, right, context=f"{context}[{index}]"
            )
        return
    assert first == second, context


def test_masked_huber_uses_equal_view_weight() -> None:
    prediction = torch.tensor([[[2.0, 0.0], [4.0, 4.0]]])
    target = torch.zeros_like(prediction)
    mask = torch.tensor([[[True, False], [True, True]]])
    loss, stats = masked_huber_reconstruction_loss(
        prediction, target, mask, delta=1.0
    )
    # view 0: 1.5; view 1: 3.5; equal-view mean: 2.5
    assert torch.isclose(loss, torch.tensor(2.5))
    assert stats["target_counts"].tolist() == [1, 2]


def test_vicreg_frozen_weighted_loss_backpropagates() -> None:
    first = torch.randn(8, 6, requires_grad=True)
    second = torch.randn(8, 6, requires_grad=True)
    loss, diagnostics = vicreg_loss(first, second)
    loss.backward()
    assert torch.isfinite(loss)
    assert torch.isfinite(first.grad).all()
    assert set(diagnostics) == {
        "loss", "invariance", "variance", "covariance"
    }


@pytest.mark.parametrize("invalid", (-1.0, float("nan"), float("inf")))
def test_optimization_rejects_invalid_vicreg_weights(invalid: float) -> None:
    with pytest.raises(ValueError, match="VICReg weights"):
        FM0OptimizationConfig(vicreg_total_weight=invalid)


def test_vicreg_connects_the_otherwise_untrained_embedding_projection() -> None:
    dataset = SyntheticFM0Dataset(
        SyntheticFM0Config(
            variant="TWIRL-FM0.1.1",
            seed=560067,
            windows_per_epoch=4,
            window_length=32,
        )
    )
    first = collate_fm0_samples(
        [dataset.sample(index, mask_view=0) for index in range(4)]
    )
    second = collate_fm0_samples(
        [dataset.sample(index, mask_view=1) for index in range(4)]
    )
    optimization = FM0OptimizationConfig(
        effective_batch_windows=4,
        max_optimizer_steps=1,
        warmup_steps=0,
        vicreg_total_weight=2.0e-4,
    )

    reconstruction_model = _tiny_model()
    reconstruction, _ = fm0_objective(
        reconstruction_model,
        first,
        second_batch=None,
        config=optimization,
    )
    reconstruction.backward()
    assert all(
        parameter.grad is None
        for parameter in reconstruction_model.embedding_projection.parameters()
    )

    repaired_model = _tiny_model()
    repaired, diagnostics = fm0_objective(
        repaired_model,
        first,
        second_batch=second,
        config=optimization,
        reconstruct_second_view=False,
    )
    repaired.backward()
    projection_gradients = [
        parameter.grad
        for parameter in repaired_model.embedding_projection.parameters()
    ]
    assert all(gradient is not None for gradient in projection_gradients)
    assert sum(float(gradient.detach().pow(2).sum()) for gradient in projection_gradients) > 0
    assert diagnostics["vicreg_variance"] > 0
    assert torch.equal(
        diagnostics["reconstruction"], diagnostics["reconstruction_first"]
    )
    assert torch.allclose(
        diagnostics["reconstruction_mean"],
        0.5
        * (
            diagnostics["reconstruction_first"]
            + diagnostics["reconstruction_second"]
        ),
    )


def test_vicreg_replay_matches_retained_graph_with_dropout() -> None:
    dataset = SyntheticFM0Dataset(
        SyntheticFM0Config(
            variant="TWIRL-FM0.1.1",
            seed=560067,
            windows_per_epoch=4,
            window_length=32,
        )
    )
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
    template = _tiny_model()
    retained_model = copy.deepcopy(template).train()
    replay_model = copy.deepcopy(template).train()

    torch.manual_seed(991)
    retained_first = []
    retained_second = []
    for first_batch, second_batch in paired_batches:
        retained_first.append(retained_model(first_batch)["z_window"])
        retained_second.append(retained_model(second_batch)["z_window"])
    retained_loss, retained_diagnostics = vicreg_loss(
        torch.cat(retained_first, dim=0).float(),
        torch.cat(retained_second, dim=0).float(),
        invariance_weight=optimization.vicreg_invariance_weight,
        variance_weight=optimization.vicreg_variance_weight,
        covariance_weight=optimization.vicreg_covariance_weight,
    )
    (optimization.vicreg_total_weight * retained_loss).backward()

    torch.manual_seed(991)
    replay_records = []
    replay_first = []
    replay_second = []
    for first_batch, second_batch in paired_batches:
        first_rng = fm0_training._capture_forward_rng_state(torch.device("cpu"))
        first_z = replay_model(first_batch)["z_window"]
        second_rng = fm0_training._capture_forward_rng_state(torch.device("cpu"))
        second_z = replay_model(second_batch)["z_window"]
        replay_first.append(first_z.detach())
        replay_second.append(second_z.detach())
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
        fm0_training._backprop_effective_batch_vicreg(
            model=replay_model,
            replay_records=replay_records,
            first_embeddings=replay_first,
            second_embeddings=replay_second,
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


def test_vicreg_uses_full_effective_batch_across_microbatches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    observed_batch_sizes: list[int] = []
    original_vicreg = fm0_training.vicreg_loss

    def recording_vicreg(first, second, **kwargs):
        observed_batch_sizes.append(int(first.shape[0]))
        return original_vicreg(first, second, **kwargs)

    monkeypatch.setattr(fm0_training, "vicreg_loss", recording_vicreg)
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=0,
        max_optimizer_steps=1,
        effective_batch_windows=64,
    )
    seed_everything(560067)
    result = run_synthetic_training(
        model=_tiny_vicreg_model(),
        dataset=SyntheticFM0Dataset(
            SyntheticFM0Config(
                variant="TWIRL-FM0.1.5",
                seed=560067,
                windows_per_epoch=64,
                window_length=32,
            )
        ),
        output_dir=tmp_path,
        run_contract={"schema_version": "test", "seed": 560067},
        optimization=optimization,
        target_step=1,
        micro_batch_windows=8,
        device="cpu",
        precision="fp32",
        use_vicreg=True,
    )
    assert observed_batch_sizes == [64]
    metrics = result["final_metrics"]
    for name in (
        "reconstruction",
        "reconstruction_first",
        "reconstruction_second",
        "reconstruction_mean",
        "vicreg",
        "vicreg_weighted",
        "vicreg_invariance",
        "vicreg_variance",
        "vicreg_covariance",
        "embedding_projection_gradient_norm",
        "total",
    ):
        assert math.isfinite(metrics[name])
    assert metrics["embedding_projection_gradient_norm"] > 0


def test_training_rejects_objective_mode_drift(tmp_path: Path) -> None:
    optimization = FM0OptimizationConfig(
        warmup_steps=0,
        max_optimizer_steps=1,
        effective_batch_windows=2,
    )
    with pytest.raises(ValueError, match="objective path disagree"):
        run_synthetic_training(
            model=_tiny_model(),
            dataset=SyntheticFM0Dataset(
                SyntheticFM0Config(
                    variant="TWIRL-FM0.1.1",
                    windows_per_epoch=2,
                    window_length=32,
                )
            ),
            output_dir=tmp_path,
            run_contract={"use_vicreg": False},
            optimization=optimization,
            target_step=1,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=True,
        )


def test_checkpoint_resume_rejects_objective_mode_drift(tmp_path: Path) -> None:
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=0,
        max_optimizer_steps=2,
        effective_batch_windows=2,
    )
    dataset_config = SyntheticFM0Config(
        variant="TWIRL-FM0.1.1",
        seed=560067,
        windows_per_epoch=4,
        window_length=32,
    )
    contract = {"schema_version": "legacy-test", "seed": 560067}
    seed_everything(560067)
    run_synthetic_training(
        model=_tiny_model(),
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=tmp_path,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=False,
    )

    seed_everything(999)
    with pytest.raises(ValueError, match="checkpoint objective state mismatch"):
        run_synthetic_training(
            model=_tiny_model(),
            dataset=SyntheticFM0Dataset(dataset_config),
            output_dir=tmp_path,
            run_contract=contract,
            optimization=optimization,
            target_step=2,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=True,
            reconstruct_second_view=False,
            resume_checkpoint=tmp_path / "checkpoint.pt",
        )


def test_checkpoint_resume_rejects_schedule_optimizer_and_dataset_drift(
    tmp_path: Path,
) -> None:
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=0,
        max_optimizer_steps=2,
        effective_batch_windows=2,
    )
    dataset_config = SyntheticFM0Config(
        variant="TWIRL-FM0.1.1",
        seed=560067,
        windows_per_epoch=4,
        window_length=32,
    )
    contract = {
        "schema_version": "test",
        "seed": 560067,
        "use_vicreg": True,
        "reconstruct_second_view": False,
        "immutable_milestone_steps": [0, 1, 2],
    }
    seed_everything(560067)
    run_synthetic_training(
        model=_tiny_model(),
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=tmp_path,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=True,
        reconstruct_second_view=False,
    )
    checkpoint = tmp_path / "checkpoint.pt"

    drifted_contract = dict(contract)
    drifted_contract["immutable_milestone_steps"] = [0, 1, 3]
    with pytest.raises(
        ValueError,
        match="immutable milestone schedule exceeds optimizer horizon",
    ):
        run_synthetic_training(
            model=_tiny_model(),
            dataset=SyntheticFM0Dataset(dataset_config),
            output_dir=tmp_path,
            run_contract=drifted_contract,
            optimization=optimization,
            target_step=2,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=True,
            reconstruct_second_view=False,
            resume_checkpoint=checkpoint,
        )

    drifted_optimization = FM0OptimizationConfig(
        learning_rate=2.0e-3,
        weight_decay=0.0,
        warmup_steps=0,
        max_optimizer_steps=2,
        effective_batch_windows=2,
    )
    with pytest.raises(ValueError, match="optimization contract mismatch"):
        run_synthetic_training(
            model=_tiny_model(),
            dataset=SyntheticFM0Dataset(dataset_config),
            output_dir=tmp_path,
            run_contract=contract,
            optimization=drifted_optimization,
            target_step=2,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=True,
            reconstruct_second_view=False,
            resume_checkpoint=checkpoint,
        )

    drifted_dataset = SyntheticFM0Config(
        variant="TWIRL-FM0.1.1",
        seed=560068,
        windows_per_epoch=4,
        window_length=32,
    )
    with pytest.raises(ValueError, match="dataset configuration mismatch"):
        run_synthetic_training(
            model=_tiny_model(),
            dataset=SyntheticFM0Dataset(drifted_dataset),
            output_dir=tmp_path,
            run_contract=contract,
            optimization=optimization,
            target_step=2,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=True,
            reconstruct_second_view=False,
            resume_checkpoint=checkpoint,
        )


def test_release_validator_strictly_loads_checkpoint(tmp_path: Path) -> None:
    run_dir = tmp_path / "valid"
    _write_structural_release(run_dir)
    result = validate_run_release(run_dir)
    assert result["checkpoint_inspected"] is True
    assert result["checkpoint_tensor_count"] > 0


@pytest.mark.parametrize("corruption", ("missing_model_key", "nonfinite_model"))
def test_release_validator_rejects_structural_checkpoint_corruption(
    tmp_path: Path, corruption: str
) -> None:
    run_dir = tmp_path / corruption
    _write_structural_release(run_dir)
    checkpoint = torch.load(run_dir / "checkpoint.pt", weights_only=False)
    first_name = next(iter(checkpoint["model_state"]))
    if corruption == "missing_model_key":
        del checkpoint["model_state"][first_name]
    else:
        value = checkpoint["model_state"][first_name].clone()
        value.reshape(-1)[0] = float("inf")
        checkpoint["model_state"][first_name] = value
    _rewrite_checkpoint_binding(run_dir, checkpoint)
    with pytest.raises(ValueError, match="model state|non-finite"):
        validate_run_release(run_dir)


@pytest.mark.parametrize("use_vicreg", (False, True))
def test_checkpoint_resume_matches_uninterrupted_training(
    tmp_path: Path, use_vicreg: bool
) -> None:
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=1,
        max_optimizer_steps=2,
        effective_batch_windows=2,
    )
    dataset_config = SyntheticFM0Config(
        variant="TWIRL-FM0.1.1",
        seed=560067,
        windows_per_epoch=8,
        window_length=32,
    )
    contract = {
        "schema_version": "test",
        "seed": 560067,
        "use_vicreg": use_vicreg,
        "reconstruct_second_view": False,
        "immutable_milestone_steps": [0, 1, 2],
    }

    full_dir = tmp_path / "full"
    full_dir.mkdir()
    seed_everything(560067)
    full_model = _tiny_model()
    run_synthetic_training(
        model=full_model,
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=full_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=2,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=use_vicreg,
        reconstruct_second_view=False,
    )

    resumed_dir = tmp_path / "resumed"
    resumed_dir.mkdir()
    seed_everything(560067)
    partial_model = _tiny_model()
    run_synthetic_training(
        model=partial_model,
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=resumed_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=use_vicreg,
        reconstruct_second_view=False,
    )
    seed_everything(999)
    resumed_model = _tiny_model()
    run_synthetic_training(
        model=resumed_model,
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=resumed_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=2,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=use_vicreg,
        reconstruct_second_view=False,
        resume_checkpoint=resumed_dir / "checkpoint_step_00000001.pt",
    )
    full_state = torch.load(full_dir / "checkpoint.pt", weights_only=False)
    resumed_state = torch.load(resumed_dir / "checkpoint.pt", weights_only=False)
    _assert_nested_equal(full_state, resumed_state)
    expected_objective = {
        "use_vicreg": use_vicreg,
        "reconstruct_second_view": False,
    }
    assert full_state["objective_state"] == expected_objective
    for directory in (full_dir, resumed_dir):
        for step in (0, 1, 2):
            milestone = directory / f"checkpoint_step_{step:08d}.pt"
            assert milestone.is_file()
            sidecar = milestone.with_name(milestone.name + ".sha256")
            assert sidecar.is_file()
            assert fm0_training._verify_sha256_sidecar(
                milestone, required=True
            )
            payload = torch.load(milestone, weights_only=False)
            assert payload["progress"]["global_step"] == step
            assert payload["objective_state"] == expected_objective


def test_checkpoint_resume_rejects_wrong_required_step(tmp_path: Path) -> None:
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=0,
        max_optimizer_steps=2,
        effective_batch_windows=2,
    )
    dataset_config = SyntheticFM0Config(
        variant="TWIRL-FM0.1.1",
        seed=560067,
        windows_per_epoch=4,
        window_length=32,
    )
    contract = {
        "schema_version": "resume-step-test",
        "seed": 560067,
        "use_vicreg": False,
        "reconstruct_second_view": False,
        "immutable_milestone_steps": [0, 1, 2],
    }
    seed_everything(560067)
    run_synthetic_training(
        model=_tiny_model(),
        dataset=SyntheticFM0Dataset(dataset_config),
        output_dir=tmp_path,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=False,
    )

    with pytest.raises(ValueError, match="required resume step"):
        run_synthetic_training(
            model=_tiny_model(),
            dataset=SyntheticFM0Dataset(dataset_config),
            output_dir=tmp_path,
            run_contract=contract,
            optimization=optimization,
            target_step=2,
            micro_batch_windows=2,
            device="cpu",
            precision="fp32",
            use_vicreg=False,
            resume_checkpoint=tmp_path / "checkpoint_step_00000001.pt",
            expected_resume_step=2,
        )


def test_milestone_checkpoint_is_hash_bound_and_non_overwriting(
    tmp_path: Path,
) -> None:
    optimization = FM0OptimizationConfig(
        learning_rate=1.0e-3,
        weight_decay=0.0,
        warmup_steps=0,
        max_optimizer_steps=1,
        effective_batch_windows=2,
    )
    contract = {
        "schema_version": "test",
        "seed": 560067,
        "use_vicreg": False,
        "reconstruct_second_view": False,
        "immutable_milestone_steps": [0, 1],
    }
    seed_everything(560067)
    result = run_synthetic_training(
        model=_tiny_model(),
        dataset=SyntheticFM0Dataset(
            SyntheticFM0Config(
                variant="TWIRL-FM0.1.1",
                seed=560067,
                windows_per_epoch=2,
                window_length=32,
            )
        ),
        output_dir=tmp_path,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=2,
        device="cpu",
        precision="fp32",
        use_vicreg=False,
    )
    assert [item["step"] for item in result["immutable_milestone_checkpoints"]] == [
        0,
        1,
    ]
    milestone = tmp_path / "checkpoint_step_00000001.pt"
    payload = torch.load(milestone, weights_only=False)
    with pytest.raises(FileExistsError, match="refusing to replace"):
        fm0_training.save_immutable_milestone_checkpoint(tmp_path, payload)

    with milestone.open("ab") as handle:
        handle.write(b"tamper")
    resume_model = _tiny_model()
    resume_optimizer, resume_scheduler = (
        fm0_training.make_optimizer_and_scheduler(resume_model, optimization)
    )
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        fm0_training.load_checkpoint(
            milestone,
            model=resume_model,
            optimizer=resume_optimizer,
            scheduler=resume_scheduler,
            expected_run_contract=contract,
            expected_optimization=optimization,
            expected_objective_state={
                "use_vicreg": False,
                "reconstruct_second_view": False,
            },
            dataset=SyntheticFM0Dataset(
                SyntheticFM0Config(
                    variant="TWIRL-FM0.1.1",
                    seed=560067,
                    windows_per_epoch=2,
                    window_length=32,
                )
            ),
        )
