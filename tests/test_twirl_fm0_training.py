from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path

import pytest

torch = pytest.importorskip("torch")

from twirl.models.fm0.dataset import (  # noqa: E402
    SyntheticFM0Config,
    SyntheticFM0Dataset,
)
from twirl.models.fm0.model import FM0ModelConfig, build_fm0_model  # noqa: E402
import twirl.models.fm0.training as fm0_training  # noqa: E402
from twirl.models.fm0.training import (  # noqa: E402
    FM0OptimizationConfig,
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
    assert math.isfinite(result["final_metrics"]["vicreg"])
    assert math.isfinite(result["final_metrics"]["total"])


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


def test_checkpoint_resume_matches_uninterrupted_training(tmp_path: Path) -> None:
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
    contract = {"schema_version": "test", "seed": 560067}

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
        use_vicreg=False,
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
        use_vicreg=False,
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
        use_vicreg=False,
        resume_checkpoint=resumed_dir / "checkpoint.pt",
    )
    full_state = torch.load(full_dir / "checkpoint.pt", weights_only=False)
    resumed_state = torch.load(resumed_dir / "checkpoint.pt", weights_only=False)
    assert full_state["loss_history"] == resumed_state["loss_history"]
    for name, value in full_state["model_state"].items():
        assert torch.equal(value, resumed_state["model_state"][name]), name
