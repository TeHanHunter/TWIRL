from __future__ import annotations

from dataclasses import asdict

import pytest

torch = pytest.importorskip("torch")

from twirl.models.fm0.dataset import (
    SyntheticFM0Config,
    SyntheticFM0Dataset,
    collate_fm0_samples,
)
from twirl.models.fm0.model import (
    FM0ModelConfig,
    architecture_for_variant,
    build_fm0_model,
    count_trainable_parameters,
    parameter_match_fraction,
)

_CADENCE_LAST_KEYS = (
    "flux",
    "flux_valid",
    "flux_error",
    "error_valid",
    "local_time_cadences",
    "delta_time_cadences",
    "time_valid",
    "segment_boundary",
    "reconstruction_mask",
    "temporal_mask",
)


def _slice_batch_context(
    batch: dict[str, torch.Tensor], length: int
) -> dict[str, torch.Tensor]:
    sliced = {name: value.clone() for name, value in batch.items()}
    for name in _CADENCE_LAST_KEYS:
        if name in sliced:
            sliced[name] = sliced[name][..., :length].clone()
    return sliced


def _tiny_model_config(architecture: str) -> FM0ModelConfig:
    return FM0ModelConfig(
        architecture=architecture,
        n_flux_views=2,
        window_length=64,
        d_model=16,
        embedding_dim=16,
        stem_kernel=5,
        dropout=0.0,
        tcn_blocks=2,
        tcn_dilation_cycle=(1, 2),
        conformer_blocks=1,
        conformer_heads=4,
        conformer_ff_multiplier=2,
        conformer_conv_kernel=7,
        patch_stride=4,
        minimum_parameters=1,
        maximum_parameters=10_000_000,
    )


def test_default_tcn_and_conformer_are_in_frozen_matched_budget() -> None:
    tcn = build_fm0_model("TWIRL-FM0.1.1")
    conformer = build_fm0_model("TWIRL-FM0.1.2")
    for model in (tcn, conformer):
        assert 8_000_000 <= count_trainable_parameters(model) <= 12_000_000
    assert parameter_match_fraction(tcn, conformer) <= 0.10


def test_fm0_2_objective_variants_reuse_the_exact_fm0_1_backbones() -> None:
    fm01_tcn = build_fm0_model("TWIRL-FM0.1.1")
    fm02_tcn = build_fm0_model("TWIRL-FM0.2.1")
    fm01_conformer = build_fm0_model("TWIRL-FM0.1.2")
    fm02_conformer = build_fm0_model("TWIRL-FM0.2.2")

    assert fm01_tcn.config == fm02_tcn.config
    assert fm01_conformer.config == fm02_conformer.config


def test_fm0_3_variant_mappings_freeze_short_cadence_preserving_geometry() -> None:
    assert architecture_for_variant("TWIRL-FM0.3.1") == "tcn"
    assert architecture_for_variant("TWIRL-FM0.3.2") == "conformer"

    tcn = build_fm0_model("TWIRL-FM0.3.1")
    conformer = build_fm0_model("TWIRL-FM0.3.2")
    assert tcn.config.architecture == "tcn"
    assert conformer.config.architecture == "conformer"
    for model in (tcn, conformer):
        assert model.config.n_flux_views == 2
        assert model.config.window_length == 128
        assert model.config.patch_stride == 1


def test_legacy_model_config_checkpoint_keys_are_unchanged() -> None:
    assert set(asdict(FM0ModelConfig(architecture="tcn", n_flux_views=2))) == {
        "architecture",
        "n_flux_views",
        "window_length",
        "d_model",
        "embedding_dim",
        "stem_kernel",
        "dropout",
        "tcn_blocks",
        "tcn_kernel",
        "tcn_dilation_cycle",
        "conformer_blocks",
        "conformer_heads",
        "conformer_ff_multiplier",
        "conformer_conv_kernel",
        "patch_stride",
        "minimum_parameters",
        "maximum_parameters",
    }


@pytest.mark.parametrize(
    ("variant", "architecture"),
    (("TWIRL-FM0.1.1", "tcn"), ("TWIRL-FM0.1.2", "conformer")),
)
def test_tiny_models_preserve_shapes_and_ignore_invalid_fill(
    variant: str, architecture: str
) -> None:
    config = FM0ModelConfig(
        architecture=architecture,
        n_flux_views=2,
        window_length=64,
        d_model=16,
        embedding_dim=16,
        stem_kernel=5,
        dropout=0.0,
        tcn_blocks=2,
        tcn_dilation_cycle=(1, 2),
        conformer_blocks=1,
        conformer_heads=4,
        conformer_ff_multiplier=2,
        conformer_conv_kernel=7,
        patch_stride=4,
        minimum_parameters=1,
        maximum_parameters=10_000_000,
    )
    model = build_fm0_model(
        variant, config_override=config, enforce_parameter_budget=False
    ).eval()
    dataset = SyntheticFM0Dataset(
        SyntheticFM0Config(variant=variant, window_length=64, windows_per_epoch=4)
    )
    batch = collate_fm0_samples([dataset[0], dataset[1]])
    with torch.no_grad():
        first = model(batch)
        changed = {name: value.clone() for name, value in batch.items()}
        changed["flux"][~changed["flux_valid"]] = 1.0e6
        second = model(changed)
    assert first["reconstruction"].shape == (2, 2, 64)
    assert first["h_cadence"].shape == (2, 64, 16)
    assert first["h_window"].shape == (2, 16)
    assert first["z_window"].shape == (2, 16)
    assert torch.equal(
        first["z_window"], model.embedding_projection(first["h_window"])
    )
    assert torch.equal(first["reconstruction"], second["reconstruction"])
    assert torch.equal(first["h_cadence"], second["h_cadence"])
    assert torch.equal(first["h_window"], second["h_window"])
    assert torch.equal(first["z_window"], second["z_window"])


def test_conformer_reconstruction_retains_subpatch_cadence_state() -> None:
    config = FM0ModelConfig(
        architecture="conformer",
        n_flux_views=2,
        window_length=64,
        d_model=16,
        embedding_dim=16,
        stem_kernel=5,
        dropout=0.0,
        conformer_blocks=1,
        conformer_heads=4,
        conformer_ff_multiplier=2,
        conformer_conv_kernel=7,
        minimum_parameters=1,
        maximum_parameters=10_000_000,
    )
    model = build_fm0_model(
        "TWIRL-FM0.1.2", config_override=config, enforce_parameter_budget=False
    ).eval()
    dataset = SyntheticFM0Dataset(
        SyntheticFM0Config(
            variant="TWIRL-FM0.1.2", window_length=64, windows_per_epoch=2
        )
    )
    batch = collate_fm0_samples([dataset[0]])
    with torch.no_grad():
        reconstruction = model(batch)["reconstruction"]
    assert not torch.equal(reconstruction[..., 0], reconstruction[..., 1])


def test_cadence_preserving_conformer_keeps_one_token_per_cadence() -> None:
    config = FM0ModelConfig(
        architecture="conformer",
        n_flux_views=2,
        window_length=64,
        d_model=16,
        embedding_dim=16,
        stem_kernel=5,
        dropout=0.0,
        conformer_blocks=1,
        conformer_heads=4,
        conformer_ff_multiplier=2,
        conformer_conv_kernel=7,
        patch_stride=1,
        minimum_parameters=1,
        maximum_parameters=10_000_000,
    )
    model = build_fm0_model(
        "TWIRL-FM0.1.2", config_override=config, enforce_parameter_budget=False
    ).eval()
    dataset = SyntheticFM0Dataset(
        SyntheticFM0Config(
            variant="TWIRL-FM0.1.2", window_length=64, windows_per_epoch=2
        )
    )
    batch = collate_fm0_samples([dataset[0]])

    with torch.no_grad():
        stem = model.stem(*model._pack_input(batch)[:2])
        token_valid = batch["time_valid"].bool()
        patched, patched_valid = model._patch_mean(stem, token_valid)
        output = model(batch)

    assert patched.shape[-1] == batch["flux"].shape[-1]
    assert patched_valid.shape[-1] == batch["flux"].shape[-1]
    assert patched is stem
    assert patched_valid is token_valid
    torch.testing.assert_close(patched, stem)
    assert torch.equal(patched_valid, batch["time_valid"].bool())
    assert output["token_valid"].shape == batch["time_valid"].shape
    assert output["h_cadence"].shape[1] == batch["flux"].shape[-1]


def test_model_config_rejects_nonpositive_patch_stride() -> None:
    with pytest.raises(ValueError, match="patch_stride must be positive"):
        FM0ModelConfig(architecture="conformer", n_flux_views=2, patch_stride=0)


@pytest.mark.parametrize(
    ("variant", "architecture"),
    (("TWIRL-FM0.1.1", "tcn"), ("TWIRL-FM0.1.2", "conformer")),
)
def test_short_context_eval_matches_an_exact_length_masked_tail(
    variant: str, architecture: str
) -> None:
    config = _tiny_model_config(architecture)
    model = build_fm0_model(
        variant, config_override=config, enforce_parameter_budget=False
    ).eval()
    dataset = SyntheticFM0Dataset(
        SyntheticFM0Config(variant=variant, window_length=64, windows_per_epoch=4)
    )
    batch = collate_fm0_samples([dataset[0], dataset[1]])
    short_length = 48
    short_batch = _slice_batch_context(batch, short_length)
    masked_tail_batch = {name: value.clone() for name, value in batch.items()}
    masked_tail_batch["flux_valid"][..., short_length:] = False
    masked_tail_batch["error_valid"][..., short_length:] = False
    masked_tail_batch["time_valid"][..., short_length:] = False
    masked_tail_batch["reconstruction_mask"][..., short_length:] = False

    with pytest.raises(
        ValueError, match="flux cadence length differs from model configuration"
    ):
        model(short_batch)

    short_output = model.forward_short_context(short_batch)
    with torch.no_grad():
        masked_tail_output = model(masked_tail_batch)

    assert short_output["reconstruction"].shape == (2, 2, short_length)
    assert short_output["token_valid"].shape[-1] == (
        short_length
        if architecture == "tcn"
        else short_length // config.patch_stride
    )
    assert all(not value.requires_grad for value in short_output.values())
    torch.testing.assert_close(
        short_output["reconstruction"],
        masked_tail_output["reconstruction"][..., :short_length],
    )
    torch.testing.assert_close(
        short_output["h_window"], masked_tail_output["h_window"]
    )
    torch.testing.assert_close(
        short_output["z_window"], masked_tail_output["z_window"]
    )


def test_short_context_forward_is_eval_only_and_checkpoint_neutral() -> None:
    config = _tiny_model_config("tcn")
    model = build_fm0_model(
        "TWIRL-FM0.1.1", config_override=config, enforce_parameter_budget=False
    )
    dataset = SyntheticFM0Dataset(
        SyntheticFM0Config(
            variant="TWIRL-FM0.1.1", window_length=64, windows_per_epoch=2
        )
    )
    batch = collate_fm0_samples([dataset[0]])
    short_batch = _slice_batch_context(batch, 32)

    with pytest.raises(RuntimeError, match="short-context forward is evaluation-only"):
        model.forward_short_context(short_batch)

    model.eval()
    full_eval_output = model.forward_short_context(batch)
    with torch.no_grad():
        full_output = model(batch)
    for name in full_output:
        torch.testing.assert_close(full_eval_output[name], full_output[name])

    oversized_batch = {
        name: (
            torch.cat((value, value[..., -1:]), dim=-1)
            if name in _CADENCE_LAST_KEYS
            else value.clone()
        )
        for name, value in batch.items()
    }
    with pytest.raises(ValueError, match="no greater than the checkpoint window"):
        model.forward_short_context(oversized_batch)

    reloaded = build_fm0_model(
        "TWIRL-FM0.1.1", config_override=config, enforce_parameter_budget=False
    )
    reloaded.load_state_dict(model.state_dict(), strict=True)
    assert tuple(reloaded.state_dict()) == tuple(model.state_dict())
