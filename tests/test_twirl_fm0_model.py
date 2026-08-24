from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from twirl.models.fm0.dataset import (  # noqa: E402
    SyntheticFM0Config,
    SyntheticFM0Dataset,
    collate_fm0_samples,
)
from twirl.models.fm0.model import (  # noqa: E402
    FM0ModelConfig,
    build_fm0_model,
    count_trainable_parameters,
    parameter_match_fraction,
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
    assert first["h_window"].shape == (2, 16)
    assert first["z_window"].shape == (2, 16)
    assert torch.equal(
        first["z_window"], model.embedding_projection(first["h_window"])
    )
    assert torch.equal(first["reconstruction"], second["reconstruction"])
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
