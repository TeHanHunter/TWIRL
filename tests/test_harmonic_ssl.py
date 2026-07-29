from __future__ import annotations

import copy
import math

import pandas as pd
import pytest

from twirl.vetting.harmonic_inputs import HARMONIC_FACTORS
from twirl.vetting.harmonic_ssl import (
    EventPreservingAugmentationConfig,
    VICRegConfig,
    _event_phase_centers,
    augment_ssl_batch,
    build_ssl_cache_identity,
    select_ssl_fold_rows,
    validate_ssl_cache_identity,
    vicreg_loss,
)


def _selection_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "review_id": "held",
                "tic": 100,
                "sector": 56,
                "fixed_split": "development",
                "cv_fold": 0,
                "is_injected_row": False,
            },
            {
                "review_id": "selected-a",
                "tic": 200,
                "sector": 57,
                "fixed_split": "development",
                "cv_fold": 1,
                "is_injected_row": False,
            },
            {
                "review_id": "selected-b",
                "tic": 200,
                "sector": 58,
                "fixed_split": "development",
                "cv_fold": 1,
                "is_injected_row": False,
            },
            {
                "review_id": "fixed-test",
                "tic": 300,
                "sector": 59,
                "fixed_split": "test",
                "cv_fold": -1,
                "is_injected_row": False,
            },
            {
                "review_id": "injected",
                "tic": 400,
                "sector": 56,
                "fixed_split": "development",
                "cv_fold": 2,
                "is_injected_row": True,
            },
            {
                "review_id": "reserved-host-earlier",
                "tic": 500,
                "sector": 62,
                "fixed_split": "development",
                "cv_fold": 3,
                "is_injected_row": False,
            },
            {
                "review_id": "reserved-host-s63",
                "tic": 500,
                "sector": 63,
                "fixed_split": "development",
                "cv_fold": 3,
                "is_injected_row": False,
            },
            {
                "review_id": "selected-c",
                "tic": 600,
                "sector": 62,
                "fixed_split": "development",
                "cv_fold": 4,
                "is_injected_row": "false",
            },
        ]
    )


def test_ssl_fold_selection_is_real_only_and_host_disjoint() -> None:
    rows = _selection_rows()
    selected, audit = select_ssl_fold_rows(rows, held_out_fold=0)

    assert selected["review_id"].tolist() == ["selected-a", "selected-b", "selected-c"]
    assert set(selected["tic"]) == {200, 600}
    assert not selected["is_injected_row"].astype(str).str.lower().eq("true").any()
    assert audit["n_selected_rows"] == 3
    assert audit["n_selected_tics"] == 2
    assert audit["primary_ssl_real_only"] is True
    assert audit["reserved_sector_exclusion_scope"] == "whole_tic"
    assert all(audit["tic_disjoint"].values())
    assert len(audit["selected_rows_sha256"]) == 64
    assert len(audit["selected_tics_sha256"]) == 64


def test_ssl_fold_selection_fails_closed_on_split_or_boolean_errors() -> None:
    rows = _selection_rows()
    rows.loc[rows["tic"].eq(200) & rows["sector"].eq(58), "cv_fold"] = 2
    with pytest.raises(ValueError, match="TIC split leakage"):
        select_ssl_fold_rows(rows, held_out_fold=0)

    rows = _selection_rows()
    rows.loc[rows["review_id"].eq("injected"), "is_injected_row"] = "maybe"
    with pytest.raises(ValueError, match="invalid boolean"):
        select_ssl_fold_rows(rows, held_out_fold=0)

    with pytest.raises(KeyError, match="missing required columns"):
        select_ssl_fold_rows(_selection_rows().drop(columns="sector"), held_out_fold=0)


def _native_batch(torch: object) -> dict[str, object]:
    batch_size = 2
    n_harmonic = len(HARMONIC_FACTORS)
    n_samples = 11
    phase = torch.linspace(-0.5, 0.5, n_samples)
    harmonic = torch.zeros(batch_size, n_harmonic, 7, n_samples)
    harmonic[:, :, 0, :] = 1.0
    harmonic[:, :, 1, :] = 2.0
    harmonic[:, :, 2, :] = 1.0
    harmonic[:, :, 3, :] = 0.2
    harmonic[:, :, 4, :] = 0.3
    harmonic[:, :, 5, :] = phase
    harmonic[:, :, 6, :] = 0.0
    harmonic_mask = torch.ones_like(harmonic, dtype=torch.bool)

    local = harmonic[:, :1].repeat(1, 14, 1, 1)
    local_mask = torch.ones_like(local, dtype=torch.bool)
    local_mask[0, 0, 0, 0] = False
    periodogram = torch.arange(4 * 13, dtype=torch.float32).reshape(1, 4, 13)
    periodogram = periodogram.repeat(batch_size, 1, 1)
    return {
        "review_id": ["a", "b"],
        "tic": torch.tensor([1, 2]),
        "period_d": torch.tensor([1.0, 2.0]),
        "harmonic_values": harmonic,
        "harmonic_mask": harmonic_mask,
        "local_values": local,
        "local_mask": local_mask,
        "periodogram_values": periodogram,
        "periodogram_mask": torch.ones_like(periodogram, dtype=torch.bool),
        "metadata": torch.ones(batch_size, 5),
    }


def test_ssl_augmentation_is_deterministic_and_does_not_mutate_input() -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    original = {
        key: value.clone() if torch.is_tensor(value) else copy.deepcopy(value)
        for key, value in batch.items()
    }
    config = EventPreservingAugmentationConfig(
        harmonic_cadence_dropout_probability=0.35,
        harmonic_flux_noise_scale=0.2,
        local_flux_noise_scale=0.1,
        periodogram_bin_dropout_probability=0.25,
    )

    first = augment_ssl_batch(
        batch,
        duration_min=torch.tensor([120.0, 180.0]),
        config=config,
        seed=56,
        view_index=0,
    )
    repeated = augment_ssl_batch(
        batch,
        duration_min=torch.tensor([120.0, 180.0]),
        config=config,
        seed=56,
        view_index=0,
    )
    second_view = augment_ssl_batch(
        batch,
        duration_min=torch.tensor([120.0, 180.0]),
        config=config,
        seed=56,
        view_index=1,
    )

    for name in (
        "harmonic_values",
        "harmonic_mask",
        "local_values",
        "local_mask",
        "periodogram_values",
        "periodogram_mask",
        "metadata",
    ):
        assert torch.equal(first[name], repeated[name])
    assert not torch.equal(first["periodogram_values"], second_view["periodogram_values"])
    assert torch.equal(first["local_mask"], batch["local_mask"])
    assert not (first["periodogram_mask"] & ~batch["periodogram_mask"]).any()
    assert torch.count_nonzero(first["metadata"]) == 0
    for values, mask in (
        (first["harmonic_values"], first["harmonic_mask"]),
        (first["local_values"], first["local_mask"]),
    ):
        coherent = mask[:, :, 0, :] & mask[:, :, 1, :]
        assert torch.allclose(
            values[:, :, 2, :][coherent],
            (values[:, :, 1, :] - values[:, :, 0, :])[coherent],
        )
    for key, value in original.items():
        if torch.is_tensor(value):
            assert torch.equal(batch[key], value)
        else:
            assert batch[key] == value


@pytest.mark.parametrize(
    ("tensor_name", "mask_name", "index"),
    (
        ("harmonic_values", "harmonic_mask", (0, 0, 0, 1)),
        ("local_values", "local_mask", (0, 0, 0, 1)),
        ("periodogram_values", "periodogram_mask", (0, 0, 1)),
        ("metadata", None, (0, 1)),
    ),
)
@pytest.mark.parametrize(
    ("invalid_value", "message"),
    (
        (float("inf"), "nonfinite model-active values"),
        (1_000_001.0, "exceeds the model-active absolute bound"),
    ),
)
def test_ssl_augmentation_rejects_invalid_active_model_values(
    tensor_name: str,
    mask_name: str | None,
    index: tuple[int, ...],
    invalid_value: float,
    message: str,
) -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    batch[tensor_name][index] = invalid_value

    with pytest.raises(ValueError, match=message):
        augment_ssl_batch(
            batch,
            duration_min=torch.tensor([120.0, 180.0]),
            config=EventPreservingAugmentationConfig(
                harmonic_cadence_dropout_probability=0.0,
                harmonic_flux_noise_scale=0.0,
                local_flux_noise_scale=0.0,
                periodogram_bin_dropout_probability=0.0,
            ),
            seed=56,
            view_index=0,
        )

    if mask_name is not None and math.isfinite(invalid_value):
        batch[mask_name][index] = False
        augment_ssl_batch(
            batch,
            duration_min=torch.tensor([120.0, 180.0]),
            config=EventPreservingAugmentationConfig(
                harmonic_cadence_dropout_probability=0.0,
                harmonic_flux_noise_scale=0.0,
                local_flux_noise_scale=0.0,
                periodogram_bin_dropout_probability=0.0,
            ),
            seed=56,
            view_index=0,
        )


def test_ssl_augmentation_rejects_masked_nonfinite_tensor_payload() -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    batch["harmonic_values"][0, 0, 0, 1] = float("nan")
    batch["harmonic_mask"][0, 0, 0, 1] = False

    with pytest.raises(ValueError, match="nonfinite tensor payload"):
        augment_ssl_batch(
            batch,
            duration_min=torch.tensor([120.0, 180.0]),
            seed=56,
            view_index=0,
        )


@pytest.mark.parametrize(
    ("channel", "invalid_value", "message"),
    (
        (3, -0.1, "error channels must be nonnegative"),
        (5, 0.5001, "phase channel must be in"),
        (6, 0.5, "quality channel must be binary"),
    ),
)
def test_ssl_augmentation_rejects_harmonic_channel_domain_violations(
    channel: int,
    invalid_value: float,
    message: str,
) -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    batch["harmonic_values"][0, 0, channel, 1] = invalid_value

    with pytest.raises(ValueError, match=message):
        augment_ssl_batch(
            batch,
            duration_min=torch.tensor([120.0, 180.0]),
            seed=56,
            view_index=0,
        )


def test_ssl_augmentation_requires_quality_bad_photometry_mask() -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    cadence = 3
    for values_name, mask_name in (
        ("harmonic_values", "harmonic_mask"),
        ("local_values", "local_mask"),
    ):
        batch[values_name][0, 0, 6, cadence] = 1.0
        batch[mask_name][0, 0, :5, cadence] = False

    config = EventPreservingAugmentationConfig(
        harmonic_cadence_dropout_probability=0.0,
        harmonic_flux_noise_scale=0.0,
        local_flux_noise_scale=0.0,
        periodogram_bin_dropout_probability=0.0,
    )
    augmented = augment_ssl_batch(
        batch,
        duration_min=torch.tensor([120.0, 180.0]),
        config=config,
        seed=56,
        view_index=0,
    )
    for values_name, mask_name in (
        ("harmonic_values", "harmonic_mask"),
        ("local_values", "local_mask"),
    ):
        assert torch.equal(
            augmented[values_name][:, :, 5:7, :],
            batch[values_name][:, :, 5:7, :],
        )
        assert torch.equal(
            augmented[mask_name][:, :, 5:7, :],
            batch[mask_name][:, :, 5:7, :],
        )

    invalid = _native_batch(torch)
    invalid["harmonic_values"][0, 0, 6, cadence] = 1.0
    with pytest.raises(ValueError, match="quality-bad cadences must mask"):
        augment_ssl_batch(
            invalid,
            duration_min=torch.tensor([120.0, 180.0]),
            seed=56,
            view_index=0,
        )


def test_ssl_augmentation_rejects_hidden_quality_coordinate_mask() -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    batch["harmonic_mask"][0, 0, 6, 1] = False

    with pytest.raises(ValueError, match="phase/quality masks must match"):
        augment_ssl_batch(
            batch,
            duration_min=torch.tensor([120.0, 180.0]),
            seed=56,
            view_index=0,
        )


def test_ssl_augmentation_rejects_post_augmentation_envelope_overflow() -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    batch["local_values"][:, :, 3:5, :] = 1_000_000.0

    with pytest.raises(
        ValueError,
        match="post-augmentation local_values exceeds",
    ):
        augment_ssl_batch(
            batch,
            duration_min=torch.tensor([120.0, 180.0]),
            config=EventPreservingAugmentationConfig(
                harmonic_cadence_dropout_probability=0.0,
                harmonic_flux_noise_scale=0.0,
                local_flux_noise_scale=10.0,
                periodogram_bin_dropout_probability=0.0,
            ),
            seed=56,
            view_index=0,
        )


def test_ssl_augmentation_protects_harmonic_event_samples() -> None:
    torch = pytest.importorskip("torch")
    batch = _native_batch(torch)
    config = EventPreservingAugmentationConfig(
        harmonic_cadence_dropout_probability=1.0,
        harmonic_flux_noise_scale=10.0,
        local_flux_noise_scale=0.0,
        periodogram_bin_dropout_probability=0.0,
        event_protection_duration_multiplier=2.0,
    )
    augmented = augment_ssl_batch(
        batch,
        duration_min=torch.tensor([144.0, 144.0]),
        config=config,
        seed=7,
        view_index=0,
    )

    phase = batch["harmonic_values"][:, :, 5, :]
    factors = torch.tensor(HARMONIC_FACTORS).view(1, -1, 1)
    half_width = (
        2.0
        * torch.tensor([144.0, 144.0]).view(-1, 1, 1)
        / (2.0 * 1440.0 * batch["period_d"].view(-1, 1, 1) * factors)
    ).clamp(max=0.25)
    protected = torch.zeros_like(phase, dtype=torch.bool)
    for view_index, factor in enumerate(HARMONIC_FACTORS):
        for center in _event_phase_centers(factor):
            distance = torch.abs(
                torch.remainder(
                    phase[:, view_index, :] - float(center) + 0.5,
                    1.0,
                )
                - 0.5
            )
            protected[:, view_index, :] |= distance <= half_width[:, view_index, :]

    protected_channels = protected.unsqueeze(2).expand_as(batch["harmonic_mask"])
    assert torch.equal(
        augmented["harmonic_mask"][protected_channels],
        batch["harmonic_mask"][protected_channels],
    )
    assert torch.equal(
        augmented["harmonic_values"][protected_channels],
        batch["harmonic_values"][protected_channels],
    )
    unprotected = (~protected).unsqueeze(2).expand_as(batch["harmonic_mask"])
    assert not augmented["harmonic_mask"][unprotected].any()
    assert torch.equal(augmented["local_mask"], batch["local_mask"])


def test_event_protection_includes_repeated_long_harmonic_aliases() -> None:
    assert _event_phase_centers(2.0) == (-0.5, -0.25, 0.0, 0.25)
    assert _event_phase_centers(3.0) == pytest.approx((
        -0.5,
        -1.0 / 3.0,
        -1.0 / 6.0,
        0.0,
        1.0 / 6.0,
        1.0 / 3.0,
    ))
    assert _event_phase_centers(4.0) == (
        -0.5,
        -0.375,
        -0.25,
        -0.125,
        0.0,
        0.125,
        0.25,
        0.375,
    )
    assert _event_phase_centers(0.5) == (0.0,)
    assert _event_phase_centers(1.0 / 3.0) == (-0.5, 0.0)


def test_vicreg_loss_reports_components_and_backpropagates() -> None:
    torch = pytest.importorskip("torch")
    torch.manual_seed(4)
    embedding_a = torch.randn(8, 6, requires_grad=True)
    embedding_b = (embedding_a.detach() + 0.1 * torch.randn(8, 6)).requires_grad_()
    loss, diagnostics = vicreg_loss(
        embedding_a,
        embedding_b,
        config=VICRegConfig(),
    )

    assert torch.isfinite(loss)
    assert set(diagnostics) >= {
        "loss",
        "invariance",
        "variance",
        "covariance",
        "std_mean_a",
        "std_mean_b",
    }
    assert torch.equal(loss, diagnostics["loss"])
    assert diagnostics["invariance"] > 0
    loss.backward()
    assert embedding_a.grad is not None
    assert embedding_b.grad is not None
    assert torch.isfinite(embedding_a.grad).all()

    collapsed = torch.zeros(4, 3)
    _, collapsed_diagnostics = vicreg_loss(collapsed, collapsed)
    assert collapsed_diagnostics["invariance"] == 0
    assert collapsed_diagnostics["variance"] > 0


def test_ssl_cache_identity_is_stable_and_strict() -> None:
    selected, _ = select_ssl_fold_rows(_selection_rows(), held_out_fold=0)
    hashes = {
        "training_table_sha256": "a" * 64,
        "native_registry_sha256": "b" * 64,
        "split_registry_sha256": "c" * 64,
    }
    identity = build_ssl_cache_identity(
        **hashes,
        selected_rows=selected,
        profile="shape_plus_periodogram_bls",
        held_out_fold=0,
        seed=56,
        model_config={"embedding_dim": 64, "metadata_dim": 31},
        code_revision="0123456789abcdef",
    )
    repeated = build_ssl_cache_identity(
        **hashes,
        selected_rows=selected.sample(frac=1.0, random_state=2),
        profile="shape_plus_periodogram_bls",
        held_out_fold=0,
        seed=56,
        model_config={"metadata_dim": 31, "embedding_dim": 64},
        code_revision="0123456789abcdef",
    )

    assert identity == repeated
    assert identity.digest() == repeated.digest()
    validate_ssl_cache_identity(identity, identity)
    validate_ssl_cache_identity(identity, identity.to_manifest())

    changed = dict(identity.to_manifest())
    changed["seed"] = 57
    with pytest.raises(ValueError, match="digest does not match"):
        validate_ssl_cache_identity(identity, changed)

    extra = dict(identity.to_manifest())
    extra["unexpected"] = True
    with pytest.raises(ValueError, match="keys differ"):
        validate_ssl_cache_identity(identity, extra)

    bad_hashes = dict(hashes)
    bad_hashes["training_table_sha256"] = "ABC"
    with pytest.raises(ValueError, match="lowercase 64-character"):
        build_ssl_cache_identity(
            **bad_hashes,
            selected_rows=selected,
            profile="shape_plus_periodogram_bls",
            held_out_fold=0,
            seed=56,
            model_config={},
            code_revision="revision",
        )
