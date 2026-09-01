from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from twirl.models.fm0.dataset import (
    SyntheticFM0Config,
    SyntheticFM0Dataset,
    evaluation_window_starts,
    prepare_model_window,
    synchronized_temporal_mask,
    variant_view_indices,
)
from twirl.models.fm0.training import _synthetic_dataset_config_contract
from twirl.models.fm0.validation import (
    RUN_CONTRACT_SCHEMA_VERSION,
    RUN_SUMMARY_SCHEMA_VERSION,
    validate_run_release,
    write_json_with_sha256,
    write_sha256_sidecar,
)


def _window(length: int = 20) -> dict[str, np.ndarray]:
    return {
        "flux": np.arange(length * 6, dtype=np.float32).reshape(length, 6),
        "flux_valid": np.ones((length, 6), dtype=bool),
        "flux_error": np.full((length, 2), 0.01, dtype=np.float32),
        "error_valid": np.ones((length, 2), dtype=bool),
        "local_time_cadences": np.arange(length, dtype=np.float32),
        "delta_time_cadences": np.r_[0.0, np.ones(length - 1)].astype(np.float32),
        "time_valid": np.ones(length, dtype=bool),
        "segment_boundary": np.r_[True, np.zeros(length - 1, dtype=bool)],
        "view_present": np.ones(6, dtype=bool),
    }


def test_evaluation_windows_keep_strict_stride_and_pad_last() -> None:
    assert evaluation_window_starts(3000) == (0, 1024, 2048)
    assert evaluation_window_starts(2048) == (0, 1024)
    assert evaluation_window_starts(0) == ()


def test_synchronized_mask_is_deterministic_and_never_targets_invalid() -> None:
    valid = np.ones((3, 100), dtype=bool)
    valid[1, 10:20] = False
    time_valid = np.ones(100, dtype=bool)
    time_valid[30:35] = False
    first_time, first = synchronized_temporal_mask(valid, time_valid, seed=7)
    second_time, second = synchronized_temporal_mask(valid, time_valid, seed=7)
    assert np.array_equal(first_time, second_time)
    assert np.array_equal(first, second)
    assert np.array_equal(first, valid & first_time[None, :])
    assert not np.any(first[:, 30:35])


def test_fm0_3_variants_keep_the_two_adp_views() -> None:
    assert variant_view_indices("TWIRL-FM0.3.1") == (2, 3)
    assert variant_view_indices("TWIRL-FM0.3.2") == (2, 3)


def test_short_synthetic_event_and_scaled_mask_are_nonempty_and_deterministic() -> None:
    config = SyntheticFM0Config(
        variant="TWIRL-FM0.3.1",
        window_length=128,
        windows_per_epoch=2,
        noise_scale=0.0,
        event_depth=0.1,
        mask_target_fraction=0.15,
        mask_span_range=(1, 4),
    )
    dataset = SyntheticFM0Dataset(config)
    first = dataset.sample(0, mask_view=0)
    repeated = dataset.sample(0, mask_view=0)
    second_mask = dataset.sample(0, mask_view=1)

    for name in first:
        assert np.array_equal(first[name], repeated[name]), name
    assert np.any(first["flux"][first["flux_valid"]] < -0.05)
    assert np.array_equal(first["flux"], second_mask["flux"])
    assert not np.array_equal(first["temporal_mask"], second_mask["temporal_mask"])

    eligible = first["time_valid"] & np.any(first["flux_valid"], axis=0)
    target = int(np.ceil(0.15 * np.count_nonzero(eligible)))
    masked = int(np.count_nonzero(first["temporal_mask"]))
    assert target <= masked <= target + 3
    assert np.array_equal(
        first["reconstruction_mask"],
        first["flux_valid"] & first["temporal_mask"][None, :],
    )


def test_legacy_synthetic_defaults_keep_checkpoint_keys_and_event_sequence() -> None:
    config = SyntheticFM0Config(
        variant="TWIRL-FM0.1.1",
        windows_per_epoch=2,
        noise_scale=0.0,
        event_depth=0.1,
    )
    assert _synthetic_dataset_config_contract(config) == {
        "variant": "TWIRL-FM0.1.1",
        "seed": 560067,
        "n_sources": 64,
        "visits_per_source": 2,
        "windows_per_epoch": 2,
        "window_length": 2048,
        "noise_scale": 0.0,
        "event_depth": 0.1,
    }

    sample = SyntheticFM0Dataset(config)[0]
    for view, expected_depth in enumerate((-0.104, -0.106)):
        event_indices = np.flatnonzero(sample["flux"][view] < -0.05)
        assert np.array_equal(event_indices, np.arange(1262, 1264))
        assert sample["flux"][view, 1262] == pytest.approx(expected_depth)


@pytest.mark.parametrize(
    "overrides",
    (
        {"mask_target_fraction": -0.01},
        {"mask_target_fraction": 1.01},
        {"mask_span_range": (0, 4)},
        {"mask_span_range": (4, 1)},
    ),
)
def test_synthetic_config_rejects_invalid_mask_policy(
    overrides: dict[str, object],
) -> None:
    with pytest.raises(ValueError, match="mask_"):
        SyntheticFM0Config(variant="TWIRL-FM0.3.1", **overrides)


def test_prepare_window_selects_frozen_views_and_pads_neutrally() -> None:
    sample = prepare_model_window(
        _window(), variant="TWIRL-FM0.1.1", mask_seed=5, window_length=32
    )
    assert sample["flux"].shape == (2, 32)
    assert np.array_equal(sample["flux"][:, :20], _window()["flux"][:, (2, 3)].T)
    assert not np.any(sample["flux_valid"][:, 20:])
    assert not np.any(sample["time_valid"][20:])
    assert np.all(sample["flux"][:, 20:] == 0)
    assert np.array_equal(
        sample["reconstruction_mask"],
        sample["flux_valid"] & sample["temporal_mask"][None, :],
    )


def test_prepare_window_fails_when_a_selected_variant_view_is_missing() -> None:
    window = _window()
    window["view_present"][3] = False
    window["flux_valid"][:, 3] = False
    with pytest.raises(ValueError, match="adp_3x3"):
        prepare_model_window(
            window, variant="TWIRL-FM0.1.1", mask_seed=5, window_length=32
        )

    window = _window()
    window["flux_valid"][:, 3] = False
    with pytest.raises(ValueError, match="no valid cadences.*adp_3x3"):
        prepare_model_window(
            window, variant="TWIRL-FM0.1.1", mask_seed=5, window_length=32
        )


def test_synthetic_run_release_validator_binds_three_artifacts(tmp_path: Path) -> None:
    contract_hash = write_json_with_sha256(
        tmp_path / "run_contract.json",
        {"schema_version": RUN_CONTRACT_SCHEMA_VERSION, "synthetic_smoke": True},
    )
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"synthetic checkpoint placeholder")
    checkpoint_hash = write_sha256_sidecar(checkpoint)
    write_json_with_sha256(
        tmp_path / "summary.json",
        {
            "schema_version": RUN_SUMMARY_SCHEMA_VERSION,
            "passed": True,
            "synthetic_only": True,
            "global_step": 1,
            "variant": "TWIRL-FM0.1.1",
            "architecture": "tcn",
            "run_contract_sha256": contract_hash,
            "checkpoint_sha256": checkpoint_hash,
        },
    )
    # This unit deliberately uses a byte placeholder so it can exercise the
    # sidecar contract in environments without the optional Torch dependency.
    result = validate_run_release(tmp_path, inspect_checkpoint=False)
    assert result["passed"] is True
    assert result["checkpoint_inspected"] is False
    checkpoint.write_bytes(b"tampered")
    try:
        validate_run_release(tmp_path, inspect_checkpoint=False)
    except ValueError as exc:
        assert "SHA256 mismatch" in str(exc)
    else:  # pragma: no cover - defensive
        raise AssertionError("tampered checkpoint passed")


def test_training_cli_exposes_wrapper_compatible_arguments() -> None:
    root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "stage5_validation" / "train_twirl_fm0.py"),
            "--help",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    for argument in (
        "--synthetic-smoke",
        "--max-steps",
        "--input-release-receipt",
        "--expected-git-sha",
    ):
        assert argument in completed.stdout
