from __future__ import annotations

import numpy as np
from pathlib import Path
import pytest
import subprocess
import sys

from twirl.models.fm0.dataset import (
    evaluation_window_starts,
    prepare_model_window,
    synchronized_temporal_mask,
)
from twirl.models.fm0.validation import (
    RUN_CONTRACT_SCHEMA_VERSION,
    RUN_SUMMARY_SCHEMA_VERSION,
    sha256_file,
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
