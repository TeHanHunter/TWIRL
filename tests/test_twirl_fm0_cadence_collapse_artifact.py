from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from twirl.models.fm0 import cadence_collapse_artifact as collapse_artifact
from twirl.models.fm0.cadence_collapse_artifact import (
    CADENCE_OBJECTIVE_IDENTITY,
    CHECKPOINT_SCHEMA_VERSION,
    PANEL_SCHEMA_VERSION,
    RESULT_SCHEMA_VERSION,
    evaluate_cadence_collapse_artifacts,
    sha256_file,
)


def _write_panel(path: Path, *, batch_size: int = 32) -> list[str]:
    rng = np.random.default_rng(560067)
    length = 128
    n_views = 2
    identities = [f"sample-{index:04d}" for index in range(batch_size)]
    delta = np.ones((batch_size, length), dtype=np.float32)
    delta[:, 0] = 0.0
    boundary = np.zeros((batch_size, length), dtype=bool)
    boundary[:, 0] = True
    flux_valid = np.ones((batch_size, n_views, length), dtype=bool)
    # A time-valid cadence with no naturally available model view must not
    # enter the collapse geometry.
    flux_valid[:, :, 10] = False
    np.savez(
        path,
        schema_version=np.asarray(PANEL_SCHEMA_VERSION),
        variant=np.asarray("TWIRL-FM0.3.1"),
        native_cadence_seconds=np.asarray(200, dtype=np.int64),
        cadence_axis_operation=np.asarray("none"),
        sample_identity=np.asarray(identities),
        flux=rng.normal(0.0, 0.01, (batch_size, n_views, length)).astype(np.float32),
        flux_valid=flux_valid,
        flux_error=np.full((batch_size, 2, length), 0.01, dtype=np.float32),
        error_valid=np.ones((batch_size, 2, length), dtype=bool),
        local_time_cadences=np.broadcast_to(
            np.arange(length, dtype=np.float32), (batch_size, length)
        ).copy(),
        delta_time_cadences=delta,
        time_valid=np.ones((batch_size, length), dtype=bool),
        segment_boundary=boundary,
        view_present=np.ones((batch_size, n_views), dtype=bool),
        reconstruction_mask=np.zeros((batch_size, n_views, length), dtype=bool),
    )
    return identities


def _checkpoint_payload(*, torch: Any, step: int) -> dict[str, Any]:
    from twirl.models.fm0.model import TWIRLFM0, FM0ModelConfig

    model_config = FM0ModelConfig(
        architecture="tcn",
        n_flux_views=2,
        window_length=128,
        d_model=32,
        embedding_dim=32,
        stem_kernel=3,
        dropout=0.0,
        tcn_blocks=1,
        tcn_kernel=3,
        tcn_dilation_cycle=(1,),
        conformer_heads=4,
        conformer_conv_kernel=3,
        patch_stride=1,
        minimum_parameters=1,
        maximum_parameters=1_000_000,
    )
    torch.manual_seed(560067)
    model = TWIRLFM0(model_config)
    run_contract = {
        "schema_version": "test",
        "variant": "TWIRL-FM0.3.1",
        "architecture": "tcn",
        "objective": CADENCE_OBJECTIVE_IDENTITY,
        "use_vicreg": True,
        "reconstruct_second_view": False,
        "mask_target_fraction": 0.15,
        "mask_span_range": [1, 4],
        "seed": 560067,
        "immutable_milestone_steps": [0, 2_000],
    }
    dataset_contract = {
        "kind": "synthetic",
        "config": {
            "variant": "TWIRL-FM0.3.1",
            "window_length": 128,
            "mask_target_fraction": 0.15,
            "mask_span_range": (1, 4),
        },
    }
    return {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "model_state": model.state_dict(),
        "optimizer_state": {},
        "scheduler_state": {},
        "scaler_state": None,
        "rng_state": {},
        "progress": {"global_step": step, "epoch": 0},
        "sampler_state": {},
        "model_config": asdict(model_config),
        "optimization_config": {},
        "dataset_contract": dataset_contract,
        "synthetic_dataset_config": dataset_contract["config"],
        "run_contract": run_contract,
        "objective_state": {
            "schema_version": "twirl_fm0_objective_state_v2",
            "identity": CADENCE_OBJECTIVE_IDENTITY,
            "use_vicreg": True,
            "reconstruct_second_view": False,
        },
        "loss_history": [],
    }


def _write_checkpoints(tmp_path: Path) -> tuple[Path, str, Path, str]:
    torch = pytest.importorskip("torch")
    step0 = tmp_path / "checkpoint_step_00000000.pt"
    step2000 = tmp_path / "checkpoint_step_00002000.pt"
    torch.save(_checkpoint_payload(torch=torch, step=0), step0)
    torch.save(_checkpoint_payload(torch=torch, step=2_000), step2000)
    return step0, sha256_file(step0), step2000, sha256_file(step2000)


@pytest.fixture
def tiny_canonical_model(monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep artifact mechanics tests small while preserving exact-config logic."""

    torch = pytest.importorskip("torch")
    from twirl.models.fm0 import model as fm0_model

    payload = _checkpoint_payload(torch=torch, step=0)
    config = fm0_model.FM0ModelConfig(**payload["model_config"])

    def build(variant: str, **_kwargs: object) -> Any:
        if variant != "TWIRL-FM0.3.1":
            raise ValueError("test fixture supports only the TCN FM0.3 arm")
        return fm0_model.TWIRLFM0(config)

    monkeypatch.setattr(fm0_model, "build_fm0_model", build)


def _run(
    panel: Path,
    panel_hash: str,
    step0: Path,
    step0_hash: str,
    step2000: Path,
    step2000_hash: str,
) -> dict[str, Any]:
    return evaluate_cadence_collapse_artifacts(
        panel_path=panel,
        panel_sha256=panel_hash,
        step0_checkpoint_path=step0,
        step0_checkpoint_sha256=step0_hash,
        step2000_checkpoint_path=step2000,
        step2000_checkpoint_sha256=step2000_hash,
    )


def test_artifact_runner_derives_tokens_for_non_authoritative_preflight(
    tmp_path: Path,
    tiny_canonical_model: None,
) -> None:
    panel = tmp_path / "panel.npz"
    identities = _write_panel(panel)
    panel_hash = sha256_file(panel)
    step0, step0_hash, step2000, step2000_hash = _write_checkpoints(tmp_path)

    result = _run(
        panel,
        panel_hash,
        step0,
        step0_hash,
        step2000,
        step2000_hash,
    )

    assert result["schema_version"] == RESULT_SCHEMA_VERSION
    assert result["status"] == "non_authoritative_preflight"
    assert result["authoritative_artifact_gate"] is False
    assert "passed" not in result
    assert result["artifact_bindings"]["panel_sha256"] == panel_hash
    assert result["artifact_bindings"]["step0_checkpoint_sha256"] == step0_hash
    assert result["artifact_bindings"]["step2000_checkpoint_sha256"] == step2000_hash
    expected_order_hash = hashlib.sha256(
        json.dumps(tuple(identities), separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    assert result["artifact_bindings"]["sample_order_sha256"] == (expected_order_hash)
    assert result["checkpoint_contract"]["exact_steps"] == [0, 2_000]
    assert result["checkpoint_contract"]["window_length"] == 128
    assert result["checkpoint_contract"]["patch_stride"] == 1
    assert result["inference"]["step0_token_shape"] == [32, 128, 32]
    assert result["inference"]["step2000_token_shape"] == [32, 128, 32]
    assert result["inference"]["identical_token_valid"] is True
    assert result["inference"]["naturally_available_token_count"] == 32 * 127
    assert result["inference"]["one_native_cadence_per_token"] is True
    assert result["inference"]["temporal_averaging_or_pooling"] is False
    assert result["gate"]["step0"]["eligible_contiguous_edges"] == 32 * 125
    assert result["criteria_satisfied"] is result["gate"]["criteria_satisfied"]
    assert "passed" not in result["gate"]
    assert result["gate"]["authoritative_artifact_gate"] is False
    assert result["formal_gate_prerequisites"] == [
        "separately_frozen_panel_receipt",
        "validated_real_run_and_input_release_contract",
        "immutable_step0_and_step2000_milestone_checkpoints",
    ]


def test_panel_substitution_is_rejected_by_actual_byte_hash(tmp_path: Path) -> None:
    panel = tmp_path / "panel.npz"
    _write_panel(panel)
    original_hash = sha256_file(panel)
    panel.write_bytes(panel.read_bytes() + b"substitution")

    with pytest.raises(ValueError, match="panel SHA-256 differs"):
        _run(
            panel,
            original_hash,
            panel,
            "a" * 64,
            panel,
            "b" * 64,
        )


def test_checkpoint_hash_substitution_and_global_step_are_rejected(
    tmp_path: Path,
    tiny_canonical_model: None,
) -> None:
    torch = pytest.importorskip("torch")
    panel = tmp_path / "panel.npz"
    _write_panel(panel)
    panel_hash = sha256_file(panel)
    step0, step0_hash, step2000, step2000_hash = _write_checkpoints(tmp_path)
    original_step2000_hash = step2000_hash
    step2000.write_bytes(step2000.read_bytes() + b"substitution")
    with pytest.raises(ValueError, match="step-2000 checkpoint SHA-256 differs"):
        _run(
            panel,
            panel_hash,
            step0,
            step0_hash,
            step2000,
            original_step2000_hash,
        )

    torch.save(_checkpoint_payload(torch=torch, step=1), step0)
    with pytest.raises(ValueError, match="must be exact step 0"):
        _run(
            panel,
            panel_hash,
            step0,
            sha256_file(step0),
            step2000,
            sha256_file(step2000),
        )


@pytest.mark.parametrize("drift", ["schema", "run", "model", "dataset", "objective"])
def test_checkpoint_schema_and_exact_contract_drift_are_rejected(
    tmp_path: Path,
    drift: str,
    tiny_canonical_model: None,
) -> None:
    torch = pytest.importorskip("torch")
    panel = tmp_path / "panel.npz"
    _write_panel(panel)
    panel_hash = sha256_file(panel)
    step0, step0_hash, step2000, _ = _write_checkpoints(tmp_path)
    payload = _checkpoint_payload(torch=torch, step=2_000)
    if drift == "schema":
        payload["schema_version"] = "wrong"
    elif drift == "run":
        payload["run_contract"]["seed"] = 560068
    elif drift == "model":
        payload["model_config"]["dropout"] = 0.1
    elif drift == "dataset":
        payload["dataset_contract"]["config"]["mask_span_range"] = (1, 3)
    else:
        payload["objective_state"]["identity"] = "window_pooled_vicreg"
    torch.save(payload, step2000)

    match = {
        "schema": "schema version differs",
        "run": "run contract differ",
        "model": "canonical build",
        "dataset": "dataset contract geometry differs",
        "objective": "objective contract differs",
    }[drift]
    with pytest.raises((TypeError, ValueError), match=match):
        _run(
            panel,
            panel_hash,
            step0,
            step0_hash,
            step2000,
            sha256_file(step2000),
        )


def test_noncanonical_fm03_model_config_is_rejected(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    panel = tmp_path / "panel.npz"
    _write_panel(panel)
    panel_hash = sha256_file(panel)
    step0, step0_hash, step2000, step2000_hash = _write_checkpoints(tmp_path)

    with pytest.raises(ValueError, match="canonical build"):
        _run(
            panel,
            panel_hash,
            step0,
            step0_hash,
            step2000,
            step2000_hash,
        )


def test_composite_dataset_contract_exposes_native_cadence_geometry() -> None:
    contract = {
        "kind": "fm0_3_composite_release",
        "variant": "TWIRL-FM0.3.1",
        "window_length": 128,
        "mask_target_fraction": 0.15,
        "mask_span_range": [1, 4],
    }
    assert collapse_artifact._dataset_geometry(contract) is contract


def test_duplicate_panel_identity_and_nonzero_mask_fail_closed(tmp_path: Path) -> None:
    panel = tmp_path / "panel.npz"
    _write_panel(panel, batch_size=2)
    with np.load(panel, allow_pickle=False) as archive:
        payload = {name: np.asarray(archive[name]).copy() for name in archive.files}
    payload["sample_identity"][:] = "duplicate"
    np.savez(panel, **payload)
    step0 = tmp_path / "step0.pt"
    step2000 = tmp_path / "step2000.pt"
    step0.write_bytes(b"step0")
    step2000.write_bytes(b"step2000")
    with pytest.raises(ValueError, match="sample identities must be unique"):
        evaluate_cadence_collapse_artifacts(
            panel_path=panel,
            panel_sha256=sha256_file(panel),
            step0_checkpoint_path=step0,
            step0_checkpoint_sha256=sha256_file(step0),
            step2000_checkpoint_path=step2000,
            step2000_checkpoint_sha256=sha256_file(step2000),
        )
