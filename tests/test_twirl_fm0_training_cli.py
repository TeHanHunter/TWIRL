from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
import yaml

from twirl.models.fm0.training import (
    CADENCE_VICREG_OBJECTIVE_IDENTITY,
    FM0_3_MASK_SPAN_RANGE,
    FM0_3_MASK_TARGET_FRACTION,
    FM0_3_WINDOW_LENGTH,
)

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "stage5_validation" / "train_twirl_fm0.py"
SPEC = importlib.util.spec_from_file_location("train_twirl_fm0_cli_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
CLI = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLI)


def _fm03_config() -> dict[str, object]:
    return {
        "schema_version": CLI.FM0_3_CONFIG_SCHEMA_VERSION,
        "objective": {
            "name": CADENCE_VICREG_OBJECTIVE_IDENTITY,
            "reconstruction": {"optimized_mask_views": ["first"]},
        },
    }


def test_cli_accepts_only_the_explicit_fm03_training_schema(tmp_path: Path) -> None:
    path = tmp_path / "fm03.yaml"
    path.write_text(yaml.safe_dump(_fm03_config()), encoding="utf-8")
    assert CLI._load_config(path)["schema_version"] == (
        "twirl_fm0_3_cadence_objective_config_v1"
    )

    with pytest.raises(ValueError, match="config schema mismatch"):
        CLI._load_config(
            ROOT / "configs" / "models" / "twirl_fm0_3_stride1_mechanics_v1.yaml"
        )


def test_cli_resolves_fm03_cadence_objective_without_window_vicreg() -> None:
    settings = CLI._objective_settings(
        _fm03_config(),
        {"objective": CADENCE_VICREG_OBJECTIVE_IDENTITY},
    )
    assert settings == (
        CADENCE_VICREG_OBJECTIVE_IDENTITY,
        True,
        False,
        CADENCE_VICREG_OBJECTIVE_IDENTITY,
    )


@pytest.mark.parametrize(
    ("variant_objective", "optimized_masks", "message"),
    [
        (
            "masked_reconstruction_plus_same_window_vicreg",
            ["first"],
            "cadence-level objective identity",
        ),
        (
            CADENCE_VICREG_OBJECTIVE_IDENTITY,
            ["first", "second"],
            "optimizes only the first mask",
        ),
    ],
)
def test_cli_rejects_fm03_objective_contract_drift(
    variant_objective: str,
    optimized_masks: list[str],
    message: str,
) -> None:
    config = _fm03_config()
    config["objective"]["name"] = variant_objective
    config["objective"]["reconstruction"]["optimized_mask_views"] = optimized_masks
    with pytest.raises(ValueError, match=message):
        CLI._objective_settings(config, {"objective": variant_objective})


def test_cli_freezes_fm03_native_cadence_dataset_geometry() -> None:
    assert CLI._dataset_geometry(CLI.FM0_3_CONFIG_SCHEMA_VERSION) == {
        "window_length": FM0_3_WINDOW_LENGTH,
        "mask_target_fraction": FM0_3_MASK_TARGET_FRACTION,
        "mask_span_range": FM0_3_MASK_SPAN_RANGE,
    }
    assert CLI._dataset_geometry("twirl_fm0_2_objective_canary_config_v1") == {}


def test_cli_keeps_fm03_invocation_stop_out_of_invariant_contract() -> None:
    first: dict[str, object] = {}
    later: dict[str, object] = {}
    for contract, target in ((first, 1), (later, 500)):
        CLI._bind_training_stop_contract(
            contract,
            schema_version=CLI.FM0_3_CONFIG_SCHEMA_VERSION,
            optimizer_horizon=20_000,
            invocation_target_step=target,
        )
    assert first == later == {
        "training_horizon_step": 20_000,
        "stop_after_step_is_execution_state_not_scientific_contract": True,
    }
    assert "target_step" not in first

    legacy: dict[str, object] = {}
    CLI._bind_training_stop_contract(
        legacy,
        schema_version="twirl_fm0_1_poc_config_v1",
        optimizer_horizon=20_000,
        invocation_target_step=7,
    )
    assert legacy == {"target_step": 7}


def test_cli_preserves_fm02_invariant_canary_stop_contract() -> None:
    contract: dict[str, object] = {}
    CLI._bind_training_stop_contract(
        contract,
        schema_version=CLI.FM0_2_CONFIG_SCHEMA_VERSION,
        optimizer_horizon=20_000,
        invocation_target_step=500,
        canary_payload={
            "immutable_milestone_steps": [0, 64, 500, 1000, 2000],
            "authorized_stop_after_steps": [64, 500, 1000, 2000],
            "fp32_synthetic_smoke_steps": 8,
        },
    )
    assert contract == {
        "training_horizon_step": 20_000,
        "immutable_milestone_steps": [0, 64, 500, 1000, 2000],
        "authorized_stop_after_steps": [64, 500, 1000, 2000],
        "synthetic_smoke_step": 8,
    }
    assert "target_step" not in contract
