from __future__ import annotations

import hashlib
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "models" / "twirl_fm0_2_s56_s64_objective_canary.yaml"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _config() -> dict:
    payload = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    assert isinstance(payload, dict)
    return payload


def test_fm0_2_freeze_authorizes_only_the_gated_canary() -> None:
    config = _config()

    assert config["schema_version"] == "twirl_fm0_2_objective_canary_config_v1"
    assert config["status"] == "frozen_gated_canary_authorized_not_started"
    assert config["authorization"] == {
        "scientific_contract_frozen": True,
        "gpu_training_authorized": True,
        "gpu_submission_requires_all_prestart_gates": True,
        "authorized_variant": "TWIRL-FM0.2.1",
        "authorized_max_stop_after_step": 2000,
        "sealed_test_access_authorized": False,
        "production_model_claim": False,
        "foundation_model_claim": False,
    }
    assert config["purpose"]["fm0_1_failure_being_addressed"] == {
        "tcn_effective_rank": 6.200319324,
        "conformer_effective_rank": 8.5792017,
        "required_effective_rank": 26,
    }

    for binding in config["evidence_bindings"].values():
        path = ROOT / binding["relative_path"]
        assert path.is_file()
        assert _sha256(path) == binding["sha256"]


def test_fm0_2_1_changes_only_the_embedding_objective() -> None:
    config = _config()
    fixed = config["held_fixed"]
    objective = config["objective"]

    assert fixed["architecture_geometry"] == "TWIRL-FM0.1.1_TCN"
    assert fixed["parameter_count"] == 8_825_602
    assert fixed["flux_views"] == ["adp_1x1", "adp_3x3"]
    assert fixed["window_length_cadences"] == 2048
    assert fixed["effective_batch_unique_windows"] == 64
    assert fixed["dedicated_ssl_projector"] is False
    assert objective["name"] == "masked_reconstruction_plus_same_window_vicreg"
    assert objective["reconstruction"]["optimized_mask_views"] == ["first"]
    assert objective["vicreg"]["acts_on"] == "z_window"
    assert objective["vicreg"]["begins_at_optimizer_step"] == 1
    assert objective["vicreg"]["hidden_weight_sweep_forbidden"] is True
    assert config["variants"]["TWIRL-FM0.2.1"]["status"] == (
        "initial_objective_canary"
    )
    assert config["variants"]["TWIRL-FM0.2.2"]["status"] == (
        "blocked_until_0_2_1_canary_passes"
    )


def test_fm0_2_vicreg_scale_and_canary_gate_are_predeclared() -> None:
    config = _config()
    vicreg = config["objective"]["vicreg"]
    gates = config["step_2000_go_no_go_gates"]

    collapsed_variance_component = 25.0 * 0.99
    weighted = vicreg["total_weight"] * collapsed_variance_component
    assert abs(weighted - 0.00495) < 1.0e-12
    assert config["canary"]["primary_go_no_go_step"] == 2000
    assert config["canary"]["extension_beyond_step_2000_authorized"] is False
    assert gates["z_window_effective_rank_minimum"] == 26
    assert gates["z_window_constant_dimensions_maximum"] == 0
    assert gates["development_masked_huber_maximum"] == 0.004864579067
    assert gates["sealed_test_access_count"] == 0
    assert config["prestart_gate_contract"]["fail_closed"] is True
    assert config["prestart_gate_contract"][
        "required_before_first_real_data_h200_submission"
    ]
    assert config["canary"]["immutable_milestone_steps"] == [
        0,
        64,
        500,
        1000,
        2000,
    ]
