from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import yaml

from twirl.models.fm0.training import CADENCE_VICREG_OBJECTIVE_IDENTITY
from twirl.models.fm0.validation import validate_frozen_authorities

ROOT = Path(__file__).resolve().parents[1]
CONFIG = (
    ROOT
    / "configs/models/twirl_fm0_3_s56_s64_s66_s77_matched_canary_v1.yaml"
)
DESIGN = ROOT / "doc/fm0_3_matched_canary_design.md"
FREEZE = (
    ROOT
    / "reports/stage5_validation/twirl_fm0_3_matched_canary_design_freeze_v1/freeze.json"
)
WRAPPER = ROOT / "scripts/orcd/slurm_twirl_fm0_3_matched_canary_h200.sbatch"
CONTROLLER = ROOT / "scripts/orcd/run_twirl_fm0_3_matched_canary_orcd.sh"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_matched_canary_config_freezes_only_the_intended_architecture_change() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    assert config["schema_version"] == "twirl_fm0_3_cadence_objective_config_v1"
    assert config["status"] == "frozen_matched_canary_authorized_not_started"
    assert config["authorization"] == {
        "scientific_contract_frozen": True,
        "gpu_training_authorized": True,
        "gpu_submission_requires_prestart_smoke": True,
        "payload_screened_evaluation_freeze_required": True,
        "authorized_variants": ["TWIRL-FM0.3.1", "TWIRL-FM0.3.2"],
        "authorized_max_stop_after_step": 2_000,
        "sealed_test_access_authorized": False,
        "production_model_claim": False,
        "foundation_model_claim": False,
    }
    assert config["composite_input"] == {
        "receipt_sha256": (
            "cc5cc3bce4c24e74bef1fbf084f407855233de7893183e6bc3486e284a2f44d9"
        ),
        "source_bindings_sha256": (
            "8cbcac99409ab89fe2dee0c36687f50122e78376206495073698489ca5424f2b"
        ),
        "role_index_sha256": (
            "abe9c616523f2486bf1b7be69dfcfda6193d534b40b3f4afc49d8ebf3e40ce5a"
        ),
    }

    cadence = config["cadence_contract"]
    assert cadence["nominal_cadence_seconds"] == 200
    assert cadence["context_length_cadences"] == 128
    assert cadence["patch_stride"] == 1
    for field in (
        "cadence_averaging",
        "cadence_merging",
        "temporal_downsampling",
        "temporal_pooling",
        "pooled_window_objective",
    ):
        assert cadence[field] is False
    assert cadence["optimized_representation"] == "h_cadence"
    assert config["data_contract"]["flux_views"] == ["adp_1x1", "adp_3x3"]
    assert config["data_contract"]["excluded_sectors"] == [65]
    assert config["data_contract"]["sealed_rows_allowed"] == 0

    variants = config["variants"]
    assert variants["TWIRL-FM0.3.1"]["architecture"] == "tcn"
    assert variants["TWIRL-FM0.3.2"]["architecture"] == "conformer"
    assert {
        payload["reconstruction_input"] for payload in variants.values()
    } == {"contextual_h_cadence_without_stem_skip"}
    assert {
        payload["objective"] for payload in variants.values()
    } == {CADENCE_VICREG_OBJECTIVE_IDENTITY}
    assert config["objective"]["name"] == CADENCE_VICREG_OBJECTIVE_IDENTITY
    assert config["objective"]["reconstruction"]["optimized_mask_views"] == [
        "first"
    ]
    assert config["objective"]["vicreg"]["total_weight"] == 0.0002
    evaluation = config["evaluation_freeze_contract"]
    assert evaluation["required_before_real_checkpoint"] is True
    assert evaluation["temporal_panel_receipt_sha256"] == (
        "78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
    )
    assert evaluation["n_identity_crops"] == 1_440
    assert evaluation["crop_length_cadences"] == 128
    assert evaluation["minimum_joint_time_and_both_adp_valid_cadences"] == 103
    assert evaluation["complete_joint_valid_support_indices_inclusive"] == [
        60,
        68,
    ]
    assert evaluation["inference"]["readout"] == "h_cadence_token_64_only"
    assert evaluation["inference"]["token_or_window_pooling"] is False
    assert evaluation["probe"] == {
        "family": "full_batch_adam_linear_logit",
        "train_only_dimensionwise_standardization": True,
        "epochs": 400,
        "learning_rate": 0.02,
        "l2_weight": 0.001,
        "seed": 560306,
        "hyperparameter_tuning": False,
    }
    assert evaluation["bls_or_periodic_search"] is False

    optimization = config["optimization"]
    assert optimization["seeds"] == [560067]
    assert optimization["effective_batch_windows"] == 64
    assert optimization["learning_rate"] == 0.0003
    assert optimization["weight_decay"] == 0.01
    assert optimization["warmup_steps"] == 1_000
    assert optimization["max_optimizer_steps"] == 20_000
    assert config["canary"] == {
        "fp32_synthetic_smoke_steps": 8,
        "bf16_real_throughput_and_restart_step": 64,
        "immutable_milestone_steps": [0, 64, 2_000],
        "authorized_stop_after_steps": [64, 2_000],
        "one_h200_per_variant": True,
    }


def test_design_freeze_binds_the_new_doc_config_and_composite() -> None:
    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    assert freeze["freeze_schema_version"] == (
        "twirl_fm0_3_design_freeze_receipt_v1"
    )
    assert freeze["implementation_status"] == (
        "authorized_matched_canary_not_started"
    )
    assert freeze["authorization"] == {
        "gpu_submission_requires_prestart_smoke": True,
        "matched_orcd_canary": True,
        "payload_screened_evaluation_freeze_required": True,
        "sealed_test_access": False,
    }
    assert freeze["design_document"]["sha256"] == _sha256(DESIGN)
    assert freeze["scientific_config"]["sha256"] == _sha256(CONFIG)
    assert freeze["evaluation_input"] == {
        "temporal_panel_receipt_sha256": (
            "78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
        )
    }
    sidecar = FREEZE.with_name("freeze.json.sha256").read_text().strip().split()
    assert sidecar == [_sha256(FREEZE), "freeze.json"]

    authorities = validate_frozen_authorities(
        design_path=DESIGN,
        config_path=CONFIG,
        freeze_receipt_path=FREEZE,
    )
    assert authorities["design_sha256"] == _sha256(DESIGN)
    assert authorities["config_sha256"] == _sha256(CONFIG)
    assert authorities["freeze_receipt_sha256"] == _sha256(FREEZE)


def test_matched_canary_scripts_are_one_h200_guarded_and_shell_valid() -> None:
    for script in (WRAPPER, CONTROLLER):
        subprocess.run(["bash", "-n", str(script)], check=True)

    wrapper = WRAPPER.read_text(encoding="utf-8")
    assert "#SBATCH --gres=gpu:h200:1" in wrapper
    assert "#SBATCH --array" not in wrapper
    assert 'case "${VARIANT}"' in wrapper
    assert 'case "${STOP_AFTER}"' in wrapper
    assert "TWIRL-FM0.3.1" in wrapper and "TWIRL-FM0.3.2" in wrapper
    assert 'readonly PRECISION="fp32"' in wrapper
    assert 'readonly PRECISION="bf16"' in wrapper
    assert "--micro-batch-windows 8" in wrapper
    assert "--design \"${DESIGN}\"" in wrapper
    assert "--input-release \"${COMPOSITE_ROOT}\"" in wrapper
    assert "--evaluation-plan \"${EVALUATION_PLAN}\"" in wrapper
    assert "validate_matched_canary_payload_plan" in wrapper
    assert "TWIRL_FM0_EVALUATION_PLAN_RECEIPT_SHA256" in wrapper
    assert "EVALUATION_PLAN_PIN" in wrapper
    assert "evaluation-plan receipt differs from campaign pin" in wrapper
    assert "differs from launch authority" in wrapper
    assert "temporal[\"receipt_sha256\"] != panel_receipt_sha" in wrapper
    assert "checkpoint_step_00000064.pt" in wrapper

    controller = CONTROLLER.read_text(encoding="utf-8")
    assert 'RUN="status"' in controller
    assert 'EXECUTE=0' in controller
    assert 'if (( ! EXECUTE ))' in controller
    for option in (
        "ControlMaster=no",
        "BatchMode=yes",
        "PasswordAuthentication=no",
        "KbdInteractiveAuthentication=no",
        "PubkeyAuthentication=no",
        "HostbasedAuthentication=no",
        "GSSAPIAuthentication=no",
        "NumberOfPasswordPrompts=0",
    ):
        assert option in controller


def test_shell_gates_fully_validate_smoke_and_step64_before_advancing() -> None:
    wrapper = WRAPPER.read_text(encoding="utf-8")
    controller = CONTROLLER.read_text(encoding="utf-8")
    for script in (wrapper, controller):
        assert "validate_run_release(root, inspect_checkpoint=True)" in script
        assert "validate_real_run_release(root, inspect_checkpoint=True)" in script
        assert 'validation.get("checkpoint_inspected") is not True' in script
        assert 'validation.get("global_step") != expected_step' in script
        assert '"design_sha256": design_sha' in script
        assert '"config_sha256": config_sha' in script
        assert '"freeze_receipt_sha256": freeze_sha' in script
        assert 'contract.get("authorities") != expected_authorities' in script
        assert "validate_matched_canary_payload_plan" in script

    assert 'validate_prior_release synthetic "${SMOKE_OUT}" 8' in wrapper
    assert 'validate_prior_release real "${OUT}" 64' in wrapper
    assert wrapper.index('validate_prior_release real "${OUT}" 64') < wrapper.index(
        '"${PYTHON}" scripts/stage5_validation/train_twirl_fm0.py'
    )

    assert "validate_prior_release synthetic \\" in controller
    assert 'validate_prior_release real "${OUTPUT}" 64' in controller
    assert "TWIRL_FM0_EVALUATION_PLAN_RECEIPT_SHA256=" in controller
    assert "evaluation_plan.sha256" in controller
    assert "set -o noclobber" in controller
    assert 'pin_evaluation_plan "${EVALUATION_RECEIPT_SHA}"' in controller
    assert "temporal[\"receipt_sha256\"] != panel_receipt_sha" in controller
    assert controller.index(
        'pin_evaluation_plan "${EVALUATION_RECEIPT_SHA}"'
    ) < controller.index("job=\\$(sbatch --parsable")
    assert "${EVALUATION_RECEIPT_SHA}" in controller
    assert r"printf '%s\\t%s\\t%s\\t%s\\t%s\\n'" in controller
    assert controller.index(
        'validate_prior_release real "${OUTPUT}" 64'
    ) < controller.index("job=\\$(sbatch --parsable")
    assert "remote \"test -s '${OUTPUT}/checkpoint_step_00000064.pt'" not in (
        controller
    )


def test_controller_submit_is_a_dry_run_without_execute() -> None:
    result = subprocess.run(
        [
            "bash",
            str(CONTROLLER),
            "--run",
            "submit",
            "--variant",
            "TWIRL-FM0.3.2",
            "--stop-after-step",
            "64",
        ],
        cwd=ROOT,
        env={**os.environ, "TWIRL_EXPECTED_GIT_SHA": "a" * 40},
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.startswith("DRY RUN: sbatch ")
    assert "twirl-fm0-032-s64" in result.stdout
    assert "No job submitted; add --execute" in result.stdout
