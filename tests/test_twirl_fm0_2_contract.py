from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "models" / "twirl_fm0_2_s56_s64_objective_canary.yaml"
CONTROLLER = ROOT / "scripts" / "orcd" / "run_twirl_fm0_2_canary_orcd.sh"
PROJECT_ROOT = "/orcd/data/mki_aryeh/001/twirl"
TRAINING_SHA = "ddf442aafb8f62966e549e2287abad3474dd556a"
POST_SHA = "80d678c4621a22258c5d0e8f0be6f7e08ffa5bf93841c315abf42e9bb006b110"
REFERENCE_SHA = "616ef9a100b7b7a3a1923f81ac19a272cab8b6f4657a4a58c2815688ee6d1191"
STEP0_REPORT = (
    ROOT / "reports" / "stage5_validation" / "twirl_fm0_2_s56_s64_step0_evaluation_v1"
)
STEP0_EVALUATION_SHA = (
    "bc75b96271487851162f542682501f30fba1bcaa40b425333ad1767c32da1fe6"
)
STEP0_RECEIPT_SHA = "5e489c8abf0ef7c70c815889701a8dfc1a7ad31eed7b5a4366158ec68d61e29c"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _config() -> dict:
    payload = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    assert isinstance(payload, dict)
    return payload


def _initialization_environment() -> dict[str, str]:
    artifact_root = (
        f"{PROJECT_ROOT}/reports/stage5_validation/"
        f"twirl_fm0_2_s56_s64_objective_canary/{TRAINING_SHA[:12]}"
    )
    run_output = f"{artifact_root}/model_runs/fm0_2_1_canary/seed560067"
    return {
        **os.environ,
        "TWIRL_EXPECTED_GIT_SHA": "1" * 40,
        "TWIRL_FM0_ARTIFACT_GIT_SHA": TRAINING_SHA,
        "TWIRL_FM0_RUN_OUTPUT": run_output,
        "TWIRL_FM0_INITIALIZATION_CHECKPOINT": (
            f"{run_output}/checkpoint_step_00000000.pt"
        ),
        "TWIRL_FM0_STEP2000_POST_VALIDATION": (
            f"{artifact_root}/validations/"
            "seed560067-real_canary-step00002000/"
            "post_validation.receipt.json"
        ),
        "TWIRL_FM0_STEP2000_POST_VALIDATION_SHA256": POST_SHA,
        "TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT": (
            f"{artifact_root}/evaluations/step_00002000/"
            "representation_health.receipt.json"
        ),
        "TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT_SHA256": REFERENCE_SHA,
    }


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
    assert config["variants"]["TWIRL-FM0.2.1"]["status"] == ("initial_objective_canary")
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


def test_fm0_2_representation_evaluation_is_cpu_only_and_fail_closed() -> None:
    controller = (
        ROOT / "scripts" / "orcd" / "run_twirl_fm0_2_canary_orcd.sh"
    ).read_text(encoding="utf-8")
    wrapper = (
        ROOT / "scripts" / "orcd" / "slurm_twirl_fm0_2_representation_health_cpu.sbatch"
    ).read_text(encoding="utf-8")

    assert "submit-representation-evaluation)" in controller
    assert "slurm_twirl_fm0_2_representation_health_cpu.sbatch" in controller
    assert "500|1000|2000" in controller
    assert "#SBATCH -c 4" in wrapper
    assert "#SBATCH --mem=16G" in wrapper
    assert "#SBATCH --exclude=node4900" in wrapper
    assert "#SBATCH --gres" not in wrapper
    assert "build_twirl_fm0_observation_sector_authority.py" in wrapper
    assert "evaluate_twirl_fm0_representation_health.py" in wrapper
    assert '--checkpoint "${CHECKPOINT}"' in wrapper
    assert '--observation-sector-authority "${AUTHORITY_CSV}"' in wrapper
    assert 'get("sealed_test_access_count") == 0' in wrapper
    assert '"scientific_go_no_go_applied": False' in wrapper


def test_fm0_2_initialization_evaluation_is_a_separate_pinned_cpu_path() -> None:
    controller = CONTROLLER.read_text(encoding="utf-8")
    wrapper = (
        ROOT
        / "scripts"
        / "orcd"
        / "slurm_twirl_fm0_2_initialization_representation_health_cpu.sbatch"
    ).read_text(encoding="utf-8")

    assert "submit-initialization-representation-evaluation)" in controller
    assert "slurm_twirl_fm0_2_initialization_representation_health_cpu.sbatch" in (
        controller
    )
    assert "checkpoint_step_00000000.pt" in controller
    assert TRAINING_SHA in controller
    assert POST_SHA in controller
    assert REFERENCE_SHA in controller
    assert 'case "${TWIRL_FM0_STOP_AFTER_STEP}" in 500|1000|2000)' in controller

    assert "#SBATCH -c 4" in wrapper
    assert "#SBATCH --mem=16G" in wrapper
    assert "#SBATCH --exclude=node4900" in wrapper
    assert "#SBATCH --gres" not in wrapper
    assert "--device cpu" in wrapper
    assert "TWIRL_fm0_2_${EXPECTED_EVALUATOR_GIT_SHA:0:12}" in wrapper
    assert "83816d07975eebe3825d76dfe7096d22b70376f5" in wrapper
    assert "2ad38165d0d89acb08970e9bdf2a07df54022bf020b6a72e310e4cb4eb3f014e" in wrapper
    assert "573335ee9e2e9f3a7cc4aedf42dec76a584f0a72cdebacbf9f0e805a81d489f8" in wrapper
    assert "92463070381486f0c6053c190d62da8d0c5c0d31be8072e2bdcd677329ac792c" in wrapper
    assert "3c802741744999fe553e71bc2dccfbef6c309a29a28c30dfaac693a798bb3718" in wrapper
    assert "bdd3e8039b312aeb21557662226a8a11b761dd1da031742be42ee9b0d1c6edc5" in wrapper
    assert 'mkdir -- "${OUTPUT_DIR}"' in wrapper
    assert 'mkdir -p "${OUTPUT_DIR}"' not in wrapper
    assert "sys.flags.optimize == 0" in wrapper
    assert "optimized Python would disable fail-closed assertions" in wrapper
    assert 'get("global_step") == 0' in wrapper
    assert 'get("sealed_test_access_count") == 0' in wrapper
    assert '"evaluation_role": "exact_same_seed_initialization"' in wrapper
    assert (
        '"random_encoder_control": "separate_matched_random_control_seed_0"' in wrapper
    )
    assert '"development_event_retention_run": False' in wrapper
    assert '"scientific_go_no_go_applied": False' in wrapper
    assert '"eligible_evidence_steps_after_receipt_validation"' in wrapper
    assert '"performed_by_this_job": False' in wrapper
    assert "${INPUT_RECEIPT_SHA256}" not in wrapper
    assert "sbatch" not in wrapper


def test_fm0_2_initialization_controller_dry_run_is_exact() -> None:
    result = subprocess.run(
        ["bash", str(CONTROLLER), "submit-initialization-representation-evaluation"],
        cwd=ROOT,
        env=_initialization_environment(),
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "slurm_twirl_fm0_2_initialization_representation_health_cpu.sbatch" in (
        result.stdout
    )
    assert "checkpoint_step_00000000.pt" in result.stdout
    assert "step00002000/post_validation.receipt.json" in result.stdout
    assert "step_00002000/representation_health.receipt.json" in result.stdout
    assert "TWIRL_FM0_STOP_AFTER_STEP" not in result.stdout


def test_fm0_2_generic_representation_controller_still_rejects_step_zero() -> None:
    result = subprocess.run(
        ["bash", str(CONTROLLER), "submit-representation-evaluation"],
        cwd=ROOT,
        env={
            **os.environ,
            "TWIRL_EXPECTED_GIT_SHA": "1" * 40,
            "TWIRL_FM0_STOP_AFTER_STEP": "0",
        },
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "Unauthorized representation milestone" in result.stderr


def test_fm0_2_step0_evaluation_receipt_is_exact_and_keeps_gates_closed() -> None:
    receipts = STEP0_REPORT / "receipts"
    evaluation_path = receipts / "step0_representation_health.json"
    receipt_path = receipts / "step0_representation_health.receipt.json"
    sidecar_path = receipts / "step0_representation_health.receipt.json.sha256"
    evaluation = json.loads(evaluation_path.read_text(encoding="utf-8"))
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))

    assert _sha256(evaluation_path) == STEP0_EVALUATION_SHA
    assert _sha256(receipt_path) == STEP0_RECEIPT_SHA
    assert sidecar_path.read_text(encoding="utf-8") == (
        f"{STEP0_RECEIPT_SHA}  representation_health.receipt.json\n"
    )
    assert receipt["passed"] is True
    assert receipt["evaluation_sha256"] == STEP0_EVALUATION_SHA
    assert receipt["checkpoint_step"] == 0
    assert receipt["checkpoint_seed"] == 560067
    assert receipt["evaluation_role"] == "exact_same_seed_initialization"
    assert receipt["output_field_role_mapping"] == {
        "trained_encoder": "exact_same_seed_initialization",
        "random_encoder_control": "separate_matched_random_control_seed_0",
        "source_paired_trained_minus_random": (
            "exact_initialization_minus_matched_random_control"
        ),
    }
    assert receipt["sealed_test_access_count"] == 0
    assert receipt["development_event_retention_run"] is False
    assert receipt["scientific_go_no_go_applied"] is False
    assert receipt["production_authorized"] is False
    assert receipt["foundation_model_claim_authorized"] is False
    assert evaluation["run"]["global_step"] == 0
    assert evaluation["data_access"]["sealed_test_access_count"] == 0


def test_fm0_2_four_point_table_contains_the_exact_same_seed_trajectory() -> None:
    table = STEP0_REPORT / "step0_step500_step1000_step2000_metrics.csv"
    with table.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert {(int(row["step"]), row["representation"]) for row in rows} == {
        (step, branch)
        for step in (0, 500, 1000, 2000)
        for branch in ("h_window", "z_window")
    }
    z_by_step = {
        int(row["step"]): row for row in rows if row["representation"] == "z_window"
    }
    assert float(z_by_step[0]["effective_rank"]) == pytest.approx(3.0274745240463248)
    assert float(z_by_step[2000]["effective_rank"]) == pytest.approx(39.39935029993325)
    assert float(z_by_step[0]["cross_sector_retrieval_clustered_mean"]) == (
        pytest.approx(0.033796296296296297)
    )
    assert float(z_by_step[2000]["cross_sector_retrieval_clustered_mean"]) == (
        pytest.approx(0.020833333333333332)
    )


@pytest.mark.parametrize(
    ("name", "value", "message"),
    [
        (
            "TWIRL_FM0_ARTIFACT_GIT_SHA",
            "2" * 40,
            "frozen FM0.2.1 training revision",
        ),
        (
            "TWIRL_FM0_INITIALIZATION_CHECKPOINT",
            f"{PROJECT_ROOT}/wrong/checkpoint_step_00000000.pt",
            "exact seed-560067 initialization",
        ),
        (
            "TWIRL_FM0_STEP2000_POST_VALIDATION",
            f"{PROJECT_ROOT}/wrong/post_validation.receipt.json",
            "not authoritative",
        ),
        (
            "TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT",
            f"{PROJECT_ROOT}/wrong/representation_health.receipt.json",
            "not authoritative",
        ),
        (
            "TWIRL_FM0_STEP2000_POST_VALIDATION_SHA256",
            "0" * 64,
            "not the frozen receipt",
        ),
        (
            "TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT_SHA256",
            "0" * 64,
            "not the frozen receipt",
        ),
    ],
)
def test_fm0_2_initialization_controller_rejects_drift_before_ssh(
    name: str,
    value: str,
    message: str,
) -> None:
    environment = _initialization_environment()
    environment[name] = value
    result = subprocess.run(
        ["bash", str(CONTROLLER), "submit-initialization-representation-evaluation"],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert message in result.stderr
    assert "ssh " not in result.stdout
