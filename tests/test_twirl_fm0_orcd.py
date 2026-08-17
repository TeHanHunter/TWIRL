from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import yaml


ROOT = Path(__file__).resolve().parents[1]
ORCD = ROOT / "scripts" / "orcd"
OP_CONFIG = ROOT / "configs" / "orcd" / "twirl_fm0_1_s56_s67_poc.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text(name: str) -> str:
    return (ORCD / name).read_text(encoding="utf-8")


def test_operational_config_binds_the_byte_exact_science_freeze() -> None:
    payload = json.loads(OP_CONFIG.read_text(encoding="utf-8"))
    bindings = payload["science_bindings"]
    for binding in bindings.values():
        assert sha256(ROOT / binding["path"]) == binding["sha256"]

    assert payload["code_binding"] == {
        "source": "TWIRL_EXPECTED_GIT_SHA",
        "format": "full_lowercase_40_hex",
        "clean_worktree_required": True,
        "detached_deployment_required": True,
    }
    assert payload["input_release_binding"]["receipt_must_be_read_only"] is True


def test_admission_reserves_stage1_before_exactly_one_fm_h200() -> None:
    payload = json.loads(OP_CONFIG.read_text(encoding="utf-8"))
    assert payload["scheduler"]["partition"] == "pg_mki_aryeh"
    assert payload["scheduler"]["configured_h200"] == 8
    assert payload["scheduler"]["gpu_array_forbidden"] is True
    assert payload["scheduler"]["multi_gpu_forbidden"] is True
    assert payload["stage1_reservation"] == {
        "job_name_regex": "^twirl-s[0-9]+-o[0-9]+-c[1-4]d[1-4]-a2v1(?:-p4)?$",
        "lanes": 4,
        "h200_per_lane": 1,
        "cpus_per_lane": 16,
        "memory_mib_per_lane": 393216,
        "live_cpu_ceiling": 78,
        "dependency_pending_reason": "Dependency",
    }
    assert payload["jobs"]["fp32_smoke_h200"]["h200"] == 1
    assert payload["jobs"]["fp32_smoke_h200"]["cpus"] == 8
    assert payload["jobs"]["fp32_smoke_h200"]["memory_mib"] == 98304
    assert payload["jobs"]["fp32_smoke_h200"]["walltime"] == "02:00:00"
    assert payload["admission"]["ttl_seconds"] == 300
    assert payload["admission"]["inspect_all_partition_jobs"] is True


def test_controller_is_noninteractive_dry_run_by_default() -> None:
    controller = ORCD / "run_twirl_fm0_1_poc_orcd.sh"
    result = subprocess.run(
        ["bash", str(controller), "probe"],
        cwd=ROOT,
        env={"PATH": "/usr/bin:/bin", "TWIRL_EXPECTED_GIT_SHA": "a" * 40},
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.startswith("+ ssh ")
    script = controller.read_text(encoding="utf-8")
    assert "DRY_RUN=1" in script
    assert "BatchMode=yes" in script
    assert "PasswordAuthentication=no" in script
    assert "KbdInteractiveAuthentication=no" in script
    assert "NumberOfPasswordPrompts=0" in script
    assert "ControlPath=${CONTROL_PATH}" in script


def test_deploy_materializes_the_frozen_receipt_in_sparse_worktrees() -> None:
    script = text("run_twirl_fm0_1_poc_orcd.sh")
    receipt = "reports/stage5_validation/twirl_fm0_1_design_freeze_v1"
    assert "sparse-checkout list" in script
    assert f"sparse-checkout add {receipt}" in script
    assert f"test -s '${{REMOTE_REPO}}/{receipt}/freeze.json'" in script


def test_environment_builder_binds_current_operational_config() -> None:
    builder = text("build_orcd_fm0_env.sh")
    assert f'EXPECTED_ORCD_CONFIG_SHA256="{sha256(OP_CONFIG)}"' in builder


def test_controller_admission_uses_all_queue_and_node_tres() -> None:
    script = text("run_twirl_fm0_1_poc_orcd.sh")
    assert 'i.get(\\"schema_version\\")' in script
    assert 'l[\\"input_receipt\\"][\\"sha256\\"]' in script
    assert 'tres(alloc_tres, "gres/gpu:h200", zero_if_absent=True)' in script
    assert 'ADMISSION_TMP="${tmp}"' in script
    assert 'ADMISSION_TMP=""' in script
    assert "squeue -p '${PARTITION}'" in script
    assert "-t RUNNING,PENDING,CONFIGURING,COMPLETING" in script
    assert "scontrol show node -o '${GPU_NODE}'" in script
    assert 'field(node_raw, "AllocTRES")' in script
    assert 'field(node_raw, "CfgTRES")' in script
    assert 'runnable Stage-1 job is pending; FM admission denied' in script
    assert "another FM0.1 GPU smoke is already live" in script
    assert '"expires_at_epoch": now+300' in script
    assert '"memory_mib": 1572864' in script
    assert "sbatch --parsable" in script
    assert '"${LAUNCH_DIR}/fp32-smoke.job"' in script
    assert 'fp32-smoke-${remote_epoch}.job' not in script


def test_registry_stage_requires_and_exports_admitted_observations() -> None:
    controller = text("run_twirl_fm0_1_poc_orcd.sh")
    wrapper = text("slurm_twirl_fm0_1_registry_cpu.sbatch")
    assert "TWIRL_FM0_OBSERVATIONS_TABLE:?set staged admitted-observations table" in controller
    assert "TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY" in controller
    assert "TWIRL_FM0_OBSERVATIONS_TABLE:?stage the admitted-observations table" in wrapper
    assert "TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY:?stage the A2v1 HDF5 source inventory" in wrapper
    assert "TWIRL_FM0_SOURCE_MANIFEST" in controller
    assert "TWIRL_FM0_SOURCE_ADAPTER" in controller
    assert "a2v1_hdf5_quality_aware_v1" in controller
    assert '--observations "${OBSERVATIONS}"' in wrapper
    assert '--a2v1-hdf5-source-inventory "${SOURCE_INVENTORY}"' in wrapper
    assert '"${OUT}/observations.csv"' in wrapper
    assert '"${OUT}/a2v1_hdf5_manifest.csv"' in wrapper


def test_slurm_wrappers_enforce_cpu_gpu_separation_and_claim_limit() -> None:
    cpu_wrappers = (
        "slurm_twirl_fm0_1_registry_cpu.sbatch",
        "slurm_twirl_fm0_1_input_validation_cpu.sbatch",
        "slurm_twirl_fm0_1_loader_smoke_cpu.sbatch",
        "slurm_twirl_fm0_1_post_validation_cpu.sbatch",
    )
    for name in cpu_wrappers:
        script = text(name)
        assert "#SBATCH -p pg_mki_aryeh" in script
        assert "#SBATCH --exclude=node4900" in script
        assert "#SBATCH --gres" not in script

    smoke = text("slurm_twirl_fm0_1_fp32_smoke_h200.sbatch")
    assert "#SBATCH --gres=gpu:h200:1" in smoke
    assert "#SBATCH -c 8" in smoke
    assert "#SBATCH --mem=96G" in smoke
    assert "#SBATCH -t 02:00:00" in smoke
    assert "#SBATCH --array" not in smoke
    assert "--precision fp32" in smoke
    assert "--synthetic-smoke" in smoke
    assert "not a real-data FM0.1 result" in smoke
    assert "TWIRL_FM0_ADMISSION_RECEIPT_SHA256" in smoke
    assert "TWIRL_FM0_INPUT_RECEIPT_SHA256" in smoke

    post = text("slurm_twirl_fm0_1_post_validation_cpu.sbatch")
    assert "validate_twirl_fm0_release.py" in post
    assert 'run_validation.get("checkpoint_inspected") is not True' in post
    assert "independent_run_validation_sha256" in post

    loader = text("slurm_twirl_fm0_1_loader_smoke_cpu.sbatch")
    assert '"${PYTHON}" -m pytest -q' in loader
    assert "tests/test_twirl_fm0_model.py" in loader
    assert "tests/test_twirl_fm0_training.py" in loader
    assert "tests/test_twirl_fm0_real_dataset.py" in loader
    assert "input_builder_git_sha" in loader

    real = text("slurm_twirl_fm0_1_real_train_h200.sbatch")
    assert "#SBATCH --gres=gpu:h200:1" in real
    assert "#SBATCH -c 8" in real
    assert "#SBATCH --mem=96G" in real
    assert "#SBATCH -t 47:30:00" in real
    assert "#SBATCH --array" not in real
    assert "--variant TWIRL-FM0.1.1" in real
    assert '--input-release "${RELEASE_ROOT}"' in real
    assert '--max-steps "${TARGET_STEP}"' in real
    assert "--synthetic-smoke" not in real
    assert "scientific_training_eligible" in real
    assert "TWIRL_FM0_ADMISSION_RECEIPT_SHA256" in real


def test_input_gate_uses_correct_local_time_and_exact_shard_allowlist() -> None:
    script = text("slurm_twirl_fm0_1_input_validation_cpu.sbatch")
    assert '"local_time_cadences"' in script
    assert '"elapsed_time_cadences"' not in script
    assert '"flux": (n, 6)' in script
    assert '"flux_error": (n, 2)' in script
    assert '"view_present": (6,)' in script
    assert '"certifies_full_campaign": False' in script
    assert '"foundation_model_claim_authorized": False' in script
    assert "--a2v1-hdf5-manifest" in script
    assert 'adapter == "a2v1_hdf5_quality_aware_v1"' in script


def test_environment_is_versioned_and_torch_is_pinned() -> None:
    env = yaml.safe_load(
        (ROOT / "configs" / "orcd" / "twirl-fm0-poc-env.yml").read_text(
            encoding="utf-8"
        )
    )
    assert "python=3.11.9" in env["dependencies"]
    assert "numpy=1.26.4" in env["dependencies"]
    assert "pytest=8.4.1" in env["dependencies"]
    builder = text("build_orcd_fm0_env.sh")
    assert 'readonly TORCH_VERSION="2.11.0"' in builder
    assert 'readonly TORCH_INDEX_URL="https://download.pytorch.org/whl/cu128"' in builder
    assert "prefix exists without a completed manifest" in builder
    assert "status --porcelain=v1 --untracked-files=all" in builder
    assert 'git -C "${REPO}" archive "${EXPECTED_SHA}"' in builder
    assert '--no-build-isolation "${BUILD_SOURCE}"' in builder
    assert '--no-build-isolation "${REPO}"' not in builder
    assert '--no-build-isolation -e "${REPO}"' not in builder


def test_all_new_shell_boundaries_parse_with_bash() -> None:
    paths = [
        ORCD / "build_orcd_fm0_env.sh",
        ORCD / "run_twirl_fm0_1_poc_orcd.sh",
        *sorted(ORCD.glob("slurm_twirl_fm0_1_*.sbatch")),
    ]
    subprocess.run(["bash", "-n", *map(str, paths)], cwd=ROOT, check=True)
