from __future__ import annotations

from pathlib import Path
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
ORCD_ROOT = REPO_ROOT / "scripts" / "orcd"

ELIGIBILITY = (
    ORCD_ROOT
    / "slurm_teacher_ssl_fullpool_native_eligibility_cpu.sbatch"
)
NATIVE = ORCD_ROOT / "slurm_teacher_ssl_fullpool_native_cpu.sbatch"
NATIVE_REGISTRY = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_native_registry_cpu.sbatch"
)
SSL_REGISTRY = ORCD_ROOT / "slurm_teacher_ssl_fullpool_registry_cpu.sbatch"
NUMERIC_AUDIT = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_numeric_audit_cpu.sbatch"
)
NUMERIC_RELEASE = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_numeric_release_cpu.sbatch"
)
SMOKE = ORCD_ROOT / "slurm_teacher_ssl_fullpool_smoke_h200.sbatch"
FOLDS = ORCD_ROOT / "slurm_teacher_ssl_fullpool_fold_h200.sbatch"
CONTROLLER = ORCD_ROOT / "run_teacher_ssl_fullpool_v2_orcd.sh"

ASSETS = (
    ELIGIBILITY,
    NATIVE,
    NATIVE_REGISTRY,
    SSL_REGISTRY,
    NUMERIC_AUDIT,
    NUMERIC_RELEASE,
    SMOKE,
    FOLDS,
    CONTROLLER,
)

V1_ROOT = (
    "teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
)
V2_ROOT = (
    "teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible"
)


@pytest.mark.parametrize("asset", ASSETS)
def test_fullpool_v2_orcd_assets_have_valid_bash_syntax(
    asset: Path,
) -> None:
    completed = subprocess.run(
        ["bash", "-n", str(asset)],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "asset",
    (
        ELIGIBILITY,
        NATIVE,
        NATIVE_REGISTRY,
        SSL_REGISTRY,
        NUMERIC_AUDIT,
        NUMERIC_RELEASE,
        SMOKE,
        FOLDS,
    ),
)
def test_v2_jobs_require_the_exact_clean_deployed_revision(
    asset: Path,
) -> None:
    wrapper = asset.read_text(encoding="utf-8")
    assert "TWIRL_EXPECTED_GIT_SHA" in wrapper
    assert 'git -C "${REPO}" rev-parse HEAD' in wrapper
    assert (
        'git -C "${REPO}" status --porcelain=v1 --untracked-files=all'
        in wrapper
    )


def test_eligibility_job_binds_complete_v1_authority_into_fresh_v2() -> None:
    wrapper = ELIGIBILITY.read_text(encoding="utf-8")

    assert "#SBATCH --exclude=node4900" in wrapper
    assert "#SBATCH --gres" not in wrapper
    assert V1_ROOT in wrapper
    assert V2_ROOT in wrapper
    assert "for sector in {56..62}; do" in wrapper
    assert "teacher_ssl_full_pool_observations.csv" in wrapper
    assert "teacher_ssl_full_pool_observations.parquet" not in wrapper
    assert "teacher_ssl_full_pool_manifest.summary.json" in wrapper
    assert "teacher_ssl_full_pool_bls.parquet" in wrapper
    assert "teacher_ssl_full_pool_bls.summary.json" in wrapper
    assert "build_teacher_ssl_fullpool_native_eligibility.py" in wrapper
    assert "--frozen-pool" in wrapper
    assert "--frozen-pool-summary" in wrapper
    assert "--bls-peaks" in wrapper
    assert "--bls-summary" in wrapper
    assert "--exclusions-out" in wrapper
    assert "--summary-out" in wrapper
    assert "--production-lock" in wrapper
    assert "Refusing to overwrite" in wrapper


@pytest.mark.parametrize(
    "asset",
    (ELIGIBILITY, NATIVE, NATIVE_REGISTRY, SSL_REGISTRY),
)
def test_v2_pool_consumers_use_the_raw_release_csv_authority(
    asset: Path,
) -> None:
    wrapper = asset.read_text(encoding="utf-8")

    assert "teacher_ssl_full_pool_observations.csv" in wrapper
    assert "teacher_ssl_full_pool_observations.parquet" not in wrapper


def test_native_v2_reuses_only_immutable_v1_raw_sources() -> None:
    wrapper = NATIVE.read_text(encoding="utf-8")

    assert V1_ROOT in wrapper
    assert V2_ROOT in wrapper
    assert (
        'RAW_DIR="${SOURCE_RUN_ROOT}/native_preparation/raw_sources"'
        in wrapper
    )
    assert "--raw-source-summary" in wrapper
    assert "--raw-export-complete" in wrapper
    assert "--raw-transfer-validation" in wrapper
    assert '${RUN}/native/s${SECTOR}/native_' in wrapper
    assert 'd["inputs"]["compact_exports"]' in wrapper
    assert (
        "teacher_v3_s56_s62_pdo_inputs_33f5855e_extract_"
        "orcdlane_20260726T0143/a2v1_cadence_reference"
        in wrapper
    )
    assert "--eligibility-exclusions" in wrapper
    assert "--eligibility-summary" in wrapper
    assert '"${RUN}" != "${EXPECTED_RUN}"' in wrapper
    assert '"${OUT_H5}" != "${RUN}/"*' in wrapper
    assert "Never impute" not in wrapper


def test_native_registry_requires_all_shards_and_their_summaries() -> None:
    wrapper = NATIVE_REGISTRY.read_text(encoding="utf-8")

    assert V2_ROOT in wrapper
    assert "for sector in {56..62}; do" in wrapper
    assert "native_shards=()" in wrapper
    assert "native_summaries=()" in wrapper
    assert 'native_shards+=("${shard}")' in wrapper
    assert 'native_summaries+=("${shard_summary}")' in wrapper
    assert "7 * 10#${N_SHARDS}" in wrapper
    assert "--eligibility-exclusions" in wrapper
    assert "--eligibility-summary" in wrapper
    assert "--native-shard-summaries" in wrapper
    assert "--release-summary-out" in wrapper
    assert "native_fullpool_release.summary.json" in wrapper


def test_final_registry_binds_every_v2_release_summary() -> None:
    wrapper = SSL_REGISTRY.read_text(encoding="utf-8")

    assert V1_ROOT in wrapper
    assert V2_ROOT in wrapper
    assert "for sector in {56..62}; do" in wrapper
    assert "--bls-peaks" in wrapper
    assert "--bls-summary" in wrapper
    assert "--eligibility-exclusions" in wrapper
    assert "--eligibility-summary" in wrapper
    assert "--native-registry" in wrapper
    assert "--native-registry-summary" in wrapper
    assert "--native-release-summary" in wrapper
    assert "native_fullpool_release.summary.json" in wrapper


def test_numeric_audit_maps_exact_112_cpu_tasks_to_seven_by_sixteen() -> None:
    wrapper = NUMERIC_AUDIT.read_text(encoding="utf-8")

    assert "#SBATCH --gres" not in wrapper
    assert "#SBATCH --array=0-111%4" in wrapper
    assert "readonly N_SHARDS=16" in wrapper
    assert "readonly N_TASKS=112" in wrapper
    assert "56 + SLURM_ARRAY_TASK_ID / N_SHARDS" in wrapper
    assert "SLURM_ARRAY_TASK_ID % N_SHARDS" in wrapper
    assert "audit_teacher_ssl_fullpool_numeric.py" in wrapper
    assert "--native-registry" in wrapper
    assert "--native-registry-summary" in wrapper
    assert "--native-release-summary" in wrapper
    assert "--registry" in wrapper
    assert "--registry-summary" in wrapper
    assert "--sector" in wrapper
    assert "--shard-index" in wrapper
    assert "--expected-code-revision" in wrapper
    assert "--report-out" in wrapper
    assert (
        "frozen/model_input_numeric_gate_v1"
        in wrapper
    )
    assert "native_${SHARD_INDEX}.numeric.json" in wrapper
    assert "already complete report=" in wrapper
    assert "recovering incomplete report=" in wrapper
    assert 'payload.get("passed") is True' in wrapper
    assert '"${REPORT_OUT}.sha256"' in wrapper


def test_numeric_release_requires_and_enumerates_all_112_reports() -> None:
    wrapper = NUMERIC_RELEASE.read_text(encoding="utf-8")

    assert "#SBATCH --gres" not in wrapper
    assert "#SBATCH --array" not in wrapper
    assert "for sector in {56..62}; do" in wrapper
    assert "for shard_index in {0..15}; do" in wrapper
    assert 'reports+=("${report}")' in wrapper
    assert '"${#reports[@]}" -ne 112' in wrapper
    assert "merge_teacher_ssl_fullpool_numeric.py" in wrapper
    assert '--shard-reports "${reports[@]}"' in wrapper
    assert "--release-out" in wrapper
    assert "model_input_numeric_gate.release.json" in wrapper
    assert "validate_numeric_gate_release" in wrapper
    assert "expected_code_revision=sys.argv[2]" in wrapper
    assert "already complete release=" in wrapper
    assert "recovering incomplete publication" in wrapper
    assert "model_input_numeric_audit.parquet" in wrapper
    assert '"${report}.sha256"' in wrapper


def test_bounded_smoke_is_one_fold_one_epoch_and_one_h200() -> None:
    wrapper = SMOKE.read_text(encoding="utf-8")

    assert "#SBATCH --gres=gpu:h200:1" in wrapper
    assert "#SBATCH --array" not in wrapper
    assert "readonly EPOCHS=1" in wrapper
    assert "readonly MAX_ROWS=4096" in wrapper
    assert "readonly BATCH_SIZE=64" in wrapper
    assert "readonly WORKERS=8" in wrapper
    assert "readonly SEED=560064" in wrapper
    assert "readonly FOLD=2" in wrapper
    assert (
        'readonly REQUIRED_OBSERVATION_ID="s0060-tic0000000722078603"'
        in wrapper
    )
    assert "--epochs \"${EPOCHS}\"" in wrapper
    assert '--max-rows "${MAX_ROWS}"' in wrapper
    assert (
        '--required-observation-id "${REQUIRED_OBSERVATION_ID}"'
        in wrapper
    )
    assert "Smoke and production output roots must be distinct." in wrapper
    assert "Smoke roots differ from the exact Teacher v4-SSL v2 contract." in wrapper
    assert "--eligibility-exclusions" in wrapper
    assert "--eligibility-summary" in wrapper
    assert "--native-registry" in wrapper
    assert "--native-registry-summary" in wrapper
    assert "--native-release-summary" in wrapper
    assert "native_model_exclusions.csv" in wrapper
    assert "native_fullpool_release.summary.json" in wrapper
    assert "model_input_numeric_gate.release.json" in wrapper
    assert "model_runs/effective_quality_mask_v1" in wrapper
    assert "--numeric-gate-release" in wrapper
    assert "validate_numeric_gate_release" in wrapper
    assert "expected_code_revision=sys.argv[2]" in wrapper
    assert "validate_teacher_ssl_fullpool_v2_smoke.py" in wrapper
    assert "--expected-code-revision" in wrapper
    assert 'fold_${FOLD}/summary.json' in wrapper
    assert '--expected-fold "${FOLD}"' in wrapper
    assert "TWIRL_SSL_FULLPOOL_SMOKE_MAX_ROWS \\" in wrapper
    assert "TWIRL_SSL_FULLPOOL_SMOKE_FOLD \\" in wrapper
    assert "TWIRL_SSL_FULLPOOL_SMOKE_FOLD:-" not in wrapper
    assert (
        "TWIRL_SSL_FULLPOOL_SMOKE_REQUIRED_OBSERVATION_ID \\"
        in wrapper
    )
    assert "Smoke contract forbids environment override" in wrapper


def test_production_folds_are_five_one_h200_jobs_capped_at_four() -> None:
    wrapper = FOLDS.read_text(encoding="utf-8")

    assert "#SBATCH --gres=gpu:h200:1" in wrapper
    assert wrapper.count("#SBATCH --gres=gpu:h200:1") == 1
    assert "#SBATCH --array=0-4%4" in wrapper
    assert V2_ROOT in wrapper
    assert "readonly EPOCHS=20" in wrapper
    assert "readonly BATCH_SIZE=64" in wrapper
    assert "readonly WORKERS=8" in wrapper
    assert "readonly SEED=560064" in wrapper
    assert "--eligibility-exclusions" in wrapper
    assert "--eligibility-summary" in wrapper
    assert "--native-registry" in wrapper
    assert "--native-registry-summary" in wrapper
    assert "--native-release-summary" in wrapper
    assert "native_model_exclusions.csv" in wrapper
    assert "native_fullpool_release.summary.json" in wrapper
    assert "model_input_numeric_gate.release.json" in wrapper
    assert "model_runs/effective_quality_mask_v1" in wrapper
    assert "--numeric-gate-release" in wrapper
    assert "validate_numeric_gate_release" in wrapper
    assert "expected_code_revision=sys.argv[2]" in wrapper
    assert "training/five_fold" in wrapper
    assert "Training output may not be published under the v1 run root." in wrapper
    assert "Fold roots differ from the exact Teacher v4-SSL v2 contract." in wrapper
    assert "TWIRL_SSL_FULLPOOL_MAX_ROWS \\" in wrapper
    assert "TWIRL_SSL_FULLPOOL_SMOKE_OUT_ROOT \\" in wrapper
    assert "Production fold contract forbids environment override" in wrapper
    assert "--max-rows" not in wrapper
    assert "smoke_args" not in wrapper
    assert "validate_teacher_ssl_fullpool_v2_smoke.py" in wrapper
    assert (
        "model_runs/effective_quality_mask_v1"
        in wrapper
    )
    assert "smoke/one_epoch/encoder_pretraining/fold_2/summary.json" in wrapper
    assert "--expected-fold 2" in wrapper


def test_v2_controller_has_explicit_nonchained_gates() -> None:
    controller = CONTROLLER.read_text(encoding="utf-8")

    for action in (
        "submit-eligibility",
        "submit-native-canary",
        "submit-native",
        "submit-native-registry",
        "submit-ssl-registry",
        "submit-numeric-audit",
        "submit-numeric-release",
        "submit-smoke",
        "submit-folds",
    ):
        assert action in controller
    assert "--dependency=" not in controller
    assert "-o PasswordAuthentication=no" in controller
    assert "-o KbdInteractiveAuthentication=no" in controller
    assert "-o NumberOfPasswordPrompts=0" in controller
    assert '-O check "${ORCD_HOST}"' in controller
    assert "worktree add --detach" in controller
    assert "codex/teacher-ssl-quality-mask" in controller
    assert "TWIRL_FULLPOOL_SHARD_INDEX=3" in controller
    assert "TIC 2019898202" in controller
    assert (
        "8e9e9c12a24d5ebc7be94b03a4e35411cd10066d62a87d921a8443b06cc188d1"
        in controller
    )
    assert (
        "6ddc8e57bb5fb938ce05389c1629221c73e0e73ac3bf40da47a2019e1a5660e6"
        in controller
    )
    assert (
        "b9f536144265e54a70bff17c782e0668ddbd96efcdf6c223ebc58f46edb7d976"
        in controller
    )
    assert "NATIVE_SHARDS=16" in controller
    assert "--array=0-15%4" in controller
    assert "for sector in {56..62}; do" in controller
    assert "verification" in controller
    assert 'n_shard_excluded_observations\\") == 1' in controller
    assert (
        'n_source_shard_observations\\") == '
        'd.get(\\"n_shard_observations\\") + 1'
        in controller
    )
    assert (
        "ddda9b053eddc744e5032ba350598d58d74cea4f4cc5cd705932ccf6e41022ab"
        in controller
    )
    assert "validate_teacher_ssl_fullpool_v2_smoke.py" in controller
    assert "--expected-code-revision" in controller
    assert "--expected-fold" in controller
    assert "readonly SMOKE_FOLD=2" in controller
    assert "Controller forbids a smoke-fold environment override." in controller
    assert "TWIRL_SSL_FULLPOOL_SMOKE_FOLD=${SMOKE_FOLD}" not in controller
    assert "TWIRL_SSL_FULLPOOL_OUT_ROOT=" not in controller
    assert "slurm_teacher_ssl_fullpool_smoke_h200.sbatch" in controller
    assert "slurm_teacher_ssl_fullpool_fold_h200.sbatch" in controller
    assert "slurm_teacher_ssl_fullpool_numeric_audit_cpu.sbatch" in controller
    assert "slurm_teacher_ssl_fullpool_numeric_release_cpu.sbatch" in controller
    assert "model_input_numeric_gate.release.json" in controller
    assert "model_runs/effective_quality_mask_v1" in controller
    assert "--numeric-gate-release" in controller
    assert "validate_numeric_gate_release" in controller
    assert "expected_code_revision=sys.argv[2]" in controller
    assert "test ! -e '${NUMERIC_SHARD_DIR}" not in controller
    assert controller.count("test ! -e '${NUMERIC_AUDIT}'") == 1
    assert controller.count("test ! -e '${NUMERIC_RELEASE}'") == 1


@pytest.mark.parametrize(
    "action",
    (
        "probe",
        "deploy",
        "submit-eligibility",
        "submit-native-canary",
        "submit-native",
        "submit-native-registry",
        "submit-ssl-registry",
        "submit-numeric-audit",
        "submit-numeric-release",
        "submit-smoke",
        "submit-folds",
        "status",
    ),
)
def test_v2_controller_actions_are_portable_dry_runs(action: str) -> None:
    completed = subprocess.run(
        [str(CONTROLLER), action],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.startswith("+ ssh ")
