from __future__ import annotations

from pathlib import Path
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
ORCD_ROOT = REPO_ROOT / "scripts" / "orcd"

FREEZE = ORCD_ROOT / "slurm_teacher_ssl_fullpool_freeze_cpu.sbatch"
BLS = ORCD_ROOT / "slurm_teacher_ssl_fullpool_bls_cpu.sbatch"
BLS_CANARY = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_bls_canary_cpu.sbatch"
)
BLS_MERGE = ORCD_ROOT / "slurm_teacher_ssl_fullpool_bls_merge_cpu.sbatch"
BLS_GLOBAL = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_bls_global_cpu.sbatch"
)
NATIVE = ORCD_ROOT / "slurm_teacher_ssl_fullpool_native_cpu.sbatch"
NATIVE_REGISTRY = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_native_registry_cpu.sbatch"
)
REGISTRY = ORCD_ROOT / "slurm_teacher_ssl_fullpool_registry_cpu.sbatch"
FOLD = ORCD_ROOT / "slurm_teacher_ssl_fullpool_fold_h200.sbatch"
PDO_RAW = (
    REPO_ROOT
    / "scripts/stage5_validation/run_teacher_ssl_fullpool_raw_export_pdo.sh"
)
CONTROLLER = ORCD_ROOT / "run_teacher_ssl_fullpool_orcd.sh"

CPU_WRAPPERS = (
    FREEZE,
    BLS,
    BLS_CANARY,
    BLS_MERGE,
    BLS_GLOBAL,
    NATIVE,
    NATIVE_REGISTRY,
    REGISTRY,
)
ALL_WRAPPERS = (*CPU_WRAPPERS, FOLD, PDO_RAW, CONTROLLER)


@pytest.mark.parametrize("wrapper_path", ALL_WRAPPERS)
def test_fullpool_orcd_wrapper_has_valid_bash_syntax(
    wrapper_path: Path,
) -> None:
    completed = subprocess.run(
        ["bash", "-n", str(wrapper_path)],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize("wrapper_path", CPU_WRAPPERS)
def test_fullpool_cpu_wrappers_request_no_gpu(wrapper_path: Path) -> None:
    wrapper = wrapper_path.read_text(encoding="utf-8")
    assert "#SBATCH -p pg_mki_aryeh" in wrapper
    assert "#SBATCH --gres" not in wrapper
    assert "nvidia-smi" not in wrapper


def test_fullpool_native_wrappers_are_real_only_and_sharded() -> None:
    native = NATIVE.read_text(encoding="utf-8")
    native_registry = NATIVE_REGISTRY.read_text(encoding="utf-8")
    raw = PDO_RAW.read_text(encoding="utf-8")

    assert "#SBATCH --exclude=node4900" in native
    assert "TWIRL_FULLPOOL_NATIVE_SHARDS:?" in native
    assert 'SHARD_INDEX="${SLURM_ARRAY_TASK_ID' in native
    assert "build_teacher_ssl_fullpool_native_shard.py" in native
    assert "--n-periods 4096" in native
    assert 'ORBITID_POLICY="reference_by_cadence"' in native

    assert "export_teacher_ssl_fullpool_raw_sources.py" in raw
    assert 'if [[ "${OUT_H5}" != /pdo/users/tehan/* ]]' in raw
    assert "--shard-index" in raw
    assert "--n-shards" in raw

    assert "TWIRL_FULLPOOL_NATIVE_SHARDS:-16" in native_registry
    assert "for sector in {56..62}; do" in native_registry
    assert 'native_shards+=("${shard}")' in native_registry
    assert "7 * 10#${N_SHARDS}" in native_registry
    assert "build_teacher_ssl_fullpool_native_registry.py" in native_registry


def test_fullpool_bls_arrays_and_s62_policy_are_locked() -> None:
    bls = BLS.read_text(encoding="utf-8")
    merge = BLS_MERGE.read_text(encoding="utf-8")

    assert "#SBATCH --array=0-55%28" in bls
    assert "N_SHARDS=8" in bls
    assert 'ORBITID_POLICY="strict"' in bls
    assert 'if [[ "${SECTOR}" == "62" ]]' in bls
    assert 'ORBITID_POLICY="reference_by_cadence"' in bls
    assert bls.count('ORBITID_POLICY="reference_by_cadence"') == 1
    assert "--n-periods 50000" in bls
    assert "--n-peaks 10" in bls
    assert "--expected-compact-lc-sha256" in bls
    assert "--expected-target-allowlist-sha256" in bls

    assert "#SBATCH --array=0-6%7" in merge
    assert "N_SHARDS=8" in merge
    assert "SECTOR=$((56 + TASK_ID))" in merge


def test_fullpool_freeze_uses_rebuilt_s56_pair() -> None:
    wrapper = FREEZE.read_text(encoding="utf-8")

    assert "s56_A2v1_adp_pair_rebuilt.manifest.json" in wrapper
    assert "s56_A2v1_adp_pair_rebuilt.h5" in wrapper


def test_fullpool_bls_canary_is_bounded_and_output_isolated() -> None:
    wrapper = BLS_CANARY.read_text(encoding="utf-8")

    assert "#SBATCH --array" not in wrapper
    assert "TWIRL_SSL_FULLPOOL_BLS_CANARY_TARGETS:-100" in wrapper
    assert "10#${N_TARGETS}" in wrapper
    assert "> 100" in wrapper
    assert "TWIRL_SSL_FULLPOOL_BLS_CANARY_SECTOR:-62" in wrapper
    assert "Canary and production BLS roots must be distinct." in wrapper
    assert 'CANARY_ROOT="${TWIRL_SSL_FULLPOOL_BLS_CANARY_ROOT' in wrapper
    assert '"selection": "first_n_sorted_tics"' in wrapper
    assert '"parent_allowlist"' in wrapper
    assert "--expected-compact-lc-sha256" in wrapper
    assert "--expected-target-allowlist-sha256" in wrapper
    assert "--n-periods 50000" in wrapper
    assert "--n-peaks 10" in wrapper
    assert "--n-shards 1" in wrapper
    assert "A2v1-fullpool-v1-canary" in wrapper


def test_global_bls_wrapper_requires_all_seven_sector_products() -> None:
    wrapper = BLS_GLOBAL.read_text(encoding="utf-8")
    assert "for sector in {56..62}; do" in wrapper
    assert 'args+=(--sector-bls "${sector}=${sector_table}")' in wrapper
    assert (
        'args+=(--sector-summary "${sector}=${sector_summary}")'
        in wrapper
    )
    assert "consolidate_teacher_ssl_full_pool_bls.py" in wrapper
    assert "--pool-summary" in wrapper
    assert "--out-parquet" in wrapper
    assert "--out-summary" in wrapper


def test_registry_wrapper_binds_every_external_authority() -> None:
    wrapper = REGISTRY.read_text(encoding="utf-8")
    assert (
        "TWIRL_SSL_FULLPOOL_NATIVE_REGISTRY:?"
        in wrapper
    )
    assert "build_teacher_ssl_fullpool_registry.py" in wrapper
    assert "--frozen-pool" in wrapper
    assert "--frozen-pool-summary" in wrapper
    assert "--bls-summary" in wrapper
    assert "--bls-peaks" in wrapper
    assert "--native-registry" in wrapper
    assert "--frozen-split-registry" in wrapper
    assert "--reserved-hosts" in wrapper
    assert "--registry-out" in wrapper
    assert "--summary-out" in wrapper


def test_fullpool_fold_wrapper_uses_four_way_one_h200_array() -> None:
    wrapper = FOLD.read_text(encoding="utf-8")
    assert "#SBATCH --array=0-4%4" in wrapper
    assert "#SBATCH --gres=gpu:h200:1" in wrapper
    assert wrapper.count("#SBATCH --gres=gpu:h200:1") == 1
    assert "#SBATCH -p pg_mki_aryeh" in wrapper
    assert "train_teacher_ssl_fullpool_fold.py" in wrapper
    assert '--fold "${SLURM_ARRAY_TASK_ID}"' in wrapper


def test_fullpool_fold_wrapper_has_fail_closed_bounded_smoke() -> None:
    wrapper = FOLD.read_text(encoding="utf-8")

    assert "TWIRL_SSL_FULLPOOL_MAX_ROWS" in wrapper
    assert "10#${MAX_ROWS}" in wrapper
    assert "> 10000" in wrapper
    assert 'if [[ "${EPOCHS}" != "1" ]]' in wrapper
    assert "TWIRL_SSL_FULLPOOL_SMOKE_OUT_ROOT" in wrapper
    assert "Smoke and production output roots must be distinct." in wrapper
    assert 'TRAIN_OUT_ROOT="${SMOKE_OUT_ROOT}"' in wrapper
    assert 'smoke_args+=(--max-rows "${MAX_ROWS}")' in wrapper
    assert '--out-root "${TRAIN_OUT_ROOT}"' in wrapper


def test_fullpool_controller_is_noninteractive_and_canary_gated() -> None:
    controller = CONTROLLER.read_text(encoding="utf-8")

    assert "-o PasswordAuthentication=no" in controller
    assert "-o KbdInteractiveAuthentication=no" in controller
    assert "-o NumberOfPasswordPrompts=0" in controller
    assert "-o ControlMaster=no" in controller
    assert '-O check "${ORCD_HOST}"' in controller
    assert "TWIRL_EXPECTED_GIT_SHA" in controller
    assert "worktree add --detach" in controller
    assert "stage-preregister" in controller
    assert "submit-preprocess" in controller
    assert "slurm_teacher_ssl_fullpool_freeze_cpu.sbatch" in controller
    assert "slurm_teacher_ssl_fullpool_bls_canary_cpu.sbatch" in controller
    assert 'cd \\"\\${repo}\\"' in controller
    assert controller.index('cd \\"\\${repo}\\"') < controller.index(
        "freeze=\\$(sbatch"
    )
    assert "--dependency=afterok:\\${freeze}" in controller
    assert "--dependency=afterok:\\${canary}" in controller
    assert "--dependency=afterok:\\${bls}" in controller
    assert "--dependency=afterok:\\${merge}" in controller
