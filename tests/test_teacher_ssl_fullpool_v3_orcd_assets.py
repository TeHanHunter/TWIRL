from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import importlib.util

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
ORCD_ROOT = REPO_ROOT / "scripts" / "orcd"
STAGE5_ROOT = REPO_ROOT / "scripts" / "stage5_validation"

CONTROLLER = ORCD_ROOT / "run_teacher_ssl_fullpool_v3_orcd.sh"
CANARY_AUDIT = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_canary_numeric_cpu.sbatch"
)
CANARY_AUDIT_PY = ORCD_ROOT / "audit_teacher_ssl_fullpool_v3_canary.py"
NATIVE = ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_native_cpu.sbatch"
NATIVE_REGISTRY = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_native_registry_cpu.sbatch"
)
SSL_REGISTRY = ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_registry_cpu.sbatch"
NUMERIC_AUDIT = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_numeric_audit_cpu.sbatch"
)
NUMERIC_RELEASE = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_numeric_release_cpu.sbatch"
)
SMOKE = ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_smoke_h200.sbatch"
FOLDS = ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_fold_h200.sbatch"
FOLD_VALIDATION = (
    ORCD_ROOT / "slurm_teacher_ssl_fullpool_v3_validate_folds_cpu.sbatch"
)

ASSETS = (
    CONTROLLER,
    CANARY_AUDIT,
    NATIVE,
    NATIVE_REGISTRY,
    SSL_REGISTRY,
    NUMERIC_AUDIT,
    NUMERIC_RELEASE,
    SMOKE,
    FOLDS,
    FOLD_VALIDATION,
)

HISTORICAL_V3_ROOT = (
    "teacher_ssl_fullpool_v3_s56_s62_a2v1_"
    "effective_quality_adp_bls_eligible"
)
V3R1_ROOT = (
    "teacher_ssl_fullpool_v3r1_s56_s62_a2v1_"
    "effective_quality_adp_btjd_bls_eligible"
)
V3R1_RUN_ROOT_ENV = "TWIRL_SSL_FULLPOOL_V3R1_RUN_ROOT"
MODEL_NAMESPACE = "model_runs/effective_quality_adp_btjd_v2"
NUMERIC_GATE_NAMESPACE = "model_input_numeric_gate_v3r1"
COMPLETION_RELEASE = "teacher_ssl_fullpool_v3r1_five_fold.complete.json"
NATIVE_CONTRACT = (
    "twirl_teacher_ssl_fullpool_real_native_v3r1_"
    "effective_quality_adp_btjd"
)
RELEASE_BINDING = (
    "teacher_ssl_fullpool_v3r1_s56_s62_a2v1_"
    "effective_quality_adp_btjd_bls_eligible"
)
BUILDER_CONTRACT = (
    "twirl_teacher_ssl_fullpool_effective_quality_adp_btjd_builder_v2"
)
DETREND_CONTRACT = "twirl_fs_adp03q_effective_quality_btjd_v2"
DETREND_TIME_CONTRACT = (
    "twirl_teacher_ssl_fullpool_adp_detrend_time_btjd_v1"
)
DETREND_TIME_DATASET = "adp_detrend_time_btjd"
WARNING_CAPTURE_POLICY = (
    "twirl_teacher_ssl_fullpool_exact_numpy_rankwarning_capture_v1"
)
EMPTY_WARNING_LEDGER_SHA256 = (
    "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945"
)


@pytest.mark.parametrize("asset", ASSETS)
def test_v3_assets_have_valid_bash_syntax(asset: Path) -> None:
    completed = subprocess.run(
        ["bash", "-n", str(asset)],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


def test_v3_canary_auditor_has_valid_python_syntax() -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "py_compile", str(CANARY_AUDIT_PY)],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "asset",
    (
        CANARY_AUDIT,
        NATIVE,
        NATIVE_REGISTRY,
        SSL_REGISTRY,
        NUMERIC_AUDIT,
        NUMERIC_RELEASE,
        SMOKE,
        FOLDS,
        FOLD_VALIDATION,
    ),
)
def test_v3r1_jobs_require_exact_clean_revision_and_fresh_root(
    asset: Path,
) -> None:
    wrapper = asset.read_text(encoding="utf-8")
    assert V3R1_ROOT in wrapper
    assert V3R1_RUN_ROOT_ENV in wrapper
    assert HISTORICAL_V3_ROOT not in wrapper
    assert "TWIRL_SSL_FULLPOOL_V3_RUN_ROOT" not in wrapper
    assert "TWIRL_EXPECTED_GIT_SHA" in wrapper
    assert 'git -C "${REPO}" rev-parse HEAD' in wrapper
    assert (
        'git -C "${REPO}" status --porcelain=v1 --untracked-files=all'
        in wrapper
    )
    assert "--resume" not in wrapper


def test_native_v3r1_uses_raw_v1_and_btjd_effective_quality_contract() -> None:
    wrapper = NATIVE.read_text(encoding="utf-8")
    assert "teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp" in wrapper
    assert V3R1_ROOT in wrapper
    assert V3R1_RUN_ROOT_ENV in wrapper
    assert 'RAW_DIR="${SOURCE_RUN_ROOT}/native_preparation/raw_sources"' in wrapper
    assert "--raw-source" in wrapper
    assert "--compact-adp-h5" in wrapper
    assert "--eligibility-exclusions" in wrapper
    assert "--eligibility-summary" in wrapper
    assert "--expected-code-revision" in wrapper
    assert "Refusing to overwrite or resume" in wrapper
    assert "native-v3r1" in wrapper
    assert HISTORICAL_V3_ROOT not in wrapper

    cli = (
        STAGE5_ROOT / "build_teacher_ssl_fullpool_native_shard.py"
    ).read_text(encoding="utf-8")
    assert '"--expected-code-revision"' in cli
    assert "expected_git_sha=args.expected_code_revision" in cli


def test_controller_is_noninteractive_explicit_and_fail_closed() -> None:
    wrapper = CONTROLLER.read_text(encoding="utf-8")
    assert 'readonly LOCKED_ORCD_HOST="tehan@orcd-login.mit.edu"' in wrapper
    assert (
        'readonly LOCKED_CONTROL_PATH="/Users/tehan/.ssh/cm/%r@%h:%p"'
        in wrapper
    )
    assert "ORCD_HOST differs from the locked noninteractive endpoint" in wrapper
    assert "ORCD_CONTROL_PATH differs from the locked user-opened socket" in wrapper
    for option in (
        "BatchMode=yes",
        "PasswordAuthentication=no",
        "KbdInteractiveAuthentication=no",
        "PubkeyAuthentication=no",
        "HostbasedAuthentication=no",
        "GSSAPIAuthentication=no",
        "NumberOfPasswordPrompts=0",
        "ControlMaster=no",
    ):
        assert option in wrapper
    for action in (
        "stage-eligibility)",
        "submit-native-canaries)",
        "submit-native-canary-audit)",
        "submit-native)",
        "submit-native-registry)",
        "submit-ssl-registry)",
        "submit-numeric-audit)",
        "submit-numeric-release)",
        "submit-smoke)",
        "submit-folds)",
        "validate-folds)",
        "verify-completion)",
    ):
        assert action in wrapper
    assert V3R1_ROOT in wrapper
    assert V3R1_RUN_ROOT_ENV in wrapper
    assert MODEL_NAMESPACE in wrapper
    assert NUMERIC_GATE_NAMESPACE in wrapper
    assert "sbatch --dependency" not in wrapper


def test_controller_locks_fresh_v3r1_identities_and_historical_v3_only_as_history(
) -> None:
    wrapper = CONTROLLER.read_text(encoding="utf-8")

    for value in (
        V3R1_ROOT,
        V3R1_RUN_ROOT_ENV,
        MODEL_NAMESPACE,
        NUMERIC_GATE_NAMESPACE,
        COMPLETION_RELEASE,
        NATIVE_CONTRACT,
        RELEASE_BINDING,
        BUILDER_CONTRACT,
        DETREND_CONTRACT,
        DETREND_TIME_CONTRACT,
        DETREND_TIME_DATASET,
        WARNING_CAPTURE_POLICY,
        EMPTY_WARNING_LEDGER_SHA256,
    ):
        assert value in wrapper

    assert (
        f'readonly HISTORICAL_V3_ROOT="${{TWIRL_ROOT}}/'
        f"reports/stage5_validation/{HISTORICAL_V3_ROOT}\""
        in wrapper
    )
    assert (
        'readonly HISTORICAL_V3_NATIVE_CONTRACT="'
        'twirl_teacher_ssl_fullpool_real_native_v3_effective_quality_adp"'
        in wrapper
    )
    assert (
        'readonly HISTORICAL_V3_RELEASE_BINDING="'
        f'{HISTORICAL_V3_ROOT}"'
        in wrapper
    )
    assert (
        f'readonly RUN_ROOT="${{TWIRL_ROOT}}/reports/stage5_validation/'
        f'{V3R1_ROOT}"'
        in wrapper
    )
    assert (
        f'readonly CLAIM_PREFIX="${{TWIRL_ROOT}}/reports/stage5_validation/.'
        f'{V3R1_ROOT}"'
        in wrapper
    )
    assert (
        f",TWIRL_SSL_FULLPOOL_V3R1_RUN_ROOT=${{RUN_ROOT}}"
        in wrapper
    )
    assert ",TWIRL_SSL_FULLPOOL_V3_RUN_ROOT=" not in wrapper


@pytest.mark.parametrize(
    ("asset", "job_name", "event_label"),
    (
        (
            CANARY_AUDIT,
            "twirl-ssl-v3r1-canary-numeric",
            "[teacher-ssl-fullpool-v3r1-canary-numeric]",
        ),
        (
            NATIVE,
            "twirl-ssl-v3r1-native",
            "[fullpool-native-v3r1]",
        ),
        (
            NATIVE_REGISTRY,
            "twirl-ssl-v3r1-native-registry",
            "[ssl-full-pool-native-registry-v3r1]",
        ),
        (
            SSL_REGISTRY,
            "twirl-ssl-v3r1-registry",
            "[ssl-full-pool-registry-v3r1]",
        ),
        (
            NUMERIC_AUDIT,
            "twirl-ssl-v3r1-numeric-audit",
            "[teacher-ssl-numeric-audit-v3r1]",
        ),
        (
            NUMERIC_RELEASE,
            "twirl-ssl-v3r1-numeric-release",
            "[teacher-ssl-numeric-release-v3r1]",
        ),
        (
            SMOKE,
            "twirl-ssl-v3r1-smoke",
            "[teacher-ssl-fullpool-v3r1-smoke]",
        ),
        (
            FOLDS,
            "twirl-ssl-v3r1-fold",
            "[teacher-ssl-fullpool-v3r1]",
        ),
        (
            FOLD_VALIDATION,
            "twirl-ssl-v3r1-validate-folds",
            "[teacher-ssl-fullpool-v3r1-fold-validation]",
        ),
    ),
)
def test_runtime_job_and_event_labels_are_unambiguously_v3r1(
    asset: Path,
    job_name: str,
    event_label: str,
) -> None:
    wrapper = asset.read_text(encoding="utf-8")
    assert f"#SBATCH --job-name={job_name}" in wrapper
    assert wrapper.count(event_label) >= 1
    assert "#SBATCH --job-name=twirl-ssl-v3-" not in wrapper


def _controller_action_block(action: str) -> str:
    wrapper = CONTROLLER.read_text(encoding="utf-8")
    start_token = f"  {action})\n"
    start = wrapper.index(start_token) + len(start_token)
    end = wrapper.find("\n    ;;", start)
    assert end > start
    return wrapper[start:end]


@pytest.mark.parametrize(
    ("action", "first_mutation"),
    (
        ("deploy", "git -C '${REMOTE_SOURCE}' fetch"),
        ("stage-eligibility", 'mkdir -p \\"\\${staging}'),
        ("submit-native-canaries", "mkdir -p '${LAUNCH_DIR}'"),
        ("submit-native-canary-audit", "sbatch --parsable"),
        ("submit-native", ": > '${NATIVE_JOBS}.submitting'"),
        ("submit-native-registry", "sbatch --parsable"),
        ("submit-ssl-registry", "sbatch --parsable"),
        ("submit-numeric-audit", "sbatch --parsable"),
        ("submit-numeric-release", "sbatch --parsable"),
        ("submit-smoke", "sbatch --parsable"),
        ("submit-folds", "sbatch --parsable"),
        ("validate-folds", "sbatch --parsable"),
    ),
)
def test_each_mutating_action_atomically_claims_before_mutation(
    action: str,
    first_mutation: str,
) -> None:
    block = _controller_action_block(action)
    claim = f"mkdir -- '${{CLAIM_PREFIX}}.{action}.claim'"
    assert block.count(claim) == 1
    assert block.index(claim) < block.index(first_mutation)


def test_stage_claims_are_unique_persistent_and_atomic(
    tmp_path: Path,
) -> None:
    wrapper = CONTROLLER.read_text(encoding="utf-8")
    assert "rm " + "'${CLAIM_PREFIX}" not in wrapper
    assert "rmdir " + "'${CLAIM_PREFIX}" not in wrapper

    claim = tmp_path / "submit-native.claim"
    first = subprocess.Popen(
        ["mkdir", str(claim)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    second = subprocess.Popen(
        ["mkdir", str(claim)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    first.wait()
    second.wait()
    assert sorted((first.returncode, second.returncode)) == [0, 1]
    assert claim.is_dir()


def test_controller_gates_canaries_and_all_native_products() -> None:
    wrapper = CONTROLLER.read_text(encoding="utf-8")
    assert "TWIRL_FULLPOOL_SECTOR=56" in wrapper
    assert "TWIRL_FULLPOOL_SHARD_INDEX=3" in wrapper
    assert "TWIRL_FULLPOOL_SECTOR=60" in wrapper
    assert "TWIRL_FULLPOOL_SHARD_INDEX=4" in wrapper
    assert 'readonly S60_CANARY_TIC="704538814"' in wrapper
    assert "targets/0000000704538814" in wrapper
    assert "native_4.model_facing_numeric.json" in wrapper
    assert "native_canary_numeric.job" in wrapper
    assert "slurm_teacher_ssl_fullpool_v3_canary_numeric_cpu.sbatch" in wrapper
    assert 'd[\\"passed\\"] is True and d[\\"n_failures\\"]==0' in wrapper
    assert 'c[\\"failed_observations\\"]==0' in wrapper
    assert "--array=0-15%4" in wrapper
    assert "total_tasks" in wrapper
    assert "[[ \\${total_tasks} -eq 112 ]]" in wrapper
    assert "[[ \\${products} -eq 112 ]]" in wrapper
    assert "'COMPLETED|0:0'" in wrapper
    assert "test ! -s" in wrapper
    assert "sha256sum -c" in wrapper
    for value in (
        NATIVE_CONTRACT,
        RELEASE_BINDING,
        BUILDER_CONTRACT,
        DETREND_CONTRACT,
        DETREND_TIME_CONTRACT,
        DETREND_TIME_DATASET,
        WARNING_CAPTURE_POLICY,
        EMPTY_WARNING_LEDGER_SHA256,
    ):
        assert value in wrapper
    assert 'd[\\"detrend_time_contract_version\\"]' in wrapper
    assert 'd[\\"detrend_time_dataset\\"]' in wrapper
    assert 'd[\\"detrend_time_system\\"]==\\"BTJD\\"' in wrapper
    assert 'd[\\"published_time_system\\"]==\\"BJD\\"' in wrapper
    assert 'd[\\"rank_warning_publication_policy\\"]==\\"require_zero\\"' in wrapper
    assert 'int(d[\\"rank_warning_count\\"])==0' in wrapper
    assert 'd[\\"rank_warning_ledger_json\\"]==\\"[]\\"' in wrapper
    assert 'd[\\"rank_warning_ledger_sha256\\"]' in wrapper
    assert 'h.attrs[\\"detrend_time_contract_version\\"]' in wrapper
    assert 'assert \\"${DETREND_TIME_DATASET}\\" in g' in wrapper
    assert "np.isfinite(btjd).all()" in wrapper
    assert "btjd.min()>=0.0" in wrapper
    assert "btjd.max()<100000.0" in wrapper
    assert "np.array_equal(bjd,btjd+2457000.0)" in wrapper
    assert (
        'int(g.attrs[\\"effective_quality_adp_small_'
        'rank_warning_count\\"])==0'
        in wrapper
    )
    assert (
        'int(g.attrs[\\"effective_quality_adp_primary_'
        'rank_warning_count\\"])==0'
        in wrapper
    )
    assert "sum(int(d[\\\"n_source_shard_observations\\\"]) for d in ds)==175366" in wrapper
    assert "sum(int(d[\\\"n_shard_observations\\\"]) for d in ds)==175347" in wrapper
    assert "sum(int(d[\\\"n_shard_excluded_observations\\\"]) for d in ds)==19" in wrapper


def test_controller_stages_eligibility_authority_byte_identically() -> None:
    wrapper = CONTROLLER.read_text(encoding="utf-8")
    assert "native_model_exclusions.csv" in wrapper
    assert "native_model_eligibility.summary.json" in wrapper
    assert "cp --" in wrapper
    assert wrapper.count("cmp -s --") >= 4
    assert "eligible_excluded_union_equals_full" in wrapper
    assert "8e9e9c12a24d5ebc7be94b03a4e35411cd10066d62a87d921a8443b06cc188d1" in wrapper
    assert "6ddc8e57bb5fb938ce05389c1629221c73e0e73ac3bf40da47a2019e1a5660e6" in wrapper
    assert "b9f536144265e54a70bff17c782e0668ddbd96efcdf6c223ebc58f46edb7d976" in wrapper


def test_numeric_smoke_and_fold_contracts_are_exact() -> None:
    numeric = NUMERIC_AUDIT.read_text(encoding="utf-8")
    release = NUMERIC_RELEASE.read_text(encoding="utf-8")
    smoke = SMOKE.read_text(encoding="utf-8")
    folds = FOLDS.read_text(encoding="utf-8")

    assert "#SBATCH --array=0-111%4" in numeric
    assert NUMERIC_GATE_NAMESPACE in numeric
    assert "model_input_numeric_gate_v2" not in numeric
    assert '"${#reports[@]}" -ne 112' in release
    assert NUMERIC_GATE_NAMESPACE in release
    assert "#SBATCH --gres=gpu:h200:1" in smoke
    assert "#SBATCH --array" not in smoke
    assert "readonly FOLD=2" in smoke
    assert "readonly MAX_ROWS=4096" in smoke
    assert "readonly EPOCHS=1" in smoke
    assert 'REQUIRED_OBSERVATION_ID="s0060-tic0000000722078603"' in smoke
    assert MODEL_NAMESPACE in smoke
    assert "validate_teacher_ssl_fullpool_v3_smoke.py" in smoke
    assert "#SBATCH --gres=gpu:h200:1" in folds
    assert "#SBATCH --array=0-4%4" in folds
    assert MODEL_NAMESPACE in folds
    assert "--resume" not in folds


def test_five_fold_completion_validation_is_independent_and_immutable() -> None:
    wrapper = FOLD_VALIDATION.read_text(encoding="utf-8")
    controller = CONTROLLER.read_text(encoding="utf-8")

    assert "#SBATCH --gres" not in wrapper
    assert "#SBATCH --array" not in wrapper
    assert V3R1_ROOT in wrapper
    assert V3R1_RUN_ROOT_ENV in wrapper
    assert MODEL_NAMESPACE in wrapper
    assert "training/five_fold" in wrapper
    assert "validate_teacher_ssl_fullpool_v3_training.py" in wrapper
    assert '--model-root "${MODEL_ROOT}"' in wrapper
    assert '--expected-code-revision "${EXPECTED_SHA}"' in wrapper
    assert '--output-release "${OUTPUT_RELEASE}"' in wrapper
    assert (
        'OUTPUT_RELEASE="${RUN_ROOT}/frozen/model/'
        f'{COMPLETION_RELEASE}"'
        in wrapper
    )
    assert "Refusing to overwrite or reuse native-v3r1 completion release" in wrapper
    assert "completion.release.json" not in wrapper
    assert (
        "twirl_teacher_ssl_fullpool_v3r1_five_fold_completion_release_v2"
        in wrapper
    )
    assert (
        "teacher_ssl_fullpool_v3r1_effective_quality_adp_btjd_"
        "five_fold_complete_v2"
        in wrapper
    )
    assert '"counts"]=={"folds":5,"completed_epochs":100}' in wrapper

    action = _controller_action_block("validate-folds")
    assert "[[ \\${#fold_states[@]} -eq 5 ]]" in action
    assert "for fold in {0..4}; do" in action
    assert "'COMPLETED|0:0'" in action
    assert "[[ \\${#errs[@]} -eq 5 ]]" in action
    assert "for err in" in action and "test ! -s" in action
    assert "fold_validation.job" in controller
    assert "slurm_teacher_ssl_fullpool_v3_validate_folds_cpu.sbatch" in action
    assert (
        f"frozen/model/{COMPLETION_RELEASE}"
        in controller
    )

    verification = _controller_action_block("verify-completion")
    assert "'COMPLETED|0:0'" in verification
    assert "twirl-ssl-v3r1-validate-folds-" in verification
    assert "test ! -s" in verification
    assert f"sha256sum -c '{COMPLETION_RELEASE}.sha256'" in verification
    assert "validate_teacher_ssl_fullpool_v3_training.py" in verification
    assert "--model-root '${MODEL_RUN_ROOT}/training/five_fold'" in verification
    assert "--output-release '${FOLD_COMPLETION_RELEASE}'" in verification
    assert '[[ \\"\\${after}\\" == \\"\\${before}\\" ]]' in verification
    assert "mkdir --" not in verification
    assert "sbatch" not in verification


def test_s60_canary_auditor_uses_only_primary_label_free_authorities() -> None:
    auditor = CANARY_AUDIT_PY.read_text(encoding="utf-8")
    wrapper = CANARY_AUDIT.read_text(encoding="utf-8")

    for api in (
        "load_frozen_ssl_full_pool",
        "load_global_full_pool_bls",
        "derive_anchor_eligibility",
        "fullpool_dataset_rows",
        "HarmonicNativeDataset",
        "collate_native_samples",
        "audit_collated_sample",
        "verify_full_pool_native_shard",
    ):
        assert api in auditor
    assert "load_fullpool_ssl_registry" not in auditor
    assert "load_full_pool_native_registry_release" not in auditor
    assert "teacher_ssl_fullpool_registry.parquet" not in auditor
    assert "native-v2" not in auditor
    assert "CANARY_TARGET_TIC = 704_538_814" in auditor
    assert '"DET_FLUX_ADP_SML"' in auditor
    assert "CANARY_SHARD_INDEX = 4" in auditor
    assert "PRODUCTION_FULL_IDENTITY_SHA256" in auditor
    assert "PRODUCTION_ELIGIBLE_IDENTITY_SHA256" in auditor
    assert "FULL_POOL_NATIVE_CONTRACT_VERSION" in auditor
    assert "FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION" in auditor
    assert "FULL_POOL_NATIVE_DETREND_TIME_DATASET" in auditor
    assert "FULL_POOL_NATIVE_BTJD_TO_BJD_OFFSET_D" in auditor
    assert "FULL_POOL_NATIVE_WARNING_CAPTURE_POLICY" in auditor
    assert "FULL_POOL_NATIVE_RANK_WARNING_PUBLICATION_POLICY" in auditor
    assert "FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_JSON" in auditor
    assert "FULL_POOL_NATIVE_EMPTY_RANK_WARNING_LEDGER_SHA256" in auditor
    assert "expected_code_revision != _code_revision()" in auditor
    assert "_publish_new(report_out, payload)" in auditor
    assert "audit_only_no_clip_no_exclusion" in auditor
    assert (
        "twirl_teacher_ssl_fullpool_v3r1_s60_shard4_"
        "model_facing_canary_v2"
        in auditor
    )
    assert '"detrend_time_system": "BTJD"' in auditor
    assert '"published_time_system": "BJD"' in auditor
    assert '"rank_warning_count": 0' in auditor
    assert "teacher_ssl_fullpool_v3r1_canary_numeric_complete" in auditor

    assert "#SBATCH --gres" not in wrapper
    assert "#SBATCH --array" not in wrapper
    assert "sector=60 shard=4 tic=704538814" in wrapper
    assert "native_4.model_facing_numeric.json" in wrapper
    assert "Refusing to overwrite or reuse" in wrapper


def test_s60_canary_over_envelope_payload_fails_closed() -> None:
    from twirl.vetting.ssl_full_pool_numeric import audit_collated_sample

    harmonic = np.zeros((1, 7, 7, 2), dtype=np.float32)
    harmonic[0, 0, 0, 0] = np.float32(5.0e36)
    batch = {
        "review_id": ["s0060-tic0000000704538814"],
        "tic": np.asarray([704_538_814], dtype=np.int64),
        "period_d": np.asarray([1.0], dtype=np.float32),
        "duration_min": np.asarray([3.0], dtype=np.float32),
        "harmonic_values": harmonic,
        "harmonic_mask": np.ones_like(harmonic, dtype=bool),
        "local_values": np.zeros((1, 14, 7, 2), dtype=np.float32),
        "local_mask": np.ones((1, 14, 7, 2), dtype=bool),
        "periodogram_values": np.zeros((1, 4, 2), dtype=np.float32),
        "periodogram_mask": np.ones((1, 4, 2), dtype=bool),
        "metadata": np.empty((1, 0), dtype=np.float32),
    }

    report = audit_collated_sample(batch)

    assert report["passed"] is False
    assert report["action"] == "audit_only_no_clip_no_exclusion"
    assert {
        str(item["code"]) for item in report["failures"]
    } >= {"numeric_envelope_exceeded", "float32_square_overflow"}


def test_s60_canary_rejects_a_drifted_native_sha_sidecar(
    tmp_path: Path,
) -> None:
    spec = importlib.util.spec_from_file_location(
        "teacher_ssl_fullpool_v3_canary_audit_test",
        CANARY_AUDIT_PY,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    native = tmp_path / "native_4.h5"
    native.write_bytes(b"current-v3r1-canary")
    sidecar = tmp_path / "native_4.h5.sha256"
    sidecar.write_text(
        f"{'0' * 64}  {native}\n",
        encoding="ascii",
    )

    with pytest.raises(ValueError, match="digest differs"):
        module._validate_sha_sidecar(
            sidecar,
            artifact_path=native,
            artifact_sha256=module.file_sha256(native),
        )


def test_controller_persists_and_gates_every_downstream_job_id() -> None:
    wrapper = CONTROLLER.read_text(encoding="utf-8")
    for manifest in (
        "native_registry.job",
        "native_canary_numeric.job",
        "ssl_registry.job",
        "numeric_audit.job",
        "numeric_release.job",
        "smoke.job",
        "folds.job",
        "fold_validation.job",
    ):
        assert manifest in wrapper
    assert wrapper.count("sacct -X -n -P") >= 5
    assert "sacct -n -P" in wrapper
    assert "errs=(" in wrapper


def test_controller_probe_defaults_to_a_portable_dry_run() -> None:
    env = os.environ.copy()
    env["TWIRL_EXPECTED_GIT_SHA"] = "a" * 40
    env["ORCD_HOST"] = "tehan@orcd-login.mit.edu"
    env["ORCD_CONTROL_PATH"] = "/Users/tehan/.ssh/cm/%r@%h:%p"
    completed = subprocess.run(
        ["bash", str(CONTROLLER), "probe"],
        cwd=REPO_ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.startswith("+ ssh ")
    assert "BatchMode=yes" in completed.stdout
    assert "tehan@orcd-login.mit.edu" in completed.stdout


@pytest.mark.parametrize(
    ("name", "value", "message"),
    (
        (
            "ORCD_HOST",
            "nobody@example.invalid",
            "ORCD_HOST differs from the locked noninteractive endpoint",
        ),
        (
            "ORCD_CONTROL_PATH",
            "/tmp/poisoned-control-socket",
            "ORCD_CONTROL_PATH differs from the locked user-opened socket",
        ),
    ),
)
def test_controller_rejects_polluted_connection_environment_before_ssh(
    name: str,
    value: str,
    message: str,
) -> None:
    env = os.environ.copy()
    env["TWIRL_EXPECTED_GIT_SHA"] = "a" * 40
    env[name] = value
    completed = subprocess.run(
        ["bash", str(CONTROLLER), "probe"],
        cwd=REPO_ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 2
    assert message in completed.stderr
    assert completed.stdout == ""
