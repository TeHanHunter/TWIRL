from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _text(name: str) -> str:
    return (ROOT / "scripts" / "orcd" / name).read_text(encoding="utf-8")


def test_v4_bls_rebuilds_uniformly_from_all_112_raw_v1_shards() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_bls_cpu.sbatch")

    assert "#SBATCH --array=0-111%28" in text
    assert "N_SHARDS=16" in text
    assert "raw_source_${SHARD_INDEX}.h5" in text
    assert "--raw-source-h5" in text
    assert "--raw-source-summary" in text
    assert "--raw-export-complete" in text
    assert "--raw-transfer-validation" in text
    assert "--frozen-pool" in text
    assert "--n-periods 50000" in text
    assert "--n-peaks 10" in text
    assert "A2v1-fullpool-v1" in text
    assert "TWIRL_SSL_FULLPOOL_EXECUTION_ALLOWLIST_ROOT" in text
    assert "--execution-allowlist" in text
    builder = (
        ROOT
        / "scripts"
        / "stage5_validation"
        / "build_s56_adp_real_bls_peaks.py"
    ).read_text(encoding="utf-8")
    assert '"cadence_time_inventory_reference_only"' in builder


def test_v4_bls_merges_exactly_16_shards_per_sector() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_bls_merge_cpu.sbatch")

    assert "#SBATCH --array=0-6%7" in text
    assert "N_SHARDS=16" in text
    assert "--n-shards \"${N_SHARDS}\"" in text


def test_v4_jobs_bind_the_exact_deployed_commit() -> None:
    for name in (
        "slurm_teacher_ssl_fullpool_v4_bls_cpu.sbatch",
        "slurm_teacher_ssl_fullpool_v4_bls_merge_cpu.sbatch",
        "slurm_teacher_ssl_fullpool_v4_bls_global_cpu.sbatch",
    ):
        text = _text(name)
        assert "TWIRL_EXPECTED_GIT_SHA" in text
        assert 'git -C "${REPO}" rev-parse HEAD' in text


def test_v4_canary_is_the_exact_34_affected_observation_set() -> None:
    root = ROOT / "configs" / "teacher_ssl_fullpool_v4_canary"
    by_sector = {
        int(path.stem[1:3]): pd.read_csv(path)["tic"].astype(int).tolist()
        for path in sorted(root.glob("s*_affected.csv"))
    }

    assert set(by_sector) == {56, 58, 61, 62}
    assert {sector: len(tics) for sector, tics in by_sector.items()} == {
        56: 7,
        58: 3,
        61: 16,
        62: 8,
    }
    assert len({(sector, tic) for sector, tics in by_sector.items() for tic in tics}) == 34


def test_v4_eligibility_is_rederived_from_corrected_global_bls() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_native_eligibility_cpu.sbatch")

    assert '${RUN_ROOT}/real_bls/global' in text
    assert "native_model_eligibility_rederived" in text
    assert "--production-lock" in text
    assert "build_teacher_ssl_fullpool_native_eligibility.py" in text


def test_v4_native_rebuilds_all_112_shards_from_raw_v1_authority() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_native_cpu.sbatch")

    assert "#SBATCH --array=0-111%28" in text
    assert 'RAW_SOURCE="${RAW_DIR}/s${SECTOR}/raw_source_${SHARD_INDEX}.h5"' in text
    assert "--raw-source-h5" in text
    assert "--compact-adp-h5" in text
    assert "--n-periods 4096" in text
    native = (
        ROOT / "src" / "twirl" / "vetting" / "ssl_full_pool_native.py"
    ).read_text(encoding="utf-8")
    assert '"cadence_time_inventory_reference_only"' in native
    assert 'output.attrs["raw_detector_mapping_authority"] = 1' in native
    assert 'output.attrs["raw_orbitid_authority"] = 1' in native
    assert 'output.attrs["raw_internal_quality_authority"] = 1' in native


def test_v4_ssl_registry_uses_corrected_bls_and_locked_eligibility() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_registry_cpu.sbatch")

    assert '${RUN_ROOT}/real_bls/global' in text
    assert "--allow-rederived-eligibility" not in text
    assert "build_teacher_ssl_fullpool_registry.py" in text
    assert "teacher_ssl_fullpool_v4r2_s56_s62" in text


def test_v4_numeric_gate_covers_all_112_native_shards() -> None:
    audit = _text("slurm_teacher_ssl_fullpool_v4_numeric_audit_cpu.sbatch")
    release = _text("slurm_teacher_ssl_fullpool_v4_numeric_release_cpu.sbatch")

    assert "#SBATCH --array=0-111%4" in audit
    assert "audit_teacher_ssl_fullpool_numeric.py" in audit
    assert "model_input_numeric_gate_v4" in audit
    assert 'if [[ "${#reports[@]}" -ne 112 ]]' in release
    assert "merge_teacher_ssl_fullpool_numeric.py" in release


def test_v4_smoke_and_four_h200_fold_plan_are_locked() -> None:
    smoke = _text("slurm_teacher_ssl_fullpool_v4_smoke_h200.sbatch")
    folds = _text("slurm_teacher_ssl_fullpool_v4_fold_h200.sbatch")

    assert "#SBATCH --gres=gpu:h200:1" in smoke
    assert "readonly MAX_ROWS=4096" in smoke
    assert "readonly EPOCHS=1" in smoke
    assert "validate_teacher_ssl_fullpool_v3_smoke.py" in smoke
    assert "#SBATCH --gres=gpu:h200:1" in folds
    assert "#SBATCH --array=0-4%4" in folds
    assert "readonly EPOCHS=20" in folds
    assert "--max-rows" not in folds


def test_v4_completion_release_requires_all_five_folds() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_validate_folds_cpu.sbatch")

    assert "for fold in {0..4}" in text
    assert "validate_teacher_ssl_fullpool_v4_training.py" in text
    assert "TWIRL_VALIDATOR_GIT_SHA" in text
    assert "TWIRL_EXPECTED_TRAINING_GIT_SHA" in text
    assert "--validator-code-revision" in text
    assert '"completed_epochs":100' in text
