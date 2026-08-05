from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
from types import ModuleType

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = (
    REPO_ROOT
    / "configs/models/teacher_ssl_v1_s56_s62_a2v1_current_adp.yaml"
)
CLI_PATH = (
    REPO_ROOT / "scripts/stage5_validation/train_teacher_ssl_v1.py"
)
SBATCH_PATH = (
    REPO_ROOT / "scripts/orcd/slurm_teacher_ssl_v1_h200.sbatch"
)


def _load_cli():
    spec = importlib.util.spec_from_file_location(
        "train_teacher_ssl_v1_contract_test",
        CLI_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_teacher_ssl_config_freezes_no_leakage_contract() -> None:
    cli = _load_cli()
    payload = cli._load_config(CONFIG_PATH)

    assert payload["model_facing_name"] == "Teacher v4-SSL"
    assert payload["encoder_name"] == "teacher_ssl_v1"
    assert payload["baseline"]["expected_summary_sha256"] == (
        "0bdcb064a7e67f2304ba58b1e79c462d"
        "aa7cc8aad5ad64d2f57a86c4dae46e99"
    )
    assert payload["scope"]["sectors"] == list(range(56, 63))
    assert payload["scope"]["development_only"] is True
    assert payload["scope"]["held_cv_fold_excluded"] is True
    assert payload["scope"]["fixed_test_model_use_forbidden"] is True
    assert payload["scope"]["prospective_sector_63_model_use_forbidden"] is True
    assert payload["scope"]["injections_forbidden_from_ssl"] is True
    assert payload["unlabeled_pool"]["uncertain_included"] is True
    assert payload["unlabeled_pool"]["wide_transit_like_included"] is True
    assert payload["unlabeled_pool"]["scalar_metadata"] == "absent"
    assert payload["ssl"]["profile"] == "shape_plus_periodogram_bls"
    assert payload["ssl"]["objective"] == "vicreg"
    assert payload["ssl"]["event_window_protected"] is True
    assert payload["ssl"]["invariant_views"] == [
        "event_protected_mask_dropout",
        "small_valid_flux_noise",
    ]
    assert payload["ssl"]["augmentation_parameters"] == {
        "harmonic_cadence_dropout_probability": 0.02,
        "periodogram_bin_dropout_probability": 0.02,
        "local_masks_preserved": True,
        "scalar_metadata_zeroed": True,
        "invalid_samples_never_unmasked": True,
    }
    assert set(payload["ssl"]["prohibited_augmentations"]) == {
        "crop",
        "phase_warp",
        "time_warp",
        "reversal",
        "smoothing",
        "flux_inversion",
        "depth_scaling",
        "aperture_swap",
        "harmonic_reassignment",
        "event_insertion",
        "event_removal",
    }
    assert payload["optimization"]["ssl_epochs"] == 20
    assert payload["optimization"]["fine_tune_epochs"] == 100
    assert payload["optimization"]["batch_size"] == 32
    assert payload["optimization"]["seed"] == 560064
    assert payload["fine_tuning"]["label_policy"] == "uncertain_masked"
    assert payload["fine_tuning"]["held_fold_used_for_early_stopping"] is True
    assert payload["fine_tuning"]["estimate_status"] == (
        "matched_development_not_untouched"
    )
    assert payload["compute"]["h200_gpus"] == 1
    assert payload["compute"]["approved_ssl_ceiling_h200_gpus"] == 4


@pytest.mark.parametrize(
    ("section", "key", "bad_value"),
    [
        ("scope", "development_only", False),
        ("scope", "held_cv_fold_excluded", False),
        ("scope", "fixed_test_model_use_forbidden", False),
        ("scope", "prospective_sector_63_model_use_forbidden", False),
        ("scope", "injections_forbidden_from_ssl", False),
        ("unlabeled_pool", "scalar_metadata", "present"),
        ("ssl", "objective", "contrastive"),
        ("ssl", "event_window_protected", False),
        ("fine_tuning", "label_policy", "uncertain_as_other"),
        ("compute", "h200_gpus", 2),
        ("compute", "approved_ssl_ceiling_h200_gpus", 5),
    ],
)
def test_teacher_ssl_config_rejects_contract_drift(
    tmp_path: Path,
    section: str,
    key: str,
    bad_value: object,
) -> None:
    cli = _load_cli()
    payload = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    payload[section][key] = bad_value
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="invalid Teacher v4-SSL config"):
        cli._load_config(path)


def test_teacher_ssl_config_rejects_unknown_keys(tmp_path: Path) -> None:
    cli = _load_cli()
    payload = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    payload["scope"]["permit_test_tics"] = True
    path = tmp_path / "unknown.yaml"
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="unknown keys"):
        cli._load_config(path)


def test_teacher_ssl_h200_wrapper_is_bounded_and_fresh() -> None:
    completed = subprocess.run(
        ["bash", "-n", str(SBATCH_PATH)],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    wrapper = SBATCH_PATH.read_text(encoding="utf-8")
    assert "#SBATCH -p pg_mki_aryeh" in wrapper
    assert "#SBATCH -c 16" in wrapper
    assert "#SBATCH --mem=192G" in wrapper
    assert "#SBATCH --gres=gpu:h200:1" in wrapper
    assert wrapper.count("--gres=gpu:h200:1") == 1
    assert "#SBATCH -t 08:00:00" in wrapper
    assert "twirl-s56-torch/bin/python" in wrapper
    assert "reports/stage5_validation/teacher_v3_s56_s62" in wrapper
    assert 'if [[ -e "${OUT_DIR}" ]]' in wrapper
    assert "Refusing to reuse non-fresh" in wrapper
    assert "--allow-cpu" not in wrapper
    assert "--training-table" in wrapper
    assert "--native-registry" in wrapper
    assert "--native-registry-summary" in wrapper
    assert "--baseline-teacher-v3-root" in wrapper
    assert "teacher_v3_s56_s62_a2v1_current_adp/full" in wrapper


def test_teacher_ssl_cli_wires_frozen_inputs_and_overrides(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cli = _load_cli()
    captured: dict[str, object] = {}
    fake_module = ModuleType("twirl.vetting.teacher_ssl_training")

    def fake_run(**kwargs):
        captured.update(kwargs)
        return {"status": "ok"}

    fake_module.run_teacher_ssl_pilot = fake_run
    monkeypatch.setitem(
        sys.modules,
        "twirl.vetting.teacher_ssl_training",
        fake_module,
    )
    result = cli.main(
        [
            "--config",
            str(CONFIG_PATH),
            "--input-root",
            str(REPO_ROOT),
            "--out-dir",
            str(tmp_path / "pilot"),
            "--epochs",
            "3",
            "--fine-tune-epochs",
            "4",
            "--batch-size",
            "5",
            "--workers",
            "0",
            "--seed",
            "7",
            "--allow-cpu",
        ]
    )

    assert result == 0
    assert Path(captured["training_table"]).name == (
        "training_table_with_splits.csv"
    )
    assert Path(captured["native_registry"]).name == (
        "native_input_registry.csv"
    )
    assert Path(captured["native_registry_summary"]).name == (
        "native_input_registry.summary.json"
    )
    assert Path(captured["baseline_teacher_v3_root"]).name == "full"
    assert captured["expected_baseline_summary_sha256"] == (
        "0bdcb064a7e67f2304ba58b1e79c462d"
        "aa7cc8aad5ad64d2f57a86c4dae46e99"
    )
    assert captured["ssl_epochs"] == 3
    assert captured["fine_tune_epochs"] == 4
    assert captured["batch_size"] == 5
    assert captured["workers"] == 0
    assert captured["seed"] == 7
    assert captured["require_cuda"] is False


def test_teacher_ssl_cli_has_valid_python_syntax() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "py_compile",
            str(CLI_PATH),
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
