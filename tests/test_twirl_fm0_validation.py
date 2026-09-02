from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from twirl.models.fm0 import validation as validation_module
from twirl.models.fm0.composite_release import (
    COMPOSITE_RELEASE_SCHEMA_VERSION,
    COMPOSITE_RELEASE_STATE,
)
from twirl.models.fm0.training import CADENCE_VICREG_OBJECTIVE_IDENTITY
from twirl.models.fm0.validation import (
    FM0_3_COMPOSITE_DATASET_KIND,
    FM0_3_TRAIN_ROLE,
    REAL_RUN_CONTRACT_SCHEMA_VERSION,
    REAL_RUN_SUMMARY_SCHEMA_VERSION,
    _validate_real_checkpoint_dataset_binding,
    sha256_file,
    validate_fm0_3_prestart_smoke,
    validate_real_run_release,
    write_json_with_sha256,
    write_sha256_sidecar,
)


def _composite_binding(root: Path) -> dict[str, object]:
    return {
        "kind": FM0_3_COMPOSITE_DATASET_KIND,
        "release_root": str(root.resolve()),
        "receipt_path": str((root / "receipt.json").resolve()),
        "receipt_sha256": sha256_file(root / "receipt.json"),
        "source_bindings_path": str((root / "source_bindings.csv").resolve()),
        "source_bindings_sha256": sha256_file(root / "source_bindings.csv"),
        "role_index_path": str((root / "role_index.csv").resolve()),
        "role_index_sha256": sha256_file(root / "role_index.csv"),
        "role": FM0_3_TRAIN_ROLE,
        "n_sources": 3,
        "n_observations": 5,
        "n_excluded_missing_required_views": 1,
    }


def _minimal_composite(root: Path) -> dict[str, object]:
    root.mkdir()
    (root / "source_bindings.csv").write_text("source\n", encoding="utf-8")
    (root / "role_index.csv").write_text("role\n", encoding="utf-8")
    (root / "receipt.json").write_text(
        json.dumps(
            {
                "schema_version": COMPOSITE_RELEASE_SCHEMA_VERSION,
                "release_state": COMPOSITE_RELEASE_STATE,
                "passed": True,
                "sources": {
                    "source_bindings_sha256": sha256_file(root / "source_bindings.csv")
                },
                "selection": {
                    "role_index_sha256": sha256_file(root / "role_index.csv")
                },
                "limits": {
                    "identity_only": True,
                    "source_shards_opened": False,
                    "shards_rewritten": False,
                    "development_rows_selected": 0,
                    "scientific_training_eligible": True,
                    "sealed_rows_selected": 0,
                    "model_training_authorized": False,
                    "real_cli_training_enabled": False,
                    "sealed_test_access_authorized": False,
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "READY").write_text(
        sha256_file(root / "receipt.json") + "\n", encoding="utf-8"
    )
    return _composite_binding(root)


def _minimal_fm03_smoke(root: Path) -> tuple[dict[str, object], str]:
    root.mkdir()
    authorities = {
        "design_sha256": "a" * 64,
        "config_sha256": "b" * 64,
        "freeze_receipt_sha256": "c" * 64,
    }
    contract = {
        "campaign_id": "twirl_fm0_3_matched_canary_test",
        "variant": "TWIRL-FM0.3.1",
        "architecture": "tcn",
        "objective": CADENCE_VICREG_OBJECTIVE_IDENTITY,
        "expected_git_sha": "d" * 40,
        "authorities": authorities,
        "synthetic_smoke": True,
        "real_data_consumed": False,
        "precision": "fp32",
        "device_request": "cuda",
        "training_horizon_step": 20_000,
        "synthetic_smoke_step": 8,
        "authorized_stop_after_steps": [64, 2_000],
        "stop_after_step_is_execution_state_not_scientific_contract": True,
        "evaluation_plan": None,
        "input_release": None,
    }
    write_json_with_sha256(root / "run_contract.json", contract)
    (root / "checkpoint.pt").write_bytes(b"trusted-smoke-checkpoint")
    write_sha256_sidecar(root / "checkpoint.pt")
    summary_sha = write_json_with_sha256(
        root / "summary.json",
        {
            "passed": True,
            "synthetic_only": True,
            "real_data_consumed": False,
            "scientific_result": False,
            "variant": "TWIRL-FM0.3.1",
            "architecture": "tcn",
            "global_step": 8,
            "requested_stop_after_step": 8,
            "precision": "fp32",
        },
    )
    return authorities, summary_sha


def test_fm03_prestart_smoke_binds_exact_same_variant_release(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "smoke"
    authorities, summary_sha = _minimal_fm03_smoke(root)
    artifacts = {
        name: sha256_file(root / name)
        for name in ("run_contract.json", "checkpoint.pt", "summary.json")
    }
    monkeypatch.setattr(
        validation_module,
        "validate_run_release",
        lambda *_args, **_kwargs: {
            "passed": True,
            "checkpoint_inspected": True,
            "variant": "TWIRL-FM0.3.1",
            "architecture": "tcn",
            "global_step": 8,
            "artifact_sha256": artifacts,
        },
    )

    binding = validate_fm0_3_prestart_smoke(
        root,
        expected_summary_sha256=summary_sha,
        expected_campaign_id="twirl_fm0_3_matched_canary_test",
        expected_variant="TWIRL-FM0.3.1",
        expected_architecture="tcn",
        expected_git_sha="d" * 40,
        expected_authorities=authorities,
    )

    assert binding["kind"] == "fm0_3_same_variant_prestart_smoke"
    assert binding["summary_sha256"] == summary_sha
    assert binding["checkpoint_sha256"] == artifacts["checkpoint.pt"]
    with pytest.raises(ValueError, match="matched canary"):
        validate_fm0_3_prestart_smoke(
            root,
            expected_summary_sha256=summary_sha,
            expected_campaign_id="twirl_fm0_3_matched_canary_test",
            expected_variant="TWIRL-FM0.3.2",
            expected_architecture="conformer",
            expected_git_sha="d" * 40,
            expected_authorities=authorities,
        )


def test_real_checkpoint_requires_exact_fm03_composite_dataset_binding() -> None:
    release = {
        "kind": FM0_3_COMPOSITE_DATASET_KIND,
        "release_root": "/composite",
        "receipt_sha256": "a" * 64,
        "source_bindings_sha256": "b" * 64,
        "role_index_sha256": "c" * 64,
        "role": FM0_3_TRAIN_ROLE,
        "n_sources": 3,
        "n_observations": 5,
        "n_excluded_missing_required_views": 1,
    }
    contract = {
        "seed": 560067,
        "input_release": release,
        "optimization": {
            "max_optimizer_steps": 20_000,
            "effective_batch_windows": 64,
        },
    }
    dataset = {
        "kind": FM0_3_COMPOSITE_DATASET_KIND,
        "composite_root": "/composite",
        "receipt_sha256": "a" * 64,
        "source_bindings_sha256": "b" * 64,
        "role_index_sha256": "c" * 64,
        "variant": "TWIRL-FM0.3.1",
        "role": FM0_3_TRAIN_ROLE,
        "seed": 560067,
        "windows_per_epoch": 1_280_000,
        "n_sources": 3,
        "n_observations": 5,
        "n_excluded_missing_required_views": 1,
    }

    _validate_real_checkpoint_dataset_binding(
        contract=contract,
        dataset_contract=dataset,
        variant="TWIRL-FM0.3.1",
        cadence_objective=True,
    )
    drifted = dict(dataset, role="temporal_holdout")
    with pytest.raises(ValueError, match="identity differs"):
        _validate_real_checkpoint_dataset_binding(
            contract=contract,
            dataset_contract=drifted,
            variant="TWIRL-FM0.3.1",
            cadence_objective=True,
        )


def test_real_release_validator_rechecks_fm03_composite_hashes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    release = _minimal_composite(tmp_path / "composite")
    deep_calls: list[dict[str, object]] = []

    def fake_deep_validate(root: Path, **kwargs: object) -> object:
        deep_calls.append({"root": root, **kwargs})
        return object()

    monkeypatch.setattr(
        "twirl.models.fm0.composite_release.validate_composite_release",
        fake_deep_validate,
    )
    evaluation_root = tmp_path / "evaluation-plan"
    evaluation_receipt = {
        "producer_git_sha": "d" * 40,
        "identity_plan": {
            "receipt_sha256": "e" * 64,
            "schedule_sha256": "f" * 64,
            "producer_git_sha": "d" * 40,
            "input_authorities": {
                "temporal_panel": {
                    "receipt_sha256": (
                        "78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
                    ),
                    "panel_sha256": "8" * 64,
                    "sector_bindings_sha256": "9" * 64,
                },
                "composite_release": {
                    "receipt_sha256": release["receipt_sha256"],
                    "source_bindings_sha256": release[
                        "source_bindings_sha256"
                    ],
                    "role_index_sha256": release["role_index_sha256"],
                }
            },
        },
        "payload_bindings": {
            "source_shard_bindings_sha256": "1" * 64,
            "crop_payload_bindings_sha256": "2" * 64,
            "n_crops_frozen": 1_440,
        },
    }
    evaluation_result = SimpleNamespace(
        root=evaluation_root.resolve(),
        receipt_path=(evaluation_root / "receipt.json").resolve(),
        receipt_sha256="3" * 64,
        schedule_path=(evaluation_root / "schedule.csv").resolve(),
        schedule_sha256="4" * 64,
        receipt=evaluation_receipt,
    )
    evaluation_calls: list[dict[str, object]] = []

    def fake_evaluation_validate(root: Path, **kwargs: object) -> object:
        evaluation_calls.append({"root": root, **kwargs})
        return evaluation_result

    monkeypatch.setattr(
        "twirl.models.fm0.matched_canary_payload_plan."
        "validate_matched_canary_payload_plan",
        fake_evaluation_validate,
    )
    smoke_binding = {
        "kind": "fm0_3_same_variant_prestart_smoke",
        "root": str((tmp_path / "smoke").resolve()),
        "summary_sha256": "5" * 64,
    }
    smoke_calls: list[dict[str, object]] = []

    def fake_smoke_validate(root: Path, **kwargs: object) -> dict[str, object]:
        smoke_calls.append({"root": root, **kwargs})
        return smoke_binding

    monkeypatch.setattr(
        validation_module,
        "validate_fm0_3_prestart_smoke",
        fake_smoke_validate,
    )
    evaluation_binding = {
        "kind": "fm0_3_payload_screened_evaluation_plan",
        "root": str(evaluation_result.root),
        "receipt_path": str(evaluation_result.receipt_path),
        "receipt_sha256": evaluation_result.receipt_sha256,
        "schedule_path": str(evaluation_result.schedule_path),
        "schedule_sha256": evaluation_result.schedule_sha256,
        "producer_git_sha": "d" * 40,
        "identity_plan_receipt_sha256": "e" * 64,
        "identity_plan_schedule_sha256": "f" * 64,
        "identity_plan_producer_git_sha": "d" * 40,
        "temporal_panel_receipt_sha256": (
            "78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
        ),
        "temporal_panel_sha256": "8" * 64,
        "temporal_panel_sector_bindings_sha256": "9" * 64,
        "source_shard_bindings_sha256": "1" * 64,
        "crop_payload_bindings_sha256": "2" * 64,
        "n_crops": 1_440,
    }
    run = tmp_path / "run"
    run.mkdir()
    contract = {
        "schema_version": REAL_RUN_CONTRACT_SCHEMA_VERSION,
        "campaign_id": "twirl_fm0_3_matched_canary_test",
        "variant": "TWIRL-FM0.3.1",
        "architecture": "tcn",
        "objective": CADENCE_VICREG_OBJECTIVE_IDENTITY,
        "synthetic_smoke": False,
        "real_data_consumed": True,
        "expected_git_sha": "d" * 40,
        "authorities": {},
        "prestart_smoke": smoke_binding,
        "input_release": release,
        "input_release_reuse": None,
        "evaluation_plan": evaluation_binding,
        "immutable_milestone_steps": [],
    }
    contract_sha = write_json_with_sha256(run / "run_contract.json", contract)
    (run / "checkpoint.pt").write_bytes(b"trusted-test-checkpoint")
    checkpoint_sha = write_sha256_sidecar(run / "checkpoint.pt")
    write_json_with_sha256(
        run / "summary.json",
        {
            "schema_version": REAL_RUN_SUMMARY_SCHEMA_VERSION,
            "passed": True,
            "synthetic_only": False,
            "real_data_consumed": True,
            "variant": "TWIRL-FM0.3.1",
            "architecture": "tcn",
            "global_step": 64,
            "run_contract_sha256": contract_sha,
            "checkpoint_sha256": checkpoint_sha,
            "immutable_milestone_checkpoints": [],
        },
    )

    result = validate_real_run_release(run, inspect_checkpoint=False)

    assert result["passed"] is True
    assert smoke_calls == [
        {
            "root": str((tmp_path / "smoke").resolve()),
            "expected_summary_sha256": "5" * 64,
            "expected_campaign_id": "twirl_fm0_3_matched_canary_test",
            "expected_variant": "TWIRL-FM0.3.1",
            "expected_architecture": "tcn",
            "expected_git_sha": "d" * 40,
            "expected_authorities": {},
            "inspect_checkpoint": False,
        }
    ]
    assert deep_calls == [
        {
            "root": (tmp_path / "composite").resolve(),
            "expected_receipt_sha256": release["receipt_sha256"],
            "expected_source_bindings_sha256": release[
                "source_bindings_sha256"
            ],
            "expected_role_index_sha256": release["role_index_sha256"],
            "require_read_only": True,
            "require_source_read_only": True,
        }
    ]
    assert evaluation_calls == [
        {
            "root": str(evaluation_result.root),
            "expected_receipt_sha256": "3" * 64,
            "require_read_only": True,
        }
    ]
    (tmp_path / "composite" / "role_index.csv").write_text("drift\n", encoding="utf-8")
    with pytest.raises(ValueError, match="changed after training"):
        validate_real_run_release(run, inspect_checkpoint=False)
