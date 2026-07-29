from __future__ import annotations

import hashlib
import importlib.util
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Callable

import h5py
import pandas as pd
import pytest

import twirl.vetting.teacher_ssl_fullpool as FULLPOOL
from twirl.vetting.harmonic_cnn import (
    HarmonicModelConfig,
    build_harmonic_cnn,
)
from twirl.vetting.harmonic_ssl import (
    EventPreservingAugmentationConfig,
    VICRegConfig,
    vicreg_loss,
)
from twirl.vetting.ssl_full_pool_numeric import (
    MODEL_INPUT_NUMERIC_RELEASE_BINDING,
    numeric_native_freshness,
)
from twirl.vetting.teacher_ssl_fullpool import (
    FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
    FULLPOOL_SSL_CHECKPOINT_SCHEMA,
    FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
    FULLPOOL_SSL_ENCODER_NAME,
    FULLPOOL_SSL_LEARNING_RATE,
    FULLPOOL_SSL_MODEL_NAMESPACE,
    FULLPOOL_SSL_RUN_CONTRACT_SCHEMA,
    FULLPOOL_SSL_RUN_ID,
    FULLPOOL_SSL_SELECTION_SCHEMA,
    FULLPOOL_SSL_SUMMARY_SCHEMA,
    FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA,
    FULLPOOL_SSL_WEIGHT_DECAY,
)


torch = pytest.importorskip("torch")
ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts"
    / "stage5_validation"
    / "validate_teacher_ssl_fullpool_v3_training.py"
)
SPEC = importlib.util.spec_from_file_location(
    "validate_teacher_ssl_fullpool_v3_training_test",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metadata(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _real_embedding_only_ssl_step(
    model: Any,
    projector: Any,
    optimizer: Any,
    *,
    steps: int,
) -> None:
    batch_size = 2
    first_batch = {
        "harmonic_values": torch.randn(batch_size, 7, 7, 13),
        "harmonic_mask": torch.ones(
            batch_size, 7, 7, 13, dtype=torch.bool
        ),
        "local_values": torch.randn(batch_size, 14, 7, 13),
        "local_mask": torch.ones(
            batch_size, 14, 7, 13, dtype=torch.bool
        ),
        "periodogram_values": torch.randn(batch_size, 4, 13),
        "periodogram_mask": torch.ones(
            batch_size, 4, 13, dtype=torch.bool
        ),
        "metadata": torch.empty(batch_size, 0),
    }
    for step in range(steps):
        second_batch = {
            name: (
                value.clone()
                if value.dtype == torch.bool
                else value + 0.01 * torch.randn_like(value)
            )
            for name, value in first_batch.items()
        }
        optimizer.zero_grad(set_to_none=True)
        first = projector(model(first_batch)["embedding"])
        second = projector(model(second_batch)["embedding"])
        loss = vicreg_loss(first, second)
        if isinstance(loss, tuple):
            loss = loss[0]
        loss.backward()
        if step == 0:
            FULLPOOL._assert_ssl_gradient_coverage(
                model,
                projector,
                optimizer,
                fold=0,
                epoch=1,
                batch_index=0,
            )
        optimizer.step()


def _fixture(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> dict[str, Any]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    revision = "1" * 40
    freshness = numeric_native_freshness(revision)
    authority_bindings: dict[str, dict[str, Any]] = {}
    for name in (
        "eligibility_exclusions",
        "eligibility_summary",
        "native_registry",
        "native_registry_summary",
        "native_release_summary",
        "registry",
        "registry_summary",
    ):
        path = tmp_path / f"{name}.authority"
        path.write_text(f"{name}\n", encoding="ascii")
        authority_bindings[name] = _metadata(path)
    numeric_gate_path = tmp_path / "numeric_gate_release.json"
    numeric_gate_payload = {
        "release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "schema_version": "numeric-test",
        "envelope_canonical_sha256": "a" * 64,
        "counts": {"eligible_observations": 10},
        "identity_hashes": {"eligible": "b" * 64},
        "code_revision": revision,
        "native_freshness": freshness,
        "authority_bindings": {
            "ssl_registry": authority_bindings["registry"],
            "ssl_registry_summary": authority_bindings[
                "registry_summary"
            ],
            "native_registry": authority_bindings["native_registry"],
            "native_registry_summary": authority_bindings[
                "native_registry_summary"
            ],
            "native_release_summary": authority_bindings[
                "native_release_summary"
            ],
        },
        "passed": True,
    }
    _write_json(numeric_gate_path, numeric_gate_payload)
    authority_bindings["numeric_gate_release"] = _metadata(
        numeric_gate_path
    )
    numeric_summary = {
        "binding": authority_bindings["numeric_gate_release"],
        "release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "schema_version": "numeric-test",
        "envelope_canonical_sha256": "a" * 64,
        "counts": {"eligible_observations": 10},
        "identity_hashes": {"eligible": "b" * 64},
        "code_revision": revision,
        "passed": True,
    }
    training_authority = {
        "schema_version": FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA,
        "production_lock_passed": True,
        "source_provenance_verified": True,
        "native_freshness": freshness,
        "numeric_gate_release": numeric_summary,
        "authority_bindings": authority_bindings,
    }
    native_path = tmp_path / "native_v3.h5"
    with h5py.File(native_path, "w") as h5:
        h5.attrs["contract_version"] = freshness["native_contract_version"]
        h5.attrs["release_binding"] = freshness["release_binding"]
        h5.attrs["expected_git_sha"] = revision
        h5.attrs["builder_contract_version"] = freshness[
            "builder_contract_version"
        ]
        h5.attrs["builder_code_sha256"] = freshness[
            "builder_code_sha256"
        ]
        h5.attrs["detrend_contract_version"] = freshness[
            "detrend_contract_version"
        ]
        h5.attrs["detrend_config_sha256"] = freshness[
            "detrend_config_sha256"
        ]
        h5.attrs["detrend_quality_source"] = "final_effective_quality"
        h5.attrs["raw_photometry_only"] = 1
        h5.attrs["compact_adp_photometry_reused"] = 0
        h5.attrs["compact_adp_flux_reused"] = 0
        h5.attrs["periodogram_n"] = freshness["periodogram_n"]
    native_record = {
        "native_h5_path": str(native_path.resolve()),
        "native_h5_sha256": _sha256(native_path),
        "native_h5_size_bytes": native_path.stat().st_size,
        "native_h5_mtime_ns": native_path.stat().st_mtime_ns,
        "native_h5_ctime_ns": native_path.stat().st_ctime_ns,
        "native_h5_device": native_path.stat().st_dev,
        "native_h5_inode": native_path.stat().st_ino,
        "native_contract_version": freshness["native_contract_version"],
        "native_release_binding": freshness["release_binding"],
        "native_expected_git_sha": revision,
        "native_builder_contract_version": freshness[
            "builder_contract_version"
        ],
        "native_builder_code_sha256": freshness["builder_code_sha256"],
        "native_detrend_contract_version": freshness[
            "detrend_contract_version"
        ],
        "native_detrend_config_sha256": freshness[
            "detrend_config_sha256"
        ],
        "native_detrend_quality_source": "final_effective_quality",
        "raw_photometry_only": True,
        "compact_adp_photometry_reused": False,
        "compact_adp_flux_reused": False,
        "periodogram_n": freshness["periodogram_n"],
        "hash_verified_now": True,
        "root_contract_verified_now": True,
        "group_identities_verified_now": True,
        "n_selected_observations": 10,
    }
    model_config = HarmonicModelConfig(metadata_dim=0)
    model = build_harmonic_cnn(
        model_config,
        profile=FULLPOOL.FULLPOOL_SSL_PROFILE,
    )
    projector = FULLPOOL._projection_head(2 * model_config.embedding_dim)
    optimizer_named_parameters = FULLPOOL._optimizer_named_parameters(
        model,
        projector,
    )
    optimizer = torch.optim.AdamW(
        [
            {
                "params": [
                    parameter
                    for _, parameter in optimizer_named_parameters
                ],
                "decoupled_weight_decay": True,
            }
        ],
        lr=FULLPOOL_SSL_LEARNING_RATE,
        betas=FULLPOOL.FULLPOOL_SSL_ADAMW_BETAS,
        eps=FULLPOOL.FULLPOOL_SSL_ADAMW_EPS,
        weight_decay=FULLPOOL_SSL_WEIGHT_DECAY,
        amsgrad=False,
        foreach=None,
        maximize=False,
        capturable=False,
        differentiable=False,
        fused=None,
    )
    _real_embedding_only_ssl_step(
        model,
        projector,
        optimizer,
        steps=20,
    )
    optimizer_config = FULLPOOL._locked_optimizer_config()
    optimizer_contract = FULLPOOL._optimizer_checkpoint_contract(
        model,
        projector,
        optimizer,
        gradient_coverage_verified=True,
        expected_step=20,
    )
    epoch_coverage = {
        "selected_observations": 10,
        "observations_per_epoch": 10,
        "optimizer_steps_per_epoch": 1,
        "singleton_observations_skipped_per_epoch": 0,
    }

    def selection_for_fold(fold: int) -> dict[str, Any]:
        return {
            "selection_schema_version": FULLPOOL_SSL_SELECTION_SCHEMA,
            "held_out_fold": fold,
            "n_registry_observations": 10,
            "n_eligible_observations": 10,
            "n_eligible_tics": 10,
            "n_held_observations": 0,
            "n_held_tics": 0,
            "n_selected_observations": 10,
            "n_selected_tics": 10,
            "max_rows": None,
            "required_observation_ids": [],
            "n_required_observations": 0,
            "required_observations_selected": True,
            "selected_rows_sha256": f"{fold + 1:x}" * 64,
            "selected_tics_sha256": f"{fold + 6:x}" * 64,
            "tic_disjoint": {
                "held_fold_tics": True,
                "fixed_test_tics": True,
                "reserved_prospective_tics": True,
            },
        }

    monkeypatch.setattr(VALIDATOR, "PRODUCTION_FULL_OBSERVATIONS", 10)
    monkeypatch.setattr(VALIDATOR, "PRODUCTION_ELIGIBLE_OBSERVATIONS", 10)
    monkeypatch.setattr(
        VALIDATOR,
        "validate_numeric_gate_release",
        lambda path, **kwargs: json.loads(
            Path(path).read_text(encoding="utf-8")
        ),
    )
    monkeypatch.setattr(
        VALIDATOR,
        "load_fullpool_ssl_registry",
        lambda **kwargs: (pd.DataFrame({"row": range(10)}), {}),
    )
    monkeypatch.setattr(
        VALIDATOR,
        "select_fullpool_ssl_fold",
        lambda registry, held_out_fold, **kwargs: (
            registry.copy(),
            selection_for_fold(held_out_fold),
        ),
    )
    monkeypatch.setattr(
        VALIDATOR,
        "_verify_native_files",
        lambda selected, **kwargs: [native_record],
    )
    model_root = (
        tmp_path
        / "model_runs"
        / FULLPOOL_SSL_MODEL_NAMESPACE
        / "training"
        / "five_fold"
    )
    for fold in range(5):
        fold_dir = model_root / "encoder_pretraining" / f"fold_{fold}"
        fold_dir.mkdir(parents=True)
        selection = selection_for_fold(fold)
        contract = {
            "schema_version": FULLPOOL_SSL_RUN_CONTRACT_SCHEMA,
            "run_id": FULLPOOL_SSL_RUN_ID,
            "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
            "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
            "numeric_release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
            "native_freshness": freshness,
            "numeric_gate_release": numeric_summary,
            "fold": fold,
            "registry_path": authority_bindings["registry"]["path"],
            "registry_sha256": authority_bindings["registry"]["sha256"],
            "registry_summary_path": authority_bindings[
                "registry_summary"
            ]["path"],
            "registry_summary_sha256": authority_bindings[
                "registry_summary"
            ]["sha256"],
            "training_authority": training_authority,
            "selection_audit": selection,
            "native_files": [native_record],
            "native_hashes_verified_now": True,
            "native_root_contracts_verified_now": True,
            "native_group_identities_verified_now": True,
            "model_config": asdict(model_config),
            "augmentation_config": asdict(
                EventPreservingAugmentationConfig()
            ),
            "vicreg_config": asdict(VICRegConfig()),
            "projection_architecture": [
                2 * model_config.embedding_dim,
                128,
                64,
            ],
            "optimizer_config": optimizer_config,
            "epoch_coverage": epoch_coverage,
            "epochs": 20,
            "batch_size": 64,
            "workers": 8,
            "seed": FULLPOOL_SSL_DEFAULT_TRAINING_SEED + 1000 * fold,
            "learning_rate": FULLPOOL_SSL_LEARNING_RATE,
            "weight_decay": FULLPOOL_SSL_WEIGHT_DECAY,
            "checkpoint_every": 1,
            "require_cuda": True,
            "max_rows": None,
            "required_observation_ids": [],
            "labels_loaded": False,
            "fixed_test_tensors_constructed": False,
            "prospective_sector_tensors_constructed": False,
            "embedding_export": False,
            "neighbor_probe": False,
            "code_revision": revision,
        }
        contract_path = fold_dir / "run_contract.json"
        _write_json(contract_path, contract)
        history_rows = [
            {
                "epoch": epoch,
                "n_observations": 10,
                "loss": 1.0 / epoch,
            }
            for epoch in range(1, 21)
        ]
        history_path = fold_dir / "history.csv"
        history_path.write_text(
            "epoch,n_observations,loss\n"
            + "".join(
                f"{row['epoch']},{row['n_observations']},{row['loss']}\n"
                for row in history_rows
            ),
            encoding="ascii",
        )
        checkpoint = {
            "schema_version": FULLPOOL_SSL_CHECKPOINT_SCHEMA,
            "run_id": FULLPOOL_SSL_RUN_ID,
            "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
            "numeric_release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
            "native_freshness": freshness,
            "numeric_gate_release": numeric_summary,
            "fold": fold,
            "epochs": 20,
            "run_contract_sha256": _sha256(contract_path),
            "selection_audit": selection,
            "model_config": asdict(model_config),
            "projection_architecture": [
                2 * model_config.embedding_dim,
                128,
                64,
            ],
            "optimizer_config": optimizer_config,
            "optimizer_contract": optimizer_contract,
            "epoch_coverage": epoch_coverage,
            "history": history_rows,
            "encoder_state_dict": model.state_dict(),
            "projection_state_dict": projector.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
        }
        checkpoint_path = fold_dir / "encoder.pt"
        torch.save(checkpoint, checkpoint_path)
        summary = {
            "schema_version": FULLPOOL_SSL_SUMMARY_SCHEMA,
            "run_id": FULLPOOL_SSL_RUN_ID,
            "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
            "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
            "numeric_release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
            "native_freshness": freshness,
            "numeric_gate_release": numeric_summary,
            "fold": fold,
            "run_contract": str(contract_path.resolve()),
            "run_contract_sha256": _sha256(contract_path),
            "checkpoint": str(checkpoint_path.resolve()),
            "checkpoint_sha256": _sha256(checkpoint_path),
            "history": str(history_path.resolve()),
            "history_sha256": _sha256(history_path),
            "completed_epochs": 20,
            "global_step": 20,
            "epoch_coverage": epoch_coverage,
            "selection_audit": selection,
            "labels_loaded": False,
            "fixed_test_status": "host_excluded_tensors_not_constructed",
            "prospective_status": "host_excluded_tensors_not_constructed",
            "automatic_production_promotion": False,
        }
        _write_json(fold_dir / "summary.json", summary)
    return {
        "model_root": model_root,
        "revision": revision,
        "native_path": native_path,
        "numeric_gate_path": numeric_gate_path,
    }


def _rewrite_contract(
    case: dict[str, Any],
    fold: int,
    mutate: Callable[[dict[str, Any]], None],
) -> None:
    fold_dir = (
        case["model_root"] / "encoder_pretraining" / f"fold_{fold}"
    )
    contract_path = fold_dir / "run_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    mutate(contract)
    _write_json(contract_path, contract)
    summary_path = fold_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["run_contract_sha256"] = _sha256(contract_path)
    _write_json(summary_path, summary)


def test_training_validator_accepts_and_revalidates_immutable_release(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    release = tmp_path / "completion.json"
    first = VALIDATOR.write_teacher_ssl_fullpool_completion_release(
        model_root=case["model_root"],
        expected_code_revision=case["revision"],
        output_path=release,
    )
    second = VALIDATOR.write_teacher_ssl_fullpool_completion_release(
        model_root=case["model_root"],
        expected_code_revision=case["revision"],
        output_path=release,
    )
    assert first["passed"] is True
    assert first["counts"] == {"folds": 5, "completed_epochs": 100}
    assert first["completion_release"]["sha256"] == second[
        "completion_release"
    ]["sha256"]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("learning_rate", 1.0e-3, "learning rate"),
        ("weight_decay", 0.0, "weight decay"),
        ("epochs", 19, "epochs"),
        ("max_rows", 4096, "smoke/reuse"),
        ("required_observation_ids", ["diagnostic"], "smoke/reuse"),
    ],
)
def test_training_validator_rejects_hyperparameter_or_smoke_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    field: str,
    value: Any,
    message: str,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    _rewrite_contract(
        case,
        0,
        lambda contract: contract.update({field: value}),
    )
    with pytest.raises(ValueError, match=message):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


def test_training_validator_rejects_missing_or_duplicate_fold(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    fold_four = case["model_root"] / "encoder_pretraining" / "fold_4"
    fold_four.rename(fold_four.with_name("fold_5"))
    with pytest.raises(ValueError, match="exactly fold 0--4"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


def test_training_validator_rejects_incomplete_history(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    fold_dir = case["model_root"] / "encoder_pretraining" / "fold_1"
    history_path = fold_dir / "history.csv"
    lines = history_path.read_text(encoding="ascii").splitlines()
    history_path.write_text("\n".join(lines[:-1]) + "\n", encoding="ascii")
    summary_path = fold_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["history_sha256"] = _sha256(history_path)
    _write_json(summary_path, summary)
    with pytest.raises(ValueError, match="exactly 20 epochs"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


@pytest.mark.parametrize(
    ("state_name", "message"),
    [
        ("encoder_state_dict", "non-finite tensor"),
        ("optimizer_state_dict", "non-finite tensor"),
    ],
)
def test_training_validator_rejects_nonfinite_model_or_optimizer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    state_name: str,
    message: str,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    fold_dir = case["model_root"] / "encoder_pretraining" / "fold_3"
    checkpoint_path = fold_dir / "encoder.pt"
    checkpoint = torch.load(
        checkpoint_path, map_location="cpu", weights_only=False
    )
    if state_name == "encoder_state_dict":
        first_name = next(iter(checkpoint[state_name]))
        checkpoint[state_name][first_name].reshape(-1)[0] = float("inf")
    else:
        checkpoint[state_name]["state"][0]["exp_avg"][0, 0] = float("nan")
    torch.save(checkpoint, checkpoint_path)
    summary_path = fold_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["checkpoint_sha256"] = _sha256(checkpoint_path)
    _write_json(summary_path, summary)
    with pytest.raises(ValueError, match=message):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("lr", 1.0e-3, "optimizer lr"),
        ("weight_decay", 0.0, "optimizer weight decay"),
    ],
)
def test_training_validator_rejects_optimizer_hyperparameter_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    field: str,
    value: float,
    message: str,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    fold_dir = case["model_root"] / "encoder_pretraining" / "fold_2"
    checkpoint_path = fold_dir / "encoder.pt"
    checkpoint = torch.load(
        checkpoint_path, map_location="cpu", weights_only=False
    )
    checkpoint["optimizer_state_dict"]["param_groups"][0][field] = value
    torch.save(checkpoint, checkpoint_path)
    summary_path = fold_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["checkpoint_sha256"] = _sha256(checkpoint_path)
    _write_json(summary_path, summary)
    with pytest.raises(ValueError, match=message):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


def test_training_validator_binds_numeric_projection_to_release_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    numeric_path = case["numeric_gate_path"]
    numeric = json.loads(numeric_path.read_text(encoding="utf-8"))
    numeric["counts"]["eligible_observations"] = 11
    _write_json(numeric_path, numeric)
    binding = _metadata(numeric_path)
    _rewrite_contract(
        case,
        0,
        lambda contract: contract["training_authority"][
            "authority_bindings"
        ].update({"numeric_gate_release": binding}),
    )
    with pytest.raises(ValueError, match="projection differs"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


def test_training_validator_calls_production_numeric_release_validator(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)

    def reject_fabricated_release(*args: Any, **kwargs: Any) -> dict[str, Any]:
        raise ValueError("production numeric release validation was called")

    monkeypatch.setattr(
        VALIDATOR,
        "validate_numeric_gate_release",
        reject_fabricated_release,
    )
    with pytest.raises(ValueError, match="production numeric release"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda contract: contract["selection_audit"].update(
                {"tic_disjoint": {}}
            ),
            "selection is not full production",
        ),
        (
            lambda contract: contract["selection_audit"].update(
                {"selected_rows_sha256": "f" * 64}
            ),
            "hash-bound registry",
        ),
        (
            lambda contract: contract.update({"native_files": []}),
            "lacks selected native files",
        ),
    ],
)
def test_training_validator_rejects_fabricated_selection_or_native_set(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: Callable[[dict[str, Any]], None],
    message: str,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    _rewrite_contract(case, 0, mutation)
    with pytest.raises(ValueError, match=message):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


def test_training_validator_rejects_epoch_observation_count_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    fold_dir = case["model_root"] / "encoder_pretraining" / "fold_0"
    history_path = fold_dir / "history.csv"
    history = history_path.read_text(encoding="ascii")
    history_path.write_text(
        history.replace("1,10,", "1,9,", 1),
        encoding="ascii",
    )
    summary_path = fold_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["history_sha256"] = _sha256(history_path)
    _write_json(summary_path, summary)
    with pytest.raises(ValueError, match="epoch/count sequence"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


def test_training_validator_rejects_optimizer_state_coverage_gap(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    fold_dir = case["model_root"] / "encoder_pretraining" / "fold_4"
    checkpoint_path = fold_dir / "encoder.pt"
    checkpoint = torch.load(
        checkpoint_path, map_location="cpu", weights_only=False
    )
    first_parameter_id = checkpoint["optimizer_state_dict"][
        "param_groups"
    ][0]["params"][0]
    del checkpoint["optimizer_state_dict"]["state"][first_parameter_id]
    torch.save(checkpoint, checkpoint_path)
    summary_path = fold_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["checkpoint_sha256"] = _sha256(checkpoint_path)
    _write_json(summary_path, summary)
    with pytest.raises(ValueError, match="optimizer state coverage"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )


def test_training_validator_rejects_stale_or_tampered_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _fixture(tmp_path, monkeypatch)
    summary_path = (
        case["model_root"]
        / "encoder_pretraining"
        / "fold_0"
        / "summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["native_freshness"]["expected_git_sha"] = "2" * 40
    _write_json(summary_path, summary)
    with pytest.raises(ValueError, match="summary authority"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )

    case = _fixture(tmp_path / "native_tamper", monkeypatch)
    with case["native_path"].open("ab") as handle:
        handle.write(b"tamper")
    with pytest.raises(ValueError, match="HDF5 binding changed"):
        VALIDATOR.validate_teacher_ssl_fullpool_training(
            model_root=case["model_root"],
            expected_code_revision=case["revision"],
        )
