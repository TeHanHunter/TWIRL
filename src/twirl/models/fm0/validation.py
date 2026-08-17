"""Fail-closed authority and artifact validation for FM0.1 smokes."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Any, Mapping


RUN_CONTRACT_SCHEMA_VERSION = "twirl_fm0_1_synthetic_run_contract_v1"
RUN_SUMMARY_SCHEMA_VERSION = "twirl_fm0_1_synthetic_run_summary_v1"
RUN_VALIDATION_SCHEMA_VERSION = "twirl_fm0_1_synthetic_run_validation_v1"
REAL_RUN_CONTRACT_SCHEMA_VERSION = "twirl_fm0_1_real_run_contract_v1"
REAL_RUN_SUMMARY_SCHEMA_VERSION = "twirl_fm0_1_real_run_summary_v1"
REAL_RUN_VALIDATION_SCHEMA_VERSION = "twirl_fm0_1_real_run_validation_v1"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON authority is not a mapping: {path}")
    return value


def write_json_with_sha256(path: str | Path, payload: Mapping[str, Any]) -> str:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    encoded = (
        json.dumps(dict(payload), indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    temporary.write_bytes(encoded)
    os.replace(temporary, destination)
    digest = hashlib.sha256(encoded).hexdigest()
    sidecar = destination.with_name(destination.name + ".sha256")
    sidecar.write_text(f"{digest}  {destination.name}\n", encoding="utf-8")
    return digest


def write_sha256_sidecar(path: str | Path) -> str:
    target = Path(path)
    digest = sha256_file(target)
    target.with_name(target.name + ".sha256").write_text(
        f"{digest}  {target.name}\n", encoding="utf-8"
    )
    return digest


def _verify_sidecar(path: Path) -> str:
    sidecar = path.with_name(path.name + ".sha256")
    if not sidecar.is_file():
        raise ValueError(f"missing SHA256 sidecar: {sidecar}")
    fields = sidecar.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2 or fields[1] != path.name:
        raise ValueError(f"malformed SHA256 sidecar: {sidecar}")
    observed = sha256_file(path)
    if fields[0] != observed:
        raise ValueError(f"SHA256 mismatch: {path}")
    return observed


def validate_frozen_authorities(
    *,
    design_path: str | Path,
    config_path: str | Path,
    freeze_receipt_path: str | Path,
) -> dict[str, str]:
    receipt_path = Path(freeze_receipt_path).resolve(strict=True)
    receipt = read_json(receipt_path)
    if receipt.get("scientific_contract_status") != "frozen":
        raise ValueError("FM0.1 scientific contract is not frozen")
    if receipt.get("implementation_status") != "authorized_not_started":
        raise ValueError("unexpected FM0.1 implementation authorization state")
    design = Path(design_path).resolve(strict=True)
    config = Path(config_path).resolve(strict=True)
    design_hash = sha256_file(design)
    config_hash = sha256_file(config)
    if design_hash != receipt.get("design_document", {}).get("sha256"):
        raise ValueError("FM0.1 design hash differs from freeze receipt")
    if config_hash != receipt.get("scientific_config", {}).get("sha256"):
        raise ValueError("FM0.1 config hash differs from freeze receipt")
    return {
        "design_sha256": design_hash,
        "config_sha256": config_hash,
        "freeze_receipt_sha256": sha256_file(receipt_path),
    }


def require_clean_git_revision(repo: str | Path, expected_sha: str) -> str:
    if len(expected_sha) != 40 or any(c not in "0123456789abcdef" for c in expected_sha):
        raise ValueError("expected Git SHA must be 40 lowercase hexadecimal characters")
    root = Path(repo).resolve(strict=True)
    observed = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "-C", str(root), "status", "--porcelain=v1", "--untracked-files=all"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if observed != expected_sha or dirty:
        raise ValueError("repository is not the exact clean expected revision")
    return observed


def _inspect_trusted_checkpoint(
    checkpoint_path: Path,
    *,
    contract: Mapping[str, Any],
    summary: Mapping[str, Any],
) -> dict[str, int]:
    """Structurally inspect a checkpoint created by this training entrypoint.

    ``torch.load(..., weights_only=False)`` can execute pickle payloads.  This
    helper is therefore intentionally private and is called only after the
    release's own three-artifact layout and SHA-256 bindings have passed.  Run
    it only on a trusted TWIRL run directory, never on an arbitrary download.
    """

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - base environment may omit Torch
        raise RuntimeError(
            "PyTorch is required to inspect an FM0 checkpoint; use "
            "inspect_checkpoint=False only in the non-Torch sidecar unit test"
        ) from exc

    from dataclasses import fields

    from .dataset import variant_view_indices
    from .model import FM0ModelConfig, TWIRLFM0, count_trainable_parameters
    from .training import (
        CHECKPOINT_SCHEMA_VERSION,
        FM0OptimizationConfig,
        make_optimizer_and_scheduler,
    )

    try:
        checkpoint = torch.load(
            checkpoint_path, map_location="cpu", weights_only=False
        )
    except Exception as exc:
        raise ValueError("FM0 checkpoint is not structurally loadable") from exc
    if not isinstance(checkpoint, Mapping):
        raise ValueError("FM0 checkpoint is not a mapping")
    required = {
        "schema_version",
        "model_state",
        "optimizer_state",
        "scheduler_state",
        "scaler_state",
        "rng_state",
        "progress",
        "sampler_state",
        "model_config",
        "optimization_config",
        "synthetic_dataset_config",
        "dataset_contract",
        "run_contract",
        "loss_history",
    }
    if not required.issubset(checkpoint):
        raise ValueError("FM0 checkpoint schema is incomplete")
    if checkpoint.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError("FM0 checkpoint schema mismatch")
    if checkpoint.get("run_contract") != dict(contract):
        raise ValueError("FM0 checkpoint run contract mismatch")

    progress = checkpoint.get("progress")
    if not isinstance(progress, Mapping):
        raise ValueError("FM0 checkpoint progress is malformed")
    checkpoint_step = progress.get("global_step")
    summary_step = summary.get("global_step")
    target_step = contract.get("target_step")
    if (
        isinstance(checkpoint_step, bool)
        or not isinstance(checkpoint_step, int)
        or checkpoint_step <= 0
        or isinstance(summary_step, bool)
        or not isinstance(summary_step, int)
        or isinstance(target_step, bool)
        or not isinstance(target_step, int)
        or checkpoint_step != summary_step
        or checkpoint_step != target_step
    ):
        raise ValueError("FM0 checkpoint global step differs from the completed run")

    model_payload = checkpoint.get("model_config")
    expected_model_fields = {field.name for field in fields(FM0ModelConfig)}
    if (
        not isinstance(model_payload, Mapping)
        or set(model_payload) != expected_model_fields
    ):
        raise ValueError("FM0 checkpoint model configuration schema mismatch")
    try:
        model_config = FM0ModelConfig(**dict(model_payload))
    except (TypeError, ValueError) as exc:
        raise ValueError("FM0 checkpoint model configuration is invalid") from exc
    variant = contract.get("variant")
    if not isinstance(variant, str):
        raise ValueError("FM0 run contract has no variant")
    if model_config.n_flux_views != len(variant_view_indices(variant)):
        raise ValueError("FM0 checkpoint model views differ from the run variant")
    if (
        model_config.architecture != contract.get("architecture")
        or model_config.architecture != summary.get("architecture")
        or variant != summary.get("variant")
    ):
        raise ValueError("FM0 checkpoint model identity differs from the release")

    model_state = checkpoint.get("model_state")
    if not isinstance(model_state, Mapping) or not model_state:
        raise ValueError("FM0 checkpoint model state is empty or malformed")
    if not all(isinstance(value, torch.Tensor) for value in model_state.values()):
        raise ValueError("FM0 checkpoint model state contains a non-tensor")
    try:
        model = TWIRLFM0(model_config)
        model.load_state_dict(model_state, strict=True)
    except (RuntimeError, TypeError, ValueError) as exc:
        raise ValueError(
            "FM0 checkpoint model state does not match its config"
        ) from exc
    parameter_count = count_trainable_parameters(model)
    if (
        isinstance(summary.get("parameter_count"), bool)
        or not isinstance(summary.get("parameter_count"), int)
        or summary.get("parameter_count") != parameter_count
    ):
        raise ValueError(
            "FM0 summary parameter count differs from the checkpoint model"
        )

    optimization = checkpoint.get("optimization_config")
    if not isinstance(optimization, Mapping) or optimization != contract.get(
        "optimization"
    ):
        raise ValueError("FM0 checkpoint optimization contract mismatch")
    effective_batch = optimization.get("effective_batch_windows")
    if (
        isinstance(effective_batch, bool)
        or not isinstance(effective_batch, int)
        or effective_batch <= 0
    ):
        raise ValueError("FM0 checkpoint effective batch is invalid")
    try:
        optimization_config = FM0OptimizationConfig(**dict(optimization))
        optimizer, scheduler = make_optimizer_and_scheduler(
            model, optimization_config
        )
        optimizer.load_state_dict(checkpoint["optimizer_state"])
        scheduler.load_state_dict(checkpoint["scheduler_state"])
    except (KeyError, TypeError, ValueError, RuntimeError) as exc:
        raise ValueError(
            "FM0 checkpoint optimizer or scheduler state is invalid"
        ) from exc
    sampler = checkpoint.get("sampler_state")
    if not isinstance(sampler, Mapping) or sampler.get("next_absolute_sample") != (
        checkpoint_step * effective_batch
    ):
        raise ValueError("FM0 checkpoint sampler state is inconsistent")

    dataset_contract = checkpoint.get("dataset_contract")
    if not isinstance(dataset_contract, Mapping):
        raise ValueError("FM0 checkpoint dataset contract is missing")
    if contract.get("synthetic_smoke"):
        dataset = checkpoint.get("synthetic_dataset_config")
        if (
            dataset_contract.get("kind") != "synthetic"
            or not isinstance(dataset, Mapping)
            or dataset_contract.get("config") != dataset
            or dataset.get("variant") != variant
            or dataset.get("seed") != contract.get("seed")
        ):
            raise ValueError("FM0 checkpoint synthetic dataset contract mismatch")
    else:
        release_binding = contract.get("input_release")
        if (
            dataset_contract.get("kind") != "fm0_input_release"
            or not isinstance(release_binding, Mapping)
            or dataset_contract.get("manifest_sha256")
            != release_binding.get("manifest_sha256")
            or dataset_contract.get("source_partition") != "poc_train"
            or dataset_contract.get("variant") != variant
            or dataset_contract.get("seed") != contract.get("seed")
        ):
            raise ValueError("FM0 checkpoint real dataset contract mismatch")
    for state_name in ("optimizer_state", "scheduler_state", "rng_state"):
        if not isinstance(checkpoint.get(state_name), Mapping):
            raise ValueError(f"FM0 checkpoint {state_name} is malformed")

    tensor_count = 0
    floating_tensor_count = 0

    def inspect_tensors(value: Any, context: str) -> None:
        nonlocal tensor_count, floating_tensor_count
        if isinstance(value, torch.Tensor):
            tensor_count += 1
            if value.is_floating_point() or value.is_complex():
                floating_tensor_count += 1
                if not bool(torch.isfinite(value).all()):
                    raise ValueError(
                        f"FM0 checkpoint {context} contains non-finite values"
                    )
            return
        if isinstance(value, Mapping):
            for key, item in value.items():
                inspect_tensors(item, f"{context}.{key}")
            return
        if isinstance(value, (list, tuple)):
            for index, item in enumerate(value):
                inspect_tensors(item, f"{context}[{index}]")

    for state_name in (
        "model_state",
        "optimizer_state",
        "scheduler_state",
        "scaler_state",
        "rng_state",
    ):
        inspect_tensors(checkpoint.get(state_name), state_name)
    if tensor_count == 0 or floating_tensor_count == 0:
        raise ValueError("FM0 checkpoint contains no learned tensor state")

    history = checkpoint.get("loss_history")
    if not isinstance(history, list) or len(history) != checkpoint_step:
        raise ValueError("FM0 checkpoint loss history is incomplete")
    for index, row in enumerate(history, start=1):
        if not isinstance(row, Mapping) or row.get("step") != float(index):
            raise ValueError("FM0 checkpoint loss history step is inconsistent")
        for name, value in row.items():
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ValueError(
                    f"FM0 checkpoint loss history {index}.{name} is non-finite"
                )
    if summary.get("final_metrics") != history[-1]:
        raise ValueError("FM0 summary metrics differ from checkpoint history")
    return {
        "tensor_count": tensor_count,
        "floating_tensor_count": floating_tensor_count,
    }


def validate_run_release(
    run_dir: str | Path,
    *,
    inspect_checkpoint: bool = True,
) -> dict[str, Any]:
    """Validate a trusted synthetic release and, by default, its checkpoint."""

    root = Path(run_dir).resolve(strict=True)
    contract_path = root / "run_contract.json"
    checkpoint_path = root / "checkpoint.pt"
    summary_path = root / "summary.json"
    hashes = {
        path.name: _verify_sidecar(path)
        for path in (contract_path, checkpoint_path, summary_path)
        if path.is_file()
    }
    if set(hashes) != {"run_contract.json", "checkpoint.pt", "summary.json"}:
        raise ValueError("FM0 synthetic run release is incomplete")
    contract = read_json(contract_path)
    summary = read_json(summary_path)
    if contract.get("schema_version") != RUN_CONTRACT_SCHEMA_VERSION:
        raise ValueError("FM0 run-contract schema mismatch")
    if summary.get("schema_version") != RUN_SUMMARY_SCHEMA_VERSION:
        raise ValueError("FM0 summary schema mismatch")
    if not contract.get("synthetic_smoke") or not summary.get("synthetic_only"):
        raise ValueError("this validator accepts only explicit synthetic smokes")
    if summary.get("run_contract_sha256") != hashes["run_contract.json"]:
        raise ValueError("summary does not bind the run contract")
    if summary.get("checkpoint_sha256") != hashes["checkpoint.pt"]:
        raise ValueError("summary does not bind the checkpoint")
    if not summary.get("passed") or int(summary.get("global_step", 0)) <= 0:
        raise ValueError("FM0 synthetic smoke did not complete")
    checkpoint_inspection = None
    if inspect_checkpoint:
        checkpoint_inspection = _inspect_trusted_checkpoint(
            checkpoint_path,
            contract=contract,
            summary=summary,
        )
    return {
        "schema_version": RUN_VALIDATION_SCHEMA_VERSION,
        "passed": True,
        "synthetic_only": True,
        "run_dir": str(root),
        "variant": summary.get("variant"),
        "architecture": summary.get("architecture"),
        "global_step": int(summary["global_step"]),
        "artifact_sha256": hashes,
        "checkpoint_inspected": bool(inspect_checkpoint),
        "checkpoint_tensor_count": (
            checkpoint_inspection["tensor_count"]
            if checkpoint_inspection is not None
            else None
        ),
    }


def validate_real_run_release(
    run_dir: str | Path,
    *,
    inspect_checkpoint: bool = True,
) -> dict[str, Any]:
    """Validate one trusted, checksum-bound real-data FM0.1 run release."""

    root = Path(run_dir).resolve(strict=True)
    contract_path = root / "run_contract.json"
    checkpoint_path = root / "checkpoint.pt"
    summary_path = root / "summary.json"
    hashes = {
        path.name: _verify_sidecar(path)
        for path in (contract_path, checkpoint_path, summary_path)
        if path.is_file()
    }
    if set(hashes) != {"run_contract.json", "checkpoint.pt", "summary.json"}:
        raise ValueError("FM0 real-data run release is incomplete")
    contract = read_json(contract_path)
    summary = read_json(summary_path)
    if contract.get("schema_version") != REAL_RUN_CONTRACT_SCHEMA_VERSION:
        raise ValueError("FM0 real run-contract schema mismatch")
    if summary.get("schema_version") != REAL_RUN_SUMMARY_SCHEMA_VERSION:
        raise ValueError("FM0 real summary schema mismatch")
    if contract.get("synthetic_smoke") or summary.get("synthetic_only"):
        raise ValueError("real-data validator rejects synthetic runs")
    if contract.get("real_data_consumed") is not True:
        raise ValueError("real-data run contract does not declare real input")
    release = contract.get("input_release")
    if not isinstance(release, Mapping):
        raise ValueError("real-data run lacks an input-release binding")
    receipt_path = Path(str(release.get("receipt_path", ""))).resolve(strict=True)
    manifest_path = Path(str(release.get("manifest_path", ""))).resolve(strict=True)
    if sha256_file(receipt_path) != release.get("receipt_sha256"):
        raise ValueError("real-data input receipt changed after training")
    if sha256_file(manifest_path) != release.get("manifest_sha256"):
        raise ValueError("real-data manifest changed after training")
    receipt = read_json(receipt_path)
    if (
        receipt.get("passed") is not True
        or receipt.get("scientific_training_eligible") is not True
        or receipt.get("release", {}).get("manifest_sha256")
        != release.get("manifest_sha256")
    ):
        raise ValueError("real-data input receipt does not authorize training")
    if summary.get("run_contract_sha256") != hashes["run_contract.json"]:
        raise ValueError("summary does not bind the real run contract")
    if summary.get("checkpoint_sha256") != hashes["checkpoint.pt"]:
        raise ValueError("summary does not bind the real checkpoint")
    if not summary.get("passed") or int(summary.get("global_step", 0)) <= 0:
        raise ValueError("FM0 real-data training did not complete")
    checkpoint_inspection = None
    if inspect_checkpoint:
        checkpoint_inspection = _inspect_trusted_checkpoint(
            checkpoint_path,
            contract=contract,
            summary=summary,
        )
    return {
        "schema_version": REAL_RUN_VALIDATION_SCHEMA_VERSION,
        "passed": True,
        "synthetic_only": False,
        "real_data_consumed": True,
        "run_dir": str(root),
        "variant": summary.get("variant"),
        "architecture": summary.get("architecture"),
        "global_step": int(summary["global_step"]),
        "artifact_sha256": hashes,
        "checkpoint_inspected": bool(inspect_checkpoint),
        "checkpoint_tensor_count": (
            checkpoint_inspection["tensor_count"]
            if checkpoint_inspection is not None
            else None
        ),
    }
