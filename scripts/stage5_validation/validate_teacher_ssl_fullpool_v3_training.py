#!/usr/bin/env python3
"""Validate and freeze the completed native-v3 five-fold SSL training run."""
from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import csv
from dataclasses import asdict
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.harmonic_cnn import (  # noqa: E402
    HarmonicModelConfig,
    build_harmonic_cnn,
    profile_branches,
)
from twirl.vetting.harmonic_ssl import (  # noqa: E402
    EventPreservingAugmentationConfig,
    VICRegConfig,
)
from twirl.vetting.ssl_full_pool_eligibility import (  # noqa: E402
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_FULL_OBSERVATIONS,
)
from twirl.vetting.ssl_full_pool_numeric import (  # noqa: E402
    MODEL_INPUT_NUMERIC_AUTHORITY_NAMES,
    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
    MODEL_INPUT_NUMERIC_RELEASE_BINDING,
    validate_numeric_gate_release,
    validate_numeric_native_freshness,
)
from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    FULLPOOL_SSL_ADAMW_BETAS,
    FULLPOOL_SSL_ADAMW_EPS,
    FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
    FULLPOOL_SSL_CHECKPOINT_SCHEMA,
    FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
    FULLPOOL_SSL_ENCODER_NAME,
    FULLPOOL_SSL_GRADIENT_COVERAGE_SCHEMA,
    FULLPOOL_SSL_LEARNING_RATE,
    FULLPOOL_SSL_MODEL_NAMESPACE,
    FULLPOOL_SSL_N_FOLDS,
    FULLPOOL_SSL_OPTIMIZER_CONTRACT_SCHEMA,
    FULLPOOL_SSL_OPTIMIZER_PARAMETER_SCOPE,
    FULLPOOL_SSL_OPTIMIZER_TYPE,
    FULLPOOL_SSL_PROFILE,
    FULLPOOL_SSL_RUN_CONTRACT_SCHEMA,
    FULLPOOL_SSL_RUN_ID,
    FULLPOOL_SSL_SELECTION_SCHEMA,
    FULLPOOL_SSL_SUMMARY_SCHEMA,
    FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA,
    FULLPOOL_SSL_WEIGHT_DECAY,
    _projection_head,
    _verify_native_files,
    load_fullpool_ssl_registry,
    select_fullpool_ssl_fold,
)


COMPLETION_RELEASE_SCHEMA = (
    "twirl_teacher_ssl_fullpool_v3_five_fold_completion_release_v1"
)
COMPLETION_RELEASE_BINDING = (
    "teacher_ssl_fullpool_v3_effective_quality_adp_five_fold_complete_v1"
)
PRODUCTION_EPOCHS = 20
PRODUCTION_BATCH_SIZE = 64
PRODUCTION_WORKERS = 8
PRODUCTION_CHECKPOINT_EVERY = 1
_SHA40 = re.compile(r"^[0-9a-f]{40}$")
_TIC_DISJOINT_KEYS = {
    "held_fold_tics",
    "fixed_test_tics",
    "reserved_prospective_tics",
}
_ADAMW_GROUP_KEYS = {
    "lr",
    "betas",
    "eps",
    "weight_decay",
    "amsgrad",
    "foreach",
    "maximize",
    "capturable",
    "differentiable",
    "fused",
    "decoupled_weight_decay",
    "params",
}
_LOCKED_OPTIMIZER_CONFIG = {
    "schema_version": FULLPOOL_SSL_OPTIMIZER_CONTRACT_SCHEMA,
    "optimizer_type": FULLPOOL_SSL_OPTIMIZER_TYPE,
    "parameter_group_count": 1,
    "defaults": {
        "lr": FULLPOOL_SSL_LEARNING_RATE,
        "betas": list(FULLPOOL_SSL_ADAMW_BETAS),
        "eps": FULLPOOL_SSL_ADAMW_EPS,
        "weight_decay": FULLPOOL_SSL_WEIGHT_DECAY,
        "amsgrad": False,
        "foreach": None,
        "maximize": False,
        "capturable": False,
        "differentiable": False,
        "fused": None,
        "decoupled_weight_decay": True,
    },
}


def _no_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"JSON contains duplicate key {key!r}")
        value[key] = item
    return value


def _reject_constant(value: str) -> None:
    raise ValueError(f"JSON contains non-finite constant {value!r}")


def _read_json(path: Path, *, context: str) -> dict[str, Any]:
    try:
        value = json.loads(
            Path(path).read_text(encoding="utf-8"),
            object_pairs_hook=_no_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{context} is not strict UTF-8 JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a JSON object")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _json_normalized(value: Any) -> Any:
    return json.loads(
        json.dumps(value, sort_keys=True, allow_nan=False)
    )


def _metadata(path: Path) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve(strict=True)
    before = resolved.stat()
    if not resolved.is_file() or before.st_size <= 0:
        raise ValueError(f"artifact is missing or empty: {resolved}")
    digest = _sha256(resolved)
    after = resolved.stat()
    if (
        int(before.st_size),
        int(before.st_mtime_ns),
        int(before.st_ctime_ns),
        int(before.st_dev),
        int(before.st_ino),
    ) != (
        int(after.st_size),
        int(after.st_mtime_ns),
        int(after.st_ctime_ns),
        int(after.st_dev),
        int(after.st_ino),
    ):
        raise RuntimeError(f"artifact changed while hashing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": int(after.st_size),
        "sha256": digest,
    }


def _require_bound_artifact(value: Any, *, context: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != {
        "path",
        "size_bytes",
        "sha256",
    }:
        raise ValueError(f"{context} is not an exact artifact binding")
    observed = _metadata(Path(str(value.get("path", ""))))
    if observed != dict(value):
        raise ValueError(f"{context} changed after training preflight")
    return observed


def _numeric_gate_projection(
    payload: Mapping[str, Any],
    *,
    binding: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "binding": dict(binding),
        "release_binding": payload.get("release_binding"),
        "schema_version": payload.get("schema_version"),
        "envelope_canonical_sha256": payload.get(
            "envelope_canonical_sha256"
        ),
        "counts": dict(payload.get("counts", {})),
        "identity_hashes": dict(payload.get("identity_hashes", {})),
        "code_revision": payload.get("code_revision"),
        "passed": payload.get("passed"),
    }


def _exact_int(value: Any, expected: int, *, context: str) -> None:
    if type(value) is not int or int(value) != int(expected):
        raise ValueError(f"{context} must equal {expected}")


def _exact_float(value: Any, expected: float, *, context: str) -> None:
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not math.isfinite(float(value))
        or float(value) != float(expected)
    ):
        raise ValueError(f"{context} must equal {expected}")


def _expected_epoch_coverage(
    n_selected_observations: int,
) -> dict[str, int]:
    selected = int(n_selected_observations)
    remainder = selected % PRODUCTION_BATCH_SIZE
    singleton_skipped = int(remainder == 1)
    return {
        "selected_observations": selected,
        "observations_per_epoch": selected - singleton_skipped,
        "optimizer_steps_per_epoch": (
            selected // PRODUCTION_BATCH_SIZE + int(remainder >= 2)
        ),
        "singleton_observations_skipped_per_epoch": singleton_skipped,
    }


def _read_history(
    path: Path,
    *,
    expected_n_observations: int,
) -> tuple[list[str], list[dict[str, float]]]:
    try:
        reader = csv.DictReader(
            io.StringIO(Path(path).read_text(encoding="utf-8"), newline="")
        )
        columns = reader.fieldnames
        raw_rows = list(reader)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise ValueError("fold history is not valid UTF-8 CSV") from exc
    if (
        not isinstance(columns, list)
        or len(columns) != len(set(columns))
        or not {"epoch", "n_observations", "loss"} <= set(columns)
    ):
        raise ValueError("fold history schema is invalid")
    if len(raw_rows) != PRODUCTION_EPOCHS:
        raise ValueError("fold history must contain exactly 20 epochs")
    rows: list[dict[str, float]] = []
    for index, raw in enumerate(raw_rows, start=1):
        row: dict[str, float] = {}
        for name in columns:
            try:
                value = float(raw[name])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"fold history {name} is not numeric") from exc
            if not math.isfinite(value):
                raise ValueError(f"fold history {name} is non-finite")
            row[name] = value
        if (
            row["epoch"] != index
            or row["n_observations"] != expected_n_observations
        ):
            raise ValueError("fold history epoch/count sequence is invalid")
        rows.append(row)
    return columns, rows


def _finite_tensor_count(value: Any, *, context: str) -> int:
    import torch

    if torch.is_tensor(value):
        if not bool(torch.isfinite(value).all().item()):
            raise ValueError(f"{context} contains a non-finite tensor")
        return 1
    if isinstance(value, Mapping):
        return sum(
            _finite_tensor_count(item, context=f"{context}.{name}")
            for name, item in value.items()
        )
    if isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    ):
        return sum(
            _finite_tensor_count(item, context=f"{context}[{index}]")
            for index, item in enumerate(value)
        )
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{context} contains a non-finite scalar")
    return 0


def _expected_parameter_manifest(
    encoder: Any,
    projector: Any,
) -> list[dict[str, Any]]:
    if getattr(encoder, "profile", None) != FULLPOOL_SSL_PROFILE:
        raise ValueError("full-pool SSL encoder profile differs")
    branches = profile_branches(FULLPOOL_SSL_PROFILE)
    prefixes: list[str] = []
    if "chronology" in branches:
        prefixes.extend(
            ("small_encoder.", "supp_encoder.", "supplement_gate.")
        )
    if {"harmonic", "single_fold"} & branches:
        prefixes.append("harmonic_encoder.")
    if "local" in branches:
        prefixes.append("local_encoder.")
    if "periodogram" in branches:
        prefixes.append("periodogram_encoder.")
    if "metadata" in branches and getattr(
        encoder, "metadata_encoder", None
    ) is not None:
        prefixes.append("metadata_encoder.")
    prefixes.extend(("branch_gate.", "fusion."))
    manifest: list[dict[str, Any]] = []
    for module_prefix, module in (
        ("encoder", encoder),
        ("projector", projector),
    ):
        for name, parameter in module.named_parameters():
            if module_prefix == "encoder" and not name.startswith(
                tuple(prefixes)
            ):
                continue
            manifest.append(
                {
                    "index": len(manifest),
                    "name": f"{module_prefix}.{name}",
                    "shape": [int(value) for value in parameter.shape],
                    "dtype": str(parameter.dtype),
                    "numel": int(parameter.numel()),
                    "requires_grad": bool(parameter.requires_grad),
                }
            )
    return manifest


def _expected_gradient_coverage(
    encoder: Any,
    projector: Any,
    manifest: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    active_names = [str(item["name"]) for item in manifest]
    active_name_set = set(active_names)
    excluded_names = [
        f"{prefix}.{name}"
        for prefix, module in (("encoder", encoder), ("projector", projector))
        for name, _ in module.named_parameters()
        if f"{prefix}.{name}" not in active_name_set
    ]
    return {
        "schema_version": FULLPOOL_SSL_GRADIENT_COVERAGE_SCHEMA,
        "verified_after_real_ssl_backward": True,
        "loss_path": "encoder.embedding->projector->vicreg",
        "active_parameter_count": len(active_names),
        "active_parameter_names_sha256": _canonical_sha256(active_names),
        "excluded_parameter_count": len(excluded_names),
        "excluded_parameter_names_sha256": _canonical_sha256(
            excluded_names
        ),
    }


def _validate_optimizer_checkpoint(
    checkpoint: Mapping[str, Any],
    *,
    fold: int,
    model_config: Mapping[str, Any],
    projection_architecture: Sequence[int],
    expected_step: int,
) -> dict[str, int]:
    import torch

    locked_checkpoint_model_config = asdict(
        HarmonicModelConfig(metadata_dim=0)
    )
    locked_contract_model_config = _json_normalized(
        locked_checkpoint_model_config
    )
    locked_projection = [
        2 * int(locked_checkpoint_model_config["embedding_dim"]),
        128,
        64,
    ]
    if (
        dict(model_config) != locked_contract_model_config
        or list(projection_architecture) != locked_projection
        or checkpoint.get("model_config") != locked_checkpoint_model_config
        or checkpoint.get("projection_architecture") != locked_projection
    ):
        raise ValueError(f"fold {fold} model architecture differs")
    encoder = build_harmonic_cnn(
        HarmonicModelConfig(**locked_checkpoint_model_config),
        profile=FULLPOOL_SSL_PROFILE,
    )
    projector = _projection_head(locked_projection[0])
    encoder_state = checkpoint.get("encoder_state_dict")
    projector_state = checkpoint.get("projection_state_dict")
    if not isinstance(encoder_state, Mapping) or not isinstance(
        projector_state, Mapping
    ):
        raise ValueError(f"fold {fold} model state is invalid")
    try:
        encoder.load_state_dict(encoder_state, strict=True)
        projector.load_state_dict(projector_state, strict=True)
    except Exception as exc:
        raise ValueError(f"fold {fold} model state schema differs") from exc
    manifest = _expected_parameter_manifest(encoder, projector)
    expected_contract = {
        **_LOCKED_OPTIMIZER_CONFIG,
        "parameter_scope": FULLPOOL_SSL_OPTIMIZER_PARAMETER_SCOPE,
        "parameter_manifest": manifest,
        "parameter_manifest_sha256": _canonical_sha256(manifest),
        "parameter_count": len(manifest),
        "parameter_numel": sum(item["numel"] for item in manifest),
        "parameter_ids": list(range(len(manifest))),
        "state_parameter_count": len(manifest),
        "state_fields": ["step", "exp_avg", "exp_avg_sq"],
        "state_step": int(expected_step),
        "state_step_dtype": "torch.float32",
        "gradient_coverage": _expected_gradient_coverage(
            encoder,
            projector,
            manifest,
        ),
    }
    if checkpoint.get("optimizer_config") != _LOCKED_OPTIMIZER_CONFIG:
        raise ValueError(f"fold {fold} optimizer config differs")
    if checkpoint.get("optimizer_contract") != expected_contract:
        raise ValueError(
            f"fold {fold} optimizer parameter contract differs"
        )

    optimizer = checkpoint.get("optimizer_state_dict")
    if (
        not isinstance(optimizer, Mapping)
        or not isinstance(optimizer.get("state"), Mapping)
        or not isinstance(optimizer.get("param_groups"), list)
        or len(optimizer["param_groups"]) != 1
        or not isinstance(optimizer["param_groups"][0], Mapping)
    ):
        raise ValueError(f"fold {fold} optimizer state is empty or invalid")
    group = optimizer["param_groups"][0]
    if set(group) != _ADAMW_GROUP_KEYS:
        raise ValueError(f"fold {fold} AdamW group field inventory differs")
    defaults = _LOCKED_OPTIMIZER_CONFIG["defaults"]
    for name, expected in defaults.items():
        observed = group.get(name)
        if name == "betas" and observed is not None:
            observed = list(observed)
        if observed != expected:
            raise ValueError(
                f"fold {fold} optimizer {name} differs from locked AdamW"
            )
    parameter_ids = group.get("params")
    if (
        not isinstance(parameter_ids, list)
        or parameter_ids != list(range(len(manifest)))
    ):
        raise ValueError(
            f"fold {fold} optimizer parameter ID coverage differs"
        )
    state = optimizer["state"]
    if set(state) != set(parameter_ids):
        raise ValueError(f"fold {fold} optimizer state coverage differs")
    for parameter_id, item in zip(parameter_ids, manifest, strict=True):
        parameter_state = state[parameter_id]
        if (
            not isinstance(parameter_state, Mapping)
            or set(parameter_state) != {"step", "exp_avg", "exp_avg_sq"}
            or not torch.is_tensor(parameter_state["step"])
            or parameter_state["step"].numel() != 1
            or parameter_state["step"].dtype != torch.float32
            or float(parameter_state["step"].item())
            != float(expected_step)
            or not torch.is_tensor(parameter_state["exp_avg"])
            or list(parameter_state["exp_avg"].shape) != item["shape"]
            or str(parameter_state["exp_avg"].dtype) != item["dtype"]
            or not bool(
                torch.isfinite(parameter_state["exp_avg"]).all().item()
            )
            or not torch.is_tensor(parameter_state["exp_avg_sq"])
            or list(parameter_state["exp_avg_sq"].shape) != item["shape"]
            or str(parameter_state["exp_avg_sq"].dtype) != item["dtype"]
            or not bool(
                torch.isfinite(parameter_state["exp_avg_sq"]).all().item()
            )
        ):
            raise ValueError(
                f"fold {fold} optimizer state tensors differ for "
                f"{item['name']}"
            )
    optimizer_count = _finite_tensor_count(
        optimizer, context=f"fold {fold} optimizer_state_dict"
    )
    if optimizer_count < len(manifest) * 2:
        raise ValueError(f"fold {fold} optimizer state tensor coverage differs")
    return {
        "encoder_state_dict": _finite_tensor_count(
            encoder_state,
            context=f"fold {fold} encoder_state_dict",
        ),
        "projection_state_dict": _finite_tensor_count(
            projector_state,
            context=f"fold {fold} projection_state_dict",
        ),
        "optimizer_state_dict": optimizer_count,
    }


def _validate_native_record(
    value: Any,
    *,
    expected_code_revision: str,
    index: int,
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"native record {index} is not a mapping")
    expected = {
        "native_contract_version": MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
            "native_contract_version"
        ],
        "native_release_binding": MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
            "release_binding"
        ],
        "native_expected_git_sha": expected_code_revision,
        "native_builder_contract_version": (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "builder_contract_version"
            ]
        ),
        "native_builder_code_sha256": MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
            "builder_code_sha256"
        ],
        "native_detrend_contract_version": (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "detrend_contract_version"
            ]
        ),
        "native_detrend_config_sha256": (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "detrend_config_sha256"
            ]
        ),
        "native_detrend_quality_source": "final_effective_quality",
        "raw_photometry_only": True,
        "compact_adp_photometry_reused": False,
        "compact_adp_flux_reused": False,
        "periodogram_n": MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
            "periodogram_n"
        ],
        "hash_verified_now": True,
        "root_contract_verified_now": True,
        "group_identities_verified_now": True,
    }
    if {name: value.get(name) for name in expected} != expected:
        raise ValueError(f"native record {index} has stale/reused provenance")
    path = Path(str(value.get("native_h5_path", ""))).resolve(strict=True)
    stat = path.stat()
    observed_stat = {
        "native_h5_size_bytes": int(stat.st_size),
        "native_h5_mtime_ns": int(stat.st_mtime_ns),
        "native_h5_ctime_ns": int(stat.st_ctime_ns),
        "native_h5_device": int(stat.st_dev),
        "native_h5_inode": int(stat.st_ino),
    }
    expected_stat = {
        name: int(value.get(name, -1)) for name in observed_stat
    }
    if (
        observed_stat != expected_stat
        or _sha256(path) != value.get("native_h5_sha256")
    ):
        raise ValueError(f"native record {index} HDF5 binding changed")
    try:
        import h5py

        with h5py.File(path, "r") as h5:
            root_expected = {
                "contract_version": expected["native_contract_version"],
                "release_binding": expected["native_release_binding"],
                "expected_git_sha": expected_code_revision,
                "builder_contract_version": expected[
                    "native_builder_contract_version"
                ],
                "builder_code_sha256": expected[
                    "native_builder_code_sha256"
                ],
                "detrend_contract_version": expected[
                    "native_detrend_contract_version"
                ],
                "detrend_config_sha256": expected[
                    "native_detrend_config_sha256"
                ],
                "detrend_quality_source": "final_effective_quality",
            }
            if {
                name: str(h5.attrs.get(name, "")) for name in root_expected
            } != root_expected:
                raise ValueError(
                    f"native record {index} HDF5 root provenance changed"
                )
            for name, expected_value in (
                ("raw_photometry_only", 1),
                ("compact_adp_photometry_reused", 0),
                ("compact_adp_flux_reused", 0),
                (
                    "periodogram_n",
                    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS["periodogram_n"],
                ),
            ):
                if int(h5.attrs.get(name, -1)) != expected_value:
                    raise ValueError(
                        f"native record {index} HDF5 root {name} changed"
                    )
    except OSError as exc:
        raise ValueError(f"native record {index} is not loadable HDF5") from exc
    return _metadata(path)


def _validate_checkpoint(
    path: Path,
    *,
    fold: int,
    contract_sha256: str,
    selection: Mapping[str, Any],
    expected_freshness: Mapping[str, Any],
    expected_numeric_release: Mapping[str, Any],
    model_config: Mapping[str, Any],
    projection_architecture: Sequence[int],
    optimizer_config: Mapping[str, Any],
    epoch_coverage: Mapping[str, Any],
    history_columns: list[str],
    history_rows: list[dict[str, float]],
) -> dict[str, int]:
    try:
        import torch

        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:
        raise ValueError(f"fold {fold} checkpoint is not loadable") from exc
    if not isinstance(checkpoint, Mapping):
        raise ValueError(f"fold {fold} checkpoint is not a mapping")
    if (
        checkpoint.get("schema_version") != FULLPOOL_SSL_CHECKPOINT_SCHEMA
        or checkpoint.get("run_id") != FULLPOOL_SSL_RUN_ID
        or checkpoint.get("checkpoint_namespace")
        != FULLPOOL_SSL_CHECKPOINT_NAMESPACE
        or checkpoint.get("numeric_release_binding")
        != MODEL_INPUT_NUMERIC_RELEASE_BINDING
        or checkpoint.get("native_freshness") != expected_freshness
        or checkpoint.get("numeric_gate_release") != expected_numeric_release
        or checkpoint.get("run_contract_sha256") != contract_sha256
        or checkpoint.get("selection_audit") != selection
        or checkpoint.get("optimizer_config") != optimizer_config
        or checkpoint.get("epoch_coverage") != epoch_coverage
    ):
        raise ValueError(f"fold {fold} checkpoint authority/schema differs")
    _exact_int(checkpoint.get("fold"), fold, context="checkpoint fold")
    _exact_int(
        checkpoint.get("epochs"),
        PRODUCTION_EPOCHS,
        context="checkpoint epochs",
    )
    checkpoint_history = checkpoint.get("history")
    if (
        not isinstance(checkpoint_history, list)
        or len(checkpoint_history) != PRODUCTION_EPOCHS
    ):
        raise ValueError(f"fold {fold} checkpoint history is incomplete")
    for expected, observed in zip(
        history_rows, checkpoint_history, strict=True
    ):
        if (
            not isinstance(observed, Mapping)
            or set(observed) != set(history_columns)
            or any(float(observed[name]) != expected[name] for name in history_columns)
        ):
            raise ValueError(f"fold {fold} checkpoint history differs from CSV")
    return _validate_optimizer_checkpoint(
        checkpoint,
        fold=fold,
        model_config=model_config,
        projection_architecture=projection_architecture,
        expected_step=(
            PRODUCTION_EPOCHS
            * int(epoch_coverage["optimizer_steps_per_epoch"])
        ),
    )


def validate_teacher_ssl_fullpool_training(
    *,
    model_root: Path,
    expected_code_revision: str,
) -> dict[str, Any]:
    """Validate exact fold 0--4 completion without trusting job exit codes."""

    if _SHA40.fullmatch(str(expected_code_revision)) is None:
        raise ValueError("expected code revision must be a lowercase 40-hex SHA")
    root = Path(model_root).expanduser().resolve(strict=True)
    expected_suffix = (
        "model_runs",
        FULLPOOL_SSL_MODEL_NAMESPACE,
        "training",
        "five_fold",
    )
    if tuple(root.parts[-len(expected_suffix) :]) != expected_suffix:
        raise ValueError("training root is outside the fresh model namespace")
    fold_root = root / "encoder_pretraining"
    observed_names = {
        path.name for path in fold_root.iterdir() if path.name.startswith("fold_")
    }
    expected_names = {f"fold_{fold}" for fold in range(FULLPOOL_SSL_N_FOLDS)}
    if observed_names != expected_names:
        raise ValueError("five-fold directory inventory is not exactly fold 0--4")

    expected_freshness = validate_numeric_native_freshness(
        {
            **MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
            "expected_git_sha": expected_code_revision,
        },
        context="completion expected native freshness",
        expected_code_revision=expected_code_revision,
    )
    common: dict[str, Any] | None = None
    fold_records: list[dict[str, Any]] = []
    seeds: set[int] = set()
    validated_numeric_releases: dict[str, dict[str, Any]] = {}
    loaded_registries: dict[
        tuple[str, str], tuple[Any, dict[str, Any]]
    ] = {}
    for fold in range(FULLPOOL_SSL_N_FOLDS):
        fold_dir = fold_root / f"fold_{fold}"
        if not fold_dir.is_dir() or fold_dir.is_symlink():
            raise ValueError(f"fold {fold} directory is invalid")
        summary_path = fold_dir / "summary.json"
        summary_binding = _metadata(summary_path)
        summary = _read_json(summary_path, context=f"fold {fold} summary")
        if (
            summary.get("schema_version") != FULLPOOL_SSL_SUMMARY_SCHEMA
            or summary.get("run_id") != FULLPOOL_SSL_RUN_ID
            or summary.get("encoder_name") != FULLPOOL_SSL_ENCODER_NAME
            or summary.get("checkpoint_namespace")
            != FULLPOOL_SSL_CHECKPOINT_NAMESPACE
            or summary.get("numeric_release_binding")
            != MODEL_INPUT_NUMERIC_RELEASE_BINDING
            or summary.get("native_freshness") != expected_freshness
        ):
            raise ValueError(f"fold {fold} summary authority/schema differs")
        _exact_int(summary.get("fold"), fold, context="summary fold")
        _exact_int(
            summary.get("completed_epochs"),
            PRODUCTION_EPOCHS,
            context="summary completed epochs",
        )
        if type(summary.get("global_step")) is not int or int(
            summary["global_step"]
        ) < 1:
            raise ValueError(f"fold {fold} summary global step is invalid")
        if (
            summary.get("labels_loaded") is not False
            or summary.get("fixed_test_status")
            != "host_excluded_tensors_not_constructed"
            or summary.get("prospective_status")
            != "host_excluded_tensors_not_constructed"
            or summary.get("automatic_production_promotion") is not False
        ):
            raise ValueError(f"fold {fold} summary isolation gates failed")

        contract_path = Path(str(summary.get("run_contract", ""))).resolve(
            strict=True
        )
        checkpoint_path = Path(str(summary.get("checkpoint", ""))).resolve(
            strict=True
        )
        history_path = Path(str(summary.get("history", ""))).resolve(
            strict=True
        )
        if (
            contract_path != fold_dir / "run_contract.json"
            or checkpoint_path != fold_dir / "encoder.pt"
            or history_path != fold_dir / "history.csv"
        ):
            raise ValueError(f"fold {fold} artifacts escape their namespace")
        contract_binding = _metadata(contract_path)
        checkpoint_binding = _metadata(checkpoint_path)
        history_binding = _metadata(history_path)
        if (
            summary.get("run_contract_sha256")
            != contract_binding["sha256"]
            or summary.get("checkpoint_sha256")
            != checkpoint_binding["sha256"]
            or summary.get("history_sha256") != history_binding["sha256"]
        ):
            raise ValueError(f"fold {fold} summary artifact binding differs")
        contract = _read_json(
            contract_path, context=f"fold {fold} run contract"
        )
        if (
            contract.get("schema_version") != FULLPOOL_SSL_RUN_CONTRACT_SCHEMA
            or contract.get("run_id") != FULLPOOL_SSL_RUN_ID
            or contract.get("encoder_name") != FULLPOOL_SSL_ENCODER_NAME
            or contract.get("checkpoint_namespace")
            != FULLPOOL_SSL_CHECKPOINT_NAMESPACE
            or contract.get("code_revision") != expected_code_revision
            or contract.get("numeric_release_binding")
            != MODEL_INPUT_NUMERIC_RELEASE_BINDING
            or contract.get("native_freshness") != expected_freshness
        ):
            raise ValueError(f"fold {fold} run contract authority/schema differs")
        for name, expected in (
            ("fold", fold),
            ("epochs", PRODUCTION_EPOCHS),
            ("batch_size", PRODUCTION_BATCH_SIZE),
            ("workers", PRODUCTION_WORKERS),
            (
                "seed",
                FULLPOOL_SSL_DEFAULT_TRAINING_SEED + 1000 * fold,
            ),
            ("checkpoint_every", PRODUCTION_CHECKPOINT_EVERY),
        ):
            _exact_int(
                contract.get(name),
                expected,
                context=f"fold {fold} {name}",
            )
        _exact_float(
            contract.get("learning_rate"),
            FULLPOOL_SSL_LEARNING_RATE,
            context=f"fold {fold} learning rate",
        )
        _exact_float(
            contract.get("weight_decay"),
            FULLPOOL_SSL_WEIGHT_DECAY,
            context=f"fold {fold} weight decay",
        )
        locked_model_config = _json_normalized(
            asdict(HarmonicModelConfig(metadata_dim=0))
        )
        locked_projection = [
            2 * int(locked_model_config["embedding_dim"]),
            128,
            64,
        ]
        if (
            contract.get("model_config") != locked_model_config
            or contract.get("augmentation_config")
            != _json_normalized(
                asdict(EventPreservingAugmentationConfig())
            )
            or contract.get("vicreg_config")
            != _json_normalized(asdict(VICRegConfig()))
            or contract.get("projection_architecture")
            != locked_projection
            or contract.get("optimizer_config")
            != _LOCKED_OPTIMIZER_CONFIG
        ):
            raise ValueError(f"fold {fold} model/optimizer config differs")
        if (
            contract.get("require_cuda") is not True
            or contract.get("max_rows") is not None
            or contract.get("required_observation_ids") != []
            or contract.get("labels_loaded") is not False
            or contract.get("fixed_test_tensors_constructed") is not False
            or contract.get("prospective_sector_tensors_constructed") is not False
            or contract.get("embedding_export") is not False
            or contract.get("neighbor_probe") is not False
            or contract.get("native_hashes_verified_now") is not True
            or contract.get("native_root_contracts_verified_now") is not True
            or contract.get("native_group_identities_verified_now") is not True
        ):
            raise ValueError(f"fold {fold} smoke/reuse/isolation gate failed")
        seed = int(contract["seed"])
        if seed in seeds:
            raise ValueError("five-fold training seeds are not distinct")
        seeds.add(seed)

        selection = contract.get("selection_audit")
        tic_disjoint = (
            selection.get("tic_disjoint")
            if isinstance(selection, Mapping)
            else None
        )
        if (
            not isinstance(selection, Mapping)
            or selection.get("selection_schema_version")
            != FULLPOOL_SSL_SELECTION_SCHEMA
            or selection.get("held_out_fold") != fold
            or selection.get("n_registry_observations")
            != PRODUCTION_FULL_OBSERVATIONS
            or selection.get("n_eligible_observations")
            != PRODUCTION_ELIGIBLE_OBSERVATIONS
            or selection.get("max_rows") is not None
            or selection.get("required_observation_ids") != []
            or selection.get("n_required_observations") != 0
            or selection.get("required_observations_selected") is not True
            or not isinstance(tic_disjoint, Mapping)
            or set(tic_disjoint) != _TIC_DISJOINT_KEYS
            or any(tic_disjoint.get(name) is not True for name in tic_disjoint)
            or type(selection.get("n_selected_observations")) is not int
            or type(selection.get("n_held_observations")) is not int
            or int(selection["n_selected_observations"])
            + int(selection["n_held_observations"])
            != PRODUCTION_ELIGIBLE_OBSERVATIONS
            or re.fullmatch(
                r"[0-9a-f]{64}",
                str(selection.get("selected_rows_sha256", "")),
            )
            is None
            or re.fullmatch(
                r"[0-9a-f]{64}",
                str(selection.get("selected_tics_sha256", "")),
            )
            is None
        ):
            raise ValueError(f"fold {fold} selection is not full production")
        authority = contract.get("training_authority")
        if (
            not isinstance(authority, Mapping)
            or authority.get("schema_version")
            != FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA
            or authority.get("production_lock_passed") is not True
            or authority.get("source_provenance_verified") is not True
            or authority.get("native_freshness") != expected_freshness
        ):
            raise ValueError(f"fold {fold} training authority differs")
        bindings = authority.get("authority_bindings")
        expected_training_bindings = {
            "eligibility_exclusions",
            "eligibility_summary",
            "native_registry",
            "native_registry_summary",
            "native_release_summary",
            "registry",
            "registry_summary",
            "numeric_gate_release",
        }
        if (
            not isinstance(bindings, Mapping)
            or set(bindings) != expected_training_bindings
        ):
            raise ValueError(f"fold {fold} lacks authority bindings")
        normalized_bindings = {
            str(name): _require_bound_artifact(
                binding, context=f"fold {fold} authority {name}"
            )
            for name, binding in bindings.items()
        }
        numeric_release_binding = normalized_bindings.get(
            "numeric_gate_release"
        )
        if numeric_release_binding is None:
            raise ValueError(
                f"fold {fold} lacks the numeric-gate release binding"
            )
        numeric_release = contract.get("numeric_gate_release")
        if (
            not isinstance(numeric_release, Mapping)
            or numeric_release.get("release_binding")
            != MODEL_INPUT_NUMERIC_RELEASE_BINDING
            or numeric_release.get("passed") is not True
            or authority.get("numeric_gate_release") != numeric_release
            or summary.get("numeric_gate_release") != numeric_release
        ):
            raise ValueError(f"fold {fold} numeric release binding differs")
        expected_numeric_authorities = {
            "ssl_registry": normalized_bindings["registry"],
            "ssl_registry_summary": normalized_bindings["registry_summary"],
            "native_registry": normalized_bindings["native_registry"],
            "native_registry_summary": normalized_bindings[
                "native_registry_summary"
            ],
            "native_release_summary": normalized_bindings[
                "native_release_summary"
            ],
        }
        if set(expected_numeric_authorities) != set(
            MODEL_INPUT_NUMERIC_AUTHORITY_NAMES
        ):
            raise RuntimeError("numeric authority mapping is incomplete")
        numeric_digest = str(numeric_release_binding["sha256"])
        if numeric_digest not in validated_numeric_releases:
            validated_numeric_releases[numeric_digest] = (
                validate_numeric_gate_release(
                    Path(str(numeric_release_binding["path"])),
                    expected_code_revision=expected_code_revision,
                    expected_authority_bindings=expected_numeric_authorities,
                )
            )
        numeric_release_payload = validated_numeric_releases[numeric_digest]
        if (
            dict(numeric_release)
            != _numeric_gate_projection(
                numeric_release_payload,
                binding=numeric_release_binding,
            )
            or numeric_release_payload.get("authority_bindings")
            != expected_numeric_authorities
            or numeric_release_payload.get("native_freshness")
            != expected_freshness
            or numeric_release_payload.get("code_revision")
            != expected_code_revision
            or numeric_release_payload.get("passed") is not True
        ):
            raise ValueError(
                f"fold {fold} numeric release projection differs from "
                "the hash-bound release file"
            )
        registry_key = (
            str(normalized_bindings["registry"]["sha256"]),
            str(normalized_bindings["registry_summary"]["sha256"]),
        )
        if registry_key not in loaded_registries:
            loaded_registries[registry_key] = load_fullpool_ssl_registry(
                registry_path=Path(
                    str(normalized_bindings["registry"]["path"])
                ),
                summary_path=Path(
                    str(normalized_bindings["registry_summary"]["path"])
                ),
            )
        registry, _registry_summary = loaded_registries[registry_key]
        selected, expected_selection = select_fullpool_ssl_fold(
            registry,
            held_out_fold=fold,
            max_rows=None,
            required_observation_ids=(),
        )
        if dict(selection) != expected_selection:
            raise ValueError(
                f"fold {fold} selection differs from the hash-bound registry"
            )
        epoch_coverage = _expected_epoch_coverage(len(selected))
        if (
            contract.get("epoch_coverage") != epoch_coverage
            or summary.get("epoch_coverage") != epoch_coverage
            or summary.get("selection_audit") != selection
        ):
            raise ValueError(f"fold {fold} epoch/selection coverage differs")
        _exact_int(
            summary.get("global_step"),
            PRODUCTION_EPOCHS
            * epoch_coverage["optimizer_steps_per_epoch"],
            context=f"fold {fold} summary global step",
        )
        native_records = contract.get("native_files")
        if not isinstance(native_records, list) or not native_records:
            raise ValueError(f"fold {fold} lacks selected native files")
        expected_native_records = _verify_native_files(
            selected,
            expected_git_sha=expected_code_revision,
        )
        if native_records != expected_native_records:
            raise ValueError(
                f"fold {fold} native files differ from registry selection"
            )
        native_bindings = [
            _validate_native_record(
                value,
                expected_code_revision=expected_code_revision,
                index=index,
            )
            for index, value in enumerate(native_records)
        ]
        if len({item["path"] for item in native_bindings}) != len(
            native_bindings
        ):
            raise ValueError(f"fold {fold} repeats a native HDF5")

        history_columns, history_rows = _read_history(
            history_path,
            expected_n_observations=epoch_coverage[
                "observations_per_epoch"
            ],
        )
        tensor_counts = _validate_checkpoint(
            checkpoint_path,
            fold=fold,
            contract_sha256=contract_binding["sha256"],
            selection=selection,
            expected_freshness=expected_freshness,
            expected_numeric_release=numeric_release,
            model_config=contract["model_config"],
            projection_architecture=contract[
                "projection_architecture"
            ],
            optimizer_config=contract["optimizer_config"],
            epoch_coverage=epoch_coverage,
            history_columns=history_columns,
            history_rows=history_rows,
        )
        common_observed = {
            "run_id": contract["run_id"],
            "encoder_name": contract["encoder_name"],
            "checkpoint_namespace": contract["checkpoint_namespace"],
            "numeric_gate_release": dict(numeric_release),
            "native_freshness": dict(expected_freshness),
            "training_authority": dict(authority),
            "authority_bindings": normalized_bindings,
            "registry_path": contract.get("registry_path"),
            "registry_sha256": contract.get("registry_sha256"),
            "registry_summary_path": contract.get("registry_summary_path"),
            "registry_summary_sha256": contract.get(
                "registry_summary_sha256"
            ),
            "model_config": contract.get("model_config"),
            "augmentation_config": contract.get("augmentation_config"),
            "vicreg_config": contract.get("vicreg_config"),
            "optimizer_config": contract.get("optimizer_config"),
            "projection_architecture": contract.get(
                "projection_architecture"
            ),
        }
        if common is None:
            common = common_observed
        elif common_observed != common:
            raise ValueError("five folds do not share one authority/run binding")
        fold_records.append(
            {
                "fold": fold,
                "seed": seed,
                "completed_epochs": PRODUCTION_EPOCHS,
                "global_step": int(summary["global_step"]),
                "selection_audit": dict(selection),
                "epoch_coverage": epoch_coverage,
                "summary": summary_binding,
                "run_contract": contract_binding,
                "history": history_binding,
                "checkpoint": checkpoint_binding,
                "tensor_counts": tensor_counts,
                "native_h5_count": len(native_bindings),
            }
        )
    assert common is not None
    return {
        "passed": True,
        "schema_version": COMPLETION_RELEASE_SCHEMA,
        "release_binding": COMPLETION_RELEASE_BINDING,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "model_namespace": FULLPOOL_SSL_MODEL_NAMESPACE,
        "expected_code_revision": expected_code_revision,
        "numeric_release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "numeric_gate_release": common["numeric_gate_release"],
        "native_freshness": expected_freshness,
        "authority_bindings": common["authority_bindings"],
        "training_hyperparameters": {
            "folds": list(range(FULLPOOL_SSL_N_FOLDS)),
            "epochs_per_fold": PRODUCTION_EPOCHS,
            "batch_size": PRODUCTION_BATCH_SIZE,
            "workers": PRODUCTION_WORKERS,
            "base_seed": FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
            "fold_seeds": sorted(seeds),
            "learning_rate": FULLPOOL_SSL_LEARNING_RATE,
            "weight_decay": FULLPOOL_SSL_WEIGHT_DECAY,
            "optimizer": common["optimizer_config"],
            "checkpoint_every": PRODUCTION_CHECKPOINT_EVERY,
            "require_cuda": True,
            "max_rows": None,
        },
        "counts": {
            "folds": FULLPOOL_SSL_N_FOLDS,
            "completed_epochs": PRODUCTION_EPOCHS * FULLPOOL_SSL_N_FOLDS,
        },
        "folds": fold_records,
    }


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publish_immutable_bytes(path: Path, payload: bytes) -> None:
    """Atomically create one immutable file, accepting only identical reuse."""

    output = Path(path)
    if output.is_symlink():
        raise RuntimeError(f"immutable output cannot be a symlink: {output}")
    if output.exists():
        if not output.is_file() or output.read_bytes() != payload:
            raise FileExistsError(
                f"immutable output exists with different bytes: {output}"
            )
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent,
        prefix=f".{output.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, output)
        except FileExistsError:
            if (
                output.is_symlink()
                or not output.is_file()
                or output.read_bytes() != payload
            ):
                raise FileExistsError(
                    "immutable output raced with different bytes: "
                    f"{output}"
                ) from None
        _fsync_directory(output.parent)
    finally:
        temporary.unlink(missing_ok=True)


def write_teacher_ssl_fullpool_completion_release(
    *,
    model_root: Path,
    expected_code_revision: str,
    output_path: Path,
) -> dict[str, Any]:
    """Publish one deterministic immutable completion release."""

    release = validate_teacher_ssl_fullpool_training(
        model_root=model_root,
        expected_code_revision=expected_code_revision,
    )
    payload = (
        json.dumps(release, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    output = Path(
        os.path.abspath(Path(output_path).expanduser())
    )
    _publish_immutable_bytes(output, payload)
    digest = _sha256(output)
    sidecar = Path(str(output) + ".sha256")
    sidecar_payload = f"{digest}  {output.name}\n".encode("ascii")
    _publish_immutable_bytes(sidecar, sidecar_payload)
    if _read_json(output, context="completion release") != release:
        raise RuntimeError("published completion release failed readback")
    if _sha256(output) != digest:
        raise RuntimeError("published completion release changed")
    return {
        **release,
        "completion_release": _metadata(output),
        "completion_release_sha256_sidecar": _metadata(sidecar),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-root", type=Path, required=True)
    parser.add_argument("--expected-code-revision", required=True)
    parser.add_argument("--output-release", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    result = write_teacher_ssl_fullpool_completion_release(
        model_root=args.model_root,
        expected_code_revision=args.expected_code_revision,
        output_path=args.output_release,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
