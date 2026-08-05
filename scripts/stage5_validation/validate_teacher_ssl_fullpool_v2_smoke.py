#!/usr/bin/env python3
"""Fail-closed validation for the Teacher v4-SSL full-pool v2 smoke."""
from __future__ import annotations

import argparse
from collections.abc import Mapping
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.ssl_full_pool_eligibility import (  # noqa: E402
    PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    PRODUCTION_EXCLUDED_OBSERVATIONS,
    PRODUCTION_FULL_IDENTITY_SHA256,
    PRODUCTION_FULL_OBSERVATIONS,
)
from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    FULL_POOL_NATIVE_CONTRACT_VERSION,
)
from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    FULLPOOL_SSL_CHECKPOINT_SCHEMA,
    FULLPOOL_SSL_RUN_CONTRACT_SCHEMA,
    FULLPOOL_SSL_RUN_ID,
    FULLPOOL_SSL_SELECTION_SCHEMA,
    FULLPOOL_SSL_SUMMARY_SCHEMA,
    FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA,
)


SMOKE_VALIDATION_SCHEMA = "twirl_teacher_ssl_fullpool_smoke_validation_v2"
SMOKE_EPOCHS = 1
SMOKE_MAX_ROWS = 4096
SMOKE_BATCH_SIZE = 64
SMOKE_WORKERS = 8
SMOKE_SEED = 560064
SMOKE_LEARNING_RATE = 3.0e-4
SMOKE_WEIGHT_DECAY = 1.0e-4
SMOKE_GLOBAL_STEPS = SMOKE_MAX_ROWS // SMOKE_BATCH_SIZE
AUTHORITY_NAMES: tuple[str, ...] = (
    "eligibility_exclusions",
    "eligibility_summary",
    "native_registry",
    "native_registry_summary",
    "native_release_summary",
    "registry",
    "registry_summary",
)


def _no_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"JSON contains duplicate key {key!r}")
        value[key] = item
    return value


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON contains non-finite constant {value!r}")


def _read_json(path: Path, *, context: str) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_no_duplicate_json_keys,
            parse_constant=_reject_json_constant,
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{context} is not strict UTF-8 JSON: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a JSON object")
    return value


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact_metadata(path: Path) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve(strict=True)
    before = resolved.stat()
    if not resolved.is_file() or before.st_size <= 0:
        raise ValueError(f"artifact is missing, empty, or not a file: {resolved}")
    digest = _file_sha256(resolved)
    after = resolved.stat()
    before_identity = (
        int(before.st_size),
        int(before.st_mtime_ns),
        int(before.st_dev),
        int(before.st_ino),
    )
    after_identity = (
        int(after.st_size),
        int(after.st_mtime_ns),
        int(after.st_dev),
        int(after.st_ino),
    )
    if before_identity != after_identity:
        raise RuntimeError(f"artifact changed while it was hashed: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": int(after.st_size),
        "sha256": digest,
    }


def _declared_metadata(value: Any, *, context: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != {
        "path",
        "size_bytes",
        "sha256",
    }:
        raise ValueError(f"{context} is not an exact artifact binding")
    path_text = value.get("path")
    digest = value.get("sha256")
    size_bytes = value.get("size_bytes")
    if (
        not isinstance(path_text, str)
        or not path_text.strip()
        or not Path(path_text).is_absolute()
    ):
        raise ValueError(f"{context} has an invalid absolute path")
    if (
        not isinstance(size_bytes, int)
        or isinstance(size_bytes, bool)
        or size_bytes <= 0
    ):
        raise ValueError(f"{context} has an invalid size")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise ValueError(f"{context} has an invalid SHA-256")
    return {
        "path": str(Path(path_text).expanduser().resolve()),
        "size_bytes": size_bytes,
        "sha256": digest,
    }


def _exact_int(value: Any, expected: int, *, context: str) -> None:
    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or value != int(expected)
    ):
        raise ValueError(f"{context} must equal {int(expected)}")


def _expected_partition() -> dict[str, dict[str, Any]]:
    return {
        "full": {
            "count": PRODUCTION_FULL_OBSERVATIONS,
            "identity_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
        },
        "eligible": {
            "count": PRODUCTION_ELIGIBLE_OBSERVATIONS,
            "identity_sha256": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
        },
        "excluded": {
            "count": PRODUCTION_EXCLUDED_OBSERVATIONS,
            "identity_sha256": PRODUCTION_EXCLUDED_IDENTITY_SHA256,
        },
    }


def _validate_partition(value: Any) -> dict[str, dict[str, Any]]:
    expected = _expected_partition()
    if not isinstance(value, Mapping) or set(value) != set(expected):
        raise ValueError("training authority has the wrong production partition")
    for name, expected_item in expected.items():
        item = value.get(name)
        if not isinstance(item, Mapping) or set(item) != {
            "count",
            "identity_sha256",
        }:
            raise ValueError(
                f"training authority production partition {name} is invalid"
            )
        _exact_int(
            item.get("count"),
            expected_item["count"],
            context=f"production partition {name} count",
        )
        if item.get("identity_sha256") != expected_item["identity_sha256"]:
            raise ValueError(
                f"production partition {name} identity SHA-256 differs"
            )
    return expected


def _absolute_declared_path(value: Any, *, context: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} path is missing")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ValueError(f"{context} path must be absolute")
    return path.resolve(strict=True)


def _verify_summary_output(
    summary: Mapping[str, Any],
    *,
    path_field: str,
    sha_field: str,
) -> dict[str, Any]:
    path = _absolute_declared_path(
        summary.get(path_field),
        context=f"smoke summary {path_field}",
    )
    metadata = _artifact_metadata(path)
    if summary.get(sha_field) != metadata["sha256"]:
        raise ValueError(f"smoke summary {path_field} SHA-256 differs")
    return metadata


def _validate_smoke_history(path: Path) -> tuple[list[str], dict[str, float]]:
    """Require exactly one complete finite 4096-row smoke epoch."""

    try:
        reader = csv.DictReader(
            io.StringIO(path.read_text(encoding="utf-8"), newline="")
        )
        columns = reader.fieldnames
        rows = list(reader)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise ValueError("smoke history is not a valid UTF-8 CSV") from exc
    if (
        not isinstance(columns, list)
        or len(columns) != len(set(columns))
        or not {"epoch", "n_observations", "loss"} <= set(columns)
    ):
        raise ValueError("smoke history has an invalid or incomplete schema")
    if len(rows) != 1:
        raise ValueError("smoke history must contain exactly one epoch")
    row = rows[0]
    try:
        epoch = int(row["epoch"])
        n_observations = int(row["n_observations"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("smoke history epoch/count fields are invalid") from exc
    if epoch != SMOKE_EPOCHS or n_observations != SMOKE_MAX_ROWS:
        raise ValueError(
            "smoke history must record epoch 1 over exactly 4096 observations"
        )
    numeric: dict[str, float] = {}
    for name in columns:
        try:
            value = float(row[name])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"smoke history {name} is not numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"smoke history {name} is non-finite")
        numeric[name] = value
    return columns, numeric


def _validate_smoke_checkpoint(
    path: Path,
    *,
    expected_fold: int,
    expected_run_contract_sha256: str,
    expected_selection: Mapping[str, Any],
    history_columns: list[str],
    history_row: Mapping[str, float],
) -> int:
    """Load the bounded checkpoint and reject non-finite learned tensors."""

    try:
        import torch

        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:
        raise ValueError("smoke checkpoint is not loadable") from exc
    if not isinstance(checkpoint, Mapping):
        raise ValueError("smoke checkpoint is not a mapping")
    if checkpoint.get("schema_version") != FULLPOOL_SSL_CHECKPOINT_SCHEMA:
        raise ValueError("smoke checkpoint has the wrong schema")
    if checkpoint.get("run_id") != FULLPOOL_SSL_RUN_ID:
        raise ValueError("smoke checkpoint has the wrong run ID")
    _exact_int(
        checkpoint.get("fold"),
        expected_fold,
        context="smoke checkpoint fold",
    )
    _exact_int(
        checkpoint.get("epochs"),
        SMOKE_EPOCHS,
        context="smoke checkpoint epochs",
    )
    if checkpoint.get("run_contract_sha256") != expected_run_contract_sha256:
        raise ValueError("smoke checkpoint has the wrong run-contract SHA-256")
    if checkpoint.get("selection_audit") != expected_selection:
        raise ValueError("smoke checkpoint selection audit differs")
    history = checkpoint.get("history")
    if not isinstance(history, list) or len(history) != 1:
        raise ValueError("smoke checkpoint history must contain one epoch")
    checkpoint_row = history[0]
    if not isinstance(checkpoint_row, Mapping) or set(checkpoint_row) != set(
        history_columns
    ):
        raise ValueError("smoke checkpoint history schema differs from CSV")
    for name in history_columns:
        try:
            value = float(checkpoint_row[name])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"smoke checkpoint history {name} is not numeric"
            ) from exc
        if not math.isfinite(value) or value != float(history_row[name]):
            raise ValueError(
                f"smoke checkpoint history {name} differs from finite CSV"
            )
    tensor_count = 0
    for state_name in ("encoder_state_dict", "projection_state_dict"):
        state = checkpoint.get(state_name)
        if not isinstance(state, Mapping) or not state:
            raise ValueError(f"smoke checkpoint {state_name} is empty or invalid")
        for tensor_name, tensor in state.items():
            if not isinstance(tensor_name, str) or not torch.is_tensor(tensor):
                raise ValueError(
                    f"smoke checkpoint {state_name} contains a non-tensor"
                )
            if not bool(torch.isfinite(tensor).all().item()):
                raise ValueError(
                    f"smoke checkpoint {state_name}.{tensor_name} is non-finite"
                )
            tensor_count += 1
    return tensor_count


def validate_teacher_ssl_fullpool_smoke(
    *,
    summary_path: Path,
    expected_code_revision: str,
    eligibility_exclusions_path: Path,
    eligibility_summary_path: Path,
    native_registry_path: Path,
    native_registry_summary_path: Path,
    native_release_summary_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    expected_fold: int = 0,
) -> dict[str, Any]:
    """Validate one bounded smoke result against all production authorities."""

    if (
        not isinstance(expected_code_revision, str)
        or len(expected_code_revision) != 40
        or any(
            character not in "0123456789abcdef"
            for character in expected_code_revision
        )
    ):
        raise ValueError("expected code revision must be a lowercase 40-hex Git SHA")
    _exact_int(expected_fold, expected_fold, context="expected fold")
    if expected_fold not in range(5):
        raise ValueError("expected fold must be in [0,4]")

    expected_authorities = {
        "eligibility_exclusions": _artifact_metadata(
            eligibility_exclusions_path
        ),
        "eligibility_summary": _artifact_metadata(eligibility_summary_path),
        "native_registry": _artifact_metadata(native_registry_path),
        "native_registry_summary": _artifact_metadata(
            native_registry_summary_path
        ),
        "native_release_summary": _artifact_metadata(
            native_release_summary_path
        ),
        "registry": _artifact_metadata(registry_path),
        "registry_summary": _artifact_metadata(registry_summary_path),
    }
    if tuple(expected_authorities) != AUTHORITY_NAMES:
        raise AssertionError("internal smoke authority order changed")

    summary_metadata = _artifact_metadata(summary_path)
    summary = _read_json(
        Path(summary_metadata["path"]),
        context="Teacher v4-SSL smoke summary",
    )
    if summary.get("schema_version") != FULLPOOL_SSL_SUMMARY_SCHEMA:
        raise ValueError("smoke summary has the wrong schema")
    if summary.get("run_id") != FULLPOOL_SSL_RUN_ID:
        raise ValueError("smoke summary has the wrong run ID")
    _exact_int(summary.get("fold"), expected_fold, context="smoke summary fold")
    _exact_int(
        summary.get("completed_epochs"),
        SMOKE_EPOCHS,
        context="smoke completed epochs",
    )
    if (
        summary.get("labels_loaded") is not False
        or summary.get("fixed_test_status")
        != "host_excluded_tensors_not_constructed"
        or summary.get("prospective_status")
        != "host_excluded_tensors_not_constructed"
        or summary.get("automatic_production_promotion") is not False
    ):
        raise ValueError(
            "smoke summary label/fixed-test/S63 isolation gates failed"
        )
    _exact_int(
        summary.get("global_step"),
        SMOKE_GLOBAL_STEPS,
        context="smoke global steps",
    )

    run_contract_path = _absolute_declared_path(
        summary.get("run_contract"),
        context="smoke run contract",
    )
    run_contract_metadata = _artifact_metadata(run_contract_path)
    if summary.get("run_contract_sha256") != run_contract_metadata["sha256"]:
        raise ValueError("smoke run-contract SHA-256 differs from its bytes")
    contract = _read_json(
        run_contract_path,
        context="Teacher v4-SSL smoke run contract",
    )
    if contract.get("schema_version") != FULLPOOL_SSL_RUN_CONTRACT_SCHEMA:
        raise ValueError("smoke run contract has the wrong schema")
    if contract.get("run_id") != FULLPOOL_SSL_RUN_ID:
        raise ValueError("smoke run contract has the wrong run ID")
    if contract.get("code_revision") != expected_code_revision:
        raise ValueError("smoke run contract has the wrong code revision")
    _exact_int(contract.get("fold"), expected_fold, context="smoke contract fold")
    _exact_int(contract.get("epochs"), SMOKE_EPOCHS, context="smoke epochs")
    _exact_int(
        contract.get("max_rows"),
        SMOKE_MAX_ROWS,
        context="smoke max rows",
    )
    _exact_int(
        contract.get("checkpoint_every"),
        1,
        context="smoke checkpoint interval",
    )
    for name, expected in (
        ("batch_size", SMOKE_BATCH_SIZE),
        ("workers", SMOKE_WORKERS),
        ("seed", SMOKE_SEED),
    ):
        _exact_int(
            contract.get(name),
            expected,
            context=f"smoke {name.replace('_', ' ')}",
        )
    for name, expected in (
        ("learning_rate", SMOKE_LEARNING_RATE),
        ("weight_decay", SMOKE_WEIGHT_DECAY),
    ):
        value = contract.get(name)
        if (
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or float(value) != expected
        ):
            raise ValueError(f"smoke {name.replace('_', ' ')} must equal {expected}")
    if contract.get("require_cuda") is not True:
        raise ValueError("smoke did not require CUDA")
    if (
        contract.get("labels_loaded") is not False
        or contract.get("fixed_test_tensors_constructed") is not False
        or contract.get("prospective_sector_tensors_constructed") is not False
        or contract.get("embedding_export") is not False
        or contract.get("neighbor_probe") is not False
    ):
        raise ValueError(
            "smoke contract label/injection/fixed-test/S63 isolation failed"
        )

    authority = contract.get("training_authority")
    if not isinstance(authority, Mapping):
        raise ValueError("smoke run contract lacks training authority")
    if authority.get("schema_version") != FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA:
        raise ValueError("smoke training authority has the wrong schema")
    if authority.get("production_lock_passed") is not True:
        raise ValueError("smoke training authority production lock failed")
    if authority.get("source_provenance_verified") is not True:
        raise ValueError("smoke training authority source provenance failed")
    partition = _validate_partition(authority.get("partition"))
    mapping_sha = authority.get("native_mapping_sha256")
    if (
        not isinstance(mapping_sha, str)
        or len(mapping_sha) != 64
        or any(character not in "0123456789abcdef" for character in mapping_sha)
    ):
        raise ValueError("smoke training authority native mapping SHA is invalid")
    declared_authorities = authority.get("authority_bindings")
    if not isinstance(declared_authorities, Mapping) or set(
        declared_authorities
    ) != set(AUTHORITY_NAMES):
        raise ValueError(
            "smoke training authority does not bind exactly seven artifacts"
        )
    for name in AUTHORITY_NAMES:
        observed = _declared_metadata(
            declared_authorities.get(name),
            context=f"training authority {name}",
        )
        if observed != expected_authorities[name]:
            raise ValueError(
                f"smoke training authority {name} metadata differs"
            )
    for name, contract_path_field, contract_sha_field in (
        ("registry", "registry_path", "registry_sha256"),
        (
            "registry_summary",
            "registry_summary_path",
            "registry_summary_sha256",
        ),
    ):
        if (
            contract.get(contract_path_field) != expected_authorities[name]["path"]
            or contract.get(contract_sha_field)
            != expected_authorities[name]["sha256"]
        ):
            raise ValueError(f"smoke contract {name} top-level binding differs")

    native_flags = (
        "native_hashes_verified_now",
        "native_root_contracts_verified_now",
        "native_group_identities_verified_now",
    )
    if any(contract.get(name) is not True for name in native_flags):
        raise ValueError("smoke did not pass all three native verification gates")
    native_files = contract.get("native_files")
    if not isinstance(native_files, list) or not native_files:
        raise ValueError("smoke run contract has no verified native files")
    selected_native_rows = 0
    native_paths: set[str] = set()
    for index, record in enumerate(native_files):
        context = f"smoke native file {index}"
        if not isinstance(record, Mapping):
            raise ValueError(f"{context} record is invalid")
        if any(
            record.get(name) is not True
            for name in (
                "hash_verified_now",
                "root_contract_verified_now",
                "group_identities_verified_now",
            )
        ):
            raise ValueError(f"{context} verification flags failed")
        if record.get("native_contract_version") != (
            FULL_POOL_NATIVE_CONTRACT_VERSION
        ):
            raise ValueError(f"{context} is not native-v2")
        path_text = record.get("native_h5_path")
        if (
            not isinstance(path_text, str)
            or not path_text
            or not Path(path_text).is_absolute()
            or path_text in native_paths
        ):
            raise ValueError(f"{context} path is invalid or duplicated")
        native_paths.add(path_text)
        sha = record.get("native_h5_sha256")
        if (
            not isinstance(sha, str)
            or len(sha) != 64
            or any(character not in "0123456789abcdef" for character in sha)
        ):
            raise ValueError(f"{context} SHA-256 is invalid")
        size_bytes = record.get("native_h5_size_bytes")
        if (
            not isinstance(size_bytes, int)
            or isinstance(size_bytes, bool)
            or size_bytes <= 0
        ):
            raise ValueError(f"{context} size is invalid")
        row_count = record.get("n_selected_observations")
        if (
            not isinstance(row_count, int)
            or isinstance(row_count, bool)
            or row_count <= 0
        ):
            raise ValueError(f"{context} selected-row count is invalid")
        selected_native_rows += row_count
    if selected_native_rows != SMOKE_MAX_ROWS:
        raise ValueError("smoke native-file rows do not total exactly 4096")

    selection = contract.get("selection_audit")
    if not isinstance(selection, Mapping):
        raise ValueError("smoke run contract lacks selection audit")
    if selection.get("selection_schema_version") != FULLPOOL_SSL_SELECTION_SCHEMA:
        raise ValueError("smoke selection audit has the wrong schema")
    for name, expected in (
        ("held_out_fold", expected_fold),
        ("n_registry_observations", PRODUCTION_FULL_OBSERVATIONS),
        ("n_eligible_observations", PRODUCTION_ELIGIBLE_OBSERVATIONS),
        ("n_selected_observations", SMOKE_MAX_ROWS),
        ("max_rows", SMOKE_MAX_ROWS),
    ):
        _exact_int(
            selection.get(name),
            expected,
            context=f"smoke selection {name}",
        )
    if selection.get("tic_disjoint") != {
        "held_fold_tics": True,
        "fixed_test_tics": True,
        "reserved_prospective_tics": True,
    }:
        raise ValueError("smoke selection leaks held, fixed-test, or S63 hosts")
    if summary.get("selection_audit") != selection:
        raise ValueError("smoke summary and run-contract selections differ")

    checkpoint_metadata = _verify_summary_output(
        summary,
        path_field="checkpoint",
        sha_field="checkpoint_sha256",
    )
    history_metadata = _verify_summary_output(
        summary,
        path_field="history",
        sha_field="history_sha256",
    )
    history_columns, history_row = _validate_smoke_history(
        Path(history_metadata["path"])
    )
    checkpoint_tensor_count = _validate_smoke_checkpoint(
        Path(checkpoint_metadata["path"]),
        expected_fold=expected_fold,
        expected_run_contract_sha256=run_contract_metadata["sha256"],
        expected_selection=selection,
        history_columns=history_columns,
        history_row=history_row,
    )
    if _artifact_metadata(Path(checkpoint_metadata["path"])) != (
        checkpoint_metadata
    ):
        raise RuntimeError("smoke checkpoint changed while it was validated")
    if _artifact_metadata(Path(history_metadata["path"])) != history_metadata:
        raise RuntimeError("smoke history changed while it was validated")
    return {
        "passed": True,
        "schema_version": SMOKE_VALIDATION_SCHEMA,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "code_revision": expected_code_revision,
        "fold": expected_fold,
        "epochs": SMOKE_EPOCHS,
        "max_rows": SMOKE_MAX_ROWS,
        "batch_size": SMOKE_BATCH_SIZE,
        "workers": SMOKE_WORKERS,
        "seed": SMOKE_SEED,
        "global_steps": SMOKE_GLOBAL_STEPS,
        "partition": partition,
        "authority_bindings": expected_authorities,
        "summary": summary_metadata,
        "run_contract": run_contract_metadata,
        "checkpoint": checkpoint_metadata,
        "history": history_metadata,
        "history_columns": history_columns,
        "checkpoint_tensor_count": checkpoint_tensor_count,
        "native_files_verified_during_smoke": len(native_files),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--expected-code-revision", required=True)
    parser.add_argument("--eligibility-exclusions", type=Path, required=True)
    parser.add_argument("--eligibility-summary", type=Path, required=True)
    parser.add_argument("--native-registry", type=Path, required=True)
    parser.add_argument("--native-registry-summary", type=Path, required=True)
    parser.add_argument("--native-release-summary", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--registry-summary", type=Path, required=True)
    parser.add_argument("--expected-fold", type=int, choices=range(5), default=0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    audit = validate_teacher_ssl_fullpool_smoke(
        summary_path=args.summary,
        expected_code_revision=args.expected_code_revision,
        eligibility_exclusions_path=args.eligibility_exclusions,
        eligibility_summary_path=args.eligibility_summary,
        native_registry_path=args.native_registry,
        native_registry_summary_path=args.native_registry_summary,
        native_release_summary_path=args.native_release_summary,
        registry_path=args.registry,
        registry_summary_path=args.registry_summary,
        expected_fold=args.expected_fold,
    )
    print(json.dumps(audit, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
