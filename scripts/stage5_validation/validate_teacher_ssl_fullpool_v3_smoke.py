#!/usr/bin/env python3
"""Fail-closed validation for the fresh effective-quality ADP H200 smoke."""
from __future__ import annotations

import argparse
from collections.abc import Mapping
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.ssl_full_pool_numeric import (  # noqa: E402
    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
    MODEL_INPUT_NUMERIC_RELEASE_BINDING,
    validate_numeric_native_freshness,
)
from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
    FULLPOOL_SSL_MODEL_NAMESPACE,
    FULLPOOL_SSL_RUN_ID,
    FULLPOOL_SSL_SMOKE_FOLD,
)


_BASE_PATH = Path(__file__).with_name(
    "validate_teacher_ssl_fullpool_v2_smoke.py"
)
_BASE_SPEC = importlib.util.spec_from_file_location(
    "validate_teacher_ssl_fullpool_v2_smoke_preserved",
    _BASE_PATH,
)
if _BASE_SPEC is None or _BASE_SPEC.loader is None:
    raise RuntimeError("cannot load the preserved v2 smoke validator")
BASE = importlib.util.module_from_spec(_BASE_SPEC)
_BASE_SPEC.loader.exec_module(BASE)


SMOKE_VALIDATION_SCHEMA = "twirl_teacher_ssl_fullpool_smoke_validation_v3"
AUTHORITY_NAMES = BASE.AUTHORITY_NAMES
SMOKE_REQUIRED_OBSERVATION_IDS = BASE.SMOKE_REQUIRED_OBSERVATION_IDS
SMOKE_EPOCHS = BASE.SMOKE_EPOCHS
SMOKE_MAX_ROWS = BASE.SMOKE_MAX_ROWS
SMOKE_BATCH_SIZE = BASE.SMOKE_BATCH_SIZE
SMOKE_WORKERS = BASE.SMOKE_WORKERS
SMOKE_BASE_SEED = BASE.SMOKE_SEED
SMOKE_SEED = SMOKE_BASE_SEED + 1000 * FULLPOOL_SSL_SMOKE_FOLD
SMOKE_LEARNING_RATE = BASE.SMOKE_LEARNING_RATE
SMOKE_WEIGHT_DECAY = BASE.SMOKE_WEIGHT_DECAY
SMOKE_GLOBAL_STEPS = BASE.SMOKE_GLOBAL_STEPS
SMOKE_FOLD = FULLPOOL_SSL_SMOKE_FOLD
FULLPOOL_SSL_CHECKPOINT_SCHEMA = BASE.FULLPOOL_SSL_CHECKPOINT_SCHEMA


def _expected_numeric_gate_summary(
    release: Mapping[str, Any],
    *,
    binding: Mapping[str, Any],
) -> dict[str, Any]:
    """Match the exact fresh numerical release embedded by training."""

    return {
        "binding": dict(binding),
        "release_binding": release.get("release_binding"),
        "schema_version": release.get("schema_version"),
        "envelope_canonical_sha256": release.get(
            "envelope_canonical_sha256"
        ),
        "counts": dict(release.get("counts", {})),
        "identity_hashes": dict(release.get("identity_hashes", {})),
        "code_revision": release.get("code_revision"),
        "passed": True,
    }


def _require_freshness(
    value: Any,
    *,
    expected_freshness: Mapping[str, Any],
    expected_code_revision: str,
    context: str,
) -> None:
    try:
        observed = validate_numeric_native_freshness(
            value,
            context=context,
            expected_code_revision=expected_code_revision,
        )
    except ValueError as exc:
        raise ValueError(
            f"{context} differs from the validated numeric-release freshness"
        ) from exc
    if observed != dict(expected_freshness):
        raise ValueError(
            f"{context} differs from the validated numeric-release freshness"
        )


def _require_model_bindings(
    value: Mapping[str, Any],
    *,
    expected_freshness: Mapping[str, Any],
    expected_code_revision: str,
    context: str,
) -> None:
    if value.get("numeric_release_binding") != (
        MODEL_INPUT_NUMERIC_RELEASE_BINDING
    ):
        raise ValueError(f"{context} has a stale numeric release binding")
    _require_freshness(
        value.get("native_freshness"),
        expected_freshness=expected_freshness,
        expected_code_revision=expected_code_revision,
        context=f"{context} native freshness",
    )


def _read_checkpoint(path: Path) -> Mapping[str, Any]:
    try:
        import torch

        value = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:
        raise ValueError("fresh smoke checkpoint is not loadable") from exc
    if not isinstance(value, Mapping):
        raise ValueError("fresh smoke checkpoint is not a mapping")
    return value


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
    numeric_gate_release_path: Path,
    expected_fold: int = FULLPOOL_SSL_SMOKE_FOLD,
) -> dict[str, Any]:
    """Validate the exact fold-2/4096-row/one-epoch fresh H200 smoke."""

    preserved = BASE._expected_numeric_gate_summary
    preserved_seed = BASE.SMOKE_SEED
    BASE._expected_numeric_gate_summary = _expected_numeric_gate_summary
    BASE.SMOKE_SEED = SMOKE_SEED
    try:
        audit = BASE.validate_teacher_ssl_fullpool_smoke(
            summary_path=summary_path,
            expected_code_revision=expected_code_revision,
            eligibility_exclusions_path=eligibility_exclusions_path,
            eligibility_summary_path=eligibility_summary_path,
            native_registry_path=native_registry_path,
            native_registry_summary_path=native_registry_summary_path,
            native_release_summary_path=native_release_summary_path,
            registry_path=registry_path,
            registry_summary_path=registry_summary_path,
            numeric_gate_release_path=numeric_gate_release_path,
            expected_fold=expected_fold,
        )
    finally:
        BASE._expected_numeric_gate_summary = preserved
        BASE.SMOKE_SEED = preserved_seed

    numeric_release = BASE._read_json(
        Path(numeric_gate_release_path).expanduser().resolve(strict=True),
        context="validated numeric-gate release",
    )
    expected_freshness = validate_numeric_native_freshness(
        numeric_release.get("native_freshness"),
        context="numeric-gate native freshness",
        expected_code_revision=expected_code_revision,
    )

    summary_resolved = Path(summary_path).expanduser().resolve(strict=True)
    expected_suffix = (
        "model_runs",
        FULLPOOL_SSL_MODEL_NAMESPACE,
        "smoke",
        "one_epoch",
        "encoder_pretraining",
        f"fold_{FULLPOOL_SSL_SMOKE_FOLD}",
        "summary.json",
    )
    if tuple(summary_resolved.parts[-len(expected_suffix) :]) != (
        expected_suffix
    ):
        raise ValueError(
            "smoke summary is outside the fresh effective_quality_adp_v1 "
            "model namespace"
        )

    summary = BASE._read_json(
        summary_resolved,
        context="fresh Teacher v4-SSL smoke summary",
    )
    _require_model_bindings(
        summary,
        expected_freshness=expected_freshness,
        expected_code_revision=expected_code_revision,
        context="smoke summary",
    )
    if (
        summary.get("run_id") != FULLPOOL_SSL_RUN_ID
        or summary.get("checkpoint_namespace")
        != FULLPOOL_SSL_CHECKPOINT_NAMESPACE
    ):
        raise ValueError("smoke summary belongs to a stale model namespace")

    contract_path = Path(str(summary.get("run_contract", ""))).resolve(
        strict=True
    )
    contract = BASE._read_json(
        contract_path,
        context="fresh Teacher v4-SSL smoke run contract",
    )
    _require_model_bindings(
        contract,
        expected_freshness=expected_freshness,
        expected_code_revision=expected_code_revision,
        context="smoke run contract",
    )
    authority = contract.get("training_authority")
    if not isinstance(authority, Mapping):
        raise ValueError("smoke run contract lacks training authority")
    _require_freshness(
        authority.get("native_freshness"),
        expected_freshness=expected_freshness,
        expected_code_revision=expected_code_revision,
        context="smoke training authority native freshness",
    )
    numeric_summary = authority.get("numeric_gate_release")
    if (
        not isinstance(numeric_summary, Mapping)
        or numeric_summary.get("release_binding")
        != MODEL_INPUT_NUMERIC_RELEASE_BINDING
    ):
        raise ValueError(
            "smoke training authority has a stale numeric release"
        )

    checkpoint_path = Path(str(summary.get("checkpoint", ""))).resolve(
        strict=True
    )
    checkpoint = _read_checkpoint(checkpoint_path)
    _require_model_bindings(
        checkpoint,
        expected_freshness=expected_freshness,
        expected_code_revision=expected_code_revision,
        context="smoke checkpoint",
    )
    if (
        checkpoint.get("run_id") != FULLPOOL_SSL_RUN_ID
        or checkpoint.get("checkpoint_namespace")
        != FULLPOOL_SSL_CHECKPOINT_NAMESPACE
    ):
        raise ValueError("smoke checkpoint belongs to a stale namespace")

    for index, record in enumerate(contract.get("native_files", [])):
        if not isinstance(record, Mapping):
            raise ValueError(f"smoke native file {index} is invalid")
        expected = {
            "native_contract_version": (
                MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                    "native_contract_version"
                ]
            ),
            "native_release_binding": (
                MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS["release_binding"]
            ),
            "native_expected_git_sha": expected_code_revision,
            "native_builder_contract_version": (
                MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                    "builder_contract_version"
                ]
            ),
            "native_builder_code_sha256": (
                MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                    "builder_code_sha256"
                ]
            ),
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
        }
        if {name: record.get(name) for name in expected} != expected:
            raise ValueError(
                f"smoke native file {index} lacks fresh ADP provenance"
            )
        native_path = Path(
            str(record.get("native_h5_path", ""))
        ).expanduser().resolve(strict=True)
        if (
            not native_path.is_file()
            or native_path.stat().st_size
            != int(record.get("native_h5_size_bytes", -1))
        ):
            raise ValueError(
                f"smoke native file {index} size binding changed"
            )
        digest = hashlib.sha256()
        with native_path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != record.get("native_h5_sha256"):
            raise ValueError(
                f"smoke native file {index} hash binding changed"
            )
        try:
            import h5py

            with h5py.File(native_path, "r") as h5:
                root_expected = {
                    "contract_version": expected[
                        "native_contract_version"
                    ],
                    "release_binding": expected[
                        "native_release_binding"
                    ],
                    "expected_git_sha": expected[
                        "native_expected_git_sha"
                    ],
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
                    "detrend_quality_source": expected[
                        "native_detrend_quality_source"
                    ],
                }
                observed_root = {
                    name: str(h5.attrs.get(name, ""))
                    for name in root_expected
                }
                if observed_root != root_expected:
                    raise ValueError(
                        f"smoke native file {index} root provenance changed"
                    )
                for name, expected_value in (
                    ("raw_photometry_only", 1),
                    ("compact_adp_photometry_reused", 0),
                    ("compact_adp_flux_reused", 0),
                    (
                        "periodogram_n",
                        MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                            "periodogram_n"
                        ],
                    ),
                ):
                    if int(h5.attrs.get(name, -1)) != expected_value:
                        raise ValueError(
                            f"smoke native file {index} root {name} changed"
                        )
        except OSError as exc:
            raise ValueError(
                f"smoke native file {index} is not loadable HDF5"
            ) from exc

    audit.update(
        {
            "schema_version": SMOKE_VALIDATION_SCHEMA,
            "run_id": FULLPOOL_SSL_RUN_ID,
            "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
            "model_namespace": FULLPOOL_SSL_MODEL_NAMESPACE,
            "numeric_release_binding": (
                MODEL_INPUT_NUMERIC_RELEASE_BINDING
            ),
            "native_freshness": dict(contract["native_freshness"]),
            "fresh_output_root_verified": True,
        }
    )
    return audit


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
    parser.add_argument("--numeric-gate-release", type=Path, required=True)
    parser.add_argument(
        "--expected-fold",
        type=int,
        choices=(FULLPOOL_SSL_SMOKE_FOLD,),
        default=FULLPOOL_SSL_SMOKE_FOLD,
    )
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
        numeric_gate_release_path=args.numeric_gate_release,
        expected_fold=args.expected_fold,
    )
    print(json.dumps(audit, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
