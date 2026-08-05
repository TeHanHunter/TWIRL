#!/usr/bin/env python3
"""Run the native-only S56--S62 Teacher v4-SSL development pilot."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import yaml


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

SCHEMA_VERSION = "twirl_teacher_ssl_training_config_v1"
RUN_NAME = "teacher_ssl_v1"
RUN_ID = "teacher_ssl_v1_s56_s62_a2v1_current_adp"
MODEL_FACING_NAME = "Teacher v4-SSL"
ENCODER_NAME = "teacher_ssl_v1"
ARCHITECTURE_VERSION = "s56_harmonic_cnn_v1"
PRIMARY_PROFILE = "shape_plus_periodogram_bls"
SECTORS = list(range(56, 63))


def _mapping(
    payload: dict[str, Any],
    name: str,
    failures: list[str],
) -> dict[str, Any]:
    value = payload.get(name)
    if not isinstance(value, dict):
        failures.append(f"{name} must be a mapping")
        return {}
    return value


def _exact_keys(
    value: dict[str, Any],
    expected: set[str],
    name: str,
    failures: list[str],
) -> None:
    actual = set(value)
    missing = sorted(expected - actual)
    unknown = sorted(actual - expected)
    if missing:
        failures.append(f"{name} lacks keys {missing}")
    if unknown:
        failures.append(f"{name} has unknown keys {unknown}")


def _expect(
    mapping: dict[str, Any],
    expected: dict[str, Any],
    name: str,
    failures: list[str],
) -> None:
    _exact_keys(mapping, set(expected), name, failures)
    for key, expected_value in expected.items():
        if mapping.get(key) != expected_value:
            failures.append(
                f"{name}.{key}={mapping.get(key)!r}, "
                f"expected {expected_value!r}"
            )


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Teacher v4-SSL config must be a YAML mapping")

    failures: list[str] = []
    top_level_keys = {
        "schema_version",
        "run_name",
        "run_id",
        "model_facing_name",
        "encoder_name",
        "architecture_version",
        "scope",
        "inputs",
        "baseline",
        "unlabeled_pool",
        "ssl",
        "optimization",
        "fine_tuning",
        "compute",
    }
    _exact_keys(payload, top_level_keys, "config", failures)
    expected_scalars = {
        "schema_version": SCHEMA_VERSION,
        "run_name": RUN_NAME,
        "run_id": RUN_ID,
        "model_facing_name": MODEL_FACING_NAME,
        "encoder_name": ENCODER_NAME,
        "architecture_version": ARCHITECTURE_VERSION,
    }
    for key, expected_value in expected_scalars.items():
        if payload.get(key) != expected_value:
            failures.append(
                f"config.{key}={payload.get(key)!r}, "
                f"expected {expected_value!r}"
            )
    _expect(
        _mapping(payload, "scope", failures),
        {
            "sectors": SECTORS,
            "native_only": True,
            "real_only": True,
            "development_only": True,
            "held_cv_fold_excluded": True,
            "fixed_test_model_use_forbidden": True,
            "prospective_sector_63_model_use_forbidden": True,
            "injections_forbidden_from_ssl": True,
            "automatic_production_promotion": False,
        },
        "scope",
        failures,
    )

    inputs = _mapping(payload, "inputs", failures)
    expected_input_keys = {
        "split_bound_training_table",
        "real_native_registry",
        "real_native_registry_summary",
        "baseline_teacher_v3_root",
    }
    _exact_keys(inputs, expected_input_keys, "inputs", failures)
    for key in sorted(expected_input_keys):
        value = inputs.get(key)
        if not isinstance(value, str) or not value.strip():
            failures.append(f"inputs.{key} must be a non-empty path")

    _expect(
        _mapping(payload, "baseline", failures),
        {
            "run_id": "teacher_v3_s56_s62_a2v1_current_adp",
            "summary_filename": "summary.json",
            "expected_summary_sha256": (
                "0bdcb064a7e67f2304ba58b1e79c462d"
                "aa7cc8aad5ad64d2f57a86c4dae46e99"
            ),
            "require_frozen_release": True,
        },
        "baseline",
        failures,
    )

    _expect(
        _mapping(payload, "unlabeled_pool", failures),
        {
            "labels_ignored_during_ssl": True,
            "uncertain_included": True,
            "wide_transit_like_included": True,
            "scalar_metadata": "absent",
            "group_column": "tic",
        },
        "unlabeled_pool",
        failures,
    )

    ssl = _mapping(payload, "ssl", failures)
    _exact_keys(
        ssl,
        {
            "profile",
            "objective",
            "augmentation_contract",
            "event_window_protected",
            "invariant_views",
            "augmentation_parameters",
            "prohibited_augmentations",
            "vicreg",
        },
        "ssl",
        failures,
    )
    expected_ssl_values = {
        "profile": PRIMARY_PROFILE,
        "objective": "vicreg",
        "augmentation_contract": "event_preserving_v1",
        "event_window_protected": True,
        "invariant_views": [
            "event_protected_mask_dropout",
            "small_valid_flux_noise",
        ],
        "prohibited_augmentations": [
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
        ],
    }
    for key, expected_value in expected_ssl_values.items():
        if ssl.get(key) != expected_value:
            failures.append(
                f"ssl.{key}={ssl.get(key)!r}, expected {expected_value!r}"
            )
    augmentation_parameters = ssl.get("augmentation_parameters")
    if not isinstance(augmentation_parameters, dict):
        failures.append("ssl.augmentation_parameters must be a mapping")
    else:
        _expect(
            augmentation_parameters,
            {
                "harmonic_cadence_dropout_probability": 0.02,
                "periodogram_bin_dropout_probability": 0.02,
                "local_masks_preserved": True,
                "scalar_metadata_zeroed": True,
                "invalid_samples_never_unmasked": True,
            },
            "ssl.augmentation_parameters",
            failures,
        )
    vicreg = ssl.get("vicreg")
    if not isinstance(vicreg, dict):
        failures.append("ssl.vicreg must be a mapping")
    else:
        _expect(
            vicreg,
            {
                "invariance_weight": 25.0,
                "variance_weight": 25.0,
                "covariance_weight": 1.0,
            },
            "ssl.vicreg",
            failures,
        )

    optimization = _mapping(payload, "optimization", failures)
    _expect(
        optimization,
        {
            "ssl_epochs": 20,
            "fine_tune_epochs": 100,
            "batch_size": 32,
            "workers": 8,
            "seed": 560064,
        },
        "optimization",
        failures,
    )
    _expect(
        _mapping(payload, "fine_tuning", failures),
        {
            "label_policy": "uncertain_masked",
            "grouped_development_oof_only": True,
            "held_fold_used_for_early_stopping": True,
            "estimate_status": "matched_development_not_untouched",
            "fixed_test_evaluation_forbidden": True,
        },
        "fine_tuning",
        failures,
    )
    _expect(
        _mapping(payload, "compute", failures),
        {
            "system": "orcd",
            "partition": "pg_mki_aryeh",
            "h200_gpus": 1,
            "approved_ssl_ceiling_h200_gpus": 4,
            "scaling_policy": (
                "one_for_native_pilot_then_parallel_folds_or_seeds"
            ),
        },
        "compute",
        failures,
    )
    if failures:
        raise ValueError(
            "invalid Teacher v4-SSL config: " + "; ".join(failures)
        )
    return payload


def _resolve_input(
    *,
    explicit: Path | None,
    config: dict[str, Any],
    key: str,
    input_root: Path,
) -> Path:
    if explicit is not None:
        return explicit.expanduser().resolve()
    value = config["inputs"][key]
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = input_root / path
    return path.resolve()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--input-root",
        type=Path,
        default=ROOT,
        help="Base for relative input paths in the YAML config.",
    )
    parser.add_argument("--training-table", type=Path)
    parser.add_argument("--native-registry", type=Path)
    parser.add_argument("--native-registry-summary", type=Path)
    parser.add_argument("--baseline-teacher-v3-root", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--epochs", type=int)
    parser.add_argument("--fine-tune-epochs", type=int)
    parser.add_argument("--batch-size", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--allow-cpu", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    config = _load_config(args.config.expanduser().resolve())
    optimization = config["optimization"]

    training_table = _resolve_input(
        explicit=args.training_table,
        config=config,
        key="split_bound_training_table",
        input_root=args.input_root.expanduser().resolve(),
    )
    native_registry = _resolve_input(
        explicit=args.native_registry,
        config=config,
        key="real_native_registry",
        input_root=args.input_root.expanduser().resolve(),
    )
    native_registry_summary = _resolve_input(
        explicit=args.native_registry_summary,
        config=config,
        key="real_native_registry_summary",
        input_root=args.input_root.expanduser().resolve(),
    )
    baseline_teacher_v3_root = _resolve_input(
        explicit=args.baseline_teacher_v3_root,
        config=config,
        key="baseline_teacher_v3_root",
        input_root=args.input_root.expanduser().resolve(),
    )
    out_dir = args.out_dir.expanduser().resolve()
    ssl_epochs = int(
        optimization["ssl_epochs"]
        if args.epochs is None
        else args.epochs
    )
    fine_tune_epochs = int(
        optimization["fine_tune_epochs"]
        if args.fine_tune_epochs is None
        else args.fine_tune_epochs
    )
    batch_size = int(
        optimization["batch_size"]
        if args.batch_size is None
        else args.batch_size
    )
    workers = int(
        optimization["workers"]
        if args.workers is None
        else args.workers
    )
    seed = int(
        optimization["seed"] if args.seed is None else args.seed
    )
    for name, value in (
        ("epochs", ssl_epochs),
        ("fine-tune-epochs", fine_tune_epochs),
        ("batch-size", batch_size),
    ):
        if value <= 0:
            raise ValueError(f"--{name} must be positive")
    if workers < 0:
        raise ValueError("--workers must be non-negative")

    start = {
        "event": "teacher_ssl_v1_start",
        "run_id": RUN_ID,
        "model_facing_name": MODEL_FACING_NAME,
        "encoder_name": ENCODER_NAME,
        "training_table": str(training_table),
        "native_registry": str(native_registry),
        "native_registry_summary": str(native_registry_summary),
        "baseline_teacher_v3_root": str(baseline_teacher_v3_root),
        "expected_baseline_summary_sha256": config["baseline"][
            "expected_summary_sha256"
        ],
        "out_dir": str(out_dir),
        "ssl_epochs": ssl_epochs,
        "fine_tune_epochs": fine_tune_epochs,
        "batch_size": batch_size,
        "workers": workers,
        "seed": seed,
        "require_cuda": not args.allow_cpu,
    }
    print(json.dumps(start, sort_keys=True), flush=True)

    from twirl.vetting.teacher_ssl_training import (  # noqa: PLC0415
        run_teacher_ssl_pilot,
    )

    summary = run_teacher_ssl_pilot(
        training_table=training_table,
        native_registry=native_registry,
        native_registry_summary=native_registry_summary,
        baseline_teacher_v3_root=baseline_teacher_v3_root,
        expected_baseline_summary_sha256=config["baseline"][
            "expected_summary_sha256"
        ],
        out_dir=out_dir,
        ssl_epochs=ssl_epochs,
        fine_tune_epochs=fine_tune_epochs,
        batch_size=batch_size,
        workers=workers,
        seed=seed,
        require_cuda=not args.allow_cpu,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
