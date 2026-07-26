#!/usr/bin/env python3
"""Train the frozen S56--S62 Teacher v3 release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_cnn import (  # noqa: E402
    MODEL_VERSION,
    HarmonicTrainConfig,
)
from twirl.vetting.teacher_v3_training import (  # noqa: E402
    TEACHER_V3_BASELINE_PROFILE,
    TEACHER_V3_PRIMARY_PROFILE,
    TEACHER_V3_RUN_ID,
    run_teacher_v3_metadata_baseline,
    run_teacher_v3_training,
)


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError("Teacher v3 config must be a YAML mapping")
    expected = {
        "run_name": "teacher_v3",
        "run_id": TEACHER_V3_RUN_ID,
        "architecture_version": MODEL_VERSION,
    }
    failures = [
        f"{name}={payload.get(name)!r}, expected {value!r}"
        for name, value in expected.items()
        if payload.get(name) != value
    ]
    model = payload.get("model", {})
    if not isinstance(model, dict):
        failures.append("model must be a mapping")
    else:
        if model.get("primary_profile") != TEACHER_V3_PRIMARY_PROFILE:
            failures.append("config primary_profile is not frozen")
        if model.get("baseline_profile") != TEACHER_V3_BASELINE_PROFILE:
            failures.append("config baseline_profile is not frozen")
    evaluation = payload.get("evaluation", {})
    if not isinstance(evaluation, dict) or evaluation.get("calibration") != (
        "one_temperature_from_concatenated_development_oof_logits"
    ):
        failures.append("config does not require pooled development OOF calibration")
    if failures:
        raise ValueError("invalid Teacher v3 config: " + "; ".join(failures))
    return payload


def _resolve_input(
    *,
    explicit: Path | None,
    config: dict[str, Any],
    config_key: str,
    input_root: Path,
    required: bool = True,
) -> Path | None:
    if explicit is not None:
        return explicit.expanduser().resolve()
    inputs = config.get("inputs", {})
    value = inputs.get(config_key, "") if isinstance(inputs, dict) else ""
    if not str(value).strip():
        if required:
            raise ValueError(f"Teacher v3 config lacks inputs.{config_key}")
        return None
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
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument(
        "--output-root",
        type=Path,
        help=(
            "Parent output directory. The run ID and mode are appended unless "
            "the path already ends with the run ID."
        ),
    )
    output.add_argument(
        "--out-dir",
        type=Path,
        help="Exact output directory; no run-ID suffix is added.",
    )
    parser.add_argument(
        "--profiles",
        default="",
        help=(
            "Optional explicit profile mode. `metadata_only` starts the "
            "native-independent baseline; otherwise omit for the frozen "
            "metadata + shape_plus_periodogram_bls run."
        ),
    )
    parser.add_argument(
        "--metadata-only-bootstrap",
        action="store_true",
        help=(
            "Train only the genuine seven-sector metadata baseline, without "
            "native HDF5 files or encoder pretraining."
        ),
    )
    parser.add_argument("--epochs", type=int)
    parser.add_argument("--pretrain-epochs", type=int)
    parser.add_argument("--batch-size", type=int)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--learning-rate", type=float)
    parser.add_argument("--weight-decay", type=float)
    parser.add_argument("--patience", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--allow-cpu", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    config = _load_config(args.config.resolve())
    requested_profiles = tuple(
        value.strip()
        for value in str(args.profiles).split(",")
        if value.strip()
    )
    metadata_only = bool(args.metadata_only_bootstrap) or requested_profiles == (
        TEACHER_V3_BASELINE_PROFILE,
    )
    if requested_profiles and not metadata_only:
        expected = {
            TEACHER_V3_BASELINE_PROFILE,
            TEACHER_V3_PRIMARY_PROFILE,
        }
        if set(requested_profiles) != expected or len(requested_profiles) != 2:
            raise ValueError(
                "--profiles must be exactly metadata_only, or the frozen pair "
                "metadata_only,shape_plus_periodogram_bls"
            )
    optimization = config.get("optimization", {})
    if not isinstance(optimization, dict):
        raise ValueError("Teacher v3 config optimization must be a mapping")

    def configured(name: str, override: Any) -> Any:
        return optimization[name] if override is None else override

    train_config = HarmonicTrainConfig(
        epochs=int(configured("epochs", args.epochs)),
        batch_size=int(configured("batch_size", args.batch_size)),
        learning_rate=float(
            configured("learning_rate", args.learning_rate)
        ),
        weight_decay=float(
            configured("weight_decay", args.weight_decay)
        ),
        patience=int(configured("patience", args.patience)),
        seed=int(configured("seed", args.seed)),
    )
    training_table = _resolve_input(
        explicit=args.training_table,
        config=config,
        config_key="split_bound_training_table",
        input_root=args.input_root.resolve(),
    )
    assert training_table is not None
    if args.out_dir is not None:
        out_dir = args.out_dir.expanduser().resolve()
    else:
        output_root = args.output_root.expanduser().resolve()
        run_root = (
            output_root
            if output_root.name == TEACHER_V3_RUN_ID
            else output_root / TEACHER_V3_RUN_ID
        )
        out_dir = run_root / (
            "metadata_only_bootstrap" if metadata_only else "full"
        )
    print(
        json.dumps(
            {
                "event": "teacher_v3_start",
                "mode": (
                    "metadata_only_bootstrap"
                    if metadata_only
                    else "fixed_primary_and_baseline"
                ),
                "run_id": TEACHER_V3_RUN_ID,
                "training_table": str(training_table),
                "out_dir": str(out_dir),
                "train_config": train_config.__dict__,
                "workers": int(args.workers),
                "require_cuda": not args.allow_cpu,
            },
            sort_keys=True,
        ),
        flush=True,
    )
    if metadata_only:
        summary = run_teacher_v3_metadata_baseline(
            training_table=training_table,
            out_dir=out_dir,
            train_config=train_config,
            workers=args.workers,
            require_cuda=not args.allow_cpu,
        )
    else:
        native_registry = _resolve_input(
            explicit=args.native_registry,
            config=config,
            config_key="real_native_registry",
            input_root=args.input_root.resolve(),
        )
        native_registry_summary = _resolve_input(
            explicit=args.native_registry_summary,
            config=config,
            config_key="real_native_registry_summary",
            input_root=args.input_root.resolve(),
        )
        assert native_registry is not None
        assert native_registry_summary is not None
        summary = run_teacher_v3_training(
            training_table=training_table,
            native_registry=native_registry,
            native_registry_summary=native_registry_summary,
            out_dir=out_dir,
            train_config=train_config,
            workers=args.workers,
            pretrain_epochs=int(
                configured(
                    "encoder_pretrain_epochs",
                    args.pretrain_epochs,
                )
            ),
            require_cuda=not args.allow_cpu,
        )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
