#!/usr/bin/env python3
"""Run an explicitly synthetic TWIRL-FM0.1 numerical smoke."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import sys

import yaml


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.dataset import SyntheticFM0Config, SyntheticFM0Dataset  # noqa: E402
from twirl.models.fm0.model import (  # noqa: E402
    architecture_for_variant,
    build_fm0_model,
    count_trainable_parameters,
)
from twirl.models.fm0.training import (  # noqa: E402
    FM0OptimizationConfig,
    run_synthetic_training,
    seed_everything,
)
from twirl.models.fm0.validation import (  # noqa: E402
    RUN_CONTRACT_SCHEMA_VERSION,
    RUN_SUMMARY_SCHEMA_VERSION,
    read_json,
    require_clean_git_revision,
    sha256_file,
    validate_frozen_authorities,
    validate_run_release,
    write_json_with_sha256,
    write_sha256_sidecar,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--design", type=Path, default=ROOT / "doc" / "foundation_model_design.md"
    )
    parser.add_argument("--freeze-receipt", type=Path, required=True)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--architecture", choices=("tcn", "conformer"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--synthetic-smoke", action="store_true")
    parser.add_argument("--input-release", type=Path)
    parser.add_argument("--input-release-receipt", type=Path)
    parser.add_argument("--expected-git-sha")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--target-step", "--max-steps", dest="target_step", type=int, default=1)
    parser.add_argument("--micro-batch-windows", type=int, default=2)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--precision", choices=("fp32", "bf16"), default="fp32")
    parser.add_argument("--resume-checkpoint", type=Path)
    return parser


def _load_config(path: Path) -> dict[str, object]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("FM0.1 config must be a YAML mapping")
    if payload.get("schema_version") != "twirl_fm0_1_poc_config_v1":
        raise ValueError("FM0.1 config schema mismatch")
    return payload


def main() -> int:
    args = _parser().parse_args()
    if not args.synthetic_smoke:
        raise SystemExit(
            "This entrypoint currently supports --synthetic-smoke only; the "
            "real input-release loader remains a separate gated integration."
        )
    if args.input_release is not None:
        raise SystemExit(
            "--input-release is intentionally unavailable in the synthetic-only trainer"
        )
    authorities = validate_frozen_authorities(
        design_path=args.design,
        config_path=args.config,
        freeze_receipt_path=args.freeze_receipt,
    )
    config = _load_config(args.config)
    variants = config.get("variants", {})
    if not isinstance(variants, dict) or args.variant not in variants:
        raise ValueError("requested variant is absent from the frozen config")
    optimization_payload = config.get("optimization", {})
    if not isinstance(optimization_payload, dict):
        raise ValueError("frozen optimization config is malformed")
    frozen_seeds = tuple(int(value) for value in optimization_payload.get("seeds", []))
    if args.seed not in frozen_seeds:
        raise ValueError(f"seed must be one of the frozen values {frozen_seeds}")
    resolved_architecture = architecture_for_variant(
        args.variant, development_winner=args.architecture
    )
    optimization = FM0OptimizationConfig(
        learning_rate=float(optimization_payload["learning_rate"]),
        weight_decay=float(optimization_payload["weight_decay"]),
        warmup_steps=int(optimization_payload["warmup_steps"]),
        max_optimizer_steps=int(optimization_payload["max_optimizer_steps"]),
        effective_batch_windows=int(optimization_payload["effective_batch_windows"]),
        huber_delta=float(config["objective"]["reconstruction"]["huber_delta_fractional_flux"]),
        vicreg_total_weight=float(config["objective"]["vicreg"]["total_weight"]),
        vicreg_invariance_weight=float(config["objective"]["vicreg"]["invariance_weight"]),
        vicreg_variance_weight=float(config["objective"]["vicreg"]["variance_weight"]),
        vicreg_covariance_weight=float(config["objective"]["vicreg"]["covariance_weight"]),
    )

    git_sha = None
    if args.expected_git_sha:
        git_sha = require_clean_git_revision(ROOT, args.expected_git_sha)
    input_receipt = None
    receipt_argument = args.input_release_receipt
    if receipt_argument is None and os.environ.get("TWIRL_FM0_INPUT_RECEIPT"):
        receipt_argument = Path(os.environ["TWIRL_FM0_INPUT_RECEIPT"])
    if receipt_argument is not None:
        receipt_path = receipt_argument.resolve(strict=True)
        expected_receipt_hash = os.environ.get("TWIRL_FM0_INPUT_RECEIPT_SHA256")
        observed_receipt_hash = sha256_file(receipt_path)
        if expected_receipt_hash and observed_receipt_hash != expected_receipt_hash:
            raise ValueError("input-release receipt differs from expected environment hash")
        input_receipt = {
            "path": str(receipt_path),
            "sha256": observed_receipt_hash,
            "consumed_by_synthetic_smoke": False,
        }
    output = args.output_dir.resolve()
    if args.resume_checkpoint is None:
        if output.exists():
            raise FileExistsError(f"refusing non-fresh output directory: {output}")
        output.mkdir(parents=True)
    elif not output.is_dir():
        raise FileNotFoundError("resume output directory does not exist")

    run_contract = {
        "schema_version": RUN_CONTRACT_SCHEMA_VERSION,
        "campaign_id": config["campaign_id"],
        "variant": args.variant,
        "architecture": resolved_architecture,
        "seed": args.seed,
        "synthetic_smoke": True,
        "real_data_consumed": False,
        "precision": args.precision,
        "device_request": args.device,
        "target_step": args.target_step,
        "micro_batch_windows": args.micro_batch_windows,
        "authorities": authorities,
        "input_release_receipt_precondition": input_receipt,
        "expected_git_sha": git_sha,
        "optimization": asdict(optimization),
    }
    contract_path = output / "run_contract.json"
    if contract_path.exists():
        if read_json(contract_path) != run_contract:
            raise ValueError("resume run contract differs from existing contract")
        contract_sha = sha256_file(contract_path)
    else:
        contract_sha = write_json_with_sha256(contract_path, run_contract)

    seed_everything(args.seed)
    model = build_fm0_model(
        args.variant,
        development_winner=args.architecture,
        enforce_parameter_budget=True,
    )
    parameter_count = count_trainable_parameters(model)
    dataset_config = SyntheticFM0Config(variant=args.variant, seed=args.seed)
    dataset = SyntheticFM0Dataset(dataset_config)
    result = run_synthetic_training(
        model=model,
        dataset=dataset,
        output_dir=output,
        run_contract=run_contract,
        optimization=optimization,
        target_step=args.target_step,
        micro_batch_windows=args.micro_batch_windows,
        device=args.device,
        precision=args.precision,
        use_vicreg=args.variant == "TWIRL-FM0.1.5",
        resume_checkpoint=args.resume_checkpoint,
    )
    checkpoint_path = output / "checkpoint.pt"
    checkpoint_sha = write_sha256_sidecar(checkpoint_path)
    summary = {
        "schema_version": RUN_SUMMARY_SCHEMA_VERSION,
        "passed": True,
        "synthetic_only": True,
        "scientific_result": False,
        "variant": args.variant,
        "architecture": model.config.architecture,
        "parameter_count": parameter_count,
        "global_step": result["global_step"],
        "final_metrics": result["final_metrics"],
        "precision": result["precision"],
        "device": result["device"],
        "elapsed_seconds_this_invocation": result["elapsed_seconds_this_invocation"],
        "run_contract_sha256": contract_sha,
        "checkpoint_sha256": checkpoint_sha,
    }
    write_json_with_sha256(output / "summary.json", summary)
    validation = validate_run_release(output)
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
