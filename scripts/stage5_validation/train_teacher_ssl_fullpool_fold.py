#!/usr/bin/env python3
"""Train or resume one fold-local encoder on the broad unlabeled LC pool."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
    FULLPOOL_SSL_MODEL_NAMESPACE,
    run_fullpool_ssl_fold,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eligibility-exclusions", type=Path, required=True)
    parser.add_argument("--eligibility-summary", type=Path, required=True)
    parser.add_argument("--native-registry", type=Path, required=True)
    parser.add_argument("--native-registry-summary", type=Path, required=True)
    parser.add_argument("--native-release-summary", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--registry-summary", type=Path, required=True)
    parser.add_argument("--numeric-gate-release", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--fold", type=int, choices=range(5), required=True)
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=FULLPOOL_SSL_DEFAULT_TRAINING_SEED)
    parser.add_argument("--learning-rate", type=float, default=3.0e-4)
    parser.add_argument("--weight-decay", type=float, default=1.0e-4)
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume only when the existing fold contract is byte-identical.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        help="Deterministic bounded smoke only; omit for the full-pool run.",
    )
    parser.add_argument(
        "--required-observation-id",
        action="append",
        default=None,
        help=(
            "Require one eligible, non-held observation in the selection; "
            "repeat the option for multiple IDs. Intended for bounded smokes."
        ),
    )
    parser.add_argument(
        "--allow-rederived-eligibility",
        action="store_true",
        help="Accept the corrected v4 BLS-derived eligibility authority.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    start = {
        "event": "teacher_ssl_fullpool_fold_start",
        "eligibility_exclusions": str(
            args.eligibility_exclusions.expanduser().resolve()
        ),
        "eligibility_summary": str(
            args.eligibility_summary.expanduser().resolve()
        ),
        "native_registry": str(args.native_registry.expanduser().resolve()),
        "native_registry_summary": str(
            args.native_registry_summary.expanduser().resolve()
        ),
        "native_release_summary": str(
            args.native_release_summary.expanduser().resolve()
        ),
        "registry": str(args.registry.expanduser().resolve()),
        "registry_summary": str(args.registry_summary.expanduser().resolve()),
        "numeric_gate_release": str(
            args.numeric_gate_release.expanduser().resolve()
        ),
        "out_root": str(args.out_root.expanduser().resolve()),
        "fold": int(args.fold),
        "epochs": int(args.epochs),
        "batch_size": int(args.batch_size),
        "workers": int(args.workers),
        "seed": int(args.seed),
        "resume": bool(args.resume),
        "require_cuda": True,
        "model_namespace": FULLPOOL_SSL_MODEL_NAMESPACE,
        "verify_native_hashes": True,
        "verify_numeric_gate_release": True,
        "max_rows": args.max_rows,
        "required_observation_ids": sorted(
            args.required_observation_id or []
        ),
    }
    print(json.dumps(start, sort_keys=True), flush=True)
    summary = run_fullpool_ssl_fold(
        eligibility_exclusions_path=args.eligibility_exclusions,
        eligibility_summary_path=args.eligibility_summary,
        native_registry_path=args.native_registry,
        native_registry_summary_path=args.native_registry_summary,
        native_release_summary_path=args.native_release_summary,
        registry_path=args.registry,
        registry_summary_path=args.registry_summary,
        numeric_gate_release_path=args.numeric_gate_release,
        out_root=args.out_root,
        fold=int(args.fold),
        epochs=int(args.epochs),
        batch_size=int(args.batch_size),
        workers=int(args.workers),
        seed=int(args.seed),
        learning_rate=float(args.learning_rate),
        weight_decay=float(args.weight_decay),
        checkpoint_every=int(args.checkpoint_every),
        resume=bool(args.resume),
        require_cuda=True,
        max_rows=args.max_rows,
        required_observation_ids=args.required_observation_id,
        eligibility_production_lock=not args.allow_rederived_eligibility,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
