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

from twirl.vetting.teacher_ssl_fullpool import (
    FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
    run_fullpool_ssl_fold,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--registry-summary", type=Path, required=True)
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
        "--allow-cpu",
        action="store_true",
        help="Testing only; production full-pool folds require CUDA.",
    )
    parser.add_argument(
        "--skip-native-hash-verification",
        action="store_true",
        help=(
            "Trust registry-declared native hashes after checking existence. "
            "Production launch should omit this option."
        ),
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        help="Deterministic bounded smoke only; omit for the full-pool run.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    start = {
        "event": "teacher_ssl_fullpool_fold_start",
        "registry": str(args.registry.expanduser().resolve()),
        "registry_summary": str(args.registry_summary.expanduser().resolve()),
        "out_root": str(args.out_root.expanduser().resolve()),
        "fold": int(args.fold),
        "epochs": int(args.epochs),
        "batch_size": int(args.batch_size),
        "workers": int(args.workers),
        "seed": int(args.seed),
        "resume": bool(args.resume),
        "require_cuda": not args.allow_cpu,
        "verify_native_hashes": not args.skip_native_hash_verification,
        "max_rows": args.max_rows,
    }
    print(json.dumps(start, sort_keys=True), flush=True)
    summary = run_fullpool_ssl_fold(
        registry_path=args.registry,
        registry_summary_path=args.registry_summary,
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
        require_cuda=not args.allow_cpu,
        verify_native_hashes=not args.skip_native_hash_verification,
        max_rows=args.max_rows,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
