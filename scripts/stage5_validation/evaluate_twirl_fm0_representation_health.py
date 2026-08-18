#!/usr/bin/env python3
"""Evaluate FM0.1 representation health on the development partition only."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.representation_health import (  # noqa: E402
    evaluate_representation_health,
)
from twirl.models.fm0.validation import (  # noqa: E402
    require_clean_git_revision,
    write_json_with_sha256,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--max-components", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--random-control-seed", type=int, default=0)
    parser.add_argument("--expected-git-sha")
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.output.exists():
        raise FileExistsError("refusing to overwrite representation-health output")
    if args.expected_git_sha:
        require_clean_git_revision(ROOT, args.expected_git_sha)
    payload = evaluate_representation_health(
        run_dir=args.run_dir,
        device=args.device,
        max_components=args.max_components,
        batch_size=args.batch_size,
        random_control_seed=args.random_control_seed,
    )
    write_json_with_sha256(args.output, payload)
    os.chmod(args.output, 0o444)
    os.chmod(args.output.with_name(args.output.name + ".sha256"), 0o444)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
