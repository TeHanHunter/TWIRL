#!/usr/bin/env python3
"""Run the non-authoritative FM0.3 native-cadence collapse preflight."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.cadence_collapse_artifact import (
    evaluate_cadence_collapse_artifacts,
)
from twirl.models.fm0.validation import (
    require_clean_git_revision,
    write_json_with_sha256,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--panel-sha256", required=True)
    parser.add_argument("--step0-checkpoint", type=Path, required=True)
    parser.add_argument("--step0-checkpoint-sha256", required=True)
    parser.add_argument("--step2000-checkpoint", type=Path, required=True)
    parser.add_argument("--step2000-checkpoint-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-git-sha", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    sidecar = args.output.with_name(args.output.name + ".sha256")
    if args.output.exists() or sidecar.exists():
        raise FileExistsError("refusing to overwrite cadence-collapse output")
    require_clean_git_revision(ROOT, args.expected_git_sha)
    payload = evaluate_cadence_collapse_artifacts(
        panel_path=args.panel,
        panel_sha256=args.panel_sha256,
        step0_checkpoint_path=args.step0_checkpoint,
        step0_checkpoint_sha256=args.step0_checkpoint_sha256,
        step2000_checkpoint_path=args.step2000_checkpoint,
        step2000_checkpoint_sha256=args.step2000_checkpoint_sha256,
    )
    payload["evaluator_git_sha"] = args.expected_git_sha
    write_json_with_sha256(args.output, payload)
    os.chmod(args.output, 0o444)
    os.chmod(sidecar, 0o444)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
