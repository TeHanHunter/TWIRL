#!/usr/bin/env python3
"""Run the CPU-only FM0.3 cadence-preserving architecture mechanics screen."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.cadence_preserving_mechanics import (
    evaluate_cadence_preserving_mechanics,
)
from twirl.models.fm0.validation import (
    require_clean_git_revision,
    write_json_with_sha256,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-git-sha", required=True)
    args = parser.parse_args()
    sidecar = args.output.with_name(args.output.name + ".sha256")
    if args.output.exists() or sidecar.exists():
        raise FileExistsError("refusing to overwrite FM0.3 mechanics output")
    require_clean_git_revision(ROOT, args.expected_git_sha)
    result = evaluate_cadence_preserving_mechanics(config_path=args.config)
    result["evaluator_git_sha"] = args.expected_git_sha
    write_json_with_sha256(args.output, result)
    os.chmod(args.output, 0o444)
    os.chmod(sidecar, 0o444)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
