#!/usr/bin/env python3
"""Run the fast raw/quality-only FM event-transfer v3 mechanics control."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.event_transfer_mechanics import (
    evaluate_event_transfer_mechanics,
)
from twirl.models.fm0.validation import (
    require_clean_git_revision,
    write_json_with_sha256,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--temporal-panel-dir", type=Path, required=True)
    parser.add_argument("--temporal-panel-receipt-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-git-sha", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    sidecar = args.output.with_name(args.output.name + ".sha256")
    if args.output.exists() or sidecar.exists():
        raise FileExistsError("refusing to overwrite event-transfer mechanics output")
    require_clean_git_revision(ROOT, args.expected_git_sha)
    payload = evaluate_event_transfer_mechanics(
        config_path=args.config,
        temporal_panel_dir=args.temporal_panel_dir,
        temporal_panel_receipt_sha256=args.temporal_panel_receipt_sha256,
    )
    payload["evaluator_git_sha"] = args.expected_git_sha
    write_json_with_sha256(args.output, payload)
    os.chmod(args.output, 0o444)
    os.chmod(sidecar, 0o444)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
