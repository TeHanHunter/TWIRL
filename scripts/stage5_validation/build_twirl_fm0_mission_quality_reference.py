#!/usr/bin/env python3
"""Build one immutable provider-neutral FM cadence-quality reference."""

from __future__ import annotations

import argparse
import json
import signal
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.mission_quality_reference import (
    build_mission_quality_reference,
)
from twirl.models.fm0.registry import FM0ContractError


def _interrupt(signum: int, _frame: object) -> None:
    raise InterruptedError(f"received {signal.Signals(signum).name}")


def _install_signal_handlers() -> None:
    signal.signal(signal.SIGTERM, _interrupt)
    signal.signal(signal.SIGINT, _interrupt)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quality-transfer-root",
        type=Path,
        required=True,
        help="Immutable S66-S77 mission-quality transfer bundle.",
    )
    parser.add_argument(
        "--quality-transfer-manifest-sha256",
        required=True,
        help="Controller-authorized SHA-256 of the transfer manifest.",
    )
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--producer-git-sha", required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = build_mission_quality_reference(
        quality_transfer_root=args.quality_transfer_root,
        expected_quality_transfer_manifest_sha256=(
            args.quality_transfer_manifest_sha256
        ),
        sector=args.sector,
        output_dir=args.output_dir,
        producer_git_sha=args.producer_git_sha,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    _install_signal_handlers()
    try:
        raise SystemExit(main())
    except (FM0ContractError, OSError, ValueError) as exc:
        print(f"[fm0-mission-quality-reference] ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
