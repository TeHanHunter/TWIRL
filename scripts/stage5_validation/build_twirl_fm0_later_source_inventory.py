#!/usr/bin/env python3
"""Build one checksum-bound, label-blind later-sector FM source inventory."""

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

from twirl.models.fm0.later_source_inventory import build_later_source_inventory
from twirl.models.fm0.registry import FM0ContractError


def _interrupt(signum: int, _frame: object) -> None:
    raise InterruptedError(f"received {signal.Signals(signum).name}")


def _install_signal_handlers() -> None:
    signal.signal(signal.SIGTERM, _interrupt)
    signal.signal(signal.SIGINT, _interrupt)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-receipt", type=Path, required=True)
    parser.add_argument("--source-receipt-sha256", required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--expected-orbit", type=int, action="append", required=True)
    parser.add_argument("--expected-target-authority-sha256", required=True)
    parser.add_argument(
        "--target-inventory",
        type=Path,
        required=True,
        help=(
            "Exact checksum-bound frozen corpus-selection CSV. Observation "
            "FITS is metadata-only and cannot publish partitioned source rows."
        ),
    )
    parser.add_argument("--target-inventory-sha256", required=True)
    parser.add_argument("--corpus-summary", type=Path, required=True)
    parser.add_argument("--corpus-summary-sha256", required=True)
    parser.add_argument("--gaia-tic-alias-authority", type=Path)
    parser.add_argument("--gaia-tic-alias-authority-sha256")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--producer-git-sha", required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if len(args.expected_orbit) != 2:
        raise FM0ContractError("--expected-orbit must be supplied exactly twice")
    summary = build_later_source_inventory(
        source_receipt_path=args.source_receipt,
        source_receipt_sha256=args.source_receipt_sha256,
        sector=args.sector,
        expected_orbits=args.expected_orbit,
        expected_target_authority_sha256=args.expected_target_authority_sha256,
        target_inventory_path=args.target_inventory,
        target_inventory_sha256=args.target_inventory_sha256,
        corpus_summary_path=args.corpus_summary,
        corpus_summary_sha256=args.corpus_summary_sha256,
        gaia_tic_alias_authority_path=args.gaia_tic_alias_authority,
        gaia_tic_alias_authority_sha256=args.gaia_tic_alias_authority_sha256,
        output_dir=args.output_dir,
        producer_git_sha=args.producer_git_sha,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    _install_signal_handlers()
    try:
        raise SystemExit(main())
    except (FM0ContractError, OSError, ValueError) as exc:
        print(f"[fm0-later-source] ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
