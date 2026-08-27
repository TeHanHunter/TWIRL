#!/usr/bin/env python3
"""Build a checksum-bound, label-blind FM0 later-sector inventory."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from twirl.models.fm0.temporal_admission import build_later_sector_inventory

_SHA256_RE = re.compile(r"[0-9a-fA-F]{64}")


def _parse_sector_receipt(value: str) -> tuple[int, Path, str]:
    """Parse ``SECTOR=SHA256=PATH`` without inferring either binding."""

    fields = value.split("=", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError(
            "sector receipt must have the form SECTOR=SHA256=PATH"
        )
    raw_sector, raw_sha256, raw_path = fields
    try:
        sector = int(raw_sector)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("receipt sector must be an integer") from exc
    if sector < 56:
        raise argparse.ArgumentTypeError(
            "FM later-sector receipts require Sector >= 56"
        )
    if _SHA256_RE.fullmatch(raw_sha256) is None:
        raise argparse.ArgumentTypeError(
            "receipt SHA-256 must contain exactly 64 hexadecimal characters"
        )
    if not raw_path:
        raise argparse.ArgumentTypeError("receipt path must not be empty")
    return sector, Path(raw_path), raw_sha256.lower()


def _sha256(value: str) -> str:
    if _SHA256_RE.fullmatch(value) is None:
        raise argparse.ArgumentTypeError(
            "SHA-256 must contain exactly 64 hexadecimal characters"
        )
    return value.lower()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Frozen later-sector admission policy YAML.",
    )
    parser.add_argument(
        "--config-sha256",
        type=_sha256,
        required=True,
        help="Expected SHA-256 of the frozen later-sector policy YAML.",
    )
    parser.add_argument(
        "--ordered-sector",
        action="append",
        type=int,
        required=True,
        help=(
            "Candidate sector in chronological order. Repeat once per sector; "
            "the order is authoritative and is never inferred."
        ),
    )
    parser.add_argument(
        "--sector-receipt",
        action="append",
        type=_parse_sector_receipt,
        required=True,
        metavar="SECTOR=SHA256=PATH",
        help=(
            "Expected SHA-256 and path for one admission receipt. Repeat to "
            "cover the ordered sectors exactly."
        ),
    )
    parser.add_argument(
        "--selected-sector",
        action="append",
        type=int,
        required=True,
        help=(
            "Sector proposed for the label-blind later panel. Repeat in the "
            "exact proposed order; this does not itself freeze a panel."
        ),
    )
    parser.add_argument(
        "--baseline-manifest",
        type=Path,
        required=True,
        help="Frozen S56-S64 observation/component manifest.",
    )
    parser.add_argument(
        "--baseline-manifest-sha256",
        type=_sha256,
        required=True,
    )
    parser.add_argument(
        "--baseline-selection",
        type=Path,
        required=True,
        help=(
            "Frozen S56-S64 label-blind corpus selection with sector and "
            "detector placement."
        ),
    )
    parser.add_argument(
        "--baseline-selection-sha256",
        type=_sha256,
        required=True,
    )
    parser.add_argument(
        "--baseline-alias-authority",
        type=Path,
        required=True,
        help="Frozen S56-S64 edge-only Gaia/TIC alias authority.",
    )
    parser.add_argument(
        "--baseline-alias-authority-sha256",
        type=_sha256,
        required=True,
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    ordered_sectors = list(args.ordered_sector)
    if len(set(ordered_sectors)) != len(ordered_sectors):
        parser.error("--ordered-sector values must be unique")

    receipts: dict[int, tuple[Path, str]] = {}
    for sector, path, digest in args.sector_receipt:
        if sector in receipts:
            parser.error(f"duplicate --sector-receipt for S{sector}")
        receipts[sector] = (path, digest)
    if set(receipts) != set(ordered_sectors):
        parser.error(
            "--sector-receipt values must cover --ordered-sector values exactly"
        )

    selected_sectors = list(args.selected_sector)
    if len(set(selected_sectors)) != len(selected_sectors):
        parser.error("--selected-sector values must be unique")
    if not set(selected_sectors).issubset(ordered_sectors):
        parser.error("--selected-sector values must be drawn from --ordered-sector")

    ordered_sector_receipts = [
        (sector, receipts[sector][0], receipts[sector][1]) for sector in ordered_sectors
    ]
    summary = build_later_sector_inventory(
        config_path=args.config,
        expected_config_sha256=args.config_sha256,
        ordered_sector_receipts=ordered_sector_receipts,
        selected_sectors=selected_sectors,
        baseline_manifest_path=args.baseline_manifest,
        baseline_manifest_sha256=args.baseline_manifest_sha256,
        baseline_selection_path=args.baseline_selection,
        baseline_selection_sha256=args.baseline_selection_sha256,
        baseline_alias_authority_path=args.baseline_alias_authority,
        baseline_alias_authority_sha256=args.baseline_alias_authority_sha256,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
