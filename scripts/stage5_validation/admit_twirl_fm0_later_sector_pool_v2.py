#!/usr/bin/env python3
"""Admit the exact deferred S66--S77 TWIRL-FM preparation pool."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from twirl.models.fm0.later_sector_admission_v2 import (
    PREPARATION_SECTORS,
    ReceiptBinding,
    SectorPreparationBindings,
    admit_preparation_pool,
)

_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


def _sha256(value: str) -> str:
    if _SHA256.fullmatch(value) is None:
        raise argparse.ArgumentTypeError("SHA-256 must contain 64 hexadecimal digits")
    return value.lower()


def _git_sha(value: str) -> str:
    if _GIT_SHA.fullmatch(value) is None:
        raise argparse.ArgumentTypeError(
            "producer Git SHA must contain 40 lowercase hexadecimal digits"
        )
    return value


def _sector_receipt(value: str) -> tuple[int, ReceiptBinding]:
    fields = value.split("=", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError(
            "receipt must have the form SECTOR=SHA256=PATH"
        )
    raw_sector, raw_sha, raw_path = fields
    try:
        sector = int(raw_sector)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("receipt sector must be an integer") from exc
    if sector not in PREPARATION_SECTORS:
        raise argparse.ArgumentTypeError("receipt sector must be in S66--S77")
    if not raw_path:
        raise argparse.ArgumentTypeError("receipt path must not be empty")
    return sector, ReceiptBinding(path=Path(raw_path), sha256=_sha256(raw_sha))


def _indexed(
    parser: argparse.ArgumentParser,
    values: list[tuple[int, ReceiptBinding]],
    *,
    label: str,
) -> dict[int, ReceiptBinding]:
    result: dict[int, ReceiptBinding] = {}
    for sector, binding in values:
        if sector in result:
            parser.error(f"duplicate {label} receipt for S{sector}")
        result[sector] = binding
    if tuple(sorted(result)) != PREPARATION_SECTORS:
        parser.error(f"{label} receipts must cover every sector S66--S77 exactly")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--config-sha256", type=_sha256, required=True)
    parser.add_argument("--exclusion-ledger", type=Path, required=True)
    parser.add_argument("--exclusion-ledger-sha256", type=_sha256, required=True)
    parser.add_argument("--producer-git-sha", type=_git_sha, required=True)
    parser.add_argument(
        "--mission-quality-reference",
        action="append",
        type=_sector_receipt,
        required=True,
        metavar="SECTOR=SHA256=PATH",
    )
    parser.add_argument(
        "--hdf5-quality",
        action="append",
        type=_sector_receipt,
        required=True,
        metavar="SECTOR=SHA256=PATH",
    )
    parser.add_argument(
        "--source-inventory",
        action="append",
        type=_sector_receipt,
        required=True,
        metavar="SECTOR=SHA256=PATH",
    )
    parser.add_argument(
        "--six-view-release",
        action="append",
        type=_sector_receipt,
        required=True,
        metavar="SECTOR=SHA256=PATH",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    mission = _indexed(
        parser, args.mission_quality_reference, label="mission-quality reference"
    )
    hdf5 = _indexed(parser, args.hdf5_quality, label="HDF5-quality")
    source = _indexed(parser, args.source_inventory, label="source-inventory")
    six_view = _indexed(parser, args.six_view_release, label="six-view release")
    bindings = [
        SectorPreparationBindings(
            sector=sector,
            mission_quality_reference=mission[sector],
            hdf5_quality=hdf5[sector],
            source_inventory=source[sector],
            six_view_release=six_view[sector],
        )
        for sector in PREPARATION_SECTORS
    ]
    receipt, digest = admit_preparation_pool(
        config_path=args.config,
        expected_config_sha256=args.config_sha256,
        exclusion_ledger_path=args.exclusion_ledger,
        expected_exclusion_ledger_sha256=args.exclusion_ledger_sha256,
        ordered_sector_bindings=bindings,
        producer_git_sha=args.producer_git_sha,
        output_path=args.output,
    )
    print(
        json.dumps(
            {
                "output": str(args.output.expanduser().resolve()),
                "receipt_sha256": digest,
                "receipt": receipt,
            },
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
