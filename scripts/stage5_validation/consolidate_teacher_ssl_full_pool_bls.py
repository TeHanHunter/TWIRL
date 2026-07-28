#!/usr/bin/env python3
"""Validate and freeze the global S56--S62 full-pool BLS table."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.ssl_full_pool_bls import write_global_full_pool_bls


def _sector_path(value: str) -> tuple[int, Path]:
    sector_text, separator, path_text = value.partition("=")
    if not separator or not sector_text.strip() or not path_text.strip():
        raise argparse.ArgumentTypeError(
            "sector path must have the form SECTOR=PATH"
        )
    try:
        sector = int(sector_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "sector path must begin with an integer"
        ) from exc
    return sector, Path(path_text)


def _unique_mapping(
    values: list[tuple[int, Path]],
    *,
    option: str,
) -> dict[int, Path]:
    output: dict[int, Path] = {}
    for sector, path in values:
        if sector in output:
            raise ValueError(f"duplicate {option} for sector {sector}")
        output[sector] = path
    return output


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-summary", type=Path, required=True)
    parser.add_argument(
        "--sector-bls",
        type=_sector_path,
        action="append",
        required=True,
        metavar="SECTOR=PARQUET",
        help="Repeat exactly once for each sector S56--S62.",
    )
    parser.add_argument(
        "--sector-summary",
        type=_sector_path,
        action="append",
        default=[],
        metavar="SECTOR=JSON",
        help=(
            "Optional merged-summary override. The default is the sector "
            "Parquet path with suffix .summary.json."
        ),
    )
    parser.add_argument("--out-parquet", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = write_global_full_pool_bls(
        pool_summary_path=args.pool_summary,
        sector_bls_paths=_unique_mapping(
            args.sector_bls,
            option="--sector-bls",
        ),
        sector_summary_paths=_unique_mapping(
            args.sector_summary,
            option="--sector-summary",
        ),
        out_parquet_path=args.out_parquet,
        out_summary_path=args.out_summary,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
