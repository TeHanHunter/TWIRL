#!/usr/bin/env python3
"""Freeze the checksum-bound S56--S62 real-light-curve pool for SSL."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.ssl_full_pool import write_ssl_full_pool


def _sector_path(value: str) -> tuple[int, Path]:
    sector_text, separator, path_text = value.partition("=")
    if not separator or not sector_text.strip() or not path_text.strip():
        raise argparse.ArgumentTypeError(
            "compact HDF5 override must have the form SECTOR=PATH"
        )
    try:
        sector = int(sector_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "compact HDF5 override sector must be an integer"
        ) from exc
    return sector, Path(path_text)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--compact-manifest",
        type=Path,
        action="append",
        required=True,
        help=(
            "One compact export *.manifest.json; repeat exactly once for "
            "each sector S56--S62."
        ),
    )
    parser.add_argument(
        "--compact-h5",
        type=_sector_path,
        action="append",
        default=[],
        metavar="SECTOR=PATH",
        help=(
            "Optional compact-HDF5 path override. By default each manifest's "
            "adjacent .h5 file is used."
        ),
    )
    parser.add_argument("--split-registry", type=Path, required=True)
    parser.add_argument("--s63-reserved-tics", type=Path, required=True)
    parser.add_argument("--s63-reserved-summary", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    overrides: dict[int, Path] = {}
    for sector, path in args.compact_h5:
        if sector in overrides:
            raise ValueError(
                f"duplicate --compact-h5 override for sector {sector}"
            )
        overrides[sector] = path
    summary = write_ssl_full_pool(
        compact_manifest_paths=args.compact_manifest,
        split_registry_path=args.split_registry,
        s63_reserved_tics_path=args.s63_reserved_tics,
        s63_reserved_summary_path=args.s63_reserved_summary,
        out_dir=args.out_dir,
        compact_h5_by_sector=overrides,
        progress=print,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
