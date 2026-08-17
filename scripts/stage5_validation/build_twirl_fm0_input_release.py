#!/usr/bin/env python3
"""Build a local, search-free TWIRL-FM0.1 six-view input release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.input_release import (
    write_a2v1_hdf5_input_release,
    write_input_release,
)
from twirl.models.fm0.registry import load_frozen_contract


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry-dir", required=True, type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--raw-manifest",
        type=Path,
        help="strict-NPZ fixture manifest; never scientific-training eligible",
    )
    source.add_argument(
        "--a2v1-hdf5-manifest",
        type=Path,
        help="checksum-bound accepted/diagnostic A2v1 HDF5 source manifest",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--freeze-receipt", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    contract = load_frozen_contract(
        design_path=args.design,
        config_path=args.config,
        freeze_receipt_path=args.freeze_receipt,
    )
    if args.raw_manifest is not None:
        summary = write_input_release(
            registry_dir=args.registry_dir,
            raw_manifest_path=args.raw_manifest,
            out_dir=args.out_dir,
            contract=contract,
        )
    else:
        summary = write_a2v1_hdf5_input_release(
            registry_dir=args.registry_dir,
            hdf5_manifest_path=args.a2v1_hdf5_manifest,
            out_dir=args.out_dir,
            contract=contract,
        )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
