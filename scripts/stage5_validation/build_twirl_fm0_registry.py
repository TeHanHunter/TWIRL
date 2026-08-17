#!/usr/bin/env python3
"""Build the local frozen-contract TWIRL-FM0.1 identity registry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.registry import (
    build_alias_registry,
    build_observation_registry,
    load_frozen_contract,
    read_rows,
    write_registry_release,
)
from twirl.models.fm0.a2v1_adapter import (
    bind_a2v1_hdf5_source_inventory,
    write_a2v1_hdf5_bound_manifest,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aliases", required=True, type=Path, help="Gaia--TIC CSV/JSON")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--observations", type=Path, help="selected A2v1 products CSV/JSON")
    source.add_argument(
        "--a2v1-hdf5-source-inventory",
        type=Path,
        help="unbound checksum inventory used to derive observations and a bound source manifest",
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
    registry = build_alias_registry(read_rows(args.aliases))
    bound_hdf5_manifest = ()
    if args.observations:
        observations = build_observation_registry(read_rows(args.observations), registry)
    elif args.a2v1_hdf5_source_inventory:
        observations, bound_hdf5_manifest = bind_a2v1_hdf5_source_inventory(
            read_rows(args.a2v1_hdf5_source_inventory), registry
        )
    else:
        observations = ()
    summary = write_registry_release(
        args.out_dir, registry, observations, contract=contract
    )
    if bound_hdf5_manifest:
        bound_path = args.out_dir / "a2v1_hdf5_manifest.csv"
        summary["a2v1_hdf5_manifest_sha256"] = write_a2v1_hdf5_bound_manifest(
            bound_path, bound_hdf5_manifest
        )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
