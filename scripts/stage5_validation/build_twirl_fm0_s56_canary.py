#!/usr/bin/env python3
"""Build the deterministic, label-blind S56 FM0.1 canary selection."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from astropy.table import Table

from twirl.models.fm0.canary import select_s56_canary, write_s56_canary_selection
from twirl.models.fm0.registry import build_alias_registry, read_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--aliases", required=True, type=Path)
    parser.add_argument("--pdo-root", required=True, type=Path)
    parser.add_argument("--per-detector", type=int, default=4)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    table = Table.read(args.observations)
    sector_rows = table[table["sector"] == 56]
    aliases = build_alias_registry(read_rows(args.aliases))
    selection, selected_aliases = select_s56_canary(
        sector_rows,
        aliases,
        pdo_root=args.pdo_root,
        per_detector=args.per_detector,
    )
    summary = write_s56_canary_selection(args.out_dir, selection, selected_aliases)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
