#!/usr/bin/env python3
"""Export one sector's compact raw/error host subset for Teacher v3."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import export_tglc_raw_sources  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--compact-adp-h5", type=Path, required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--orbits", type=int, nargs="+", required=True)
    args = parser.parse_args()
    if args.sector < 56 or args.sector > 62:
        raise SystemExit("Teacher-v3 raw export is bounded to sectors 56-62")
    if len(set(args.orbits)) != len(args.orbits) or any(
        orbit <= 0 for orbit in args.orbits
    ):
        raise SystemExit("--orbits must contain distinct positive orbit IDs")
    summary = export_tglc_raw_sources(
        training_table=args.training_table,
        raw_root=args.raw_root,
        out_h5=args.out_h5,
        orbits=tuple(args.orbits),
        compact_adp_h5=args.compact_adp_h5,
    )
    summary["sector"] = int(args.sector)
    summary["orbits"] = [int(value) for value in args.orbits]
    summary_path = args.out_h5.with_suffix(".summary.json")
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
