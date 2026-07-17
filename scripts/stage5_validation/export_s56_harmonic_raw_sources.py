#!/usr/bin/env python3
"""Export the compact S56 raw-flux/error host subset from the PDO TGLC tree."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import export_tglc_raw_sources  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--orbits", default="119,120")
    parser.add_argument("--compact-adp-h5", type=Path)
    args = parser.parse_args()
    orbits = tuple(int(value) for value in args.orbits.split(",") if value.strip())
    summary = export_tglc_raw_sources(
        training_table=args.training_table,
        raw_root=args.raw_root,
        out_h5=args.out_h5,
        orbits=orbits,
        compact_adp_h5=args.compact_adp_h5,
    )
    summary_path = args.out_h5.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
