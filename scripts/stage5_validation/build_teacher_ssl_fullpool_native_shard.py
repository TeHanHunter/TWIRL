#!/usr/bin/env python3
"""Build one checksum-bound real-only full-pool native HDF5 shard."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    build_full_pool_native_shard,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--frozen-pool", type=Path, required=True)
    parser.add_argument("--frozen-pool-summary", type=Path, required=True)
    parser.add_argument("--sector-allowlist", type=Path, required=True)
    parser.add_argument("--raw-source-h5", type=Path, required=True)
    parser.add_argument("--compact-adp-h5", type=Path, required=True)
    parser.add_argument("--cadence-reference-table", type=Path, required=True)
    parser.add_argument(
        "--cadence-reference-manifest", type=Path, required=True
    )
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument("--n-periods", type=int, default=4096)
    parser.add_argument(
        "--orbitid-policy",
        choices=("strict", "reference_by_cadence"),
        required=True,
    )
    args = parser.parse_args()
    summary = build_full_pool_native_shard(
        sector=args.sector,
        pool_path=args.frozen_pool,
        pool_summary_path=args.frozen_pool_summary,
        allowlist_path=args.sector_allowlist,
        raw_source_h5=args.raw_source_h5,
        compact_adp_h5=args.compact_adp_h5,
        cadence_reference_table=args.cadence_reference_table,
        cadence_reference_manifest=args.cadence_reference_manifest,
        out_h5=args.out_h5,
        shard_index=args.shard_index,
        n_shards=args.n_shards,
        n_periods=args.n_periods,
        orbitid_policy=args.orbitid_policy,
    )
    summary_path = args.out_h5.with_suffix(".summary.json")
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
