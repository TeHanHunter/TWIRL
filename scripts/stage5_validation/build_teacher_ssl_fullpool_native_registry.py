#!/usr/bin/env python3
"""Freeze exact full-pool native shards into an observation-keyed registry."""
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
    write_full_pool_native_registry,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-pool", type=Path, required=True)
    parser.add_argument("--frozen-pool-summary", type=Path, required=True)
    parser.add_argument("--eligibility-exclusions", type=Path, required=True)
    parser.add_argument("--eligibility-summary", type=Path, required=True)
    parser.add_argument(
        "--native-shards", type=Path, nargs="+", required=True
    )
    parser.add_argument(
        "--native-shard-summaries", type=Path, nargs="+", required=True
    )
    parser.add_argument("--source-out", type=Path, required=True)
    parser.add_argument("--registry-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--release-summary-out", type=Path, required=True)
    parser.add_argument(
        "--allow-rederived-eligibility",
        action="store_true",
        help=(
            "Accept a hash-bound corrected BLS-derived eligibility summary "
            "without imposing the superseded v2 19-row production lock."
        ),
    )
    args = parser.parse_args()
    result = write_full_pool_native_registry(
        pool_path=args.frozen_pool,
        pool_summary_path=args.frozen_pool_summary,
        eligibility_exclusions_path=args.eligibility_exclusions,
        eligibility_summary_path=args.eligibility_summary,
        native_shard_paths=args.native_shards,
        native_shard_summary_paths=args.native_shard_summaries,
        source_path=args.source_out,
        registry_path=args.registry_out,
        summary_path=args.summary_out,
        release_summary_path=args.release_summary_out,
        eligibility_production_lock=not args.allow_rederived_eligibility,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
