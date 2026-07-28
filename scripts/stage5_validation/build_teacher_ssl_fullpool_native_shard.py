#!/usr/bin/env python3
"""Build one checksum-bound real-only full-pool native HDF5 shard."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    build_full_pool_native_shard,
)


def _publish_immutable_json(path: Path, value: object) -> None:
    """Fsync and atomically publish JSON without overwriting another writer."""

    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--frozen-pool", type=Path, required=True)
    parser.add_argument("--frozen-pool-summary", type=Path, required=True)
    parser.add_argument("--eligibility-exclusions", type=Path, required=True)
    parser.add_argument("--eligibility-summary", type=Path, required=True)
    parser.add_argument("--sector-allowlist", type=Path, required=True)
    parser.add_argument("--raw-source-h5", type=Path, required=True)
    parser.add_argument("--raw-source-summary", type=Path, required=True)
    parser.add_argument("--raw-export-complete", type=Path, required=True)
    parser.add_argument(
        "--raw-transfer-validation",
        type=Path,
        required=True,
    )
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
        eligibility_exclusions_path=args.eligibility_exclusions,
        eligibility_summary_path=args.eligibility_summary,
        allowlist_path=args.sector_allowlist,
        raw_source_h5=args.raw_source_h5,
        raw_source_summary_path=args.raw_source_summary,
        raw_export_complete_path=args.raw_export_complete,
        raw_transfer_validation_path=args.raw_transfer_validation,
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
    _publish_immutable_json(summary_path, summary)
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
