#!/usr/bin/env python3
"""Seed Teacher-v2 native-input verification cache from merge summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.teacher_v2_training import (
    seed_teacher_v2_native_verification_cache,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--native-summary",
        action="append",
        nargs=2,
        required=True,
        metavar=("NATIVE_H5", "SUMMARY_JSON"),
        type=Path,
    )
    parser.add_argument("--out-cache", type=Path, required=True)
    args = parser.parse_args()

    cache_path = seed_teacher_v2_native_verification_cache(
        args.out_cache,
        [(native, summary) for native, summary in args.native_summary],
    )
    payload = json.loads(cache_path.read_text())
    print(
        json.dumps(
            {
                "cache_path": str(cache_path.resolve()),
                "n_entries": len(payload.get("entries", {})),
                "schema": payload.get("schema", ""),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
