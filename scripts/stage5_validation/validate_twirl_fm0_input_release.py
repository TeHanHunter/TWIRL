#!/usr/bin/env python3
"""Independently validate every shard in an immutable FM0.1 input release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.input_release import validate_scientific_input_release
from twirl.models.fm0.registry import load_frozen_contract


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--freeze-receipt", type=Path)
    args = parser.parse_args()
    contract = load_frozen_contract(
        design_path=args.design,
        config_path=args.config,
        freeze_receipt_path=args.freeze_receipt,
    )
    result = validate_scientific_input_release(
        args.release_dir, contract=contract, workers=args.workers
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
