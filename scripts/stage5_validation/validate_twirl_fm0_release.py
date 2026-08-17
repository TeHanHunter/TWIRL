#!/usr/bin/env python3
"""Independently validate one synthetic or real TWIRL-FM0.1 run release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.validation import (  # noqa: E402
    REAL_RUN_CONTRACT_SCHEMA_VERSION,
    RUN_CONTRACT_SCHEMA_VERSION,
    read_json,
    validate_real_run_release,
    validate_run_release,
    write_json_with_sha256,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    contract = read_json(args.run_dir / "run_contract.json")
    schema = contract.get("schema_version")
    if schema == REAL_RUN_CONTRACT_SCHEMA_VERSION:
        result = validate_real_run_release(args.run_dir)
    elif schema == RUN_CONTRACT_SCHEMA_VERSION:
        result = validate_run_release(args.run_dir)
    else:
        raise ValueError(f"unsupported FM0 run-contract schema: {schema!r}")
    if args.output is not None:
        write_json_with_sha256(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
