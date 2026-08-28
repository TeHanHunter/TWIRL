#!/usr/bin/env python3
"""Freeze Phase-A authorities or run final S66--S77 admission-v2."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from twirl.models.fm0.preparation_authority import (
    admit_from_preparation_authorities,
    freeze_phase_a_authorities,
    load_phase_a_authority_record,
    load_preparation_authority_map,
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


def _sha256(value: str) -> str:
    if _SHA256.fullmatch(value) is None:
        raise argparse.ArgumentTypeError("SHA-256 must be lowercase 64-hex")
    return value


def _git_sha(value: str) -> str:
    if _GIT_SHA.fullmatch(value) is None:
        raise argparse.ArgumentTypeError("Git SHA must be lowercase 40-hex")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    for action in ("freeze-phase-a", "validate-phase-a", "admit"):
        child = subparsers.add_parser(action)
        child.add_argument("--authority-map", type=Path, required=True)
        child.add_argument("--authority-map-sha256", type=_sha256, required=True)
        child.add_argument("--producer-git-sha", type=_git_sha, required=True)
        if action != "validate-phase-a":
            child.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    authority = load_preparation_authority_map(
        args.authority_map,
        expected_sha256=args.authority_map_sha256,
    )
    if args.action == "freeze-phase-a":
        result, digest = freeze_phase_a_authorities(
            authority,
            producer_git_sha=args.producer_git_sha,
            output_path=args.output,
        )
    elif args.action == "validate-phase-a":
        result, path, digest = load_phase_a_authority_record(
            authority,
            expected_producer_git_sha=args.producer_git_sha,
        )
        result = {**result, "validated_path": str(path)}
    else:
        result, digest = admit_from_preparation_authorities(
            authority,
            producer_git_sha=args.producer_git_sha,
            output_path=args.output,
        )
    print(
        json.dumps(
            {"action": args.action, "sha256": digest, "result": result},
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
