#!/usr/bin/env python3
"""Build, validate, or attach a Teacher-v1 native-input registry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION
from twirl.vetting.teacher_native_registry import (
    validate_native_input_registry_path,
    write_native_input_registry,
    write_native_registry_attachment,
)


def _add_expected_contract(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--expected-native-contract",
        default=RAW_PAIR_CONTRACT_VERSION,
        help=(
            "Required HDF5 root contract_version. Pass an empty string only "
            "to audit and retain heterogeneous nonblank contracts."
        ),
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser(
        "build",
        help="Build a unique (sector,TIC) registry from a CSV/Parquet table.",
    )
    build.add_argument("--source-table", type=Path, required=True)
    build.add_argument("--out-registry", type=Path, required=True)
    build.add_argument("--out-summary", type=Path, required=True)
    build.add_argument(
        "--path-base",
        type=Path,
        help="Base for relative native_h5_path values (default: source parent).",
    )
    _add_expected_contract(build)

    validate = subparsers.add_parser(
        "validate",
        help="Validate registry identities, files, hashes, contracts, and groups.",
    )
    validate.add_argument("--registry", type=Path, required=True)
    validate.add_argument(
        "--summary",
        type=Path,
        help="Optional frozen registry summary whose hashes/counts must match.",
    )
    _add_expected_contract(validate)

    attach = subparsers.add_parser(
        "attach",
        help="Attach a registry to a corpus using the exact (sector,TIC) key.",
    )
    attach.add_argument("--corpus", type=Path, required=True)
    attach.add_argument("--registry", type=Path, required=True)
    attach.add_argument(
        "--registry-summary",
        type=Path,
        help="Optional frozen registry summary to verify before attachment.",
    )
    attach.add_argument("--out-table", type=Path, required=True)
    attach.add_argument("--out-summary", type=Path, required=True)
    _add_expected_contract(attach)
    return parser


def _expected_contract(value: str) -> str | None:
    text = str(value).strip()
    return text or None


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    expected_contract = _expected_contract(args.expected_native_contract)
    if args.command == "build":
        result = write_native_input_registry(
            source_path=args.source_table,
            registry_path=args.out_registry,
            summary_path=args.out_summary,
            path_base=args.path_base,
            expected_contract_version=expected_contract,
        )
    elif args.command == "validate":
        result = validate_native_input_registry_path(
            registry_path=args.registry,
            summary_path=args.summary,
            expected_contract_version=expected_contract,
        )
    elif args.command == "attach":
        result = write_native_registry_attachment(
            corpus_path=args.corpus,
            registry_path=args.registry,
            registry_summary_path=args.registry_summary,
            output_path=args.out_table,
            summary_path=args.out_summary,
            expected_contract_version=expected_contract,
        )
    else:  # pragma: no cover - argparse enforces the subcommand.
        raise AssertionError(args.command)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
